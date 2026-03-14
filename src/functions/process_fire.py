# run_per_fire_parallel.py
import os
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Tuple, List, Sequence
import random
import yaml

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio import windows
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling, transform_bounds
from shapely.geometry import Polygon
from pyproj import Transformer

# -------------------------
# CONFIG (edit these)
# -------------------------

def load_config(filepath: Path) -> dict:
    with open(filepath, "r") as f:
        return yaml.safe_load(f)

# Load config relative to this file (adjust if needed)
SETTINGS = load_config(Path(__file__).parent.parent.parent / "config.yml")

PEAT_RASTER  = SETTINGS["python"]["peatland_raster"]
DNBR_DIR  = SETTINGS["python"]["dnbr_rasters"]
PROGRESSION_PATH = SETTINGS["python"]["progression_2_output"]
LANDSCAPE_PATH = SETTINGS["python"]["progression"]
OUT_DIR = SETTINGS["python"]["new_ua_processed"]  # new folder for Parquet outputs

DNBR_THRESHOLD = 0.1
ALL_TOUCHED = False             # rasterize option
WRITE_BATCH_SIZE = 1_000_000    # rows per Parquet shard
MAX_PROCESSES = max(1, os.cpu_count() - 1)   # cap as you like

# (Optional) Arrow/Parquet perf and GDAL cache
os.environ.setdefault("GDAL_CACHEMAX", "2048")
os.environ.setdefault("GDAL_NUM_THREADS", "ALL_CPUS")
# If installed: faster vector I/O
try:
    gpd.options.io_engine = "pyogrio"
except Exception:
    pass


# -------------------------
# UTILS
# -------------------------
def dnbr_path_for_fire(fire_id: str) -> Optional[str]:
    p = os.path.join(DNBR_DIR, f"{fire_id}_dnbr.tif")
    return p if os.path.exists(p) else None


def get_available_fire_ids(dnbr_dir: str) -> List[str]:
    return sorted({
        f[:-9] for f in os.listdir(dnbr_dir)
        if f.lower().endswith("_dnbr.tif")
    })


def intersect_bounds(a: Tuple[float, float, float, float],
                     b: Tuple[float, float, float, float]) -> Optional[Tuple[float, float, float, float]]:
    minx = max(a[0], b[0]); miny = max(a[1], b[1])
    maxx = min(a[2], b[2]); maxy = min(a[3], b[3])
    return None if (minx >= maxx or miny >= maxy) else (minx, miny, maxx, maxy)


def build_peat_window(peat_src, geom) -> Tuple[windows.Window, rasterio.Affine]:
    """Tight peat window covering geom (geom must be in peat CRS)."""
    bounds = gpd.GeoSeries([geom], crs=peat_src.crs).total_bounds
    win = windows.from_bounds(*bounds, transform=peat_src.transform,
                              width=peat_src.width, height=peat_src.height)
    pad = 1
    win = windows.Window(
        max(0, int(np.floor(win.col_off)) - pad),
        max(0, int(np.floor(win.row_off)) - pad),
        min(peat_src.width  - int(np.floor(win.col_off)) + pad, int(np.ceil(win.width))  + 2*pad),
        min(peat_src.height - int(np.floor(win.row_off)) + pad, int(np.ceil(win.height)) + 2*pad),
    )
    win = windows.Window(int(win.col_off), int(win.row_off), int(win.width), int(win.height))
    return win, windows.transform(win, peat_src.transform)


def rasterize_mask(geom, shape, transform) -> np.ndarray:
    if geom is None or geom.is_empty:
        return np.zeros(shape, dtype=np.uint8)
    return rasterize(
        [(geom, 1)],
        out_shape=shape,
        transform=transform,
        fill=0,
        dtype=np.uint8,
        all_touched=ALL_TOUCHED
    )


# -------------------------
# WORKER (one fire per process)
# -------------------------
from rasterio.warp import calculate_default_transform, transform_bounds

def process_one_fire(
    fire_id: str,
    progression_path: str,
    landscape_path: str,
    out_dir: str,
    dnbr_threshold: float = 0.1,
) -> str:
    """
    Process all polygons for a single fire ID:
      * Align work to peat grid
      * Reproject dNBR per peat-window
      * Write Parquet shards per CLUSTERID

    Adds:
      - area_factor_m2_per_pixel : scalar area (m²) of one pixel when the peat window
                                   grid is projected to EPSG:3347 (Lambert). 
                                   area(m²) = sum(used == 1) * area_factor_m2_per_pixel
    """
    dnbr_path = dnbr_path_for_fire(fire_id)
    if not dnbr_path:
        return f"skip: no dNBR for FireID={fire_id}"

    lambert_stats_can_epsg = 3347  # same as your R code

    # Open peat and dNBR once in this process
    with rasterio.open(PEAT_RASTER) as peat_src, rasterio.open(dnbr_path) as dnbr_ds:
        peat_crs = peat_src.crs
        peat_transform = peat_src.transform
        peat_nodata = peat_src.nodata

        # Read vectors, projected to peat CRS
        prog = gpd.read_file(progression_path).to_crs(peat_crs)
        landscape = gpd.read_file(landscape_path).to_crs(peat_crs)

        # Filter to this fire
        if "K_FireID" not in prog.columns:
            return f"error: progression missing K_FireID"
        prog["K_FireID_str"] = prog["K_FireID"].astype(str)
        prog = prog[prog["K_FireID_str"] == fire_id].copy()
        if prog.empty:
            return f"skip: no polygons in progression for FireID={fire_id}"

        # Identify the K_UniqueID column name if present (case variations)
        kunique_col = None
        for cand in ("K_UniqueID", "K_UNIQUEID", "K_uniqueid"):
            if cand in prog.columns:
                kunique_col = cand
                break

        # Dissolve landscape by CLUSTERID (ensure one geometry per id)
        if "CLUSTERID" not in prog.columns:
            return f"error: progression missing CLUSTERID"
        if "CLUSTERID" not in landscape.columns:
            return f"error: landscape missing CLUSTERID"

        landscape = landscape.dissolve(by="CLUSTERID").reset_index()
        land_map = {int(r.CLUSTERID): r.geometry for _, r in landscape.iterrows()}

        # Loop over polygons for this fire
        for _, row in prog.iterrows():
            cluster_id = int(row["CLUSTERID"])
            kunique_id = int(row[kunique_col]) if kunique_col and pd.notnull(row[kunique_col]) else None
            date_val = row.get("DATE", None)
            poly_geom = row.geometry
            land_geom = land_map.get(cluster_id)
            if land_geom is None:
                warnings.warn(f"[{fire_id}] missing landscape geom for CLUSTERID={cluster_id}; skip")
                continue

            # ---- Build peat window for landscape geom (availability universe)
            win, win_transform = build_peat_window(peat_src, land_geom)
            if win.width <= 0 or win.height <= 0:
                warnings.warn(f"[{fire_id}] empty peat window for CLUSTERID={cluster_id}; skip")
                continue
            out_shape = (int(win.height), int(win.width))

            # Peat class window (categorical; no resample)
            lc_arr = peat_src.read(1, window=win)

            # ---- Reproject the minimal dNBR source window into peat window grid
            # 1) bounds in peat CRS for dst window
            win_bounds = windows.bounds(win, transform=peat_transform)  # (minx, miny, maxx, maxy)
            # 2) transform dst bounds -> dNBR CRS to compute source window
            try:
                src_bounds = transform_bounds(peat_crs, dnbr_ds.crs, *win_bounds, densify_pts=21)
            except Exception as e:
                warnings.warn(f"[{fire_id}] bounds transform failed for CLUSTERID={cluster_id}: {e}")
                continue

            src_win = windows.from_bounds(*src_bounds, transform=dnbr_ds.transform,
                                          width=dnbr_ds.width, height=dnbr_ds.height)
            # Clip to dataset
            src_col_off = max(0, int(np.floor(src_win.col_off)))
            src_row_off = max(0, int(np.floor(src_win.row_off)))
            src_w = min(dnbr_ds.width  - src_col_off, int(np.ceil(src_win.width)))
            src_h = min(dnbr_ds.height - src_row_off, int(np.ceil(src_win.height)))
            if src_w <= 0 or src_h <= 0:
                warnings.warn(f"[{fire_id}] no dNBR overlap for CLUSTERID={cluster_id}; skip")
                continue
            src_win = windows.Window(src_col_off, src_row_off, src_w, src_h)
            src_win_transform = windows.transform(src_win, dnbr_ds.transform)

            # Read dNBR source window
            dnbr_src_arr = dnbr_ds.read(1, window=src_win, boundless=True)

            # Reproject to peat window grid
            dnbr_aligned = np.zeros(out_shape, dtype=np.float32)
            reproject(
                source=dnbr_src_arr,
                destination=dnbr_aligned,
                src_transform=src_win_transform,
                src_crs=dnbr_ds.crs,
                dst_transform=win_transform,
                dst_crs=peat_crs,
                resampling=Resampling.bilinear,  # continuous index
                num_threads=2
            )

            # ---- Compute area factor (m²/pixel) by projecting the peat window grid to EPSG:3347
            # This mirrors the R logic: project raster, then area = count * (res_x * res_y)
            dst_crs = rasterio.crs.CRS.from_epsg(lambert_stats_can_epsg)
            # We keep the same output array size (out_shape) so factor is uniform
            dst_transform, dst_width, dst_height = calculate_default_transform(
                peat_crs, dst_crs,
                width=out_shape[1], height=out_shape[0],
                left=win_bounds[0], bottom=win_bounds[1], right=win_bounds[2], top=win_bounds[3]
            )
            # Pixel area in projected CRS
            area_factor_m2_per_pixel = abs(dst_transform.a * dst_transform.e)

            # ---- Rasterize masks in peat window
            avail_mask = rasterize_mask(land_geom, out_shape, win_transform)  # 0/1
            poly_mask  = rasterize_mask(poly_geom,  out_shape, win_transform)  # 0/1

            burn_mask = (dnbr_aligned > dnbr_threshold).astype(np.uint8)
            used_mask = (burn_mask & (poly_mask == 1)).astype(np.uint8)

            # Keep pixels with valid lc_class that are in availability OR in the progression polygon
            valid_lc = (lc_arr != peat_nodata) if (peat_nodata is not None and not np.isnan(peat_nodata)) else np.isfinite(lc_arr)
            keep_mask = valid_lc & ((avail_mask == 1) | (poly_mask == 1))

            rows_idx, cols_idx = np.where(keep_mask)
            if rows_idx.size == 0:
                warnings.warn(f"[{fire_id}] no valid pixels for CLUSTERID={cluster_id}; skip")
                continue

            # Coordinates (peat CRS)
            xs, ys = rasterio.transform.xy(win_transform, rows_idx, cols_idx, offset="center")
            xs = np.asarray(xs); ys = np.asarray(ys)

            # Values
            lc_vals = lc_arr[rows_idx, cols_idx].astype(np.int32, copy=False)
            used_vals = used_mask[rows_idx, cols_idx].astype(np.uint8, copy=False)
            in_landscape_vals = (avail_mask[rows_idx, cols_idx] == 1).astype(np.uint8)

            # Write Parquet shards
            out_fire_cluster = Path(out_dir) / f"FireID={fire_id}" / f"CLUSTERID={cluster_id}"
            out_fire_cluster.mkdir(parents=True, exist_ok=True)

            n = len(rows_idx)
            start = 0
            part = 0
            while start < n:
                end = min(start + WRITE_BATCH_SIZE, n)
                df = pd.DataFrame({
                    "K_FireID": fire_id,
                    "CLUSTERID": cluster_id,
                    "K_UniqueID": kunique_id,
                    "DATE": date_val,
                    "x": xs[start:end],
                    "y": ys[start:end],
                    "lc_class": lc_vals[start:end],
                    "used": used_vals[start:end],
                    "in_landscape": in_landscape_vals[start:end],
                    "area_factor_m2_per_pixel": area_factor_m2_per_pixel,  # <-- NEW (scalar)
                })
                df.to_parquet(out_fire_cluster / f"part-{part:04d}.parquet", index=False)
                part += 1
                start = end

        return f"ok: FireID={fire_id}"
# -------------------------
# DRIVER (parallel per fire)
# -------------------------
from typing import Optional, Sequence
import random

def run_all_fires_parallel(
    progression_path: str,
    landscape_path: str,
    out_dir: str,
    max_processes: int = MAX_PROCESSES,
    dnbr_threshold: float = DNBR_THRESHOLD,
    n: Optional[int] = None,                 # <-- NEW: limit number of fires
    random_sample: bool = False,             # <-- NEW: choose fires randomly
    fires: Optional[Sequence[str]] = None,   # <-- NEW: explicit list of fire IDs to run
    dry_run: bool = False,                   # <-- NEW: list and exit
):
    """
    Discover available fires (from *_dnbr.tif) and run one process per fire.
    Use `n` to limit the number of fires for testing.

    Args:
        progression_path: path to progression (GeoJSON or Shapefile)
        landscape_path: path to landscape polygons (shapefile)
        out_dir: output directory for Parquet
        max_processes: number of worker processes
        dnbr_threshold: dNBR threshold for burned pixels
        n: if provided, run only N fires
        random_sample: if True, sample N fires randomly (requires n)
        fires: if provided, use exactly these fire IDs (strings matching <FireID>_dnbr.tif)
        dry_run: if True, print selected fires and return without processing
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # 1) Determine candidate fire IDs
    if fires is not None:
        # use user-specified list
        fire_ids = [str(fid) for fid in fires]
    else:
        # discover from dNBR filenames
        discovered = get_available_fire_ids(DNBR_DIR)  # e.g., ["66", "1234", ...]
        if not discovered:
            raise RuntimeError("No *_dnbr.tif found in DNBR_DIR")
        fire_ids = discovered

    # 2) Filter progression to keep only those with DNBR available
    prog = gpd.read_file(progression_path)
    if "K_FireID" not in prog.columns:
        raise ValueError("progression is missing K_FireID")
    prog["K_FireID_str"] = prog["K_FireID"].astype(str)

    fire_ids = sorted(set(prog["K_FireID_str"]) & set(fire_ids))
    if not fire_ids:
        raise RuntimeError("No overlap between progression K_FireID and dNBR filenames")

    # 3) Reduce to N fires if requested
    if n is not None:
        if n <= 0:
            raise ValueError("n must be a positive integer")
        if random_sample:
            random.seed(42)  # make deterministic if you like
            fire_ids = random.sample(fire_ids, k=min(n, len(fire_ids)))
        else:
            fire_ids = fire_ids[:n]

    # 4) Dry run?
    if dry_run:
        print(f"[DRY RUN] Would process {len(fire_ids)} fire(s): {fire_ids}")
        return

    # 5) Fan out: one process per fire
    with ProcessPoolExecutor(max_workers=max_processes) as ex:
        futures = [
            ex.submit(process_one_fire, fid, progression_path, landscape_path, out_dir, dnbr_threshold)
            for fid in fire_ids
        ]
        for fut in as_completed(futures):
            try:
                print(fut.result())
            except Exception as e:
                print("error:", e)


if __name__ == "__main__":
    run_all_fires_parallel(
        progression_path=PROGRESSION_PATH,
        landscape_path=LANDSCAPE_PATH,
        out_dir=OUT_DIR,
        max_processes=8,
        dnbr_threshold=0.1
    )
