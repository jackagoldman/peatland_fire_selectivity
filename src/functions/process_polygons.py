import os
import re
import warnings
from typing import Optional, Union, Iterable
from pathlib import Path


import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from shapely.geometry import Point
import yaml


# ---------- Configuration helpers ----------

def load_config(filepath: Path) -> dict:
    with open(filepath, "r") as f:
        return yaml.safe_load(f)

# Load config relative to this file (adjust if needed)
SETTINGS = load_config(Path(__file__).parent.parent.parent / "config.yml")

#----------- Main function to create prog_sub ----------
def make_prog_sub(
    progression_path: str,
    peat_raster_path: str,
    dnbr_folder: str,
    area_quantile: float = 0.20,      # keep polygons >= this quantile (top 80% by area)
    buffer_distance_m: float = 200.0, # buffer in meters
    working_crs: str = "EPSG:3347",   # NAD83 / LCC Canada
    fix_invalid_geoms: bool = True,
) -> gpd.GeoDataFrame:
    """
    Rebuild `filtered_polygons` and then subset to polygons with available dNBR rasters,
    mirroring the R steps from NickPel.

    Returns a GeoDataFrame in `working_crs` (default EPSG:3347) with:
        - filtered + buffered polygons
        - only polygons intersecting peat raster coverage (via centroid test)
        - only polygons whose K_FireID has a matching <FireID>_dnbr.tif in dnbr_folder
    """
    # ---------- 1) Load inputs ----------
    gdf = gpd.read_file(progression_path)
    if gdf.empty:
        raise ValueError("Progression shapefile contains no features.")

    # Make sure required attributes exist
    required_cols = ["K_FireID"]      # CLUSTERID is used later in your pipeline; warn if missing
    for col in required_cols:
        if col not in gdf.columns:
            raise ValueError(f"Required column '{col}' not found in {progression_path}")

    if "CLUSTERID" not in gdf.columns:
        warnings.warn("Column 'CLUSTERID' not found. Your later selectivity step will need it.")

    # Work in projected CRS for correct area/buffer ops
    gdf = gdf.to_crs(working_crs)

    # Optionally clean invalid geometries (self-intersections, etc.)
    if fix_invalid_geoms:
        # Classic fix: zero-width buffer
        gdf["geometry"] = gdf.buffer(0)

    # Drop empty/invalid
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.is_valid].copy()
    if gdf.empty:
        raise ValueError("All geometries are empty/invalid after cleaning.")

    # ---------- 2) Keep largest 80% by area ----------
    # Compute area in m² in working CRS
    gdf["AREA_CALC"] = gdf.geometry.area.astype("float64")
    thresh = gdf["AREA_CALC"].quantile(area_quantile)
    gdf = gdf[gdf["AREA_CALC"] >= thresh].copy()

    if gdf.empty:
        raise ValueError("No polygons remain after area-quantile filtering.")

    # ---------- 3) Buffer by 200 m (if requested) ----------
    if buffer_distance_m and buffer_distance_m > 0:
        gdf["geometry"] = gdf.geometry.buffer(buffer_distance_m)

    # Drop any new invalids created by buffer
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.is_valid].copy()
    if gdf.empty:
        raise ValueError("No polygons remain after buffering and validity checks.")

    # ---------- 4) Keep polygons within peat raster coverage (centroid sample) ----------
    with rasterio.open(peat_raster_path) as src:
        raster_crs = src.crs
        nodata = src.nodata

        # Compute centroids in raster CRS
        centroids = gdf.geometry.centroid
        centroids = centroids.to_crs(raster_crs)
        coords = [(pt.x, pt.y) for pt in centroids]

        # Sample raster at centroid locations
        # src.sample yields arrays of band values; assume single-band peat map here
        vals = [next(src.sample([xy]))[0] for xy in coords]

    vals = np.array(vals)
    # Create a mask of "valid peat value": not nan and not nodata (if nodata is set)
    if nodata is not None and not np.isnan(nodata):
        valid_mask = (~np.isnan(vals)) & (vals != nodata)
    else:
        valid_mask = ~np.isnan(vals)

    gdf = gdf.loc[valid_mask].copy()
    if gdf.empty:
        raise ValueError("No polygons overlap the peat raster (centroid test).")

    # ---------- 5) Subset to polygons with matching dNBR rasters ----------
    if not os.path.isdir(dnbr_folder):
        raise FileNotFoundError(f"dNBR folder not found: {dnbr_folder}")

    # Collect files that look like <FireID>_dnbr.tif (case-insensitive)
    dnbr_files = [f for f in os.listdir(dnbr_folder) if re.search(r"dnbr\.tif$", f, flags=re.IGNORECASE)]

    # Extract FireIDs from filenames
    available_ids = set()
    for fname in dnbr_files:
        m = re.match(r"(.+?)_dnbr\.tif$", fname, flags=re.IGNORECASE)
        if m:
            available_ids.add(m.group(1))

    if not available_ids:
        warnings.warn("No <FireID>_dnbr.tif files found. Returning empty subset.")
        return gdf.iloc[0:0].copy()

    # Compare as strings to avoid numeric/string mismatches
    gdf["K_FireID_str"] = gdf["K_FireID"].astype(str)
    prog_sub = gdf[gdf["K_FireID_str"].isin(available_ids)].copy()

    # ---------- Sanity checks ----------
    if prog_sub.empty:
        warnings.warn("No polygons matched available dNBR rasters after filtering.")
    else:
        print(f"Polygons after all filters: {len(prog_sub)}")
        print(f"Unique K_FireID in subset: {prog_sub['K_FireID_str'].nunique()}")

    # Keep in working CRS (EPSG:3347) for downstream area operations
    return prog_sub


