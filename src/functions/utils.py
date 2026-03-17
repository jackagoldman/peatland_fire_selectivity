
from __future__ import annotations

import os
import re
import warnings
from typing import Optional, Union, Iterable, List
from pathlib import Path
import polars as pl


import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from shapely.geometry import Point
import yaml
import shutil


# ---------- Configuration helpers ----------

def load_config(filepath: Path) -> dict:
    with open(filepath, "r") as f:
        return yaml.safe_load(f)

# Load config relative to this file (adjust if needed)
SETTINGS = load_config(Path(__file__).parent.parent.parent / "config.yml")


# -------- Collect cluster ids to parquet



def collect_cluster_parquets(root_dir):
    """
    Scans FireID=*/CLUSTERID=*/ folders inside `root_dir`,
    copies each .parquet file into root_dir/ua_clusterid/,
    and renames them as ua_<clusterid>.parquet.
    """
    root = Path(root_dir)
    out = root.parent / "ua_clusterid"
    out.mkdir(exist_ok=True)

    for fire_dir in root.glob("FireID=*"):
        for cluster_dir in fire_dir.glob("CLUSTERID=*"):
            cluster_id = cluster_dir.name.split("=")[1]

            # find parquet file
            for pq in cluster_dir.glob("*.parquet"):
                new_name = f"ua_{cluster_id}.parquet"
                shutil.copy2(pq, out / new_name)
                break  # use first parquet if multiple


#----------- calculate area of polygons from clusters  ----------

def cluster_areas(progression_path: str,
                  landscape_path: str,
                  out_path: str,
                  id_col: str = "CLUSTERID",
                  area_crs: str = "EPSG:3347") -> pd.DataFrame:
    """
    Compute per-CLUSTERID polygon areas from two shapefiles and return a
    merged DataFrame with columns:
        CLUSTERID, poly_area (progression), buffer_area (landscape)

    Parameters
    ----------
    progression_path : str
        Path to the progression shapefile (polygons).
    landscape_path : str
        Path to the landscape shapefile (polygons).
    out_path : str
        Path to write the output Parquet file with area calculations.
    id_col : str, default "CLUSTERID"
        The ID column present in both shapefiles.
    area_crs : str, default "EPSG:6933"
        Equal-area CRS used for area calculation if input is geographic.

    Returns
    -------
    pandas.DataFrame
        Columns: [CLUSTERID, poly_area, buffer_area]; areas in square meters.
    """

    def _prep_and_area(path, id_col, area_crs):
        gdf = gpd.read_file(path)

        if id_col not in gdf.columns:
            raise KeyError(f"'{id_col}' not found in {path} columns: {list(gdf.columns)}")

        # Drop rows without an ID or geometry
        gdf = gdf.dropna(subset=[id_col, gdf.geometry.name]).copy()

        if gdf.empty:
            return pd.DataFrame(columns=[id_col, "area_m2"]).astype({id_col: gdf.dtypes.get(id_col, "object")})

        # Ensure we have a CRS; if missing, you may need to set it manually before calling this function
        if gdf.crs is None:
            raise ValueError(f"CRS is missing in {path}. Set the CRS before area calculation.")

        # If geographic (degrees), project to equal-area for correct area computation
        if not gdf.crs.is_projected:
            gdf = gdf.to_crs(area_crs)

        # Compute area (m²) for polygons; non-polygons will yield 0
        gdf["area_m2"] = gdf.geometry.area.fillna(0)

        # Sum by ID
        return (gdf.groupby(id_col, dropna=False, as_index=False)["area_m2"]
                    .sum())

    prog = _prep_and_area(progression_path, id_col, area_crs).rename(columns={"area_m2": "poly_area"})
    land = _prep_and_area(landscape_path,  id_col, area_crs).rename(columns={"area_m2": "buffer_area"})

    # Outer join to keep all IDs present in either file
    out = pd.merge(prog, land, on=id_col, how="outer")

    # Replace NaNs with 0 for areas where an ID is missing in one source
    out[["poly_area", "buffer_area"]] = out[["poly_area", "buffer_area"]].fillna(0)

    # write out as parquet
    out.to_parquet(out_path, index=False)

    # Optional: sort by ID
    return out.sort_values(by=id_col).reset_index(drop=True)

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# ------ calculate area of burned pixels from cluster using factor column in consolidated fire files ----------
def cluster_burned_area(
    fire_parquet_path: str,
    out_path: str,
    id_col: str = "CLUSTERID",
    factor_col: str = "area_factor_m2_per_pixel",
    *,
    output_format: str = "int",      # "float" | "int" | "string" | "decimal128"
    decimal_scale: int = 0           # used when output_format="decimal128"
) -> pd.DataFrame:
    """
    Calculate burned area for each cluster by summing the area factor column
    across all pixels in the consolidated fire Parquet file.

    Parameters
    ----------
    fire_parquet_path : str
        Path to the consolidated fire Parquet file containing pixel-level data.
    id_col : str, default "CLUSTERID"
        The column name in the Parquet file that identifies clusters.
    factor_col : str, default "area_factor_m2_per_pixel"
        The column name that contains the area factor (m² per pixel).
    output_format : {"float","int","string","decimal128"}, default "int"
        How to store burned_area_m2 in the output parquet.
        - "float": keep float64 (may display with scientific notation in some tools)
        - "int": round to nearest meter and store as nullable Int64 (no sci-notation)
        - "string": store as string formatting (e.g., "1234567")
        - "decimal128": store as Decimal128 with the given scale (fixed point)
    decimal_scale : int, default 0
        Scale for Decimal128 (e.g., 0 for whole meters, 2 for centimeters).

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns [CLUSTERID, burned_area_m2] (in-memory dtypes reflect output_format).
    """
    # Read only needed columns
    df = pd.read_parquet(fire_parquet_path, columns=[id_col, factor_col])

    if id_col not in df.columns:
        raise KeyError(f"'{id_col}' not found in Parquet columns: {list(df.columns)}")
    if factor_col not in df.columns:
        raise KeyError(f"'{factor_col}' not found in Parquet columns: {list(df.columns)}")

    # Sum area factors by cluster ID
    burned_area = (
        df.groupby(id_col, dropna=False)[factor_col]
          .sum()
          .reset_index()
          .rename(columns={factor_col: "burned_area_m2"})
    )

    # ---- Normalize type / formatting for output ----
    fmt = output_format.lower()
    if fmt == "float":
        # Ensure float64
        burned_area["burned_area_m2"] = burned_area["burned_area_m2"].astype("float64")

        # Write as standard parquet via pandas (will be float)
        burned_area.to_parquet(out_path, index=False)

    elif fmt == "int":
        # Round to whole meters & use nullable Int64 to avoid float sci-notation
        burned_area["burned_area_m2"] = (
            burned_area["burned_area_m2"].round().astype("Int64")
        )
        burned_area.to_parquet(out_path, index=False)

    elif fmt == "string":
        # Convert to strings explicitly to avoid any scientific formatting downstream
        # You can apply custom formatting, e.g., commas:
        # burned_area["burned_area_m2"] = burned_area["burned_area_m2"].map(lambda x: f"{x:,.0f}")
        burned_area["burned_area_m2"] = burned_area["burned_area_m2"].map(
            lambda x: "" if pd.isna(x) else f"{x:.0f}"
        )
        burned_area.to_parquet(out_path, index=False)

    elif fmt == "decimal128":
        # Use pyarrow to write Decimal128 (fixed point) to parquet
        # Decide precision based on maximum digits; default to a safe high precision like 38
        # scale = decimal_scale (0 means whole meters)
        pa_tbl = pa.Table.from_pandas(burned_area, preserve_index=False)

        # Compute a Decimal128 type
        # Precision must be >= digits (before + after decimal); 38 is a common safe max for Arrow
        dec_type = pa.decimal128(38, decimal_scale)

        # Cast the burned_area_m2 column: first to string->decimal or float->decimal via scaling
        # Simplest: cast via string to avoid rounding surprises
        col = pa_tbl["burned_area_m2"]

        # If decimal_scale == 0 and col is float, formatting to string with no sci-notation helps
        col_as_str = pa.compute.cast(col, pa.large_string())
        col_dec = pa.compute.cast(col_as_str, dec_type)

        pa_tbl = pa_tbl.set_column(
            pa_tbl.column_names.index("burned_area_m2"),
            "burned_area_m2",
            col_dec
        )

        pq.write_table(pa_tbl, out_path, compression="zstd")

    else:
        raise ValueError("output_format must be one of: 'float', 'int', 'string', 'decimal128'")

    # Return a sorted DataFrame for convenience
    return burned_area.sort_values(by=id_col).reset_index(drop=True)




# merge cluster areas and burned area back to individual fire files 

def _align_key_dtype(lf: pl.LazyFrame, key: str, target_dtype: pl.DataType) -> pl.LazyFrame:
    """Cast join key to target dtype when present."""
    return lf.with_columns(pl.col(key).cast(target_dtype)) if key in lf.columns else lf


def merge_fire_file_inplace(
    fire_parquet_path: Union[str, Path],
    progression_landscape_path: Union[str, Path],
    burned_area_path: Union[str, Path],
    *,
    cluster_col: str = "CLUSTERID",
    poly_area_col: str = "poly_area",
    buffer_area_col: str = "buffer_area",
    burned_area_col: str = "burned_area_m2",
) -> None:
    """
    Read one K_FireID_*.parquet and LEFT JOIN:
      1) [poly_area, buffer_area] from progression_landscape_areas_by_cluster.parquet on CLUSTERID
      2) [burned_area_m2] from burned_area.parquet on CLUSTERID

    Overwrites the original file in place (same path, same name).
    - Left joins preserve the row count of the fire file; unmatched keys yield nulls.
    - The join key dtype is aligned to the fire file schema to avoid join errors.
    """
    fire_parquet_path = Path(fire_parquet_path)
    progression_landscape_path = Path(progression_landscape_path)
    burned_area_path = Path(burned_area_path)

    # Lazy scans
    lf_fire = pl.scan_parquet(fire_parquet_path)

    lf_prog_land = (
        pl.scan_parquet(progression_landscape_path)
          .select([cluster_col, poly_area_col, buffer_area_col])
    )

    lf_burned = (
        pl.scan_parquet(burned_area_path)
          .select([cluster_col, burned_area_col])
    )

    # Align join key dtypes to the fire file’s dtype (fallback Int64)
    fire_dtype = lf_fire.schema.get(cluster_col, pl.Int64)
    lf_fire      = _align_key_dtype(lf_fire, cluster_col, fire_dtype)
    lf_prog_land = _align_key_dtype(lf_prog_land, cluster_col, fire_dtype)
    lf_burned    = _align_key_dtype(lf_burned, cluster_col, fire_dtype)

    # Left joins preserve fire row count
    lf_out = (
        lf_fire.join(lf_prog_land, on=cluster_col, how="left")
               .join(lf_burned,    on=cluster_col, how="left")
    )

    # Collect and overwrite in place
    out_df = lf_out.collect()
    out_df.write_parquet(fire_parquet_path)  # overwrite same file
    print(f"[ok] Overwrote: {fire_parquet_path} (rows={out_df.height}, cols={out_df.width})")


def merge_fire_dir_inplace(
    fire_dir: Union[str, Path],
    progression_landscape_path: Union[str, Path],
    burned_area_path: Union[str, Path],
    *,
    recursive: bool = False,
    cluster_col: str = "CLUSTERID",
    poly_area_col: str = "poly_area",
    buffer_area_col: str = "buffer_area",
    burned_area_col: str = "burned_area_m2",
) -> List[Path]:
    """
    Process all K_FireID_*.parquet files under fire_dir (optionally recursively),
    performing the joins and overwriting each file in place.

    Returns the list of updated file paths.
    """
    fire_dir = Path(fire_dir)
    pattern = "**/K_FireID_*.parquet" if recursive else "K_FireID_*.parquet"
    files = sorted(fire_dir.glob(pattern))
    if not files:
        print(f"[warn] No files matching {pattern} in {fire_dir}")
        return []

    updated: List[Path] = []
    for f in files:
        try:
            merge_fire_file_inplace(
                fire_parquet_path=f,
                progression_landscape_path=progression_landscape_path,
                burned_area_path=burned_area_path,
                cluster_col=cluster_col,
                poly_area_col=poly_area_col,
                buffer_area_col=buffer_area_col,
                burned_area_col=burned_area_col,
            )
            updated.append(f)
        except Exception as e:
            print(f"[error] Failed to update {f}: {e}")

    print(f"[done] updated={len(updated)} in {fire_dir}")
    return updated