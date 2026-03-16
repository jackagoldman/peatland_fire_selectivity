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
