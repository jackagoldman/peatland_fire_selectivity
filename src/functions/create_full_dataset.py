#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import yaml
import polars as pl
import geopandas as gpd
from datetime import date




# ---------- Read in all selectivity files using polars and combine into one dataset  ----------
def concat_dataframes(
    input_dir: Path,
    output_dir: Path | None = None,
    output_filename: str = "full_dataset.parquet",
    save_output: bool = True
) -> pl.DataFrame:

    all_files = list(input_dir.glob("*.parquet"))
    if not all_files:
        print(f"No parquet files found in {input_dir}")
        return None

    datasets = [pl.read_parquet(f) for f in all_files]
    full_dataset = pl.concat(datasets, how="vertical_relaxed")

    if save_output:
        if output_dir is None:
            raise ValueError("output_dir must be provided when save_output=True")
        output_path = output_dir / output_filename
        full_dataset.write_parquet(output_path)
        print(f"Full dataset created at {output_path}")

    return full_dataset

# ---------- get east vs west per progression centroid ----------
def get_east_west_from_centroids(
    input_file: Path,
    output_dir: Path | None = None,
    save_output: bool = True,
    output_filename: str = "progression_east_west.parquet"
) -> pl.DataFrame:

    prog_gdf = gpd.read_file(input_file).to_crs(4326)
    prog_gdf["x"] = prog_gdf.centroid.x

    ont_border = -95.15
    prog_gdf["region"] = prog_gdf["x"].apply(lambda x: "east" if x > ont_border else "west")

    result_df = pl.from_pandas(prog_gdf[["CLUSTERID", "region"]])

    result_df = result_df.with_columns(
    pl.col("CLUSTERID").cast(pl.Int64))



    if save_output:
        if output_dir is None:
            raise ValueError("output_dir must be provided when save_output=True")
        out_path = output_dir / output_filename
        result_df.write_parquet(out_path)
        print(f"East/west file written to {out_path}")

    return result_df


 # ---------- filter out land class when availability is below 5% ----------
 # proportion available is av

def filter_low_availability_classes(
    input_df: pl.DataFrame,
    output_dir: Path | None = None,
    threshold: float = 0.05,
    save_output: bool = True,
    output_filename: str = "filtered_full_dataset.parquet"
) -> pl.DataFrame:

    filtered_df = input_df.filter(pl.col("av") >= threshold)

    if save_output:
        if output_dir is None:
            raise ValueError("output_dir must be provided when save_output=True")
        out_path = output_dir / output_filename
        filtered_df.write_parquet(out_path)
        print(f"Filtered dataset created at {out_path}")

    return filtered_df

# ---------- add FWI data ----------
def merge_fwi_data(
    input_df: pl.DataFrame,
    fwi_file: Path,
    output_dir: Path | None = None,
    save_output: bool = True,
    output_filename: str = "full_dataset_with_fwi.parquet"
) -> pl.DataFrame:

    fwi_df = pl.read_csv(fwi_file)

    merged_df = input_df.join(fwi_df, on="K_UniqueID", how="left")


    def cast_fwi_to_numeric(df: pl.DataFrame) -> pl.DataFrame:
        numeric_cols = ["DC.mean", "ISI.mean", "BUI.mean", "FWI.mean"]
        # cast numeric FWI columns and write them to uniquely named columns
        # to avoid duplicate column-name collisions after joins
        return df.with_columns([
            pl.col(c).cast(pl.Float64, strict=False).alias(f"{c}_num")
            for c in numeric_cols
        ])
    merged_df = cast_fwi_to_numeric(merged_df)


    if save_output:
        if output_dir is None:
            raise ValueError("output_dir must be provided when save_output=True")
        out_path = output_dir / output_filename
        merged_df.write_parquet(out_path)
        print(f"FWI dataset created at {out_path}")

    return merged_df

# ---------- classify FWI into bins based on values ----------
"""
DC bins: 

low: 0-140
Moderate: 141-240
High: 241-340
Extreme: > 340

ISI bins:

Low: 0-2.2
Moderate: 2.3-5.0
High: 5.1-10.0
Extreme: > 10.0

BUI bins:

Low: 0-20
Moderate: 21-36
High: 37-60
Extreme: > 60

FWI bins:

Low: 0-3
Moderate: 4-10
High: 11-22
Extreme: > 22

"""
def classify_fwi_bins(
    input_df: pl.DataFrame,
    output_dir: Path | None = None,
    save_output: bool = True,
    output_filename: str = "full_dataset_with_fwi_bins_only.parquet"
) -> pl.DataFrame:


    # verify numeric columns exist
    print({c: input_df.schema.get(f"{c}_num") for c in ["DC.mean", "ISI.mean", "BUI.mean", "FWI.mean"]})

    classified_df = input_df.with_columns([

        # ---------- DC ----------
        pl.when(pl.col("DC.mean_num") <= 140).then(pl.lit("low"))
        .when(pl.col("DC.mean_num") <= 240).then(pl.lit("moderate"))
        .when(pl.col("DC.mean_num") <= 340).then(pl.lit("high"))
        .otherwise(pl.lit("extreme"))
        .alias("DC_bin"),

        # ---------- ISI ----------
        pl.when(pl.col("ISI.mean_num") <= 2.2).then(pl.lit("low"))
        .when(pl.col("ISI.mean_num") <= 5.0).then(pl.lit("moderate"))
        .when(pl.col("ISI.mean_num") <= 10.0).then(pl.lit("high"))
        .otherwise(pl.lit("extreme"))
        .alias("ISI_bin"),

        # ---------- BUI ----------
        pl.when(pl.col("BUI.mean_num") <= 20).then(pl.lit("low"))
        .when(pl.col("BUI.mean_num") <= 36).then(pl.lit("moderate"))
        .when(pl.col("BUI.mean_num") <= 60).then(pl.lit("high"))
        .otherwise(pl.lit("extreme"))
        .alias("BUI_bin"),

        # ---------- FWI ----------
        pl.when(pl.col("FWI.mean_num") <= 3).then(pl.lit("low"))
        .when(pl.col("FWI.mean_num") <= 10).then(pl.lit("moderate"))
        .when(pl.col("FWI.mean_num") <= 22).then(pl.lit("high"))
        .otherwise(pl.lit("extreme"))
        .alias("FWI_bin"),
    ])

    if save_output:
        if output_dir is None:
            raise ValueError("output_dir must be provided when save_output=True")
        out_path = output_dir / output_filename
        classified_df.write_parquet(out_path)
        print(f"FWI bin dataset written to {out_path}")

    return classified_df


# ---------- classify FWI based on quartile bins ----------

""" 
Create quartile bins for each FWI component based on the distribution of values in the dataset.
First, create q_bin_{fwi component} columns for each FWI component based on the quartiles of the distribution of values in the dataset (no na in quartile calc)
Using q_bin{fwi component} columns, create custom bin labels based on quartile bounrdaries: 
Q1: Min value - 0.25 quantile
Q2: 0.25 - 0.5 quantile
Q3: 0.5 - 0.75 quantile
Q4: 0.75 quantile - Max value
Assign custom bin labels to q_bin_{fwi_component} columns:

Then, write the resulting dataset to a new parquet file.
"""
def classify_fwi_quartile_bins(
    input_df: pl.DataFrame,
    output_dir: Path | None = None,
    save_output: bool = True,
    output_filename: str = "full_dataset_with_fwi_quartile_bins.parquet"
) -> pl.DataFrame:

    # use numeric columns created in merge_fwi_data
    cols = ["DC.mean_num", "ISI.mean_num", "BUI.mean_num", "FWI.mean_num"]

    def add_quartile_bins(df: pl.DataFrame, col: str) -> pl.DataFrame:
        # Compute quartiles on non-null values and alias them to avoid
        # duplicate output names in the select projection.
        q_df = (
            df.filter(pl.col(col).is_not_null())
            .select(
                pl.col(col).quantile(0.00).alias("q0"),
                pl.col(col).quantile(0.25).alias("q25"),
                pl.col(col).quantile(0.50).alias("q50"),
                pl.col(col).quantile(0.75).alias("q75"),
                pl.col(col).quantile(1.00).alias("q100"),
            )
        )

        q_dict = q_df.to_dict(as_series=False)
        # extract single values (quantile returns a single-element list)
        q0 = q_dict.get("q0", [None])[0]
        q25 = q_dict.get("q25", [None])[0]
        q50 = q_dict.get("q50", [None])[0]
        q75 = q_dict.get("q75", [None])[0]
        q100 = q_dict.get("q100", [None])[0]

        # If any quantile is None (e.g. all values null), skip binning
        if None in (q0, q25, q50, q75, q100):
            return df.with_columns(
                pl.lit(None).alias(f"{col}_bin"),
            )

        # Bin numeric quartiles 1–4
        # ensure strictly increasing bin edges (avoid equal edges)
        bins = [float(q0), float(q25), float(q50), float(q75), float(q100)]
        for i in range(1, len(bins)):
            if bins[i] <= bins[i - 1]:
                bins[i] = bins[i - 1] + 1e-9

        # Use when/then chain instead of pl.cut for compatibility
        df = df.with_columns(
            pl.when(pl.col(col) <= bins[1]).then(pl.lit("1"))
            .when(pl.col(col) <= bins[2]).then(pl.lit("2"))
            .when(pl.col(col) <= bins[3]).then(pl.lit("3"))
            .otherwise(pl.lit("4"))
            .alias(f"{col}_bin")
        )

        # Assign human‑readable quartile labels
        # create readable labels; remove trailing '_num' from column name for display
        disp_name = col.replace("_num", "")
        labels = {
            "1": f"{disp_name} Q1: {round(q0)}–{round(q25)}",
            "2": f"{disp_name} Q2: {round(q25)}–{round(q50)}",
            "3": f"{disp_name} Q3: {round(q50)}–{round(q75)}",
            "4": f"{disp_name} Q4: {round(q75)}–{round(q100)}",
        }

        df = df.with_columns(
            pl.when(pl.col(f"{col}_bin") == "1").then(pl.lit(labels["1"]))
            .when(pl.col(f"{col}_bin") == "2").then(pl.lit(labels["2"]))
            .when(pl.col(f"{col}_bin") == "3").then(pl.lit(labels["3"]))
            .when(pl.col(f"{col}_bin") == "4").then(pl.lit(labels["4"]))
            .otherwise(None)
            .alias(f"{col}_Quartile_Label")
        )
        return df

    # Apply binning to each FWI component
    df_out = input_df
    for col in cols:
        df_out = add_quartile_bins(df_out, col)

    # Save if requested
    if save_output:
        if output_dir is None:
            raise ValueError("output_dir must be provided when save_output=True")
        out_path = output_dir / output_filename
        df_out.write_parquet(out_path)
        print(f"Quartile FWI dataset written to {out_path}")

    return df_out


# ---------- Canopy and class bin classification ----------
""" create two new columns: 
canopy_bin
class_bin

canopy levels = open, forested, treed
if lc_class_name contains open then canopy_bin = "open"
if lc_class_name contains forested then canopy_bin = "forested"
if lc_class_name contains treed then canopy_bin = "treed"

class levels = bog, rich fen, poor fen, permafrost, upland
if lc_class_name contains "bog" then class_bin = "bog"
if lc_class_name contains "rich fen" then class_bin = "rich fen"
if lc_class_name contains "poor fen" then class_bin = "poor fen"
if lc_class_name contains "permafrost" then class_bin = "permafrost"
if lc_class_name contains "upland" then class_bin = "upland"

then factor the columns and write the resulting dataset to a new parquet file.
"""
def classify_canopy_and_class_bins(
    input_df: pl.DataFrame,
    output_dir: Path | None = None,
    save_output: bool = True,
    output_filename: str = "full_dataset_with_canopy_and_class_bins.parquet"
) -> pl.DataFrame:

    df_out = (
        input_df
        # canopy
        .with_columns(
            pl.when(pl.col("lc_class_name").str.contains("Open")).then(pl.lit("open"))
            .when(pl.col("lc_class_name").str.contains("Forested")).then(pl.lit("forested"))
            .when(pl.col("lc_class_name").str.contains("Treed")).then(pl.lit("treed"))
            .otherwise(None)
            .alias("canopy_bin")
        )
        # class grouping
        .with_columns(
            pl.when(pl.col("lc_class_name").str.contains("Bog")).then(pl.lit("bog"))
            .when(pl.col("lc_class_name").str.contains("Rich Fen")).then(pl.lit("rich fen"))
            .when(pl.col("lc_class_name").str.contains("Poor Fen")).then(pl.lit("poor fen"))
            .when(pl.col("lc_class_name").str.contains("PPC")).then(pl.lit("permafrost"))
            .when(pl.col("lc_class_name").str.contains("Upland")).then(pl.lit("upland"))
            .otherwise(None)
            .alias("class_bin")
        )
        .with_columns(
            pl.col("canopy_bin").cast(pl.Categorical),
            pl.col("class_bin").cast(pl.Categorical),
        )
    )

    # Save if needed
    if save_output:
        if output_dir is None:
            raise ValueError("output_dir must be provided when save_output=True")
        out_path = output_dir / output_filename
        df_out.write_parquet(out_path)
        print(f"Canopy/class bin dataset written to {out_path}")

    return df_out

#----- classify season based on Date column ------
"""
if Date<20230616 = Early (May to mid-June) else Late (mid-June to August)
"""
def classify_season_bins(
    input_df: pl.DataFrame,
    output_dir: Path | None = None,
    save_output: bool = True,
    output_filename: str = "fire_selectivity_full_dataset"
) -> pl.DataFrame:

    # Cast `Date` to integer if possible and compare to integer threshold
    date_threshold = 20230616
    date_expr = pl.col("Date").cast(pl.Int64, strict=False)

    df_out = input_df.with_columns(
        pl.when(date_expr < pl.lit(date_threshold)).then(pl.lit("Early"))
        .otherwise(pl.lit("Late"))
        .alias("season_bin")
    ).with_columns(
        pl.col("season_bin").cast(pl.Categorical)
    )

    if save_output:
        if output_dir is None:
            raise ValueError("output_dir must be provided when save_output=True")
        date_str = date.today().strftime("%Y%m%d")
        out_path = output_dir / f"{output_filename}_{date_str}.parquet"
        df_out.write_parquet(out_path)
        print(f"Season‑binned dataset written to {out_path}")

    return df_out

