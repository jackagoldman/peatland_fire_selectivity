#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import yaml
import polars as pl
import geopandas as gpd



# ---------- Configuration helpers ----------

def load_config(filepath: Path) -> dict:
    with open(filepath, "r") as f:
        return yaml.safe_load(f)

# Load config relative to this file (adjust if needed)
SETTINGS = load_config(Path(__file__).parent.parent / "config.yml")

# Confirm keys exist
in_rel = SETTINGS["python"]["kunique_selectivity_output"]       # e.g., "data/ua_fire_combined/"
out_rel = SETTINGS["python"]["full_dataset_output"]  # e.g., "data/full_dataset/"
prog_rel = SETTINGS["python"]["progression"]  # e.g., "data/progression/"
fwi_rel = SETTINGS["python"]["fwi_input"]  # e.g., "data/fwi/fwi.parquet"


# ---------- Read in all selectivity files using polars and combine into one dataset  ----------
def create_full_dataset(input_dir: Path, output_dir: Path, output_filename: str = "full_dataset.parquet") -> pl.DataFrame:
    # Read all parquet files in the input directory
    all_files = list(input_dir.glob("*.parquet"))
    if not all_files:
        print(f"No parquet files found in {input_dir}")
        return

    # Load and concatenate all datasets
    datasets = [pl.read_parquet(file) for file in all_files]
    full_dataset = pl.concat(datasets, how="vertical_relaxed")

    # Write the combined dataset to a new parquet file
    output_path = output_dir / output_filename
    full_dataset.write_parquet(output_path)
    print(f"Full dataset created at {output_path}")
    return full_dataset

# ---------- get east vs west per progression centroid ----------
def get_east_west_from_centroids(input_file: Path, output_dir: Path) -> pl.DataFrame:
    # read in the progression shapefile
    prog_gdf = gpd.read_file(input_file).to_crs(4326)  # ensure it's in WGS84
    # get the longitude of the centroid for each progression polygon
    prog_gdf['x'] = prog_gdf.centroid.x
    # boundary of east vs west is Ontario's western border
    ontario_western_border = -95.15
    # assign east/west based on the longitude of the centroid
    prog_gdf['region'] = prog_gdf['x'].apply(lambda x: 'east' if x > ontario_western_border else 'west')
    # create a dataframe with kunqiue_id, CLUSTERID and region
    result_df = pl.from_pandas(prog_gdf[['CLUSTERID', 'region']])

    # write the result to a parquet file
    output_path = output_dir / "progression_east_west.parquet"
    result_df.write_parquet(output_path)

    return result_df


 # ---------- filter out land class when availability is below 5% ----------
 # proportion available is pi_avail

def filter_low_availability_classes(
    input_file: pl.DataFrame, 
    output_dir: Path, 
    threshold: float = 0.05
    ) -> pl.DataFrame:
    # filter out rows where pi_avail is below the threshold
    filtered_df = input_file.filter(pl.col("pi_avail") >= threshold)
    # write the filtered dataset to a new parquet file
    output_path = output_dir / "filtered_full_dataset.parquet"
    filtered_df.write_parquet(output_path)
    print(f"Filtered dataset created at {output_path}")
    return filtered_df

# ---------- add FWI data ----------
def merge_fwi_data(
    input_file: pl.DataFrame,
    fwi_file: Path,
    output_dir: Path
    ) -> pl.DataFrame:
    # read in the FWI data
    fwi_df = pl.read_csv(fwi_file)
    # merge the FWI data with the input dataset on CLUSTERID
    merged_df = input_file.join(fwi_df, on=["CLUSTERID"], how="left")
    # write the merged dataset to a new parquet file
    output_path = output_dir / "full_dataset_with_fwi.parquet"
    merged_df.write_parquet(output_path)
    print(f"Full dataset with FWI created at {output_path}")
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
    input_file: pl.DataFrame,
    output_dir: Path
    ) -> pl.DataFrame:
    # classify FWI components into bins based on the specified thresholds
    classified_df = input_file.with_columns(
        pl.when(pl.col("DC.mean") <= 140).then("low")
        .when((pl.col("DC.mean") > 140) & (pl.col("DC.mean") <= 240)).then("moderate")
        .when((pl.col("DC.mean") > 240) & (pl.col("DC.mean") <= 340)).then("high")
        .otherwise("extreme").alias("DC_bin"),
        pl.when(pl.col("ISI.mean") <= 2.2).then("low")
        .when((pl.col("ISI.mean") > 2.2) & (pl.col("ISI.mean") <= 5.0)).then("moderate")
        .when((pl.col("ISI.mean") > 5.0) & (pl.col("ISI.mean") <= 10.0)).then("high")
        .otherwise("extreme").alias("ISI_bin"),
        pl.when(pl.col("BUI.mean") <= 20).then("low")
        .when((pl.col("BUI.mean") > 20) & (pl.col("BUI.mean") <= 36)).then("moderate")
        .when((pl.col("BUI.mean") > 36) & (pl.col("BUI") <= 60)).then("high")
        .otherwise("extreme").alias("BUI_bin"),
        pl.when(pl.col("FWI.mean") <= 3).then("low")
        .when((pl.col("FWI.mean") > 3) & (pl.col("FWI.mean") <= 10)).then("moderate")
        .when((pl.col("FWI.mean") > 10) & (pl.col("FWI.mean") <= 22)).then("high")
        .otherwise("extreme").alias("FWI_bin"),
    )
    # write the classified dataset to a new parquet file
    output_path = output_dir / "full_dataset_with_fwi_bins_only.parquet"
    classified_df.write_parquet(output_path)
    print(f"Full dataset with FWI bins created at {output_path}")
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
    input_file: pl.DataFrame,
    output_dir: Path
    ) -> pl.DataFrame:

    cols = ["DC.mean", "ISI.mean", "BUI.mean", "FWI.mean"]

    def add_quartile_bins(df, col):
        # compute quartile boundaries
        q = (
            df.select(pl.col(col).quantile([0.25, 0.5, 0.75], interpolation="nearest").to_series().to_list())
        )
        q0, q25, q50, q75, q100 = q

        # create bin dumns (1-4)
        df = df.with_columns(
            pl.cut(
                pl.col(col),
                bins=[q0,q25, q50, q75, q100],
                labels=["1", "2", "3", "4"],
                include_lowest=True
            ).alis(f"{col}_bin")
        )

        labels = {
            "1": f"Q1: {round(q0)} - {round(q25)}",
            "2": f"Q2: {round(q25)} - {round(q50)}",
            "3": f"Q3: {round(q50)} - {round(q75)}",
            "4": f"Q4: {round(q75)} - {round(q100)}"        
            }
        # map numeric bins to custom labels
        df = df.with_columns(
            pl.col(f"{col}_bin").map_dict(labels).alias(f"{col}_Quartile_Label")
        )
        return df
    # apply the function to each FWI component
    for col in cols:
        input_file = add_quartile_bins(input_file, col)

    # write the classified dataset to a new parquet file
    output_path = output_dir / "full_dataset_with_fwi_quartile_bins.parquet"
    input_file.write_parquet(output_path)
    print(f"Full dataset with FWI quartile bins created at {output_path}")
    return input_file


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
    input_file: pl.DataFrame,
    output_dir: Path
    ) -> pl.DataFrame:
    # create canopy_bin column
    input_file = input_file.with_columns(
        pl.when(pl.col("lc_class_name").str.contains("open", case=False)).then("open")
        .when(pl.col("lc_class_name").str.contains("forested", case=False)).then("forested")
        .when(pl.col("lc_class_name").str.contains("treed", case=False)).then("treed")
        .otherwise(None).alias("canopy_bin")
    )
    # create class_bin column
    input_file = input_file.with_columns(
        pl.when(pl.col("lc_class_name").str.contains("bog", case=False)).then("bog")
        .when(pl.col("lc_class_name").str.contains("rich fen", case=False)).then("rich fen")
        .when(pl.col("lc_class_name").str.contains("poor fen", case=False)).then("poor fen")
        .when(pl.col("lc_class_name").str.contains("permafrost", case=False)).then("permafrost")
        .when(pl.col("lc_class_name").str.contains("upland", case=False)).then("upland")
        .otherwise(None).alias("class_bin")
    )
    # factor the columns
    input_file = input_file.with_columns(
        pl.col("canopy_bin").cast(pl.Categorical),
        pl.col("class_bin").cast(pl.Categorical)
    )
    # write the classified dataset to a new parquet file
    output_path = output_dir / "full_dataset_with_canopy_and_class_bins.parquet"
    input_file.write_parquet(output_path)
    print(f"Full dataset with canopy and class bins created at {output_path}")
    return input_file


#----- classify season based on Date column ------
"""
if Date<20230616 = Early (May to mid-June) else Late (mid-June to August)
"""
def classify_season_bins(
    input_file: pl.DataFrame,
    output_dir: Path
    ) -> pl.DataFrame:
    input_file = input_file.with_columns(
        pl.when(pl.col("Date") < pl.lit("20230616")).then("Early")
        .otherwise("Late").alias("season_bin")
    )
    # factor the column
    input_file = input_file.with_columns(
        pl.col("season_bin").cast(pl.Categorical)
    )
    # write the classified dataset to a new parquet file
    output_path = output_dir / "full_dataset_with_season_bins.parquet"
    input_file.write_parquet(output_path)
    print(f"Full dataset with season bins created at {output_path}")
    return input_file


# ---------- Main function ----------
def main() -> int:
    # Load configuration and resolve paths
    input_dir = Path(in_rel)
    output_dir = Path(out_rel)
    prog_file = Path(prog_rel)
    fwi_file = Path(fwi_rel)

    print(f"Input directory: {input_dir}\nOutput directory: {output_dir}\nProgression file: {prog_file}\nFWI file: {fwi_file}")

    # Create full dataset
    full_dataset = create_full_dataset(input_dir, output_dir)
    print("Full dataset created. Number of rows:" , full_dataset.shape[0])

    # Get east vs west classification from progression centroids
    east_west_df = get_east_west_from_centroids(prog_file, output_dir)
    print("East/west classification created. Number of rows:" , east_west_df.shape[0])

    # Merge east/west classification with full dataset
    full_dataset = full_dataset.join(east_west_df, on="CLUSTERID", how="left")
    print("Full dataset merged with east/west classification. Number of rows:" , full_dataset.shape[0])

    # Filter out low availability classes
    filtered_dataset = filter_low_availability_classes(full_dataset, output_dir)
    
    # Merge FWI data
    merged_dataset = merge_fwi_data(filtered_dataset, fwi_file, output_dir)

    # Classify FWI into bins
    classified_fwi_bins = classify_fwi_bins(merged_dataset, output_dir)

    # Classify FWI into quartile bins
    classified_fwi_quartiles = classify_fwi_quartile_bins(classified_fwi_bins, output_dir)

    # Classify canopy and class bins
    classified_canopy_class = classify_canopy_and_class_bins(classified_fwi_quartiles, output_dir)

    # Classify season bins
    classified_season = classify_season_bins(classified_canopy_class, output_dir)


if __name__ == "__main__":
       sys.exit(main())
