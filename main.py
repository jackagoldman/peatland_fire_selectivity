# MAIN WORKFLOW SCRIPT
# This script orchestrates the entire workflow for processing fire progression data, calculating burned areas, 
# computing selectivity indices, and creating a consolidated dataset for analysis.
#  It calls functions from the src/functions directory to perform each step of the workflow.


# please see file_descriptions.txt for detailed descriptions of each function and workflow step, 
# as well as the config.yml file for configuration details.


# Import necessary libraries
import os
from pathlib import Path
from datetime import datetime




# Bring in functions from the src/functions directory
from src.functions.utils import SETTINGS, cluster_areas, cluster_burned_area, collect_cluster_parquets
from src.functions.process_fire import run_all_fires_parallel
from src.functions.process_polygons import make_prog_sub
from src.functions.calculate_ratios import calculate_ratios
from src.functions.convert_lc_class import add_lc_class_name_to_parquet
from src.functions.selectivity_indices import calculate_selectivity
from src.functions.create_full_dataset import concat_dataframes, get_east_west_from_centroids, filter_low_availability_classes, merge_fwi_data, classify_fwi_bins, classify_fwi_quartile_bins, classify_canopy_and_class_bins, classify_season_bins
from src.functions.summary_statistics import summarize_lands_abundance, find_peatland_fires, find_peatland_dominant_fires, chunk_fires_three_bins


# Define main function to run the workflow

def main():
    print("Run the full workflow")

    now = datetime.now()
    print(f"Workflow started at {now}")

# STEP 1:
    # ----------- process the progressions to get ready for processing
    print("Step 1: processing progression shapefile to subset to polygons with peat and dNBR coverage, and buffer them")
    # progression_path = SETTINGS["python"]["progression_2"]  
    # peat_raster_path = SETTINGS["python"]["peatland_raster"]
    # dnbr_folder = SETTINGS["python"]["dnbr_rasters"]

    # prog_sub = make_prog_sub(
    #     progression_path=progression_path,
    #     peat_raster_path=peat_raster_path,
    #     dnbr_folder=dnbr_folder,
    #     area_quantile=0.20,
    #     buffer_distance_m=200.0,
    #     working_crs="EPSG:3347",
    #     fix_invalid_geoms=True,
    # )

    # # Optional: save the result
    # if not prog_sub.empty:
    #     out_path = SETTINGS["python"]["progression_2_output"]
    #     prog_sub.to_file(out_path, driver="GeoJSON")
    #     print(f"Saved prog_sub → {out_path}")

# STEP 2: 
    # Step 2.1 ------- process the fires in parrallel
    print("Step 2.1: Process fire files in parallel using worker function")
    # PROGRESSION_PATH = SETTINGS["python"]["progression_2_output"]
    # LANDSCAPE_PATH = SETTINGS["python"]["progression"]
    # OUT_DIR = SETTINGS["python"]["new_ua_processed"]  # new folder for Parquet outputs

    # run_all_fires_parallel(
    #         progression_path=PROGRESSION_PATH,
    #         landscape_path=LANDSCAPE_PATH,
    #         out_dir=OUT_DIR,
    #         max_processes=8,
    #         dnbr_threshold=0.1
    #     )
    # Step 2.2 ------- collect the cluster parquets into a single directory for consolidation step
    print("Step 2.2: Collect cluster Parquet files into a single directory for consolidation")
    #
    # collect_cluster_parquets(SETTINGS['python']['new_ua_processed'])

# STEP 3:
    # ----- 3 Get burned area of polygons via cluster areas 
    print("Step 3: Calculate cluster areas from progression and landscape shapefiles")
    # check if outpath exists, if it does, skip this step (since it can be time consuming and we don't want to overwrite)
    if os.path.exists(SETTINGS["python"]["cluster_areas"]):
        print(f"Cluster areas file already exists at {SETTINGS['python']['cluster_areas']}, skipping area calculation.")
    if not os.path.exists(SETTINGS["python"]["cluster_areas"]):
        cluster_areas(progression_path=SETTINGS["python"]["progression_2_output"],
                    landscape_path=SETTINGS["python"]["progression"],
                    out_path=SETTINGS["python"]["cluster_areas"],
                    id_col="CLUSTERID",
                    area_crs="EPSG:3347")
    
# STEP 4:
    # 4.1 conver lc class
    # add_lc_class_name_to_parquet(SETTINGS["python"]["by_cluster_ids"])
    # 4.2 calculate burn and available ratios
    # calculate_ratios(in_dir=SETTINGS["python"]["by_cluster_ids"], out_dir=SETTINGS["python"]["ratios_output"])


# STEP 5: 
    # ------ Calculate selectivity indices
    #calculate_selectivity(in_dir=SETTINGS["python"]["ratios_output"], out_dir=SETTINGS["python"]["selectivity_output"])

# STEP 6:
    # ------ Consolidate all selectivity outputs into a single Parquet for analysis

    # ------ build full datatset using combined_df
    # Confirm keys exist
    in_rel = SETTINGS["python"]["selectivity_output"]       # e.g., "data/ua_fire_combined/"
    out_rel = SETTINGS["python"]["full_dataset_output"]  # e.g., "data/full_dataset/"
    prog_rel = SETTINGS["python"]["progression"]  # e.g., "data/progression/"
    fwi_rel = SETTINGS["python"]["fwi_input"]  # e.g., "data/fwi/fwi.parquet"

    # Load configuration and resolve paths
    input_dir = Path(in_rel)
    output_dir = Path(out_rel)
    prog_file = Path(prog_rel)
    fwi_file = Path(fwi_rel)

    print(f"Input directory: {input_dir}\nOutput directory: {output_dir}\nProgression file: {prog_file}\nFWI file: {fwi_file}")

    # Create full dataset
    full_dataset = concat_dataframes(input_dir, save_output=False)
    print("Full dataset created. Number of rows:" , full_dataset.shape[0])

    # Get east vs west classification from progression centroids
    east_west_df = get_east_west_from_centroids(prog_file, save_output=False)
    print("East/west classification created. Number of rows:" , east_west_df.shape[0])

    # using east vs west centroids, get the ecozone classification for each fire 
    # read in boreal NA ecozones shapefile 
    # make new column called ecozone which is the ecozone classification for each fire (based on the centroid of the fire progression)


    # Merge east/west classification with full dataset
    full_dataset = full_dataset.join(east_west_df, on="CLUSTERID", how="left")
    print("Full dataset merged with east/west classification. Number of rows:" , full_dataset.shape[0])

    # Filter out low availability classes
    filtered_dataset = filter_low_availability_classes(full_dataset, save_output=False)
    print("Low availability classes filtered. Number of rows:" , filtered_dataset.shape[0])

    # Merge FWI data
    merged_dataset = merge_fwi_data(filtered_dataset, fwi_file, save_output=False)
    print("FWI data merged. Number of rows:" , merged_dataset.shape[0])

    # Classify FWI into bins
    classified_fwi_bins = classify_fwi_bins(merged_dataset, save_output=False)
    print("FWI classified into bins. Number of rows:" , classified_fwi_bins.shape[0])

    # Classify FWI into quartile bins
    classified_fwi_quartiles = classify_fwi_quartile_bins(classified_fwi_bins, save_output=False)
    print("FWI classified into quartiles. Number of rows:" , classified_fwi_quartiles.shape[0])

    # Classify canopy and class bins
    classified_canopy_class = classify_canopy_and_class_bins(classified_fwi_quartiles, save_output=False)
    print("Canopy and class bins classified. Number of rows:" , classified_canopy_class.shape[0])

    # Classify season bins
    classified_season = classify_season_bins(classified_canopy_class, output_dir)
    print("Season bins classified and full dataset saved. Number of rows:" , classified_season.shape[0])
    date_str = datetime.today().strftime("%Y%m%d")
    print(f"Full dataset saved to {output_dir / f'fire_selectivity_full_dataset_{date_str}.parquet'}")



##-------------------------------
#  Summary statistics
##-------------------------------

# summary statistics calculated as per Dan's request
# 1. landscape abundance vs burn frequency
# 2. burn rate of peatland fires vs adjacent uplands
# 3. what fire shows a pereference for burning in peatlands? preference vs indifference vs avoidance





if __name__ == "__main__":
    main()
   

