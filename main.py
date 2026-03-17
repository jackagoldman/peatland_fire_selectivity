import os
import re
import warnings
from typing import Optional, Union, Iterable
from pathlib import Path


from src.functions.consolidate_fires import consolidate_all_fires, consolidate_one_fire, concatenate_fires
from src.functions.utils import SETTINGS, cluster_areas, cluster_burned_area, merge_fire_dir_inplace
from src.functions.process_fire import run_all_fires_parallel

def main():
    print("Run the full workflow")
# STEP 2: 

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


# STEP 3:
    # ----- 3.1 Get burned area of polygons via cluster areas 
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
    # ------- Execute consolidation of fire files after worker processing is complete.
    print("Step 4: Consolidate fire files from worker outputs into final Parquet files")

    # get file paths from config
    OUT_DIR = SETTINGS["python"]["new_ua_processed"]
    FINAL_DIR = SETTINGS["python"]["new_ua_consolidated"]

    # if FINAL_DIR has 644 files in it, we can assume consolidation is done and skip this step (since it can be time consuming)
    if os.path.exists(FINAL_DIR) and len(os.listdir(FINAL_DIR)) == 644:
        print(f"Final consolidated files already exist in {FINAL_DIR}, skipping consolidation step.")
    if len(os.listdir(FINAL_DIR)) < 644:
        consolidate_all_fires(
            out_dir=OUT_DIR,
            dest_dir=FINAL_DIR,
            # Optional filters:
            # fires=["12345", "67890"],
            # n=10, random_sample=True,
            # Per-chunk dedup (subset of columns):
            # drop_duplicates_on=["x", "y", "lc_class", "used"],
        )


# STEP 5: 
    print("Step 5: Calculate burned area of pixels from consolidated fire files using area factor column")

    if os.path.exists(SETTINGS["python"]["burned_areas"]):
        print(f"Burned area file already exists at {SETTINGS['python']['burned_areas']}, skipping burned area calculation.")
    # -----  Calculate area of burned pixels
    cluster_burned_area(
        fire_parquet_path=SETTINGS["python"]["new_ua_consolidated"],
        out_path=SETTINGS['python']['burned_areas'],
        id_col="CLUSTERID",
        factor_col="area_factor_m2_per_pixel",
    ) 


# Step 6: 
  # join area to individual parquet files
  
    merge_fire_dir_inplace(
                fire_dir=SETTINGS["python"]["new_ua_consolidated"],
                progression_landscape_path=SETTINGS["python"]["cluster_areas"],
                burned_area_path=SETTINGS['python']['burned_areas'],
                recursive=False,               # set True to recurse subdirectories
                cluster_col="CLUSTERID",
                poly_area_col="poly_area",
                buffer_area_col="buffer_area",
                burned_area_col="burned_area_m2",
            )




if __name__ == "__main__":
    main()
   

