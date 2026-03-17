import os
import re
import warnings
from typing import Optional, Union, Iterable
from pathlib import Path

# Bring in functions from the src/functions directory
from src.functions.utils import SETTINGS, cluster_areas, cluster_burned_area, merge_fire_dir_inplace, collect_cluster_parquets
from src.functions.process_fire import run_all_fires_parallel
from src.functions.process_polygons import make_prog_sub
from src.functions.calculate_ratios import calculate_ratios
from src.functions.convert_lc_class import add_lc_class_name_to_parquet

def main():
    print("Run the full workflow")

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
    collect_cluster_parquets(SETTINGS['python']['new_ua_processed'])

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
    add_lc_class_name_to_parquet(SETTINGS["python"]["by_cluster_ids"])
    # 4.2 calculate burn and available ratios
    calculate_ratios(in_dir=SETTINGS["python"]["by_cluster_ids"], out_dir=SETTINGS["python"]["ratios_output"])


# STEP 5: 
   # ------ Calculate selectivity indices


if __name__ == "__main__":
    main()
   

