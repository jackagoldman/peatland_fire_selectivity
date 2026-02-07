# Script to run fire processing in parallel
# This script processes multiple fires concurrently using parallel computing
# to extract burned/unburned pixel data for peatland fire selectivity analysis.

# get defaults
config <- config::get(file = "config.yml")

# Source the main processing function
# set working directory to the script's directory
project_root <- config$proj_dir

# get functions
source(file.path(project_root, "src/functions/process_fire.R"))# Script to run fire processing in parallel

# Load required libraries
library(raster)  # For raster operations
library(sf)      # For spatial data handling
library(foreach) # For parallel loops
library(doParallel) # For parallel backend

# Set file paths - UPDATE THESE WITH YOUR ACTUAL FILE PATHS
dnbr_path <- config$dnbr_path      # Path to DNBR raster files
peatland_path <-config$peatland_path # Path to merged peatland-canopy land cover raster
progression_path <- config$progression_path # Path to fire progression shapefile

# Read in fire progression shapefile
prog_poly <- st_read(progression_path)                    # Load fire progression polygons

# Get list of DNBR files from the dnbr_path folder
file_list <- list.files(dnbr_path)                        # List all files in DNBR folder
dnbr_files <- file_list[grep("dnbr", file_list)]          # Filter for files containing "dnbr"

# Extract unique fire IDs from DNBR filenames (remove "_dnbr.tif" suffix)
fire_ids <- unique(gsub("_dnbr", "", gsub(".tif", "", dnbr_files)))

# Filter shapefile to only include fires with available DNBR data
prog_poly <- prog_poly[prog_poly$K_FireID %in% fire_ids, ]

# Set up parallel processing
num_clusters <- min(10, detectCores() - 1)               # Use up to 20 cores or available-1
cl <- makeCluster(num_clusters, outfile = config$worker_log)   # Create cluster
registerDoParallel(cl)                                    # Register parallel backend

# Run parallel processing for each fire
# Note: Update 'prog_poly' in the foreach call if using a subset for testing
results <- foreach(i = 1:nrow(prog_poly), .packages = c("raster", "sf")) %dopar% {
  res <- process_fire(i, prog_poly, dnbr_path, peatland_path)
  arrow::write_parquet(res, paste0("data/used_available_per_fire/used_available_", res$K_UniqueID[1], ".parquet"), compression = "snappy")
  NULL
}

# Stop the parallel cluster
stopCluster(cl)

