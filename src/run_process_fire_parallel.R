# Script to run fire processing in parallel
# This script processes multiple fires concurrently using parallel computing
# to extract burned/unburned pixel data for peatland fire selectivity analysis.

# Source the main processing function
source("src/functions/process_fire.R")

# Load required libraries
library(raster)  # For raster operations
library(sf)      # For spatial data handling
library(foreach) # For parallel loops
library(doParallel) # For parallel backend

# Set file paths - UPDATE THESE WITH YOUR ACTUAL FILE PATHS
dnbr_path <- ("")          # Path to folder containing DNBR .tif files (e.g., "/path/to/dnbr/")
peatland_path <- ("")      # Path to merged peatland-canopy raster file (e.g., "/path/to/landcover.tif")
progression_path <- ("")   # Path to fire progression shapefile (e.g., "/path/to/progression.shp")

# Read in raster data
peatland_data <- raster(peatland_path)                    # Load merged peatland-canopy land cover raster

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
num_clusters <- min(20, detectCores() - 1)               # Use up to 20 cores or available-1
cl <- makeCluster(num_clusters)                           # Create cluster
registerDoParallel(cl)                                    # Register parallel backend

# Run parallel processing for each fire
# Note: Update 'prog_poly' in the foreach call if using a subset for testing
results <- foreach(i = 1:nrow(prog_poly), .packages = c("raster", "sf")) %dopar% {
  process_fire(i, prog_poly, dnbr_path, peatland_data)
}

# Stop the parallel cluster
stopCluster(cl)

# Combine all fire results into a single data frame
final_df <- do.call(rbind, results)

# Write results to CSV file
write.csv(final_df, "results/hsf_fire_selectivity_data.csv", row.names = FALSE)