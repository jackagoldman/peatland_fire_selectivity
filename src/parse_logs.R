# Parse fires from error logs to identify which fires need to be reprocessed
# This script reads error logs from the parallel processing run to extract fire IDs that encountered errors during processing. It then saves these fire IDs to a CSV file for further analysis or reprocessing.

# Load required libraries
library(dplyr)  # For data manipulation
library(stringr) # For string operations


# Set path to the error log file generated during parallel processing
config <- config::get(file = "config.yml")

# Read the error log file
error_log_path <- config$error_log
error_log <- readLines(error_log_path)

# extract all values before : in the error log lines
ku_ids <- str_extract(error_log, "^[^:]+") |> 
  na.omit() |> 
  unique()

# remove "Error processing fire index " from the fire IDs
ku_ids <- str_replace(ku_ids, "Error processing fire index ", "")

# get progression shapefile to get fire IDs
prog_poly <- sf::st_read(config$progression_path)

# filter K_UniqueID in the progression shapefile to only include those in the ku_ids list
prog_poly <- prog_poly %>%
  filter(K_UniqueID %in% ku_ids)

# save as text file
writeLines(ku_ids, "data/fires_to_reprocess.txt")
