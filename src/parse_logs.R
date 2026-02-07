# Parse fires from error logs to identify which fires need to be reprocessed
# This script reads error logs from the parallel processing run to extract fire IDs that encountered errors during processing. It then saves these fire IDs to a CSV file for further analysis or reprocessing.

# Load required libraries
library(dplyr)  # For data manipulation
library(stringr) # For string operations


# Set path to the error log file generated during parallel processing
config <- config::get(file = "config.yml")

# get progression path
progression_path <- config::get("progression_path", config = "local")

# Read the error log file
error_log_path <- config$error_log
error_log <- readLines(error_log_path)

# extract all values before : in the error log lines
ku_ids <- str_extract(error_log, "^[^:]+") |> 
  na.omit() |> 
  unique()

# remove "Error processing fire index " from the fire IDs
ku_ids <- str_replace(ku_ids, "Error processing fire index ", "")

#cealn ku_ids and convert to integer
ku_ids <- stringr::str_trim(ku_ids)
ku_ids <- as.numeric(ku_ids)     

# get progression shapefile to get fire IDs
prog_poly <- sf::st_read(progression_path)


# lengths differ 54 ku_ids to 36 in prog_poly, likely some ku_ids not in the shapefile?
prog_poly_ku_ids <- prog_poly$K_UniqueID

# length of prog_poly_ku_ids
length(prog_poly_ku_ids) 

# duplicates in prog_poly_ku_ids
sum(duplicated(prog_poly_ku_ids)) # no duplicates

# check which ku_ids are not in the progression shapefile
ku_ids_not_in_shapefile <- ku_ids[!ku_ids %in% prog_poly_ku_ids]
 length(ku_ids_not_in_shapefile) # 18 



# seems like there are 18 missing, but these are not in the selectivity stats which is odd.
# How could these of been processed or missing? check the logs and get the K_UniqueID from processing logs and make a list to check against
path_to_processing_log <- config$processing_logs

# loop through each processing log file and extract K_UniqueIDs
processing_log_files <- list.files(path_to_processing_log, full.names = TRUE)
processed_ku_ids <- numeric(0)
for (log_file in processing_log_files) {
  log_lines <- readr::read_tsv(log_file, quote = "\"", trim_ws = TRUE)
  log_lines$K_UniqueID <- as.numeric(log_lines$K_UniqueID) # convert to numeric
  ku_ids_in_log <- log_lines$K_UniqueID
 # append
  processed_ku_ids <- c(processed_ku_ids, ku_ids_in_log)
}
processed_ku_ids

# check which ku_ids from the error log are in the processed ku_ids list
ku_ids_in_processed <- ku_ids[ku_ids %in% processed_ku_ids]

# how many duplicates in the processed ku_ids?
num_duplicates <- sum(duplicated(processed_ku_ids)) # no duplicates

# check processed in the progression shapefile
prog_poly_processed <- prog_poly |> 
  dplyr::filter(K_UniqueID %in% processed_ku_ids)

# check which progression poly ids arent in the processed ku ids
not_processed_ku_ids_in_shapefile <- prog_poly_ku_ids[!prog_poly_ku_ids%in% processed_ku_ids ]


# are these in the error log?
not_processed_ku_ids_in_shapefile_in_error_log <- not_processed_ku_ids_in_shapefile[not_processed_ku_ids_in_shapefile %in% ku_ids]
# error log is wrong!!?

#filter prog poly by to only include not processed
not_processed_poly <- prog_poly |> 
  dplyr::filter(K_UniqueID %in% not_processed_ku_ids_in_shapefile)

nrow(not_processed_poly)

# save these as fires to reprocess shapefile
sf::st_write(not_processed_poly, "data/fire_progression_polygons/fires_to_reprocess/progression_polygons_to_reprocess.shp")
