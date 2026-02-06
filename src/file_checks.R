###### Script to check files to see which are missing ######
library(dplyr)
library(purrr)
#get path to koreen folder which should have all processed fires
km_path <- "W:/jack/NickPelletier - do not delete/selectivity stats km"

# my processed files
jg_path <- "W:/jack/per_fire2"


# get file list
files_km <- list.files(path = km_path, pattern = ".csv$")
files_jg <- list.files(path = jg_path, pattern = ".csv$")


# get pattern to match for the numbers 
km_list <- gsub("[^0-9]", "", files_km)

# get pattern to match for the numbers for jg
jg_list <- gsub("[^0-9]", "", files_jg)

# get the fires that were not processed (next time error log should be more robust)
not_processed <- km_list[!(km_list %in% jg_list)]
# - weird that there are over 2000, it might be that the KUnique is different from clsuterid

# check the progression shp
progression_path <- ("G:/Fire_Selectivity/NickPelletier - do not delete/fire polygons 2023/landscape_processed_polygons_km_oct18.shp")   # Path to fire progression shapefile (e.g., "/path/to/progression.shp")
prog_poly <- sf::st_read(progression_path)

# check str
str(prog_poly)

#loop through rows clusterID, K_FireID and K_UniqueID to see which row matches all in km_list
id_cols <- prog_poly %>% sf::st_drop_geometry() %>% select(c(CLUSTERID, K_FireID, K_UniqueID))

matches <- id_cols %>% 
  summarise(across(everything(), ~sum(.x %in% km_list, na.rm = TRUE)))

print(matches)
# i used K_FireID, she used K_UniqueID. Seems like I will have to redo.

# read in 1 file from files_jg to check columns
first_file_jg <- files_jg[[1]]

# build file path 
path1 <- file.path(paste0(km_path, "/selectivity_stats_331.csv"))

# read the csv
first_file <- read.csv(path1)

#check column names
colnames(first_file)


