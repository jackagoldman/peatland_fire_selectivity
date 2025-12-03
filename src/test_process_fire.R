library(raster)
library(sf)


# set working directory to the script's directory
project_root <- getwd()

# get functions
source(file.path(project_root, "src/functions/process_fire.R"))


# set paths
dnbr_path = ("G:/Fire_Selectivity/NickPelletier - do not delete/dNBR rasters/")
peatland_path = ("G:/Fire_Selectivity/NickPelletier - do not delete/Peat Map Pontone/PeatlandMap8b_2023_07_17.tif")
progression_path = ("G:/Fire_Selectivity/NickPelletier - do not delete/fire polygons 2023/landscape_processed_polygons_km_oct18.shp")

# read in rasters
peatland_data <- raster(peatland_path)
# read in fire progression shapefile
fire_prog <- st_read(progression_path)

#
print(dim(peatland_data))
print(dim(fire_prog))

# Get files from dnbr_raster folder
file_list <- list.files(dnbr_path)
#filter files with ".shp" in the name
dnbr_files = file_list[grep("dnbr", file_list)]

# extract fire IDs from filenames
fire_ids <- unique(gsub("_dnbr", "", gsub(".tif", "", dnbr_files))) 

# Filter shapefile
prog_poly <- fire_prog[fire_prog$K_FireID %in% fire_ids, ]

# for testing, filter to a single fire
prog_poly = prog_poly %>% filter(NFIREID == "191")

# read in dnbr file
dnbr_test = raster::raster(paste0(dnbr_path, dnbr_files[[1]]))


results = list()
for (i in 1:nrow(prog_poly)){
  cat("Processing fire", i, "\n")
  results[[i]] = process_fire(i, prog_poly, dnbr_path, peatland_data)
}

# View Test Results
cat("\nTest Results Summary:\n")
cat("Number of fires processed:", length(results), "\n")
for (i in 1:length(results)) {
  df <- results[[i]] 
  cat("Fire", i, "- Pixels:", nrow(df), "| Burned:", sum(df$Used == 1), "| Unburned:", sum (df$Used == 0), "\n")
  if (nrow(df) > 0) {
  cat("Sample rows:\n")
  print(head(df, 3))
  cat("\n")
  }
}

# save csv for further inspection
test_output <- do.call(rbind, results)
write.csv(test_output, "test/test_process_fire_results.csv", row.names = FALSE)
cat("Full results saved to test_process_fire_results.csv\n")


