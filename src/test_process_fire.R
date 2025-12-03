library(raster)
library(sf)

# set paths
dnbr_path = ("G:/Fire_Selectivity/NickPelletier - do not delete/dNBR rasters/")
peatland_path = ("G:/Fire_Selectivity/NickPelletier - do not delete/Peat Map Pontone/PeatlandMap8b_2023_07_17.tif")
canopy_path = ("E:/Jack/data/peatland_fire_selectivity/scanfi/scanfi_canopy_data_classified_8bit.tif")  
progression_path = ("G:/Fire_Selectivity/NickPelletier - do not delete/polygons 2023/landscape_processed_polygons_km_oct18.shp")
# read in rasters
peatland_data <- raster(peatland_path)
canopy_data_classified <- raster(canopy_path)  
# read in fire progression shapefile
fire_prog <- st_read(progression_path)

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
  results[[i]] = process_fire(i, prog_poly, canopy_data_classified, peatland_data)
}