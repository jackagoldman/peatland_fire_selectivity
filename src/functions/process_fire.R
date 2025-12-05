# Function to process each fire

#' Process a single fire to extract burned and unburned pixel data
#'
#' This function processes a fire by its index in the progression polygon dataset,
#' extracting land cover information for burned and unburned pixels within the fire boundary.
#' It uses DNBR and pre-merged peatland-canopy rasters to identify and classify pixels.
#' Requires global variable `dnbr_path` set to the folder containing DNBR .tif files.
#'
#' @param i Integer index of the fire in prog_poly
#' @param prog_poly sf data frame of fire progression polygons with columns like K_FireID, CLUSTERID, etc.
#' @param peatland_path Character string path to the pre-merged peatland-canopy raster file
#' @return A data frame with rows for each pixel (burned and unburned), containing:
#'   - Used: 1 for burned, 0 for unburned
#'   - Lc_class: Land cover class (1-17 from merged peatland-canopy)
#'   - X, Y: Longitude and latitude of pixel centroid (EPSG:4326)
#'   - P_size: Pixel area in square units
#'   - Fire_ID, K_UniqueID, CLUSTERID, AREA, C_AREA: Fire metadata
#' @import sf raster 
#' @export

process_fire <- function(i, prog_poly, dnbr_path , peatland_path) {
  require(sf)
  require(raster)
  
  # Create logs directory if it doesn't exist
  dir.create("logs", showWarnings = FALSE)
  
  # Initialize error and fire logs
  error_log_file <- "logs/error_log.txt"
  fire_log_file <- paste0("logs/fire_log_", i, ".txt")
  
  tryCatch({
    poly <- prog_poly[i, ]
  
  # Transform polygon to CRS 4326 for centroids later
  poly_4326 <- st_transform(poly, 4326)
  
  k_fireid <- poly$K_FireID # for matchin dnbr files
  cluster_id <- poly$CLUSTERID # so we only work on one cluster at a time for a given fire
  
  # Read corresponding DNBR raster
  dnbr_raster <- raster(paste0(dnbr_path, k_fireid, "_dnbr.tif"))
  
  # Load peatland raster
  peatland_data <- raster(peatland_path)
  
  #raster crs - retrieve CRS for later use in transforming points
  raster_crs = crs(dnbr_raster)
  
  # Crop peatland raster to DNBR extent - ensures only relevant area is processed
  peatland_raster = crop(peatland_data, dnbr_raster)
  
  # Discard the full peatland raster to free memory
  peatland_data <- NULL
  gc()  # Force garbage collection
  
  #dnbr
  dnbr_raster_res = resample(dnbr_raster, peatland_raster)


  # Find burned pixels (> 0.1)
  burned <- dnbr_raster > 0.1
  
  #mask landc cover raster to burned areas
  masked_peatland_raster = peatland_raster
  masked_peatland_raster[!burned] = NA # contains land cover only in burned areas
  
  # make sure crs for poly is the same as peatland_raster
  poly <- st_transform(poly, crs(peatland_raster))

  # rasterize polygon for spatial masking
  polygon_raster = rasterize(poly, peatland_raster)# conversts the fire progression polygon into a reask mask matching the grid of peatland_raster. Pixels inside the polygon are set to a value of 1 and outside toNA or 0. Create a spatial max for fire boundary
  
  #align extents to prevent cropping errors
  polygon_raster = crop(polygon_raster, st_bbox(masked_peatland_raster))
  masked_peatland_raster = crop(masked_peatland_raster, st_bbox(polygon_raster)) # crops both rasters to each others bounding box to ensure identical extents, avoids errors from mismatched raster dimensions during masking, as mask requires aligned inputs.
  
  #final mask to the polygon boundary
  burned_peatland_raster = mask(masked_peatland_raster, polygon_raster) # applis the polygon mask. setting pixels outside the polygon to NA. The results is a raster with alndcover values only for burned pixels inside the fire polygon.

  # Find unburned pixels (< 0.1)
  unburned <- dnbr_raster < 0.1 #unburned mask
  
  #apply the mask for unburned pixels
  masked_peatland_raster_ub = peatland_raster
  masked_peatland_raster_ub[!unburned] = NA
  
  # cross crop using rasterized polygon from above
  polygon_raster_ub = crop(polygon_raster, st_bbox(masked_peatland_raster_ub))
  masked_peatland_raster_ub = crop(masked_peatland_raster_ub, st_bbox(polygon_raster_ub))
  
  #mask the raster with the polygon
  unburned_peatland_raster = mask(masked_peatland_raster_ub, polygon_raster_ub) # unburned peatland raster with landcover values only for unburned pixels inside the fire polygon.
  
  
  # Get centroids for burned pixels
  if (cellStats(burned_peatland_raster, 'countNA', na.rm = TRUE) > 0) {
    burned_pts <- rasterToPoints(burned_peatland_raster)
    coords_burned_raw <- burned_pts[, 1:2]
    lc_burned <- burned_pts[, 3]
    burned_pts_sf <- st_as_sf(data.frame(x = coords_burned_raw[, 1], y = coords_burned_raw[, 2]), 
                              coords = c("x", "y"), crs = raster_crs)
    burned_pts_4326 <- st_transform(burned_pts_sf, 4326)
    coords_burned <- st_coordinates(burned_pts_4326)
    burned_df <- data.frame(
      Used = 1,
      Lc_class = lc_burned,
      X = coords_burned[, 1],  # Longitude
      Y = coords_burned[, 2],  # Latitude
      Fire_ID = k_fireid,
      K_UniqueID = poly$K_UniqueID,
      CLUSTERID = cluster_id,
      AREA = poly$AREA, 
      C_AREA = poly$C_AREA
    )
  } else {
    burned_df <- data.frame(
      Used = integer(),
      Lc_class = integer(),
      X = numeric(),
      Y = numeric(),
      Fire_ID = character(),
      K_UniqueID = character(),
      CLUSTERID = character(),
      AREA = numeric(),
      C_AREA = numeric()
    )
  }
  
  # Get centroids for unburned pixels
  if (cellStats(unburned_peatland_raster, 'countNA', na.rm = TRUE) > 0) {
    unburned_pts <- rasterToPoints(unburned_peatland_raster)
    coords_unburned_raw <- unburned_pts[, 1:2]
    lc_unburned <- unburned_pts[, 3]
    unburned_pts_sf <- st_as_sf(data.frame(x = coords_unburned_raw[, 1], y = coords_unburned_raw[, 2]), 
                                coords = c("x", "y"), crs = raster_crs)
    unburned_pts_4326 <- st_transform(unburned_pts_sf, 4326)
    coords_unburned <- st_coordinates(unburned_pts_4326)
    unburned_df <- data.frame(
      Used = 0,
      Lc_class = lc_unburned,
      X = coords_unburned[, 1],  # Longitude
      Y = coords_unburned[, 2],  # Latitude
      Fire_ID = k_fireid,
      K_UniqueID = poly$K_UniqueID,
      CLUSTERID = cluster_id, 
      AREA = poly$AREA, 
      C_AREA = poly$C_AREA
    )
  } else {
    unburned_df <- data.frame(
      Used = integer(),
      Lc_class = integer(),
      X = numeric(),
      Y = numeric(),
      Fire_ID = character(),
      K_UniqueID = character(),
      CLUSTERID = character(),
      AREA = numeric(),
      C_AREA = numeric()
    )
   
  }

  # Combine burned and unburned for this fire
  fire_df <- rbind(burned_df, unburned_df)
  
  # get total pixels, burned and unburned
  total_pixels <- nrow(fire_df)
  burned_pixels <- sum(fire_df$Used ==1)
  unburned_pixels <- sum(fire_df$Used ==0)
  
  # add to df 
  fire_df$P_num <- ifelse(fire_df$Used == 1, burned_pixels, unburned_pixels)
  fire_df$P_total <- total_pixels
  
  # Log fire processing results
  fire_log <- data.frame(
    Fire_ID = k_fireid,
    Burned_Pixels = burned_pixels,
    Unburned_Pixels = unburned_pixels,
    Total_Pixels = total_pixels
  )
  write.table(fire_log, file = fire_log_file, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  return(fire_df)
  }, error = function(e) {
    # Log error
    error_msg <- paste("Error processing fire index", i, ":", e$message)
    write(error_msg, file = error_log_file, append = TRUE)
    # Return empty data frame
    return(data.frame(
      Used = integer(),
      Lc_class = integer(),
      X = numeric(),
      Y = numeric(),
      Fire_ID = character(),
      K_UniqueID = character(),
      CLUSTERID = character(),
      AREA = numeric(),
      C_AREA = numeric(),
      P_num = integer(),
      P_total = integer()
    ))
  })
}
