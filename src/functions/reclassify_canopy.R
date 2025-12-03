#' Reclassify Canopy Data into Three Classes
#'
#' This function reclassifies canopy closure percentage data into three categories:
#' 1 (Open: 0-9%), 2 (Treed: 10-24%), 3 (Forested: 25%+).
#'
#' @param canopy_data RasterLayer of canopy closure percentages
#' @return RasterLayer with reclassified canopy classes (1-3)
#' @import raster

reclassify_canopy <- function(canopy_data) {
  reclass_df <- c(0, 9, 1,    # Open
                  10, 24, 2,  # Treed
                  25, Inf, 3) # Forested

  reclass_m <- matrix(reclass_df,
                      ncol = 3,
                      byrow = TRUE)

  canopy_data_classified <- reclassify(canopy_data, reclass_m)

  return(canopy_data_classified)
}

