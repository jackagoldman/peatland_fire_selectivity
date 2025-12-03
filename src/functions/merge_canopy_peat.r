# merged canopy and peatland function
library(raster)

#' Merge Peatland and Canopy Rasters into Unified Land Cover Classification
#'
#' This function combines peatland type and canopy closure rasters into a single
#' land cover classification raster. It reclassifies based on combinations of
#' peatland types (1-9) and canopy classes (1-3), producing 17 land cover classes.
#'
#' @param PeatlandMap RasterLayer of peatland types (1-9)
#' @param CanopyClosure RasterLayer of canopy closure classes (1-3)
#' @return RasterLayer with unified land cover classes (1-17):
#'   - 1-3: Bog types (Open, Treed, Forested)
#'   - 4-6: Rich Fen types
#'   - 7-9: Poor Fen types
#'   - 10-12: PPC (Peat Plateau Complex) types
#'   - 13: Mineral Wetlands
#'   - 14: Water
#'   - 15: Upland
#'   - 16: Agriculture
#'   - 17: Urban
#' @import raster

# Function to merge canopy and peat maps
MergeCanopyPeat <- function(PeatlandMap, CanopyClosure) {
  require(raster)
  r1 <- PeatlandMap
  r2 <- CanopyClosure
  
  r3 <- raster(r1)
  
  # Fill based on conditions
  r3[r1 == 1 & r2 == 1] <- 1  # Open Bog
  r3[r1 == 1 & r2 == 2] <- 2  # Treed Bog
  r3[r1 == 1 & r2 == 3] <- 3  # Forested Bog
  
  r3[r1 == 2 & r2 == 1] <- 4  # Open Rich Fen
  r3[r1 == 2 & r2 == 2] <- 5  # Treed Rich Fen
  r3[r1 == 2 & r2 == 3] <- 6  # Forested Rich Fen
  
  r3[r1 == 3 & r2 == 1] <- 7  # Open Poor Fen
  r3[r1 == 3 & r2 == 2] <- 8  # Treed Poor Fen
  r3[r1 == 3 & r2 == 3] <- 9  # Forested Poor Fen
  
  r3[r1 == 4 & r2 == 1] <- 10 # Open PPC 
  r3[r1 == 4 & r2 == 2] <- 11 # Treed PPC
  r3[r1 == 4 & r2 == 3] <- 12 # Forested PPC
  
  r3[r1 == 5] <- 13 # Mineral Wetlands
  r3[r1 == 6] <- 14 # Water
  r3[r1 == 7] <- 15 # Upland
  r3[r1 == 8] <- 16 # Agriculture
  r3[r1 == 9] <- 17 # Urban
  
  return(r3)
}


