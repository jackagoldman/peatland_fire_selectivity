
#steps: 
# merge all parquet files into one big dataset
# convert lc class column to class names:
# 1 = Open Bog
# 2 = Treed Bog
# 3 = Forested Bog
# 4 = Open Rich Fen
# 5 = Treed Rich Fen
# 6 = Forested Rich Fen
# 7 = Open Poor Fen
# 8 = Treed Poor Fen
# 9 = Forested Poor Fen
# 10 = Open PPC 
# 11 = Treed PPC
# 12 = Forested PPC
# 13 = Mineral Wetlands
# 14 = Water
# 15  = Upland
# 16 = Agriculture
# 17 = Urban
  
# get defaults
config <- config::get(file = "config.yml", config = "defaults")

# get the output of combined by fire used available
combined_path <- config$combined_output

# for all files in the combined output folder, read them in and merge into one big dataset
library(arrow)
library(dplyr)
library(stringr)

output_dir = config$full_dataset_output
print("Starting merge process...")

# Use Arrow's open_dataset for out-of-memory processing
ds <- open_dataset(combined_path, format = "parquet")



# convert lc class column to class names
lc_class_names <- c(
  "1" = "Open Bog",
  "2" = "Treed Bog",
  "3" = "Forested Bog",
  "4" = "Open Rich Fen",
  "5" = "Treed Rich Fen",
  "6" = "Forested Rich Fen",
  "7" = "Open Poor Fen",
  "8" = "Treed Poor Fen",
  "9" = "Forested Poor Fen",
  "10" = "Open PPC",
  "11" = "Treed PPC",
  "12" = "Forested PPC",
  "13" = "Mineral Wetlands",
  "14" = "Water",
  "15" = "Upland",
  "16" = "Agriculture",
  "17" = "Urban"
)

# Create lookup table
lc_lookup <- tibble(
  Lc_class = as.integer(names(lc_class_names)),
  lc_class_name = unname(lc_class_names)
)

# Use Arrow join (no collect)
combined <- ds %>%
  mutate(Lc_class = as.character(Lc_class)) %>%
  left_join(lc_lookup, by = "Lc_class")

# write output
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write_parquet(combined, file.path(output_dir, "full_dataset.parquet"), compression = "snappy")

# write out a files of just the X and Y columns 
write_parquet(select(combined, X, Y), file.path(output_dir, "full_dataset_coords.parquet"), compression = "snappy")