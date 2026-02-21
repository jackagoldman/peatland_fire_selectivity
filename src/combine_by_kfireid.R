# combine all files by k_fireID

# Steps:
#- Get df k_fireID and corresponding kuID with each row being a list of kuids?
#- Write a function that: filters through each row or k_fireID, read all the corresponding used_available files that have the kuid in kfireid. Bind them to a single df and write to parquet.
#- in binded per fire file, make a row that is step_id which corresponds to k_unique_id.
#- find which k_unique_id is which progression step (first to final)
#- 


# Make a df which has a column for k_fireID and each corresponding K_UniqueID
# get progression path
prog_path <- config::get("progression_path")

# read in progression polygons
progs <- sf::st_read(prog_path)

# get the required columns
progs_ids <- dplyr::select(progs, c("K_FireID", "K_UniqueID")


# get path to files
ua_path <– config::get("used_available")

#set folder path
ua_folder <- file.path(ua_path)

library(arrow)
library(dplyr)
library(stringr)

df = prog_ids
input_dir = ua_path 
output_dir = config::get("combined_output")

merge_parquet_by_fire <- function(df, input_dir, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # get all parquet files
  files <- list.files(input_dir, pattern = "\\.parquet$", full.names = TRUE)

  # loop over fire IDs
  for (fire in unique(df$K_fireID)) {

    # all UIDs for this fire
    uids <- df$K_UniqueID[df$K_fireID == fire]

    # find parquet files where the filename contains one of the UIDs
    matched_files <- files[
      sapply(files, function(f) any(str_detect(basename(f), paste0("\\b", uids, "\\b"))))
    ]

    if (length(matched_files) == 0) {
      message("No files found for fire: ", fire)
      next
    }

    # read and merge
    tables <- lapply(matched_files, read_parquet)
    combined <- bind_rows(tables)

    # write output
    write_parquet(combined, file.path(output_dir, paste0("fire_", fire, ".parquet")))
    message("Wrote fire_", fire)
  }
}