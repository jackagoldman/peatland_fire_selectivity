#####################################################################
# FIRE SELECTIVITY PLOTS
# Jack A. Goldman
#####################################################################

# required libraries
library(tidyverse)
library(ggplot2)
library(sf)
library(googledrive)

googledrive::drive_auth(scopes = c("https://www.googleapis.com/auth/drive", "https://www.googleapis.com/auth/drive.readonly"))
shared_drives <- googledrive::drive_find(q = "name = 'peatlands/fire selectivity'", type = "folder")
print(shared_drives)# set file path for working directory
# List folders in the shared folder
folders <- googledrive::drive_ls(shared_drives$id, type = "folder")
print(folders)
# Create summary_stats directory if it doesn't exist
# Get the ID of the summary_stats folder for later use
summary_stats_id <- folders$id[folders$name == "summary_stats"]


pathname = "/Users/jgoldman/Library/CloudStorage/GoogleDrive-jandrewgoldman@gmail.com/My Drive/3_POSTDOC/projects/selectivity/"

# read in selectivity data
selectivity = read_csv(paste0(pathname, 'Dataset_selectivity_Nov11_KM_notLog_ALL.csv'))


# get landscape polygons, and join to na boreal ecozones
# read in landscape polygons
landscape_polygons <- sf::st_read("data/landscape_processed_polygons_km_oct18.shp") %>%
  st_transform(crs = 4326) # ensure it's in WGS84 lat-long

# read in ecozone polygons
ecozone_polygons <- sf::st_read("/Users/jgoldman/Library/CloudStorage/GoogleDrive-jandrewgoldman@gmail.com/My Drive/3_POSTDOC/projects/selectivity/Ecozones/ecozones.shp") |>
  st_transform(crs = 4326) # ensure it's in WGS84 lat-long

# perform spatial join to assign ecozone to each landscape polygon
landscape_with_ecozones <- st_join(landscape_polygons, dplyr::select(ecozone_polygons, ZONE_NAME), left = TRUE)

# save the joined dataset for later use
st_write(landscape_with_ecozones,  "data/landscape_polygons_with_ecozones.shp")

#  get zone name and K_UniqueID and ZONE_NAME for each polygon
landscape_ecozone_info <- landscape_with_ecozones %>%
  st_set_geometry(NULL) %>% # drop geometry for easier handling
  select(K_UniqueID, ZONE_NAME)

# join landscape ecozone info to selectivity dataset by K_UniqueID
selectivity_with_ecozones <- selectivity %>%
  left_join(landscape_ecozone_info, by = "K_UniqueID")

# save the joined dataset for later use
write_csv(selectivity_with_ecozones,  "data/Dataset_selectivity_km_nov11_all_notlog_with_ecozones.csv")

# read in sel_eco
sel_eco <- read_csv("data/Dataset_selectivity_km_nov11_all_notlog_with_ecozones.csv")

####### MEAN SELECTIVITY BY CLASS WITH ARROWS ########
my_colors <- c("yellow","orange", "red", "darkred")

selectivity.mean<- selectivity %>% group_by(variable) %>% 
  summarise(mean=mean(jacobs,na.rm=T),
            median=median(jacobs,na.rm=T),
            sd=sd(jacobs,na.rm=T),
            p10=quantile(jacobs,probs=0.10,na.rm=T),
            p90=quantile(jacobs,probs=0.90,na.rm=T))

#write
temp_file <- tempfile(fileext = ".csv")
write_csv(selectivity.mean, temp_file)
googledrive::drive_upload(media = temp_file, path = summary_stats_id, name = "jacobs_summary_overall.csv")

# percent area of peatland (variable) type per ecozone (ZONE_NAME) 

area_ecozone <- selectivity_with_ecozones |>  group_by(ZONE_NAME, variable) |> 
  summarise(mean_area = mean(av.area, na.rm = TRUE),
            total_area = sum(av.area, na.rm = TRUE)) |> 
  ungroup() |> 
  group_by(ZONE_NAME) |> 
  mutate(percent_area = mean_area / sum(mean_area, na.rm = TRUE) * 100) |> 
    mutate(across(c(mean_area, total_area, percent_area), ~round(.x, 2)))

# write table of percent area by variable and ecozone and save to summar_stats_id
temp_file <- tempfile(fileext = ".csv")
write_csv(area_ecozone, temp_file)
googledrive::drive_upload(media = temp_file, path = summary_stats_id, name = "percent_area_by_variable_and_ecozone.csv")

# for each ecozone, make a plot and facet wrap them by percent area of each peatland type in that ecozone (variable) and color by variable
p_ecozone <- ggplot(area_ecozone, aes(x = variable, y = percent_area, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  labs(x = "", y = "Percent Area", fill = "Peatland Type") +
  scale_fill_brewer(palette = "Set3") +  # or scale_fill_viridis_d()
  theme_bw() +
  facet_wrap(~ZONE_NAME, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# write table of percent area by variable and ecozone and save to summar_stats_id
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = p_ecozone, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = summary_stats_id, name = "test_percent_area_by_variable_and_ecozone.png")


# mean bui and mean isi by ecozone
fwi_ecozone <- selectivity_with_ecozones |> 
  group_by(ZONE_NAME) |> 
  summarise(
    median_bui = median(BUI.mean, na.rm = TRUE),
    median_isi = median(ISI.mean, na.rm = TRUE),
    median_fwi = median(FWI.mean, na.rm = TRUE),
    median_dc = median(DC.mean, na.rm = TRUE),
    median_dmc = median(DMC.mean, na.rm = TRUE)
  ) |> 
  ungroup() |> 
  mutate(across(where(is.numeric), ~round(.x, 2)))

# write 
temp_file <- tempfile(fileext = ".csv")
write_csv(fwi_ecozone, temp_file)
googledrive::drive_upload(media = temp_file, path = summary_stats_id, name = "median_fwi_bui_isi_by_ecozone.csv")

# convert date column and get median day of burning per ecozone
selectivity_with_ecozones <- selectivity_with_ecozones |> 
  mutate(date = as.Date(Date, format = "%Y-%m-%d")) |> 
  group_by(ZONE_NAME) |> 
  summarise(median_burn_day = median(as.numeric(format(date, "%j")), na.rm = TRUE)) |> 
  ungroup() |> 
  mutate(median_burn_day = round(median_burn_day, 0))
# write
temp_file <- tempfile(fileext = ".csv")
write_csv(selectivity_with_ecozones, temp_file)
googledrive::drive_upload(media = temp_file, path = summary_stats_id, name = "median_burn_day_by_ecozone.csv")

# make a new dir called jacobs plots in google drive folders
# Create a new directory called "jacobs plots" in the Google Drive shared folder
jacobs_plots_id <- googledrive::drive_mkdir("jacobs_plots", path = shared_drives$id)


##### Plot with arrows indicating increasing/decreasing selectivity from "upland" category
sel_plot <- selectivity.mean %>%
  filter(variable != "total_peat" & variable != "mineral") %>%
  arrange(mean) %>%
  mutate(x = row_number())

# find upland index (case-insensitive)
upland_idx <- which(tolower(sel_plot$variable) == "upland")
if (length(upland_idx) == 0) stop("No 'upland' category found in variable column. Adjust the name or inspect sel_plot$variable.")

# compute vertical position just above the upper CI for the upland category
y_range <- max(sel_plot$p90, na.rm = TRUE) - min(sel_plot$p10, na.rm = TRUE)
y_pos_above <- sel_plot$p90[upland_idx] + y_range * 0.05  # adjust 0.05 to move arrows higher/lower

# arrow geometry (tweak lengths as desired)
left_start  <- upland_idx - 0.25
left_end    <- upland_idx - 2
right_start <- upland_idx + 0.25
right_end   <- upland_idx + 2
# clamp to plot range
n <- nrow(sel_plot)
left_end  <- pmax(1, left_end)
right_end <- pmin(n, right_end)

# base plot using numeric x but label axis with category names
p <- ggplot(sel_plot, aes(x = x, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = p10, ymax = p90), width = 0.2) +
  scale_x_continuous(breaks = sel_plot$x, labels = sel_plot$variable) +
  xlab("") +
  ylab("jacobs's index") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

p

# add two arrows with heads pointing away from the upland category,
# positioned above the CI (use inherit.aes = FALSE and pass constants directly)
left_mid  <- (left_start + left_end) / 2
right_mid <- (right_start + right_end) / 2

# vertical offset "one space" above arrows: use a small fraction of the plot y-range
label_offset <- y_range * 0.02
label_y <- y_pos_above + label_offset

p2 <- p +
  geom_segment(x = left_start, xend = left_end, y = y_pos_above, yend = y_pos_above,
               inherit.aes = FALSE,
               arrow = arrow(length = grid::unit(0.3, "cm")),
               lineend = "round",
               color = "black") +
  geom_segment(x = right_start, xend = right_end, y = y_pos_above, yend = y_pos_above,
               inherit.aes = FALSE,
               arrow = arrow(length = grid::unit(0.3, "cm")),
               lineend = "round",
               color = "black") +
  # place labels centered on the arrow lines, one "space" above the arrow
  geom_label(x = right_mid, y = label_y, label = "Increasing",
             inherit.aes = FALSE,
             fontface = "bold", hjust = 0.5, vjust = 0,
             fill = "white", label.size = 0, size = 3) +
  geom_label(x = left_mid, y = label_y, label = "Decreasing",
             inherit.aes = FALSE,
             fontface = "bold", hjust = 0.5, vjust = 0,
             fill = "white", label.size = 0, size = 3)

#save to google drive jacobs
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = p, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "mean_selectivity_by_class.png")

## same plot but with "total_" categories removed
sel_plot2 <- selectivity.mean %>%
  filter(!str_detect(variable, "^total_")) %>%
  arrange(mean) %>%
  mutate(x = row_number())      

# remove underscore from variable names for better x-axis labels
sel_plot2 <- sel_plot2 |> 
  mutate(variable = str_replace_all(variable, "_", " "))
# find upland index (case-insensitive)
upland_idx2 <- which(tolower(sel_plot2$variable) == "upland")
if (length(upland_idx2) == 0) stop("No 'upland' category found in variable column. Adjust the name or inspect sel_plot2$variable.")         
# compute vertical position just above the upper CI for the upland category
y_range2 <- max(sel_plot2$p90, na.rm = TRUE) - min(sel_plot2$p10, na.rm = TRUE)
y_pos_above2 <- sel_plot2$p90[upland_idx2] + y_range2 * 0.05  # adjust 0.05 to move arrows higher/lower
# arrow geometry (tweak lengths as desired)
left_start2  <- upland_idx2 - 0.25
left_end2    <- upland_idx2 - 2
right_start2 <- upland_idx2 + 0.25
right_end2   <- upland_idx2 + 2
# clamp to plot range
n2 <- nrow(sel_plot2)
left_end2  <- pmax(1, left_end2)
right_end2 <- pmin(n2, right_end2)
# base plot using numeric x but label axis with category names
p2 <- ggplot(sel_plot2, aes(x = x, y = mean)) +
    geom_point() +      
    geom_errorbar(aes(ymin = p10, ymax = p90), width = 0.2) +
    scale_x_continuous(breaks = sel_plot2$x, labels = sel_plot2$variable) +
    xlab("") +
    ylab("jacobs's index") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# add two arrows with heads pointing away from the upland category,
# positioned above the CI (use inherit.aes = FALSE and pass constants directly)
left_mid2  <- (left_start2 + left_end2) / 2
right_mid2 <- (right_start2 + right_end2) / 2
# vertical offset "one space" above arrows: use a small fraction of the plot y-range
label_offset2 <- y_range2 * 0.02
label_y2 <- y_pos_above2 + label_offset2
p2 = p2 +
    geom_segment(x = left_start2, xend = left_end2, y = y_pos_above2, yend = y_pos_above2,
                 inherit.aes = FALSE,
                 arrow = arrow(length = grid::unit(0.3, "cm")),
                 lineend = "round",
                 color = "black") +     
    geom_segment(x = right_start2, xend = right_end2, y = y_pos_above2, yend = y_pos_above2,
                 inherit.aes = FALSE,
                 arrow = arrow(length = grid::unit(0.3, "cm")),
                 lineend = "round",
                 color = "black") +     
    # place labels centered on the arrow lines, one "space" above the arrow
    geom_label(x = right_mid2, y = label_y2, label = "Increasing",
               inherit.aes = FALSE,
               fontface = "bold", hjust = 0.5, vjust = 0,
               fill = "white", label.size = 0, size = 3) +     
    geom_label(x = left_mid2, y = label_y2, label = "Decreasing",
               inherit.aes = FALSE,
               fontface = "bold", hjust = 0.5, vjust = 0,
               fill = "white", label.size = 0, size = 3)        

## ggsave
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = p2, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "mean_selectivity_by_class.png")


########### BARPLOTS OF MEAN SELECTIVITY BY FWI AND BUI QUARTILES ###########
## Expand Y-axis to make sample sizes more visible
my_colors <- c("yellow","orange", "red", "darkred")

# FWI
selectivity.mean<- selectivity %>% filter(!is.na(class.bin)) %>% group_by(Quartile_Bins.FWI,class.bin,canopy.bin) %>% 
  summarise(mean=mean(jacobs,na.rm=T),
            median=median(jacobs,na.rm=T),
            sd=sd(jacobs,na.rm=T),
            p05=quantile(jacobs,probs=0.05,na.rm=T),
            p95=quantile(jacobs,probs=0.95,na.rm=T),
            count=n(), 
            .groups = "drop") %>% # drop returns an ungrouped tibble, 
  # ensure "total" appears as the last (right-most) column in the facet
  mutate(canopy.bin = forcats::fct_relevel(as.factor(canopy.bin), "total", after = Inf))
 



fwi.mean.barplot = ggplot(data = selectivity.mean %>% filter(!is.na(Quartile_Bins.FWI)),
       aes(x = Quartile_Bins.FWI, y = mean, fill = Quartile_Bins.FWI)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  geom_text(aes(label = count), vjust = -0.5, color = "black", size = 3) +  # Adding text labels
  labs(x = "", y = "jacobs's index", fill = "FWI quartiles") +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  facet_grid(class.bin ~ canopy.bin, scales = "fixed") +
  theme(axis.text.x = element_text(angle = 90)) +
  # increase vertical spacing between y-axis ticks: fewer, prettier breaks and extra top padding
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     expand = expansion(mult = c(0, 0.15)))


# save
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = fwi.mean.barplot, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "mean_selectivity_by_fwi.png")


##BUI

selectivity.mean<- selectivity %>% filter(!is.na(class.bin)) %>% group_by(Quartile_Bins.BUI,class.bin,canopy.bin) %>% 
  summarise(mean=mean(jacobs,na.rm=T),
            median=median(jacobs,na.rm=T),
            sd=sd(jacobs,na.rm=T),
            p05=quantile(jacobs,probs=0.05,na.rm=T),
            p95=quantile(jacobs,probs=0.95,na.rm=T),
            count=n(), 
            .groups = "drop") %>% # drop returns an ungrouped tibble, 
  # ensure "total" appears as the last (right-most) column in the facet
  mutate(canopy.bin = forcats::fct_relevel(as.factor(canopy.bin), "total", after = Inf))


bui.median.barplot = ggplot(data=selectivity.mean %>% filter(!is.na(Quartile_Bins.BUI)),
       aes(x=Quartile_Bins.BUI,y=median,fill=Quartile_Bins.BUI))+
  geom_bar(stat='identity',alpha=0.9)+
    geom_text(aes(label = count), vjust = -0.5, color = "black", size = 3) +  # Adding text labels
  labs(x="",
       y="jacobs's index",
       fill="BUI quartiles")+
    scale_fill_manual(values=my_colors)+
  theme_bw()+
  facet_grid(class.bin~canopy.bin)+
  theme(axis.text.x = element_text(angle = 90))+
  # increase vertical spacing between y-axis ticks: fewer, prettier breaks and extra top padding
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     expand = expansion(mult = c(0, 0.15)))

bui.median.barplot

temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = bui.median.barplot, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "mean_selectivity_by_bui.png")


##### EXPLORING MEDIAN SELECTIVITY == 0 ISSUE #####
# why are upland and permafrost median selectivity 0?
# look at the median selectivity for only permafrost sites
selectivity.permafrost <- selectivity.mean %>%
  filter(class.bin == "permafrost") 

# median permafrost is all 0 
  head(selectivity.permafrost)

# look at selectivity values for upland sites
selectivity.upland <- selectivity.mean %>%
  filter(class.bin == "upland")

# median upland is all 0
head(selectivity.upland) # whats with the quartile bin 5? , 145 that are binned in NA

### mean jacobs for fwi quatiles structure 2, free_x

selectivity.mean <- selectivity %>%
  group_by(Quartile_Bins.FWI,class.bin,canopy.bin) %>%
  summarise(
    mean = mean(jacobs, na.rm = TRUE),
    median = median(jacobs, na.rm = TRUE),
    sd = sd(jacobs, na.rm = TRUE),
    p05 = quantile(jacobs, probs = 0.05, na.rm = TRUE),
    p95 = quantile(jacobs, probs = 0.95, na.rm = TRUE),
    count = n()
  )


# Create the barplot with counts
p <- ggplot(data = selectivity.mean %>% filter(!is.na(Quartile_Bins.FWI))  %>% filter(!is.na(class.bin)),
            aes(x = as.factor(Quartile_Bins.FWI), y = mean,fill=Quartile_Bins.FWI)) +
  geom_bar(stat = "identity",alpha=0.9) +
  geom_text(aes(label = count), vjust = -0.5, color = "black", size = 3) +
  labs(x="",
       y="jacobs's index",
       fill="FWI quartiles")+ 
  theme_bw() +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(class.bin~canopy.bin, scales = "free_x",ncol=4)+
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + # remove x-axis text
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     expand = expansion(mult = c(0, 0.15)))

p

# save
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = p, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "mean_selectivity_by_fwi_free_x_with_counts.png")



####### East vs West boxplots

##### Original sample size location
# Calculate counts for each combination of 'class.bin', 'Quartile_Bins.BUI', and 'season.bin'
count_data <- selectivity %>% filter(!is.na(class.bin)) %>%
  group_by(class.bin, Quartile_Bins.BUI, region) %>%
  summarise(count = n(),
            y_pos=quantile(jacobs,0.95,na.rm=T))

# Plot the boxplot with counts
y_min <- min(selectivity$jacobs, na.rm = TRUE)
y_max <- max(selectivity$jacobs, na.rm = TRUE)
y_pad <- 0.05 * (y_max - y_min)
y_limits <- c(y_min - y_pad, y_max + y_pad)

# Plot the boxplot with counts placed below the zero line
ggplot(data = selectivity %>% filter(!is.na(Quartile_Bins.BUI)) %>% filter(!is.na(class.bin)),
       aes(fill = region, x = Quartile_Bins.BUI)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.6, aes(y = jacobs)) +
  geom_text(data = count_data %>% filter(!is.na(Quartile_Bins.BUI)),
            aes(label = count, color = region, y = y_below),
            vjust = 1,
            position = position_dodge(width = 0.75),
            size = 3) +
  labs(x = "", y = "jacobs's index", fill = "Region") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1", guide = FALSE) +  # suppress duplicate legend for text
  facet_wrap(~class.bin, scales = "fixed", ncol = 5) +
  scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(angle = 90))


##### Adjusted sample size location
y_min <- min(selectivity$jacobs, na.rm = TRUE)
y_max <- max(selectivity$jacobs, na.rm = TRUE)
y_pad <- 0.05 * (y_max - y_min)
y_limits <- c(y_min - y_pad, y_max + y_pad)
y_range <- diff(y_limits)

# Calculate counts for each combination of 'class.bin', 'Quartile_Bins.BUI', and 'region'
# and position those sample-size labels a little below the zero line
count_data <- selectivity %>% filter(!is.na(class.bin)) %>%
  group_by(class.bin, Quartile_Bins.BUI, region) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(y_below = 0 - 0.04 * y_range)   # adjust 0.04 to move labels further/closer to zero

# Plot the boxplot with counts placed below the zero line
jacobs.east.west = ggplot(data = selectivity %>% filter(!is.na(Quartile_Bins.BUI)) %>% filter(!is.na(class.bin)),
       aes(fill = region, x = Quartile_Bins.BUI)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.6, aes(y = jacobs)) +
  geom_text(data = count_data %>% filter(!is.na(Quartile_Bins.BUI)),
            aes(label = count, color = region, y = y_below),
            vjust = 1,
            position = position_dodge(width = 0.75), angle = 90,
            size = 3) +
  labs(x = "", y = "jacobs's index", fill = "Region") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1", guide = FALSE) +  # suppress duplicate legend for text
  facet_wrap(~class.bin, scales = "fixed", ncol = 5) +
  scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(angle = 90))

#save
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = jacobs.east.west, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "jacobs_by_class_east_west_boxplot.png")


########### Quartile bin BUI count as alpha in barplot #############
selectivity.mean <- selectivity %>%
  group_by(Quartile_Bins.BUI,class.bin,canopy.bin, region) %>%
  summarise(
    mean = mean(jacobs, na.rm = TRUE),
    median = median(jacobs, na.rm = TRUE),
    sd = sd(jacobs, na.rm = TRUE),
    p05 = quantile(jacobs, probs = 0.05, na.rm = TRUE),
    p95 = quantile(jacobs, probs = 0.95, na.rm = TRUE),
    count = n()
  )

selectivity.mean <- selectivity.mean %>%
  mutate(canopy.bin = forcats::fct_relevel(as.factor(canopy.bin), "total", after = Inf))

selectivity.mean$alpha=round(selectivity.mean$count/max(selectivity.mean$count),3)+0.5

y_min <- min(selectivity.mean$mean, na.rm = TRUE)
y_max <- max(selectivity.mean$mean, na.rm = TRUE)
if (all(c("p05", "p95") %in% names(selectivity.mean))) {
  y_min <- min(y_min, min(selectivity.mean$p05, na.rm = TRUE))
  y_max <- max(y_max, max(selectivity.mean$p95, na.rm = TRUE))
}
y_max = 0.10
y_pad <- 0.05 * (y_max - y_min)
y_limits <- c(y_min - y_pad, y_max + y_pad)

### BUI East-West barplot with counts as alpha
bui.east.west.alpha = ggplot(data = selectivity.mean %>% filter(!is.na(Quartile_Bins.BUI))  %>% filter(!is.na(class.bin)),
       aes(x = Quartile_Bins.BUI, y = mean, fill = region)) +
  geom_bar(aes(alpha = count), stat = "identity", position = "dodge") +
  labs(x = "", y = "jacobs's index", fill = "Region") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  # use fixed scales so y is the same across facets
  facet_grid(class.bin ~ canopy.bin, scales = "fixed") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 5))+ # enforce common y-limits and prettier breaks
  # show alpha legend as 6 categories from 0 to 1250
  scale_alpha_continuous(
    range = c(0.3, 1),
    limits = c(0, 1250),
    breaks = seq(0, 1250, length.out = 6),
    labels = as.integer(seq(0, 1250, length.out = 6)),
    na.value = 0
  )

# save
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = bui.east.west.alpha, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "jacobs_by_class_east_west_boxplot.png")


ggsave(paste0(pathname, "figures/bui_east_west_barplot_alpha_counts.png"), plot = bui.east.west.alpha, width = 12, height = 8)  
  
########## seasonal ############

selectivity.mean<- selectivity %>% group_by(Quartile_Bins.BUI,class.bin,canopy.bin,season.bin) %>% 
  summarise(mean=mean(jacobs,na.rm=T),
            median=median(jacobs,na.rm=T),
            sd=sd(jacobs,na.rm=T),
            p05=quantile(jacobs,probs=0.05,na.rm=T),
            p95=quantile(jacobs,probs=0.95,na.rm=T),
            count=n())

selectivity.mean <- selectivity.mean %>%
  mutate(canopy.bin = forcats::fct_relevel(as.factor(canopy.bin), "total", after = Inf))

selectivity.mean$alpha=round(selectivity.mean$count/max(selectivity.mean$count),3)+0.5

y_min <- min(selectivity.mean$mean, na.rm = TRUE)
y_max <- max(selectivity.mean$mean, na.rm = TRUE)
if (all(c("p05", "p95") %in% names(selectivity.mean))) {
  y_min <- min(y_min, min(selectivity.mean$p05, na.rm = TRUE))
  y_max <- max(y_max, max(selectivity.mean$p95, na.rm = TRUE))
}
y_max = 0.15
y_pad <- 0.05 * (y_max - y_min)
y_limits <- c(y_min - y_pad, y_max + y_pad)

## spring-fall barplot with counts as alpha
bui.spring.fall.alpha = ggplot(data = selectivity.mean %>% filter(!is.na(Quartile_Bins.BUI))  %>% filter(!is.na(class.bin)),
       aes(x = Quartile_Bins.BUI, y = mean, fill = season.bin)) +
  geom_bar(aes(alpha = count), stat = "identity", position = "dodge") +
  labs(x = "", y = "jacobs's index", fill = "Season") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  # use fixed scales so y is the same across facets
  facet_grid(class.bin ~ canopy.bin, scales = "fixed") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 5))+ # enforce common y-limits and prettier breaks
  # show alpha legend as 6 categories from 0 to 1250
  scale_alpha_continuous(
    range = c(0.3, 1),
    limits = c(0, 1250),
    breaks = seq(0, 1250, length.out = 6),
    labels = as.integer(seq(0, 1250, length.out = 6)),
    na.value = 0
  )
#save
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = bui.spring.fall.alpha, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "bui_early_late_barplot_alpha_counts.png")


## simple barplot with counts ##

# get scale of y xis
y_min <- min(selectivity.mean$mean, na.rm = TRUE)
y_max <- max(selectivity.mean$mean, na.rm = TRUE)
if (all(c("p05", "p95") %in% names(selectivity.mean))) {
  y_min <- min(y_min, min(selectivity.mean$p05, na.rm = TRUE))
  y_max <- max(y_max, max(selectivity.mean$p95, na.rm = TRUE))
}

y_max = 0.13 # seems like limit? NA must have limit closer to 0.25
y_pad <- 0.05 * (y_max - y_min)
y_limits <- c(y_min - y_pad, y_max + y_pad)

# early-late simple with counts
bui.early.late.counts = ggplot(data=selectivity.mean %>% filter(!is.na(Quartile_Bins.BUI))  %>% filter(!is.na(class.bin)),
       aes(x=Quartile_Bins.BUI,y=mean,fill=season.bin))+
       geom_bar(stat = "identity", alpha = 0.9,
           position = position_dodge(width = 0.9), width = 0.8) +
  # place count labels on top of each bar (use same dodge width as geom_bar)
  geom_text(aes(label = count, y = mean),
            position = position_dodge(width = 0.9),
            vjust = -0.4,
            color = "black",
            size = 2.5) +
  labs(x = "",
       y = "jacobs's index",
       fill = "Period") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(class.bin ~ canopy.bin, scales = "free") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 5)) # enforce common y-limits and prettier breaks

#save
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = bui.early.late.counts, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "bui_early_late_barplot_simple_counts.png")




# filter just western ecoregion #####
selectivity.mean<- selectivity %>% filter(region == "West") %>% group_by(Quartile_Bins.BUI,class.bin,canopy.bin,season.bin) %>% 
  summarise(mean=mean(jacobs,na.rm=T),
            median=median(jacobs,na.rm=T),
            sd=sd(jacobs,na.rm=T),
            p05=quantile(jacobs,probs=0.05,na.rm=T),
            p95=quantile(jacobs,probs=0.95,na.rm=T),
            count=n())
selectivity.mean <- selectivity.mean %>%
  mutate(canopy.bin = forcats::fct_relevel(as.factor(canopy.bin), "total", after = Inf))

#get scale of y xis
y_min <- min(selectivity.mean$mean, na.rm = TRUE)
y_max <- max(selectivity.mean$mean, na.rm = TRUE)
if (all(c("p05", "p95") %in% names(selectivity.mean))) {
  y_min <- min(y_min, min(selectivity.mean$p05, na.rm = TRUE))
  y_max <- max(y_max, max(selectivity.mean$p95, na.rm = TRUE))
}

y_max = 0.13 # seems like limit? NA must have limit closer to 0.25
y_pad <- 0.05 * (y_max - y_min)
y_limits <- c(y_min - y_pad, y_max + y_pad)

# plot
bui.early.late.west = ggplot(data=selectivity.mean %>% filter(!is.na(Quartile_Bins.BUI))  %>% filter(!is.na(class.bin)),
       aes(x=Quartile_Bins.BUI,y=mean,fill=season.bin))+
       geom_bar(stat = "identity", alpha = 0.9,
           position = position_dodge(width = 0.9), width = 0.8) +
  # place count labels on top of each bar (use same dodge width as geom_bar)
  geom_text(aes(label = count, y = mean),
            position = position_dodge(width = 0.9),
            vjust = -0.4,
            color = "black",
            size = 2.5) +
  labs(x = "",
       y = "jacobs's index",
       fill = "Period") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(class.bin ~ canopy.bin, scales = "free") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 5)) # enforce common y-limits and prettier breaks

# save
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = bui.early.late.west, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "bui_early_late_west.png")

ggsave(paste0(pathname, "figures/bui_early_late_west.png"), plot = bui.early.late.west, width = 12, height = 8)

################ mean jacobss for landcover bui col isi row ##########
selectivity.mean<- selectivity %>% group_by(class.bin,Quartile_Bins.ISI,Quartile_Bins.BUI) %>% 
  summarise(mean=mean(jacobs,na.rm=T),
            median=median(jacobs,na.rm=T),
            sd=sd(jacobs,na.rm=T),
            p05=quantile(jacobs,probs=0.05,na.rm=T),
            p95=quantile(jacobs,probs=0.95,na.rm=T),
            count=n())


mean.jacobs.bui.isi = ggplot(data=selectivity.mean %>% filter(!is.na(Quartile_Bins.ISI)&!is.na(Quartile_Bins.BUI)&!is.na(class.bin)&class.bin!='total'),
       aes(x=as.factor(reorder(class.bin, mean)),y=mean,fill=class.bin))+
  geom_bar(stat = "identity",alpha=0.9)+
    geom_text(aes(label = count), vjust = -0.5, color = "black", size = 3) +  # Adding text labels
  labs(x="BUI quartile",
       y="Mean jacobs's Index",
       fill="Land cover")+
  theme_bw()+
  geom_hline(yintercept = 0)+
  scale_fill_brewer(palette="Set1")+
  facet_grid(Quartile_Bins.ISI~Quartile_Bins.BUI)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + scale_y_continuous(
      sec.axis = sec_axis(
        trans = ~./4,
        name = "ISI quartile")) + 
        theme(
    axis.ticks.y.right = element_blank(), # Remove ticks on the right y-axis
    axis.text.y.right = element_blank()   # Remove labels on the right y-axis
  )

# save
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = mean.jacobs.bui.isi, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "mean_jacobs_bui_isi.png")

### spring-fall boxplot with counts
selectivity.mean<- selectivity %>% filter(!is.na(Quartile_Bins.BUI)) %>% filter(!is.na(class.bin)) %>%group_by(Quartile_Bins.BUI,class.bin,canopy.bin,season.bin) %>% 
  summarise(mean=mean(jacobs,na.rm=T),
            median=median(jacobs,na.rm=T),
            sd=sd(jacobs,na.rm=T),
            p05=quantile(jacobs,probs=0.05,na.rm=T),
            p95=quantile(jacobs,probs=0.95,na.rm=T),
            count=n())

# Calculate counts for each combination of 'class.bin', 'Quartile_Bins.BUI', and 'season.bin'
count_data <- selectivity %>% filter(!is.na(class.bin)) %>%
  group_by(class.bin, Quartile_Bins.BUI, season.bin) %>%
  summarise(count = n(),
            y_pos=quantile(jacobs,0.95,na.rm=T))  %>%
  mutate(y_below = 0 - 0.04 * y_range)   # adjust 0.04 to move labels further/closer to zero

y_min <- min(selectivity$jacobs, na.rm = TRUE)
y_max <- max(selectivity$jacobs, na.rm = TRUE)
y_pad <- 0.05 * (y_max - y_min)

y_max = 1
y_limits <- c(y_min - y_pad, y_max + y_pad)
y_range <- diff(y_limits)

# Plot the boxplot with counts
jacobs.early.late.boxplot = ggplot(data = selectivity %>%filter(!is.na(Quartile_Bins.BUI))  %>% filter(!is.na(class.bin)), aes(fill = season.bin, x = Quartile_Bins.BUI)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.6,aes(y=jacobs)) +
  geom_text(data = count_data %>% filter(!is.na(Quartile_Bins.BUI)), 
            aes(label = count, color=season.bin,y=y_below),
            vjust=1,
            angle=90,
            position=position_dodge(width=0.75),
            size = 3) +
  labs(x = "", y = "jacobs's index", fill = "Season") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1", guide = FALSE) +
  facet_wrap(~class.bin, scales = "free", ncol = 5) +
    scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(angle = 90))

# save
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = jacobs.early.late.boxplot, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "jacobs.early.late.boxplot_with_counts.png")
ggsave(paste0(pathname, "figures/jacobs_early_late_boxplot_with_counts.png"), plot = jacobs.early.late.boxplot, width = 12, height = 6)

### warning message on above plot 
#Warning message:
#Removed 274 rows containing non-finite outside the scale range(`stat_boxplot()`). 
###

## Progression polygon shapefile
progression_data <- st_read(paste0(pathname, "progression_polygons/progression_ClipedToFilledNBAC_withNFIREID.shp"))
extent_data <- st_read(paste0(pathname, "fire_extent/nbac_union_dissolve2_buffer2000m_SpatialJoin_dissolve.shp"))

extent_data<- st_transform(extent_data,crs=3347)
progression_data<-st_transform(progression_data,crs=3347)

#out<- extent_data %>% filter(MaxDate<20230810)  #cut this out because he does another removal of dates below
out.index<-unique(extent_data$NFIREID)

#extent_data <- extent_data %>% filter(MaxDate<20230810) #duplicate? 

lookup_table_extent <- st_drop_geometry(extent_data)
lookup_out <- st_drop_geometry(extent_data)
lookup_progression <-st_drop_geometry(progression_data)


#Sort the progression polygon to keep only the top 80% in size
#vfind the 0.8 quantile of area

threshold<- quantile(lookup_progression$UpdateArea,0.2)

prog.dat <- progression_data %>% filter(UpdateArea>= threshold)


#Buffer the progression polygons by 2000 m

prog.dat <- prog.dat %>% st_buffer(2000)  #changed from 200 m?  

#Remove fire polygons outside of Pontone's map boundaries

prog.temp<-st_transform(prog.dat,crs=4326)

prog.centroids<- st_centroid(prog.temp)

ecozones.data<-st_read(paste0(pathname, "Ecozones/ecozones.shp"))

centroid_ecozones<- st_join(prog.centroids,ecozones.data%>%st_transform(crs=st_crs(prog.centroids)))%>%
  select(K_UniqueID,NFIREID,ZONE_NAME)%>%
  st_drop_geometry()

selectivity<-selectivity%>%merge(centroid_ecozones,by.x=c('K_UniqueID','NFIREID.x'), by.y=c('K_UniqueID','NFIREID'))

## save selectivity with ecozone info
write.csv(selectivity,paste0(pathname,"selectivity_with_ecozone_info.csv"))

#. by ecozone
########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NEXT UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########


selectivity.mean<- selectivity_with_ecozones %>% group_by(Quartile_Bins.BUI,class.bin,ZONE_NAME) %>% 
  summarise(mean=mean(jacobs,na.rm=T),
            median=median(jacobs,na.rm=T),
            sd=sd(jacobs,na.rm=T),
            p05=quantile(jacobs,probs=0.05,na.rm=T),
            p95=quantile(jacobs,probs=0.95,na.rm=T),
            count=n())


# Calculate counts for each combination of 'class.bin', 'Quartile_Bins.BUI', and 'ZONE_NAME'
plot_df <- selectivity_with_ecozones %>%
  filter(!is.na(Quartile_Bins.BUI), !is.na(class.bin), !is.na(ZONE_NAME)) %>%
  filter(!ZONE_NAME %in% c("Atlantic Maritime", "Boreal Cordillera", "Prairie", "Southern Arctic"))

# common y-limits with a small padding
y_min <- min(plot_df$jacobs, na.rm = TRUE)
y_max <- max(plot_df$jacobs, na.rm = TRUE)
y_pad <- 0.05 * (y_max - y_min)
y_limits <- c(y_min - y_pad, y_max + y_pad)
y_range <- diff(y_limits)

# target y for count labels (absolute on jacobs scale)
y_target <- 0.5

# counts + label positions: place at y = 0.5 when that's inside the plot limits,
# otherwise fall back to previous clamped position
# y - target option
count_data <- plot_df %>%
  group_by(class.bin, Quartile_Bins.BUI, ZONE_NAME) %>%
  summarise(count = n(), y_val = quantile(jacobs, 0.90, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_pos = if (y_target >= y_limits[1] && y_target <= y_limits[2]) {
                   y_target
                 } else {
                   pmin(y_val, y_limits[2] - 0.02 * y_range)
                 })

# base option: y_pos based on quantile
# counts + label positions (clamped inside y-limits)
count_data <- plot_df %>%
  group_by(class.bin, Quartile_Bins.BUI, ZONE_NAME) %>%
  summarise(count = n(), y_pos = quantile(jacobs, 0.90, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_pos = pmin(y_pos, y_limits[2] - 0.02 * y_range))


# boxplot with same y scale across facets (geom_text uses y_pos)
box.plot.eco <- ggplot(plot_df, aes(x = Quartile_Bins.BUI, y = jacobs)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.6) +
  geom_text(data = count_data, aes(x = Quartile_Bins.BUI, label = count, y = y_pos),
            vjust = -0.5, position = position_dodge(width = 0.75), size = 3) +
  labs(x = "", y = "jacobs's index", fill = "Period") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(ZONE_NAME ~ class.bin, scales = "fixed") +
  scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(angle = 90))

box.plot.eco

temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = box.plot.eco, width = 10, height = 6)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "bui_quartile_by_class_per_ecozone.png")
########### Sping and Summer East and West ######

#########

is.na(count_data)


# make a bar plot about the total number of fires by class.bin and canopy bin with counts on the y axis
dcount <- selectivity %>%
  filter(!is.na(class.bin)) %>%
  group_by(class.bin, canopy.bin) %>%
  summarise(total_fires = n(), .groups = "drop") %>%
  ungroup()

# use position_dodge2(preserve = "single") so single bars (e.g. upland/total) keep the same width as individual dodged bars
pd <- position_dodge2(width = 0.9, preserve = "single")

total_fires_plot = ggplot(data = dcount,
                          aes(x = class.bin, y = total_fires, fill = canopy.bin)) +
  geom_bar(stat = "identity", position = pd, width = 0.8) +
  labs(x = "Land Cover Type", y = "Total Progression Days", fill = "Forest Cover Type") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  geom_text(aes(label = total_fires), vjust = -0.5, color = "black", size = 3,
            position = pd) +  # use same dodge for labels
  theme(axis.text.x = element_text(angle = 90))

total_fires_plot

tempfile <- tempfile(fileext = ".png")
ggsave(tempfile, plot = total_fires_plot, width = 10, height = 6)
googledrive::drive_upload(media = tempfile, path = jacobs_plots_id, name = "total_progressions_by_class_canopy.png")

#
selectivity.mean<- selectivity %>% group_by(Quartile_Bins.FWI,class.bin,canopy.bin) %>% 
  summarise(mean=mean(jacobs,na.rm=T),
            median=median(jacobs,na.rm=T),
            sd=sd(jacobs,na.rm=T),
            p05=quantile(jacobs,probs=0.05,na.rm=T),
            p95=quantile(jacobs,probs=0.95,na.rm=T),
            count=n())

violin_fwi_class = ggplot(data=selectivity %>% filter(!is.na(Quartile_Bins.FWI) & variable != "water" & variable != 'mineral' & is.finite(jacobs)),
       aes(y=jacobs,fill=Quartile_Bins.FWI,x=Quartile_Bins.FWI))+
  geom_violin(alpha=0.5)+
    #geom_text(aes(label = count), vjust = -0.5, color = "black", size = 3) +  # Adding text labels
  labs(x="",
       y="jacobs's index",
       fill="FWI quartiles")+
  theme_bw()+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values=my_colors)+
  facet_grid(~class.bin,scales="free")+
  theme(axis.text.x = element_text(angle = 90, size = 20),
        axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22))

temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, plot = violin_fwi_class, width = 20, height = 20)
googledrive::drive_upload(media = temp_file, path = jacobs_plots_id, name = "violin_fwi_class.png")


# for each fire progression, calculate where fires burned at the same rate as adjacent updates
# from selectivity, calculate burn rate
# br = burned area/ avaiable area

selectivity_with_ecozones <- selectivity_with_ecozones |> 
  mutate(br = burned.area / av.area)

# find progressions where br is equal to uplands
# get just the totals
burn_rates_by_fire <- selectivity_with_ecozones |> 
  filter(canopy.bin == "total") |> 
  group_by(NFIREID.x) |>  
  summarise(ba_peatland = sum(burned.area, na.rm = TRUE),
            ba_upland = sum(burned.area[class.bin == "upland"], na.rm = TRUE),
            ba_peatland = (ba_peatland - ba_upland), 
            av_peatland = sum(av.area[class.bin == "peatland"], na.rm = TRUE),
            av_upland = sum(av.area[class.bin == "upland"], na.rm = TRUE),
            av_peatland = (av_peatland - av_upland),
            br_peatland = ba_peatland / av_peatland,
            br_upland = ba_upland / av_upland,
            .groups = "drop") 

# find fires where br_peatland is approximately equal to br_upland (e.g. within 10% of each other)
similar_burn_rates <- burn_rates_by_fire |> 
  filter(abs(br_peatland - br_upland) / ((br_peatland + br_upland) / 2) < 0.1)

# save
tempfile <- tempfile(fileext = ".csv")
write.csv(similar_burn_rates, tempfile, row.names = FALSE)
googledrive::drive_upload(media = tempfile, path = summary_stats_id, name = "similar_burn_rates_by_fire_all_peat_within10pct.csv")


burn_rates_by_progression <- selectivity_with_ecozones |> 
  filter(canopy.bin == "total") |> 
  group_by(K_UniqueID) |>  
  summarise(ba_peatland = sum(burned.area, na.rm = TRUE),
            ba_upland = sum(burned.area[class.bin == "upland"], na.rm = TRUE),
            ba_peatland = (ba_peatland - ba_upland), 
            av_peatland = sum(av.area[class.bin == "peatland"], na.rm = TRUE),
            av_upland = sum(av.area[class.bin == "upland"], na.rm = TRUE),
            av_peatland = (av_peatland - av_upland),
            br_peatland = ba_peatland / av_peatland,
            br_upland = ba_upland / av_upland,
            .groups = "drop") 

# find fires where br_peatland is approximately equal to br_upland (e.g. within 10% of each other)
similar_burn_rates_prog <- burn_rates_by_progression |> 
  filter(abs(br_peatland - br_upland) / ((br_peatland + br_upland) / 2) < 0.1)

# save
tempfile <- tempfile(fileext = ".csv")
write.csv(similar_burn_rates_prog, tempfile, row.names = FALSE)
googledrive::drive_upload(media = tempfile, path = summary_stats_id, name = "similar_burn_rates_by_progression_all_peat_within10pct.csv")

#what fires showed a preference for burning in peatlands?  Let's chunk these fires into those three bins: avoidance, indifference, preference.
# based on jacobs index
selectivity_with_ecozones <- selectivity_with_ecozones |> 
  mutate(preference_bin = case_when(
    jacobs < -0.33 ~ "avoidance",
    jacobs > 0.33 ~ "preference",
    TRUE ~ "indifference"
  ))

#