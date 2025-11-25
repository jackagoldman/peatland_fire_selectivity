#############################
# Fire Selectivity median vs high selectivity comparison per landcover and assocaited drivers
# Jack A. Goldman
################################
# required libraries
library(tidyverse)
library(ggplot2)
library(lqmm)
library(purrr)
library(patchwork)

# reasoning

# fire selectivity (Chesson), these show the "typical" selectivity (50%) 
#and the "high" selectivity threshold (90%) across fires for a landcover type.
# Reasoning and Rationale for Using 50% and 90%
#Why these specific quartiles? They capture different parts of the selectivity distribution:
#50% (Median): Represents the "average" or most common fire behavior. 
#It's robust to outliers and shows what happens in typical fires. For fire selectivity, it answers: "How selective are most fires for this landcover?" 
#If the median is low (e.g., near 0), most fires are neutral or avoidant.
#90%: Focuses on the upper tail—the rare, highly selective fires. 
#It answers: "What drives the most extreme selectivity?" 
#This is useful for risk assessment, as it highlights conditions leading to strong fire preference (e.g., in dry weather).
#Rationale in fire ecology: Fires aren't uniform; some are mild, others extreme. 
#Analyzing median (50%) vs. high-end (90%) reveals if factors like weather (ISI) affect all fires similarly or only the selective ones. 
#If selectivity varies by quantile, it suggests heterogeneity (e.g., dry conditions amplify selectivity only for certain fires).
#Why not just the mean? Means can be skewed by outliers; quartiles are more reliable for non-normal data like selectivity scores.


#set file path for working directory
pathname = "/Users/jgoldman/Library/CloudStorage/GoogleDrive-jandrewgoldman@gmail.com/My Drive/3_POSTDOC/projects/selectivity/"

# read in selectivity data
selectivity = read_csv(paste0(pathname, 'Dataset_selectivity_Nov11_KM_notLog_ALL.csv'))

# turn landcover variable into a factor
selectivity$fvariable = as.factor(selectivity$variable)


# what is NFIREID and KUNIQUE ID
selectivity$NFIREID.x
selectivity$K_UniqueID


# prepare df_lq
df_lq <- selectivity %>%
  filter(!is.na(variable),
         !is.na(K_UniqueID), !is.na(NFIREID.x),
         is.finite(Chesson), is.finite(BUI.mean)) %>%
  mutate(fvariable = as.factor(variable)) %>%
  droplevels()

# minimum obs per landcover to attempt lqmm
min_n <- 10

# keep only sufficiently sampled landcovers
keep_lvls <- df_lq %>% count(fvariable) %>% filter(n >= min_n) %>% pull(fvariable)

# calculate and compare median and 90th quantile selectivity per landcover type
selectivity %>%
  group_by(fvariable) %>%
  summarise(
    median_sel = quantile(Chesson, 0.5, na.rm = TRUE),
    q90_sel = quantile(Chesson, 0.9, na.rm = TRUE)
  ) %>%
  arrange(desc(median_sel))

# For each landcover class, fit LQMM at taus 0.5 and 0.9
all_mods <- list()
for (lvl in keep_lvls) {
  cat("\n\nFitting models for landcover:", lvl, "\n")
  
  taus <- c(0.5, 0.9)
  mods <- mclapply(taus, function(tau) {
    lqmm::lqmm(fixed = Chesson ~ BUI.mean + ISI.mean,
               random = ~ 1 | NFIREID.x/K_UniqueID,
               group = NFIREID.x,
               data = df_lq[df_lq$fvariable == lvl, ],
               nK = 10,
               tau = tau,
               control = lqmmControl(method = "df", startQR = TRUE))
  }, mc.cores = 2)
  
  # Store results by landcover name
  all_mods[[lvl]] <- mods
  
  cat("Summary for tau = 0.5:\n")
  print(summary(mods[[1]]))
  cat("\nSummary for tau = 0.9:\n")
  print(summary(mods[[2]]))
}
# all_mods is a list where each element is named by fvariable, containing a list of two models (tau 0.1 and 0.9)

summary(all_mods[["forested_bog"]][[1]])

# save all models to RDS in /results
path_wd = getwd()
saveRDS(all_mods, file = paste0(path_wd, "/results/lqmm_median_vs_high_selectivity_models.rds"))

# Extract fixed effects into a table
fixed_effects_table <- data.frame()
taus <- c(0.5, 0.9)

for (lvl in names(all_mods)) {
  for (i in 1:2) {
    mod <- all_mods[[lvl]][[i]]
    tau_val <- taus[i]
    
    # Get the summary table
    summ <- summary(mod)
    ttable <- summ$tTable
    
    # Convert to data frame
    df_temp <- as.data.frame(ttable)
    df_temp$Coefficient <- rownames(ttable)
    df_temp$Landcover <- lvl
    df_temp$Tau <- tau_val
    
    # Rename columns to match request
    colnames(df_temp) <- c("Value", "Std.Error", "CI_Lower", "CI_Upper", "P", "Coefficient", "Landcover", "Tau")
    
    # Reorder columns
    df_temp <- df_temp[, c("Landcover", "Tau", "Coefficient", "Value", "Std.Error", "CI_Lower", "CI_Upper", "P")]
    
    # Append to main table
    fixed_effects_table <- rbind(fixed_effects_table, df_temp)
  }
}

# View the table
print(fixed_effects_table)

# save to csv in /results
write.csv(fixed_effects_table, file = paste0(path_wd, "/results/lqmm_median_vs_high_selectivity_fixed_effects.csv"), row.names = FALSE)


# Create individual ggplot plots for each landcover showing slopes for tau 0.5 and 0.9
plot_list <- list()

for (lvl in names(all_mods)) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  
  # Compute mean BUI.mean for fixing
  mean_BUI <- mean(dat$BUI.mean, na.rm = TRUE)
  
  # Prepare data for lines
  line_data <- data.frame()
  for (i in 1:2) {
    mod <- all_mods[[lvl]][[i]]
    tau_val <- taus[i]
    
    if (inherits(mod, "error") || is.null(mod)) next
    
    coefs <- coef(mod)
    if (any(is.na(coefs))) next
    
    intercept_adj <- coefs[1] + coefs[2] * mean_BUI
    slope <- coefs[3]
    
    temp <- data.frame(Tau = tau_val, Intercept = intercept_adj, Slope = slope)
    line_data <- rbind(line_data, temp)
  }
  
  # Create plot
  # Create line data for geom_line
  x_range <- seq(min(dat$ISI.mean, na.rm = TRUE), max(dat$ISI.mean, na.rm = TRUE), length = 100)
  line_pred <- data.frame()
  for (i in 1:nrow(line_data)) {
    tau_val <- line_data$Tau[i]
    int <- line_data$Intercept[i]
    sl <- line_data$Slope[i]
    y_vals <- int + sl * x_range
    temp <- data.frame(x = x_range, y = y_vals, Tau = tau_val)
    line_pred <- rbind(line_pred, temp)
  }
  
  p <- ggplot(dat, aes(x = ISI.mean, y = Chesson)) +
    geom_point(alpha = 0.6) +
    geom_line(data = line_pred, aes(x = x, y = y, color = factor(Tau)), linewidth = 1) +
    scale_color_manual(values = c("0.5" = "red", "0.9" = "blue"), name = "Tau", 
                       guide = guide_legend(direction = "horizontal")) +
    labs(title = (lvl),
         x = "ISI.mean", y = "Chesson's Index") +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "bottom")
  
  # Add to list
  plot_list[[lvl]] <- p
}

# Example: print one plot
print(plot_list[["forested_bog"]])

levels(df_lq$fvariable)

# Create a grid of all plots except "mineral"
plots_to_grid <- plot_list[names(plot_list) != "mineral"]
if (length(plots_to_grid) > 0) {
  grid_plot <- wrap_plots(plots_to_grid)
  print(grid_plot)
}

# To save all plots, you can loop: for (name in names(plot_list)) { ggsave(filename = paste0(name, "_slopes.png"), plot = plot_list[[name]]) }


