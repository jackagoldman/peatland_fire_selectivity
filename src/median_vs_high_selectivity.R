#############################
# Fire Selectivity median vs high selectivity comparison per landcover and assocaited drivers
# Jack A. Goldman
################################
# required libraries
library(tidyverse)
library(ggplot2)
library(lqmm)
library(purrr)
library(cowplot)
library(paletteer)
# reasoning

# Fire selectivity (Chesson) at 0.5 and 0.9 quantiles show the "typical" selectivity (50%) 
# and the "high" selectivity threshold (90%) across fires for a landcover type.
# Reasoning and Rationale for Using 50% and 90% - 
# These specific quartiles capture different parts of the selectivity distribution:
# 50% (Median): Represents the "average" or most common fire selectivity behaviour. 
# It's robust to outliers and shows what happens in typical fires. 
# For fire selectivity, it answers: "How selective are most fires for this landcover?" 
# If the median is low (e.g., near 0), most fires are neutral or avoidant.
# 90%: Focuses on the upper tail which are the rare, highly selective fires. 
# It answers: "What drives the most extreme selectivity?" 
# This is useful for risk assessment, as it highlights conditions leading to strong fire preference (e.g., in dry weather).
# Rationale in fire ecology: Fires aren't uniform; some are mild, others extreme. 
# Analyzing median (50%) vs. high-end (90%) reveals if factors like weather (ISI) affect all fires similarly or only the selective ones. 
# If selectivity varies by quantile, it suggests heterogeneity (e.g., dry conditions amplify selectivity only for certain fires).
# Why not just the mean? Means can be skewed by outliers; quartiles are more reliable for non-normal data like selectivity scores.


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

# check one model for tau = 0.9
summary(all_mods[["open_bog"]][[2]])

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

# ISI plots
# Create individual ggplot plots for each landcover showing slopes for tau 0.5 and 0.9
plot_list <- list()

for (lvl in names(all_mods)) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  
  if (nrow(dat) == 0) next
  
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
  # make title name by removing underscors from lvl
  title_name <- gsub("_", " ", lvl)
  
  # Get p-values for ISI.mean
  p_vals <- fixed_effects_table %>%
    filter(Landcover == lvl, Coefficient == "ISI.mean") %>%
    select(Tau, P) %>%
    deframe()
  
  p_text <- paste0("italic('T 0.5 p = ", round(p_vals["0.5"], 3), "')\nitalic('T 0.9 p = ", round(p_vals["0.9"], 3), "')")
  
  p <- ggplot(dat, aes(x = ISI.mean, y = Chesson)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_line(data = line_pred, aes(x = x, y = y, color = factor(Tau)), linewidth = 3) +
    scale_color_manual(values = c("0.5" = "#BF6860FF", "0.9" = "#37738DFF"), name = "Tau", 
                       guide = guide_legend(direction = "horizontal")) +
    labs(title = (title_name),
         x = "Mean ISI", y = "Chesson's Index") +
    annotate("text", x = max(dat$ISI.mean, na.rm = TRUE), y = max(dat$Chesson, na.rm = TRUE), 
             label = paste0("italic('T 0.5 p = ", ifelse(p_vals["0.5"] < 0.0001, "< 0.0001", round(p_vals["0.5"], 3)), "')"), hjust = 1, vjust = 1, size = 10, parse = TRUE) +
    annotate("text", x = max(dat$ISI.mean, na.rm = TRUE), y = max(dat$Chesson, na.rm = TRUE) - 0.08 * diff(range(dat$Chesson, na.rm = TRUE)), 
             label = paste0("italic('T 0.9 p = ", ifelse(p_vals["0.9"] < 0.0001, "< 0.0001", round(p_vals["0.9"], 3)), "')"), hjust = 1, vjust = 1, size = 10, parse = TRUE) +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 3),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          plot.margin = margin(10, 10, 10, 10))
  
  # Add to list
  plot_list[[lvl]] <- p
}

# print one plot
print(plot_list[["forested_bog"]])

# save ISI to /plots/ISI
for (name in names(plot_list)) { ggsave(filename = paste0("plots/ISI/", name, "_t5_t9.png"), plot = plot_list[[name]]) }


levels(df_lq$fvariable)

# Create a grid of all plots except "mineral" and all the start with total_
plot_list_sub <- plot_list[names(plot_list) != "mineral" & !grepl("^total_", names(plot_list))]



# Combine data from all landcovers into one dataframe
combined_data <- do.call(rbind, lapply(names(plot_list_sub), function(lvl) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  dat$landcover <- lvl
  dat
}))

# Prepare line predictions for all landcovers
line_data_all <- data.frame()
for (lvl in names(plot_list_sub)) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  mean_BUI <- mean(dat$BUI.mean, na.rm = TRUE)
  x_range <- seq(min(dat$ISI.mean, na.rm = TRUE), max(dat$ISI.mean, na.rm = TRUE), length = 100)
  
  for (i in 1:2) {
    mod <- all_mods[[lvl]][[i]]
    tau_val <- taus[i]
    if (inherits(mod, "error") || is.null(mod)) next
    coefs <- coef(mod)
    if (any(is.na(coefs))) next
    intercept_adj <- coefs[1] + coefs[2] * mean_BUI
    slope <- coefs[3]
    y_vals <- intercept_adj + slope * x_range
    temp <- data.frame(x = x_range, y = y_vals, Tau = tau_val, landcover = lvl)
    line_data_all <- rbind(line_data_all, temp)
  }
}

# Prepare annotation data for p-values
annotation_data <- data.frame()
for (lvl in names(plot_list_sub)) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  p_vals <- fixed_effects_table %>%
    filter(Landcover == lvl, Coefficient == "ISI.mean") %>%
    select(Tau, P) %>%
    deframe()
  
  # For tau 0.5
  p_05 <- ifelse(p_vals["0.5"] < 0.0001, "< 0.0001", round(p_vals["0.5"], 3))
  ann1 <- data.frame(landcover = lvl, 
                     x = max(dat$ISI.mean, na.rm = TRUE), 
                     y = max(dat$Chesson, na.rm = TRUE), 
                     label = paste0("italic('T 0.5 p = ", p_05, "')"))
  
  # For tau 0.9
  p_09 <- ifelse(p_vals["0.9"] < 0.0001, "< 0.0001", round(p_vals["0.9"], 3))
  ann2 <- data.frame(landcover = lvl, 
                     x = max(dat$ISI.mean, na.rm = TRUE), 
                     y = max(dat$Chesson, na.rm = TRUE) - 0.08 * diff(range(dat$Chesson, na.rm = TRUE)), 
                     label = paste0("italic('T 0.9 p = ", p_09, "')"))
  
  annotation_data <- rbind(annotation_data, ann1, ann2)
}

# Create faceted plot
faceted_plot <- ggplot(combined_data, aes(x = ISI.mean, y = Chesson)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_line(data = line_data_all, aes(x = x, y = y, color = factor(Tau)), linewidth = 3) +
  geom_text(data = annotation_data, aes(x = x, y = y, label = label), hjust = 1, vjust = 1, size = 3, parse = TRUE, inherit.aes = FALSE) +
  scale_color_manual(values = c("0.5" = "#BF6860FF", "0.9" = "#37738DFF"), name = "Tau") +
  labs(x = "Mean ISI", y = "Chesson's Index") +
  facet_wrap(~ landcover, scales = "free") +  # Adjust scales as needed
  theme_bw() +
  theme(strip.text = element_text(size = 12))

# Print the faceted plot (without saving)
print(faceted_plot)

# save faceted plot
ggsave(filename = paste0(path_wd, "/plots/ISI/faceted_ISI_median_vs_high_selectivity.png"), plot = faceted_plot, width = 12, height = 20, dpi = 300)


# Recreate the same plotting workflow but for BUI instead of ISI

# BUI plots
# Create individual ggplot plots for each landcover showing slopes for tau 0.5 and 0.9
plot_list_BUI <- list()

for (lvl in names(all_mods)) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  
  if (nrow(dat) == 0) next
  
  # Compute mean ISI.mean for fixing
  mean_ISI <- mean(dat$ISI.mean, na.rm = TRUE)
  
  # Prepare data for lines
  line_data <- data.frame()
  for (i in 1:2) {
    mod <- all_mods[[lvl]][[i]]
    tau_val <- taus[i]
    
    if (inherits(mod, "error") || is.null(mod)) next
    
    coefs <- coef(mod)
    if (any(is.na(coefs))) next
    
    intercept_adj <- coefs[1] + coefs[3] * mean_ISI
    slope <- coefs[2]
    
    temp <- data.frame(Tau = tau_val, Intercept = intercept_adj, Slope = slope)
    line_data <- rbind(line_data, temp)
  }
  
  # Create plot
  # Create line data for geom_line
  x_range <- seq(min(dat$BUI.mean, na.rm = TRUE), max(dat$BUI.mean, na.rm = TRUE), length = 100)
  line_pred <- data.frame()
  for (i in 1:nrow(line_data)) {
    tau_val <- line_data$Tau[i]
    int <- line_data$Intercept[i]
    sl <- line_data$Slope[i]
    y_vals <- int + sl * x_range
    temp <- data.frame(x = x_range, y = y_vals, Tau = tau_val)
    line_pred <- rbind(line_pred, temp)
  }
  # make title name by removing underscors from lvl
  title_name <- gsub("_", " ", lvl)
  
  # Get p-values for BUI.mean
  p_vals <- fixed_effects_table %>%
    filter(Landcover == lvl, Coefficient == "BUI.mean") %>%
    select(Tau, P) %>%
    deframe()
  
  p_text <- paste0("italic('T 0.5 p = ", round(p_vals["0.5"], 3), "')\nitalic('T 0.9 p = ", round(p_vals["0.9"], 3), "')")
  
  p <- ggplot(dat, aes(x = BUI.mean, y = Chesson)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_line(data = line_pred, aes(x = x, y = y, color = factor(Tau)), linewidth = 3) +
    scale_color_manual(values = c("0.5" = "#BF6860FF", "0.9" = "#37738DFF"), name = "Tau", 
                       guide = guide_legend(direction = "horizontal")) +
    labs(title = (title_name),
         x = "Mean BUI", y = "Chesson's Index") +
    annotate("text", x = max(dat$BUI.mean, na.rm = TRUE), y = max(dat$Chesson, na.rm = TRUE), 
             label = paste0("italic('T 0.5 p = ", ifelse(p_vals["0.5"] < 0.0001, "< 0.0001", round(p_vals["0.5"], 3)), "')"), hjust = 1, vjust = 1, size = 10, parse = TRUE) +
    annotate("text", x = max(dat$BUI.mean, na.rm = TRUE), y = max(dat$Chesson, na.rm = TRUE) - 0.08 * diff(range(dat$Chesson, na.rm = TRUE)), 
             label = paste0("italic('T 0.9 p = ", ifelse(p_vals["0.9"] < 0.0001, "< 0.0001", round(p_vals["0.9"], 3)), "')"), hjust = 1, vjust = 1, size = 10, parse = TRUE) +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 3),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          plot.margin = margin(10, 10, 10, 10))
  
  # Add to list
  plot_list_BUI[[lvl]] <- p
}

# print one plot
print(plot_list_BUI[["forested_bog"]])

# save BUI to /plots/BUI
for (name in names(plot_list_BUI)) { ggsave(filename = paste0("plots/BUI/", name, "_t5_t9.png"), plot = plot_list_BUI[[name]]) }


levels(df_lq$fvariable)

# Create a grid of all plots except "mineral" and all the start with total_
plot_list_sub_BUI <- plot_list_BUI[names(plot_list_BUI) != "mineral" & !grepl("^total_", names(plot_list_BUI))]

n_plots <- length(plot_list_sub_BUI)
ncol <- 3
nrow <- ceiling(n_plots / ncol)
panel_plot_BUI <- plot_grid(plotlist = plot_list_sub_BUI, 
                        ncol = ncol,
                        rel_widths = rep(1, ncol),
                        rel_heights = rep(1, nrow),
                        label_size = 12,
                        align = 'none')

# Print the panel plot
print(panel_plot_BUI)

# Combine data from all landcovers into one dataframe
combined_data_BUI <- do.call(rbind, lapply(names(plot_list_sub_BUI), function(lvl) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  dat$landcover <- lvl
  dat
}))

# Prepare line predictions for all landcovers
line_data_all_BUI <- data.frame()
for (lvl in names(plot_list_sub_BUI)) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  mean_ISI <- mean(dat$ISI.mean, na.rm = TRUE)
  x_range <- seq(min(dat$BUI.mean, na.rm = TRUE), max(dat$BUI.mean, na.rm = TRUE), length = 100)
  
  for (i in 1:2) {
    mod <- all_mods[[lvl]][[i]]
    tau_val <- taus[i]
    if (inherits(mod, "error") || is.null(mod)) next
    coefs <- coef(mod)
    if (any(is.na(coefs))) next
    intercept_adj <- coefs[1] + coefs[3] * mean_ISI
    slope <- coefs[2]
    y_vals <- intercept_adj + slope * x_range
    temp <- data.frame(x = x_range, y = y_vals, Tau = tau_val, landcover = lvl)
    line_data_all_BUI <- rbind(line_data_all_BUI, temp)
  }
}

# Prepare annotation data for p-values
annotation_data_BUI <- data.frame()
global_max_BUI <- max(df_lq$BUI.mean, na.rm = TRUE)
for (lvl in names(plot_list_sub_BUI)) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  p_vals <- fixed_effects_table %>%
    filter(Landcover == lvl, Coefficient == "BUI.mean") %>%
    select(Tau, P) %>%
    deframe()
  
  # For tau 0.5
  p_05 <- ifelse(p_vals["0.5"] < 0.0001, "< 0.0001", round(p_vals["0.5"], 3))
  ann1 <- data.frame(landcover = lvl, 
                     x = global_max_BUI, 
                     y = 0.85, 
                     label = paste0("italic('T 0.5 p = ", p_05, "')"))
  
  # For tau 0.9
  p_09 <- ifelse(p_vals["0.9"] < 0.0001, "< 0.0001", round(p_vals["0.9"], 3))
  ann2 <- data.frame(landcover = lvl, 
                     x = global_max_BUI, 
                     y = 0.85 - 0.08 * 0.85, 
                     label = paste0("italic('T 0.9 p = ", p_09, "')"))
  
  annotation_data_BUI <- rbind(annotation_data_BUI, ann1, ann2)
}

# Create faceted plot
faceted_plot_BUI <- ggplot(combined_data_BUI, aes(x = BUI.mean, y = Chesson)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_line(data = line_data_all_BUI, aes(x = x, y = y, color = factor(Tau)), linewidth = 3) +
  geom_text(data = annotation_data_BUI, aes(x = x, y = y, label = label), hjust = 1, vjust = 1, size = 3, parse = TRUE, inherit.aes = FALSE) +
  scale_color_manual(values = c("0.5" = "#BF6860FF", "0.9" = "#37738DFF"), name = "Tau") +
  labs(x = "Mean BUI", y = "Chesson's Index") +
  facet_wrap(~ landcover, scales = "fixed") +
  coord_cartesian(ylim = c(0, 0.85)) +
  theme_bw() +
  theme(strip.text = element_text(size = 12))

# Print the faceted plot (without saving)
print(faceted_plot_BUI)

# save faceted plot
ggsave(filename = paste0(path_wd, "/plots/BUI/faceted_BUI_median_vs_high_selectivity.png"), plot = faceted_plot_BUI, width = 12, height = 20, dpi = 300)


# get residuals and compute pinball loss for all landcovers at tau = 0.5 and 0.9
pinball_df <- data.frame()
for (lvl in names(all_mods)) {
  dat <- df_lq[df_lq$fvariable == lvl, ]
  for (i in 1:2) {
    tau_val <- taus[i]
    mod <- all_mods[[lvl]][[i]]
    if (inherits(mod, "error") || is.null(mod)) next
    
    p_pred <- predict(mod)
    obs <- dat$Chesson
    
    resid <- obs - p_pred
    
    # Compute pinball loss
    loss <- ifelse(obs >= p_pred,
                   tau_val * abs(obs - p_pred),
                   (1 - tau_val) * abs(obs - p_pred))
    
    avg_loss <- mean(loss, na.rm = TRUE)
    
    temp_df <- data.frame(Landcover = lvl, Tau = tau_val, Avg_Pinball_Loss = avg_loss)
    pinball_df <- rbind(pinball_df, temp_df)
  }
}

print(pinball_df)

#save to csv in /results
write.csv(pinball_df, file = paste0(path_wd, "/results/lqmm_median_vs_high_selectivity_pinball_loss.csv"), row.names = FALSE)

