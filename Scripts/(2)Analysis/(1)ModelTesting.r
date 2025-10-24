
setwd("G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\Plot_Level_Datasets")

con_all_results = read.csv("VariableVotes_Conifer.csv")
shrub_all_results = read.csv("G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\Plot_Level_Datasets\\VariableVotes_Shrubs.csv")
oak_all_results = read.csv("G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\Plot_Level_Datasets\\VariableVotes_Oak.csv")

####################### PLOT METRICS ###################################


con_all_results <- arrange(con_all_results, desc(significant_count))

sub_variable = con_all_results[-c(2,3,5,6,7,10:16),]

# Order variables
top_sig_summary <- sub_variable %>%
  arrange(desc(significant_count)) %>%
  mutate(sig_vars = factor(sig_vars, levels = rev(sig_vars)))

top_sig_summary = top_sig_summary[c(1:10),]

# Build a long-format dataframe for clean aesthetics + legend
plot_df <- top_sig_summary %>%
  select(sig_vars, mean_aic, mean_auc, mean_kappa) %>%
  mutate(
    auc_scaled = mean_auc * 100,
    kappa_scaled = mean_kappa * 100
  ) %>%
  pivot_longer(
    cols = c(mean_aic, auc_scaled, kappa_scaled),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    mean_aic = "AIC",
                    auc_scaled = "AUC",
                    kappa_scaled = "Kappa")
  )

# Define color and shape
metric_colors <- c("AUC" = "darkgreen", "Kappa" = "purple", "AIC" = "red")
metric_shapes <- c("AUC" = 16, "Kappa" = 17, "AIC" = 15)

# Final plot
c <- ggplot(plot_df, aes(x = value, y = sig_vars, color = metric, shape = metric)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = metric_colors) +
  scale_shape_manual(values = metric_shapes) +
  scale_x_continuous(
    name = "AIC",
    sec.axis = sec_axis(~./100, name = "AUC / Kappa")
  ) +
  labs(
    y = "Variable",
    color = " ",
    shape = " ",
    title = "         Conifer"
  ) +
  theme_minimal(base_size = 16) +  # Increase base font size
  theme(
    axis.title.y = element_text(margin = margin(r = 15), size = 16),
    axis.title.x.top = element_text(color = "black", size = 16),
    axis.title.x.bottom = element_text(color = "red", size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

c
#---------------------     OAK      ---------------------


oak_all_results <- arrange(oak_all_results, desc(significant_count))

sub_variable = oak_all_results[-c(3:22, 24:29),]

# Order variables
top_sig_summary <- sub_variable %>%
  arrange(desc(significant_count)) %>%
  mutate(sig_vars = factor(sig_vars, levels = rev(sig_vars)))

top_sig_summary = top_sig_summary[c(1:10),]

# Build a long-format dataframe for clean aesthetics + legend
plot_df <- top_sig_summary %>%
  select(sig_vars, mean_aic, mean_auc, mean_kappa) %>%
  mutate(
    auc_scaled = mean_auc * 100,
    kappa_scaled = mean_kappa * 100
  ) %>%
  pivot_longer(
    cols = c(mean_aic, auc_scaled, kappa_scaled),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    mean_aic = "AIC",
                    auc_scaled = "AUC",
                    kappa_scaled = "Kappa")
  )

# Define color and shape
metric_colors <- c("AUC" = "darkgreen", "Kappa" = "purple", "AIC" = "red")
metric_shapes <- c("AUC" = 16, "Kappa" = 17, "AIC" = 15)

# Final plot
o <- ggplot(plot_df, aes(x = value, y = sig_vars, color = metric, shape = metric)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = metric_colors) +
  scale_shape_manual(values = metric_shapes) +
  scale_x_continuous(
    name = "AIC",
    sec.axis = sec_axis(~./100, name = "AUC / Kappa")
  ) +
  labs(
    y = "Variable",
    color = " ",
    shape = " ",
    title = "         Oak"
  ) +
  theme_minimal(base_size = 16) +  # Increase base font size
  theme(
    axis.title.y = element_text(margin = margin(r = 15), size = 16),
    axis.title.x.top = element_text(color = "black", size = 16),
    axis.title.x.bottom = element_text(color = "red", size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

o


################ SHRUB

shrub_all_results <- arrange(shrub_all_results, desc(significant_count))

sub_variable = shrub_all_results[-c(3:22, 24:29),]

# Order variables
top_sig_summary <- sub_variable %>%
  arrange(desc(significant_count)) %>%
  mutate(sig_vars = factor(sig_vars, levels = rev(sig_vars)))

top_sig_summary = top_sig_summary[c(1:10),]

# Build a long-format dataframe for clean aesthetics + legend
plot_df <- top_sig_summary %>%
  select(sig_vars, mean_aic, mean_auc, mean_kappa) %>%
  mutate(
    auc_scaled = mean_auc * 100,
    kappa_scaled = mean_kappa * 100
  ) %>%
  pivot_longer(
    cols = c(mean_aic, auc_scaled, kappa_scaled),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    mean_aic = "AIC",
                    auc_scaled = "AUC",
                    kappa_scaled = "Kappa")
  )

# Define color and shape
metric_colors <- c("AUC" = "darkgreen", "Kappa" = "purple", "AIC" = "red")
metric_shapes <- c("AUC" = 16, "Kappa" = 17, "AIC" = 15)

# Final plot
s <- ggplot(plot_df, aes(x = value, y = sig_vars, color = metric, shape = metric)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = metric_colors) +
  scale_shape_manual(values = metric_shapes) +
  scale_x_continuous(
    name = "AIC",
    sec.axis = sec_axis(~./100, name = "AUC / Kappa")
  ) +
  labs(
    y = "Variable",
    color = " ",
    shape = " ",
    title = "         Shrub"
  ) +
  theme_minimal(base_size = 16) +  # Increase base font size
  theme(
    axis.title.y = element_text(margin = margin(r = 15), size = 16),
    axis.title.x.top = element_text(color = "black", size = 16),
    axis.title.x.bottom = element_text(color = "red", size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

s

#--------------------------- ALL TOGETHER

library(patchwork)

# Combine all three plots side by side
c + o + s +
  plot_layout(ncol = 3) +
  plot_annotation(title = "Top Significant Variables: AIC, AUC & Kappa")







