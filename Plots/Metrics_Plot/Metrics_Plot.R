setwd("./Plots/Metrics_Plot")
library(tidyverse)
library(RColorBrewer)
# Create plot df
plot_df_mp <- data.frame(
  Model = c("HiDRA", "HiDRA", "HiDRA", "HiDRA", "PaccMann"),
  DR = c("256-bit Morgan Fingerprints", "1024-bit Morgan Fingerprints", 
         "2048-bit Morgan Fingerprints", "PubChem Fingerprints", "SMILES"),
  Split = c("MP", "MP", "MP","MP","MP"),
  RMSE = c(0.985, 0.982, 0.976, 0.974, 1.137),
  PCC = c(0.935, 0.934, 0.935 ,0.935, 0.913)
)

plot_df_mp <- plot_df_mp %>%
  mutate(`Scaled RMSE` = RMSE / max(RMSE))

# Melt the data frame to long format for ggplot2
plot_df_mp_long <- plot_df_mp %>%
  select(Model, DR, `Scaled RMSE`, PCC) %>%
  pivot_longer(cols = c(`Scaled RMSE`, PCC), names_to = "Metric", values_to = "Value")

plot_df_mp_long <- plot_df_mp_long %>%
  mutate(DR_Model = paste(Model, DR, sep = " "))

plot_df_mc <- data.frame(
  Model = c("PaccMann"),
  DR = c("SMILES"),
  Split = c("MC"),
  RMSE = c(1.316),
  PCC = c(0.878)
)

plot_df_md <- data.frame(
  Model = c("HiDRA", "HiDRA", "SRMF", "SRMF", "SRMF", "ADRML", "ADRML", "ADRML"),
  DR = c("512-bit Morgan Fingerprints", "PubChem Fingerprints", 
         "256-bit Morgan Fingerprints", "1024-bit Morgan Fingerprints", 
         "PubChem Fingerprints", "256-bit Morgan Fingerprints", 
         "512-bit Morgan Fingerprints", "1024-bit Morgan Fingerprints"),
  Split = c("MD", "MD", "MD","MD", "MD", "MD", "MD","MD"),
  RMSE = c(2.475, 2.402, 3.457, 3.524, 3.048, 3.483, 3.537, 3.539),
  PCC = c(0.407, 0.449, 0.183, 0.212, 0.317, 0.310, 0.314, 0.324)
)
plot_df_md <- plot_df_md %>%
  mutate(`Scaled RMSE` = RMSE / max(RMSE))

# Melt the data frame to long format for ggplot2
plot_df_md_long <- plot_df_md %>%
  select(Model, DR, `Scaled RMSE`, PCC) %>%
  pivot_longer(cols = c(`Scaled RMSE`, PCC), names_to = "Metric", values_to = "Value")

plot_df_md_long <- plot_df_md_long %>%
  mutate(DR_Model = paste(Model, DR, sep = " "))

# Combine unique 'DR_Model' values from both datasets
unique_dr_models <- unique(c(unique(plot_df_mp_long$DR_Model), unique(plot_df_md_long$DR_Model)))

# Define a color palette for the combined set
dr_model_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique_dr_models))
names(dr_model_colors) <- unique_dr_models


# Create the bar plot for scaled RMSE and PCC

metrics_mp <- ggplot(plot_df_mp_long, aes(x = Metric, y = Value, fill = DR_Model)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  labs(x = "Metric", y = "Performance", fill = "Method", title = "Performance Metrics") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),  # Adjust size of x-axis text
    axis.text.y = element_text(size = 14),             # Adjust size of y-axis text
    axis.title.x = element_text(size = 16),            # Adjust size of x-axis label
    axis.title.y = element_text(size = 16),            # Adjust size of y-axis label
    plot.title = element_text(size = 18, hjust = 0.5), # Adjust size of title and center it
    legend.title = element_text(size = 16),            # Adjust size of legend title
    legend.text = element_text(size = 14)              # Adjust size of legend text
  ) +
  scale_fill_manual(values = dr_model_colors)

# Print the plot
print(metrics_mp)

# Create the bar plot for scaled RMSE and PCC
metrics_md <- ggplot(plot_df_md_long, aes(x = Metric, y = Value, fill = DR_Model)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  labs(x = "Metric", y = "Performance", fill = "Method", title = "Performance Metrics") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),  # Adjust size of x-axis text
    axis.text.y = element_text(size = 14),             # Adjust size of y-axis text
    axis.title.x = element_text(size = 16),            # Adjust size of x-axis label
    axis.title.y = element_text(size = 16),            # Adjust size of y-axis label
    plot.title = element_text(size = 18, hjust = 0.5), # Adjust size of title and center it
    legend.title = element_text(size = 16),            # Adjust size of legend title
    legend.text = element_text(size = 14)              # Adjust size of legend text
  ) +
  scale_fill_manual(values = dr_model_colors)

# Print the plot
print(metrics_md)

ggsave("MetricsPlot_mp.pdf", metrics_mp, width = 10, height = 5, dpi = 500)
ggsave("MetricsPlot_md.pdf", metrics_md, width = 10, height = 5, dpi = 500)
