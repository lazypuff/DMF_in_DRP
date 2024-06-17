setwd("./Plots/RR_plots")
library(tidyverse)
library(ggthemes)
library(patchwork)
# Read in Data
SR_df_new <- read.csv("SRMF_RRplots_new.csv")
SR_df_new <- SR_df_new[,-1]

AD_df_new <- read.csv("ADRML_RRplots_new.csv")
AD_df_new <- AD_df_new[,-1]

pdsp_df_new <- read.csv("PathDSP_RRplots_new_shuffled.csv")
pdsp_df_new <- pdsp_df_new[,-1]

Hi_df_new <- read.csv("HiDRA_RRplots_new_shuffled.csv")
Hi_df_new <- Hi_df_new[,-1]

PM_df_new <- read.csv("PM_RRplots_new_shuffled.csv")
PM_df_new <- PM_df_new[,-1]

plot_df_new <- rbind(SR_df_new,AD_df_new,pdsp_df_new,Hi_df_new, PM_df_new)
# data processing code
df_process_onlyWD <- function(df,split_set) {
  plot_df_processed <- df %>%
    separate(model, into = c("model_name", "cond", "split"), sep = "_") %>%
    # Filter for 'split' == 'md', 'cond' != 'nd', and 'm' == 10
    filter(split == split_set) %>%
    rename(Setting = cond) %>%
    mutate(model_name = recode(model_name,
                               'AD' = 'ADRML',
                               'Hi' = 'HiDRA',
                               'pdsp' = 'PathDSP',
                               'PM' = 'PaccMann',
                               'SR' = 'SRMF'),
           Setting = recode(Setting,
                            'pbfp' = 'PubChem Fingerprints',
                            '256b' = '256-bit Morgan Fingerprints',
                            '512b' = '512-bit Morgan Fingerprints',
                            '1024b' = '1024-bit Morgan Fingerprints',
                            's' = 'SMILES')
    ) %>% 
    mutate(model_setting = paste(model_name, Setting, sep = " "))

  return(plot_df_processed)
}

### create the only WD for the best method in each split
# mask-pairs
p_1 <- df_process_onlyWD(plot_df_new, "mp") %>% 
  filter(k <= 5, model_setting == 'HiDRA PubChem Fingerprints') %>%
  ggplot(aes(x = k, y = RR1)) + 
  geom_line(size = 1.5, alpha=0.8, color = "#449B75") +
  labs(x = "The Number of Top Performing Drugs", y = expression(RR[1])) + 
  scale_x_continuous(breaks = 1:10) +  # Set x-axis breaks
  theme_fivethirtyeight() + 
  theme(axis.title = element_text(), legend.position = "none") + 
  scale_color_brewer(palette = "Set1") + 
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),  # Adjust size of x-axis text
    axis.text.y = element_text(size = 14),             # Adjust size of y-axis text
    axis.title.x = element_text(size = 16),            # Adjust size of x-axis label
    axis.title.y = element_text(size = 16),            # Adjust size of y-axis label
    plot.title = element_text(size = 18, hjust = 0.5), # Adjust size of title and center it
    legend.title = element_text(size = 16),            # Adjust size of legend title
    legend.text = element_text(size = 14)              # Adjust size of legend text
  )

p_2 <- df_process_onlyWD(plot_df_new, "mp") %>% 
  filter(model_setting == 'HiDRA PubChem Fingerprints') %>%
  ggplot(aes(x = k, y = RR2)) + 
  geom_line(size = 1.5, alpha=0.8, color = "#449B75") +
  labs(x = "The Number of Top Performing Drugs", y = expression(RR[2])) + 
  scale_x_continuous(breaks = 1:10) +  # Set x-axis breaks
  theme_fivethirtyeight() + 
  theme(axis.title = element_text(), legend.position = "none") + 
  scale_color_brewer(palette = "Set1") + 
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),  # Adjust size of x-axis text
    axis.text.y = element_text(size = 14),             # Adjust size of y-axis text
    axis.title.x = element_text(size = 16),            # Adjust size of x-axis label
    axis.title.y = element_text(size = 16),            # Adjust size of y-axis label
    plot.title = element_text(size = 18, hjust = 0.5), # Adjust size of title and center it
    legend.title = element_text(size = 16),            # Adjust size of legend title
    legend.text = element_text(size = 14)              # Adjust size of legend text
  )

combined_plot_mp <- (p_1 + p_2 + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 20),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"))  # Adjust the size here
# mask-cells
p_1 <- df_process_onlyWD(plot_df_new, "mc") %>% 
  filter(k <= 5, model_setting == 'PaccMann SMILES') %>%
  ggplot(aes(x = k, y = RR1)) + 
  geom_line(size = 1.5, alpha=0.8, color = "#6B886D") +
  labs(x = "The Number of Top Performing Drugs", y = expression(RR[1])) + 
  scale_x_continuous(breaks = 1:10) +  # Set x-axis breaks
  theme_fivethirtyeight() + 
  theme(axis.title = element_text(), legend.position = "none") + 
  scale_color_brewer(palette = "Set1") + 
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),  # Adjust size of x-axis text
    axis.text.y = element_text(size = 14),             # Adjust size of y-axis text
    axis.title.x = element_text(size = 16),            # Adjust size of x-axis label
    axis.title.y = element_text(size = 16),            # Adjust size of y-axis label
    plot.title = element_text(size = 18, hjust = 0.5), # Adjust size of title and center it
    legend.title = element_text(size = 16),            # Adjust size of legend title
    legend.text = element_text(size = 14)              # Adjust size of legend text
  )

p_2 <- df_process_onlyWD(plot_df_new, "mc") %>% 
  filter(model_setting == 'PaccMann SMILES') %>%
  ggplot(aes(x = k, y = RR2)) + 
  geom_line(size = 1.5, alpha=0.8, color = "#6B886D") +
  labs(x = "The Number of Top Performing Drugs", y = expression(RR[2])) + 
  scale_x_continuous(breaks = 1:10) +  # Set x-axis breaks
  theme_fivethirtyeight() + 
  theme(axis.title = element_text(), legend.position = "none") + 
  scale_color_brewer(palette = "Set1") + 
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),  # Adjust size of x-axis text
    axis.text.y = element_text(size = 14),             # Adjust size of y-axis text
    axis.title.x = element_text(size = 16),            # Adjust size of x-axis label
    axis.title.y = element_text(size = 16),            # Adjust size of y-axis label
    plot.title = element_text(size = 18, hjust = 0.5), # Adjust size of title and center it
    legend.title = element_text(size = 16),            # Adjust size of legend title
    legend.text = element_text(size = 14)              # Adjust size of legend text
  )

combined_plot_mc <- (p_1 + p_2 + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 20),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"))  # Adjust the size here
# mask-drugs
p_1 <- df_process_onlyWD(plot_df_new, "md") %>% 
  filter(k <= 5, model_setting == 'HiDRA PubChem Fingerprints') %>%
  ggplot(aes(x = k, y = RR1)) + 
  geom_line(size = 1.5, alpha=0.8, color = "#449B75") +
  labs(x = "The Number of Top Performing Drugs", y = expression(RR[1])) + 
  scale_x_continuous(breaks = 1:10) +  # Set x-axis breaks
  theme_fivethirtyeight() + 
  theme(axis.title = element_text(), legend.position = "none") + 
  scale_color_brewer(palette = "Set1") + 
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),  # Adjust size of x-axis text
    axis.text.y = element_text(size = 14),             # Adjust size of y-axis text
    axis.title.x = element_text(size = 16),            # Adjust size of x-axis label
    axis.title.y = element_text(size = 16),            # Adjust size of y-axis label
    plot.title = element_text(size = 18, hjust = 0.5), # Adjust size of title and center it
    legend.title = element_text(size = 16),            # Adjust size of legend title
    legend.text = element_text(size = 14)              # Adjust size of legend text
  )

p_2 <- df_process_onlyWD(plot_df_new, "md") %>% 
  filter(model_setting == 'HiDRA PubChem Fingerprints') %>%
  ggplot(aes(x = k, y = RR2)) + 
  geom_line(size = 1.5, alpha=0.8, color = "#449B75") +
  labs(x = "The Number of Top Performing Drugs", y = expression(RR[2])) + 
  scale_x_continuous(breaks = 1:10) +  # Set x-axis breaks
  theme_fivethirtyeight() + 
  theme(axis.title = element_text(), legend.position = "none") + 
  scale_color_brewer(palette = "Set1") + 
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),  # Adjust size of x-axis text
    axis.text.y = element_text(size = 14),             # Adjust size of y-axis text
    axis.title.x = element_text(size = 16),            # Adjust size of x-axis label
    axis.title.y = element_text(size = 16),            # Adjust size of y-axis label
    plot.title = element_text(size = 18, hjust = 0.5), # Adjust size of title and center it
    legend.title = element_text(size = 16),            # Adjust size of legend title
    legend.text = element_text(size = 14)              # Adjust size of legend text
  )

combined_plot_md <- (p_1 + p_2 + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 20),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"))  # Adjust the size here

ggsave("combined_bestRRplot_mp.pdf", combined_plot_mp, width = 10, height = 5, dpi = 500)
ggsave("combined_bestRRplot_mc.pdf", combined_plot_mc, width = 10, height = 5, dpi = 500)
ggsave("combined_bestRRplot_md.pdf", combined_plot_md, width = 10, height = 5, dpi = 500)