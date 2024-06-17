library(tidyverse)
library(gridExtra)
# mask pair
PM_s_mp <- read.csv("./PaccMann/analysis/compact_results/smiles_mp.csv")
PM_s_mp <- PM_s_mp[,-1]
PM_nd_mp <- read.csv("./PaccMann/analysis/compact_results/nd2_mp.csv")
PM_nd_mp <- PM_nd_mp[,-1]

# mask cell
PM_s_mc <- read.csv("./PaccMann/analysis/compact_results/smiles_mc.csv")
PM_s_mc <- PM_s_mc[,-1]
PM_nd_mc <- read.csv("./PaccMann/analysis/compact_results/nd2_mc.csv")
PM_nd_mc <- PM_nd_mc[,-1]

# mask drug
PM_s_md <- read.csv("./PaccMann/analysis/compact_results/smiles_md.csv")
PM_s_md <- PM_s_md[,-1]
PM_nd_md <- read.csv("./PaccMann/analysis/compact_results/nd2_md.csv")
PM_nd_md <- PM_nd_md[,-1]

# Calculate the RR for PaccMann Results
getRR_PM <- function(df, select_n, top_n) {
  # Step 1: Find the drugs with the smallest predicted resp for each cell
  smallest_resp_per_cell <- df %>%
    group_by(cell_line) %>%
    slice_min(order_by = IC50, n = top_n) %>%
    ungroup()
  
  # Step 2: Identify the top n drugs with the smallest IC50 for each cell
  top_predictions_per_cell <- df %>%
    group_by(cell_line) %>%
    slice_min(order_by = prediction, n = select_n) %>%
    ungroup()
  
  # Step 3: Check if the drug with the smallest resp is among the top 10 drugs with the smallest prediction for each cell
  overlap_check <- smallest_resp_per_cell %>%
    left_join(top_predictions_per_cell, by = "cell_line", suffix = c("_smallest_resp", "_top_pred"), relationship = "many-to-many") %>%
    mutate(overlap = drug_smallest_resp == drug_top_pred)
  
  # Now, filter or examine overlap_check to see where overlaps occur
  overlap_summary <- overlap_check %>%
    group_by(cell_line) %>%
    summarise(number_in_top_10 = sum(overlap))
  
  return(
    list(mean(overlap_summary$number_in_top_10 >= 1),
         (mean(overlap_summary$number_in_top_10) / select_n))
  )
}


##### create plot df for PaccMann new
create_plot_df_new <- function(df, top_n = 10) {
  df_name <- deparse(substitute(df))
  plot_df <- data.frame(model = character(),
                        k = integer(),
                        m = integer(),
                        RR1 = numeric(),
                        RR2 = numeric(),
                        stringsAsFactors = FALSE)
  
  for (i in 1:top_n) {
    # Call getRR for each i, with df as the model, i as select_n, and 10 as top_n
    RR1_value <- unlist(getRR_PM(df, 1, i))[1]
    RR2_value <- unlist(getRR_PM(df, i, i))[2]
    
    # Append the result to plot_df
    plot_df <- rbind(plot_df, data.frame(model = df_name, k = i, m = i, RR1 = RR1_value,
                                         RR2 = RR2_value))
  }
  return(plot_df)
}

# add mp
plot_df_new <- create_plot_df_new(PM_s_mp,10)
plot_df_new <- rbind(plot_df_new, create_plot_df_new(PM_nd_mp,10))
# add mc
plot_df_new <- rbind(plot_df_new, create_plot_df_new(PM_s_mc,10))
plot_df_new <- rbind(plot_df_new, create_plot_df_new(PM_nd_mc,10))
# add md
plot_df_new <- rbind(plot_df_new, create_plot_df_new(PM_s_md,10))
plot_df_new <- rbind(plot_df_new, create_plot_df_new(PM_nd_md,10))

# output PaccMann plot df
write.csv(plot_df_new, file = "./Plots/RR_Plot/PM_RRplots_new_shuffled.csv")
