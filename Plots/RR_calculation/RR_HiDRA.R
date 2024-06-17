library(data.table)
library(tidyverse)
library(gridExtra)

read_in_dat_Hi <- function(path){
  df <- fread(path)
  df <- df[,-1]
  colnames(df) <- c("drug", "cell", "prediction", "resp", "fold")
  return(df)
}
# best mp
Hi_pbfp_mp <- read_in_dat_Hi("./HiDRA/analysis/compact_results/pbfp_mp.csv")
Hi_256b_mp <- read_in_dat_Hi("./HiDRA/analysis/compact_results/256b_mp.csv")
Hi_1024b_mp <- read_in_dat_Hi("./HiDRA/analysis/compact_results/1024b_mp.csv")
# nd mp
Hi_nd_mp <- read_in_dat_Hi("./HiDRA/analysis/compact_results/pbfpnd_mp.csv")
Hi_pbfp_mc <- read_in_dat_Hi("./HiDRA/analysis/compact_results/pbfp_mc.csv")
Hi_nd_mc <- read_in_dat_Hi("./HiDRA/analysis/compact_results/pbfpnd3_mc.csv")
# best md
Hi_pbfp_md <- read_in_dat_Hi("./HiDRA/analysis/compact_results/pbfp_md.csv")
Hi_512b_md <- read_in_dat_Hi("./HiDRA/analysis/compact_results/512b_md.csv")
Hi_nd_md <- read_in_dat_Hi("./HiDRA/analysis/compact_results/pbfpnd3_md.csv")

##### Calculate RR for HiDRA
####### First Step: convert all results into a prediction and response matrix format(if needed)
####### the mat should be columns represent drugs and rows represent cells (312 rows and 144 cols in this case).


####### Second Step: the getRR function, reads in the prediction and response matrix and 
####### output the Recomendation Rate.

####### Third Step: plot the RR over top n best drugs

getRR_Hi <- function(df, select_n, top_n) {
  # Step 1: Find the drugs with the smallest predicted resp for each cell
  smallest_resp_per_cell <- df %>%
    group_by(cell) %>%
    slice_min(order_by = resp, n = top_n) %>%
    ungroup()
  
  # Step 2: Identify the top n drugs with the smallest IC50 for each cell
  top_predictions_per_cell <- df %>%
    group_by(cell) %>%
    slice_min(order_by = prediction, n = select_n) %>%
    ungroup()
  
  # Step 3: Check if the drug with the smallest resp is among the top 10 drugs with the smallest prediction for each cell
  overlap_check <- smallest_resp_per_cell %>%
    left_join(top_predictions_per_cell, by = "cell", suffix = c("_smallest_resp", "_top_pred"), relationship = "many-to-many") %>%
    mutate(overlap = drug_smallest_resp == drug_top_pred)
  
  # Now, filter or examine overlap_check to see where overlaps occur
  overlap_summary <- overlap_check %>%
    group_by(cell) %>%
    summarise(number_in_top_10 = sum(overlap))
  
  return(
    list(mean(overlap_summary$number_in_top_10 >= 1),
         (mean(overlap_summary$number_in_top_10) / select_n))
  )
}

### create plot_df
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
    RR1_value <- unlist(getRR_Hi(df, 1, i))[1]
    RR2_value <- unlist(getRR_Hi(df, i, i))[2]
    
    # Append the result to plot_df
    plot_df <- rbind(plot_df, data.frame(model = df_name, k = i, m = i, RR1 = RR1_value,
                                         RR2 = RR2_value))
  }
  return(plot_df)
}

## create plot df for HiDRA new
# add mp
plot_df_new <- create_plot_df_new(Hi_pbfp_mp,10)
plot_df_new <- rbind(plot_df_new, create_plot_df_new(Hi_nd_mp,10))
# add mc
plot_df_new <- rbind(plot_df_new, create_plot_df_new(Hi_pbfp_mc,10))
plot_df_new <- rbind(plot_df_new, create_plot_df_new(Hi_nd_mc,10))
# add md
plot_df_new <- rbind(plot_df_new, create_plot_df_new(Hi_pbfp_md,10))
plot_df_new <- rbind(plot_df_new, create_plot_df_new(Hi_nd_md,10))
## output the plot_df_new
write.csv(plot_df_new, file = "./Plots/RR_Plot/HiDRA_RRplots_new_shuffled.csv")
