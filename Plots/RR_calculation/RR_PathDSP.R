library(data.table)
library(tidyverse)
library(gridExtra)

# best mp
pdsp_1024b_mp <- fread("./PathDSP/results/1024b_mask_comb/output_oGDSC_IF_1024b_mask_comb.FNN.cv_10.Prediction.txt")
pdsp_nd_mp <- fread("./PathDSP/results/1024b_mask_comb/output_oGDSC_IF_1024bnds3_mask_comb.FNN.cv_10.Prediction.txt")
# best mc
pdsp_pbfp_mc <- fread("./PathDSP/results/pbfp_mask_cell/output_oGDSC_IF_pbfp_mask_cell.FNN.cv_10.Prediction.txt")
pdsp_nd_mc <- fread("./PathDSP/results/pbfp_mask_cell/output_oGDSC_IF_pbfpnds3_mask_cell.FNN.cv_10.Prediction.txt")
# best md
pdsp_512b_md <- fread("./PathDSP/results/512b_mask_drug/output_oGDSC_IF_512b_mask_drug.FNN.cv_10.Prediction.txt")
pdsp_nd_md <- fread("./PathDSP/results/512b_mask_drug/output_oGDSC_IF_512bnds1_mask_drug.FNN.cv_10.Prediction.txt")

####### First Step: convert all results into a prediction and response matrix format(if needed)
####### the mat should be columns represent drugs and rows represent cells (312 rows and 144 cols in this case).


####### Second Step: the getRR function, reads in the prediction and response matrix and 
####### output the Recomendation Rate.

####### Third Step: plot the RR over top n best drugs


getRR_PathDSP <- function(df, select_n, top_n) {
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

### create plot_df new
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
    RR1_value <- unlist(getRR_PathDSP(df, 1, i))[1]
    RR2_value <- unlist(getRR_PathDSP(df, i, i))[2]
    
    # Append the result to plot_df
    plot_df <- rbind(plot_df, data.frame(model = df_name, k = i, m = i, RR1 = RR1_value,
                                         RR2 = RR2_value))
  }
  return(plot_df)
}


## create plot df for PathDSP new Mar 30 2024, shuffled null-drug
# add mp
plot_df_new <- create_plot_df_new(pdsp_1024b_mp,10)
plot_df_new <- rbind(plot_df_new, create_plot_df_new(pdsp_nd_mp,10))
# add mc
plot_df_new <- rbind(plot_df_new, create_plot_df_new(pdsp_pbfp_mc,10))
plot_df_new <- rbind(plot_df_new, create_plot_df_new(pdsp_nd_mc,10))
# add md
plot_df_new <- rbind(plot_df_new, create_plot_df_new(pdsp_512b_md,10))
plot_df_new <- rbind(plot_df_new, create_plot_df_new(pdsp_nd_md,10))
## output the plot_df_new
write.csv(plot_df_new, file = "./Plots/RR_Plot/PathDSP_RRplots_new_shuffled.csv")
