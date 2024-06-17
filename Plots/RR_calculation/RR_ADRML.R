library(tidyverse)
library(gridExtra)

# the original response matrix
resp_mat <- read.csv("./ADRML/data/resp_ogdsc.csv")

# best and nd, mp
read_in_data <- function(path){
  df <- read.csv(path)
  df <- data.frame(SANGER_MODEL_ID = NA, df)
  colnames(df) <- colnames(resp_mat)
  df$SANGER_MODEL_ID <- resp_mat$SANGER_MODEL_ID
  df <- replace(df, df == 0, NA)
  return(df)
}
# mp
AD_256b_mp <- read_in_data("./ADRML/results/predMat_IF_Ogdsc_256b_mask_comb.csv")
AD_nd_mp <- read_in_data("./ADRML/results/predMat_IF_Ogdsc_nd_mask_comb.csv")
# mc
AD_512b_mc <- read_in_data("./ADRML/results/predMat_IF_Ogdsc_512b_mask_cell.csv")
AD_nd_mc <- read_in_data("./ADRML/results/predMat_IF_Ogdsc_nd_mask_cell.csv")
# md
AD_256b_md <- read_in_data("./ADRML/results/predMat_IF_Ogdsc_256b_mask_drug.csv")
AD_512b_md <- read_in_data("./ADRML/results/predMat_IF_Ogdsc_512b_mask_drug.csv")
AD_1024b_md <- read_in_data("./ADRML/results/predMat_IF_Ogdsc_1024b_mask_drug.csv")
AD_nd_md <- read_in_data("./ADRML/results/predMat_IF_Ogdsc_nd_mask_drug.csv")

## Calculate the RR for ADRML
getRR_sim <- function(real_df, pred_df, select_n, top_n){
  results <- sapply(1:nrow(real_df), function(i) {
    # Extract drug names (assuming the same across both data frames)
    drug_names <- names(real_df)[-1]
    
    # Extract and convert the row to a numeric vector, excluding the first column (SANGER_MODEL_ID)
    real_values_vector <- as.numeric(real_df[i, -1])
    
    # Now sort this numeric vector to get the indices of the top n smallest real values
    top_n_indices <- order(real_values_vector)[1:top_n]
    
    # Get the names of the top 10 drugs with the smallest real values for the current cell line
    top_n_real <- drug_names[top_n_indices]
    
    # Convert pred_df row to a numeric vector to match operations
    pred_values_vector <- as.numeric(pred_df[i, -1])
    
    select_n_indices <- order(pred_values_vector)[1:select_n]
    # Identify the drug with the smallest predicted value
    select_n_pred_drug <- drug_names[select_n_indices]
    
    # Check if the drug with the smallest predicted value is in the list of top 10 real
    select_n_pred_drug %in% top_n_real
  })
  if(is.matrix(results)) {
    RR1 <- mean(apply(results, 2, function(x) any(x == TRUE)))
  } else if(is.vector(results)) {
    RR1 <- mean(results)  # Will be 1 if any TRUE, 0 otherwise
  } else {
    warning("results is not a logical vector or matrix.")
  }
  RR2 <- mean(results)
  return(list(RR1,RR2))
}

#### new plot_df, RR1 fixed select_n = 1 and top_n(m) drugs are varied
#### RR2, select_n(k) = top_n(m)
### create plot_df
create_plot_df_new <- function(real_df, pred_df, n = 10) {
  df_name <- deparse(substitute(pred_df))
  plot_df <- data.frame(model = character(),
                        k = integer(),
                        m = integer(),
                        RR1 = numeric(),
                        RR2 = numeric(),
                        stringsAsFactors = FALSE)
  
  for (i in 1:n) {
    # Call getRR for each i, with pred_df as the model, i as select_n, 
    RR1_value <- unlist(getRR_sim(real_df, pred_df, select_n = 1, top_n = i))[1]
    RR2_value <- unlist(getRR_sim(real_df, pred_df, select_n = i, top_n = i))[2]
    
    # Append the result to plot_df
    plot_df <- rbind(plot_df, data.frame(model = df_name, k = i, m = i, RR1 = RR1_value,
                                         RR2 = RR2_value))
  }
  return(plot_df)
}
# add mp
plot_df_new <- create_plot_df_new(resp_mat,AD_256b_mp,10)
plot_df_new <- rbind(plot_df_new, create_plot_df_new(resp_mat,AD_nd_mp,10))
# add mc
plot_df_new <- rbind(plot_df_new, create_plot_df_new(resp_mat,AD_512b_mc,10))
plot_df_new <- rbind(plot_df_new, create_plot_df_new(resp_mat,AD_nd_mc,10))
# add md
plot_df_new <- rbind(plot_df_new, create_plot_df_new(resp_mat,AD_256b_md,10))
plot_df_new <- rbind(plot_df_new, create_plot_df_new(resp_mat,AD_nd_md,10))
# output
write.csv(plot_df_new, file = "./Plots/RR_Plot/ADRML_RRplots_new.csv")
