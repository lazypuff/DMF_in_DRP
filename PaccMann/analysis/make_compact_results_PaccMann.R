library(tidyverse)

getDF_PaccMann <- function(split_stra, cond){
  combined_df <- data.frame()
  for (i in 1:10) {
    # Generate the filename dynamically
    filename <- sprintf("./PaccMann/paccmann_result/random/gdsc1_%s_%s_%d.csv", split_stra, cond,i)
    
    # Read the current CSV file
    current_df <- read.csv(filename)
    
    # Add a new column to current_df to record the source file index
    
    start_col <- which(names(current_df) == "drug")
    end_col <- which(names(current_df) == "prediction")
    
    subset_current_df <- current_df[, start_col:end_col]
    subset_current_df$fold <- i
    
    combined_df <- rbind(combined_df, subset_current_df)
  }
  return(combined_df)
}
#######

# mask pairs
nd_mp <- getDF_PaccMann("mskcomb", "random")

# mask cell
nd_mc <- getDF_PaccMann("maskcell", "random")

# mask drug
nd_md <- getDF_PaccMann("maskdrug", "random")

write.csv(nd_mp, file = "./PaccMann/analysis/compact_results/nd_mp.csv")
write.csv(nd_mc, file = "./PaccMann/analysis/compact_results/nd_mc.csv")
write.csv(nd_md, file = "./PaccMann/analysis/compact_results/nd_md.csv")


## Mar 26, the results of original and two shuffled
library(tidyverse)

# split_stra = 'maskcell', 'maskcomb' or 'maskdrug'
# cond = 'result' or 'random'
getDF_PaccMann1 <- function(folder, split_stra, cond){
  combined_df <- data.frame()
  for (i in 1:10) {
    # Generate the filename dynamically
    filename <- sprintf("./PaccMann/paccmann_result/random1/%s/gdsc1_%s_%s_%d.csv", folder,split_stra, cond,i)
    
    # Read the current CSV file
    current_df <- read.csv(filename)
    
    # Add a new column to current_df to record the source file index
    
    start_col <- which(names(current_df) == "drug")
    end_col <- which(names(current_df) == "prediction")
    
    subset_current_df <- current_df[, start_col:end_col]
    subset_current_df$fold <- i
    
    combined_df <- rbind(combined_df, subset_current_df)
  }
  return(combined_df)
}
#######

# mask pairs
s_mp <- getDF_PaccMann1("mask_comb","maskcomb", "result")
nd_mp <- getDF_PaccMann1("mask_comb","mskcomb", "random")

# mask cell
s_mc <- getDF_PaccMann1("mask_cell","maskcell", "result")
nd_mc <- getDF_PaccMann1("mask_cell","maskcell", "random")

# mask drug
s_md <- getDF_PaccMann("mask_drug","maskdrug", "result")
nd_md <- getDF_PaccMann("mask_drug","maskdrug", "random")

write.csv(s_mp, file = "./PaccMann/analysis/compact_results/smiles_mp.csv")
write.csv(nd_mp, file = "./PaccMann/analysis/compact_results/nd2_mp.csv")
write.csv(s_mc, file = "./PaccMann/analysis/compact_results/smiles_mc.csv")
write.csv(nd_mc, file = "./PaccMann/analysis/compact_results/nd2_mc.csv")
write.csv(s_md, file = "./PaccMann/analysis/compact_results/smiles_md.csv")
write.csv(nd_md, file = "./PaccMann/analysis/compact_results/nd2_md.csv")


# split_stra = 'maskcell', 'maskcomb' or 'maskdrug'
# cond = 'result' or 'random'
getDF_PaccMann2 <- function(folder, split_stra, cond){
  combined_df <- data.frame()
  for (i in 1:10) {
    # Generate the filename dynamically
    filename <- sprintf("./PaccMann/paccmann_result/random2/%s/gdsc1_%s_%s_%d.csv", folder,split_stra, cond,i)
    
    # Read the current CSV file
    current_df <- read.csv(filename)
    
    # Add a new column to current_df to record the source file index
    
    start_col <- which(names(current_df) == "drug")
    end_col <- which(names(current_df) == "prediction")
    
    subset_current_df <- current_df[, start_col:end_col]
    subset_current_df$fold <- i
    
    combined_df <- rbind(combined_df, subset_current_df)
  }
  return(combined_df)
}
#######

# mask pairs
nd_mp <- getDF_PaccMann2("mask_comb2","mskcomb", "random2")

# mask cell
nd_mc <- getDF_PaccMann2("mask_cell2","maskcell", "random2")

# mask drug
nd_md <- getDF_PaccMann2("mask_drug2","maskdrug", "random2")

write.csv(nd_mp, file = "./PaccMann/analysis/compact_results/nd3_mp.csv")
write.csv(nd_mc, file = "./PaccMann/analysis/compact_results/nd3_mc.csv")
write.csv(nd_md, file = "./PaccMann/analysis/compact_results/nd3_md.csv")
