import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr
from scipy.stats import wilcoxon

### old gdsc2 dataset results
# mask comb
og_256b_mask_comb = pd.read_csv('./SRMF/results/results_oldgdsc_256b_mask_comb.csv',header=None)
og_512b_mask_comb = pd.read_csv('./SRMF/results/results_oldgdsc_512b_mask_comb.csv',header=None)
og_1024b_mask_comb = pd.read_csv('./SRMF/results/results_oldgdsc_1024b_mask_comb.csv',header=None)
og_pbfp_mask_comb = pd.read_csv('./SRMF/results/results_oldgdsc_pbfp_mask_comb.csv',header=None)
og_nd_mask_comb = pd.read_csv('./SRMF/results/results_oldgdsc_nd_mask_comb.csv',header=None)
# mask cell
og_256b_mask_cell = pd.read_csv('./SRMF/results/results_oldgdsc_256b_mask_cell.csv',header=None)
og_512b_mask_cell = pd.read_csv('./SRMF/results/results_oldgdsc_512b_mask_cell.csv',header=None)
og_1024b_mask_cell = pd.read_csv('./SRMF/results/results_oldgdsc_1024b_mask_cell.csv',header=None)
og_pbfp_mask_cell = pd.read_csv('./SRMF/results/results_oldgdsc_pbfp_mask_cell.csv',header=None)
og_nd_mask_cell = pd.read_csv('./SRMF/results/results_oldgdsc_nd_mask_cell.csv',header=None)
# mask drug
og_256b_mask_drug = pd.read_csv('./SRMF/results/results_oldgdsc_256b_mask_drug.csv',header=None)
og_512b_mask_drug = pd.read_csv('./SRMF/results/results_oldgdsc_512b_mask_drug.csv',header=None)
og_1024b_mask_drug = pd.read_csv('./SRMF/results/results_oldgdsc_1024b_mask_drug.csv',header=None)
og_pbfp_mask_drug = pd.read_csv('./SRMF/results/results_oldgdsc_pbfp_mask_drug.csv',header=None)
og_nd_mask_drug = pd.read_csv('./SRMF/results/results_oldgdsc_nd_mask_drug.csv',header=None)



def getstatistics(metrics_data):
    mean_rmse = metrics_data[4].mean()
    std_rmse = metrics_data[4].std()

    mean_pcc = metrics_data[5].mean()
    std_pcc = metrics_data[5].std()

    mean_r2 = metrics_data[6].mean()
    std_r2 = metrics_data[6].std()

    median_rmse = np.median(metrics_data[4])
    q1_rmse = np.percentile(metrics_data[4], 25)
    q3_rmse = np.percentile(metrics_data[4], 75)
    iqr_rmse = q3_rmse - q1_rmse

    # Calculate median and IQR for Pearson Correlation Coefficient
    median_pcc = np.median(metrics_data[5])
    q1_pcc = np.percentile(metrics_data[5], 25)
    q3_pcc = np.percentile(metrics_data[5], 75)
    iqr_pcc = q3_pcc - q1_pcc

    # Calculate median and IQR for R2
    median_r2 = np.median(metrics_data[6])
    q1_r2 = np.percentile(metrics_data[6], 25)
    q3_r2 = np.percentile(metrics_data[6], 75)
    iqr_r2 = q3_r2 - q1_r2

    # Print the calculated means and standard deviations
    print(f"Mean RMSE: {mean_rmse:.4f} (Std: {std_rmse:.4f})",f"Median: {median_rmse:.4f}",f"IQR: {iqr_rmse:.4f}")
    print(f"Mean Pearson Correlation: {mean_pcc:.4f} (Std: {std_pcc:.4f})",f"Median: {median_pcc:.4f}",f"IQR: {iqr_pcc:.4f}")
    print(f"Mean R2: {mean_r2:.4f} (Std: {std_r2:.4f})",f"Median: {median_r2:.4f}",f"IQR: {iqr_r2:.4f}")

def get_pairedWilcoxStat(metrics1, metrics2):
    metric_diff = metrics1 - metrics2
    rmse_stat, rmse_p_value = wilcoxon(metric_diff[4])

    # Perform paired Wilcoxon signed-rank test for Pearson Correlation Coefficient
    pcc_stat, pcc_p_value = wilcoxon(metric_diff[5])

    # Perform paired Wilcoxon signed-rank test for R2
    r2_stat, r2_p_value = wilcoxon(metric_diff[6])

    print("Paired Wilcoxon signed-rank test results:")
    print(f"RMSE - Statistic: {rmse_stat:.4f}, p-value: {rmse_p_value:.4f}")
    print(f"Pearson Correlation - Statistic: {pcc_stat:.4f}, p-value: {pcc_p_value:.4f}")
    print(f"R2 - Statistic: {r2_stat:.4f}, p-value: {r2_p_value:.4f}")

### old gdsc2 811 trained performance statistics
# mask comb
getstatistics(og_256b_mask_comb)
getstatistics(og_512b_mask_comb)
getstatistics(og_1024b_mask_comb)
getstatistics(og_pbfp_mask_comb)
getstatistics(og_nd_mask_comb)
# mask cell
getstatistics(og_256b_mask_cell)
getstatistics(og_512b_mask_cell)
getstatistics(og_1024b_mask_cell)
getstatistics(og_pbfp_mask_cell)
getstatistics(og_nd_mask_cell)
# mask drug
getstatistics(og_256b_mask_drug)
getstatistics(og_512b_mask_drug)
getstatistics(og_1024b_mask_drug)
getstatistics(og_pbfp_mask_drug)
getstatistics(og_nd_mask_drug)



### now check the drug effect
get_pairedWilcoxStat(og_256b_mask_comb,og_nd_mask_comb)
get_pairedWilcoxStat(og_512b_mask_comb,og_nd_mask_comb)
get_pairedWilcoxStat(og_1024b_mask_comb,og_nd_mask_comb)
get_pairedWilcoxStat(og_pbfp_mask_comb,og_nd_mask_comb)


get_pairedWilcoxStat(og_256b_mask_cell,og_nd_mask_cell)
get_pairedWilcoxStat(og_512b_mask_cell,og_nd_mask_cell)
get_pairedWilcoxStat(og_1024b_mask_cell,og_nd_mask_cell)
get_pairedWilcoxStat(og_pbfp_mask_cell,og_nd_mask_cell)


get_pairedWilcoxStat(og_256b_mask_drug,og_nd_mask_drug)
get_pairedWilcoxStat(og_512b_mask_drug,og_nd_mask_drug)
get_pairedWilcoxStat(og_1024b_mask_drug,og_nd_mask_drug)
get_pairedWilcoxStat(og_pbfp_mask_drug,og_nd_mask_drug)

