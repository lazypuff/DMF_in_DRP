import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr
from scipy.stats import wilcoxon


og_256b_mask_comb = pd.read_csv('./ADRML/results/IF_Ogdsc_256b_mask_comb.csv')
og_256b_mask_cell = pd.read_csv('./ADRML/results/IF_Ogdsc_256b_mask_cell.csv')
og_256b_mask_drug = pd.read_csv('./ADRML/results/IF_Ogdsc_256b_mask_drug.csv')

og_512b_mask_comb = pd.read_csv('./ADRML/results/IF_Ogdsc_512b_mask_comb.csv')
og_512b_mask_cell = pd.read_csv('./ADRML/results/IF_Ogdsc_512b_mask_cell.csv')
og_512b_mask_drug = pd.read_csv('./ADRML/results/IF_Ogdsc_512b_mask_drug.csv')

og_1024b_mask_comb = pd.read_csv('./ADRML/results/IF_Ogdsc_1024b_mask_comb.csv')
og_1024b_mask_cell = pd.read_csv('./ADRML/results/IF_Ogdsc_1024b_mask_cell.csv')
og_1024b_mask_drug = pd.read_csv('./ADRML/results/IF_Ogdsc_1024b_mask_drug.csv')

og_2048b_mask_comb = pd.read_csv('./ADRML/results/IF_Ogdsc_2048b_mask_comb.csv')
og_2048b_mask_cell = pd.read_csv('./ADRML/results/IF_Ogdsc_2048b_mask_cell.csv')
og_2048b_mask_drug = pd.read_csv('./ADRML/results/IF_Ogdsc_2048b_mask_drug.csv')

og_nd_mask_comb = pd.read_csv('./ADRML/results/IF_Ogdsc_nd_mask_comb.csv')
og_nd_mask_cell = pd.read_csv('./ADRML/results/IF_Ogdsc_nd_mask_cell.csv')
og_nd_mask_drug = pd.read_csv('./ADRML/results/IF_Ogdsc_nd_mask_drug.csv')


def getstatistics(metrics_data):
    mean_rmse = metrics_data['rmse'].mean()
    std_rmse = metrics_data['rmse'].std()

    mean_pcc = metrics_data['PCC'].mean()
    std_pcc = metrics_data['PCC'].std()

    mean_r2 = metrics_data['R2'].mean()
    std_r2 = metrics_data['R2'].std()

    median_rmse = np.median(metrics_data['rmse'])
    q1_rmse = np.percentile(metrics_data['rmse'], 25)
    q3_rmse = np.percentile(metrics_data['rmse'], 75)
    iqr_rmse = q3_rmse - q1_rmse

    # Calculate median and IQR for Pearson Correlation Coefficient
    median_pcc = np.median(metrics_data['PCC'])
    q1_pcc = np.percentile(metrics_data['PCC'], 25)
    q3_pcc = np.percentile(metrics_data['PCC'], 75)
    iqr_pcc = q3_pcc - q1_pcc

    # Calculate median and IQR for R2
    median_r2 = np.median(metrics_data['R2'])
    q1_r2 = np.percentile(metrics_data['R2'], 25)
    q3_r2 = np.percentile(metrics_data['R2'], 75)
    iqr_r2 = q3_r2 - q1_r2

    # Print the calculated means and standard deviations
    print(f"Mean RMSE: {mean_rmse:.4f} (Std: {std_rmse:.4f})",f"Median: {median_rmse:.4f}",f"IQR: {iqr_rmse:.4f}")
    print(f"Mean Pearson Correlation: {mean_pcc:.4f} (Std: {std_pcc:.4f})",f"Median: {median_pcc:.4f}",f"IQR: {iqr_pcc:.4f}")
    print(f"Mean R2: {mean_r2:.4f} (Std: {std_r2:.4f})",f"Median: {median_r2:.4f}",f"IQR: {iqr_r2:.4f}")

def get_pairedWilcoxStat(metrics1, metrics2):
    metric_diff = metrics1 - metrics2
    rmse_stat, rmse_p_value = wilcoxon(metric_diff['rmse'])

    # Perform paired Wilcoxon signed-rank test for Pearson Correlation Coefficient
    pcc_stat, pcc_p_value = wilcoxon(metric_diff['PCC'])

    # Perform paired Wilcoxon signed-rank test for R2
    r2_stat, r2_p_value = wilcoxon(metric_diff['R2'])

    print("Paired Wilcoxon signed-rank test results:")
    print(f"RMSE - Statistic: {rmse_stat:.4f}, p-value: {rmse_p_value:.4f}")
    print(f"Pearson Correlation - Statistic: {pcc_stat:.4f}, p-value: {pcc_p_value:.4f}")
    print(f"R2 - Statistic: {r2_stat:.4f}, p-value: {r2_p_value:.4f}")

getstatistics(og_256b_mask_comb)
getstatistics(og_256b_mask_cell)
getstatistics(og_256b_mask_drug)

getstatistics(og_512b_mask_comb)
getstatistics(og_512b_mask_cell)
getstatistics(og_512b_mask_drug)

getstatistics(og_1024b_mask_comb)
getstatistics(og_1024b_mask_cell)
getstatistics(og_1024b_mask_drug)

getstatistics(og_nd_mask_comb)
getstatistics(og_nd_mask_cell)
getstatistics(og_nd_mask_drug)

getstatistics(og_2048b_mask_comb)
getstatistics(og_2048b_mask_cell)
getstatistics(og_2048b_mask_drug)

get_pairedWilcoxStat(og_256b_mask_comb,og_nd_mask_comb)
get_pairedWilcoxStat(og_512b_mask_comb,og_nd_mask_comb)
get_pairedWilcoxStat(og_1024b_mask_comb,og_nd_mask_comb)
get_pairedWilcoxStat(og_2048b_mask_comb,og_nd_mask_comb)

get_pairedWilcoxStat(og_256b_mask_cell,og_nd_mask_cell)
get_pairedWilcoxStat(og_512b_mask_cell,og_nd_mask_cell)
get_pairedWilcoxStat(og_1024b_mask_cell,og_nd_mask_cell)
get_pairedWilcoxStat(og_2048b_mask_cell,og_nd_mask_cell)

get_pairedWilcoxStat(og_512b_mask_cell,og_1024b_mask_cell)

get_pairedWilcoxStat(og_256b_mask_drug,og_nd_mask_drug)
get_pairedWilcoxStat(og_512b_mask_drug,og_nd_mask_drug)
get_pairedWilcoxStat(og_1024b_mask_drug,og_nd_mask_drug)
get_pairedWilcoxStat(og_2048b_mask_drug,og_nd_mask_drug)

get_pairedWilcoxStat(og_256b_mask_drug,og_512b_mask_drug)
