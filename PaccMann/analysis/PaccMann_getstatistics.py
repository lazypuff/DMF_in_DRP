import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr
from scipy.stats import wilcoxon


def getperform(dataset):
    rmse_list = []
    pearson_corr_list = []
    r2_list = []

    # Group the DataFrame by 'fold' column
    grouped = dataset.groupby('fold')

    # Iterate over each group (fold)
    for fold, group_data in grouped:
        true_values = group_data['IC50']
        predicted_values = group_data['prediction']

        # Calculate RMSE
        rmse = mean_squared_error(true_values, predicted_values, squared=False)
        rmse_list.append(rmse)

        # Calculate Pearson Correlation Coefficient
        pearson_corr, _ = pearsonr(true_values, predicted_values)
        pearson_corr_list.append(pearson_corr)

        # Calculate R2 Score
        r2 = r2_score(true_values, predicted_values)
        r2_list.append(r2)

    # Create a DataFrame to store the metrics for each fold
    metrics_df = pd.DataFrame({
        'fold': grouped.groups.keys(),
        'RMSE': rmse_list,
        'Pearson_Corr': pearson_corr_list,
        'R2': r2_list
    })
    return metrics_df


def getstatistics(metrics_data):
    mean_rmse = metrics_data['RMSE'].mean()
    std_rmse = metrics_data['RMSE'].std()

    mean_pcc = metrics_data['Pearson_Corr'].mean()
    std_pcc = metrics_data['Pearson_Corr'].std()

    mean_r2 = metrics_data['R2'].mean()
    std_r2 = metrics_data['R2'].std()

    median_rmse = np.median(metrics_data['RMSE'])
    q1_rmse = np.percentile(metrics_data['RMSE'], 25)
    q3_rmse = np.percentile(metrics_data['RMSE'], 75)
    iqr_rmse = q3_rmse - q1_rmse

    # Calculate median and IQR for Pearson Correlation Coefficient
    median_pcc = np.median(metrics_data['Pearson_Corr'])
    q1_pcc = np.percentile(metrics_data['Pearson_Corr'], 25)
    q3_pcc = np.percentile(metrics_data['Pearson_Corr'], 75)
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
    rmse_stat, rmse_p_value = wilcoxon(metric_diff['RMSE'])

    # Perform paired Wilcoxon signed-rank test for Pearson Correlation Coefficient
    pcc_stat, pcc_p_value = wilcoxon(metric_diff['Pearson_Corr'])

    # Perform paired Wilcoxon signed-rank test for R2
    r2_stat, r2_p_value = wilcoxon(metric_diff['R2'])

    print("Paired Wilcoxon signed-rank test results:")
    print(f"RMSE - Statistic: {rmse_stat:.4f}, p-value: {rmse_p_value:.4f}")
    print(f"Pearson Correlation - Statistic: {pcc_stat:.4f}, p-value: {pcc_p_value:.4f}")
    print(f"R2 - Statistic: {r2_stat:.4f}, p-value: {r2_p_value:.4f}")


def read_data_paccmann(path):
    df = pd.read_csv(path)
    df = df.iloc[:,1:]
    return df

## authentic smiles
og_mask_comb_s = read_data_paccmann("./PaccMann/analysis/compact_results/smiles_mp.csv")
og_mask_cell_s = read_data_paccmann("./PaccMann/analysis/compact_results/smiles_mc.csv")
og_mask_drug_s = read_data_paccmann("./PaccMann/analysis/compact_results/smiles_md.csv")
## shuffled smiles 1
og_mask_comb_nd1 = read_data_paccmann("./PaccMann/analysis/compact_results/nd_mp.csv")
og_mask_cell_nd1 = read_data_paccmann("./PaccMann/analysis/compact_results/nd_mc.csv")
og_mask_drug_nd1 = read_data_paccmann("./PaccMann/analysis/compact_results/nd_md.csv")
## shuffled smiles 2
og_mask_comb_nd2 = read_data_paccmann("./PaccMann/analysis/compact_results/nd2_mp.csv")
og_mask_cell_nd2 = read_data_paccmann("./PaccMann/analysis/compact_results/nd2_mc.csv")
og_mask_drug_nd2 = read_data_paccmann("./PaccMann/analysis/compact_results/nd2_md.csv")
## shuffled smiles 3
og_mask_comb_nd3 = read_data_paccmann("./PaccMann/analysis/compact_results/nd3_mp.csv")
og_mask_cell_nd3 = read_data_paccmann("./PaccMann/analysis/compact_results/nd3_mc.csv")
og_mask_drug_nd3 = read_data_paccmann("./PaccMann/analysis/compact_results/nd3_md.csv")


metrics_mp_s = getperform(og_mask_comb_s)
metrics_mc_s = getperform(og_mask_cell_s)
metrics_md_s = getperform(og_mask_drug_s)

metrics_mp_nd1 = getperform(og_mask_comb_nd1)
metrics_mc_nd1 = getperform(og_mask_cell_nd1)
metrics_md_nd1 = getperform(og_mask_drug_nd1)

metrics_mp_nd2 = getperform(og_mask_comb_nd2)
metrics_mc_nd2 = getperform(og_mask_cell_nd2)
metrics_md_nd2 = getperform(og_mask_drug_nd2)

metrics_mp_nd3 = getperform(og_mask_comb_nd3)
metrics_mc_nd3 = getperform(og_mask_cell_nd3)
metrics_md_nd3 = getperform(og_mask_drug_nd3)

## get metrics statsitics
# authentic smiles
getstatistics(metrics_mp_s)
getstatistics(metrics_mc_s)
getstatistics(metrics_md_s)
# shuffled smiles 1
getstatistics(metrics_mp_nd1)
getstatistics(metrics_mc_nd1)
getstatistics(metrics_md_nd1)
# shuffled smiles 2
getstatistics(metrics_mp_nd2)
getstatistics(metrics_mc_nd2)
getstatistics(metrics_md_nd2)
# shuffled smiles 3
getstatistics(metrics_mp_nd3)
getstatistics(metrics_mc_nd3)
getstatistics(metrics_md_nd3)

### paired wilcoxon rank test
# with nd1
get_pairedWilcoxStat(metrics_mp_s, metrics_mp_nd1)
get_pairedWilcoxStat(metrics_mc_s, metrics_mc_nd1)
get_pairedWilcoxStat(metrics_md_s, metrics_md_nd1)

# with nd2
get_pairedWilcoxStat(metrics_mp_s, metrics_mp_nd2)
get_pairedWilcoxStat(metrics_mc_s, metrics_mc_nd2)
get_pairedWilcoxStat(metrics_md_s, metrics_md_nd2)

# with nd3
get_pairedWilcoxStat(metrics_mp_s, metrics_mp_nd3)
get_pairedWilcoxStat(metrics_mc_s, metrics_mc_nd3)
get_pairedWilcoxStat(metrics_md_s, metrics_md_nd3)
