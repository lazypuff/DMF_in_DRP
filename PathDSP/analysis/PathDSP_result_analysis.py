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
        true_values = group_data['resp']
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

## get the paired Wilcoxon test
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


### OLD GDSC2
og_mask_comb_256b = pd.read_csv("./PathDSP/results/256b_mask_comb/output_oGDSC_IF_256b_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_256b = pd.read_csv("./PathDSP/results/256b_mask_cell/output_oGDSC_IF_256b_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_256b = pd.read_csv("./PathDSP/results/256b_mask_drug/output_oGDSC_IF_256b_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_512b = pd.read_csv("./PathDSP/results/512b_mask_comb/output_oGDSC_IF_512b_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_512b = pd.read_csv("./PathDSP/results/512b_mask_cell/output_oGDSC_IF_512b_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_512b = pd.read_csv("./PathDSP/results/512b_mask_drug/output_oGDSC_IF_512b_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_1024b = pd.read_csv("./PathDSP/results/1024b_mask_comb/output_oGDSC_IF_1024b_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_1024b = pd.read_csv("./PathDSP/results/1024b_mask_cell/output_oGDSC_IF_1024b_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_1024b = pd.read_csv("./PathDSP/results/1024b_mask_drug/output_oGDSC_IF_1024b_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_pbfp = pd.read_csv("./PathDSP/results/pbfp_mask_comb/output_oGDSC_IF_pbfp_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_pbfp = pd.read_csv("./PathDSP/results/pbfp_mask_cell/output_oGDSC_IF_pbfp_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_pbfp = pd.read_csv("./PathDSP/results/pbfp_mask_drug/output_oGDSC_IF_pbfp_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')

og_metric_mask_comb_256b = getperform(og_mask_comb_256b)
og_metric_mask_cell_256b = getperform(og_mask_cell_256b)
og_metric_mask_drug_256b = getperform(og_mask_drug_256b)

og_metric_mask_comb_512b = getperform(og_mask_comb_512b)
og_metric_mask_cell_512b = getperform(og_mask_cell_512b)
og_metric_mask_drug_512b = getperform(og_mask_drug_512b)

og_metric_mask_comb_1024b = getperform(og_mask_comb_1024b)
og_metric_mask_cell_1024b = getperform(og_mask_cell_1024b)
og_metric_mask_drug_1024b = getperform(og_mask_drug_1024b)

og_metric_mask_comb_pbfp = getperform(og_mask_comb_pbfp)
og_metric_mask_cell_pbfp = getperform(og_mask_cell_pbfp)
og_metric_mask_drug_pbfp = getperform(og_mask_drug_pbfp)

getstatistics(og_metric_mask_comb_256b)
getstatistics(og_metric_mask_comb_512b)
getstatistics(og_metric_mask_comb_1024b)
getstatistics(og_metric_mask_comb_pbfp)

getstatistics(og_metric_mask_cell_256b)
getstatistics(og_metric_mask_cell_512b)
getstatistics(og_metric_mask_cell_1024b)
getstatistics(og_metric_mask_cell_pbfp)

getstatistics(og_metric_mask_drug_256b)
getstatistics(og_metric_mask_drug_512b)
getstatistics(og_metric_mask_drug_1024b)
getstatistics(og_metric_mask_drug_pbfp)


### OLD GDSC2 PERMUTATION TYPE, GROUP 1
og_mask_comb_256bnds1 = pd.read_csv("./PathDSP/results/256b_mask_comb/output_oGDSC_IF_256bnds1_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_256bnds1 = pd.read_csv("./PathDSP/results/256b_mask_cell/output_oGDSC_IF_256bnds1_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_256bnds1 = pd.read_csv("./PathDSP/results/256b_mask_drug/output_oGDSC_IF_256bnds1_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_512bnds1 = pd.read_csv("./PathDSP/results/512b_mask_comb/output_oGDSC_IF_512bnds1_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_512bnds1 = pd.read_csv("./PathDSP/results/512b_mask_cell/output_oGDSC_IF_512bnds1_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_512bnds1 = pd.read_csv("./PathDSP/results/512b_mask_drug/output_oGDSC_IF_512bnds1_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_1024bnds1 = pd.read_csv("./PathDSP/results/1024b_mask_comb/output_oGDSC_IF_1024bnds1_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_1024bnds1 = pd.read_csv("./PathDSP/results/1024b_mask_cell/output_oGDSC_IF_1024bnds1_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_1024bnds1 = pd.read_csv("./PathDSP/results/1024b_mask_drug/output_oGDSC_IF_1024bnds1_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_pbfpnds1 = pd.read_csv("./PathDSP/results/pbfp_mask_comb/output_oGDSC_IF_pbfpnds1_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_pbfpnds1 = pd.read_csv("./PathDSP/results/pbfp_mask_cell/output_oGDSC_IF_pbfpnds1_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_pbfpnds1 = pd.read_csv("./PathDSP/results/pbfp_mask_drug/output_oGDSC_IF_pbfpnds1_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')

og_metric_mask_comb_256bnds1 = getperform(og_mask_comb_256bnds1)
og_metric_mask_cell_256bnds1 = getperform(og_mask_cell_256bnds1)
og_metric_mask_drug_256bnds1 = getperform(og_mask_drug_256bnds1)

og_metric_mask_comb_512bnds1 = getperform(og_mask_comb_512bnds1)
og_metric_mask_cell_512bnds1 = getperform(og_mask_cell_512bnds1)
og_metric_mask_drug_512bnds1 = getperform(og_mask_drug_512bnds1)

og_metric_mask_comb_1024bnds1 = getperform(og_mask_comb_1024bnds1)
og_metric_mask_cell_1024bnds1 = getperform(og_mask_cell_1024bnds1)
og_metric_mask_drug_1024bnds1 = getperform(og_mask_drug_1024bnds1)

og_metric_mask_comb_pbfpnds1 = getperform(og_mask_comb_pbfpnds1)
og_metric_mask_cell_pbfpnds1 = getperform(og_mask_cell_pbfpnds1)
og_metric_mask_drug_pbfpnds1 = getperform(og_mask_drug_pbfpnds1)

getstatistics(og_metric_mask_comb_256bnds1)
getstatistics(og_metric_mask_comb_512bnds1)
getstatistics(og_metric_mask_comb_1024bnds1)
getstatistics(og_metric_mask_comb_pbfpnds1)

getstatistics(og_metric_mask_cell_256bnds1)
getstatistics(og_metric_mask_cell_512bnds1)
getstatistics(og_metric_mask_cell_1024bnds1)
getstatistics(og_metric_mask_cell_pbfpnds1)

getstatistics(og_metric_mask_drug_256bnds1)
getstatistics(og_metric_mask_drug_512bnds1)
getstatistics(og_metric_mask_drug_1024bnds1)
getstatistics(og_metric_mask_drug_pbfpnds1)

### OLD GDSC2 PERMUTATION TYPE, GROUP 2
og_mask_comb_256bnds2 = pd.read_csv("./PathDSP/results/256b_mask_comb/output_oGDSC_IF_256bnds2_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_256bnds2 = pd.read_csv("./PathDSP/results/256b_mask_cell/ogdsc_Jun12/output_oGDSC_IF_256bnds2_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_256bnds2 = pd.read_csv("./PathDSP/results/256b_mask_drug/output_oGDSC_IF_256bnds2_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_512bnds2 = pd.read_csv("./PathDSP/results/512b_mask_comb/output_oGDSC_IF_512bnds2_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_512bnds2 = pd.read_csv("./PathDSP/results/512b_mask_cell/output_oGDSC_IF_512bnds2_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_512bnds2 = pd.read_csv("./PathDSP/results/512b_mask_drug/output_oGDSC_IF_512bnds2_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_1024bnds2 = pd.read_csv("./PathDSP/results/1024b_mask_comb/output_oGDSC_IF_1024bnds2_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_1024bnds2 = pd.read_csv("./PathDSP/results/1024b_mask_cell/output_oGDSC_IF_1024bnds2_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_1024bnds2 = pd.read_csv("./PathDSP/results/1024b_mask_drug/output_oGDSC_IF_1024bnds2_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_pbfpnds2 = pd.read_csv("./PathDSP/results/pbfp_mask_comb/output_oGDSC_IF_pbfpnds2_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_pbfpnds2 = pd.read_csv("./PathDSP/results/pbfp_mask_cell/output_oGDSC_IF_pbfpnds2_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_pbfpnds2 = pd.read_csv("./PathDSP/results/pbfp_mask_drug/output_oGDSC_IF_pbfpnds2_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')

og_metric_mask_comb_256bnds2 = getperform(og_mask_comb_256bnds2)
og_metric_mask_cell_256bnds2 = getperform(og_mask_cell_256bnds2)
og_metric_mask_drug_256bnds2 = getperform(og_mask_drug_256bnds2)

og_metric_mask_comb_512bnds2 = getperform(og_mask_comb_512bnds2)
og_metric_mask_cell_512bnds2 = getperform(og_mask_cell_512bnds2)
og_metric_mask_drug_512bnds2 = getperform(og_mask_drug_512bnds2)

og_metric_mask_comb_1024bnds2 = getperform(og_mask_comb_1024bnds2)
og_metric_mask_cell_1024bnds2 = getperform(og_mask_cell_1024bnds2)
og_metric_mask_drug_1024bnds2 = getperform(og_mask_drug_1024bnds2)

og_metric_mask_comb_pbfpnds2 = getperform(og_mask_comb_pbfpnds2)
og_metric_mask_cell_pbfpnds2 = getperform(og_mask_cell_pbfpnds2)
og_metric_mask_drug_pbfpnds2 = getperform(og_mask_drug_pbfpnds2)


getstatistics(og_metric_mask_comb_256bnds2)
getstatistics(og_metric_mask_comb_512bnds2)
getstatistics(og_metric_mask_comb_1024bnds2)
getstatistics(og_metric_mask_comb_pbfpnds2)

getstatistics(og_metric_mask_cell_256bnds2)
getstatistics(og_metric_mask_cell_512bnds2)
getstatistics(og_metric_mask_cell_1024bnds2)
getstatistics(og_metric_mask_cell_pbfpnds2)

getstatistics(og_metric_mask_drug_256bnds2)
getstatistics(og_metric_mask_drug_512bnds2)
getstatistics(og_metric_mask_drug_1024bnds2)
getstatistics(og_metric_mask_drug_pbfpnds2)

### OLD GDSC2 PERMUTATION TYPE, GROUP 3
og_mask_comb_256bnds3 = pd.read_csv("./PathDSP/results/256b_mask_comb/output_oGDSC_IF_256bnds3_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_256bnds3 = pd.read_csv("./PathDSP/results/256b_mask_cell/output_oGDSC_IF_256bnds3_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_256bnds3 = pd.read_csv("./PathDSP/results/256b_mask_drug/output_oGDSC_IF_256bnds3_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_512bnds3 = pd.read_csv("./PathDSP/results/512b_mask_comb/output_oGDSC_IF_512bnds3_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_512bnds3 = pd.read_csv("./PathDSP/results/512b_mask_cell/output_oGDSC_IF_512bnds3_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_512bnds3 = pd.read_csv("./PathDSP/results/512b_mask_drug/output_oGDSC_IF_512bnds3_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_1024bnds3 = pd.read_csv("./PathDSP/results/1024b_mask_comb/output_oGDSC_IF_1024bnds3_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_1024bnds3 = pd.read_csv("./PathDSP/results/1024b_mask_cell/output_oGDSC_IF_1024bnds3_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_1024bnds3 = pd.read_csv("./PathDSP/results/1024b_mask_drug/output_oGDSC_IF_1024bnds3_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_comb_pbfpnds3 = pd.read_csv("./PathDSP/results/pbfp_mask_comb/output_oGDSC_IF_pbfpnds3_mask_comb.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_cell_pbfpnds3 = pd.read_csv("./PathDSP/results/pbfp_mask_cell/output_oGDSC_IF_pbfpnds3_mask_cell.FNN.cv_10.Prediction.txt",
                             sep='\t')
og_mask_drug_pbfpnds3 = pd.read_csv("./PathDSP/results/pbfp_mask_drug/output_oGDSC_IF_pbfpnds3_mask_drug.FNN.cv_10.Prediction.txt",
                             sep='\t')

og_metric_mask_comb_256bnds3 = getperform(og_mask_comb_256bnds3)
og_metric_mask_cell_256bnds3 = getperform(og_mask_cell_256bnds3)
og_metric_mask_drug_256bnds3 = getperform(og_mask_drug_256bnds3)

og_metric_mask_comb_512bnds3 = getperform(og_mask_comb_512bnds3)
og_metric_mask_cell_512bnds3 = getperform(og_mask_cell_512bnds3)
og_metric_mask_drug_512bnds3 = getperform(og_mask_drug_512bnds3)

og_metric_mask_comb_1024bnds3 = getperform(og_mask_comb_1024bnds3)
og_metric_mask_cell_1024bnds3 = getperform(og_mask_cell_1024bnds3)
og_metric_mask_drug_1024bnds3 = getperform(og_mask_drug_1024bnds3)

og_metric_mask_comb_pbfpnds3 = getperform(og_mask_comb_pbfpnds3)
og_metric_mask_cell_pbfpnds3 = getperform(og_mask_cell_pbfpnds3)
og_metric_mask_drug_pbfpnds3 = getperform(og_mask_drug_pbfpnds3)


getstatistics(og_metric_mask_comb_256bnds3)
getstatistics(og_metric_mask_comb_512bnds3)
getstatistics(og_metric_mask_comb_1024bnds3)
getstatistics(og_metric_mask_comb_pbfpnds3)

getstatistics(og_metric_mask_cell_256bnds3)
getstatistics(og_metric_mask_cell_512bnds3)
getstatistics(og_metric_mask_cell_1024bnds3)
getstatistics(og_metric_mask_cell_pbfpnds3)

getstatistics(og_metric_mask_drug_256bnds3)
getstatistics(og_metric_mask_drug_512bnds3)
getstatistics(og_metric_mask_drug_1024bnds3)
getstatistics(og_metric_mask_drug_pbfpnds3)

### Paired Wilcoxon test
#### mask comb
get_pairedWilcoxStat(og_metric_mask_comb_256b,og_metric_mask_comb_256bnds1)
get_pairedWilcoxStat(og_metric_mask_comb_512b,og_metric_mask_comb_512bnds1)
get_pairedWilcoxStat(og_metric_mask_comb_1024b,og_metric_mask_comb_1024bnds1)
get_pairedWilcoxStat(og_metric_mask_comb_pbfp,og_metric_mask_comb_pbfpnds1)
# check which is better
# diff_mask_comb_pbfpnds1 = og_metric_mask_comb_pbfpnds1 - og_metric_mask_comb_pbfp

get_pairedWilcoxStat(og_metric_mask_comb_256b,og_metric_mask_comb_256bnds2)
get_pairedWilcoxStat(og_metric_mask_comb_512b,og_metric_mask_comb_512bnds2)
get_pairedWilcoxStat(og_metric_mask_comb_1024b,og_metric_mask_comb_1024bnds2)
get_pairedWilcoxStat(og_metric_mask_comb_pbfp,og_metric_mask_comb_pbfpnds2)
# check which is better
# diff_mask_comb_pbfpnds2 = og_metric_mask_comb_pbfpnds2 - og_metric_mask_comb_pbfp

get_pairedWilcoxStat(og_metric_mask_comb_256b,og_metric_mask_comb_256bnds3)
get_pairedWilcoxStat(og_metric_mask_comb_512b,og_metric_mask_comb_512bnds3)
get_pairedWilcoxStat(og_metric_mask_comb_1024b,og_metric_mask_comb_1024bnds3)
get_pairedWilcoxStat(og_metric_mask_comb_pbfp,og_metric_mask_comb_pbfpnds3)

#### mask cell
# check which is better
# diff_mask_comb_pbfp = og_metric_mask_cell_pbfp - og_metric_mask_cell_nd

get_pairedWilcoxStat(og_metric_mask_cell_256b,og_metric_mask_cell_256bnds1)
get_pairedWilcoxStat(og_metric_mask_cell_512b,og_metric_mask_cell_512bnds1)
get_pairedWilcoxStat(og_metric_mask_cell_1024b,og_metric_mask_cell_1024bnds1)
get_pairedWilcoxStat(og_metric_mask_cell_pbfp,og_metric_mask_cell_pbfpnds1)

get_pairedWilcoxStat(og_metric_mask_cell_256b,og_metric_mask_cell_256bnds2)
get_pairedWilcoxStat(og_metric_mask_cell_512b,og_metric_mask_cell_512bnds2)
get_pairedWilcoxStat(og_metric_mask_cell_1024b,og_metric_mask_cell_1024bnds2)
get_pairedWilcoxStat(og_metric_mask_cell_pbfp,og_metric_mask_cell_pbfpnds2)

get_pairedWilcoxStat(og_metric_mask_cell_256b,og_metric_mask_cell_256bnds3)
get_pairedWilcoxStat(og_metric_mask_cell_512b,og_metric_mask_cell_512bnds3)
get_pairedWilcoxStat(og_metric_mask_cell_1024b,og_metric_mask_cell_1024bnds3)
get_pairedWilcoxStat(og_metric_mask_cell_pbfp,og_metric_mask_cell_pbfpnds3)


#### mask drug
# check which is better
# diff_mask_comb_pbfp = og_metric_mask_cell_pbfp - og_metric_mask_cell_nd

get_pairedWilcoxStat(og_metric_mask_drug_256b,og_metric_mask_drug_256bnds1)
get_pairedWilcoxStat(og_metric_mask_drug_512b,og_metric_mask_drug_512bnds1)
get_pairedWilcoxStat(og_metric_mask_drug_1024b,og_metric_mask_drug_1024bnds1)
get_pairedWilcoxStat(og_metric_mask_drug_pbfp,og_metric_mask_drug_pbfpnds1)

get_pairedWilcoxStat(og_metric_mask_drug_256b,og_metric_mask_drug_256bnds2)
get_pairedWilcoxStat(og_metric_mask_drug_512b,og_metric_mask_drug_512bnds2)
get_pairedWilcoxStat(og_metric_mask_drug_1024b,og_metric_mask_drug_1024bnds2)
get_pairedWilcoxStat(og_metric_mask_drug_pbfp,og_metric_mask_drug_pbfpnds2)

get_pairedWilcoxStat(og_metric_mask_drug_256b,og_metric_mask_drug_256bnds3)
get_pairedWilcoxStat(og_metric_mask_drug_512b,og_metric_mask_drug_512bnds3)
get_pairedWilcoxStat(og_metric_mask_drug_1024b,og_metric_mask_drug_1024bnds3)
get_pairedWilcoxStat(og_metric_mask_drug_pbfp,og_metric_mask_drug_pbfpnds3)



