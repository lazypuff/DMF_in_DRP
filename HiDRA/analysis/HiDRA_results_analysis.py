import os
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr
from scipy.stats import wilcoxon

def getdf(folder_path_prediction, folder_path_true,which_dataset,suffix_exp_set,split):
    if split == 'mc':
        suffix_split = "_mc"
    elif split == 'md':
        suffix_split = "_md"
    else:
        suffix_split = ""

    folder_path = folder_path_prediction  # Replace with the actual path to your folder

    # Loop through the range of numbers from 1 to 10
    data_frames = []

    # Loop through the range of numbers from 1 to 10
    for i in range(1, 11):
        csv_file = f'{which_dataset}_{suffix_exp_set}_cv{i}_pred{suffix_split}.csv'
        file_path = os.path.join(folder_path, csv_file)

        # Check if the file exists before attempting to read
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            df['fold'] = i  # Add a new 'fold' column with the current fold number
            data_frames.append(df)
        else:
            print(f"File {csv_file} not found.")

    # Concatenate all data frames into one
    big_data_frame = pd.concat(data_frames, ignore_index=True)

    # read in the true values
    folder_path_true = folder_path_true  # Replace with the actual path to your folder

    # Loop through the range of numbers from 1 to 10
    data_frames_true = []

    # Loop through the range of numbers from 1 to 10
    for i in range(1, 11):
        csv_file = f'pred_cv_{i}.csv'
        file_path = os.path.join(folder_path_true, csv_file)

        # Check if the file exists before attempting to read
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            df['fold'] = i  # Add a new 'fold' column with the current fold number
            data_frames_true.append(df)
        else:
            print(f"File {csv_file} not found.")

    # Concatenate all data frames into one
    big_data_frame_true = pd.concat(data_frames_true, ignore_index=True)

    final_df = pd.merge(big_data_frame, big_data_frame_true, on=["Drug name", "Sanger ID"])

    final_df = final_df[["Drug name", "Sanger ID", "result", "IC50", "fold_y"]]
    return final_df



def getperform(dataset):
    rmse_list = []
    pearson_corr_list = []
    r2_list = []

    # Group the DataFrame by 'fold' column
    grouped = dataset.groupby('fold_y')

    # Iterate over each group (fold)
    for fold, group_data in grouped:
        true_values = group_data['IC50']
        predicted_values = group_data['result']

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

# metrics_df = getperform(final_df)

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

# getstatistics(metrics_df)

ogdsc_256b_mask_comb = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_comb",
              'OGDSC','256b','')
ogdsc_256b_mask_cell = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_cell",
              'OGDSC','256b','mc')
ogdsc_256b_mask_drug = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_drug",
              'OGDSC','256b','md')
ogdsc_256bnd_mask_comb = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_comb",
              'OGDSC','256bnd','')
ogdsc_256bnd_mask_cell = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_cell",
              'OGDSC','256bnd','mc')
ogdsc_256bnd_mask_drug = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_drug",
              'OGDSC','256bnd','md')
ogdsc_256bnd2_mask_comb = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_comb",
              'OGDSC','256bnd','')
ogdsc_256bnd2_mask_cell = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_cell",
              'OGDSC','256bnd','mc')
ogdsc_256bnd2_mask_drug = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_drug",
              'OGDSC','256bnd','md')
ogdsc_256bnd3_mask_comb = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_comb",
              'OGDSC','256bnd','')
ogdsc_256bnd3_mask_cell = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_cell",
              'OGDSC','256bnd','mc')
ogdsc_256bnd3_mask_drug = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_drug",
              'OGDSC','256bnd','md')


ogdsc_512b_mask_comb = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_comb",
              'OGDSC','512b','')
ogdsc_512b_mask_cell = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_cell",
              'OGDSC','512b','mc')
ogdsc_512b_mask_drug = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_drug",
              'OGDSC','512b','md')
ogdsc_512bnd_mask_comb = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_comb",
              'OGDSC','512bnd','')
ogdsc_512bnd_mask_cell = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_cell",
              'OGDSC','512bnd','mc')
ogdsc_512bnd_mask_drug = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_drug",
              'OGDSC','512bnd','md')

ogdsc_512bnd2_mask_comb = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_comb",
              'OGDSC','512bnd','')
ogdsc_512bnd2_mask_cell = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_cell",
              'OGDSC','512bnd','mc')
ogdsc_512bnd2_mask_drug = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_drug",
              'OGDSC','512bnd','md')
ogdsc_512bnd3_mask_comb = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_comb",
              'OGDSC','512bnd','')
ogdsc_512bnd3_mask_cell = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_cell",
              'OGDSC','512bnd','mc')
ogdsc_512bnd3_mask_drug = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_drug",
              'OGDSC','512bnd','md')

ogdsc_1024b_mask_comb = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_comb",
              'OGDSC','1024b','')
ogdsc_1024b_mask_cell = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_cell",
              'OGDSC','1024b','mc')
ogdsc_1024b_mask_drug = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_drug",
              'OGDSC','1024b','md')
ogdsc_1024bnd_mask_comb = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_comb",
              'OGDSC','1024bnd','')
ogdsc_1024bnd_mask_cell = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_cell",
              'OGDSC','1024bnd','mc')
ogdsc_1024bnd_mask_drug = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_drug",
              'OGDSC','1024bnd','md')
ogdsc_1024bnd2_mask_comb = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_comb",
              'OGDSC','1024bnd','')
ogdsc_1024bnd2_mask_cell = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_cell",
              'OGDSC','1024bnd','mc')
ogdsc_1024bnd2_mask_drug = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_drug",
              'OGDSC','1024bnd','md')
ogdsc_1024bnd3_mask_comb = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_comb",
              'OGDSC','1024bnd','')
ogdsc_1024bnd3_mask_cell = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_cell",
              'OGDSC','1024bnd','mc')
ogdsc_1024bnd3_mask_drug = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_drug",
              'OGDSC','1024bnd','md')

ogdsc_pbfp_mask_comb = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_comb",
              'OGDSC','pbfp','')
ogdsc_pbfp_mask_cell = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_cell",
              'OGDSC','pbfp','mc')
ogdsc_pbfp_mask_drug = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_drug",
              'OGDSC','pbfp','md')
ogdsc_pbfpnd_mask_comb = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_comb",
              'OGDSC','pbfpnd','')
ogdsc_pbfpnd_mask_cell = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_cell",
              'OGDSC','pbfpnd','mc')
ogdsc_pbfpnd_mask_drug = getdf("./HiDRA/pred_results",
              "./HiDRA/training_mask_drug",
              'OGDSC','pbfpnd','md')
ogdsc_pbfpnd2_mask_comb = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_comb",
              'OGDSC','pbfpnd','')
ogdsc_pbfpnd2_mask_cell = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_cell",
              'OGDSC','pbfpnd','mc')
ogdsc_pbfpnd2_mask_drug = getdf("./HiDRA/pred_results_nd2",
              "./HiDRA/training_mask_drug",
              'OGDSC','pbfpnd','md')
ogdsc_pbfpnd3_mask_comb = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_comb",
              'OGDSC','pbfpnd','')
ogdsc_pbfpnd3_mask_cell = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_cell",
              'OGDSC','pbfpnd','mc')
ogdsc_pbfpnd3_mask_drug = getdf("./HiDRA/pred_results_nd3",
              "./HiDRA/training_mask_drug",
              'OGDSC','pbfpnd','md')

ogdsc_2048b_mask_comb = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_comb",
              'OGDSC','2048b','')
ogdsc_2048b_mask_cell = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_cell",
              'OGDSC','2048b','mc')
ogdsc_2048b_mask_drug = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_drug",
              'OGDSC','2048b','md')
ogdsc_2048bnd_mask_comb = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_comb",
              'OGDSC','2048bnd','')
ogdsc_2048bnd_mask_cell = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_cell",
              'OGDSC','2048bnd','mc')
ogdsc_2048bnd_mask_drug = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_drug",
              'OGDSC','2048bnd','md')
ogdsc_2048bnd2_mask_comb = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_comb",
              'OGDSC','2048bnd2','')
ogdsc_2048bnd2_mask_cell = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_cell",
              'OGDSC','2048bnd2','mc')
ogdsc_2048bnd2_mask_drug = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_drug",
              'OGDSC','2048bnd2','md')
ogdsc_2048bnd3_mask_comb = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_comb",
              'OGDSC','2048bnd3','')
ogdsc_2048bnd3_mask_cell = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_cell",
              'OGDSC','2048bnd3','mc')
ogdsc_2048bnd3_mask_drug = getdf("./HiDRA/pred_results_2048b",
              "./HiDRA/training_mask_drug",
              'OGDSC','2048bnd3','md')


### get performance
metrics_ogdsc_256b_mask_comb = getperform(ogdsc_256b_mask_comb)
metrics_ogdsc_256b_mask_cell = getperform(ogdsc_256b_mask_cell)
metrics_ogdsc_256b_mask_drug = getperform(ogdsc_256b_mask_drug)
metrics_ogdsc_256bnd_mask_comb = getperform(ogdsc_256bnd_mask_comb)
metrics_ogdsc_256bnd_mask_cell = getperform(ogdsc_256bnd_mask_cell)
metrics_ogdsc_256bnd_mask_drug = getperform(ogdsc_256bnd_mask_drug)
metrics_ogdsc_256bnd2_mask_comb = getperform(ogdsc_256bnd2_mask_comb)
metrics_ogdsc_256bnd2_mask_cell = getperform(ogdsc_256bnd2_mask_cell)
metrics_ogdsc_256bnd2_mask_drug = getperform(ogdsc_256bnd2_mask_drug)
metrics_ogdsc_256bnd3_mask_comb = getperform(ogdsc_256bnd3_mask_comb)
metrics_ogdsc_256bnd3_mask_cell = getperform(ogdsc_256bnd3_mask_cell)
metrics_ogdsc_256bnd3_mask_drug = getperform(ogdsc_256bnd3_mask_drug)

metrics_ogdsc_512b_mask_comb = getperform(ogdsc_512b_mask_comb)
metrics_ogdsc_512b_mask_cell = getperform(ogdsc_512b_mask_cell)
metrics_ogdsc_512b_mask_drug = getperform(ogdsc_512b_mask_drug)
metrics_ogdsc_512bnd_mask_comb = getperform(ogdsc_512bnd_mask_comb)
metrics_ogdsc_512bnd_mask_cell = getperform(ogdsc_512bnd_mask_cell)
metrics_ogdsc_512bnd_mask_drug = getperform(ogdsc_512bnd_mask_drug)
metrics_ogdsc_512bnd2_mask_comb = getperform(ogdsc_512bnd2_mask_comb)
metrics_ogdsc_512bnd2_mask_cell = getperform(ogdsc_512bnd2_mask_cell)
metrics_ogdsc_512bnd2_mask_drug = getperform(ogdsc_512bnd2_mask_drug)
metrics_ogdsc_512bnd3_mask_comb = getperform(ogdsc_512bnd3_mask_comb)
metrics_ogdsc_512bnd3_mask_cell = getperform(ogdsc_512bnd3_mask_cell)
metrics_ogdsc_512bnd3_mask_drug = getperform(ogdsc_512bnd3_mask_drug)

metrics_ogdsc_1024b_mask_comb = getperform(ogdsc_1024b_mask_comb)
metrics_ogdsc_1024b_mask_cell = getperform(ogdsc_1024b_mask_cell)
metrics_ogdsc_1024b_mask_drug = getperform(ogdsc_1024b_mask_drug)
metrics_ogdsc_1024bnd_mask_comb = getperform(ogdsc_1024bnd_mask_comb)
metrics_ogdsc_1024bnd_mask_cell = getperform(ogdsc_1024bnd_mask_cell)
metrics_ogdsc_1024bnd_mask_drug = getperform(ogdsc_1024bnd_mask_drug)
metrics_ogdsc_1024bnd2_mask_comb = getperform(ogdsc_1024bnd2_mask_comb)
metrics_ogdsc_1024bnd2_mask_cell = getperform(ogdsc_1024bnd2_mask_cell)
metrics_ogdsc_1024bnd2_mask_drug = getperform(ogdsc_1024bnd2_mask_drug)
metrics_ogdsc_1024bnd3_mask_comb = getperform(ogdsc_1024bnd3_mask_comb)
metrics_ogdsc_1024bnd3_mask_cell = getperform(ogdsc_1024bnd3_mask_cell)
metrics_ogdsc_1024bnd3_mask_drug = getperform(ogdsc_1024bnd3_mask_drug)

metrics_ogdsc_pbfp_mask_comb = getperform(ogdsc_pbfp_mask_comb)
metrics_ogdsc_pbfp_mask_cell = getperform(ogdsc_pbfp_mask_cell)
metrics_ogdsc_pbfp_mask_drug = getperform(ogdsc_pbfp_mask_drug)
metrics_ogdsc_pbfpnd_mask_comb = getperform(ogdsc_pbfpnd_mask_comb)
metrics_ogdsc_pbfpnd_mask_cell = getperform(ogdsc_pbfpnd_mask_cell)
metrics_ogdsc_pbfpnd_mask_drug = getperform(ogdsc_pbfpnd_mask_drug)
metrics_ogdsc_pbfpnd2_mask_comb = getperform(ogdsc_pbfpnd2_mask_comb)
metrics_ogdsc_pbfpnd2_mask_cell = getperform(ogdsc_pbfpnd2_mask_cell)
metrics_ogdsc_pbfpnd2_mask_drug = getperform(ogdsc_pbfpnd2_mask_drug)
metrics_ogdsc_pbfpnd3_mask_comb = getperform(ogdsc_pbfpnd3_mask_comb)
metrics_ogdsc_pbfpnd3_mask_cell = getperform(ogdsc_pbfpnd3_mask_cell)
metrics_ogdsc_pbfpnd3_mask_drug = getperform(ogdsc_pbfpnd3_mask_drug)

metrics_ogdsc_2048b_mask_comb = getperform(ogdsc_2048b_mask_comb)
metrics_ogdsc_2048b_mask_cell = getperform(ogdsc_2048b_mask_cell)
metrics_ogdsc_2048b_mask_drug = getperform(ogdsc_2048b_mask_drug)
metrics_ogdsc_2048bnd_mask_comb = getperform(ogdsc_2048bnd_mask_comb)
metrics_ogdsc_2048bnd_mask_cell = getperform(ogdsc_2048bnd_mask_cell)
metrics_ogdsc_2048bnd_mask_drug = getperform(ogdsc_2048bnd_mask_drug)
metrics_ogdsc_2048bnd2_mask_comb = getperform(ogdsc_2048bnd2_mask_comb)
metrics_ogdsc_2048bnd2_mask_cell = getperform(ogdsc_2048bnd2_mask_cell)
metrics_ogdsc_2048bnd2_mask_drug = getperform(ogdsc_2048bnd2_mask_drug)
metrics_ogdsc_2048bnd3_mask_comb = getperform(ogdsc_2048bnd3_mask_comb)
metrics_ogdsc_2048bnd3_mask_cell = getperform(ogdsc_2048bnd3_mask_cell)
metrics_ogdsc_2048bnd3_mask_drug = getperform(ogdsc_2048bnd3_mask_drug)

getstatistics(metrics_ogdsc_256b_mask_comb)
getstatistics(metrics_ogdsc_256b_mask_cell)
getstatistics(metrics_ogdsc_256b_mask_drug)
getstatistics(metrics_ogdsc_256bnd_mask_comb)
getstatistics(metrics_ogdsc_256bnd_mask_cell)
getstatistics(metrics_ogdsc_256bnd_mask_drug)
getstatistics(metrics_ogdsc_256bnd2_mask_comb)
getstatistics(metrics_ogdsc_256bnd2_mask_cell)
getstatistics(metrics_ogdsc_256bnd2_mask_drug)
getstatistics(metrics_ogdsc_256bnd3_mask_comb)
getstatistics(metrics_ogdsc_256bnd3_mask_cell)
getstatistics(metrics_ogdsc_256bnd3_mask_drug)

getstatistics(metrics_ogdsc_512b_mask_comb)
getstatistics(metrics_ogdsc_512b_mask_cell)
getstatistics(metrics_ogdsc_512b_mask_drug)
getstatistics(metrics_ogdsc_512bnd_mask_comb)
getstatistics(metrics_ogdsc_512bnd_mask_cell)
getstatistics(metrics_ogdsc_512bnd_mask_drug)
getstatistics(metrics_ogdsc_512bnd2_mask_comb)
getstatistics(metrics_ogdsc_512bnd2_mask_cell)
getstatistics(metrics_ogdsc_512bnd2_mask_drug)
getstatistics(metrics_ogdsc_512bnd3_mask_comb)
getstatistics(metrics_ogdsc_512bnd3_mask_cell)
getstatistics(metrics_ogdsc_512bnd3_mask_drug)

getstatistics(metrics_ogdsc_1024b_mask_comb)
getstatistics(metrics_ogdsc_1024b_mask_cell)
getstatistics(metrics_ogdsc_1024b_mask_drug)
getstatistics(metrics_ogdsc_1024bnd2_mask_comb)
getstatistics(metrics_ogdsc_1024bnd2_mask_cell)
getstatistics(metrics_ogdsc_1024bnd2_mask_drug)
getstatistics(metrics_ogdsc_1024bnd3_mask_comb)
getstatistics(metrics_ogdsc_1024bnd3_mask_cell)
getstatistics(metrics_ogdsc_1024bnd3_mask_drug)

getstatistics(metrics_ogdsc_pbfp_mask_comb)
getstatistics(metrics_ogdsc_pbfp_mask_cell)
getstatistics(metrics_ogdsc_pbfp_mask_drug)
getstatistics(metrics_ogdsc_pbfpnd_mask_comb)
getstatistics(metrics_ogdsc_pbfpnd_mask_cell)
getstatistics(metrics_ogdsc_pbfpnd_mask_drug)
getstatistics(metrics_ogdsc_pbfpnd2_mask_comb)
getstatistics(metrics_ogdsc_pbfpnd2_mask_cell)
getstatistics(metrics_ogdsc_pbfpnd2_mask_drug)
getstatistics(metrics_ogdsc_pbfpnd3_mask_comb)
getstatistics(metrics_ogdsc_pbfpnd3_mask_cell)
getstatistics(metrics_ogdsc_pbfpnd3_mask_drug)

getstatistics(metrics_ogdsc_2048b_mask_comb)
getstatistics(metrics_ogdsc_2048b_mask_cell)
getstatistics(metrics_ogdsc_2048b_mask_drug)
getstatistics(metrics_ogdsc_2048bnd_mask_comb)
getstatistics(metrics_ogdsc_2048bnd_mask_cell)
getstatistics(metrics_ogdsc_2048bnd_mask_drug)
getstatistics(metrics_ogdsc_2048bnd2_mask_comb)
getstatistics(metrics_ogdsc_2048bnd2_mask_cell)
getstatistics(metrics_ogdsc_2048bnd2_mask_drug)
getstatistics(metrics_ogdsc_2048bnd3_mask_comb)
getstatistics(metrics_ogdsc_2048bnd3_mask_cell)
getstatistics(metrics_ogdsc_2048bnd3_mask_drug)

### get p-value from paired wilcoxon test
get_pairedWilcoxStat(metrics_ogdsc_256b_mask_comb,metrics_ogdsc_256bnd_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_256b_mask_cell,metrics_ogdsc_256bnd_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_256b_mask_drug,metrics_ogdsc_256bnd_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_256b_mask_comb,metrics_ogdsc_256bnd2_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_256b_mask_cell,metrics_ogdsc_256bnd2_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_256b_mask_drug,metrics_ogdsc_256bnd2_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_256b_mask_comb,metrics_ogdsc_256bnd3_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_256b_mask_cell,metrics_ogdsc_256bnd3_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_256b_mask_drug,metrics_ogdsc_256bnd3_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_512b_mask_comb,metrics_ogdsc_512bnd_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_512b_mask_cell,metrics_ogdsc_512bnd_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_512b_mask_drug,metrics_ogdsc_512bnd_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_512b_mask_comb,metrics_ogdsc_512bnd2_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_512b_mask_cell,metrics_ogdsc_512bnd2_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_512b_mask_drug,metrics_ogdsc_512bnd2_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_512b_mask_comb,metrics_ogdsc_512bnd3_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_512b_mask_cell,metrics_ogdsc_512bnd3_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_512b_mask_drug,metrics_ogdsc_512bnd3_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_1024b_mask_comb,metrics_ogdsc_1024bnd_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_1024b_mask_cell,metrics_ogdsc_1024bnd_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_1024b_mask_drug,metrics_ogdsc_1024bnd_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_1024b_mask_comb,metrics_ogdsc_1024bnd2_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_1024b_mask_cell,metrics_ogdsc_1024bnd2_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_1024b_mask_drug,metrics_ogdsc_1024bnd2_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_1024b_mask_comb,metrics_ogdsc_1024bnd3_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_1024b_mask_cell,metrics_ogdsc_1024bnd3_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_1024b_mask_drug,metrics_ogdsc_1024bnd3_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_comb,metrics_ogdsc_pbfpnd_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_cell,metrics_ogdsc_pbfpnd_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_drug,metrics_ogdsc_pbfpnd_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_comb,metrics_ogdsc_pbfpnd2_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_cell,metrics_ogdsc_pbfpnd2_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_drug,metrics_ogdsc_pbfpnd2_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_comb,metrics_ogdsc_pbfpnd3_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_cell,metrics_ogdsc_pbfpnd3_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_drug,metrics_ogdsc_pbfpnd3_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_comb,metrics_ogdsc_1024b_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_cell,metrics_ogdsc_1024b_mask_cell) # 0.8457, 0.8457, 0.8457
get_pairedWilcoxStat(metrics_ogdsc_pbfp_mask_drug,metrics_ogdsc_1024b_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_2048b_mask_comb,metrics_ogdsc_2048bnd_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_2048b_mask_cell,metrics_ogdsc_2048bnd_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_2048b_mask_drug,metrics_ogdsc_2048bnd_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_2048b_mask_comb,metrics_ogdsc_2048bnd2_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_2048b_mask_cell,metrics_ogdsc_2048bnd2_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_2048b_mask_drug,metrics_ogdsc_2048bnd2_mask_drug)

get_pairedWilcoxStat(metrics_ogdsc_2048b_mask_comb,metrics_ogdsc_2048bnd3_mask_comb)
get_pairedWilcoxStat(metrics_ogdsc_2048b_mask_cell,metrics_ogdsc_2048bnd3_mask_cell)
get_pairedWilcoxStat(metrics_ogdsc_2048b_mask_drug,metrics_ogdsc_2048bnd3_mask_drug)

### output all these compact data frames for further plots.
# pbfp
ogdsc_pbfpnd_mask_comb.to_csv("./HiDRA/analysis/compact_results/pbfpnd_mp.csv")
ogdsc_pbfpnd_mask_cell.to_csv("./HiDRA/analysis/compact_results/pbfpnd_mc.csv")
ogdsc_pbfpnd_mask_drug.to_csv("./HiDRA/analysis/compact_results/pbfpnd_md.csv")
ogdsc_pbfpnd2_mask_comb.to_csv("./HiDRA/analysis/compact_results/pbfpnd2_mp.csv")
ogdsc_pbfpnd2_mask_cell.to_csv("./HiDRA/analysis/compact_results/pbfpnd2_mc.csv")
ogdsc_pbfpnd2_mask_drug.to_csv("./HiDRA/analysis/compact_results/pbfpnd2_md.csv")
ogdsc_pbfpnd3_mask_comb.to_csv("./HiDRA/analysis/compact_results/pbfpnd3_mp.csv")
ogdsc_pbfpnd3_mask_cell.to_csv("./HiDRA/analysis/compact_results/pbfpnd3_mc.csv")
ogdsc_pbfpnd3_mask_drug.to_csv("./HiDRA/analysis/compact_results/pbfpnd3_md.csv")

ogdsc_pbfp_mask_comb.to_csv("./HiDRA/analysis/compact_results/pbfp_mp.csv")
ogdsc_pbfp_mask_cell.to_csv("./HiDRA/analysis/compact_results/pbfp_mc.csv")
ogdsc_pbfp_mask_drug.to_csv("./HiDRA/analysis/compact_results/pbfp_md.csv")

# 256b
ogdsc_256bnd_mask_comb.to_csv("./HiDRA/analysis/compact_results/256bnd_mp.csv")
ogdsc_256bnd_mask_cell.to_csv("./HiDRA/analysis/compact_results/256bnd_mc.csv")
ogdsc_256bnd_mask_drug.to_csv("./HiDRA/analysis/compact_results/256bnd_md.csv")
ogdsc_256bnd2_mask_comb.to_csv("./HiDRA/analysis/compact_results/256bnd2_mp.csv")
ogdsc_256bnd2_mask_cell.to_csv("./HiDRA/analysis/compact_results/256bnd2_mc.csv")
ogdsc_256bnd2_mask_drug.to_csv("./HiDRA/analysis/compact_results/256bnd2_md.csv")
ogdsc_256bnd3_mask_comb.to_csv("./HiDRA/analysis/compact_results/256bnd3_mp.csv")
ogdsc_256bnd3_mask_cell.to_csv("./HiDRA/analysis/compact_results/256bnd3_mc.csv")
ogdsc_256bnd3_mask_drug.to_csv("./HiDRA/analysis/compact_results/256bnd3_md.csv")

ogdsc_256b_mask_comb.to_csv("./HiDRA/analysis/compact_results/256b_mp.csv")
ogdsc_256b_mask_cell.to_csv("./HiDRA/analysis/compact_results/256b_mc.csv")
ogdsc_256b_mask_drug.to_csv("./HiDRA/analysis/compact_results/256b_md.csv")

# 512b
ogdsc_512bnd_mask_comb.to_csv("./HiDRA/analysis/compact_results/512bnd_mp.csv")
ogdsc_512bnd_mask_cell.to_csv("./HiDRA/analysis/compact_results/512bnd_mc.csv")
ogdsc_512bnd_mask_drug.to_csv("./HiDRA/analysis/compact_results/512bnd_md.csv")
ogdsc_512bnd2_mask_comb.to_csv("./HiDRA/analysis/compact_results/512bnd2_mp.csv")
ogdsc_512bnd2_mask_cell.to_csv("./HiDRA/analysis/compact_results/512bnd2_mc.csv")
ogdsc_512bnd2_mask_drug.to_csv("./HiDRA/analysis/compact_results/512bnd2_md.csv")
ogdsc_512bnd3_mask_comb.to_csv("./HiDRA/analysis/compact_results/512bnd3_mp.csv")
ogdsc_512bnd3_mask_cell.to_csv("./HiDRA/analysis/compact_results/512bnd3_mc.csv")
ogdsc_512bnd3_mask_drug.to_csv("./HiDRA/analysis/compact_results/512bnd3_md.csv")

ogdsc_512b_mask_comb.to_csv("./HiDRA/analysis/compact_results/512b_mp.csv")
ogdsc_512b_mask_cell.to_csv("./HiDRA/analysis/compact_results/512b_mc.csv")
ogdsc_512b_mask_drug.to_csv("./HiDRA/analysis/compact_results/512b_md.csv")

# 1024b
ogdsc_1024bnd_mask_comb.to_csv("./HiDRA/analysis/compact_results/1024bnd_mp.csv")
ogdsc_1024bnd_mask_cell.to_csv("./HiDRA/analysis/compact_results/1024bnd_mc.csv")
ogdsc_1024bnd_mask_drug.to_csv("./HiDRA/analysis/compact_results/1024bnd_md.csv")
ogdsc_1024bnd2_mask_comb.to_csv("./HiDRA/analysis/compact_results/1024bnd2_mp.csv")
ogdsc_1024bnd2_mask_cell.to_csv("./HiDRA/analysis/compact_results/1024bnd2_mc.csv")
ogdsc_1024bnd2_mask_drug.to_csv("./HiDRA/analysis/compact_results/1024bnd2_md.csv")
ogdsc_1024bnd3_mask_comb.to_csv("./HiDRA/analysis/compact_results/1024bnd3_mp.csv")
ogdsc_1024bnd3_mask_cell.to_csv("./HiDRA/analysis/compact_results/1024bnd3_mc.csv")
ogdsc_1024bnd3_mask_drug.to_csv("./HiDRA/analysis/compact_results/1024bnd3_md.csv")

ogdsc_1024b_mask_comb.to_csv("./HiDRA/analysis/compact_results/1024b_mp.csv")
ogdsc_1024b_mask_cell.to_csv("./HiDRA/analysis/compact_results/1024b_mc.csv")
ogdsc_1024b_mask_drug.to_csv("./HiDRA/analysis/compact_results/1024b_md.csv")

# 2048b
ogdsc_2048bnd_mask_comb.to_csv("./HiDRA/analysis/compact_results/2048bnd_mp.csv")
ogdsc_2048bnd_mask_cell.to_csv("./HiDRA/analysis/compact_results/2048bnd_mc.csv")
ogdsc_2048bnd_mask_drug.to_csv("./HiDRA/analysis/compact_results/2048bnd_md.csv")
ogdsc_2048bnd2_mask_comb.to_csv("./HiDRA/analysis/compact_results/2048bnd2_mp.csv")
ogdsc_2048bnd2_mask_cell.to_csv("./HiDRA/analysis/compact_results/2048bnd2_mc.csv")
ogdsc_2048bnd2_mask_drug.to_csv("./HiDRA/analysis/compact_results/2048bnd2_md.csv")
ogdsc_2048bnd3_mask_comb.to_csv("./HiDRA/analysis/compact_results/2048bnd3_mp.csv")
ogdsc_2048bnd3_mask_cell.to_csv("./HiDRA/analysis/compact_results/2048bnd3_mc.csv")
ogdsc_2048bnd3_mask_drug.to_csv("./HiDRA/analysis/compact_results/2048bnd3_md.csv")

ogdsc_2048b_mask_comb.to_csv("./HiDRA/analysis/compact_results/2048b_mp.csv")
ogdsc_2048b_mask_cell.to_csv("./HiDRA/analysis/compact_results/2048b_mc.csv")
ogdsc_2048b_mask_drug.to_csv("./HiDRA/analysis/compact_results/2048b_md.csv")

