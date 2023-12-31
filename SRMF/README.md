# SRMF

Source code: <https://github.com/linwang1982/SRMF>.

## Platform
MatLab(R2023a, 64-bit)

## Details about each file and folder
### All data are contained in GitHub
Codes, Data and Results are all under the folder.

#### Codes
alg_update.m, CMF.m, compute_loss.m: Three source codes of SRMF about the algorithm.

new_code_10f_xxxx_mask_yyyy.m (xxxx indicates the type of drug molecular fingerprints and yyyy indicates the splitting strategy): main codes that get the results of three metrics and hyper-parameters.

#### Data
identity_gdsc.csv, simD_gdsc_xxxx.csv (xxxx indicates the type of drug molecular fingerprints): drug similarity matrices from different drug molecular fingerprints.

simC_ogdsc_gexp.csv: cell line similarity matrix from gene expression profiles.

resp_ogdsc.csv: the response matrix for drug-cell pairs in logIC50.

srmf_valid_info_mask_xxxx_og (xxxx indicates the splitting strategy): the identical fold infomation used in other studies (here is the information for splitting train and test set rather than splitting train and validation set, sorry for the typo here)

#### Results
results_oldgdsc_xxxx_mask_yyyy.csv (xxxx indicates the type of drug molecular fingerprints and yyyy indicates the splitting strategy): in each .csv file, there are 7 columns, which are fold ID, lambda_l, lambda_d, lambda_c, RMSE, PCC, R2 respectively. 

lambda_l, lambda_d, lambda_c are three hyper-parameters that used in SRMF, and please see the original paper for details about these three hyper-parameters. lambda_d and lambda_c are deliberately tuned not to be 0 as 0 means no contribution from drug/cell similarity matrix.
## Example
This example is to get the result named as "results_oldgdsc_1024b_mask_drug.csv"

matlab -nodesktop -nosplash -singleCompThread -r new_code_10f_1024b_mask_drug -logfile new_code_10f_1024b_md.out

Use "matlab -help" in terminal to check the best options for your working machine to run the codes.


