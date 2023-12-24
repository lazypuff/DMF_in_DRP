# ADRML

Source code: <https://github.com/fahmadimoughari/AdrML>.

## Platform

Python 3.7.14

## Details about each file and folder
Please download everything below before running code. Thank you!

### All things are contained in GitHub
#### data
contain all the data used in ADRML in this study.

fold_info_mask_xxxx.csv: csv files that record the identical fold information, where xxxx indicates different splitting strategies.

resp_ogdsc.csv: the drug-cell pairs' response.

simC_ogdsc_gexp.csv: the similarity matrix from gene expression profiles for every pair of cells, calculated by Jaccard Index.

simD_(o)gdsc_xxxx.csv: the similarity matrices from different types of drug molecular fingerprints for every pair of drugs, calculated by Jaccard Index. xxxx indicates the type of drug molecular fingerprints.

#### code
ADRML_IFv2_xx.py: the main code used to train and predict, the output is the results of three metrics of each fold. xx indicates the splitting strategy, and left blank means mask pair/combination.

manifoldv2_mask_xxxx.py: the code plugged in the main code, where xxxx indicates the splitting strategy.

#### results
IF_Ogdsc_xxxx_mask_yyyy.csv: the results of three metrics for every fold in different configurations. xxxx indicates the drug molecular fingerprints type and yyyy indicates the splitting strategy.

## Example
This example is to train the ADRML using 256-bit Morgan fingerprints under mask cell splitting.

python3 ./code/ADRML_IFv2_mc.py response_dirc='./data/resp_ogdsc.csv' simC_dirc='./data/simC_ogdsc_gexp.csv' simD_dirc='./data/simD_gdsc_256b.csv' fold_info='./data/fold_info_mask_cell.csv' dim=0.7 miu=8 lambda=4 CV=10 repeat=1 out_name='IF_Ogdsc_256b_mask_cell.csv'

response_dirc: the input response in LogIC50 for each drug-cell pair;

simC_dirc: the input cell line similarity matrix;

simD_dirc: the input drug similarity matrix;

fold_info: the input fold info for desired splitting strategy;

dim, miu, lambda: the hyper-parameters in ADRML, use dim=0.7, miu=8, lambda=4 as authors suggested so;

CV, repeat: please use CV=10(10-fold CV) and repeat=1(no repeatement) for replicating this study;

out_name: desired name for the output file.





