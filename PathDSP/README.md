# PathDSP

Source code: <https://github.com/TangYiChing/PathDSP>.

## Platform

Python 3.9.7

## Environment

numpy 1.23.5

torch 2.5.1

shap 0.46.0


## Details about each file and folder
Please download everything below before running code. Thank you!
### Additional Data
Due to limit file size requested by GitHub, some input date needs to be preprocessed to retrieve. Please see the "preprocess" folder for details.

After preprocessing, there are 12 .txt files named as: "xxxx_wd.txt", "xxxx_nds1.txt", "xxxx_nds2.txt" or "xxxx_nds3.txt", where xxxx indicates the configuration (256b for 256-bit Morgan fingerprints; 512b for 512-bit Morgan fingerprints; 1024b for 1024-bit Morgan fingerprints; pbfp for PubChem fingerprints).
### Contained in GitHub

#### preprocess
This folder contains all codes that preprocess data and the output dataset needs to be used with the files under the "data" folder

#### code
This folder contain all the code needed to run PathDSP in this study. The main code is: FNNv2.py.

#### data
fold_info_mask_xxxx.csv: the fold information for different splitting strategies, where xxxx indicates the splitting strategy.
#### results
12 folders under this folder, which are named as xxxx_mask_yyyy, where xxxx indicates the drug molecular fingerprints type and yyyy indicates the splitting strategy.

Under each sub-folder, Four files are here. They are prediction results under different experimental settings. each file records the prediction and true value of IC50 for each drug-cell pair across every fold.

OutputName_b.Prediction.txt: The prediction results for the "With-Drug" setting for PathDSP.

OutputName_nds1.Prediction.txt: The prediction results for the "Null-Drug" setting 1 for PathDSP.

OutputName_nds2.Prediction.txt: The prediction results for the "Null-Drug" setting 2 for PathDSP.

OutputName_nds3.Prediction.txt: The prediction results for the "Null-Drug" setting 3 for PathDSP.

These "Null-Drug" settings differ with each other due to different shuffling of the drug molecular fingerprints. Three replicates are designated to measure the variability of the "Null-Drug" setting.

#### analysis

PathDSP_result_analysis.py: the python code to calculate the statistics shown in the paper. Please use it interactively.


## Example
Running the main code will give back the training history, predictions for each pair and best model in PyTorch format.

This example is to get the results(training history, predictions and best model) for third type of shuffled 256-bit Morgan fingerprints under mask drug splitting strategy.

python ./PathDSP/code/FNNv2.py -i 256b_s3.txt -c 10 -o ./PathDSP/output_oGDSC_IF_256bnds3_mask_drug -f fold_info_mask_drug.csv

-i: the input data, should be retrived from the section "preprocess"; -c: number of CV, please use 10 for the purpose of replicating our study; -o: the desired output name; -f: the fold information used in training.




