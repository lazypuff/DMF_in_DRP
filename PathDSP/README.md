# PathDSP

Source code: <https://github.com/TangYiChing/PathDSP>.

## Platform

Python 3.9.7

## Details about each file and folder
Please download everything below before running code. Thank you!
### Additional Data
The data needed for training and predicting PathDSP in this study can be downloaded from:

In this link, there are five .txt files named as: "pdsp_gdsc_common_xxxx.txt" or "ogdsc_xxxx.txt", where xxxx indicates the configuration (xxxx is left blank for 256-bit Morgan fingerprints).
### Contained in GitHub
#### code
This folder contain all the code needed to run PathDSP in this study. The main code is: FNNv2.py.

#### data
fold_info_mask_xxxx.csv: the fold information for different splitting strategies, where xxxx indicates the splitting strategy.
#### results
15 folders under this folder, which are named as xxxx_mask_yyyy, where xxxx indicates the drug molecular fingerprints type and yyyy indicates the splitting strategy.

Under each sub-folder, three files are here. 

OutputName.Loss.txt: the file that records the training history;

OutputName.Prediction.txt: the file that records the prediction and true value of each drug-cell pair for every fold;

OutputName.best_model.pt: best model after training in PyTorch format.
## Example
Running the main code will give back the training history, predictions for each pair and best model in PyTorch format.

This example is to get the results(training history, predictions and best model) for 1024-bit Morgan fingerprints under mask drug splitting strategy.

python ./code/FNNv2.py -i ./data/ogdsc_1024b.txt -c 10 -o ./results/1024b_mask_drug/output_oGDSC_IF_1024b_mask_drug -f ./data/fold_info_mask_drug.csv

-i: the input data, should be downloaded from the section "Additional Data"; -c: number of CV, please use 10 for the purpose of replicating our study; -o: the desired output name; -f: the fold information used in training.




