# PaccMann

Source code: <https://github.com/PaccMann/paccmann_predictor>.

## Platform

Python 3.9.7

## Details about each file and folder

The Python code for training, predicting, and summary results can be found in the code folder. The author of Paccmann provides instructions and examples, minor changes are adapted to fit the GDSC dataset. The gdsc_old folder contains data and preprocessing of the 10 times cross-validation data. SMILES were scraped from PubChem, see the gdsc_old.ipynb for details. The response data were not uploaded due to size limits. 

The code folder contains the source code of PaccMann, minor changes was made to run the model on GDSC dataset. 

### Additional Data

The training data are obtained from the GDSC Database at https://www.cancerrxgene.org/. Processing codes are included. Training and testing data in 3 settings, mask cell, mask drug and mask combination. Smiles and gene expression data.

### contained in GitHub


The command folder contains the bash commands used to train and test the tasks in 10-fold cross-validation.
Running the training command will generate 10 paralleled training processes and automatically save the best parameter in the home directory. Then run the test command to test the model performance. 

The result folder contains the prediction result of 10-fold cvs. As well as the model performance merit calculation.

Please download everything below before running code.

#### code

#### commands

#### paccmann_results

Three folders under this folder, which are "random", "random1" and "random2". The results of "With-Drug" experimental settings are under the folder "random1". Results of three differently shuffled "Null-Drug" settings are under each folder.

#### analysis

compact_results: contained the results for one experiment that combines results across all CV fold. The compact results can be retrieved by running the below code: "make_compact_results_PaccMann.R".

make_compact_results_PaccMann.R: to retrieve the compact results.

PaccMann_getstatistics.py: the python code to calculate the statistics shown in the paper. Please use it interactively.


