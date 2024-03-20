# PaccMann

Source code: <https://github.com/PaccMann/paccmann_predictor>.

## Platform

Python 3.9.7

## Details about each file and folder
The Python code for training, predicting, and summary results can be found in the code folder. The author of Paccmann provides instructions and examples, minor changes are adapted to fit the GDSC dataset. The gdsc_old folder contains data and preprocessing of the 10 times cross-validation data. SMILES were scraped from PubChem, see the gdsc_old.ipynb for details. The response data were not uploaded due to size limits. 

The command folder contains the bash commands used to train and test the tasks in 10-fold cross-validation.
Running the training command will generate 10 paralleled training processes and automatically save the best parameter in the home directory. Then run the test command to test the model performance. 

The result folder contains the prediction result of 10-fold cvs. As well as the model performance merit calculation.

Please download everything below before running code.

### Additional Data
Training and testing data in 3 settings, mask cell, mask drug and mask combination. Smiles and gene expression data.



