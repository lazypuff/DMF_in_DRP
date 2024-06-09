# HiDRA

The source code for HiDRA can be downloaded from <https://pubs.acs.org/doi/10.1021/acs.jcim.1c00706?goto=supporting-info>, in the supporting information section.

## Platform

Python 3.9.7

## Details about each file and folder

Please download everything below before running code. Thank you. \### Additional Data expression.csv, GDSC_response.csv and geneset.gmt are in the additional folder that stores data, which can be downloaded from:

expression.csv: The original gene expression provided by GDSC.

### Contained in GitHub

#### code

HiDRA_predict.py: The code of prediction for HiDRA's output model.

HiDRA_training_xxxx.py: Codes for training HiDRA model using different types of drug molecular fingerprints (HiDRA_training.py is used for 512-bit Morgan fingerprints).

#### data

GDSC_response.csv: The original drug response of each drug-cell pair in LogIC50 (Though the colname is "IC50", it is LogIC50).

geneset.gmt:The KEGG pathway data used in HiDRA.

input_dir_oGDSC_xxxx: Folders that contain required gene data for cell lines and fingerprints for drugs in each configuration.

training_mask_xxxx: Folders that contain training data and validation data for each splitting strategy.

pred_mask_xxxx: Folders that contain test data for each splitting strategy.

#### results

pred_results: results of predictions on the test data. (In this folder, the naming system is OGDSC_x\_cv#\_pred_y.csv, where x indicates the type of drug molecular fingerprints, y indicates the splitting strategy and \# indicates integer from 1 to 10. For eaxmple, "OGDSC_256b_cv3_pred.csv" means this is the prediction results for mask-pairs splitting, using 256-bit Morgan fingerprints for fold 3; "OGDSC_256bnd_cv6_pred_md.csv" means this is the prediction results for mask-drugs splitting, using random shuffling 256-bit Morgan fingerprints for fold 6.)

pred_results_nd2: results of predictions on the test data, for the second different random shuffling "Null-Drug" setting. Naming system is the same as above.

pred_results_nd3: results of predictions on the test data, for the third different random shuffling "Null-Drug" setting. Naming system is the same as above.

## Example

This example is to train the HiDRA using PubChem molecular fingerprints using fold 6 under mask drug splitting strategy.

### Training

python HiDRA_training_pbfp.py -t training_mask_drug/train_cv_6.csv -v training_mask_drug/valid_cv_6.csv -i input_dir_oGDSC_pbfp -e 30 -o OGDSC_pbfp_cv6_md.hdf5

-t: training data; -v: validation data; -i: input folder that contain corresponding gene information for cells and drug molecular fingerprints for drugs, used to differentiate splitting strategy; -e: number of epoch; -o: name for output model.

### Predicting

python HiDRA_predict.py -m OGDSC_pbfp_cv6_md.hdf5 -p pred_mask_drug/pred_set_6.csv -i input_dir_oGDSC_pbfp -o pred_results/OGDSC_pbfp_cv6_pred_md.csv

-m: the model using for predictions, which is the output from training; -p: test data; -i: input folder that contain corresponding gene information for cells and drug molecular fingerprints for drugs; -o: name for the output prediction.
