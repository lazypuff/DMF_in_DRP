# PathDSP Preprocessing

Source code: <https://github.com/TangYiChing/PathDSP>.

## Platform

Python 3.9.7

## Original data from PathDSP
Please download from here: https://zenodo.org/records/7532963.

## data contained in this folder
GDSC_response.csv: used in "Preprocessing_pdsp.py" to match records with all other models.

ogdsc_xxxx.csv: the fingerprints used to construct the PathDSP input dataset, where xxxx refers to 256-bit, 512-bit, 1024-bit Morgan fingerprints or PubChem fingerprints.

## Steps
### Step 1
Download the original data from PathDSP, which should be called as: "input.txt" from the link: https://zenodo.org/records/7532963.

### Step 2
Run the "Preprocessing_pdsp.py" first, to get the "With-Drug" setting input such as 256b_wd.txt.

### Step 3
Run the "pdsp_shuffle_by_drugs.py" to get the shuffled fingerprints input dataset such as 256b_s1.txt.

Note that there are two shuffle ways for experimental use. "pdsp_shuffle_by_drugs.py" is to randomly assign the fingerprints to drugs, while "PathDSP_shuffle.py" is to shuffle the fingerprints within each drug, which means creating a fingerprint that destroys the structure of fingerprints for each drug.

The shuffle manner in "PathDSP_shuffle.py" is used to monitor if there are any performance improvement purely due to the artifact from deep learning.


