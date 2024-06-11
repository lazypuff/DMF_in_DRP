# HiDRA Preprocessing

The source code for HiDRA can be downloaded from <https://pubs.acs.org/doi/10.1021/acs.jcim.1c00706?goto=supporting-info>, in the supporting information section.

## Platform

Python 3.9.7

## Raw data needed to run the preprocessing code
Please download from here: https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html, to retrieve the expression data from GDSC, which should be called as "Cell_line_RMA_proc_basalExp.txt".

Please download from here: https://www.cancerrxgene.org/downloads/bulk_download, to get the the GDSC1 and GDSC2 Drug Screening IC50s datasets, which should be called as "GDSC1_fitted_dose_response_ddmmyy.xlsx" and "GDSC2_fitted_dose_response_ddmmyy.xlsx"

## data contained in the parent folder
GDSC_response.csv: used in "Preprocessing_hiDRA.py" to match cell line records with all other models.

## Steps
### Step 1
Please download from here: https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html, to retrieve the expression data from GDSC, which should be called as "Cell_line_RMA_proc_basalExp.txt".

Please also download from here: https://www.cancerrxgene.org/downloads/bulk_download, to get the the GDSC1 and GDSC2 Drug Screening IC50s datasets, which should be called as "GDSC1_fitted_dose_response_ddmmyy.xlsx" and "GDSC2_fitted_dose_response_ddmmyy.xlsx". Put these three files under this folder.

### Step 2
Run the "Preprocessing_HiDRA.py", to get the required dataset "expression.csv", and put it in the parent folder "HiDRA" with all other data together.

