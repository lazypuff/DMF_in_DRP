#### edited by XIAO Meisheng, Jun 11.
import pandas as pd
from scipy.stats import zscore

# from GDSC data portal(https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html)
# download the 'RMA normalised basal expression profiles for all the cell-lines'
expression_df=pd.read_csv('./HiDRA/preprocess/Cell_line_RMA_proc_basalExp.txt'
                          ,sep='\t',index_col=0)
expression_df=expression_df[expression_df.columns[1:]]
use_cols = expression_df.columns.to_list()
number_only_cols = [item[5:] for item in use_cols]
expression_df.columns = number_only_cols

# Loading cell line information
# Cell line information file is 'Annotated list of cell-lines'
# from GDSC data portal(https://www.cancerrxgene.org/downloads/bulk_download)
# download the GDSC1 and GDSC2 Drug Screening IC50s datasets.
gdsc2_resp_raw = pd.read_excel('./HiDRA/preprocess/GDSC2_fitted_dose_response_24Jul22.xlsx')
gdsc1_resp_raw = pd.read_excel('./HiDRA/preprocess/GDSC1_fitted_dose_response_24Jul22.xlsx')
GDSC_response_raw = pd.concat([gdsc1_resp_raw,gdsc2_resp_raw]) # 575197 records
# pick out the common cell line info
cellline_information = GDSC_response_raw[['COSMIC_ID','CELL_LINE_NAME','SANGER_MODEL_ID']]
cellline_information = cellline_information.drop_duplicates()
cellline_information.columns=['COSMIC identifier','Cell line name','Sanger ID']
## read in the ccl_list in the our common dataset
resp_data = pd.read_csv('/Users/meishengxiao/PycharmProjects/pythonProject/Identical_Fold/Preprocessing_Code/PathDSP/GDSC_response.csv')
ccl_list = pd.unique(resp_data['Sanger ID'])
ccl_list.sort()
cellline_information = cellline_information[cellline_information['Sanger ID'].isin(ccl_list)]


#Excluding cell lines whose expression values are not valid
cellline_list=cellline_information['COSMIC identifier']
cellline_list=[str(x) for x in cellline_list]
cosmic_list=expression_df.columns
isin_list=[(cosmic in cellline_list) for cosmic in cosmic_list]
expression_df=expression_df.loc[:,isin_list]

#Excluding expressions that are not gene
expression_list=expression_df.index
expression_list=[str(x) for x in expression_list]
isin_list=[x!='nan' for x in expression_list]
expression_df=expression_df.loc[isin_list]

#Converting COSMIC identifier into Sanger ID
cellline_name_dic={}
for idx,x in cellline_information.iterrows():
    cellline_name_dic[str(x['COSMIC identifier'])]=x['Sanger ID']
cosmic_list=expression_df.columns
cellline_new_col=[cellline_name_dic[cosmic] for cosmic in cosmic_list]
expression_df.columns=cellline_new_col

#Transform expression values into z-score
expression_df=expression_df.apply(zscore)

expression_df.index=expression_df.index.rename('Gene_Symbol')

#Output
expression_df.to_csv('./HiDRA/expression.csv')
