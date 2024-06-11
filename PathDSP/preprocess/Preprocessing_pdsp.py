# Meisheng Jun 10
import pandas as pd
import numpy as np

orig_input = pd.read_csv("./PathDSP/preprocess/input.txt", delimiter='\t')
ogdsc_resp = pd.read_csv("./PathDSP/preprocess/GDSC_response.csv")

# pick out the common drugs and cells
drug_names = ogdsc_resp['Drug name'].unique()
sanger_ids = ogdsc_resp['Sanger ID'].unique()

# Merge orig_input with ogdsc_resp on 'drug' and 'Drug name', and 'cell' and 'Sanger ID'
filtered_input = orig_input.merge(ogdsc_resp, left_on=['drug', 'cell'], right_on=['Drug name', 'Sanger ID'], how='inner')
filtered_input = filtered_input.drop(columns=['Origin_idx', 'Drug name', 'Sanger ID', 'resp'])
filtered_input = filtered_input.rename(columns={'IC50': 'resp'})



# Now add fingerprints info.
## read in fingerprints info
mfp256 = pd.read_csv("./PathDSP/preprocess/ogdsc_256mfp.csv",index_col=0)
mfp512 = pd.read_csv("./PathDSP/preprocess/ogdsc_512mfp.csv",index_col=0)
mfp1024 = pd.read_csv("./PathDSP/preprocess/ogdsc_1024mfp.csv",index_col=0)
pbfp = pd.read_csv("./PathDSP/preprocess/ogdsc_pbfp.csv",index_col=0)

def make_input_pdsp(filter_df, fp):
    # Extract the 'drug' column values from filtered_input
    drug_column = filter_df['drug']

    # Create a DataFrame that maps 'drug' names to the 'fp' values
    chem_values = fp.loc[drug_column].reset_index(drop=True)

    # Determine the number of columns in chem_values and rename them accordingly
    num_columns = chem_values.shape[1]
    chem_columns = [f'CHEM_mBit_{i}' for i in range(num_columns)]
    chem_values.columns = chem_columns

    # Identify original 'CHEM_mBit_x' columns
    original_chem_columns = [col for col in filter_df.columns if col.startswith('CHEM_mBit_')]

    # Drop original 'CHEM_mBit_x' columns
    filter_df = filter_df.drop(columns=original_chem_columns)

    # Concatenate the new columns to filter_df
    filter_df = pd.concat([filter_df.reset_index(drop=True), chem_values], axis=1)

    # Reorder the columns to place the new 'CHEM_mBit_x' columns right after the 'cell' column
    cell_index = filter_df.columns.get_loc('cell')
    reordered_columns = (
        filter_df.columns[:cell_index + 1].tolist() +
        chem_columns +
        filter_df.columns[cell_index + 1:-len(chem_columns)].tolist()
    )
    output = filter_df[reordered_columns]

    return output

pdsp256b = make_input_pdsp(filtered_input,mfp256)
pdsp512b = make_input_pdsp(filtered_input,mfp512)
pdsp1024b = make_input_pdsp(filtered_input,mfp1024)
pdsppbfp = make_input_pdsp(filtered_input,pbfp)

# output
pdsp256b.to_csv("./PathDSP/data/256b_wd.txt",sep='\t', index=False)
pdsp512b.to_csv("./PathDSP/data/512b_wd.txt",sep='\t', index=False)
pdsp1024b.to_csv("./PathDSP/data/1024b_wd.txt",sep='\t', index=False)
pdsppbfp.to_csv("./PathDSP/data/pbfp_wd.txt",sep='\t', index=False)
