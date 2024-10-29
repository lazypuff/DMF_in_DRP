import pandas as pd

# read in pathdsp original dataset
pdsp_256b = pd.read_csv('./PathDSP/data/256b_wd.txt',sep='\t')
pdsp_512b = pd.read_csv('./PathDSP/data/512b_wd.txt',sep='\t')
pdsp_1024b = pd.read_csv('./PathDSP/data/1024b_wd.txt',sep='\t')
pdsp_2048b = pd.read_csv('./PathDSP/data/2048b_wd.txt',sep='\t')
pdsp_pbfp = pd.read_csv('./PathDSP/data/pbfp_wd.txt',sep='\t')


def shuffle_by_drugs(pdsp_df,seed):
    # read in the morgan fps
    chem_mbit_columns = [col for col in pdsp_df.columns if col.startswith('CHEM_mBit_')]
    # Add the 'drug' column to the list
    columns_to_select = ['drug'] + chem_mbit_columns

    morganf_pr = pdsp_df[columns_to_select]
    morganf_pr_unique = morganf_pr.drop_duplicates()
    morganf_pr_unique.set_index('drug', inplace=True)
    # shuffle
    morganf_pr_shuffle = morganf_pr_unique.sample(frac=1, random_state=seed)
    morganf_pr_shuffle.index = morganf_pr_unique.index

    pdsp_nd = pdsp_df.drop([col for col in pdsp_df.columns if col.startswith('CHEM_mBit_')], axis=1)
    pdsp_mfp_shuffled = morganf_pr_shuffle.copy()

    ### merge dataset
    pdsp_data = pdsp_nd.merge(pdsp_mfp_shuffled, how='left', left_on='drug', right_index=True)
    # reorder
    cell_index = pdsp_data.columns.get_loc('cell')

    # Get the columns starting with 'CHEM_mBit_#'
    chem_columns = [col for col in pdsp_data.columns if col.startswith('CHEM_mBit_')]

    # Reorder the DataFrame columns
    new_order = list(pdsp_data.columns[:cell_index + 1]) + chem_columns + list(pdsp_data.columns[cell_index + 1:])

    # Create a new DataFrame with the reordered columns
    pdsp_data_reordered = pdsp_data[new_order]

    # Get the index of the 'resp' column
    resp_index = pdsp_data_reordered.columns.get_loc('resp')

    # Keep all columns up to and including 'resp'
    pdsp_data_reordered = pdsp_data_reordered.iloc[:, :resp_index + 1]

    # output the preprocessed dataset
    return pdsp_data_reordered

# 1024b
pdsp1024b_s1 = shuffle_by_drugs(pdsp_1024b,66)
pdsp1024b_s2 = shuffle_by_drugs(pdsp_1024b,77)
pdsp1024b_s3 = shuffle_by_drugs(pdsp_1024b,88)

pdsp1024b_s1.to_csv('./PathDSP/data/1024b_nds1sd.txt', sep='\t', index=False)
pdsp1024b_s2.to_csv('./PathDSP/data/1024b_nds2sd.txt', sep='\t', index=False)
pdsp1024b_s3.to_csv('./PathDSP/data/1024b_nds3sd.txt', sep='\t', index=False)

# 256b
pdsp256b_s1 = shuffle_by_drugs(pdsp_256b,66)
pdsp256b_s2 = shuffle_by_drugs(pdsp_256b,77)
pdsp256b_s3 = shuffle_by_drugs(pdsp_256b,88)

pdsp256b_s1.to_csv('./PathDSP/data/256b_nds1sd.txt', sep='\t', index=False)
pdsp256b_s2.to_csv('./PathDSP/data/256b_nds2sd.txt', sep='\t', index=False)
pdsp256b_s3.to_csv('./PathDSP/data/256b_nds3sd.txt', sep='\t', index=False)

# 512b
pdsp512b_s1 = shuffle_by_drugs(pdsp_512b,66)
pdsp512b_s2 = shuffle_by_drugs(pdsp_512b,77)
pdsp512b_s3 = shuffle_by_drugs(pdsp_512b,88)

pdsp512b_s1.to_csv('./PathDSP/data/512b_nds1sd.txt', sep='\t', index=False)
pdsp512b_s2.to_csv('./PathDSP/data/512b_nds2sd.txt', sep='\t', index=False)
pdsp512b_s3.to_csv('./PathDSP/data/512b_nds3sd.txt', sep='\t', index=False)

# 2048b
pdsp2048b_s1 = shuffle_by_drugs(pdsp_2048b,66)
pdsp2048b_s2 = shuffle_by_drugs(pdsp_2048b,77)
pdsp2048b_s3 = shuffle_by_drugs(pdsp_2048b,88)

pdsp2048b_s1.to_csv('./PathDSP/data/2048b_nds1sd.txt', sep='\t', index=False)
pdsp2048b_s2.to_csv('./PathDSP/data/2048b_nds2sd.txt', sep='\t', index=False)
pdsp2048b_s3.to_csv('./PathDSP/data/2048b_nds3sd.txt', sep='\t', index=False)

# pbfp
pdsppbfp_s1 = shuffle_by_drugs(pdsp_pbfp,66)
pdsppbfp_s2 = shuffle_by_drugs(pdsp_pbfp,77)
pdsppbfp_s3 = shuffle_by_drugs(pdsp_pbfp,88)

pdsppbfp_s1.to_csv('./PathDSP/data/pbfp_nds1sd.txt', sep='\t', index=False)
pdsppbfp_s2.to_csv('./PathDSP/data/pbfp_nds2sd.txt', sep='\t', index=False)
pdsppbfp_s3.to_csv('./PathDSP/data/pbfp_nds3sd.txt', sep='\t', index=False)
