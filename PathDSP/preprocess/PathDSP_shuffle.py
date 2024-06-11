import pandas as pd
import random


pdsp_ogdsc_256b = pd.read_csv("./PathDSP/data/256b_wd.txt", delimiter='\t')
pdsp_ogdsc_512b = pd.read_csv("/PathDSP/data/512b_wd.txt", delimiter='\t')
pdsp_ogdsc_1024b = pd.read_csv("/PathDSP/data/1024b_wd.txt", delimiter='\t')
pdsp_ogdsc_pbfp = pd.read_csv("/PathDSP/data/pbfp_wd.txt", delimiter='\t')


def shuffle_chem_columns_by_drug(df, base_seed=66):
    """
    Shuffles 'CHEM_mBit_x' columns in the DataFrame such that all rows with the same 'drug'
    have the same shuffled order of these columns.

    Parameters:
    - df: pandas.DataFrame containing the data.
    - base_seed: An integer seed for reproducibility.

    Returns:
    - A pandas.DataFrame with the 'CHEM_mBit_x' columns shuffled consistently within each 'drug' group.
    """
    random.seed(base_seed)  # Seed for reproducibility

    # Identify the 'CHEM_mBit_' columns to shuffle
    chem_columns = [col for col in df.columns if 'CHEM_mBit_' in col]

    # Copy the DataFrame to avoid changing the original data
    shuffled_df = df.copy()

    # Initialize a dictionary to hold the shuffled columns for each drug
    shuffled_columns_by_drug = {}

    # Iterate through each unique drug value
    for drug in df['drug'].unique():
        # Shuffle the columns for the first occurrence and store the order
        if drug not in shuffled_columns_by_drug:
            shuffled_order = chem_columns.copy()
            random.shuffle(shuffled_order)
            shuffled_columns_by_drug[drug] = shuffled_order

        # Retrieve the shuffled column order for the current drug
        shuffled_order = shuffled_columns_by_drug[drug]

        # Apply the shuffled order to all rows with the current drug value
        for i, col in enumerate(chem_columns):
            shuffled_df.loc[df['drug'] == drug, col] = df.loc[df['drug'] == drug, shuffled_order[i]].values

    return shuffled_df


def check_equal(df1,df2):
    dataframes_are_identical = df1.equals(df2)
    # Step 1: Identify the 'CHEM_mBit_' columns
    chem_columns = [col for col in df1.columns if 'CHEM_mBit_' in col]

    # Step 2: Drop the 'CHEM_mBit_' columns to compare the rest of the DataFrame
    df_original_without_chem = df1.drop(columns=chem_columns)
    df_shuffled_without_chem = df2.drop(columns=chem_columns)

    # Step 3: Compare the DataFrames without the 'CHEM_mBit_' columns
    # This will return True if all the remaining data is identical, False otherwise
    dataframes_are_identical_drop = df_original_without_chem.equals(df_shuffled_without_chem)

    # Output the result
    print("DataFrames are identical (including 'CHEM_mBit_' columns):", dataframes_are_identical)
    print("DataFrames are identical (excluding 'CHEM_mBit_' columns):", dataframes_are_identical_drop)

### 256b
pdsp_ogdsc_256b_shuffle1 = shuffle_chem_columns_by_drug(pdsp_ogdsc_256b, base_seed=66)
check_equal(pdsp_ogdsc_256b,pdsp_ogdsc_256b_shuffle1)
pdsp_ogdsc_256b_shuffle2 = shuffle_chem_columns_by_drug(pdsp_ogdsc_256b, base_seed=77)
check_equal(pdsp_ogdsc_256b,pdsp_ogdsc_256b_shuffle2)
pdsp_ogdsc_256b_shuffle3 = shuffle_chem_columns_by_drug(pdsp_ogdsc_256b, base_seed=88)
check_equal(pdsp_ogdsc_256b,pdsp_ogdsc_256b_shuffle3)
### 512b
pdsp_ogdsc_512b_shuffle1 = shuffle_chem_columns_by_drug(pdsp_ogdsc_512b, base_seed=66)
check_equal(pdsp_ogdsc_512b,pdsp_ogdsc_512b_shuffle1)
pdsp_ogdsc_512b_shuffle2 = shuffle_chem_columns_by_drug(pdsp_ogdsc_512b, base_seed=77)
check_equal(pdsp_ogdsc_512b,pdsp_ogdsc_512b_shuffle2)
pdsp_ogdsc_512b_shuffle3 = shuffle_chem_columns_by_drug(pdsp_ogdsc_512b, base_seed=88)
check_equal(pdsp_ogdsc_512b,pdsp_ogdsc_512b_shuffle3)
### 1024b
pdsp_ogdsc_1024b_shuffle1 = shuffle_chem_columns_by_drug(pdsp_ogdsc_1024b, base_seed=66)
check_equal(pdsp_ogdsc_1024b,pdsp_ogdsc_1024b_shuffle1)
pdsp_ogdsc_1024b_shuffle2 = shuffle_chem_columns_by_drug(pdsp_ogdsc_1024b, base_seed=77)
check_equal(pdsp_ogdsc_1024b,pdsp_ogdsc_1024b_shuffle2)
pdsp_ogdsc_1024b_shuffle3 = shuffle_chem_columns_by_drug(pdsp_ogdsc_1024b, base_seed=88)
check_equal(pdsp_ogdsc_1024b,pdsp_ogdsc_1024b_shuffle3)
### PubChem
pdsp_ogdsc_pbfp_shuffle1 = shuffle_chem_columns_by_drug(pdsp_ogdsc_pbfp, base_seed=66)
check_equal(pdsp_ogdsc_pbfp,pdsp_ogdsc_pbfp_shuffle1)
pdsp_ogdsc_pbfp_shuffle2 = shuffle_chem_columns_by_drug(pdsp_ogdsc_pbfp, base_seed=77)
check_equal(pdsp_ogdsc_pbfp,pdsp_ogdsc_pbfp_shuffle2)
pdsp_ogdsc_pbfp_shuffle3 = shuffle_chem_columns_by_drug(pdsp_ogdsc_pbfp, base_seed=88)
check_equal(pdsp_ogdsc_pbfp,pdsp_ogdsc_pbfp_shuffle3)

## output data frames
### 256b
pdsp_ogdsc_256b_shuffle1.to_csv("./PathDSP/data/256b_s1.txt",sep='\t', index=False)

pdsp_ogdsc_256b_shuffle2.to_csv("./PathDSP/data/256b_s2.txt",sep='\t', index=False)

pdsp_ogdsc_256b_shuffle3.to_csv("./PathDSP/data/256b_s3.txt",sep='\t', index=False)
### 512b
pdsp_ogdsc_512b_shuffle1.to_csv("./PathDSP/data/512b_s1.txt",sep='\t', index=False)

pdsp_ogdsc_512b_shuffle2.to_csv("./PathDSP/data/512b_s2.txt",sep='\t', index=False)

pdsp_ogdsc_512b_shuffle3.to_csv("./PathDSP/data/512b_s3.txt",sep='\t', index=False)
### 1024b
pdsp_ogdsc_1024b_shuffle1.to_csv("./PathDSP/data/1024b_s1.txt",sep='\t', index=False)

pdsp_ogdsc_1024b_shuffle2.to_csv("./PathDSP/data/1024b_s2.txt",sep='\t', index=False)

pdsp_ogdsc_1024b_shuffle3.to_csv("./PathDSP/data/1024b_s3.txt",sep='\t', index=False)
### pbfp
pdsp_ogdsc_pbfp_shuffle1.to_csv("./PathDSP/data/pbfp_s1.txt",sep='\t', index=False)

pdsp_ogdsc_pbfp_shuffle2.to_csv("./PathDSP/data/pbfp_s2.txt",sep='\t', index=False)

pdsp_ogdsc_pbfp_shuffle3.to_csv("./PathDSP/data/pbfp_s3.txt",sep='\t', index=False)

