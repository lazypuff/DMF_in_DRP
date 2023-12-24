# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 05:09:27 2020

@author: fahma

edited by XIAO Meisheng AUG 2 2023, the loose splitting(mask combination) strategy 10-CV using Identical Folds

sbatch --mail-type ALL --mail-user meisheng@email.unc.edu -p general -N 1 --mem=16g -n 1 -t 1-
--wrap="python3 /nas/longleaf/home/meisheng/FP_Project/code_running/ADRML/mask_drug_code/ADRML_mask_drug_gdsc_512raw.py
response_dirc='/nas/longleaf/home/meisheng/FP_Project/code_running/ADRML/Similarity_based_use_data/resp_gdsc_adrml.csv'
simC_dirc='/nas/longleaf/home/meisheng/FP_Project/code_running/ADRML/Similarity_based_use_data/simC_gdsc_raw.csv'
simD_dirc='/nas/longleaf/home/meisheng/FP_Project/code_running/ADRML/Similarity_based_use_data/simD_gdsc_512b.csv'
fold_info=''
dim=0.7 miu=8 lambda=4 CV=10 repeat=1
out_name='.csv'"
"""

# Importing required libraries
import sys

# Add a directory to the Python path
sys.path.append('/nas/longleaf/home/meisheng/FP_Project/code_running/ADRML/IF_code')
import copy
import numpy as np
import random
from sklearn.metrics import mean_squared_error,r2_score#, ndcg_score
import math
from manifoldv2_mask_combination import manifold_learning
import os
import pandas as pd

os.chdir('/nas/longleaf/home/meisheng/FP_Project/code_running/ADRML/IF_code')

def normalise_sim(similarity_matrix):
       """ This function aims to normalize the similarity matrix
       using symmetric normalized Laplacian
       The input must be a square matrix
       """
       
       similarity_matrix = np.matrix(similarity_matrix)
       
       for round in range(200):
            summ = np.sum(similarity_matrix, axis=1)
            a = np.matrix(summ)
            D = np.diag(a.A1) # the diagonal matrix
            D1 = np.linalg.pinv(np.sqrt(D)); 
            similarity_matrix = D1 * similarity_matrix * D1;
    
       return similarity_matrix

# CHECK
def modelEvaluation(real_matrix,predict_matrix,testPosition): 
       """ This function computes the evaluation criteria
       
       real_matrix: is a matrix with cell lines in rows, drugs in columns,
       and real IC50 in its elemnts
       
       predict_matrix: has the same size as the real matrix 
       with the predicted IC50 values
       
       testPosition: is a vecoto, containing the pairs of (i,j) indices of 
       cell line-drug pairs that were considered as the test samples 
       in cross validation
       """
       
       real_pred=[]
       real_labels=[]
       predicted_probability=[]
        
       # gather the test position values in real_matrix and predict_matrix 
       # into vectors
       for i in range(0,len(testPosition)):
           real_pred.append(real_matrix[testPosition[i][0], testPosition[i][1]])
           predicted_probability.append(predict_matrix[testPosition[i][0], testPosition[i][1]])

       real_labels = np.array(real_labels)
       
       # computing evaluation criteria
       mse = mean_squared_error(real_pred, predicted_probability)
       rmse = math.sqrt(mse)
       R2 = r2_score(real_pred, predicted_probability)
       pearson = np.corrcoef(real_pred, predicted_probability)[0, 1]
       results = [mse, rmse, R2, pearson]

       return results


def runapp(response_file, simC_name, simD_name, fold_infos, percent, miu, landa, CV_num, repetition, output_name):
    
    """ This function runs the cross validtion
    
    response_file is the address and name of real IC50 file
    SimC_name and SimD_name are the address and name of similarity matrices
    of cell lines and drugs, respectively
    percent is the rank of latent matrix
    miu and landa are two model hyperparametrs.
    miu is the regularization coeeficient for latent matrices
    landa controls the similarity conservation while manifold learning
    CV_num is the number of folds in cross validation
    repetition is the number of repeting the cross validation
    """
    #-----------------------------------------------------------
    
    
    #reading IC50 file
    R = np.loadtxt(response_file, dtype=float, delimiter=",", skiprows=1, usecols=range(1,21))
    R_info = pd.read_csv(response_file, index_col=0)
        
    # reading similarity matrices
    simD = np.loadtxt(simD_name, dtype=float, delimiter=",", skiprows=1, usecols=range(1,21))
    simC = np.loadtxt(simC_name, dtype=float, delimiter=",", skiprows=1, usecols=range(1,101))

    # Read in the identical fold information
    folds_info = pd.read_csv(fold_infos)

    # initilaizing evaluation criteria
    mse = []
    rmse = []
    R2 = []
    pear = []
    
    # repeting the cross valiation
    for rep in range(repetition):
        print('*********repetition:' + str(rep+1) + "**********\n")
        #running the cross validation
        for CV in range(0, CV_num):
            n_fold = CV + 1
            print('*********round:' + str(n_fold) + "**********\n")
            # seleting test positions
            fold_valid = folds_info[(folds_info['fold number'] == n_fold) & (folds_info['indicator'] == 'valid')]

            # Step 2: Get the 'pair information' from the filtered rows and convert to a set
            fold_valid_pairs = set(fold_valid['pair information'])

            # Step 3 and 4: Find the matching positions in the 'R_info' DataFrame
            matching_positions = []
            for sanger_id in R_info.index:
                for drug_name in R_info.columns:
                    pair_info = f'{drug_name}/{sanger_id}'
                    if pair_info in fold_valid_pairs:
                        row_index = R_info.index.get_loc(sanger_id)
                        col_index = R_info.columns.get_loc(drug_name)
                        matching_positions.append((row_index, col_index))

            testPosition = np.array(matching_positions)

            train_IC = copy.deepcopy(R)
            
            # set the IC50 values for the test positions to zero
            for i in range(0, len(testPosition)):
                train_IC[testPosition[i, 0], testPosition[i, 1]] = 0
            testPosition = list(testPosition)
          
            # initialize the latent matrices
            
            N = len(train_IC)
            M = len(train_IC[0])
            dim = min(N, M)
            K  = int(round (percent * dim))
            P = np.random.rand(N ,K)
            Q = np.random.rand(M, K)
            
            # call manifold learning
            predict_matrix1, A1,B1 = manifold_learning(train_IC, P, Q, K, simD, simC, landa, miu)
            predict_matrix2, A2,B2 = manifold_learning(train_IC.T, B1, A1, K, simC, simD, landa, miu)
            predict_matrix = 0.5 * (predict_matrix1 + predict_matrix2.T)
            
            # evaluate the model #CHECK
            results = modelEvaluation(R, predict_matrix, testPosition)
            mse.append(results[0])
            rmse.append(results[1])
            R2.append(results[2])
            pear.append(results[3])

    print(mse, ' ', rmse, ' ', R2, ' ', pear, ' ')
    dict = {'mse': mse, 'rmse': rmse, 'R2': R2, 'PCC': pear}
    df = pd.DataFrame(dict)
    df.to_csv(output_name)

def main():
    # get the options from user
    for arg in sys.argv[1:]:
      (key,val) = arg.rstrip().split('=')
      if key == 'response_dirc':
          response_file=val
      elif key=='simC_dirc':
          simC_name=val
      elif key=='simD_dirc':
          simD_name=val
      elif key=='fold_info':
          fold_infos=val
      elif key=='dim':
          percent=float(val)
      elif key=='miu':
          miu=float(val)
      elif key=='lambda':
          landa=float(val)
      elif key=='CV':
          CV_num=int(val)
      elif key=='repeat':
          repetition=int(val)
      elif key=='out_name':
          output_name=val
          
    # call the method
    runapp(response_file, simC_name, simD_name, fold_infos, percent, miu, landa, CV_num, repetition, output_name)
    
main()
