"""
Predict with the HiDRA model

Requirement:
    model.hdf: Pre-trained model file
    predict.csv: The list of cancer-drug pairs that will be predicted. Consists of idx (int), Drug name, Sanger ID.
    input_dir: The directory that includes input files
"""
#Import basic packages
import numpy as np
import pandas as pd
import csv

import os
import argparse

#Import keras modules
import tensorflow.compat.v1 as tf
import keras.backend as K
import keras
import keras.layers
from keras.layers import Layer 
import keras.initializers
from keras.models import Model, Sequential,load_model
from keras.layers import Input, Dense, Dropout, BatchNormalization, Activation, Multiply, multiply,dot
from keras.layers import Concatenate,concatenate
from keras.optimizers import Adam
from keras.utils import plot_model

#Fix the random seed
np.random.seed(5)

#Using 20% of GPU memory only
def get_session(gpu_fraction=0.1):

    num_threads = 1
    gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=gpu_fraction)

    if num_threads:
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options, intra_op_parallelism_threads=num_threads))
    else:
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))

def read_files(Predict_list,Input_directory):
    
    #Read predict list
    GDSC_without_duplicated=pd.read_csv(Predict_list,index_col=0)
    GDSC_without_duplicated.columns=['Drug name','Sanger ID']
    
    #Read input files for HiDRA model
    X_origin=[pd.read_csv(Input_directory+'/'+str(x)+'.csv',index_col=0) for x in range(186)]
    X_origin=[df.loc[GDSC_without_duplicated['Sanger ID']] for df in X_origin]
    Drug=pd.read_csv(Input_directory+'/'+'drug.csv',index_col=0)
    X_origin.append(Drug.loc[GDSC_without_duplicated['Drug name']])
    
    return X_origin


def main():
    K.set_session(get_session())

    #Reading argument 
    parser=argparse.ArgumentParser(description='HiDRA:Hierarchical Network for Drug Response Prediction with Attention-Predict')
    
    #Options
    parser.add_argument('-m',type=str,help='The model file')
    parser.add_argument('-p',type=str,help='The prediction list')
    parser.add_argument('-i',type=str,help='The input directory')
    parser.add_argument('-o',type=str,help='The output file path that prediction result be stored')
    
    args=parser.parse_args()
   
    #Read input files from predict list 
    model_input=read_files(args.p,args.i)
    #Load model in hdf5 file format
    model=load_model(args.m,compile=False)
    #Predict
    result=model.predict(model_input)
    result=[y[0] for y in result]
    predict_list=pd.read_csv(args.p)
    predict_list['result']=result
    #Save the predict results to the output directory
    predict_list.to_csv(args.o)
    
if __name__=="__main__":
    main()
