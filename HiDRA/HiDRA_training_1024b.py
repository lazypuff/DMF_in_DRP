"""
Editted by Meisheng Xiao Aug 3 2023.
Training new HiDRA model
Requirement:
    expression.csv: Expression of all genes for all cell lines
    geneset.gmt: Gene set description file. Consists of Gene set name, Source, Gene set member 1, Gene set member 2, ...
    Training.csv: Training pair list. Consists of idx, Drug name, Sanger ID, IC50 value for that pair.
    Validation.csv: Validation pair list. Consists of idx, Drug name, Sanger ID, IC50 value for that pair.
    input_dir: The directory that includes input files.
"""

#Import basic packages
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
import seaborn as sns

import os
import argparse

#Import keras modules
import keras
import tensorflow.compat.v1 as tf
import keras.backend as K
import keras.layers
from keras.layers import Layer 
import keras.initializers
from keras.models import Model, Sequential,load_model
from keras.layers import Input, Dense, Dropout, BatchNormalization, Activation, Multiply, multiply,dot
from keras.layers import Concatenate,concatenate
from keras.optimizers import Adam
from keras.utils import plot_model

#Import rdkit
import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols

#Fix the random seed
np.random.seed(5)


#Get setssion to control GPU memory usage (only 20% of GPU memory will be used)
def get_session(gpu_fraction=0.1):

    num_threads = 1
    gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=gpu_fraction)

    if num_threads:
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options, intra_op_parallelism_threads=num_threads))
    else:
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))

#Read Input file for given pair list
def read_files(GDSC_df,Input_directory):
    #Read Gene expression file
    GeneExpression_with_Symbol=pd.read_csv("expression.csv")
    GeneExpression_with_Symbol.index=GeneExpression_with_Symbol.Gene_Symbol

    #Read Gene set file
    GeneSet_List=[]
    with open("geneset.gmt") as f:
        reader = csv.reader(f)
        data = list(list(rec) for rec in csv.reader(f, delimiter='\t')) #reads csv into a list of lists
        for row in data:
            GeneSet_List.append(row)

    #Gene set processing
    GeneSet_Dic={}
    for GeneSet in GeneSet_List:
        GeneSet_Dic[GeneSet[0]]=GeneSet[2:]
        
    GeneSet_Dic_withoutNA={}
    for GeneSet in GeneSet_Dic:
        inter1 = GeneSet_Dic[GeneSet]
        newList = []
        for element in inter1:
            if element in GeneExpression_with_Symbol.index:
                newList.append(element)
        GeneSet_Dic_withoutNA[GeneSet]=GeneExpression_with_Symbol['SIDM00003'][newList].dropna().index.values
        
    #Read GDSC's IC50 information
    GDSC_without_duplicated=GDSC_df
    GDSC_without_duplicated.columns=['Drug name','Sanger ID']

    #Read input files for HiDRA model
    X_origin=[pd.read_csv(Input_directory+'/'+str(x)+'.csv',index_col=0) for x in range(len(GeneSet_Dic_withoutNA))] #Read all inputs for cell lines
    X_origin=[df.loc[GDSC_without_duplicated['Sanger ID']] for df in X_origin] #Pick up data for cell lines in given pairs
    Drug=pd.read_csv(Input_directory+'/'+'drug.csv',index_col=0) #Read drug input
    X_origin.append(Drug.loc[GDSC_without_duplicated['Drug name']]) #Pick up data for drugs in given pairs

    return X_origin


#Make new HiDRA model
def Making_Model():
    #Read Gene expression and Gene set for model making
    #They are same with the source codes in the function 'Read_feils'
    #Read Gene expression file
    GeneExpression_with_Symbol=pd.read_csv("expression.csv")
    GeneExpression_with_Symbol.index=GeneExpression_with_Symbol.Gene_Symbol


    #Read Gene set file
    GeneSet_List=[]
    with open("geneset.gmt") as f:
        reader = csv.reader(f)
        data = list(list(rec) for rec in csv.reader(f, delimiter='\t')) #reads csv into a list of lists
        for row in data:
            GeneSet_List.append(row)

    GeneSet_Dic={}
    for GeneSet in GeneSet_List:
        GeneSet_Dic[GeneSet[0]]=GeneSet[2:]

    GeneSet_Dic_withoutNA={}
    for GeneSet in GeneSet_Dic:
        inter1 = GeneSet_Dic[GeneSet]
        newList = []
        for element in inter1:
            if element in GeneExpression_with_Symbol.index:
                newList.append(element)
        GeneSet_Dic_withoutNA[GeneSet]=GeneExpression_with_Symbol['SIDM00003'][newList].dropna().index.values
        
    
    #HiDRA model with keras
    #Drug-level network
    Drug_feature_length=1024
    Drug_Input=Input((Drug_feature_length,), dtype='float32', name='Drug_Input')
    
    Drug_Dense1=Dense(256, name='Drug_Dense_1')(Drug_Input)
    Drug_Dense1=BatchNormalization(name='Drug_Batch_1')(Drug_Dense1)
    Drug_Dense1=Activation('relu', name='Drug_RELU_1')(Drug_Dense1)

    Drug_Dense2=Dense(128, name='Drug_Dense_2')(Drug_Dense1)
    Drug_Dense2=BatchNormalization(name='Drug_Batch_2')(Drug_Dense2)
    Drug_Dense2=Activation('relu', name='Drug_RELU_2')(Drug_Dense2)
    
    #Drug network that will be used to attention network in the Gene-level network and Pathway-level network
    Drug_Dense_New1=Dense(128, name='Drug_Dense_New1')(Drug_Input)
    Drug_Dense_New1=BatchNormalization(name='Drug_Batch_New1')(Drug_Dense_New1)
    Drug_Dense_New1=Activation('relu', name='Drug_RELU_New1')(Drug_Dense_New1)

    Drug_Dense_New2=Dense(32, name='Drug_Dense_New2')(Drug_Dense_New1)
    Drug_Dense_New2=BatchNormalization(name='Drug_Batch_New2')(Drug_Dense_New2)
    Drug_Dense_New2=Activation('relu', name='Drug_RELU_New2')(Drug_Dense_New2)

    #Gene-level network
    GeneSet_Model=[]
    GeneSet_Input=[]
    
    #Making networks whose number of node is same with the number of member gene in each pathway    
    for GeneSet in GeneSet_Dic_withoutNA.keys():
        Gene_Input=Input(shape=(len(GeneSet_Dic_withoutNA[GeneSet]),),dtype='float32', name=GeneSet+'_Input')
        Drug_effected_Model_for_Attention=[Gene_Input]
        #Drug also affects to the Gene-level network attention mechanism
        Drug_Dense_Geneset=Dense(int(len(GeneSet_Dic_withoutNA[GeneSet])/4)+1,dtype='float32',name=GeneSet+'_Drug')(Drug_Dense_New2)
        Drug_Dense_Geneset=BatchNormalization(name=GeneSet+'_Drug_Batch')(Drug_Dense_Geneset)
        Drug_Dense_Geneset=Activation('relu', name=GeneSet+'Drug_RELU')(Drug_Dense_Geneset)
        Drug_effected_Model_for_Attention.append(Drug_Dense_Geneset) #Drug feature to attention layer
 
        Gene_Concat=concatenate(Drug_effected_Model_for_Attention,axis=1,name=GeneSet+'_Concat')
        #Gene-level attention network
        Gene_Attention = Dense(len(GeneSet_Dic_withoutNA[GeneSet]), activation='tanh', name=GeneSet+'_Attention_Dense')(Gene_Concat)
        Gene_Attention=Activation(activation='softmax', name=GeneSet+'_Attention_Softmax')(Gene_Attention)
        Attention_Dot=dot([Gene_Input,Gene_Attention],axes=1,name=GeneSet+'_Dot')
        Attention_Dot=BatchNormalization(name=GeneSet+'_BatchNormalized')(Attention_Dot)
        Attention_Dot=Activation('relu',name=GeneSet+'_RELU')(Attention_Dot)
        
	#Append the list of Gene-level network (attach new pathway)
        GeneSet_Model.append(Attention_Dot)
        GeneSet_Input.append(Gene_Input)

    Drug_effected_Model_for_Attention=GeneSet_Model.copy()
    
    #Pathway-level network
    Drug_Dense_Sample=Dense(int(len(GeneSet_Dic_withoutNA)/16)+1,dtype='float32',name='Sample_Drug_Dense')(Drug_Dense_New2)
    Drug_Dense_Sample=BatchNormalization(name=GeneSet+'Sample_Drug_Batch')(Drug_Dense_Sample)
    Drug_Dense_Sample=Activation('relu', name='Sample_Drug_ReLU')(Drug_Dense_Sample)    #Drug feature to attention layer
    Drug_effected_Model_for_Attention.append(Drug_Dense_Sample)
    GeneSet_Concat=concatenate(GeneSet_Model,axis=1, name='GeneSet_Concatenate')
    Drug_effected_Concat=concatenate(Drug_effected_Model_for_Attention,axis=1, name='Drug_effected_Concatenate')
    #Pathway-level attention
    Sample_Attention=Dense(len(GeneSet_Dic_withoutNA.keys()),activation='tanh', name='Sample_Attention_Dense')(Drug_effected_Concat)
    Sample_Attention=Activation(activation='softmax', name='Sample_Attention_Softmax')(Sample_Attention)
    Sample_Multiplied=multiply([GeneSet_Concat,Sample_Attention], name='Sample_Attention_Multiplied')
    Sample_Multiplied=BatchNormalization(name='Sample_Attention_BatchNormalized')(Sample_Multiplied)
    Sample_Multiplied=Activation('relu',name='Sample_Attention_Relu')(Sample_Multiplied)
    
    #Making input list
    Input_for_model=[]
    for GeneSet_f in GeneSet_Input:
        Input_for_model.append(GeneSet_f)
    Input_for_model.append(Drug_Input)
    
    #Concatenate two networks: Pathway-level network, Drug-level network 
    Total_model=[Sample_Multiplied,Drug_Dense2]
    Model_Concat=concatenate(Total_model,axis=1, name='Total_Concatenate')

    #Response prediction network
    Concated=Dense(128, name='Total_Dense')(Model_Concat)
    Concated=BatchNormalization(name='Total_BatchNormalized')(Concated)
    Concated=Activation(activation='relu', name='Total_RELU')(Concated)

    Final=Dense(1, name='Output')(Concated)
    model=Model(inputs=Input_for_model,outputs=Final)
    
    return model


def main():
    K.set_session(get_session())

   
    #Reading argument 
    parser=argparse.ArgumentParser(description='HiDRA:Hierarchical Network for Drug Response Prediction with Attention-Training')
    
    #Options
    parser.add_argument('-t',type=str,help='The training pair list')
    parser.add_argument('-v',type=str,help='The validation pair list')
    parser.add_argument('-i',type=str,help='The input directory')
    parser.add_argument('-e',type=str,help='The epoch in the training process')
    parser.add_argument('-o',type=str,help='The output path that model file be stored')
    
    args=parser.parse_args()
    
    #Read Training and Validation files
    training_pair=pd.read_csv(args.t)
    validation_pair=pd.read_csv(args.v)
    training_label=training_pair[['IC50']]
    validation_label=validation_pair[['IC50']]
    training_df=training_pair[['Drug name','Sanger ID']]
    validation_df=validation_pair[['Drug name','Sanger ID']]
    training_input=read_files(training_df,args.i)
    validation_input=read_files(validation_df,args.i)
    
    #Training
    epoch=int(args.e)
    model=Making_Model()
    model.compile(loss='mean_squared_error',optimizer='adam')
    hist=model.fit(training_input, training_label,shuffle=True,epochs=epoch,batch_size=1024,verbose=1,validation_data=(validation_input,validation_label))
    model.save(args.o) #Save the model to the output directory
    
    
    
if __name__=="__main__":
    main()
