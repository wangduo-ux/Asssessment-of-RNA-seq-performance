import scipy.stats
import os
import pandas as pd
import numpy as np

##################################################################################################################################
### calculating pearson correoaltion coefficient between ERCC concentrations and absolute expression of ERCC gene in your data ###
##################################################################################################################################

def ERCC_accuracy(expr_data, ercc_control):
    ERCC_data = expr_data[expr_data['ensembl_id'].str.contains("SPIKEIN")].copy()
    ERCC_data['mean'] = ERCC_data[["M8_1", "M8_2", "M8_3"]].mean(axis=1)
    ERCC_data['mean_2'] = ERCC_data[["D6_1", "D6_2", "D6_3"]].mean(axis=1)
    ERCC_data_M8 = ERCC_data[["ensembl_id", 'mean']].copy()
    ERCC_data_D6 = ERCC_data[["ensembl_id", 'mean_2']].copy()
    ERCC_data_M8['ensembl_id'] = "M8_" + ERCC_data_M8['ensembl_id']
    ERCC_data_D6['ensembl_id'] = "D6_" + ERCC_data_D6['ensembl_id']
    ERCC_data_D6 = ERCC_data_D6.rename(columns={'mean_2': 'mean'})
    ERCC_expr = pd.concat([ERCC_data_M8, ERCC_data_D6], axis=0)
    ERCC_expr['mean'] = np.log2(ERCC_expr['mean'] + 0.01)
    
    ercc_control = ercc_control.rename(columns={'SPIKEIN_ID': 'ensembl_id'})
    mix1 = ercc_control[['SPIKEIN_ID', 'Mix_1']].copy()
    mix1['ensembl_id'] = "M8_" + mix1['SPIKEIN_ID']
    mix2 = ercc_control[['SPIKEIN_ID', 'Mix_2']].copy()
    mix2['ensembl_id'] = "D6_" + mix2['SPIKEIN_ID']
    mix2 = mix2.rename(columns={'Mix_2': 'Mix_1'})
    ercc_control_expr = pd.concat([mix1, mix2], axis=0)
    ercc_control_expr['Mix_1'] = np.log2(ercc_control_expr['Mix_1'] + 0.01)
    
    combine = pd.merge(ERCC_expr, ercc_control_expr, on="ensembl_id", how="inner")
    pearson, p = scipy.stats.pearsonr(combine['Mix_1'], combine['mean'])
    return pearson
      
### start ###

### 1. read a expression data for all 24 samples, tab separated   ###
###    Column names must be "ensembl_id","A_1","A_2","A_3","B_1","B_2","B_3","M8_1", "M8_2", "M8_3", "F7_1", "F7_2", "F7_3","D5_1", "D5_2","D5_3", "T1_1", "T1_2","T1_3","T2_1", "T2_2", "T2_3", "D6_1", "D6_2","D6_3" ###
expr_data = pd.read_table("/your path/", header=0, sep = '\t',index_col=False)


### 2. read the ERCC control data   ###
ercc_control = pd.read_table('/your path/ercc_control_data.txt',index_col=False,header = 0)


### 3. calculating rmse between expected ERCC ratio and relative expression of ERCC gene in your data ###
pearson = ERCC_accuracy(ERCC_data, ercc_data)
    
     
