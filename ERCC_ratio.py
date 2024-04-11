import pandas as pd
import numpy as np
import os


############################################
### calculating relative gene expression ###
############################################

def ratio(df):
  df.iloc[:, 1:7] = np.log2(df.iloc[:, 1:7] + 0.01)
  df['mean'] = df["D6_1", "D6_2","D6_3"]].mean(axis=1)
  df.iloc[:, 1:7] = df.iloc[:, 1:7].sub(df['mean'], axis=0)
  df['M8'] = df[["M8_1","M8_2","M8_3"]].mean(axis=1)
  M8 = df.loc[:, ["ensembl_id", "M8"]]
  return M8


######################################################################################################
### calculating rmse between expected ERCC ratio and relative expression of ERCC gene in your data ###
######################################################################################################

def ERCC_accuracy(expr_data,ercc_control):
  expr_data = expr_data[["ensembl_id","M8_1","M8_2","M8_3","D6_1", "D6_2","D6_3"]]
  relative_data = ratio(expr_data)
  combine = pd.merge(ercc_control,relative_data,on="ensembl_id",how = "inner")
  rmse = np.sqrt(((combine['M8']-combine['expected']) ** 2).mean())
  return rmse  

### start ###

### 1. read a expression data for all 24 samples, tab separated   ###
###    Column names must be "ensembl_id","A_1","A_2","A_3","B_1","B_2","B_3","M8_1", "M8_2", "M8_3", "F7_1", "F7_2", "F7_3","D5_1", "D5_2","D5_3", "T1_1", "T1_2","T1_3","T2_1", "T2_2", "T2_3", "D6_1", "D6_2","D6_3" ###
expr_data = pd.read_table("/your path/", header=0, sep = '\t',index_col=False)
ERCC_data = expr_data[expr_data['ensembl_id'].str.contains("SPIKEIN")]

### 2. read the ERCC control data   ###
ercc_control = pd.read_table('/your path/ercc_control_data.txt',index_col=False,header = 0)
ercc_control = ercc_control.rename(columns={'SPIKEIN_ID': 'ensembl_id'})
ercc_control['expected'] = np.log2(ercc_control['expected'])

### 3. calculating rmse between expected ERCC ratio and relative expression of ERCC gene in your data ###
rmse = ERCC_accuracy(ERCC_data, ercc_data)
