import pandas as pd
import numpy as np
from openpyxl import Workbook
from scipy.stats import spearmanr


#######################################################
### obtain gene abundance from maqc TaqMan datasets ###
#######################################################

def absolute_truth_2(truth):
  A = truth.loc[:, ["ensembl_id", "Sample_A_Mean"]]
  B = truth.loc[:, ["ensembl_id", "Sample_B_Mean"]]
  A = A.rename(columns={'Sample_A_Mean': 'abundance'})
  B = B.rename(columns={'Sample_B_Mean': 'abundance'})
  A['ensembl_id'] = 'A_' + A['ensembl_id'].astype(str)
  B['ensembl_id'] = 'B_' + B['ensembl_id'].astype(str)
  truth = pd.concat([A, B], axis=0)
  truth['abundance'] = np.log2(truth['abundance']+ 0.01)
  return truth

##########################################################
### obtain gene abundance from quartet TaqMan datasets ###
##########################################################

def absolute_truth_1(truth):
  test = truth.loc[:, ["ensembl_id", "test_mean"]]
  D6 = truth.loc[:, ["ensembl_id", "D6_mean"]]
  test = test.rename(columns={'test_mean': 'abundance'})
  D6 = D6.rename(columns={'D6_mean': 'abundance'})
  test['ensembl_id'] = test['Sample'] + "_" + test['ensembl_id'].astype(str)
  D6['ensembl_id'] = 'D6_' + D6['ensembl_id'].astype(str)
  truth = pd.concat([test, D6], axis=0)
  truth['abundance'] = np.log2(truth['abundance']+0.01)
  return truth

#########################################################################################
### calculating correlation coefficient between your data and two reference datasets ###
#########################################################################################

def absolute_accuracy(expr_data,truth_1,truth_2):
  data = {}
  for prefix in ["A", "M8", "D5", "F7", "D6"]:
    data[f"{prefix}"] = expr_data[["ensembl_id",f"{prefix}_1",f"{prefix}_2",f"{prefix}_3"]]
    data[f"{prefix}"].loc[:, 'mean'] = data[f"{prefix}"].iloc[:, 1:4].mean(axis=1)
    data[f"{prefix}"] = data[f"{prefix}"][['ensembl_id','mean']]
    data[f"{prefix}"].loc[:, 'ensembl_id'] = prefix + "_" + data[f"{prefix}"].loc[:, 'ensembl_id']
  expr_all = pd.concat(data, axis=0) 
  expr_all = expr_all.reset_index(drop=True)
  expr_all['mean'] = np.log2(expr_all['mean']+0.01)
  quartet_reference = absolute_truth_1(truth_1)
  maqc_reference = absolute_truth_2(truth_2)
  combine = pd.merge(quartet_reference,expr_all,on="ensembl_id",how = "inner")
  cor = combine['abundance'].corr(combine['mean'])
  combine_2 = pd.merge(maqc_reference,expr_all,on="ensembl_id",how = "inner")
  cor_2 = combine_2['abundance'].corr(combine_2['mean'])
  accuracy_list = [cor, cor_2]
  return accuracy_list


###   start ###

### 1. read a expression data for all 24 samples, tab separated   ###
###    Column names must be "ensembl_id","A_1","A_2","A_3","B_1","B_2","B_3","M8_1", "M8_2", "M8_3", "F7_1", "F7_2", "F7_3","D5_1", "D5_2","D5_3", "T1_1", "T1_2","T1_3","T2_1", "T2_2", "T2_3", "D6_1", "D6_2","D6_3" ###
expr_data = pd.read_table("/your path/", header=0, sep = '\t',index_col=False)

### 2. Read the Quartet TaqMan datasets ###

truth_1 = pd.read_table("/your path/quartet_taqman_datasets.txt", header=0, sep = '\t',index_col=False)

### 3. Read the MAQC TaqMan datasets    ###

truth_2 = pd.read_table("/your path/maqc_taqman_datasets.txt", header=0, sep = '\t',index_col=False)

### 4. calculating Pearson correlation coefficients between your data and reference datasets ###
Pearson_correlation = absolute_accuracy(expr_data)  
