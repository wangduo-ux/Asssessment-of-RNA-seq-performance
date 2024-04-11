import pandas as pd
import numpy as np
from openpyxl import Workbook

#######################################
### calculating relative expression ###
#######################################

def relative_expression(expr_data,lab):
  expr_data = expr_data.assign(mean=expr_data.iloc[:, 1:7].mean(axis=1))
  expr_data['mean'] = expr_data.iloc[:, 4:7].mean(axis=1)
  expr_data.iloc[:, 1:4] = (expr_data.iloc[:, 1:4].sub(expr_data.iloc[:, -1], axis=0)
  expr_data['mean'] = expr_data.iloc[:, 1:4].mean(axis=1)
  data = expr_data.loc[:, ["ensembl_id", "mean"]]
  return data

#######################################
### extract data for each sample    ###
#######################################

def processdata(expr_data):
  #data1 = df[["Gene ID",f"{lab}202201", f"{lab}202202",f"{lab}202203",f"{lab}202204", f"{lab}202205",f"{lab}202206"]]
  expr_data.iloc[:, 1:25] = np.log2(expr_data.iloc[:, 1:25] + 0.01)
  A_B = expr_data[["ensembl_id","A_1", "A_2", "A_3","B_1", "B_2","B_3"]]
  m8_d6 = expr_data[["ensembl_id","M8_1", "M8_2", "M8_3","D6_1", "D6_2","D6_3"]]
  f7_d6 = expr_data[["ensembl_id","D5_1", "D5_2", "D5_3","D6_1", "D6_2","D6_3"]]
  d5_d6 = expr_data[["ensembl_id","F7_1", "F7_2", "F7_3","D6_1", "D6_2","D6_3"]]
  m8_d6 = relative_expression(m8_d6)
  f7_d6 = relative_expression(f7_d6)
  d5_d6 = relative_expression(d5_d6)
  m8_d6['ensembl_id'] = 'M8/D6_' + m8_d6['ensembl_id'].astype(str)
  f7_d6['ensembl_id'] = 'F7/D6_' + f7_d6['ensembl_id'].astype(str)
  d5_d6['ensembl_id'] = 'D5/D6_' + d5_d6['ensembl_id'].astype(str)
  combine = pd.concat([m8_d6, f7_d6, d5_d6,A_B], axis=0)
  combined_df = combined_df.rename(columns={'mean': 'FC'})
  return combined_df

#######################################
###         calculating rmse        ###
#######################################

def accuracy(expr_data):
  ratio_expr = processdata(expr_data)
  combine_1 = pd.merge(truth,ratio_expr,on="ensembl_id",how = "inner")
  combine_2 = pd.merge(truth_2,ratio_expr,on="ensembl_id",how = "inner")
  combine_3 = pd.merge(truth_3,ratio_expr,on="ensembl_id",how = "inner")
  rmse = np.sqrt(((combine['log2FC'] - combine['FC']) ** 2).mean())
  rmse_2 = np.sqrt(((combine_2['log2FC'] - combine_2['FC']) ** 2).mean())
  rmse_3 = np.sqrt(((combine_3['log2FC'] - combine_3['FC']) ** 2).mean())
  rmse_list = [rmse,rmse_2,rmse_3]
  return rmse_list




### start ###

 
### 1. read a expression data for all 24 samples, tab separated   ###
###    Column names must be "ensembl_id","A_1","A_2","A_3","B_1","B_2","B_3","M8_1", "M8_2", "M8_3", "F7_1", "F7_2", "F7_3","D5_1", "D5_2","D5_3", "T1_1", "T1_2","T1_3","T2_1", "T2_2", "T2_3", "D6_1", "D6_2","D6_3" ###
expr_data = pd.read_table("/your path/", header=0, sep = '\t',index_col=False)

### 2. Read the Quartet reference datasets ###

truth = pd.read_table("/your path/quartet_reference_datsets.txt", header=0, sep = '\t',index_col=False)

### 3. Read the  Quartet TaqMan datasets   ###

truth_2 = pd.read_table("/your path/quartet_taqman_datasets.txt", header=0, sep = '\t',index_col=False)

### 4. Read the MAQC TaqMan datasets       ###

truth_3 = pd.read_table("/your path/maqc_taqman_datasets.txt", header=0, sep = '\t',index_col=False)

### 5. calculating RMSE between your data and reference datasets ###
rmse = accuracy(expr_data)
