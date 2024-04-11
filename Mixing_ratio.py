import pandas as pd
import numpy as np
import os
import scipy.stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

#######################################
### calculating relative expression ###
#######################################

def relative_expression(expr_data):
  expr_data['mean'] = expr_data.iloc[:, 10:13].mean(axis=1)
  expr_data.iloc[:, 1:10] = (expr_data.iloc[:, 1:10].sub(expr_data['mean'], axis=0)
  expr_data['m8'] = expr_data.iloc[:, 1:4].mean(axis=1)
  expr_data['t1'] = expr_data.iloc[:, 4:7].mean(axis=1)
  expr_data['t2'] = expr_data.iloc[:, 7:10].mean(axis=1)
  data = expr_data[["ensembl_id", 'm8', 't1' , 't2']]
  return data

#######################################
### extract data for each sample    ###
#######################################

def ratio_accuracy(expr_data):
  expr_data.iloc[:, 1:25] = np.log2(expr_data.iloc[:, 1:25] + 0.01)
  expr_sample = expr_data[["ensembl_id","M8_1", "M8_2", "M8_3","T1_1", "T1_2", "T1_3","T2_1", "T2_2", "T2_3","D6_1", "D6_2","D6_3"]]
  expr_relative = relative_expression(expr_sample)
  expr_relative = expr_relative.reset_index(drop=True)
  expr_relative['expected_T1'] = expr_relative['m8'].apply(lambda x: np.log2(0.245182377 + 0.745119679 * 2**x))
  expr_relative['expected_T2'] = expr_relative['m8'].apply(lambda x: np.log2(0.739998388 + 0.240258211 * 2**x))
  rmse1 = np.sqrt(((expr_relative['t1'] -  expr_relative['expected_T1']) ** 2).mean())
  rmse2 = np.sqrt(((expr_relative['t2'] -  expr_relative['expected_T2']) ** 2).mean())
  rmse_list = [rmse1,rmse2]
  return rmse_list


###   start ###

### 1. read a expression data for all 24 samples, tab separated   ###
###    Column names must be "ensembl_id","A_1","A_2","A_3","B_1","B_2","B_3","M8_1", "M8_2", "M8_3", "F7_1", "F7_2", "F7_3","D5_1", "D5_2","D5_3", "T1_1", "T1_2","T1_3","T2_1", "T2_2", "T2_3", "D6_1", "D6_2","D6_3" ###
expr_data = pd.read_table("/your path/", header=0, sep = '\t',index_col=False)

### 2. calculating RMSE between expected relative expression and observed relative expression in T1/D6 or T2/D6 ###
RMSE = ratio_accuracy(expr_data)  
