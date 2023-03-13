
# -*- coding: utf-8 -*-
"""
"""
import itertools
from itertools import combinations
from scipy import stats
from scipy.ndimage.filters import gaussian_filter1d
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

################################################################################################################################
## FIGURE 3 showing average spatial population synchrony on te breeding ground for nonmigrants vs migrants under 
# different environmental autocorrelation simulations. 
################################################################################################################################

# NO Mig Route: 
abs_output_0 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0_Migroutes0_Nit1000_Carryover0.csv', header = None)
distances_1 = pd.read_csv('Midpoint_Distances_Gridsize150_Thin2.csv', header = None)
abs_output_0 = abs_output_0.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
distances_1 = distances_1.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
distlist_1 = distances_1.values.tolist()
flatdist_1 = (np.array(distlist_1)).flatten()
correlations_0 = {}
columns_0= abs_output_0.columns.tolist()
for col_a, col_b in itertools.combinations(columns_0, 2):
     correlations_0[col_a, '__' , col_b] = stats.pearsonr(abs_output_0.loc[:, col_a], abs_output_0.loc[:, col_b])
result_0 = pd.DataFrame.from_dict(correlations_0, orient='index')
result_0.columns = ['PCC', 'p-value']
result_0['distances'] = np.array(flatdist_1)
average_corr_0k = result_0.groupby('distances')['PCC'].mean()
average_corr_0k.drop(index=average_corr_0k.index[0], axis=0, inplace=True)


# 1 Mig Route: K
abs_output_1mig = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr1_Migroutes1_Nit1000_Carryover1.csv', header = None)
abs_output_1mig = abs_output_1mig.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_1mig = {}
columns_1mig = abs_output_1mig.columns.tolist()
for col_a, col_b in itertools.combinations(columns_1mig, 2):
     correlations_1mig[col_a, '__' , col_b] = stats.pearsonr(abs_output_1mig.loc[:, col_a], abs_output_1mig.loc[:, col_b])
result_1mig = pd.DataFrame.from_dict(correlations_1mig, orient='index')
result_1mig.columns = ['PCC', 'p-value']
result_1mig['distances'] = np.array(flatdist_1)
average_corr_1migk = result_1mig.groupby('distances')['PCC'].mean()
average_corr_1migk.drop(index=average_corr_1migk.index[0], axis=0, inplace=True)



###### TWO MIGRATION ROUTES ################################################################################################
# 2 Mig Route: 0 correlation, r
#abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.5_F2.0_Corr0_Migroutes2_Nit1000_Carryover1.csv', header = None)
#abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_3 = {}
#columns_3 = abs_output_3.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_3, 2):
#     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
#result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
#result_3.columns = ['PCC', 'p-value']
#result_3['distances'] = np.array(flatdist_1)
#average_corr_r20 = result_3.groupby('distances')['PCC'].mean()
#average_corr_r20.drop(index=average_corr_r20.index[0], axis=0, inplace=True)
#average_corr_smoothed_k3 = gaussian_filter1d(average_corr_k3, sigma=2)


# 2 Mig Route: 0 correlation, K
#abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0_Migroutes2_Nit1000_Carryover1.csv', header = None)
#abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_3 = {}
#columns_3 = abs_output_3.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_3, 2):
#     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
#result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
#result_3.columns = ['PCC', 'p-value']
#result_3['distances'] = np.array(flatdist_1)
#average_corr_k20 = result_3.groupby('distances')['PCC'].mean()
#average_corr_k20.drop(index=average_corr_k20.index[0], axis=0, inplace=True)
#average_corr_smoothed_k3 = gaussian_filter1d(average_corr_k3, sigma=2)


# 2 Mig Route: 0.25 correlation,  r
#abs_output_4 = pd.read_csv('CorrWithin0.75/Output_S0.5_F2.0_Corr0.25_Migroutes2_Nit1000_Carryover1.csv', header = None)
#abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_4 = {}
#columns_4 = abs_output_4.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_4, 2):
#     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
#result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
#result_4.columns = ['PCC', 'p-value']
#result_4['distances'] = np.array(flatdist_1)
#average_corr_r2025 = result_4.groupby('distances')['PCC'].mean()
#average_corr_r2025.drop(index=average_corr_r2025.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.25 correlation, K
#abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.25_Migroutes2_Nit1000_Carryover1.csv', header = None)
#abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_3 = {}
#columns_3 = abs_output_3.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_3, 2):
#     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
#result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
#result_3.columns = ['PCC', 'p-value']
#result_3['distances'] = np.array(flatdist_1)
#average_corr_k2025 = result_3.groupby('distances')['PCC'].mean()
#average_corr_k2025.drop(index=average_corr_k2025.index[0], axis=0, inplace=True)
#average_corr_smoothed_k3 = gaussian_filter1d(average_corr_k3, sigma=2)


# 2 Mig Route: 0.75 correlation, r 
#abs_output_5 = pd.read_csv('CorrWithin0.75/Output_S0.5_F2.0_Corr0.75_Migroutes2_Nit1000_Carryover1.csv', header = None)
#abs_output_5 = abs_output_5.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_5 = {}
#columns_5 = abs_output_5.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_5, 2):
#     correlations_5[col_a, '__' , col_b] = stats.pearsonr(abs_output_5.loc[:, col_a], abs_output_5.loc[:, col_b])
#result_5 = pd.DataFrame.from_dict(correlations_5, orient='index')
#result_5.columns = ['PCC', 'p-value']
#result_5['distances'] = np.array(flatdist_1)
#average_corr_r2075 = result_5.groupby('distances')['PCC'].mean()
#average_corr_r2075.drop(index=average_corr_r2075.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.75 correlation, K
#abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.75_Migroutes2_Nit1000_Carryover1.csv', header = None)
#abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_3 = {}
#columns_3 = abs_output_3.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_3, 2):
#     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
#result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
#result_3.columns = ['PCC', 'p-value']
#result_3['distances'] = np.array(flatdist_1)
#average_corr_k2075 = result_3.groupby('distances')['PCC'].mean()
#average_corr_k2075.drop(index=average_corr_k2075.index[0], axis=0, inplace=True)
#average_corr_smoothed_k3 = gaussian_filter1d(average_corr_k3, sigma=2)


# 2 Mig Route: 1 correlation, r
#abs_output_11 = pd.read_csv('CorrWithin0.75/Output_S0.5_F2.0_Corr1_Migroutes2_Nit1000_Carryover1.csv', header = None)
#abs_output_11 = abs_output_11.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_11 = {}
#columns_11 = abs_output_11.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_11, 2):
#     correlations_11[col_a, '__' , col_b] = stats.pearsonr(abs_output_11.loc[:, col_a], abs_output_11.loc[:, col_b])
#result_11 = pd.DataFrame.from_dict(correlations_11, orient='index')
#result_11.columns = ['PCC', 'p-value']
#result_11['distances'] = np.array(flatdist_1)
#average_corr_r21 = result_11.groupby('distances')['PCC'].mean()
#average_corr_r21.drop(index=average_corr_r21.index[0], axis=0, inplace=True)


# 2 Mig Route: 1 correlation, K
#abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr1_Migroutes2_Nit1000_Carryover1.csv', header = None)
#abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_3 = {}
#columns_3 = abs_output_3.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_3, 2):
#     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
#result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
#result_3.columns = ['PCC', 'p-value']
#result_3['distances'] = np.array(flatdist_1)
#average_corr_k21 = result_3.groupby('distances')['PCC'].mean()
#average_corr_k21.drop(index=average_corr_k21.index[0], axis=0, inplace=True)
#average_corr_smoothed_k3 = gaussian_filter1d(average_corr_k3, sigma=2)

########### 4 MIGRATION ROUTES ############################################################################################################

# 4 Mig Route: 0 correlation, r
abs_output_6 = pd.read_csv('CorrWithin0.75/Output_S0.5_F2.0_Corr0_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_6 = abs_output_6.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_6 = {}
columns_6 = abs_output_6.columns.tolist()
for col_a, col_b in itertools.combinations(columns_6, 2):
     correlations_6[col_a, '__' , col_b] = stats.pearsonr(abs_output_6.loc[:, col_a], abs_output_6.loc[:, col_b])
result_6 = pd.DataFrame.from_dict(correlations_6, orient='index')
result_6.columns = ['PCC', 'p-value']
result_6['distances'] = np.array(flatdist_1)
average_corr_r40 = result_6.groupby('distances')['PCC'].mean()
average_corr_r40.drop(index=average_corr_r40.index[0], axis=0, inplace=True)


# 4 Mig Route: 0 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k40 = result_3.groupby('distances')['PCC'].mean()
average_corr_k40.drop(index=average_corr_k40.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.25 correlation, r
abs_output_8 = pd.read_csv('CorrWithin0.75/Output_S0.5_F2.0_Corr0.25_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_8 = {}
columns_8 = abs_output_8.columns.tolist()
for col_a, col_b in itertools.combinations(columns_8, 2):
     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
result_8.columns = ['PCC', 'p-value']
result_8['distances'] = np.array(flatdist_1)
average_corr_r4025 = result_8.groupby('distances')['PCC'].mean()
average_corr_r4025.drop(index=average_corr_r4025.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.25 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.25_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k4025 = result_3.groupby('distances')['PCC'].mean()
average_corr_k4025.drop(index=average_corr_k4025.index[0], axis=0, inplace=True)

# 4 Mig Route: 0.5 correlation, r
abs_output_8 = pd.read_csv('CorrWithin0.75/Output_S0.5_F2.0_Corr0.5_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_8 = {}
columns_8 = abs_output_8.columns.tolist()
for col_a, col_b in itertools.combinations(columns_8, 2):
     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
result_8.columns = ['PCC', 'p-value']
result_8['distances'] = np.array(flatdist_1)
average_corr_r4050 = result_8.groupby('distances')['PCC'].mean()
average_corr_r4050.drop(index=average_corr_r4025.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.5 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.5_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k4050 = result_3.groupby('distances')['PCC'].mean()
average_corr_k4050.drop(index=average_corr_k4025.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.75 correlation, r
abs_output_9 = pd.read_csv('CorrWithin0.75/Output_S0.5_F2.0_Corr0.75_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_9 = abs_output_9.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_9 = {}
columns_9 = abs_output_9.columns.tolist()
for col_a, col_b in itertools.combinations(columns_9, 2):
     correlations_9[col_a, '__' , col_b] = stats.pearsonr(abs_output_9.loc[:, col_a], abs_output_9.loc[:, col_b])
result_9 = pd.DataFrame.from_dict(correlations_9, orient='index')
result_9.columns = ['PCC', 'p-value']
result_9['distances'] = np.array(flatdist_1)
average_corr_r4075 = result_9.groupby('distances')['PCC'].mean()
average_corr_r4075.drop(index=average_corr_r4075.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.75 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.75_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k4075 = result_3.groupby('distances')['PCC'].mean()
average_corr_k4075.drop(index=average_corr_k4075.index[0], axis=0, inplace=True)


# 4 Mig Route: 1 correlation, r
abs_output_7 = pd.read_csv('CorrWithin0.75/Output_S0.5_F2.0_Corr1_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_7 = abs_output_7.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_7 = {}
columns_7 = abs_output_7.columns.tolist()
for col_a, col_b in itertools.combinations(columns_7, 2):
     correlations_7[col_a, '__' , col_b] = stats.pearsonr(abs_output_7.loc[:, col_a], abs_output_7.loc[:, col_b])
result_7 = pd.DataFrame.from_dict(correlations_7, orient='index')
result_7.columns = ['PCC', 'p-value']
result_7['distances'] = np.array(flatdist_1)
average_corr_r41 = result_7.groupby('distances')['PCC'].mean()
average_corr_r41.drop(index=average_corr_r41.index[0], axis=0, inplace=True)

# 4 Mig Route: 1 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr1_Migroutes4_Nit1000_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k41 = result_3.groupby('distances')['PCC'].mean()
average_corr_k41.drop(index=average_corr_k41.index[0], axis=0, inplace=True)

## No Migration K selected Species
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0_Migroutes0_Nit1000_Carryover0.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k0 = result_3.groupby('distances')['PCC'].mean()
average_corr_k0.drop(index=average_corr_k0.index[0], axis=0, inplace=True)



## Plotting 4 migration routes comparing autocorrelation, K selected species
plt.plot(average_corr_k0, linestyle="solid", color="grey", label = "No migration")
#plt.plot(average_corr_1migk, color="black", linestyle="solid", linewidth=2, label = "1 nonbreeding ground")
plt.plot(average_corr_k40, linestyle="dashed", color="darkgreen", linewidth=2, label = "0 corr") ## 
plt.plot(average_corr_k4025, linestyle="dotted",  color ="orange",  linewidth=2, label = "0.25 corr") ## 
plt.plot(average_corr_k4050, linestyle="dashdot",  color ="blue",  linewidth=2, label = "0.50 corr") ## 
plt.plot(average_corr_k4075, linestyle="dashdot",  color = "slategrey", linewidth=2, label = "0.75 corr")  ##
plt.plot(average_corr_k41, linestyle="--", color="purple", label = "1.0 corr") ## 


## Plotting r vs K simple plot


plt.plot(average_corr_r4025, linestyle="solid",  color = "black", linewidth=2, label = "r-selected species, 0.25 corr")  ##
plt.plot(average_corr_k4025, linestyle="dotted", color="orange", label = "K-selected species, 0.25 corr") ## 

plt.plot(average_corr_r4075, linestyle="solid",  color = "grey", linewidth=2, label = "r-selected species, 0.75 corr")  ##
plt.plot(average_corr_k4075, linestyle="dotted", color="red", label = "K-selected species, 0.75 corr") ## 




plt.title('')
plt.xlabel('Distance')
plt.xlim([3, 50])
plt.ylim([0.75, 1])
plt.ylim(ymin=0.)
plt.ylabel('Average Correlation')
plt.legend()
plt.legend(loc=(1.04, 0), frameon=False, title="")
plt.show()

plt.rcParams.update({'font.size': 14})


#############################################################################################################
########## FIGURE 4 #########################################################################################
#############################################################################################################


######## RANDOM MIGRATION #######################################
#######################################################################################
## No Migration K selected Species
## Import distance file
distances_1 = pd.read_csv('Midpoint_Distances_Gridsize150_Thin2.csv', header = None)
distances_1 = distances_1.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
distlist_1 = distances_1.values.tolist()
flatdist_1 = (np.array(distlist_1)).flatten()


abs_output_3 = pd.read_csv('CorrWithin0.75_Fig4/Output_S0.9_F0.25_Corr0_Migroutes0_Nit500_Carryover0.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k0 = result_3.groupby('distances')['PCC'].mean()
average_corr_k0.drop(index=average_corr_k0.index[0], axis=0, inplace=True)



# 2 Mig Route: 0 correlation, K
#abs_output_3 = pd.read_csv('CorrWithin0.75_RandomMig/##Output_S0.9_F0.23_Corr0_Migroutes2_Nit1000_DemStoch0_Carryover1_Random.csv', header = None)
#abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_3 = {}
#columns_3 = abs_output_3.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_3, 2):
#     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
#result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
#result_3.columns = ['PCC', 'p-value']
#result_3['distances'] = np.array(flatdist_1)
#average_corr_k20_RANDOM = result_3.groupby('distances')['PCC'].mean()
#average_corr_k20_RANDOM.drop(index=average_corr_k20_RANDOM.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.25 correlation, K
#abs_output_4 = pd.read_csv('Random_Grid_Carryover1/Output_S0.9_F0.23_Corr0.25_Migroutes2_Nit500_Carryover1.csv', header = None)
#abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_4 = {}
#columns_4 = abs_output_4.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_4, 2):
#     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
#result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
#result_4.columns = ['PCC', 'p-value']
#result_4['distances'] = np.array(flatdist_1)
#average_corr_k2025_RANDOM = result_4.groupby('distances')['PCC'].mean()
#average_corr_k2025_RANDOM.drop(index=average_corr_k2025_RANDOM.index[0], axis=0, inplace=True)

# 2 Mig Route: 0.5 correlation, K
#abs_output_4 = pd.read_csv('Random_Grid_Carryover1/#Output_S0.9_F0.23_Corr0.5_Migroutes2_Nit1000_DemStoch0_Carryover1_Random.csv', header = None)
#abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_4 = {}
#columns_4 = abs_output_4.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_4, 2):
#     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
#result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
#result_4.columns = ['PCC', 'p-value']
#result_4['distances'] = np.array(flatdist_1)
#average_corr_k2050_RANDOM = result_4.groupby('distances')['PCC'].mean()
#average_corr_k2050_RANDOM.drop(index=average_corr_k2050_RANDOM.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.75 correlation,  K
#abs_output_5 = pd.read_csv('Random_Grid_Carryover1/Output_S0.9_F0.25_Corr0.75_Migroutes2_Nit500_Carryover1.csv', header = None)
#abs_output_5 = abs_output_5.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_5 = {}
#columns_5 = abs_output_5.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_5, 2):
#     correlations_5[col_a, '__' , col_b] = stats.pearsonr(abs_output_5.loc[:, col_a], abs_output_5.loc[:, col_b])
#result_5 = pd.DataFrame.from_dict(correlations_5, orient='index')
#result_5.columns = ['PCC', 'p-value']
#result_5['distances'] = np.array(flatdist_1)
#average_corr_k2075_UNEVEN = result_5.groupby('distances')['PCC'].mean()
#average_corr_k2075_UNEVEN.drop(index=average_corr_k2075_RANDOM.index[0], axis=0, inplace=True)

# 2 Mig Route: 1 correlation, K
#abs_output_11 = pd.read_csv('Random_Grid_Carryover1/#Output_S0.9_F0.23_Corr1_Migroutes2_Nit1000_DemStoch0_Carryover1_Random.csv', header = None)
#abs_output_11 = abs_output_11.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_11 = {}
#columns_11 = abs_output_11.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_11, 2):
#     correlations_11[col_a, '__' , col_b] = stats.pearsonr(abs_output_11.loc[:, col_a], abs_output_11.loc[:, #col_b])
#result_11 = pd.DataFrame.from_dict(correlations_11, orient='index')
#result_11.columns = ['PCC', 'p-value']
#result_11['distances'] = np.array(flatdist_1)
#average_corr_k21_RANDOM = result_11.groupby('distances')['PCC'].mean()
#average_corr_k21_RANDOM.drop(index=average_corr_k21_RANDOM.index[0], axis=0, inplace=True)

# 4 Mig Route: 0 correlation, K
abs_output_6 = pd.read_csv('CorrWithin0.75_RandomMig/Output_S0.9_F0.25_Corr0_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_6 = abs_output_6.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_6 = {}
columns_6 = abs_output_6.columns.tolist()
for col_a, col_b in itertools.combinations(columns_6, 2):
     correlations_6[col_a, '__' , col_b] = stats.pearsonr(abs_output_6.loc[:, col_a], abs_output_6.loc[:, col_b])
result_6 = pd.DataFrame.from_dict(correlations_6, orient='index')
result_6.columns = ['PCC', 'p-value']
result_6['distances'] = np.array(flatdist_1)
average_corr_k40_RANDOM = result_6.groupby('distances')['PCC'].mean()
average_corr_k40_RANDOM.drop(index=average_corr_k40_RANDOM.index[0], axis=0, inplace=True)

# 4 Mig Route: 0.25 correlation, K
#abs_output_8 = pd.read_csv('Random_Grid_Carryover1/Output_S0.9_F0.25_Corr0.25_Migroutes4_Nit500_Carryover1.csv', header = None)
#abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_8 = {}
#columns_8 = abs_output_8.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_8, 2):
#     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
#result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
#result_8.columns = ['PCC', 'p-value']
#result_8['distances'] = np.array(flatdist_1)
#average_corr_k4025_UNEVEN = result_8.groupby('distances')['PCC'].mean()
#average_corr_k4025_UNEVEN.drop(index=average_corr_k4025_UNEVEN.index[0], axis=0, inplace=True)

# 4 Mig Route: 0.5 correlation, K
abs_output_8 = pd.read_csv('CorrWithin0.75_RandomMig/Output_S0.9_F0.25_Corr0.5_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_8 = {}
columns_8 = abs_output_8.columns.tolist()
for col_a, col_b in itertools.combinations(columns_8, 2):
     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
result_8.columns = ['PCC', 'p-value']
result_8['distances'] = np.array(flatdist_1)
average_corr_k4050_RANDOM = result_8.groupby('distances')['PCC'].mean()
average_corr_k4050_RANDOM.drop(index=average_corr_k4050_RANDOM.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.75 correlation, K
#abs_output_9 = pd.read_csv('Random_Grid_Carryover1/Output_S0.9_F0.25_Corr0.75_Migroutes4_Nit500_Carryover1.csv', header = None)
#abs_output_9 = abs_output_9.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_9 = {}
#columns_9 = abs_output_9.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_9, 2):
#     correlations_9[col_a, '__' , col_b] = stats.pearsonr(abs_output_9.loc[:, col_a], abs_output_9.loc[:, col_b])
#result_9 = pd.DataFrame.from_dict(correlations_9, orient='index')
#result_9.columns = ['PCC', 'p-value']
#result_9['distances'] = np.array(flatdist_1)
#average_corr_k4075_UNEVEN = result_9.groupby('distances')['PCC'].mean()
#average_corr_k4075_UNEVEN.drop(index=average_corr_k4075_UNEVEN.index[0], axis=0, inplace=True)

# 4 Mig Route: 1 correlation, K
abs_output_7 = pd.read_csv('CorrWithin0.75_RandomMig/Output_S0.9_F0.25_Corr1_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_7 = abs_output_7.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_7 = {}
columns_7 = abs_output_7.columns.tolist()
for col_a, col_b in itertools.combinations(columns_7, 2):
     correlations_7[col_a, '__' , col_b] = stats.pearsonr(abs_output_7.loc[:, col_a], abs_output_7.loc[:, col_b])
result_7 = pd.DataFrame.from_dict(correlations_7, orient='index')
result_7.columns = ['PCC', 'p-value']
result_7['distances'] = np.array(flatdist_1)
average_corr_k41_RANDOM = result_7.groupby('distances')['PCC'].mean()
average_corr_k41_RANDOM.drop(index=average_corr_k41_RANDOM.index[0], axis=0, inplace=True)


######## PROXIMITY GRID MIGRATION ############# #######################################
#######################################################################################

# 2 Mig Route: 0 correlation, K
abs_output_6 = pd.read_csv('CorrWithin0.75_Fig4/Output_S0.9_F0.25_Corr0_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_6 = abs_output_6.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_6 = {}
columns_6 = abs_output_6.columns.tolist()
for col_a, col_b in itertools.combinations(columns_6, 2):
     correlations_6[col_a, '__' , col_b] = stats.pearsonr(abs_output_6.loc[:, col_a], abs_output_6.loc[:, col_b])
result_6 = pd.DataFrame.from_dict(correlations_6, orient='index')
result_6.columns = ['PCC', 'p-value']
result_6['distances'] = np.array(flatdist_1)
average_corr_k20 = result_6.groupby('distances')['PCC'].mean()
average_corr_k20.drop(index=average_corr_k20.index[0], axis=0, inplace=True)



# 2 Mig Route: 0.50 correlation, K
abs_output_9 = pd.read_csv('CorrWithin0.75_Fig4/Output_S0.9_F0.25_Corr0.5_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_9 = abs_output_9.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_9 = {}
columns_9 = abs_output_9.columns.tolist()
for col_a, col_b in itertools.combinations(columns_9, 2):
     correlations_9[col_a, '__' , col_b] = stats.pearsonr(abs_output_9.loc[:, col_a], abs_output_9.loc[:, col_b])
result_9 = pd.DataFrame.from_dict(correlations_9, orient='index')
result_9.columns = ['PCC', 'p-value']
result_9['distances'] = np.array(flatdist_1)
average_corr_k2050 = result_9.groupby('distances')['PCC'].mean()
average_corr_k2050.drop(index=average_corr_k2050.index[0], axis=0, inplace=True)



# 2 Mig Route: 1 correlation, K
abs_output_7 = pd.read_csv('CorrWithin0.75_Fig4/Output_S0.9_F0.25_Corr1_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_7 = abs_output_7.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_7 = {}
columns_7 = abs_output_7.columns.tolist()
for col_a, col_b in itertools.combinations(columns_7, 2):
     correlations_7[col_a, '__' , col_b] = stats.pearsonr(abs_output_7.loc[:, col_a], abs_output_7.loc[:, col_b])
result_7 = pd.DataFrame.from_dict(correlations_7, orient='index')
result_7.columns = ['PCC', 'p-value']
result_7['distances'] = np.array(flatdist_1)
average_corr_k21 = result_7.groupby('distances')['PCC'].mean()
average_corr_k21.drop(index=average_corr_k21.index[0], axis=0, inplace=True)



# 4 Mig Route: 0 correlation, K
abs_output_6 = pd.read_csv('CorrWithin0.75_Fig4/Output_S0.9_F0.25_Corr0_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_6 = abs_output_6.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_6 = {}
columns_6 = abs_output_6.columns.tolist()
for col_a, col_b in itertools.combinations(columns_6, 2):
     correlations_6[col_a, '__' , col_b] = stats.pearsonr(abs_output_6.loc[:, col_a], abs_output_6.loc[:, col_b])
result_6 = pd.DataFrame.from_dict(correlations_6, orient='index')
result_6.columns = ['PCC', 'p-value']
result_6['distances'] = np.array(flatdist_1)
average_corr_k40 = result_6.groupby('distances')['PCC'].mean()
average_corr_k40.drop(index=average_corr_k40.index[0], axis=0, inplace=True)



# 4 Mig Route: 0.50 correlation, K
abs_output_9 = pd.read_csv('CorrWithin0.75_Fig4/Output_S0.9_F0.25_Corr0.5_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_9 = abs_output_9.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_9 = {}
columns_9 = abs_output_9.columns.tolist()
for col_a, col_b in itertools.combinations(columns_9, 2):
     correlations_9[col_a, '__' , col_b] = stats.pearsonr(abs_output_9.loc[:, col_a], abs_output_9.loc[:, col_b])
result_9 = pd.DataFrame.from_dict(correlations_9, orient='index')
result_9.columns = ['PCC', 'p-value']
result_9['distances'] = np.array(flatdist_1)
average_corr_k4050 = result_9.groupby('distances')['PCC'].mean()
average_corr_k4050.drop(index=average_corr_k4050.index[0], axis=0, inplace=True)


# 4 Mig Route: 1 correlation, K
abs_output_7 = pd.read_csv('CorrWithin0.75_Fig4/Output_S0.9_F0.25_Corr1_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_7 = abs_output_7.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_7 = {}
columns_7 = abs_output_7.columns.tolist()
for col_a, col_b in itertools.combinations(columns_7, 2):
     correlations_7[col_a, '__' , col_b] = stats.pearsonr(abs_output_7.loc[:, col_a], abs_output_7.loc[:, col_b])
result_7 = pd.DataFrame.from_dict(correlations_7, orient='index')
result_7.columns = ['PCC', 'p-value']
result_7['distances'] = np.array(flatdist_1)
average_corr_k41 = result_7.groupby('distances')['PCC'].mean()
average_corr_k41.drop(index=average_corr_k41.index[0], axis=0, inplace=True)



######################################################################################################################
## Plotting FIGURE 4A: Proximity vs Random Migration
    
    
plt.plot(average_corr_k0, linestyle="dotted", color="grey", linewidth=2, label="No migration") ## No migraiton, K selected species.

plt.plot(average_corr_k40, linestyle="solid", color="darkgreen", linewidth=2, label="0 corr, Proximity migration") ## Proximity migration, 4 mig routes, 0 corr
plt.plot(average_corr_k40_RANDOM, linestyle="dashed", color="darkgreen", linewidth=2, label="0 corr, Random migration") ## Random migration, 4 mig routes, 0 corr

plt.plot(average_corr_k4050, linestyle="solid",  color ="orange",  linewidth=2, label="0.5 corr, Proximity migration") ## Proximity migration, 4 mig routes, 0.5 corr
plt.plot(average_corr_k4050_RANDOM, linestyle="dashed",  color="orange", linewidth=2, label="0.5 corr, Random migration") ## Random migration, 4 mig routes, 0.5 corr

plt.plot(average_corr_k41, linestyle="solid", color="purple", linewidth=2, label="1 corr, Proximity migration") ## Proximity migration, 4 mig routes, 1 corr
plt.plot(average_corr_k41_RANDOM,linestyle="dashed", color="purple", linewidth=2, label="1 corr, Random migration") ## Random migration, 4 mig routes, 1 corr



plt.title('')
plt.xlabel('Distance')
plt.xlim([3, 50])
plt.ylim([0, 1])
plt.ylim(ymin=0.)
plt.ylabel('Average Correlation')
plt.legend()
plt.legend(loc=(1.04, 0), frameon=False, title="")
plt.show()

plt.rcParams.update({'font.size': 14})




## Plotting FIGURE 4B: 2 mig routes vs 4 mig routes (comparing number and size of mig routes)

plt.plot(average_corr_k0, linestyle="dotted", color="grey", linewidth=2, label="No migration") ## No migration, K selected species.

plt.plot(average_corr_k20, linestyle="dotted", color="darkgreen", linewidth=2, label="2 nonbreeding grounds") ## K selected species, 2 mig routes, 0 corr. (proximity mig)
plt.plot(average_corr_k40, color="darkgreen", linewidth=2)  ## K selected species, 4 mig routes, 0 corr. (proximity mig)

plt.plot(average_corr_k2050, linestyle="dotted",  color ="orange",  linewidth=2, label="2 nonbreeding grounds")  ## K selected species, 2 mig routes, 0.5 corr. (proximity mig)
plt.plot(average_corr_k4050, color="orange", linewidth=2)  ## K selected species, 4 mig routes, 0.5 corr. (proximity mig)


plt.plot(average_corr_k21, linestyle="dotted", color="purple", linewidth=2, label="2 nonbreeding grounds") ## K selected species, 2 mig routes, 1 corr. (proximity mig)
plt.plot(average_corr_k41, color="purple", linewidth=2)   ## K selected species, 4 mig routes, 1 corr. (proximity mig)


plt.title('')
plt.xlabel('Distance')
plt.xlim([3, 50])
plt.ylim([0, 1])
plt.ylim(ymin=0.)
plt.ylabel('Average Correlation')
plt.legend()
plt.legend(loc=(1.04, 0), frameon=False, title="")
plt.show()

plt.rcParams.update({'font.size': 14})

#############################################################################################################
########## FIGURE 5 #########################################################################################
#############################################################################################################

# -*- coding: utf-8 -*-
"""
######## FIGURE FOR UNEVEN BREEDING GROUND POPULATIONS #########


Created on Mon Sep  5 11:05:10 2022

@author: ellencm

"""
import itertools
from itertools import combinations
from scipy import stats
from scipy.ndimage.filters import gaussian_filter1d
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


## Set working directory:
import os
if os.environ['HOMEPATH'] == "\\Users\\ellencm":
  os.chdir('//home.ansatt.ntnu.no/ellencm/Desktop/SIMOUTPUT_AUGUST/')


## No Migration K selected Species
## Import distance file
distances_1 = pd.read_csv('Midpoint_Distances_Gridsize150_Thin2.csv', header = None)
distances_1 = distances_1.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
distlist_1 = distances_1.values.tolist()
flatdist_1 = (np.array(distlist_1)).flatten()


######## UNEVEN SIZED BREEDING GROUND MIGRATION #######################################
#######################################################################################
# 2 Mig Route: 0 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75_TEST/Output_S0.9_F0.25_Corr0_Migroutes2_Nit500_Carryover1_UNEVEN.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k20_UNEVEN = result_3.groupby('distances')['PCC'].mean()
average_corr_k20_UNEVEN.drop(index=average_corr_k20_UNEVEN.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.25 correlation, K
abs_output_4 = pd.read_csv('CorrWithin0.75_UnevenPopSizes/Output_S0.9_F0.25_Corr0.25_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_4 = {}
columns_4 = abs_output_4.columns.tolist()
for col_a, col_b in itertools.combinations(columns_4, 2):
     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
result_4.columns = ['PCC', 'p-value']
result_4['distances'] = np.array(flatdist_1)
average_corr_k2025_UNEVEN = result_4.groupby('distances')['PCC'].mean()
average_corr_k2025_UNEVEN.drop(index=average_corr_k2025_UNEVEN.index[0], axis=0, inplace=True)

# 2 Mig Route: 0.5 correlation, K
abs_output_4 = pd.read_csv('CorrWithin0.75_TEST/Output_S0.9_F0.25_Corr0.5_Migroutes2_Nit500_Carryover1_UNEVEN.csv', header = None)
abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_4 = {}
columns_4 = abs_output_4.columns.tolist()
for col_a, col_b in itertools.combinations(columns_4, 2):
     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
result_4.columns = ['PCC', 'p-value']
result_4['distances'] = np.array(flatdist_1)
average_corr_k2050_UNEVEN = result_4.groupby('distances')['PCC'].mean()
average_corr_k2050_UNEVEN.drop(index=average_corr_k2050_UNEVEN.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.75 correlation,  K
#abs_output_5 = pd.read_csv('CorrWithin0.75_TEST/#Output_S0.9_F0.25_Corr0.75_Migroutes2_Nit500_Carryover1_UNEVEN.csv', header = None)
#abs_output_5 = abs_output_5.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_5 = {}
#columns_5 = abs_output_5.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_5, 2):
#     correlations_5[col_a, '__' , col_b] = stats.pearsonr(abs_output_5.loc[:, col_a], abs_output_5.loc[:, col_b])
#result_5 = pd.DataFrame.from_dict(correlations_5, orient='index')
#result_5.columns = ['PCC', 'p-value']
#result_5['distances'] = np.array(flatdist_1)
#average_corr_k2075_UNEVEN = result_5.groupby('distances')['PCC'].mean()
#average_corr_k2075_UNEVEN.drop(index=average_corr_k2075_UNEVEN.index[0], axis=0, inplace=True)

# 2 Mig Route: 1 correlation, K
abs_output_11 = pd.read_csv('CorrWithin0.75_TEST/Output_S0.9_F0.25_Corr1_Migroutes2_Nit500_Carryover1_UNEVEN.csv', header = None)
abs_output_11 = abs_output_11.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_11 = {}
columns_11 = abs_output_11.columns.tolist()
for col_a, col_b in itertools.combinations(columns_11, 2):
     correlations_11[col_a, '__' , col_b] = stats.pearsonr(abs_output_11.loc[:, col_a], abs_output_11.loc[:, col_b])
result_11 = pd.DataFrame.from_dict(correlations_11, orient='index')
result_11.columns = ['PCC', 'p-value']
result_11['distances'] = np.array(flatdist_1)
average_corr_k21_UNEVEN = result_11.groupby('distances')['PCC'].mean()
average_corr_k21_UNEVEN.drop(index=average_corr_k21_UNEVEN.index[0], axis=0, inplace=True)

# 4 Mig Route: 0 correlation, K
abs_output_6 = pd.read_csv('CorrWithin0.75_TEST/Output_S0.9_F0.25_Corr0_Migroutes4_Nit500_Carryover1_UNEVEN.csv', header = None)
abs_output_6 = abs_output_6.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_6 = {}
columns_6 = abs_output_6.columns.tolist()
for col_a, col_b in itertools.combinations(columns_6, 2):
     correlations_6[col_a, '__' , col_b] = stats.pearsonr(abs_output_6.loc[:, col_a], abs_output_6.loc[:, col_b])
result_6 = pd.DataFrame.from_dict(correlations_6, orient='index')
result_6.columns = ['PCC', 'p-value']
result_6['distances'] = np.array(flatdist_1)
average_corr_k40_UNEVEN = result_6.groupby('distances')['PCC'].mean()
average_corr_k40_UNEVEN.drop(index=average_corr_k40_UNEVEN.index[0], axis=0, inplace=True)

# 4 Mig Route: 0.25 correlation, K
abs_output_8 = pd.read_csv('CorrWithin0.75_UnevenPopSizes/Output_S0.9_F0.25_Corr0.25_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_8 = {}
columns_8 = abs_output_8.columns.tolist()
for col_a, col_b in itertools.combinations(columns_8, 2):
     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
result_8.columns = ['PCC', 'p-value']
result_8['distances'] = np.array(flatdist_1)
average_corr_k4025_UNEVEN = result_8.groupby('distances')['PCC'].mean()
average_corr_k4025_UNEVEN.drop(index=average_corr_k4025_UNEVEN.index[0], axis=0, inplace=True)



# 4 Mig Route: 0.5 correlation, K
abs_output_8 = pd.read_csv('CorrWithin0.75_TEST/Output_S0.9_F0.25_Corr0.5_Migroutes4_Nit500_Carryover1_UNEVEN.csv', header = None)
abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_8 = {}
columns_8 = abs_output_8.columns.tolist()
for col_a, col_b in itertools.combinations(columns_8, 2):
     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
result_8.columns = ['PCC', 'p-value']
result_8['distances'] = np.array(flatdist_1)
average_corr_k4050_UNEVEN = result_8.groupby('distances')['PCC'].mean()
average_corr_k4050_UNEVEN.drop(index=average_corr_k4050_UNEVEN.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.75 correlation, K
#abs_output_9 = pd.read_csv('CorrWithin0.75_TEST/#Output_S0.9_F0.25_Corr0.75_Migroutes4_Nit500_Carryover1_UNEVEN.csv', header = None)
# = abs_output_9.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_9 = {}
#columns_9 = abs_output_9.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_9, 2):
#     correlations_9[col_a, '__' , col_b] = stats.pearsonr(abs_output_9.loc[:, col_a], abs_output_9.loc[:, col_b])
#result_9 = pd.DataFrame.from_dict(correlations_9, orient='index')
#result_9.columns = ['PCC', 'p-value']
#result_9['distances'] = np.array(flatdist_1)
#average_corr_k4075_UNEVEN = result_9.groupby('distances')['PCC'].mean()
#average_corr_k4075_UNEVEN.drop(index=average_corr_k4075_UNEVEN.index[0], axis=0, inplace=True)

# 4 Mig Route: 1 correlation, K
abs_output_7 = pd.read_csv('CorrWithin0.75_TEST/Output_S0.9_F0.25_Corr1_Migroutes4_Nit500_Carryover1_UNEVEN.csv', header = None)
abs_output_7 = abs_output_7.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_7 = {}
columns_7 = abs_output_7.columns.tolist()
for col_a, col_b in itertools.combinations(columns_7, 2):
     correlations_7[col_a, '__' , col_b] = stats.pearsonr(abs_output_7.loc[:, col_a], abs_output_7.loc[:, col_b])
result_7 = pd.DataFrame.from_dict(correlations_7, orient='index')
result_7.columns = ['PCC', 'p-value']
result_7['distances'] = np.array(flatdist_1)
average_corr_k41_UNEVEN = result_7.groupby('distances')['PCC'].mean()
average_corr_k41_UNEVEN.drop(index=average_corr_k41_UNEVEN.index[0], axis=0, inplace=True)







### EVEN SIZED MIGRATION ############################################### ####################################################################################################

### NO MIGRATION, K SELECTED 
abs_output_1k = pd.read_csv('CorrWithin0.75_Fig5/Output_S0.9_F0.25_Corr0_Migroutes0_Nit500_Carryover0.csv', header = None)
distances_1 = pd.read_csv('Midpoint_Distances_Gridsize150_Thin2.csv', header = None)
abs_output_1k = abs_output_1k.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
distances_1 = distances_1.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
distlist_1 = distances_1.values.tolist()
flatdist_1 = (np.array(distlist_1)).flatten()
correlations_1k = {}
columns_1k = abs_output_1k.columns.tolist()
for col_a, col_b in itertools.combinations(columns_1k, 2):
     correlations_1k[col_a, '__' , col_b] = stats.pearsonr(abs_output_1k.loc[:, col_a], abs_output_1k.loc[:, col_b])
result_1k = pd.DataFrame.from_dict(correlations_1k, orient='index')
result_1k.columns = ['PCC', 'p-value']
result_1k['distances'] = np.array(flatdist_1)
average_corr_k0 = result_1k.groupby('distances')['PCC'].mean()
average_corr_k0.drop(index=average_corr_k0.index[0], axis=0, inplace=True)



# 2 Mig Route: 0 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k20 = result_3.groupby('distances')['PCC'].mean()
average_corr_k20.drop(index=average_corr_k20.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.25 correlation, K
#abs_output_4 = pd.read_csv('CorrWithin0.75_Fig5/Output_S0.9_F0.25_Corr0.25_Migroutes2_Nit500_Carryover1.csv', #header = None)
#abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_4 = {}
#columns_4 = abs_output_4.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_4, 2):
#     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
#result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
#result_4.columns = ['PCC', 'p-value']
#result_4['distances'] = np.array(flatdist_1)
#average_corr_k2025 = result_4.groupby('distances')['PCC'].mean()
#average_corr_k2025.drop(index=average_corr_k2025.index[0], axis=0, inplace=True)



# 2 Mig Route: 0.5 correlation, K
abs_output_4 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.5_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_4 = {}
columns_4 = abs_output_4.columns.tolist()
for col_a, col_b in itertools.combinations(columns_4, 2):
     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
result_4.columns = ['PCC', 'p-value']
result_4['distances'] = np.array(flatdist_1)
average_corr_k2050 = result_4.groupby('distances')['PCC'].mean()
average_corr_k2050.drop(index=average_corr_k2050.index[0], axis=0, inplace=True)



# 2 Mig Route: 0.75 correlation, K
#abs_output_5 = pd.read_csv('CorrWithin0.75_Fig5/Output_S0.9_F0.25_Corr0.75_Migroutes2_Nit500_Carryover1.csv', #header = None)
#abs_output_5 = abs_output_5.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
#correlations_5 = {}
#columns_5 = abs_output_5.columns.tolist()
#for col_a, col_b in itertools.combinations(columns_5, 2):
#     correlations_5[col_a, '__' , col_b] = stats.pearsonr(abs_output_5.loc[:, col_a], abs_output_5.loc[:, col_b])
#result_5 = pd.DataFrame.from_dict(correlations_5, orient='index')
#result_5.columns = ['PCC', 'p-value']
#result_5['distances'] = np.array(flatdist_1)
#average_corr_k2075 = result_5.groupby('distances')['PCC'].mean()
#average_corr_k2075.drop(index=average_corr_k2075.index[0], axis=0, inplace=True)


# 2 Mig Route: 1 correlation, K
abs_output_11 = pd.read_csv('CorrWithin0.75_TEST/Output_S0.9_F0.25_Corr1_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_11 = abs_output_11.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_11 = {}
columns_11 = abs_output_11.columns.tolist()
for col_a, col_b in itertools.combinations(columns_11, 2):
     correlations_11[col_a, '__' , col_b] = stats.pearsonr(abs_output_11.loc[:, col_a], abs_output_11.loc[:, col_b])
result_11 = pd.DataFrame.from_dict(correlations_11, orient='index')
result_11.columns = ['PCC', 'p-value']
result_11['distances'] = np.array(flatdist_1)
average_corr_k21 = result_11.groupby('distances')['PCC'].mean()
average_corr_k21.drop(index=average_corr_k21.index[0], axis=0, inplace=True)


# 4 Mig Route: 0 correlation, K
abs_output_6 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_6 = abs_output_6.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_6 = {}
columns_6 = abs_output_6.columns.tolist()
for col_a, col_b in itertools.combinations(columns_6, 2):
     correlations_6[col_a, '__' , col_b] = stats.pearsonr(abs_output_6.loc[:, col_a], abs_output_6.loc[:, col_b])
result_6 = pd.DataFrame.from_dict(correlations_6, orient='index')
result_6.columns = ['PCC', 'p-value']
result_6['distances'] = np.array(flatdist_1)
average_corr_k40 = result_6.groupby('distances')['PCC'].mean()
average_corr_k40.drop(index=average_corr_k40.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.25 correlation, K
abs_output_8 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.25_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_8 = {}
columns_8 = abs_output_8.columns.tolist()
for col_a, col_b in itertools.combinations(columns_8, 2):
     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
result_8.columns = ['PCC', 'p-value']
result_8['distances'] = np.array(flatdist_1)
average_corr_k4025 = result_8.groupby('distances')['PCC'].mean()
average_corr_k4025.drop(index=average_corr_k4025.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.5 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.5_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k4050 = result_3.groupby('distances')['PCC'].mean()
average_corr_k4050.drop(index=average_corr_k4025.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.75 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.75_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k4075 = result_3.groupby('distances')['PCC'].mean()
average_corr_k4075.drop(index=average_corr_k4075.index[0], axis=0, inplace=True)


# 4 Mig Route: 1 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr1_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k41 = result_3.groupby('distances')['PCC'].mean()
average_corr_k41.drop(index=average_corr_k41.index[0], axis=0, inplace=True)


### EVEN NONBREEDING GROUND AREAS, r selected species  ###################################################################################
######################################################################################################################
## Plotting Even vs Uneven Sized Breeding Ground Dispersal:
###FIGURE 5A 2 migraiton routes  
plt.plot(average_corr_k21_UNEVEN, linewidth=2, linestyle="dashed", color="black", label="Even Sized")
plt.plot(average_corr_k21_UNEVEN, linewidth=2, color="black", label="Uneven Sized") ## Uneven breeding ground, 2 mig routes, 1 corr


plt.plot(average_corr_k20, linestyle="dashed", linewidth=3, color="green", label="") ##
plt.plot(average_corr_k20_UNEVEN, color="green", linewidth=3, label="0 corr") ##

    
#plt.plot(average_corr_k2025, linestyle="dashed",  color ="orange",  linewidth=2, label="0.5 corr, random mig") ## E
#plt.plot(average_corr_k2025_UNEVEN, color="orange", linewidth=2, label="0.5 corr, proximity mig") ##


plt.plot(average_corr_k2050, linewidth=3, linestyle="dashed", color="slategrey",  label="") ##
plt.plot(average_corr_k4050_UNEVEN,linewidth=3,  color="slategrey", label="0.5 corr") ##


plt.plot(average_corr_k21, linewidth=3, linestyle="dashed", color="purple", label="") ## Even breeding ground, 2 mig routes, 1 corr
plt.plot(average_corr_k21_UNEVEN,linewidth=3,  color="purple", label="1 corr") ## Uneven breeding ground, 2 mig routes, 1 corr



plt.title('')
plt.xlabel('Distance')
plt.xlim([3, 50])
plt.ylim([0, 1])
plt.ylim(ymin=0)
plt.ylabel('Average Correlation')
plt.legend()
plt.legend(loc=(0.04, 0), frameon=False, title="", fontsize=12, labelspacing=0)
plt.rcParams.update({'font.size': 14})
plt.savefig('FIGURE5A_PAPERI.png', dpi = 600)
plt.close()



###FIGURE 5B 4 migraiton routes  
plt.plot(average_corr_k41_UNEVEN, linewidth=2, linestyle="dashed", color="black", label="Even Sized")
plt.plot(average_corr_k41_UNEVEN, linewidth=2, color="black", label="Uneven Sized") ## Uneven breeding ground, 2 mig routes, 1 corr


plt.plot(average_corr_k40, linestyle="dashed", color="darkgreen", linewidth=3) ## Even breeding ground, 4 mig routes, 0 corr
plt.plot(average_corr_k40_UNEVEN, color="darkgreen", linewidth=3, label="0 corr") ## Uneven breeding ground, 4 mig routes, 0 corr


#plt.plot(average_corr_k4025, linestyle="dashed",  color ="orange",  linewidth=2) ## Even breeding ground, 4 mig, 0.25 corr
#plt.plot(average_corr_k4025_UNEVEN, color="orange", linewidth=2, label="Corr 0.25") ## Uneven breeding ground, 4 mig routes, 0.25 corr


plt.plot(average_corr_k4050, linestyle="dashed",  color = "slategrey", linewidth=3)  ## Even breeding ground, 4 mig, 0.75 cor
plt.plot(average_corr_k4050_UNEVEN, color="slategrey", linewidth=3, label="0.5 corr") ## Uneven breeding ground, 4 mig routes, 0.75 corr

plt.plot(average_corr_k41, linestyle="dashed", color="purple", linewidth=3) ## Even breeding ground, 4 mig routes, 1 corr
plt.plot(average_corr_k41_UNEVEN, color="purple", linewidth=3, label="1 corr") ## Uneven breeding ground, 4 mig routes, 1 corr



plt.title('')
plt.xlabel('Distance')
plt.xlim([3, 50])
plt.ylim([0, 1])
plt.ylim(ymin=0)
plt.ylabel('Average Correlation')
plt.legend()
plt.legend(loc=(0.64, 0.67), frameon=False, title="", fontsize=12, labelspacing=0)
plt.rcParams.update({'font.size': 14})
plt.savefig('FIGURE5B_PAPERI.png', dpi = 600)
plt.close()

#############################################################################################################
########## FIGURE 6 #########################################################################################
#############################################################################################################


## No Migration K selected Species
## Import distance file
distances_1 = pd.read_csv('Midpoint_Distances_Gridsize150_Thin2.csv', header = None)
distances_1 = distances_1.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
distlist_1 = distances_1.values.tolist()
flatdist_1 = (np.array(distlist_1)).flatten()


######## WITHIN NONBREEDING CORR = 1 ##############################################################################
#######################################################################################
# 2 Mig Route: 0 correlation, K
abs_output_3 = pd.read_csv('CorrWithin1/Output_S0.9_F0.25_Corr0_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k20_withincorr1 = result_3.groupby('distances')['PCC'].mean()
average_corr_k20_withincorr1.drop(index=average_corr_k20_withincorr1.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.25 correlation, K
abs_output_4 = pd.read_csv('CorrWithin1/Output_S0.9_F0.25_Corr0.25_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_4 = {}
columns_4 = abs_output_4.columns.tolist()
for col_a, col_b in itertools.combinations(columns_4, 2):
     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
result_4.columns = ['PCC', 'p-value']
result_4['distances'] = np.array(flatdist_1)
average_corr_k2025_withincorr1 = result_4.groupby('distances')['PCC'].mean()
average_corr_k2025_withincorr1.drop(index=average_corr_k2025_withincorr1.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.75 correlation,  K
abs_output_5 = pd.read_csv('CorrWithin1/Output_S0.9_F0.25_Corr0.75_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_5 = abs_output_5.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_5 = {}
columns_5 = abs_output_5.columns.tolist()
for col_a, col_b in itertools.combinations(columns_5, 2):
     correlations_5[col_a, '__' , col_b] = stats.pearsonr(abs_output_5.loc[:, col_a], abs_output_5.loc[:, col_b])
result_5 = pd.DataFrame.from_dict(correlations_5, orient='index')
result_5.columns = ['PCC', 'p-value']
result_5['distances'] = np.array(flatdist_1)
average_corr_k2075_withincorr1 = result_5.groupby('distances')['PCC'].mean()
average_corr_k2075_withincorr1.drop(index=average_corr_k2075_withincorr1.index[0], axis=0, inplace=True)

# 2 Mig Route: 1 correlation, K
abs_output_11 = pd.read_csv('CorrWithin1/Output_S0.9_F0.25_Corr1_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_11 = abs_output_11.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_11 = {}
columns_11 = abs_output_11.columns.tolist()
for col_a, col_b in itertools.combinations(columns_11, 2):
     correlations_11[col_a, '__' , col_b] = stats.pearsonr(abs_output_11.loc[:, col_a], abs_output_11.loc[:, col_b])
result_11 = pd.DataFrame.from_dict(correlations_11, orient='index')
result_11.columns = ['PCC', 'p-value']
result_11['distances'] = np.array(flatdist_1)
average_corr_k21_withincorr1 = result_11.groupby('distances')['PCC'].mean()
average_corr_k21_withincorr1.drop(index=average_corr_k21_withincorr1.index[0], axis=0, inplace=True)


# 4 Mig Route: 0 correlation, K
abs_output_6 = pd.read_csv('CorrWithin1/Output_S0.9_F0.25_Corr0_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_6 = abs_output_6.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_6 = {}
columns_6 = abs_output_6.columns.tolist()
for col_a, col_b in itertools.combinations(columns_6, 2):
     correlations_6[col_a, '__' , col_b] = stats.pearsonr(abs_output_6.loc[:, col_a], abs_output_6.loc[:, col_b])
result_6 = pd.DataFrame.from_dict(correlations_6, orient='index')
result_6.columns = ['PCC', 'p-value']
result_6['distances'] = np.array(flatdist_1)
average_corr_k40_withincorr1 = result_6.groupby('distances')['PCC'].mean()
average_corr_k40_withincorr1.drop(index=average_corr_k40_withincorr1.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.25 correlation, K
abs_output_8 = pd.read_csv('CorrWithin1/Output_S0.9_F0.25_Corr0.25_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_8 = {}
columns_8 = abs_output_8.columns.tolist()
for col_a, col_b in itertools.combinations(columns_8, 2):
     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
result_8.columns = ['PCC', 'p-value']
result_8['distances'] = np.array(flatdist_1)
average_corr_k4025_withincorr1 = result_8.groupby('distances')['PCC'].mean()
average_corr_k4025_withincorr1.drop(index=average_corr_k4025_withincorr1.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.75 correlation, K
abs_output_9 = pd.read_csv('CorrWithin1/Output_S0.9_F0.25_Corr0.75_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_9 = abs_output_9.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_9 = {}
columns_9 = abs_output_9.columns.tolist()
for col_a, col_b in itertools.combinations(columns_9, 2):
     correlations_9[col_a, '__' , col_b] = stats.pearsonr(abs_output_9.loc[:, col_a], abs_output_9.loc[:, col_b])
result_9 = pd.DataFrame.from_dict(correlations_9, orient='index')
result_9.columns = ['PCC', 'p-value']
result_9['distances'] = np.array(flatdist_1)
average_corr_k4075_withincorr1 = result_9.groupby('distances')['PCC'].mean()
average_corr_k4075_withincorr1.drop(index=average_corr_k4075_withincorr1.index[0], axis=0, inplace=True)

# 4 Mig Route: 1 correlation, K
abs_output_7 = pd.read_csv('CorrWithin1/Output_S0.9_F0.25_Corr1_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_7 = abs_output_7.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_7 = {}
columns_7 = abs_output_7.columns.tolist()
for col_a, col_b in itertools.combinations(columns_7, 2):
     correlations_7[col_a, '__' , col_b] = stats.pearsonr(abs_output_7.loc[:, col_a], abs_output_7.loc[:, col_b])
result_7 = pd.DataFrame.from_dict(correlations_7, orient='index')
result_7.columns = ['PCC', 'p-value']
result_7['distances'] = np.array(flatdist_1)
average_corr_k41_withincorr1 = result_7.groupby('distances')['PCC'].mean()
average_corr_k41_withincorr1.drop(index=average_corr_k41_withincorr1.index[0], axis=0, inplace=True)





######## WITHIN NONBREEDING CORR = 0.75 ##############################################################################
#######################################################################################
# 2 Mig Route: 0 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k20_withincorr75 = result_3.groupby('distances')['PCC'].mean()
average_corr_k20_withincorr75.drop(index=average_corr_k20_withincorr75.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.25 correlation, K
abs_output_4 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.25_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_4 = {}
columns_4 = abs_output_4.columns.tolist()
for col_a, col_b in itertools.combinations(columns_4, 2):
     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
result_4.columns = ['PCC', 'p-value']
result_4['distances'] = np.array(flatdist_1)
average_corr_k2025_withincorr75 = result_4.groupby('distances')['PCC'].mean()
average_corr_k2025_withincorr75.drop(index=average_corr_k2025_withincorr75.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.75 correlation,  K
abs_output_5 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.75_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_5 = abs_output_5.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_5 = {}
columns_5 = abs_output_5.columns.tolist()
for col_a, col_b in itertools.combinations(columns_5, 2):
     correlations_5[col_a, '__' , col_b] = stats.pearsonr(abs_output_5.loc[:, col_a], abs_output_5.loc[:, col_b])
result_5 = pd.DataFrame.from_dict(correlations_5, orient='index')
result_5.columns = ['PCC', 'p-value']
result_5['distances'] = np.array(flatdist_1)
average_corr_k2075_withincorr75 = result_5.groupby('distances')['PCC'].mean()
average_corr_k2075_withincorr75.drop(index=average_corr_k2075_withincorr75.index[0], axis=0, inplace=True)

# 2 Mig Route: 1 correlation, K
abs_output_11 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr1_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_11 = abs_output_11.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_11 = {}
columns_11 = abs_output_11.columns.tolist()
for col_a, col_b in itertools.combinations(columns_11, 2):
     correlations_11[col_a, '__' , col_b] = stats.pearsonr(abs_output_11.loc[:, col_a], abs_output_11.loc[:, col_b])
result_11 = pd.DataFrame.from_dict(correlations_11, orient='index')
result_11.columns = ['PCC', 'p-value']
result_11['distances'] = np.array(flatdist_1)
average_corr_k21_withincorr75 = result_11.groupby('distances')['PCC'].mean()
average_corr_k21_withincorr75.drop(index=average_corr_k21_withincorr75.index[0], axis=0, inplace=True)


# 4 Mig Route: 0 correlation, K
abs_output_6 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_6 = abs_output_6.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_6 = {}
columns_6 = abs_output_6.columns.tolist()
for col_a, col_b in itertools.combinations(columns_6, 2):
     correlations_6[col_a, '__' , col_b] = stats.pearsonr(abs_output_6.loc[:, col_a], abs_output_6.loc[:, col_b])
result_6 = pd.DataFrame.from_dict(correlations_6, orient='index')
result_6.columns = ['PCC', 'p-value']
result_6['distances'] = np.array(flatdist_1)
average_corr_k40_withincorr75 = result_6.groupby('distances')['PCC'].mean()
average_corr_k40_withincorr75.drop(index=average_corr_k40_withincorr75.index[0], axis=0, inplace=True)

# 4 Mig Route: 0.25 correlation, K
abs_output_8 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.25_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_8 = {}
columns_8 = abs_output_8.columns.tolist()
for col_a, col_b in itertools.combinations(columns_8, 2):
     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
result_8.columns = ['PCC', 'p-value']
result_8['distances'] = np.array(flatdist_1)
average_corr_k4025_withincorr75 = result_8.groupby('distances')['PCC'].mean()
average_corr_k4025_withincorr75.drop(index=average_corr_k4025_withincorr75.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.75 correlation, K
abs_output_9 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr0.75_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_9 = abs_output_9.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_9 = {}
columns_9 = abs_output_9.columns.tolist()
for col_a, col_b in itertools.combinations(columns_9, 2):
     correlations_9[col_a, '__' , col_b] = stats.pearsonr(abs_output_9.loc[:, col_a], abs_output_9.loc[:, col_b])
result_9 = pd.DataFrame.from_dict(correlations_9, orient='index')
result_9.columns = ['PCC', 'p-value']
result_9['distances'] = np.array(flatdist_1)
average_corr_k4075_withincorr75 = result_9.groupby('distances')['PCC'].mean()
average_corr_k4075_withincorr75.drop(index=average_corr_k4075_withincorr75.index[0], axis=0, inplace=True)

# 4 Mig Route: 1 correlation, K
abs_output_7 = pd.read_csv('CorrWithin0.75/Output_S0.9_F0.25_Corr1_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_7 = abs_output_7.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_7 = {}
columns_7 = abs_output_7.columns.tolist()
for col_a, col_b in itertools.combinations(columns_7, 2):
     correlations_7[col_a, '__' , col_b] = stats.pearsonr(abs_output_7.loc[:, col_a], abs_output_7.loc[:, col_b])
result_7 = pd.DataFrame.from_dict(correlations_7, orient='index')
result_7.columns = ['PCC', 'p-value']
result_7['distances'] = np.array(flatdist_1)
average_corr_k41_withincorr75 = result_7.groupby('distances')['PCC'].mean()
average_corr_k41_withincorr75.drop(index=average_corr_k41_withincorr75.index[0], axis=0, inplace=True)



######## WITHIN NONBREEDING CORR = 0 ##################################################
#######################################################################################
# 2 Mig Route: 0 correlation, K
abs_output_3 = pd.read_csv('CorrWithin0/Output_S0.9_F0.25_Corr0_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_3 = abs_output_3.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_3 = {}
columns_3 = abs_output_3.columns.tolist()
for col_a, col_b in itertools.combinations(columns_3, 2):
     correlations_3[col_a, '__' , col_b] = stats.pearsonr(abs_output_3.loc[:, col_a], abs_output_3.loc[:, col_b])
result_3 = pd.DataFrame.from_dict(correlations_3, orient='index')
result_3.columns = ['PCC', 'p-value']
result_3['distances'] = np.array(flatdist_1)
average_corr_k20_withincorr0 = result_3.groupby('distances')['PCC'].mean()
average_corr_k20_withincorr0.drop(index=average_corr_k20_withincorr0.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.25 correlation, K
abs_output_4 = pd.read_csv('CorrWithin0/Output_S0.9_F0.25_Corr0.25_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_4 = abs_output_4.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_4 = {}
columns_4 = abs_output_4.columns.tolist()
for col_a, col_b in itertools.combinations(columns_4, 2):
     correlations_4[col_a, '__' , col_b] = stats.pearsonr(abs_output_4.loc[:, col_a], abs_output_4.loc[:, col_b])
result_4 = pd.DataFrame.from_dict(correlations_4, orient='index')
result_4.columns = ['PCC', 'p-value']
result_4['distances'] = np.array(flatdist_1)
average_corr_k2025_withincorr0 = result_4.groupby('distances')['PCC'].mean()
average_corr_k2025_withincorr0.drop(index=average_corr_k2025_withincorr0.index[0], axis=0, inplace=True)


# 2 Mig Route: 0.75 correlation,  K
abs_output_5 = pd.read_csv('CorrWithin0/Output_S0.9_F0.25_Corr0.75_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_5 = abs_output_5.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_5 = {}
columns_5 = abs_output_5.columns.tolist()
for col_a, col_b in itertools.combinations(columns_5, 2):
     correlations_5[col_a, '__' , col_b] = stats.pearsonr(abs_output_5.loc[:, col_a], abs_output_5.loc[:, col_b])
result_5 = pd.DataFrame.from_dict(correlations_5, orient='index')
result_5.columns = ['PCC', 'p-value']
result_5['distances'] = np.array(flatdist_1)
average_corr_k2075_withincorr0 = result_5.groupby('distances')['PCC'].mean()
average_corr_k2075_withincorr0.drop(index=average_corr_k2075_withincorr0.index[0], axis=0, inplace=True)

# 2 Mig Route: 1 correlation, K
abs_output_11 = pd.read_csv('CorrWithin0/Output_S0.9_F0.25_Corr1_Migroutes2_Nit500_Carryover1.csv', header = None)
abs_output_11 = abs_output_11.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_11 = {}
columns_11 = abs_output_11.columns.tolist()
for col_a, col_b in itertools.combinations(columns_11, 2):
     correlations_11[col_a, '__' , col_b] = stats.pearsonr(abs_output_11.loc[:, col_a], abs_output_11.loc[:, col_b])
result_11 = pd.DataFrame.from_dict(correlations_11, orient='index')
result_11.columns = ['PCC', 'p-value']
result_11['distances'] = np.array(flatdist_1)
average_corr_k21_withincorr0 = result_11.groupby('distances')['PCC'].mean()
average_corr_k21_withincorr0.drop(index=average_corr_k21_withincorr0.index[0], axis=0, inplace=True)


# 4 Mig Route: 0 correlation, K
abs_output_6 = pd.read_csv('CorrWithin0/Output_S0.9_F0.25_Corr0_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_6 = abs_output_6.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_6 = {}
columns_6 = abs_output_6.columns.tolist()
for col_a, col_b in itertools.combinations(columns_6, 2):
     correlations_6[col_a, '__' , col_b] = stats.pearsonr(abs_output_6.loc[:, col_a], abs_output_6.loc[:, col_b])
result_6 = pd.DataFrame.from_dict(correlations_6, orient='index')
result_6.columns = ['PCC', 'p-value']
result_6['distances'] = np.array(flatdist_1)
average_corr_k40_withincorr0 = result_6.groupby('distances')['PCC'].mean()
average_corr_k40_withincorr0.drop(index=average_corr_k40_withincorr0.index[0], axis=0, inplace=True)

# 4 Mig Route: 0.25 correlation, K
abs_output_8 = pd.read_csv('CorrWithin0/Output_S0.9_F0.25_Corr0.25_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_8 = abs_output_8.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_8 = {}
columns_8 = abs_output_8.columns.tolist()
for col_a, col_b in itertools.combinations(columns_8, 2):
     correlations_8[col_a, '__' , col_b] = stats.pearsonr(abs_output_8.loc[:, col_a], abs_output_8.loc[:, col_b])
result_8 = pd.DataFrame.from_dict(correlations_8, orient='index')
result_8.columns = ['PCC', 'p-value']
result_8['distances'] = np.array(flatdist_1)
average_corr_k4025_withincorr0 = result_8.groupby('distances')['PCC'].mean()
average_corr_k4025_withincorr0.drop(index=average_corr_k4025_withincorr0.index[0], axis=0, inplace=True)


# 4 Mig Route: 0.75 correlation, K
abs_output_9 = pd.read_csv('CorrWithin0/Output_S0.9_F0.25_Corr0.75_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_9 = abs_output_9.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_9 = {}
columns_9 = abs_output_9.columns.tolist()
for col_a, col_b in itertools.combinations(columns_9, 2):
     correlations_9[col_a, '__' , col_b] = stats.pearsonr(abs_output_9.loc[:, col_a], abs_output_9.loc[:, col_b])
result_9 = pd.DataFrame.from_dict(correlations_9, orient='index')
result_9.columns = ['PCC', 'p-value']
result_9['distances'] = np.array(flatdist_1)
average_corr_k4075_withincorr0 = result_9.groupby('distances')['PCC'].mean()
average_corr_k4075_withincorr0.drop(index=average_corr_k4075_withincorr0.index[0], axis=0, inplace=True)

# 4 Mig Route: 1 correlation, K
abs_output_7 = pd.read_csv('CorrWithin0/Output_S0.9_F0.25_Corr1_Migroutes4_Nit500_Carryover1.csv', header = None)
abs_output_7 = abs_output_7.astype(str).replace({"\[":"", "\]":""}, regex=True).astype(float)
correlations_7 = {}
columns_7 = abs_output_7.columns.tolist()
for col_a, col_b in itertools.combinations(columns_7, 2):
     correlations_7[col_a, '__' , col_b] = stats.pearsonr(abs_output_7.loc[:, col_a], abs_output_7.loc[:, col_b])
result_7 = pd.DataFrame.from_dict(correlations_7, orient='index')
result_7.columns = ['PCC', 'p-value']
result_7['distances'] = np.array(flatdist_1)
average_corr_k41_withincorr0 = result_7.groupby('distances')['PCC'].mean()
average_corr_k41_withincorr0.drop(index=average_corr_k41_withincorr0.index[0], axis=0, inplace=True)






###FIGURE 6 2 migraiton routes  
x = np.arange(0,63,1)

#plt.plot(average_corr_k20_withincorr1, linestyle="solid", linewidth=3, color="red", label="") ##
#plt.plot(average_corr_k2025_withincorr1, linestyle="dashed", linewidth=1, color="pink", label="") ##
#plt.plot(average_corr_k2075_withincorr1, linestyle="dashed", linewidth=1, color="pink", label="") ##
#plt.plot(average_corr_k21_withincorr1, linestyle="solid", linewidth=3, color="red", label="") ##
y1 = pd.DataFrame(average_corr_k20_withincorr1)
y2 = pd.DataFrame(average_corr_k2025_withincorr1)
y3 = pd.DataFrame(average_corr_k2075_withincorr1)
y4 = pd.DataFrame(average_corr_k21_withincorr1)

plt.fill_between(x, y1['PCC'], y4['PCC'], color='red',alpha=0.2, label="1 corr within")

#plt.plot(average_corr_k20_withincorr75, linestyle="solid", linewidth=3, color="grey", label="") ##
#plt.plot(average_corr_k2025_withincorr75, linestyle="dashed", linewidth=1, color="grey", label="") ##
#plt.plot(average_corr_k2075_withincorr75, linestyle="dashed", linewidth=1, color="grey", label="") ##
#plt.plot(average_corr_k21_withincorr75, linestyle="solid", linewidth=3, color="grey", label="") ##
y5 = pd.DataFrame(average_corr_k20_withincorr75)
y6 = pd.DataFrame(average_corr_k2025_withincorr75)
y7 = pd.DataFrame(average_corr_k2075_withincorr75)
y8 = pd.DataFrame(average_corr_k21_withincorr75)

plt.fill_between(range(len(x)), y5['PCC'], y8['PCC'],color='grey',alpha=0.2, label="0.75 corr within")


#plt.plot(average_corr_k20_withincorr0, linestyle="solid", linewidth=3, color="blue", label="") ##
#plt.plot(average_corr_k2025_withincorr0, linestyle="dashed", linewidth=1, color="lightblue", label="") ##
#plt.plot(average_corr_k2075_withincorr0, linestyle="dashed", linewidth=1, color="lightblue", label="") ##
#plt.plot(average_corr_k21_withincorr0, linestyle="solid", linewidth=3, color="blue", label="test") ##
y9 = pd.DataFrame(average_corr_k20_withincorr0)
y10 = pd.DataFrame(average_corr_k2025_withincorr0)
y11 = pd.DataFrame(average_corr_k2075_withincorr0)
y12 = pd.DataFrame(average_corr_k21_withincorr0)

plt.fill_between(range(len(x)), y9['PCC'], y12['PCC'],color='orange',alpha=0.2, label="0 corr within")



plt.title('')
plt.xlabel('Distance')
plt.xlim([3, 50])
plt.ylim([0, 1])
plt.ylim(ymin=0)
plt.ylabel('Average Correlation')
plt.legend()
plt.legend(loc=(1, 0.8), frameon=False, title="", fontsize=12, labelspacing=0)
plt.rcParams.update({'font.size': 14})


plt.savefig('FIGURE6_PAPERI.png', dpi = 800)
plt.close()


