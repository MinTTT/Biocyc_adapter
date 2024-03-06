# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
import sys
# […]

# Libs
import pandas as pd
import numpy as np  # Or any other
from matplotlib import pyplot as plt
# […]
#%%
# Own modules
data_path = r'D:\OneDrive\zjw_data\colony\data\omic_data\2023116-from SA-RNAseq different growth rate-TPM-withMultiPars.xlsx'

df = pd.read_excel(data_path, sheet_name='Sheet1-GR已更正')
#%%
sample_list = ['MOPS_20mM_mannose_1', 'MOPS_20mM_mannose_2', 'MOPS_20mM_mannose_3', 'MOPS_5mM_mannose_1', 'MOPS_5mM_mannose_2', 'MOPS_5mM_mannose_3', 'MOPS_fructose_1', 'MOPS_fructose_2', 'MOPS_fructose_3', 'MOPS_gluconate_1', 'MOPS_gluconate_2', 'MOPS_gluconate_3', 'MOPS_glucose_1', 'MOPS_glucose_2', 'MOPS_glucose_3', 'MOPS_glucose_AUCG_1', 'MOPS_glucose_AUCG_2', 'MOPS_glucose_AUCG_3', 'MOPS_glucose_caa_1', 'MOPS_glucose_caa_2', 'MOPS_glucose_caa_3', 'MOPS_glucose_EZ_1', 'MOPS_glucose_EZ_2', 'MOPS_glucose_EZ_3', 'MOPS_glycerol_1', 'MOPS_glycerol_2', 'MOPS_glycerol_3', 'MOPS_glycerol_AUCG_1', 'MOPS_glycerol_AUCG_2', 'MOPS_glycerol_AUCG_3', 'MOPS_glycerol_CAA_1', 'MOPS_glycerol_CAA_2', 'MOPS_glycerol_CAA_3', 'MOPS_glycerol_EZ_1', 'MOPS_glycerol_EZ_2', 'MOPS_glycerol_EZ_3', 'MOPS_sodium_acetate_1', 'MOPS_sodium_acetate_2', 'MOPS_sodium_acetate_3', 'MOPS_succinate_1', 'MOPS_succinate_2', 'MOPS_succinate_3', 'RDM_gluconate_1', 'RDM_gluconate_2', 'RDM_gluconate_3', 'RDM_glucose_1', 'RDM_glucose_2', 'RDM_glucose_3', 'RDM_glycerol_1', 'RDM_glycerol_2', 'RDM_glycerol_3', 'RDM_20mM mannose_1', 'RDM_20mM mannose_2', 'RDM_20mM mannose_3']

conditions_set = list(set(['_'.join(smp.split('_')[:-1]) for smp in sample_list]))

avg_data = np.empty((len(df), len(conditions_set)))
std_data = np.empty((len(df), len(conditions_set)))
cv_data = np.empty((len(df), len(conditions_set)))

for i, condition in enumerate(conditions_set):
    data_temp = []
    for j in range(5):
        sample = condition + f'_{j}'
        if sample in sample_list:
            data_temp.append(df[sample].to_numpy().reshape(-1, 1))
    data_temp = np.hstack(data_temp)
    avg_data[:, i] = np.mean(data_temp, axis=1)
    std_data[:, i] = np.std(data_temp, axis=1)
    cv_data[:, i] = std_data[:, i] /  avg_data[:, i]

fig, ax = plt.subplots()
img = ax.imshow(cv_data, interpolation='nearest', cmap='jet', aspect='auto', vmin=0, vmax=.3)
plt.colorbar(img)
fig.tight_layout()
fig.show()
