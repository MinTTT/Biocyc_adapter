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
# […]

# Own modules
#%%


data_dir = r'Y:\chupan\fulab_zc_1\seq_data\202312_RNA_seq'

output_fold = [folder for folder in os.listdir(data_dir) if folder.split('_')[-1] == 'output']

statistic_dict = {}
htseq_dict = {}
for fold in output_fold:
    path = os.path.join(data_dir, fold)
    for file in os.listdir(path):
        if file.endswith('csv'):
            if file.split('.')[-2] == 'expression_statistic':
                statistic_dict[file.split('.')[0]] = pd.read_csv(os.path.join(path, file))
                htseq_dict[file.split('.')[0]] = pd.read_csv(os.path.join(path, 'htseq.cds_counts.tsv'), delimiter='\t')


#%%
# effective reads number
effective_summary = []
for key, data in statistic_dict.items():
    effective_reads = data['htseq_counts'].sum()
    effective_summary.append([key, effective_reads])
effective_summary = pd.DataFrame(effective_summary)

#%%
# all htseq-TPM
gene_number = len(statistic_dict['1126_2L'])
sample_number = len(statistic_dict.keys())
data_all = np.zeros((gene_number, sample_number))
i = 0
for key, data in statistic_dict.items():
    data_all[:, i] = data['htseq_TPM'].to_numpy()
    i += 1
#%%
# all TPM
gene_number = len(statistic_dict['1126_2L'])
sample_number = len(statistic_dict.keys())
data_all = np.zeros((gene_number, sample_number))
i = 0
for key, data in statistic_dict.items():
    data_all[:, i] = data['TPM'].to_numpy()
    i += 1