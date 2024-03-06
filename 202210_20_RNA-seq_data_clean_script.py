
"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
import io
# Built-in/Generic Imports
import os
# […]

# Libs
import pandas as pd
import numpy as np  # Or any other
from matplotlib import pyplot as plt

# […]

# Own modules
# %%


data_dir = r'Y:\chupan\fulab_zc_1\seq_data\2022_20_RNA_seq'

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

gene_df = statistic_dict[list(statistic_dict.keys())[0]][['gene', 'ECOCYC']]

#%%
# all htseq-TPM
gene_number = len(statistic_dict[list(statistic_dict.keys())[0]])
gene_list = statistic_dict[list(statistic_dict.keys())[0]]['gene'].to_list()
sample_number = len(statistic_dict.keys())
sample_list = list(statistic_dict.keys())
sample_list.sort(key=lambda val: int(val.split('_')[0]))
htTPM_data_all = np.zeros((gene_number, sample_number))


for i, key in enumerate(sample_list):
    htTPM_data_all[:, i] = statistic_dict[key]['htseq_TPM'].to_numpy()

htTPM_df = gene_df.copy(deep=True)
htTPM_df[sample_list] = htTPM_data_all

# %%
# all TPM
gene_number = len(statistic_dict[list(statistic_dict.keys())[0]])
sample_number = len(statistic_dict.keys())
TPM_data_all = np.zeros((gene_number, sample_number))
i = 0
for key, data in statistic_dict.items():
    TPM_data_all[:, i] = data['TPM'].to_numpy()
    i += 1

# %%
# all reads number
gene_number = len(statistic_dict[list(statistic_dict.keys())[0]])
sample_number = len(statistic_dict.keys())
reads_data_all = np.zeros((gene_number, sample_number))
i = 0
sample_name_list = []
for key, data in statistic_dict.items():
    reads_data_all[:, i] = data['htseq_counts'].to_numpy()
    sample_name_list.append(key)
    i += 1
# %% sample info

sample_info = """
sample_name	sample	strain	media	Condition	growth_rate	growth_rate_std
4041_L03	4041	NH3.24	RDM Glucose	8	1.58760375	0.022712675
4041_L04	4041	NH3.24	RDM Glucose	8	1.58760375	0.022712675
4042_L03	4042	NH3.23	RDM Glucose	1	1.52532975	0.028251071
4042_L04	4042	NH3.23	RDM Glucose	1	1.52532975	0.028251071
4043_L03	4043	NH3.23	RDM Glucose	1	1.52532975	0.028251071
4043_L04	4043	NH3.23	RDM Glucose	1	1.52532975	0.028251071
4044_L03	4044	NH3.23	RDM Glucose	1	1.52532975	0.028251071
4044_L04	4044	NH3.23	RDM Glucose	1	1.52532975	0.028251071
4045_L03	4045	NH3.24	RDM Glucose	8	1.58760375	0.022712675
4045_L04	4045	NH3.24	RDM Glucose	8	1.58760375	0.022712675
4046_L03	4046	NH3.24	RDM Glucose	8	1.58760375	0.022712675
4046_L04	4046	NH3.24	RDM Glucose	8	1.58760375	0.022712675
4047_L03	4047	NH3.23	RDM Glycerol	2	1.217752	0.015268
4047_L04	4047	NH3.23	RDM Glycerol	2	1.217752	0.015268
4048_L03	4048	NH3.23	RDM Glycerol	2	1.217752	0.015268
4048_L04	4048	NH3.23	RDM Glycerol	2	1.217752	0.015268
4049_L03	4049	NH3.23	RDM Glycerol	2	1.217752	0.015268
4049_L04	4049	NH3.23	RDM Glycerol	2	1.217752	0.015268
4082_L03	4082	NH3.23	MOPS Glucose	4	0.8843265	0.011667561
4082_L04	4082	NH3.23	MOPS Glucose	4	0.8843265	0.011667561
4083_L03	4083	NH3.23	MOPS Glucose	4	0.8843265	0.011667561
4083_L04	4083	NH3.23	MOPS Glucose	4	0.8843265	0.011667561
4084_L03	4084	NH3.24	MOPS Glucose	11	0.9680055	0.01376398
4084_L04	4084	NH3.24	MOPS Glucose	11	0.9680055	0.01376398
4086_L03	4086	NH3.24	MOPS Glucose	11	0.9680055	0.01376398
4086_L04	4086	NH3.24	MOPS Glucose	11	0.9680055	0.01376398
4087_L03	4087	NH3.23	MOPS Acetate	6	0.40619375	0.011502244
4087_L04	4087	NH3.23	MOPS Acetate	6	0.40619375	0.011502244
4088_L03	4088	NH3.23	MOPS Acetate	6	0.40619375	0.011502244
4088_L04	4088	NH3.23	MOPS Acetate	6	0.40619375	0.011502244
4089_L03	4089	NH3.23	MOPS Acetate	6	0.40619375	0.011502244
4089_L04	4089	NH3.23	MOPS Acetate	6	0.40619375	0.011502244
4111_L03	4111	NH3.23	MOPS CAA Glycerol	3	1.126825333	0.016655681
4111_L04	4111	NH3.23	MOPS CAA Glycerol	3	1.126825333	0.016655681
4112_L03	4112	NH3.23	MOPS CAA Glycerol	3	1.126825333	0.016655681
4112_L04	4112	NH3.23	MOPS CAA Glycerol	3	1.126825333	0.016655681
4113_L03	4113	NH3.23	MOPS CAA Glycerol	3	1.126825333	0.016655681
4113_L04	4113	NH3.23	MOPS CAA Glycerol	3	1.126825333	0.016655681
4114_L03	4114	NH3.24	MOPS CAA Glycerol	10	1.150488	0.024666912
4114_L04	4114	NH3.24	MOPS CAA Glycerol	10	1.150488	0.024666912
4115_L03	4115	NH3.24	MOPS CAA Glycerol	10	1.150488	0.024666912
4115_L04	4115	NH3.24	MOPS CAA Glycerol	10	1.150488	0.024666912
4116_L03	4116	NH3.24	MOPS CAA Glycerol	10	1.150488	0.024666912
4116_L04	4116	NH3.24	MOPS CAA Glycerol	10	1.150488	0.024666912
4117_L03	4117	NH3.23	MOPS Glycerol	5	0.624194333	0.010284432
4117_L04	4117	NH3.23	MOPS Glycerol	5	0.624194333	0.010284432
4118_L03	4118	NH3.23	MOPS Glycerol	5	0.624194333	0.010284432
4118_L04	4118	NH3.23	MOPS Glycerol	5	0.624194333	0.010284432
4119_L03	4119	NH3.23	MOPS Glycerol	5	0.624194333	0.010284432
4119_L04	4119	NH3.23	MOPS Glycerol	5	0.624194333	0.010284432
40411_L03	40411	NH3.24	RDM Glycerol	9	1.281985	0.01831
40411_L04	40411	NH3.24	RDM Glycerol	9	1.281985	0.01831
40412_L03	40412	NH3.24	RDM Glycerol	9	1.281985	0.01831
40412_L04	40412	NH3.24	RDM Glycerol	9	1.281985	0.01831
40810_L03	40810	NH3.24	MOPS Acetate	13	0.38972975	0.004989396
40810_L04	40810	NH3.24	MOPS Acetate	13	0.38972975	0.004989396
40811_L03	40811	NH3.24	MOPS Acetate	13	0.38972975	0.004989396
40811_L04	40811	NH3.24	MOPS Acetate	13	0.38972975	0.004989396
40812_L03	40812	NH3.24	MOPS Acetate	13	0.38972975	0.004989396
40812_L04	40812	NH3.24	MOPS Acetate	13	0.38972975	0.004989396
40813_L03	40813	NH3.23	MOPS Glutamate Glucose	7	0.18176	0.005815633
40813_L04	40813	NH3.23	MOPS Glutamate Glucose	7	0.18176	0.005815633
40815_L03	40815	NH3.23	MOPS Glutamate Glucose	7	0.18176	0.005815633
40815_L04	40815	NH3.23	MOPS Glutamate Glucose	7	0.18176	0.005815633
41110_L03	41110	NH3.24	MOPS Glycerol	12	0.67794775	0.006927106
41110_L04	41110	NH3.24	MOPS Glycerol	12	0.67794775	0.006927106
41111_L03	41111	NH3.24	MOPS Glycerol	12	0.67794775	0.006927106
41111_L04	41111	NH3.24	MOPS Glycerol	12	0.67794775	0.006927106
41112_L03	41112	NH3.24	MOPS Glycerol	12	0.67794775	0.006927106
41112_L04	41112	NH3.24	MOPS Glycerol	12	0.67794775	0.006927106
"""
import io
sample_info = pd.read_table(io.StringIO(sample_info), delimiter='\t')

condition_set = list(set(sample_info['Condition'].tolist()))
condition_set.sort()

avg_data = np.empty((gene_number, len(condition_set)))
std_data = np.empty((gene_number, len(condition_set)))
cv_data = np.empty((gene_number, len(condition_set)))
for s_i, cdtion in enumerate(condition_set):
    sample_name = sample_info[sample_info['Condition'] == cdtion]['sample_name']

    sample_data = htTPM_df.loc[:, sample_name]
    sample_avg = np.mean(sample_data, axis=1)
    sample_std = np.std(sample_data, axis=1)
    avg_data[:, s_i] = sample_avg
    std_data[:, s_i] = sample_std
    cv_data[:, s_i] = sample_std / sample_avg

avg_zscore = (avg_data - np.average(avg_data, axis=1).reshape(-1, 1)) / np.std(avg_data, axis=1).reshape(-1, 1)
avg_zscore_df = pd.DataFrame(data=np.hstack((gene_df['gene'].to_numpy().reshape(-1, 1), avg_zscore)),
                             columns=['gene'] + condition_set)
avg_tpm_df = pd.DataFrame(data=np.hstack((gene_df['gene'].to_numpy().reshape(-1, 1), avg_data)),
                          columns=['gene'] + condition_set)

std_tpm_df = pd.DataFrame(data=np.hstack((gene_df['gene'].to_numpy().reshape(-1, 1), std_data)),
                          columns=['gene'] + condition_set)

# plot cv
fig, ax = plt.subplots()
img = ax.imshow(cv_data, interpolation='nearest', cmap='jet', aspect='auto', vmin=0, vmax=.3)
plt.colorbar(img)
fig.tight_layout()
fig.show()
# %%  statistic genes in toggle
condition_set = list(set(sample_info['Condition'].tolist()))

condition_set.sort()

toggle_gene_list = ['LacI', 'tetR', 'mCherry', 'aph%283%27%29-II %28or nptII%29', 'GFP']
R_filter = gene_df['gene'].isin(toggle_gene_list)
R_hist_num = np.sum(R_filter)
R_sector = htTPM_data_all[R_filter]

R_data = np.empty((R_hist_num, len(condition_set)))
for s_i, sample in enumerate(condition_set):
    sample_mask = sample_info['Condition'] == sample
    sample_data = R_sector[:, sample_mask]
    sample_avg = np.mean(sample_data, axis=1)
    R_data[:, s_i] = sample_avg

R_total = np.sum(R_data, axis=0)
R_total = pd.DataFrame(data=np.array([condition_set, R_total, R_total / 1e6]).T,
                       columns=['condition', 'TPM_total', 'psi'])

#%%
import json
with open(r'./exported_data/GOterm_hui_MSB_2015.json') as jfile:
    hui_GO = json.load(jfile)

condition_set = list(set(sample_info['sample'].tolist()))
condition_set.sort()
toggle_gene_list = hui_GO['R']['name']
# toggle_gene_list = ['LacI', 'tetR', 'mCherry', 'aph%283%27%29-II %28or nptII%29', 'GFP']
# toggle_gene_list = ['LacI',  'aph%283%27%29-II %28or nptII%29', 'GFP']
# toggle_gene_list = ['tetR', 'mCherry']
R_filter = gene_df['gene'].isin(toggle_gene_list)
R_hist_num = np.sum(R_filter)
R_sector = htTPM_df.loc[R_filter, :]

R_data = np.empty((R_hist_num, len(condition_set)))
for s_i, sample in enumerate(condition_set):
    samples = sample_info[sample_info['sample'] == sample]['sample_name']

    sample_data = R_sector.loc[:, samples]
    sample_avg = np.mean(sample_data, axis=1)
    R_data[:, s_i] = sample_avg

R_total = np.sum(R_data, axis=0)
R_total = pd.DataFrame(data=np.array([condition_set, R_total, R_total / 1e6]).T,
                       columns=['sample', 'TPM_total', 'psi'])

#%%
for i in range(13):
    print(sample_info[sample_info['Condition'] == i+1]['media'].tolist()[0])