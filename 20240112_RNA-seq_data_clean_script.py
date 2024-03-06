# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
# […]

# Libs
import pandas as pd
import numpy as np  # Or any other

# […]

# Own modules
# %%


data_dir = r'Z:\home_port\fulab\data'

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

# %%
# effective reads number
effective_summary = []
for key, data in statistic_dict.items():
    effective_reads = data['htseq_counts'].sum()
    effective_summary.append([key, effective_reads])
effective_summary = pd.DataFrame(effective_summary)

"""
Lib1221_10
Lib1221_11
Lib1221_12
Lib1221_15
Lib1221_1
Lib1221_2
Lib1221_3
Lib1221_4
Lib1221_5
Lib1221_6
Lib1221_7
Lib1221_9
SRR17762338
"""

# %%
# all htseq-TPM
gene_number = len(statistic_dict[list(statistic_dict.keys())[0]])
sample_number = len(statistic_dict.keys())
htTPM_data_all = np.zeros((gene_number, sample_number))
i = 0
for key, data in statistic_dict.items():
    htTPM_data_all[:, i] = data['htseq_TPM'].to_numpy()
    i += 1
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

# %%
# from scipy.stats import ttest_rel
# from statsmodels.stats.multitest import fdrcorrection
# t_test = ttest_rel(htTPM_data_all[:, [4, 5]],
#                     htTPM_data_all[:, [6, 7]], axis=1, )
# p_value = t_test.pvalue
#
# q_value = fdrcorrection(p_value)

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
# from statsmodels.stats.multitest import fdrcorrection, multipletests
from pthr_go_annots import get_request_obj
import json
import matplotlib.pyplot as plt

# %% test Red green init in steady state
counsts_df = pd.DataFrame(data=reads_data_all[:, [4, 5, 6, 7]].T.astype(int))
meta = pd.DataFrame(data=[['s1', 'A'], ['s2', 'A'], ['s3', 'B'], ['s4', 'B']], columns=['sample', 'C'])
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    counts=counsts_df,
    metadata=meta,
    design_factors='C',
    refit_cooks=True,
    inference=inference,
)
dds.deseq2()

stat_res = DeseqStats(dds, inference=inference)
stat_res.summary()
stat_summary = stat_res.results_df

# %% test Red green init in steady state
counsts_df = pd.DataFrame(data=reads_data_all[:, [4, 5, 9, 10]].T.astype(int))
meta = pd.DataFrame(data=[['s1', 'A'], ['s2', 'A'], ['s3', 'B'], ['s4', 'B']], columns=['sample', 'C'])
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    counts=counsts_df,
    metadata=meta,
    design_factors='C',
    refit_cooks=True,
    inference=inference, )
dds.deseq2()

stat_res = DeseqStats(dds, inference=inference)
stat_res.summary()
stat_summary = stat_res.results_df

stat_summary['-log10FDR'] = - np.log10(stat_summary['padj'])

# %% all data
group = """A
C
B
B
1
1
2
2
3
3
4
A
O"""
group = group.split('\n')
meta = pd.DataFrame(data=dict(sample=range(len(group)), C=group))
counsts_df = pd.DataFrame(data=reads_data_all.T.astype(int))
inference = DefaultInference(n_cpus=12)
dds = DeseqDataSet(
    counts=counsts_df,
    metadata=meta,
    design_factors='C',
    refit_cooks=True,
    inference=inference,
)
dds.deseq2()

stat_res = DeseqStats(dds, inference=inference, contrast=['C', 'A', '1'])
stat_res.summary()
stat_summary = stat_res.results_df
stat_summary['-log10FDR'] = - np.log10(stat_summary['padj'])

thre_logFDR = -np.log10(0.01)
thre_log2fold = np.log2(2)

gene_df = statistic_dict[list(statistic_dict.keys())[0]][['gene', 'ECOCYC']]
stat_summary[['gene', 'ECOCYC']] = gene_df.to_numpy()
up_mask = np.logical_and(stat_summary['log2FoldChange'] >= thre_log2fold, stat_summary['-log10FDR'] >= thre_logFDR)
down_mask = np.logical_and(stat_summary['log2FoldChange'] <= -thre_log2fold, stat_summary['-log10FDR'] >= thre_logFDR)

Avs1_up_name_list = gene_df.loc[up_mask.to_list(), :]
Avs1_down_name_list = gene_df.loc[down_mask.to_list(), :]

# %%

stat_res = DeseqStats(dds, inference=inference, contrast=['C', 'B', '1'])
stat_res.summary()
stat_summary = stat_res.results_df
stat_summary['-log10FDR'] = - np.log10(stat_summary['padj'])

thre_logFDR = -np.log10(0.01)
thre_log2fold = np.log2(2)

gene_df = statistic_dict[list(statistic_dict.keys())[0]][['gene', 'ECOCYC']]
stat_summary[['gene', 'ECOCYC']] = gene_df.to_numpy()
up_mask = np.logical_and(stat_summary['log2FoldChange'] >= thre_log2fold, stat_summary['-log10FDR'] >= thre_logFDR)
down_mask = np.logical_and(stat_summary['log2FoldChange'] <= -thre_log2fold, stat_summary['-log10FDR'] >= thre_logFDR)

Bvs1_up_name_list = gene_df.loc[up_mask.to_list(), :]
Bvs1_down_name_list = gene_df.loc[down_mask.to_list(), :]

# %%

stat_res = DeseqStats(dds, inference=inference, contrast=['C', 'C', '1'])
stat_res.summary()
stat_summary = stat_res.results_df
stat_summary['-log10FDR'] = - np.log10(stat_summary['padj'])

thre_logFDR = -np.log10(0.01)
thre_log2fold = np.log2(2)

gene_df = statistic_dict[list(statistic_dict.keys())[0]][['gene', 'ECOCYC']]
stat_summary[['gene', 'ECOCYC']] = gene_df.to_numpy()
up_mask = np.logical_and(stat_summary['log2FoldChange'] >= thre_log2fold, stat_summary['-log10FDR'] >= thre_logFDR)
down_mask = np.logical_and(stat_summary['log2FoldChange'] <= -thre_log2fold, stat_summary['-log10FDR'] >= thre_logFDR)

Cvs1_up_name_list = gene_df.loc[up_mask.to_list(), :]
Cvs1_down_name_list = gene_df.loc[down_mask.to_list(), :]

# %%

stat_res = DeseqStats(dds, inference=inference, contrast=['C', '3', '1'])
stat_res.summary()
stat_summary = stat_res.results_df
stat_summary['-log10FDR'] = - np.log10(stat_summary['padj'])

thre_logFDR = -np.log10(0.01)
thre_log2fold = np.log2(2)

gene_df = statistic_dict[list(statistic_dict.keys())[0]][['gene', 'ECOCYC']]
stat_summary[['gene', 'ECOCYC']] = gene_df.to_numpy()
up_mask = np.logical_and(stat_summary['log2FoldChange'] >= thre_log2fold, stat_summary['-log10FDR'] >= thre_logFDR)
down_mask = np.logical_and(stat_summary['log2FoldChange'] <= -thre_log2fold, stat_summary['-log10FDR'] >= thre_logFDR)

C3vs1_up_name_list = gene_df.loc[up_mask.to_list(), :]
C3vs1_down_name_list = gene_df.loc[down_mask.to_list(), :]

# %% load growth law GO

with open(r'./exported_data/GOterm_hui_MSB_2015.json') as jfile:
    hui_GO = json.load(jfile)

# %%

condition_set = list(set(meta['C'].tolist()))
condition_set.sort()

avg_data = np.empty((gene_number, len(condition_set)))
for s_i, sample in enumerate(condition_set):
    sample_mask = meta['C'] == sample
    sample_data = htTPM_data_all[:, sample_mask]
    sample_avg = np.mean(sample_data, axis=1)
    avg_data[:, s_i] = sample_avg

avg_zscore = (avg_data - np.average(avg_data, axis=1).reshape(-1, 1)) / np.std(avg_data, axis=1).reshape(-1, 1)
avg_zscore_df = pd.DataFrame(data=np.hstack((gene_df['gene'].to_numpy().reshape(-1, 1), avg_zscore)))

# %% abundance of gene for cluster and sort the heat map
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering

condition_set = list(set(meta['C'].tolist()))
condition_set.sort()
R_filter = gene_df['gene'].isin(hui_GO['A']['name'])
R_hist_num = np.sum(R_filter)
R_sector = htTPM_data_all[R_filter]

R_data = np.empty((R_hist_num, len(condition_set)))
for s_i, sample in enumerate(condition_set):
    sample_mask = meta['C'] == sample
    sample_data = R_sector[:, sample_mask]
    sample_avg = np.mean(sample_data, axis=1)
    R_data[:, s_i] = sample_avg

R_total = np.sum(R_data, axis=0)
R_total = pd.DataFrame(data=np.array([condition_set, R_total, R_total / 1e6]).T,
                       columns=['condition', 'TPM_total', 'psi'])
R_total = R_total.T

R_zscore = (R_data - np.average(R_data, axis=1).reshape(-1, 1)) / np.std(R_data, axis=1).reshape(-1, 1)
R_zscore_df = pd.DataFrame(data=np.hstack((gene_df[R_filter]['gene'].to_numpy().reshape(-1, 1), R_zscore)))

model = AgglomerativeClustering(n_clusters=8, compute_distances=True)
model = model.fit(R_zscore)
sort_index = np.argsort(model.labels_)
sorted_gens = R_zscore_df[0][sort_index]
sorted_zscores = R_zscore[sort_index, :]
sorted_tpm = R_data[sort_index]
sorted_labels = model.labels_[sort_index]

fig1, ax1 = plt.subplots(figsize=(10, 10))
ax1.imshow(sorted_zscores, interpolation='nearest', cmap='jet', aspect='auto')
fig1.show()

# show subgroup
label = 7
mask = sorted_labels == label
sub_sorted_genes = sorted_gens[mask]
sub_sorted_zscores = sorted_zscores[mask, :]
sub_sorted_tpm = sorted_tpm[mask, :]
fig2, ax2 = plt.subplots(figsize=(10, 10))
ax2.imshow(sub_sorted_zscores, interpolation='nearest', cmap='jet', aspect='auto')
fig2.show()
# %% get subset

# gene_set = ['aspC', 'gltB', 'gltD', 'glnA']
# gene_set = ['gadA', 'gadB'] # glutamate util
# gene_set = ['alaA'] # alanine syn
# aa synthesis
gene_set = ['alaA', 'argA', 'argB', 'argC', 'argD', 'argE', 'argF', 'argH', 'argG', 'argI', 'aspC', 'asnA', 'asnB',
            'cysK', 'cysE', 'aspC', 'gltB', 'gltD', 'glnA', 'glyA', 'ltaE', 'hisB', 'hisC', 'hisD', 'hisH', 'hisF',
            'hisA', 'hisI', 'hisG', 'ilvI', 'ilvH', 'ilvN', 'ilvB', 'ilvM', 'ilvC', 'ilvD', 'ilvE', 'leuB', 'leuC',
            'leuD', 'tdcB', 'ilvA', 'ilvI', 'ilvH', 'ilvN', 'ilvB', 'ilvM', 'ilvC', 'ilvD', 'ilvE', 'leuA', 'leuB',
            'leuC', 'leuD', 'tdcB', 'ilvA', 'metL', 'thrA', 'lysC', 'asd', 'dapA', 'dapB', 'dapD', 'dapE', 'dapF',
            'lysA', 'argD', 'metE', 'metH', 'metA', 'metB', 'metC', 'metL', 'malY', 'thrA', 'proA', 'proB', 'proC',
            'serA', 'serB', 'serC', 'glyA', 'metL', 'ltaE', 'thrA', 'thrB', 'thrC', 'aroC', 'aroA', 'aroL', 'aroK',
            'aroE', 'aroD', 'aroB', 'aroG', 'aroH', 'aroF', 'trpA', 'trpB', 'trpC', 'trpD', 'trpE', 'pheA', 'tyrB',
            'tyrA', 'aroC', 'aroA', 'aroL', 'aroK', 'aroE', 'aroD', 'aroB', 'aroG', 'aroH', 'aroF', 'ilvI', 'ilvH',
            'ilvN', 'ilvB', 'ilvM', 'ilvC', 'ilvD', 'ilvE', 'tdcB', 'ilvA', 'pheA', 'tyrB', 'tyrA', 'aroC', 'aroA',
            'aroL', 'aroK', 'aroE', 'aroD', 'aroB', 'aroG', 'aroH', 'aroF']
# gene_set = ['argA', 'argB', 'argC', 'argD', 'argE', 'argF', 'argH', 'argG', 'argI']
# gene_set = ['satP', 'caiT', 'thiQ', 'thiP', 'thiB', 'setA', 'aroP', 'fhuA', 'fhuC', 'fhuD', 'fhuB', 'btuF', 'metQ', 'metI', 'metN', 'mmuP', 'betT', 'codB', 'mhpT', 'tauA', 'tauB', 'tauC', 'sbmA', 'brnQ', 'ampG', 'amtB', 'acrB', 'acrA', 'emrE', 'pheP', 'fepA', 'fepC', 'fepG', 'fepD', 'entS', 'fepB', 'cstA', 'citT', 'dcuC', 'gltL', 'gltK', 'gltJ', 'gltI', 'potE', 'dtpD', 'pnuC', 'zitB', 'acrZ', 'modA', 'modB', 'modC', 'fiu', 'glnQ', 'glnP', 'glnH', 'rhtA', 'gsiA', 'gsiB', 'gsiC', 'gsiD', 'mdfA', 'potF', 'potG', 'potH', 'potI', 'artJ', 'artM', 'artQ', 'artP', 'lysO', 'cydC', 'cydD', 'ssuB', 'ssuC', 'ssuA', 'rutG', 'putA', 'putP', 'mdtH', 'fhuE', 'potD', 'potC', 'potB', 'potA', 'dauA', 'narK', 'oppA', 'oppB', 'oppC', 'oppD', 'oppF', 'sapF', 'sapD', 'sapC', 'sapB', 'puuP', 'mppA', 'abgT', 'zntB', 'narU', 'yddG', 'gadC', 'lsrA', 'lsrC', 'lsrD', 'lsrB', 'eamA', 'mdtI', 'mdtJ', 'tqsA', 'uidB', 'dtpA', 'mdtK', 'btuD', 'btuC', 'tcyP', 'nimT', 'leuE', 'mntP', 'znuA', 'znuC', 'znuB', 'tyrP', 'tcyN', 'tcyL', 'tcyJ', 'shiA', 'yeeO', 'plaP', 'mdtA', 'mdtB', 'mdtC', 'mdtD', 'rcnA', 'yehW', 'yehX', 'yehY', 'osmF', 'mglC', 'mglA', 'mglB', 'cirA', 'lysP', 'setB', 'bcr', 'yojI', 'glpT', 'hisP', 'hisM', 'hisQ', 'hisJ', 'argT', 'dsdX', 'mntH', 'nupC', 'xapB', 'cysZ', 'cysA', 'cysW', 'cysU', 'cysP', 'murP', 'acrD', 'eamB', 'kgtP', 'gabP', 'alaE', 'proV', 'proW', 'proX', 'ygaZ', 'ygaH', 'emrA', 'emrB', 'sdaC', 'fucP', 'lplT', 'xanQ', 'ghxQ', 'uacT', 'argO', 'galP', 'nupG', 'glcA', 'pitB', 'zupT', 'ttdT', 'sstT', 'exuT', 'tdcC', 'garP', 'mtr', 'nanT', 'panF', 'acrE', 'acrF', 'nirC', 'feoB', 'gntT', 'gntU', 'ugpC', 'ugpE', 'ugpA', 'ugpB', 'livF', 'livG', 'livM', 'livH', 'livK', 'livJ', 'yhhQ', 'nikA', 'nikB', 'nikC', 'nikD', 'nikE', 'pitA', 'dtpB', 'arsB', 'mdtE', 'mdtF', 'dctA', 'dppF', 'dppD', 'dppC', 'dppB', 'dppA', 'xylF', 'xylG', 'xylH', 'lldP', 'gltS', 'xanP', 'nepI', 'adeQ', 'uhpT', 'emrD', 'tnaB', 'adeP', 'pstB', 'pstA', 'pstC', 'pstS', 'rbsA', 'rbsC', 'rbsB', 'corA', 'rhtC', 'rhtB', 'bioP', 'yihO', 'rhaT', 'kdgT', 'fieF', 'sbp', 'btuB', 'xylE', 'malG', 'malF', 'malE', 'malK', 'ghxP', 'actP', 'gltP', 'mdtP', 'mdtO', 'mdtN', 'alsC', 'alsA', 'alsB', 'phnE', 'phnD', 'phnC', 'proP', 'adiC', 'melB', 'dcuB', 'dtpC', 'cadB', 'dcuA', 'yjeH', 'gdx', 'cycA', 'ytfQ', 'idnT', 'fecE', 'fecD', 'fecC', 'fecB', 'fecA', 'gntP', 'mdtM', 'btsT', 'lgoT']
condition_set = list(set(meta['C'].tolist()))
condition_set.sort()
R_filter = gene_df['gene'].isin(gene_set)
R_hist_num = np.sum(R_filter)
R_sector = htTPM_data_all[R_filter]

R_data = np.empty((R_hist_num, len(condition_set)))
for s_i, sample in enumerate(condition_set):
    sample_mask = meta['C'] == sample
    sample_data = R_sector[:, sample_mask]
    sample_avg = np.mean(sample_data, axis=1)
    R_data[:, s_i] = sample_avg

R_total = np.sum(R_data, axis=0)
R_total = pd.DataFrame(data=np.array([condition_set, R_total, R_total / 1e6]).T,
                       columns=['condition', 'TPM_total', 'psi'])
R_total = R_total.T
R_zscore = (R_data - np.average(R_data, axis=1).reshape(-1, 1)) / np.std(R_data, axis=1).reshape(-1, 1)
R_zscore_df = pd.DataFrame(data=np.hstack((gene_df[R_filter]['gene'].to_numpy().reshape(-1, 1), R_zscore)))
R_tpm_df = pd.DataFrame(data=np.hstack((gene_df[R_filter]['gene'].to_numpy().reshape(-1, 1), R_data)))

fig1, ax1 = plt.subplots(figsize=(10, 10))
ax1.imshow(R_zscore, interpolation='nearest', cmap='coolwarm', aspect='auto')
fig1.show()

model = AgglomerativeClustering(n_clusters=8, compute_distances=True)
model = model.fit(R_zscore)
sort_index = np.argsort(model.labels_)
sorted_gens = R_zscore_df[0][sort_index]
sorted_zscores = R_zscore[sort_index, :]
sorted_tpm = R_data[sort_index]
sorted_labels = model.labels_[sort_index]

fig1, ax1 = plt.subplots(figsize=(10, 10))
ax1.imshow(sorted_zscores, interpolation='nearest', cmap='jet', aspect='auto')
fig1.show()

# show subgroup
label = 2
mask = sorted_labels == label
sub_sorted_genes = sorted_gens[mask]
sub_sorted_zscores = sorted_zscores[mask, :]
sub_sorted_tpm = sorted_tpm[mask, :]
fig2, ax2 = plt.subplots(figsize=(10, 10))
ax2.imshow(sub_sorted_zscores, interpolation='nearest', cmap='jet', aspect='auto')
fig2.show()

# %%  statistic genes in toggle
condition_set = list(set(meta['C'].tolist()))
condition_set.sort()

toggle_gene_list = ['LacI', 'tetR', 'mCherry', 'aph%283%27%29-II %28or nptII%29', 'GFP']
R_filter = gene_df['gene'].isin(toggle_gene_list)
R_hist_num = np.sum(R_filter)
R_sector = htTPM_data_all[R_filter]

R_data = np.empty((R_hist_num, len(condition_set)))
for s_i, sample in enumerate(condition_set):
    sample_mask = meta['C'] == sample
    sample_data = R_sector[:, sample_mask]
    sample_avg = np.mean(sample_data, axis=1)
    R_data[:, s_i] = sample_avg

R_total = np.sum(R_data, axis=0)
R_total = pd.DataFrame(data=np.array([condition_set, R_total, R_total / 1e6]).T,
                       columns=['condition', 'TPM_total', 'psi'])

# %%
# run GOterm enrichment in web server panther
enrich_pars = """
{
  "organism": "83333",
  "refOrganism": "83333",
  "annotDataSet": "GO:0008150",
  "enrichmentTestType": "FISHER",
  "correction": "FDR"
}"""
request = get_request_obj('enrich')
request_parameters = json.loads(enrich_pars)
request_parameters['geneInputList'] = ','.join(sub_sorted_genes.to_list())
request.parameters = request_parameters
response = request.call_service()
response.print_results()

# %% query gene list
info_pars = """
{
  "organism": "83333"
}
"""
request = get_request_obj('geneinfo')
request_parameters = json.loads(info_pars)
request_parameters['geneInputList'] = ','.join(['purH'])
request.parameters = request_parameters
response = request.call_service()
response.print_results()
