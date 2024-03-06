# -*- coding: utf-8 -*-

"""
 This script compares two sets of RNA-seq, miniBac data and the steady states of E.coli strains in

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
import matplotlib.pyplot as plt
import json

# %%
# load data, and clean data (remove genes without transcription)
data_1 = r'D:\OneDrive\zjw_data\colony\data\omic_data\202210_20_RNA-seq.xlsx'
data_2 = r'D:\OneDrive\zjw_data\colony\data\omic_data\20240112_miniBac-seq.xlsx'

tpm_1_df = pd.read_excel(data_1, sheet_name='TPM_mean_all_conditions')
tpm_2_df = pd.read_excel(data_2, sheet_name='conditions_TPM_average')
tpm_1_info = pd.read_excel(data_1, sheet_name='sample_info')
tpm_2_info = pd.read_excel(data_2, sheet_name='sample_info')

all_tpm = np.hstack([tpm_1_df.to_numpy()[:, 1:].astype(np.float64),
                     tpm_2_df.to_numpy()[:, 1:].astype(np.float64)])
all_zscore = (all_tpm - np.mean(all_tpm, axis=1).reshape(-1, 1)) / np.std(all_tpm, axis=1).reshape(-1, 1)
nan_mask = np.any(~np.isnan(all_zscore), axis=1)

tpm_1_cds_name = tpm_1_df.columns[1:].to_numpy()
tpm_1_media_name = np.array([tpm_1_info[tpm_1_info['Condition'] == i]['media'].tolist()[0] for i in tpm_1_cds_name])
tpm_2_cds_name = tpm_2_df.columns[1:].to_numpy()
all_condition = np.hstack((tpm_1_media_name, tpm_2_cds_name))

all_genes_list = tpm_1_df['gene'].to_numpy()
cleaned_z_score = all_zscore[nan_mask, :]
cleaned_gene_list = all_genes_list[nan_mask]

# %% decomposition via PCA
from sklearn.decomposition import PCA

model = PCA(svd_solver='full', n_components=3)
model.fit(cleaned_z_score[5:, :].T)

tpm_1_pca = model.transform(cleaned_z_score[5:, :].T)
# tpm_1_pca = model.transform(cleaned_z_score[5:, :len(tpm_1_df.columns)-1].T)
# minibac_pca = model.transform(cleaned_z_score[5:, len(tpm_1_df.columns)-1:].T)

fig1, ax1 = plt.subplots(1, 1)
ax1.scatter(tpm_1_pca[:len(tpm_1_cds_name), 0], tpm_1_pca[:len(tpm_1_cds_name), 1])
ax1.scatter(tpm_1_pca[len(tpm_1_cds_name):, 0], tpm_1_pca[len(tpm_1_cds_name):, 1],
            label=tpm_2_df.columns.to_list()[1:])
# plt.legend()
fig1.show()

# %% decomposition via PCA
from sklearn.decomposition import PCA

mask = np.logical_and(tpm_1_cds_name != 7, tpm_1_cds_name != 3)
model = PCA(svd_solver='full', n_components=3)
model.fit(cleaned_z_score[5:, :len(tpm_1_cds_name)][:, mask].T)

tpm_1_pca = model.transform(cleaned_z_score[5:, :].T)
# tpm_1_pca = model.transform(cleaned_z_score[5:, :len(tpm_1_df.columns)-1].T)
# minibac_pca = model.transform(cleaned_z_score[5:, len(tpm_1_df.columns)-1:].T)

fig1, ax1 = plt.subplots(1, 1)
ax1.scatter(tpm_1_pca[:len(tpm_1_cds_name), 0], tpm_1_pca[:len(tpm_1_cds_name), 1])
ax1.scatter(tpm_1_pca[len(tpm_1_cds_name):, 0], tpm_1_pca[len(tpm_1_cds_name):, 1],
            label=tpm_2_df.columns.to_list()[1:])
# plt.legend()
fig1.show()
# %% decomposition via t-SNE
from sklearn.manifold import TSNE

model = TSNE(n_components=2)
# model.fit(cleaned_z_score[5:, :].T)
tpm_1_tSNE = model.fit_transform(cleaned_z_score[5:, :].T)
# minibac_tSNE = model.fit_transform(cleaned_z_score[5:, len(tpm_1_df.columns)-1:].T)

fig1, ax1 = plt.subplots(1, 1)
ax1.scatter(tpm_1_tSNE[:len(tpm_1_df.columns) - 1, 0], tpm_1_tSNE[:len(tpm_1_df.columns) - 1, 1])
ax1.scatter(tpm_1_tSNE[len(tpm_1_df.columns) - 1:, 0], tpm_1_tSNE[len(tpm_1_df.columns) - 1:, 1],
            label=tpm_2_df.columns.to_list()[1:])
# plt.legend()
fig1.show()

# %%
with open(r'./exported_data/GOterm_hui_MSB_2015.json') as jfile:
    hui_GO = json.load(jfile)

# %%
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering

# condition_set = samples_list
# condition_set.sort()
condition_set = all_condition

# genes_list = hui_GO['R']['name']
# ribiosome, Wu 2023
wu_classification = dict(
    # Translational proteins
    ribosomal_proteins=['rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplI', 'rplJ', 'rplK', 'rplL', 'rplM', 'rplN',
                        'rplO', 'rplP', 'rplQ', 'rplR', 'rplS', 'rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY', 'rpmA',
                        'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG' 'rpmH', 'rpmI', 'rpmJ', 'rpsA', 'rpsB', 'rpsC',
                        'rpsD', 'rpsE', 'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK', 'rpsL', 'rpsM', 'rpsN', 'rpsO',
                        'rpsP', 'rpsQ', 'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra'],
    affiliated_translational_apparatus=['arfA', 'arfB', 'efp', 'frr', 'fusA', 'infA', 'infB', 'infC', 'lepA', 'prfA',
                                        'prfB', 'prfC', 'tsf', 'tufA', 'tufB'],
    tRNA_synthesis=['alaS', 'argS', 'asnS', 'aspS', 'cysS', 'glnS', 'gltX', 'glyQ', 'glyS', 'hisS', 'ileS', 'leuS',
                    'lysS', 'lysU', 'metG', 'pheS', 'pheT', 'proS', 'serS', 'thrS', 'trpS', 'tyrS', 'valS'],
    # Transporters
    aa_transport=['abgT', 'alaE', 'ansP', 'argO', 'argT', 'aroP', 'artI', 'artJ', 'artM', 'artP', 'artQ', 'brnQ',
                  'cadB',
                  'cstA', 'cycA', 'dppA', 'dppB', 'dppC', 'dppD', 'dppF', 'dtpA', 'dtpB', 'dtpC', 'dtpD', 'eamA',
                  'eamB',
                  'frlA', 'gabP', 'gadC', 'glnH', 'glnP', 'glnQ', 'gltI', 'gltJ', 'gltK', 'gltL', 'gltP', 'gltS',
                  'hisJ',
                  'leuE', 'livF', 'livG', 'livH', 'livJ', 'livM', 'lysP', 'metQ', 'mmuP', 'mtr', 'oppA', 'oppB', 'oppC',
                  'oppD', 'oppF', 'pheP', 'plaP', 'potA', 'potB', 'potC', 'potD', 'potE', 'potF', 'potG', 'potH',
                  'potI',
                  'proV', 'proW', 'proX', 'proY', 'putP', 'puuP', 'rhtA', 'rhtB', 'rhtC', 'sgrR', 'sstT', 'tcyJ',
                  'tcyL',
                  'tcyN', 'tcyP', 'tdcC', 'tnaB', 'tyrP'],
    carbon_catabolism=['aceA', 'aceB', 'acs', 'alsA', 'alsB', 'alsC', 'araF', 'araG', 'araH', 'crr', 'fruA', 'fruB',
                       'fucK', 'fucP', 'galE', 'galF', 'galK', 'galM', 'galP', 'galR', 'galS', 'galT', 'galU', 'gatA',
                       'gatB', 'gatC', 'lacY', 'lacZ', 'malE', 'malF', 'malG', 'malK', 'malM', 'malP', 'malT', 'malX',
                       'manX', 'manY', 'manZ', 'melB', 'mglA', 'mglB', 'mglC', 'mtlA', 'mtlD', 'nagE', 'ptsG', 'ptsH',
                       'ptsI', 'rbsA', 'rbsB', 'rbsC', 'rbsD', 'rbsK', 'srlA', 'srlB', 'srlE', 'treB', 'ugpA', 'ugpB',
                       'ugpC', 'ugpE', 'ulaA', 'ulaB', 'ulaC', 'xylA', 'xylB', 'xylE', 'xylF', 'xylG', 'xylH'],
    glycerol_uptake=['glpA', 'glpB', 'glpC', 'glpD', 'glpF', 'glpK'],
    outer_membrane_porin=['nmpC', 'ompA', 'ompC', 'ompF', 'ompG', 'ompN', 'phoE'],
    # Nucleotide biosynthesis
    nucleotide_biosynthesis=['adk', 'cmk', 'deoB', 'dut', 'gmk', 'guaA', 'guaB', 'ndk', 'nrdA', 'nrdB', 'nrdE', 'nrdF',
                             'prs', 'purA', 'purB', 'purC', 'purD', 'purE', 'purF', 'purH', 'purK', 'purL', 'purM',
                             'purN', 'purT', 'pyrB', 'pyrC', 'pyrD', 'pyrE', 'pyrF', 'pyrG', 'pyrH', 'pyrI', 'pyrL',
                             'thyA', 'tmk', 'trxA'],
    # Central carbon metabolism and energy
    TCA=['acnA', 'acnB', 'fumA', 'fumB', 'fumC', 'fumD', 'fumE', 'gltA', 'icd', 'lpd*', 'sdhA', 'sdhB', 'sdhC', 'sdhD',
         'sucA', 'sucB', 'sucC', 'sucD'],
    glycolysis_gluco_neo_genesis=['eno', 'fbaA', 'fbaB', 'fbp', 'gapA', 'glpX', 'gpmA', 'gpmM', 'maeA', 'maeB', 'mdh',
                                  'mqo', 'pck', 'pfkA', 'pfkB', 'pgi', 'pgk', 'ppsA', 'pykA', 'pykF', 'tpiA', 'ybhA'],
    fermentation_to_acetate=['aceE', 'aceF', 'ackA', 'lpd', 'pta'],
    ATPase=['atpA', 'atpB', 'atpC', 'atpD', 'atpE', 'atpF', 'atpG', 'atpH', 'atpI'],
    # AA biosynthesis
    arg_group=['argA', 'argB', 'argC', 'argD', 'argE', 'argF', 'argG', 'argH', 'argI'],
    aro_group=['aroA', 'aroB', 'aroC', 'aroD', 'aroE', 'aroF', 'aroG', 'aroH', 'aroK', 'aroL'],
    cyc_group=['cysC', 'cysD', 'cysE', 'cysH', 'cysI', 'cysJ', 'cysK', 'cysM', 'cysN'],
    glt_group=['gdhA', 'gltB', 'gltD'],
    his_group=['hisA', 'hisB', 'hisC', 'hisD', 'hisF', 'hisG', 'hisH', 'hisI'],
    ilv_group=['ilvA', 'ilvB', 'ilvC', 'ilvD', 'ilvE', 'ilvG_1', 'ilvH', 'ilvI', 'ilvM', 'ilvN'],
    leu_group=['leuA', 'leuB', 'leuC', 'leuD'],
    lys_group=['asd', 'dapA', 'dapB', 'dapD', 'dapE', 'dapF', 'lysA', 'lysC'],
    met_group=['metA', 'metB', 'metC', 'metE', 'metH', 'metL'],
    other_group=['alaA', 'alaC', 'asnA', 'asnB', 'aspC', 'glnA', 'glyA', 'proA', 'proB', 'proC'],
    phetyr_group=['pheA', 'tyrA', 'tyrB'],
    ser_group=['serA', 'serB', 'serC'],
    thr_group=['thrA', 'thrB', 'thrC'],
    trp_group=['trpA', 'trpB', 'trpC', 'trpD', 'trpE'])
wu_classification_structure = dict(
    translational_proteins=['ribosomal_proteins', 'affiliated_translational_apparatus', 'tRNA_synthesis'],
    transporters=['aa_transport', 'carbon_catabolism', 'glycerol_uptake', 'outer_membrane_porin'],
    nucleotide_biosynthesis=['nucleotide_biosynthesis'],
    central_carbon_metabolism_and_energy=['TCA', 'glycolysis_gluco_neo_genesis', 'fermentation_to_acetate',
                                          'ATPase'],
    aa_biosynthesis=['arg_group', 'aro_group', 'cyc_group', 'cyc_group', 'glt_group', 'his_group', 'ilv_group',
                     'leu_group', 'lys_group', 'met_group', 'other_group', 'phetyr_group', 'ser_group', 'thr_group',
                     'trp_group'])
for key, genes_list in wu_classification.items():
    genes_list = genes_list
    gene_filter = tpm_1_df['gene'].isin(genes_list)
    gene_hist_num = np.sum(gene_filter)
    genes_sector = all_genes_list[gene_filter]
    tpm_sector = all_tpm[gene_filter, :]
    zscore_sector = all_zscore[gene_filter, :]

    genes_total = np.sum(tpm_sector, axis=0)
    genes_total = pd.DataFrame(data=np.array([condition_set, genes_total, genes_total / 1e6]).T,
                               columns=['condition', 'TPM_total', 'psi'])
    genes_total_T = genes_total.T

    fig1, ax1 = plt.subplots(figsize=(20, 10))
    ax1.bar(genes_total['condition'], genes_total['psi'])
    ax1.set_title(key)
    fig1.show()
#%%
model = AgglomerativeClustering(n_clusters=6, compute_distances=True)
model = model.fit(zscore_sector)
sort_index = np.argsort(model.labels_)
sorted_gens = genes_sector[sort_index]
sorted_zscores = zscore_sector[sort_index, :]
sorted_tpm = tpm_sector[sort_index]
sorted_labels = model.labels_[sort_index]

fig1, ax1 = plt.subplots(figsize=(10, 10))
ax1.imshow(sorted_zscores, interpolation='nearest', cmap='jet', aspect='auto')
fig1.show()

# show subgroup
# label = 5
# mask = sorted_labels == label
# sub_sorted_genes = sorted_gens[mask]
# sub_sorted_zscores = sorted_zscores[mask, :]
# sub_sorted_tpm = sorted_tpm[mask, :]
# fig2, ax2 = plt.subplots(figsize=(10, 10))
# ax2.imshow(sub_sorted_zscores, interpolation='nearest', cmap='jet', aspect='auto')
# fig2.show()

# %% 2.a analysis by each sample not each condition

data_1 = r'D:\OneDrive\zjw_data\colony\data\omic_data\202210_20_RNA-seq.xlsx'
data_2 = r'D:\OneDrive\zjw_data\colony\data\omic_data\20240112_miniBac-seq.xlsx'

tpm_1_df = pd.read_excel(data_1, sheet_name='TPM_all_samples')
tpm_2_df = pd.read_excel(data_2, sheet_name='samples_RNA-seq_TPM')
tpm_1_info = pd.read_excel(data_1, sheet_name='sample_info')
tpm_2_info = pd.read_excel(data_2, sheet_name='sample_info')

tpm_1_data = tpm_1_df.to_numpy()[:, 3:].astype(np.float64)
tpm_2_data = tpm_2_df.to_numpy()[:, 2:].astype(np.float64)
print(tpm_1_data.shape, tpm_2_data.shape)

all_tpm = np.hstack((tpm_1_data, tpm_2_data))

all_zscores = (all_tpm - np.mean(all_tpm, axis=1).reshape(-1, 1)) / np.std(all_tpm, axis=1).reshape(-1, 1)
nan_mask = np.any(~np.isnan(all_zscores), axis=1)

tpm_1_sample_name = tpm_1_df.columns[3:].to_numpy()
tpm_2_sample_name = tpm_2_df.columns[2:].to_numpy()
tpm_1_condition = np.array(
    [tpm_1_info[tpm_1_info['sample_name'] == i]['media'].to_list()[0] for i in tpm_1_sample_name])
tpm_2_condition = np.array([tpm_2_info[tpm_2_info['sample'] == i]['condition'].to_list()[0] for i in tpm_2_sample_name])
all_condition = np.hstack((tpm_1_condition, tpm_2_condition))

genes_list = tpm_1_df['gene'].to_numpy()
cleaned_z_score = all_zscores[nan_mask, :]
cleaned_genes_list = genes_list[nan_mask]
# %% 2.b analysis by each sample not each condition decomposition via PCA
from sklearn.decomposition import PCA

model = PCA(svd_solver='full', n_components=3)
model.fit(cleaned_z_score[5:, :len(tpm_1_sample_name)].T)  # remove toggle genes

eigen = model.components_
zscore_pca = model.transform(cleaned_z_score[5:, :].T)
# tpm_1_pca = model.transform(cleaned_z_score[5:, :len(tpm_1_df.columns)-1].T)
# minibac_pca = model.transform(cleaned_z_score[5:, len(tpm_1_df.columns)-1:].T)

fig1, ax1 = plt.subplots(1, 1)
ax1.scatter(zscore_pca[:len(tpm_1_sample_name), 0], zscore_pca[:len(tpm_1_sample_name), 1])
ax1.scatter(zscore_pca[len(tpm_1_sample_name):, 0], zscore_pca[len(tpm_1_sample_name):, 1])
# plt.legend()
fig1.show()

cp_1_gene_order = np.take(genes_list[5:], np.argsort(eigen[0, :])[::-1])
cp_2_gene_order = np.take(genes_list[5:], np.argsort(eigen[1, :])[::-1])
cp_3_gene_order = np.take(genes_list[5:], np.argsort(eigen[2, :])[::-1])

# %% 2.c PCA without 'MOPS Glutamate Glucose'
mask = tpm_1_condition != 'MOPS Glutamate Glucose'
model = PCA(svd_solver='full', n_components=3)
model.fit(cleaned_z_score[5:, :len(tpm_1_sample_name)][:, mask].T)  # remove toggle genes

eigen = model.components_
zscore_pca = model.transform(cleaned_z_score[5:, :].T)
# tpm_1_pca = model.transform(cleaned_z_score[5:, :len(tpm_1_df.columns)-1].T)
# minibac_pca = model.transform(cleaned_z_score[5:, len(tpm_1_df.columns)-1:].T)

fig1, ax1 = plt.subplots(1, 1)
ax1.scatter(zscore_pca[:len(tpm_1_sample_name), 0], zscore_pca[:len(tpm_1_sample_name), 1])
ax1.scatter(zscore_pca[len(tpm_1_sample_name):, 0], zscore_pca[len(tpm_1_sample_name):, 1])
# plt.legend()
fig1.show()
