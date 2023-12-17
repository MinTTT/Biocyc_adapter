# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
import sys

import numpy as np
# […]

# Libs
import pandas as pd
# import numpy as np  # Or any other
from Bio import SeqIO



ecocyc_gene_info = pd.read_csv(r'./exported_data/Ecoli_gene_info.csv')
all_gene_info = pd.read_csv(r'./exported_data/Ecoli_all_genes.csv')
all_features = pd.read_csv(r'./exported_data/Ecoli_features.NC_000913.3.csv')
gb_table = SeqIO.read(r'./exported_data/NC_000913.3.gb', 'genbank')

#%%

new_data = pd.read_csv(r'./transcriptiome_data/GSE195463_Raw_gene_counts_matrix.txt',
                       delimiter='\t', encoding='utf-16')


locus_list = new_data['Gene'].tolist()

eco_number = [None] * len(locus_list)
no_hit_index = []
for g_i, g_locus in enumerate(locus_list):

    g_id = all_features[all_features['locus_tag'] == g_locus]['genes_id'].tolist()
    if g_id:
        eco_number[g_i] = g_id[0]
    else:
        no_hit_index.append(g_i)
        print(g_locus)

hit_fet = {}
no_hit_index_r2 = []
for g_i in no_hit_index:
    # g_name = terry_gene_list[g_i]
    g_locus = locus_list[g_i]
    hit = False
    for fet in gb_table.features:
        for key, value in fet.qualifiers.items():
            # if g_name in value:
            #     hit_fet[g_i] = fet
            #     hit = True
            #     print(g_name, value)
            #     break
            if g_locus in value:
                # if key == 'locus_tag':
                hit_fet[g_i] = fet
                hit = True
                print(g_locus, value, key)
                break
        if hit:
            break
    if not hit:
        no_hit_index_r2.append(g_i)
        # pass

no_hit_index_r3 = []
for g_i in no_hit_index_r2:


    g_locus = locus_list[g_i]
    hit = False
    for recd_i, rcd in ecocyc_gene_info.iterrows():
        if isinstance(rcd['synonym'], str):
            syn = rcd['synonym'].split('; ')
            if g_locus in syn:
                eco_number[g_i] = rcd['id'].split(':')[-1]
                hit = True
                break

    if hit is False:
        print(g_locus)
        no_hit_index_r3.append(g_i)

# for g_i in no_hit_index_r2:
#     g_name = terry_gene_list[g_i]
#     g_locus = terry_locus_list[g_i]
#     print(g_name, g_locus)


for index, hit in hit_fet.items():
    db = hit.qualifiers['db_xref']
    for id in db:
        if id.split(':')[0] == 'ECOCYC':
            eco_number[index] = id.split(':')[-1]


region = ecocyc_gene_info['region'].tolist()
length = [None] * len(region)
for rg_i, rg in enumerate(region):
    if 'None' not in rg:
        tg = rg[1:-1].split(', ')
        pos = [int(loc) for loc in tg]

        length[rg_i] = pos[1] - pos[0] + 1


gene_length = [None] * len(eco_number)
gene = [None] * len(eco_number)
for cyc_i, cyc in enumerate(eco_number):
    if cyc:
        query = 'ECOLI:' + cyc
        rec = ecocyc_gene_info[ecocyc_gene_info['id'] == query]
        gene_length[cyc_i] = length[rec.index.tolist()[0]]
        gene[cyc_i] = rec['name'].tolist()[0]

    # for recd_i, rcd in ecocyc_gene_info.iterrows():
    #     if rcd['id'].split(':')[-1] == cyc:
    #         gene_length[cyc_i] = length[recd_i]
    #         break



new_data_col = new_data.columns.to_list()
new_data = new_data.to_numpy()
new_data_extended = np.hstack([np.array(eco_number).reshape(-1, 1),
                               np.array(gene).reshape(-1, 1),
                               np.array(gene_length).reshape(-1, 1),
                               new_data])
new_data_extended_pd = pd.DataFrame(data=new_data_extended, columns=['ECOCYC', 'gene name', 'gene length (nt)'] + new_data_col)
new_data_extended_pd.to_excel(r'./transcriptiome_data/jintao_2022_env.micro._data.xlsx')


#%%

sample_name = ['Fig2bc_Commercial',
       'Fig2bc_Customized', 'Fig3a_replicate1_5ab_5ng',
       'Fig3a_replicate2_5ab_5ng', 'Fig3bcde_ourmethod', 'Fig3bc_TruSeq',
       'Fig4ab_1ng_replicate1', 'Fig4ab_1ng_replicate2',
       'Fig4ab_100pg_replicate1', 'Fig4ab_100pg_replicate2',
       'Fig4ab_10pg_replicate1', 'Fig4ab_10pg_replicate2',
       'Fig4ab_1pg_replicate1', 'Fig4ab_1pg_replicate2',
       'Fig4ab_100fg_replicate1', 'Fig4ab_100fg_replicate2', 'Fig5c_RTprimer1',
       'Fig5c_RTprimer2', 'Fig5c_RTprimer3', 'Fig5c_RTprimer4',
       'Fig5c_RTprimer5', 'Fig5c_RTprimer6', 'Fig5c_RTprimer7',
       'Fig5c_RTprimer8', 'Fig5c_RTprimer9', 'Fig5c_RTprimer10',
       'Fig5c_RTprimer11', 'Fig5c_RTprimer12', 'Fig5c_RTprimer13',
       'Fig5c_RTprimer14', 'Fig5c_RTprimer15', 'Fig5c_RTprimer16',
       'Fig5c_RTprimer17', 'Fig5c_RTprimer18', 'Fig5c_RTprimer19',
       'Fig5c_RTprimer20', 'Fig5c_RTprimer21', 'Fig5c_RTprimer22',
       'Fig5c_RTprimer23', 'Fig5c_RTprimer24', 'Fig6_neg_replicate1',
       'Fig6_neg_replicate2', 'Fig6_pos_replicate1', 'Fig6_pos_replicate2']

total_reads = new_data_extended_pd[sample_name].sum(axis=0)

ratio = new_data_extended_pd[sample_name] / total_reads

reads_over_length = new_data_extended_pd[sample_name].astype(np.float64).to_numpy() / new_data_extended_pd['gene length (nt)'].astype(np.float64).to_numpy().reshape(-1, 1)
weighted_reads = np.nansum(reads_over_length, axis=0)
psi = reads_over_length / weighted_reads

data_head = new_data_extended_pd[['ECOCYC', 'gene name', 'gene length (nt)', 'Gene',]].to_numpy()

psi_df = pd.DataFrame(np.hstack([data_head, psi]), columns=['ECOCYC', 'gene name', 'gene length (nt)', 'Gene',] + sample_name)
psi_df.to_excel(r'./transcriptiome_data/jintao_2022_env.micro._psi.xlsx', index=False)
