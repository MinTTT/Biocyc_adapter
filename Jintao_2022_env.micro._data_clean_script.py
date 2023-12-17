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


# for g_i in no_hit_index_r2:
#     g_name = terry_gene_list[g_i]
#     g_locus = terry_locus_list[g_i]
#     print(g_name, g_locus)


for index, hit in hit_fet.items():
    db = hit.qualifiers['db_xref']
    for id in db:
        if id.split(':')[0] == 'ECOCYC':
            eco_number[index] = id.split(':')[-1]


new_data_col = new_data.columns.to_list()
new_data = new_data.to_numpy()
new_data_extended = np.hstack([np.array(eco_number).reshape(-1, 1), new_data])
new_data_extended_pd = pd.DataFrame(data=new_data_extended, columns=['ECOCYC'] + new_data_col)
new_data_extended_pd.to_excel(r'./transcriptiome_data/new_sci_data.xlsx')


