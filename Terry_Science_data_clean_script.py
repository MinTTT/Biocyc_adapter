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
# […]
def get_simple_name(name):
    if '_' in name:
        name = name.split('_')[0]
    return name

# Own modules

terry_sci_data = pd.read_excel(r'./transcriptiome_data/science.abk2066_table_s3.xlsx', '2 - RNAseq-ss (fractions)')
terry_gene_list = terry_sci_data['gene'].to_list()
terry_locus_list = terry_sci_data['locus'].to_list()
terry_gene_set = set(terry_sci_data['gene'].to_list())
terry_locus_set = set(terry_sci_data['locus'].to_list())



ecocyc_gene_info = pd.read_csv(r'./exported_data/Ecoli_gene_info.csv')
all_gene_info = pd.read_csv(r'./exported_data/Ecoli_all_genes.csv')
all_features = pd.read_csv(r'./exported_data\Ecoli_features.NC_000913.3.csv')
gb_table = SeqIO.read(r'.\exported_data\NC_000913.3.gb', 'genbank')


eco_number = [None] * len(terry_gene_list)
no_hit_index = []
for g_i, g_name in enumerate(terry_gene_list):
    g_locus = terry_locus_list
    g_n = get_simple_name(g_name)  # some gene have more than one copy
    g_id = all_features[all_features['genes_name'] == g_n]['genes_id'].tolist()
    if g_id:
        eco_number[g_i] = g_id[0]
    else:
        no_hit_index.append(g_i)
        print(f'{g_name}: no record.')

hit_fet = {}
no_hit_index_r2 = []
for g_i in no_hit_index:
    g_name = terry_gene_list[g_i]
    g_n = get_simple_name(g_name)
    g_locus = terry_locus_list[g_i]
    hit = False
    for fet in gb_table.features:
        for key, value in fet.qualifiers.items():
            if g_n in value:
                hit_fet[g_i] = fet
                hit = True
                print(g_n, value)
                break
            elif g_locus in value:
                # if key == 'locus_tag':
                hit_fet[g_i] = fet
                hit = True
                print('find gene record via locus tag:', g_locus, value, key)
                break
        if hit:
            break
    if not hit:
        no_hit_index_r2.append(g_i)
        # pass

# find gene in ecocyc record:
no_hit_index_r3 = []
for g_i in no_hit_index_r2:
    g_name = terry_gene_list[g_i]
    g_n = get_simple_name(g_name)
    g_locus = terry_locus_list[g_i]
    hit = False
    for recd_i, rcd in ecocyc_gene_info.iterrows():
        if isinstance(rcd['synonym'], str):
            syn = rcd['synonym'].split('; ')
            if g_locus in syn:
                eco_number[g_i] = rcd['id'].split(':')[-1]
                hit = True
                break
            elif g_n in syn:
                eco_number[g_i] = rcd['id'].split(':')[-1]
                hit = True
                break
    if hit is False:
        print(g_name, g_n, g_locus)
        no_hit_index_r3.append(g_i)




for index, hit in hit_fet.items():
    db = hit.qualifiers['db_xref']
    for id in db:
        if id.split(':')[0] == 'ECOCYC':
            eco_number[index] = id.split(':')[-1]


terry_data_col = terry_sci_data.columns.to_list()
terry_data = terry_sci_data.to_numpy()
terry_data_extended = np.hstack([np.array(eco_number).reshape(-1, 1), terry_data])
terry_data_extended_pd = pd.DataFrame(data=terry_data_extended, columns=['ECOCYC'] + terry_data_col)
terry_data_extended_pd.to_excel(r'./transcriptiome_data/terry_sci_data.xlsx')



# ecocyc_gene_set = set(ecocyc_gene_info['name'].to_list())
# ecocyc_locus_set = set(ecocyc_gene_info['locus_tag'].to_list())
# terry_only_set = terry_locus_set - ecocyc_locus_set
# ecocyc_only_set = ecocyc_locus_set - terry_locus_set
# terry_only_name = terry_sci_data[[True if locus in terry_only_set else False for locus in terry_sci_data['locus']]][
#     'gene']
#
# ecocyc_only_table = all_features[[True if locus in ecocyc_only_set else False for locus in all_features['locus_tag']]]
#
# terry_only_table = all_features[[True if locus in terry_only_set else False
#                                  for locus in all_features['locus_tag']]]
# terry_only_table_by_name = all_features[[True if name in terry_only_name else False
#                                          for name in all_features['genes_name']]]
