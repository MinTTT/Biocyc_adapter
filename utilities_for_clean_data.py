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
from Bio import SeqIO
import json


# Own modules
def get_simple_name(name):
    if '_' in name:
        name = name.split('_')[0]
    return name


def name2ecoid(gene_list=None, locus_tag=None):
    ecocyc_gene_info = pd.read_csv(r'./exported_data/Ecoli_gene_info.csv')
    # all_gene_info = pd.read_csv(r'./exported_data/Ecoli_all_genes.csv')
    all_features = pd.read_csv(r'.\exported_data\Ecoli_features.csv')
    gb_table = SeqIO.read(r'.\exported_data\U00096.3.gb', 'genbank')

    if locus_tag is None:
        locus_tag = [None] * len(gene_list)

    eco_number = [None] * len(gene_list)
    no_hit_index = []
    for g_i, g_name in enumerate(gene_list):
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
        g_name = gene_list[g_i]
        g_n = get_simple_name(g_name)
        g_locus = locus_tag[g_i]
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
        g_name = gene_list[g_i]
        g_n = get_simple_name(g_name)
        g_locus = locus_tag[g_i]
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

    return eco_number


# %%
if __name__ == '__main__':
    go_term = pd.read_excel('.\exported_data\GOterm_hui_MSB_2015.xlsx')
    gene_name_list = go_term['Gene'].tolist()
    eco_number = name2ecoid(gene_name_list)
    go_term_pd = pd.DataFrame(np.hstack([np.array(eco_number).reshape(-1, 1), go_term.values]),
                              columns=['id'] + list(go_term.columns))
    none_rcd_mask = go_term_pd[[True if term is None else False for term in go_term_pd['id'].tolist()]]
    cleaned_rcd_mask = go_term_pd[[False if term is None else True for term in go_term_pd['id'].tolist()]]

    go_term_list = set(cleaned_rcd_mask['Sector assigned'].tolist())
    GO_dict = {}
    for term in go_term_list:
        GO_id = cleaned_rcd_mask[cleaned_rcd_mask['Sector assigned'] == term].id.tolist()
        GO_Gene = cleaned_rcd_mask[cleaned_rcd_mask['Sector assigned'] == term]['Gene'].tolist()
        GO_dict[term] = {'name': GO_Gene, 'id': GO_id}

    GO_json = json.dumps(GO_dict, indent=4)

    with open('./exported_data/GOterm_hui_MSB_2015.json', 'w') as f:
        f.write(GO_json)

    terry_data = pd.read_excel(r'.\transcriptiome_data\terry_sci_data.xlsx')

    sample_list = ['c4',
                   'c4_1', 'c3_1', 'c3', 'c2', 'c2_1', 'c1', 'c1_1', 'c5', 'c0_1', 'a4',
                   'a3', 'a4_1', 'a4_1.1', 'a2', 'a2_1', 'a1', 'a1_1', 'r4_1', 'r3_1',
                   'r5', 'r4', 'r2_1', 'r3', 'r2', 'r1_1', 'r1', 'r0', 'r0_1']
    transcriptiome_psi = terry_data[sample_list].to_numpy()

    length = terry_data['gene length (nt)'].to_numpy().astype(np.float64)


    classified_GO = {}
    psi_stat = []

    sector = []
    for key, valus in GO_dict.items():
        sub_data = transcriptiome_psi[terry_data['ECOCYC'].isin(valus['id'])]
        # psi_sub_data = psi[terry_data['ECOCYC'].isin(valus['id'])]        classified_GO[key] = sub_data
        sector.append(key)
        psi_sector = sub_data.sum(axis=0).reshape(-1, 1)
        psi_stat.append(psi_sector)


    psi_stat = pd.DataFrame(data=np.hstack(psi_stat), columns=sector, index=sample_list)


#%%
str = "cysC cysD cysE cysH cysI cysJ cysK cysM cysN"
gene_list = str.split(' ')
eco_number = name2ecoid(gene_list)
print(eco_number)
for gene in gene_list:
    print(f'''"{gene}", ''')
