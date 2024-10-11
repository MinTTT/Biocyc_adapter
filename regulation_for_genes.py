# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
#%%
# Built-in/Generic Imports
import os
import sys
# […]

# Libs
import pandas as pd
import numpy as np  # Or any other
# […]

# Own modules
from fetch_ecocyc_data import BioCycQuery
import re

#%%
# fetch data from ecocyc
ecocyc = BioCycQuery()
# get all genes from table


fime_path = r'.\exported_data\genome_browser_features.xlsx'

tables_dic = pd.read_excel(fime_path, sheet_name=None)
print(list(tables_dic.keys()))

table_names = ['Replicon', 'Transcription-Units', 'GENES', 'PROMOTERS', 'Terminators', 'Attenuators',
               'Transcription-Factor-Binding-Sites', 'mRNA-Binding-Sites', 'Ribosome-Binding-Sites',
               'Extragenic-Sites', 'IS-Elements']
genes_table = tables_dic['GENES']
genes_id_list = genes_table['*UNIQUE-ID']
# ['*UNIQUE-ID', 'ACCESSION-1', 'NAME', 'SYNONYMS', 'PRODUCT-NAME',
#        'SWISS-PROT-ID', 'START-BASE', 'END-BASE', 'DIRECTION', 'GENE-TYPE',
#        'TU-COLOR']

#%%
genes_xmls = ecocyc.query_by_IDs(list(genes_id_list))
#%% find genes are belong to which TU
# ===========
# component
#   |____ Promoter
#           |____ binds-sigma-factor
#           |____ regulated-by
# -----------
genesTU_list = []
for gene_xml in genes_xmls:
    # gene_xml = genes_xmls[0]
    component_of_list = gene_xml.find_all('Transcription-Unit')

    TU_list = [component_of.attrs['frameid'] for component_of in component_of_list
               if component_of.attrs['frameid'].startswith('TU')]
    genesTU_list.append(TU_list)

#%%
genesPromoters_list = []
genesSigma_list = []
genesRegulator_activate_list = []
genesRegulator_inhibit_list = []

for gene_index, gene_xml in enumerate(genes_xmls):
    TU_list = genesTU_list[gene_index]
    promoter_list = []
    sigma_list = []
    regulated_by_list = []
    if TU_list:
        # get TU xmls
        TU_xmls = ecocyc.query_by_IDs(TU_list)
        # find promoters
        for TU_xml in TU_xmls:
            promoters_xml = TU_xml.find_all('Promoter')
            promoters_id = [promoter.attrs['frameid'] for promoter in promoters_xml]
            promoter_list.append(promoters_id)
            for promt_xml in promoters_xml:
                sigma_xml = promt_xml.find('binds-sigma-factor')
                sigma_id = sigma_xml.find().attrs['frameid'] if sigma_xml else None
                sigma_list.append(sigma_id)

        TUs_info = ecocyc.query_info_by_IDs(TU_list)
        regulators_name = []
        regulators_activate_name = []
        regulators_inhibit_name = []
        for TU_info in TUs_info:
            Tu_name = TU_info[list(TU_info.keys())[0]]['names'].strip('<B>Operon:</B>')
            try:
                Tu_regulators_inhibit = TU_info[list(TU_info.keys())[0]]['inhibitors']
                # match the characters in parentheses: \(([^)]+)\), and remove <[^>]*>
                Tu_regulators_inhibit = re.findall(r'\(([^)]+)\)', Tu_regulators_inhibit)
                Tu_regulators_inhibit = [re.sub(r'<[^>]*>', '', _inhibit)
                                         for _inhibit in Tu_regulators_inhibit]
                Tu_regulators_inhibit = '; '.join(Tu_regulators_inhibit)
            except KeyError:
                Tu_regulators_inhibit = None
            try:
                Tu_regulators_activate = TU_info[list(TU_info.keys())[0]]['activators']
                Tu_regulators_activate = re.findall(r'\(([^)]+)\)', Tu_regulators_activate)
                Tu_regulators_activate = [re.sub(r'<[^>]*>', '', _activate)
                                         for _activate in Tu_regulators_activate]
                Tu_regulators_activate = '; '.join(Tu_regulators_activate)
            except KeyError:
                Tu_regulators_activate = None
            regulators_activate_name.append(Tu_regulators_activate)
            regulators_inhibit_name.append(Tu_regulators_inhibit)


    genesSigma_list.append(sigma_list)
    genesPromoters_list.append(promoter_list)
    genesRegulator_activate_list.append(regulators_activate_name)
    genesRegulator_inhibit_list.append(regulators_inhibit_name)

#%%
genes_mat = genes_table.to_numpy()
genes_mat = np.hstack([genes_mat, np.array(genesTU_list, dtype=object).reshape(-1, 1),
                       np.array(genesPromoters_list, dtype=object).reshape(-1, 1),
                       np.array(genesSigma_list, dtype=object).reshape(-1, 1),
                       np.array(genesRegulator_activate_list, dtype=object).reshape(-1, 1),
                       np.array(genesRegulator_inhibit_list, dtype=object).reshape(-1, 1)])

genes_pd = pd.DataFrame(genes_mat,
                        columns=list(genes_table.columns) +
                                ['TU', 'Promoters', 'Sigma', 'Regulator_activate', 'Regulator_inhibit'])

# save genes_pd to excel and csv
genes_pd.to_excel(r'./exported_data/Ecoli_genes_info.xlsx', index=False)
genes_pd.to_csv(r'./exported_data/Ecoli_genes_info.csv', index=False)