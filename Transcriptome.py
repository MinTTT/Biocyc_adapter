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

terry_sci_data = pd.read_excel(r'./transcriptiome_data/science.abk2066_table_s3.xlsx', '2 - RNAseq-ss (fractions)')

terry_gene_set = set(terry_sci_data['gene'].to_list())
terry_locus_set = set(terry_sci_data['locus'].to_list())

ecocyc_gene_info = pd.read_csv(r'./exported_data/Ecoli_gene_info.csv')
all_gene_info = pd.read_csv(r'./exported_data/Ecoli_all_genes.csv')

all_features = pd.read_csv(r'./exported_data/Ecoli_features.csv')

ecocyc_gene_set = set(ecocyc_gene_info['name'].to_list())
ecocyc_locus_set = set(ecocyc_gene_info['locus_tag'].to_list())
terry_only_set = terry_locus_set - ecocyc_locus_set
ecocyc_only_set = ecocyc_locus_set - terry_locus_set
terry_only_name = terry_sci_data[ [True if locus in terry_only_set else False for locus in terry_sci_data['locus']]]['gene']



ecocyc_only_table = all_features[[True if locus in ecocyc_only_set else False for locus in all_features['locus_tag']] ]

terry_only_table = all_features[[True if locus in terry_only_set else False
                                 for locus in all_features['locus_tag']] ]
terry_only_table_by_name = all_features[[True if name in terry_only_name else False
                                         for name in all_features['genes_name']] ]
