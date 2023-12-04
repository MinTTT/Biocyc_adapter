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
ecocyc_gene_set = set(ecocyc_gene_info['name'].to_list())
ecocyc_locus_set = set(ecocyc_gene_info['locus_tag'].to_list())
terry_only = terry_locus_set - ecocyc_locus_set
