# %%
import os
import time

import requests
from requests.exceptions import SSLError
import numpy
import scipy
import pandas as pd
from bs4 import BeautifulSoup

command_lib = {
    'kb-version': 'kb-version?orgid=ECOLI',
    'web_server': 'https://websvc.biocyc.org'
}

# ====================
# Find Tables in biocyc: https://biocyc.org/gene-search.shtml
# %%


account = 'pan.chu@siat.ac.cn'
psw = 'nXqxWN7DTtnCaDP'

webServer = 'https://websvc.biocyc.org'
biocyc_web = 'https://biocyc.org'
s = requests.Session()  # create session
# Post login credentials to session:
s.post('https://websvc.biocyc.org/credentials/login/', data={'email': account, 'password': psw})
# Issue web service request:
# %%
r = s.get('https://websvc.biocyc.org/kb-version?orgid=ECOLI')

genes = s.get(f'{webServer}/ECOLI/class-instances?object=Genes')
genes_soup = BeautifulSoup(genes.text, "lxml")
tables = genes_soup.find_all('table')
for table in tables:
    if 'class' in list(table.attrs.keys()):
        if table.attrs['class'][0] == 'sortableSAQPoutputTable':
            genes_web_table = table
genes_dict = dict(
    genes_name=[],
    genes_id=[],
    genes_href=[],
)

for gene_row in genes_web_table.find_all('a'):
    gene_name = gene_row.text
    gene_href = gene_row.attrs['href']
    gene_id = gene_href.split('=')[-1]
    genes_dict['genes_name'].append(gene_name)
    genes_dict['genes_href'].append(biocyc_web + gene_href)
    genes_dict['genes_id'].append(gene_id)

genes_table = pd.DataFrame(data=genes_dict)
genes_table.to_csv(r'./exported_data/Ecoli_Genes.csv')
genes_table.to_excel(r'./exported_data/Ecoli_Genes.xlsx')

# %% Retrieve all gene info
import time
from tqdm import tqdm
from threading import Thread

info_dict = {}


def request_info(g_id, save_dict):
    retrieve_url = f'https://websvc.biocyc.org/getxml?id=ECOLI:{g_id}&detail=full'

    while True:
        try:
            save_dict[g_id] = s.get(retrieve_url).text
            break
        except requests.exceptions.SSLError:
            pass

    return None


for g_id in tqdm(genes_dict['genes_id']):
    thread_gene_retrieve = Thread(target=request_info, args=(g_id, info_dict))
    thread_gene_retrieve.start()

# back up gene info
    for g_id, g_info in info_dict.items():

        g_name = genes_table[genes_table['genes_id'] == g_id]['genes_name'].tolist()[0]
        with open(f'./exported_data/Genes_info_xml/{g_id}_{g_name}.xml', 'w') as g_xml:

            g_xml.write(g_info)


#%% load gene xlm files

genes_table = pd.read_csv(r'.\exported_data\Ecoli_Genes.csv')
genes_id = genes_table['genes_id'].to_list()
info_dict = {}

genes_file = os.scandir(r'./exported_data/Genes_info_xml/')
for gene_xml in genes_file:
    if gene_xml.name.split('.')[-1] == 'xml':
        file_name = gene_xml.name
        gene_id = file_name.split('_')[0]
        if gene_id in genes_id:
            with open(f'./exported_data/Genes_info_xml/{file_name}', 'r') as g_xml:
                info_dict[gene_id] = g_xml.read()


gene_soup = BeautifulSoup(info_dict['EG10001'], 'lxml')


