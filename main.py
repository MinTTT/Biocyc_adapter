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
# https://biocyc.org/getxml?id=ECOLI:PM0-1352&detail=low
# https://biocyc.org/getxml?ECOLI:RPOH-MONOMER
# %%


account = 'pan.chu@siat.ac.cn'
psw = 'nXqxWN7DTtnCaDP'

webServer = 'https://websvc.biocyc.org'
biocyc_web = 'https://biocyc.org'
s = requests.Session()  # create session
# Post login credentials to session:
s.post('https://websvc.biocyc.org/credentials/login/', data={'email': account, 'password': psw})
# Issue web service request:
r = s.get('https://websvc.biocyc.org/kb-version?orgid=ECOLI')
# %% get all gene
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
# %% get all RNAs
genes = s.get(f'{webServer}/ECOLI/class-instances?object=RNAs')
genes_soup = BeautifulSoup(genes.text, "xml")
tables = genes_soup.find_all('table')
for table in tables:
    if 'class' in list(table.attrs.keys()):
        if table.attrs['class'][0] == 'sortableSAQPoutputTable':
            genes_web_table = table
genes_dict = dict(
    genes_name=[],
    genes_id=[],
    genes_href=[],
    RNA_name=[],
    RNA_id=[],
    RNA_href=[],
)

for gene_row in genes_web_table.find_all('tr'):
    # gene_row = gene_row.find_all('td')
    td = gene_row.find_all('a')
    if len(td) != 0:
        if len(td) == 1:
            obj = td[0]
            gene_name = None
            gene_href = None
            gene_id = None

        else:
            obj, gene = td
            gene_name = gene.text
            gene_href = biocyc_web + gene.attrs['href']
            gene_id = gene_href.split('=')[-1]
        object_name = obj.text
        object_href = biocyc_web + obj.attrs['href']
        object_id = object_href.split('=')[-1]

        genes_dict['genes_name'].append(gene_name)
        genes_dict['genes_href'].append(gene_href)
        genes_dict['genes_id'].append(gene_id)
        genes_dict['RNA_name'].append(object_name)
        genes_dict['RNA_id'].append(object_id)
        genes_dict['RNA_href'].append(object_href)


genes_table = pd.DataFrame(data=genes_dict)
genes_table.to_csv(r'./exported_data/Ecoli_RNAs.csv')
genes_table.to_excel(r'./exported_data/Ecoli_RNAs.xlsx')

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

# %% load gene xlm files

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

gene_soup = BeautifulSoup(info_dict['EG10001'], 'xml')


def get_gene_info(gene_xml):
    """
    product
    transcription-direction
    synonym: <synonym>
    region: <left_end_position, right_end_position>
    name: <common-name>
    Parameters
    ----------
    gene_xml :

    Returns
    -------

    """
    gene_soup = BeautifulSoup(gene_xml, 'xml')
    product = gene_soup.find('product').findChild().name
    transcription_direction = gene_soup.find('transcription-direction').text
    left_end_position = int(gene_soup.find('left-end-position').text)
    right_end_position = int(gene_soup.find('right-end-position').text)
    region = (left_end_position, right_end_position)
    synonym = gene_soup.find('synonym').text
    info_dict = {
        'transcription-direction': transcription_direction,
        'product': product,
        'region': region,
        'synonym': synonym
    }
    return info_dict


get_gene_info(info_dict['EG10001'])
