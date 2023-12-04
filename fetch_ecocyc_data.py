# %%
import os
import time

import pylab as p
import requests
from requests.exceptions import SSLError
import numpy
import scipy
import pandas as pd
from bs4 import BeautifulSoup
from tqdm import tqdm
from Bio import SeqIO

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

gb_data = SeqIO.read(r'./exported_data/U00096.3.gb', 'genbank')

fet_dict = dict(genes_name=[], genes_id=[], product=[], locus_tag=[])
fet_key = dict(genes_name='gene', genes_id='db_xref', product='product', locus_tag='locus_tag')
for fet in gb_data.features:
    fet_ecocyc_id = None
    if 'gene' in fet.qualifiers.keys():

        if 'db_xref' in fet.qualifiers.keys():
            fet_db = fet.qualifiers['db_xref']
            fet_db_type = [db.split(':')[0] for db in fet_db]
            fet_db_number = [db.split(':')[1] for db in fet_db]

            for db_i, db_t in enumerate(fet_db_type):
                if db_t == 'ECOCYC':
                    fet_ecocyc_id = fet_db_number[db_i]

        fet_dict['genes_id'].append(fet_ecocyc_id)

        for key, val in fet_key.items():
            if val != 'db_xref':
                if val in fet.qualifiers.keys():
                    fet_dict[key].append(fet.qualifiers[val][0])
                else:
                    fet_dict[key].append(None)

fet_dict = pd.DataFrame(data=fet_dict)
fet_dict.to_csv(r'./exported_data/Ecoli_features.csv')
fet_dict.to_excel(r'./exported_data/Ecoli_features.xlsx')

# concat features and genes to all_genes
genes_id = genes_table['genes_id']
features_id = fet_dict['genes_id']

genes_id_set = set(genes_id)
features_id_set = set(features_id)

fet_only_id = list(features_id_set - genes_id_set)

fet_only_table = fet_dict[[True if id in fet_only_id else False for id in fet_dict['genes_id']]]
all_genes_dict = dict()
for key in ['genes_name', 'genes_id', 'genes_href']:
    all_genes_dict[key] = genes_table[key].to_list()
for id in fet_only_id:
    if id is not None:
        fet = fet_dict[fet_dict['genes_id'] == id]
        all_genes_dict['genes_id'].append(fet['genes_id'].values[0])
        all_genes_dict['genes_name'].append(fet['genes_name'].values[0])
        href = f'https://biocyc.org/gene?orgid=ECOLI&id={id}'
        all_genes_dict['genes_href'].append(href)
all_genes_table = pd.DataFrame(data=all_genes_dict)
all_genes_table.to_csv(r'./exported_data/Ecoli_all_genes.csv')
all_genes_table.to_excel(r'./exported_data/Ecoli_all_genes.xlsx')
# %% get all RNAs
genes = s.get(f'{webServer}/ECOLI/class-instances?object=RNAs')
genes_soup = BeautifulSoup(genes.text, "lxml")
tables = genes_soup.find_all('table')
for table in tables:
    if 'class' in list(table.attrs.keys()):
        if table.attrs['class'][0] == 'sortableSAQPoutputTable':
            RNAs_web_table = table
RNAs_dict = dict(
    genes_name=[],
    genes_id=[],
    genes_href=[],
    RNA_name=[],
    RNA_id=[],
    RNA_href=[],
)

for gene_row in RNAs_web_table.find_all('tr'):
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

        RNAs_dict['genes_name'].append(gene_name)
        RNAs_dict['genes_href'].append(gene_href)
        RNAs_dict['genes_id'].append(gene_id)
        RNAs_dict['RNA_name'].append(object_name)
        RNAs_dict['RNA_id'].append(object_id)
        RNAs_dict['RNA_href'].append(object_href)

genes_table = pd.DataFrame(data=RNAs_dict)
genes_table.to_csv(r'./exported_data/Ecoli_RNAs.csv')
genes_table.to_excel(r'./exported_data/Ecoli_RNAs.xlsx')

# %% Retrieve all gene info

from tqdm import tqdm
from threading import Thread
from joblib import delayed, Parallel
import time
import random
info_dict = {}

s = requests.Session()  # create session
# Post login credentials to session:
s.post('https://websvc.biocyc.org/credentials/login/', data={'email': account, 'password': psw})


def request_info(g_id, save_dict, session):
    retrieve_url = f'https://websvc.biocyc.org/getxml?id=ECOLI:{g_id}&detail=full'

    while True:
        try:
            xml_gene = session.get(retrieve_url).text
            save_dict[g_id] = xml_gene
            break
        except requests.exceptions.SSLError:
            print(g_id, 'try next time')
        except requests.exceptions.ConnectionError:
            print(g_id, 'try next time')
        except requests.exceptions.ChunkedEncodingError:
            print(g_id, 'try next time')
        time.sleep(random.random())
    return None


genes_dict = pd.read_excel(r'./exported_data/Ecoli_Genes.xlsx')
features_dict = pd.read_excel(r'./exported_data/Ecoli_features.xlsx')
genes_id = genes_dict['genes_id']
features_id = features_dict['genes_id']

genes_id_set = set(genes_id)
features_id_set = set(features_id)

all_id = list(set(genes_id_set | features_id_set))

# all_thread = []
# for g_id in tqdm(all_id):
#     thread_gene_retrieve = Thread(target=request_info, args=(g_id, info_dict, s))
#     all_thread.append(thread_gene_retrieve)

process_id = all_id.copy()
all_thread = Parallel(n_jobs=24, require='sharedmem')(delayed(request_info)(g_id, info_dict, s)
                                                     for g_id in tqdm(all_id))

print(f'Reqired {len(all_id)}, Success {len(info_dict.keys())}')
# back up gene info

for id_i, id in enumerate(all_id):
    g_name = None
    try:
        g_name = genes_dict[genes_dict['genes_id'] == id]['genes_name'].tolist()[0]
    except:
        g_name = features_dict[features_dict['genes_id'] == id]['genes_name'].tolist()[0]
    print(id, g_name, id_i)
    with open(f'./exported_data/Genes_info_xml/{id}_{g_name}.xml', 'w', encoding='utf-8') as g_xml:
        g_xml.write(info_dict[id])

# %% load gene xlm files and extract information

genes_table = pd.read_csv(r'./exported_data/Ecoli_all_genes.csv')
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
    try:
        product = gene_soup.find('product').findChild().name
    except AttributeError:
        try:
            product = gene_soup.find('product').text
        except AttributeError:
            product = None
        # some gene have no product or pseudo product.\


    try:
        transcription_direction = gene_soup.find('transcription-direction').text
    except:
        transcription_direction = None
    try:
        left_end_position = int(gene_soup.find('left-end-position').text)
    except:
        left_end_position = None
    try:
        right_end_position = int(gene_soup.find('right-end-position').text)
    except:
        right_end_position = None
    try:
        name = gene_soup.find('common-name').text
    except:
        name = None
    try:
        id = gene_soup.find('Gene').attrs['ID']
    except:
        try:
            id = gene_soup.find('Gene-Fragment').attrs['ID']
        except:
            id = None

    region = (left_end_position, right_end_position)
    synonym = gene_soup.find_all('synonym')
    if len(synonym) >= 1:
        synonym = '; '.join([s.text for s in synonym])
    else:
        synonym = None
    regulator = gene_soup.find('regulated-by')
    if regulator:
        regulation = regulator.find_all('Regulation')
        regulation = '; '.join([r.attrs['frameid'] for r in regulation])
    else:
        regulation = None

    component = gene_soup.find('component-of')
    if component:
        component = component.find_all('Transcription-Unit')
        component = '; '.join([c.attrs['frameid'] for c in component])
    else:
        component = None
    locus_tage = gene_soup.find('accession-1')
    if locus_tage:
        locus_tage = locus_tage.text
    else:
        locus_tage = None

    info_dict = {
        'id': id,
        'name': name,
        'transcription-direction': transcription_direction,
        'product': product,
        'region': region,
        'synonym': synonym,
        'regulation': regulation,
        'component': component,
        'locus_tag': locus_tage
    }
    return info_dict


# test get_gene_info() function
# get_gene_info(info_dict['EG10001'])
#


all_gene_id = list(info_dict.keys())
gene_id_number = len(all_gene_id)
genes_dict = {
    'id': [None] * gene_id_number,
    'name': [None] * gene_id_number,
    'transcription-direction': [None] * gene_id_number,
    'product': [None] * gene_id_number,
    'region': [None] * gene_id_number,
    'synonym': [None] * gene_id_number,
    'regulation': [None] * gene_id_number,
    'component': [None] * gene_id_number,
    'locus_tag': [None] * gene_id_number
}
for id_i, id in enumerate(tqdm(all_gene_id)):
    for prop_key, prop in get_gene_info(info_dict[id]).items():
        genes_dict[prop_key][id_i] = prop

genes_info_table = pd.DataFrame(data=genes_dict)
genes_info_table.to_csv(r'./exported_data/Ecoli_gene_info.csv')
genes_info_table.to_excel(r'./exported_data/Ecoli_gene_info.xlsx')

# %% get all DNA/RNA sites data
# https://ecocyc.org/ECOLI/search-query?type=Sites&siteTypes=ATTENUATORS%20CRYPTIC-PROPHAGES%20MRNA-SITES%20ORIGIN%20OTHER-EXTRAGENIC%20PHAGE-SITES%20PROMOTERS%20RECOMB-SITES%20REPEATS%20RIBOSWITCHES%20TERMINATORS%20TF-SITES%20TUS%20TRANSPOSONS
DNA_RNA_url = "https://ecocyc.org/ECOLI/search-query?type=Sites&siteTypes=ATTENUATORS%20CRYPTIC-PROPHAGES%20MRNA-SITES%20ORIGIN%20OTHER-EXTRAGENIC%20PHAGE-SITES%20PROMOTERS%20RECOMB-SITES%20REPEATS%20RIBOSWITCHES%20TERMINATORS%20TF-SITES%20TUS%20TRANSPOSONS"

data = s.get(DNA_RNA_url).text
xml_data = BeautifulSoup(data, 'html').find_all('table')
xml_table = None
for xml in xml_data:
    if 'class' in list(xml.attrs.keys()):
        if xml.attrs['class'][0] == 'sortableSAQPoutputTable':
            xml_table = xml
sites = []
for row in xml_table.find_all('tr')[1:]:
    hints = row.find_all('td')
    if len(hints) > 0:

        hit_name = '; '.join([reg.text for reg in hints[0].find_all('a')])
        hit_id = '; '.join([reg.attrs['href'].split('=')[-1] for reg in hints[0].find_all('a')])
        hit_type = hints[1].text
        try:
            left = int(hints[2].text)
            right = int(hints[3].text)
            region = (left, right)
        except ValueError:
            region = None
        strand = hints[4].text
        regulation = hints[5].find_all('a')
        if len(regulation) > 0:
            regulation_name = '; '.join([reg.text for reg in regulation])
            regulation_id = '; '.join([reg.attrs['href'].split('=')[-1] for reg in regulation])
        else:
            regulation_name = None
            regulation_id = None

        sites.append(dict(site_name=hit_name, site_id=hit_id, site_type=hit_type, region=region,
                          strand=strand, transcript_name=regulation_name, transcript_id=regulation_id))

sites_number = len(sites)
all_sites = dict(site_name=[None] * sites_number, site_id=[None] * sites_number, site_type=[None] * sites_number,
                 region=[None] * sites_number,
                 strand=[None] * sites_number, transcript_name=[None] * sites_number,
                 transcript_id=[None] * sites_number)
for s_i, site in enumerate(sites):
    for key, val in site.items():
        all_sites[key][s_i] = val

DNA_RNA_sites_df = pd.DataFrame(data=all_sites)
DNA_RNA_sites_df.to_csv(r'./exported_data/Ecoli_sites_info.csv')
DNA_RNA_sites_df.to_excel(r'./exported_data/Ecoli_sites_info.xlsx')
