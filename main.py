# %%
import requests
import numpy
import scipy
import pandas as pd
from bs4 import BeautifulSoup

command_lib = {
    'kb-version': 'kb-version?orgid=ECOLI',
    'web_server': 'https://websvc.biocyc.org'
}
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
    genes_dict['genes_href'].append(biocyc_web+gene_href)
    genes_dict['genes_id'].append(gene_id)

genes_table = pd.DataFrame(data=genes_dict)
genes_table.to_csv(r'./exported_data/Ecoli_Genes.csv')
genes_table.to_excel(r'./exported_data/Ecoli_Genes.xlsx')