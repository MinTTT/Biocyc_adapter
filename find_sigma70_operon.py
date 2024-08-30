

#%%
import numpy as np
import pandas as pd
import requests
from requests.exceptions import SSLError
from bs4 import BeautifulSoup



account = 'pan.chu@siat.ac.cn'
psw = 'nXqxWN7DTtnCaDP'

webServer = 'https://websvc.biocyc.org'
biocyc_web = 'https://biocyc.org'
s = requests.Session()  # create session
# Post login credentials to session:
s.post('https://websvc.biocyc.org/credentials/login/', data={'email': account, 'password': psw})
# Issue web service request:
r = s.get('https://websvc.biocyc.org/kb-version?orgid=ECOLI')
s = requests.Session()  # create session

#%%
# read all tables in the file
fime_path = r'.\exported_data\genome_browser_features.xlsx'

tables_dic = pd.read_excel(fime_path, sheet_name=None)
print(list(tables_dic.keys()))

table_names = ['Replicon', 'Transcription-Units', 'GENES', 'PROMOTERS', 'Terminators', 'Attenuators', 'Transcription-Factor-Binding-Sites', 'mRNA-Binding-Sites', 'Ribosome-Binding-Sites', 'Extragenic-Sites', 'IS-Elements']

sigma70_table = tables_dic['PROMOTERS']
sigma70_table = sigma70_table[sigma70_table['SIGMA-FACTOR'] == '70']
sigma70_table = sigma70_table[sigma70_table['HIGH-QUALITY-EVIDENCE?'] == 'T']
# %%
# promoter details

promoter_unregulated = []
promoter_regulated = []
xml_files_dir = './exported_data/Promoter_xml'
xml_list = os.listdir(xml_files_dir)
xml_name = [xml_file.strip('.xml') for xml_file in xml_list ]

for promoter_id in sigma70_table['*UNIQUE-ID']:
    # check if the promoter was qureied before
    if promoter_id in xml_name:
        data = open(os.path.join(xml_files_dir, promoter_id + '.xml'), 'r').read()
    else:
        quar = F'https://biocyc.org/getxml?ECOLI:{promoter_id}&detail=full' 
        try:
            data = s.get(quar).text
            with open(os.path.join(xml_files_dir, promoter_id + '.xml'), 'w') as file:
                file.write(data)
        except SSLError:
            print(f"SSL error in {promoter_id}")
            continue
    id_soup = BeautifulSoup(data, 'xml')
    regulation = id_soup.find_all('regulator')
    try:
        regulation_id = regulation[0].find_all()[0].attrs['frameid']
        print(f"{promoter_id} is regulated by {regulation[0].find_all()[0].attrs['frameid']}")
    except IndexError:
        promoter_unregulated.append(promoter_id)
        print(f"{promoter_id} is not regulated")

    belonging = id_soup.find_all('component-of')
    try:
        belonging_id = belonging[0].find('Transcription-Unit').attrs['frameid']
        print(f"{promoter_id} is component of {regulation[0].find_all()[0].attrs['frameid']}")
    except IndexError:
        pass
    # belonging_id = belonging[0].find('Transcription-Unit').attrs['frameid']
    # print(f"{promoter_id} is regulated by {regulation[0].find_all()[0].attrs['frameid']}")  

# promoter_id = sigma70_table.iloc[0, 0]
# promoter_id

# quar = 'https://biocyc.org/getxml?ECOLI:' + promoter_id
# data = s.get(quar).text
# id_soup = BeautifulSoup(data, 'xml')

#%%
regulation
# %%
sigma70_table
# %%
