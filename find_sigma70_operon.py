
#%%
import numpy as np
import pandas as pd
import requests
from requests.exceptions import SSLError
from bs4 import BeautifulSoup
import os


account = 'pan.chu@siat.ac.cn'
psw = 'nXqxWN7DTtnCaDP'

webServer = 'https://websvc.biocyc.org'
biocyc_web = 'https://biocyc.org'
s = requests.Session()  # create session
# Post login credentials to session:
s.post('https://websvc.biocyc.org/credentials/login/', data={'email': account, 'password': psw})
# Issue web service request:
r = s.get('https://websvc.biocyc.org/kb-version?orgid=ECOLI')
print(r.text)

#%%
# read all tables in the file
fime_path = r'.\exported_data\genome_browser_features.xlsx'

tables_dic = pd.read_excel(fime_path, sheet_name=None)
print(list(tables_dic.keys()))

table_names = ['Replicon', 'Transcription-Units', 'GENES', 'PROMOTERS', 'Terminators', 'Attenuators', 'Transcription-Factor-Binding-Sites', 'mRNA-Binding-Sites', 'Ribosome-Binding-Sites', 'Extragenic-Sites', 'IS-Elements']

sigma_type = '70' # type can be set  {'19', '24', '28', '32', '54', '70', 'S', nan}
sigma70_table = tables_dic['PROMOTERS']
sigma70_table = sigma70_table[sigma70_table['SIGMA-FACTOR'] == sigma_type]
sigma70_table = sigma70_table[sigma70_table['HIGH-QUALITY-EVIDENCE?'] == 'T']
# %%
# promoter details
promoter_id_list = []
regulation_list = []
promoter_unregulated = []
TU_id = []
Tu_genes = []
xml_files_dir = './exported_data/Promoter_xml'
xml_list = os.listdir(xml_files_dir)
xml_name = [xml_file.strip('.xml') for xml_file in xml_list ]
tu_xml_list = os.listdir('./exported_data/Transcription_Unit_xml')
tu_xml_name = [xml_file.strip('.xml') for xml_file in tu_xml_list]
table_length = len(sigma70_table)
for pro_i, promoter_id in enumerate(sigma70_table['*UNIQUE-ID']):
    # check if the promoter was qureied before
    print(f"{pro_i}/{table_length}: {promoter_id}")
    promoter_id_list.append(promoter_id)
    if promoter_id in xml_name:
        data = open(os.path.join(xml_files_dir, promoter_id + '.xml'), 'r').read()
    else:
        quar = f'https://websvc.biocyc.org/getxml?ECOLI:{promoter_id}&detail=full'
        print(quar)
        try:
            data = s.get(quar).text
            with open(os.path.join(xml_files_dir, promoter_id + '.xml'), 'w') as file:
                file.write(data)
        except SSLError:
            print(f"SSL error in {promoter_id}")
            continue
    id_soup = BeautifulSoup(data, 'xml')
    regulation = id_soup.find_all('regulator')
    if len(regulation) == 0:
        regulation_list.append([])
    else:
        regulation_id = [reg.attrs['frameid'] for reg in regulation[0].find_all()]
        print(f"{promoter_id} is regulated by {' '.join(regulation_id)}")
        regulation_list.append(regulation_id)


    # try:
    #     regulation_id = regulation[0].find_all()[0].attrs['frameid']
    #     print(f"{promoter_id} is regulated by {regulation[0].find_all()[0].attrs['frameid']}")
    #     regulation_list.append(regulation_id)
    # except IndexError:
    #     promoter_unregulated.append(promoter_id)
    #     print(f"{promoter_id} is not regulated")
    #     regulation_list.append(None)

    belonging = id_soup.find_all('component-of')
    try:
        belonging_id = belonging[0].find('Transcription-Unit')
        if belonging_id:
            belonging_id = belonging_id.attrs['frameid']
        else:
            belonging_id = None
        print(f"{promoter_id} is component of {belonging_id}")

    except IndexError:
        belonging_id = None
    TU_id.append(belonging_id)
    if belonging_id is not None:
        if belonging_id in tu_xml_name:
            TU_data = open(os.path.join('./exported_data/Transcription_Unit_xml', belonging_id + '.xml'), 'r').read()
        else:
            quar = f'https://websvc.biocyc.org/getxml?ECOLI:{belonging_id}&detail=full'
            try:
                TU_data = s.get(quar).text
                with open(os.path.join('./exported_data/Transcription_Unit_xml', belonging_id + '.xml'), 'w') as file:
                    file.write(data)
            except SSLError:
                print(f"SSL error in {belonging_id}")
                continue
        TU_soup = BeautifulSoup(TU_data, 'xml')
        # find all genes
        genes = TU_soup.find_all('Gene')
        genes_id = [gene.attrs['frameid'] for gene in genes]
        Tu_genes.append(genes_id)
    else:
        Tu_genes.append([])
    # belonging_id = belonging[0].find('Transcription-Unit').attrs['frameid']
    # print(f"{promoter_id} is regulated by {regulation[0].find_all()[0].attrs['frameid']}")  
#%%
# all genes in tables
all_genes = []
for genes in Tu_genes:
    all_genes += genes
all_genes = list(set(all_genes))

#
# promoter_id_list = []
# regulation_list = []
# promoter_unregulated = []
# TU_id = []
# Tu_genes = []
# find genes' promoter
gene_promoter_list = []
for gene in all_genes:
    _promoter_list = []
    _Tu_list = []
    _regulation_list = []
    for genes_i, genes in enumerate(Tu_genes):
        if gene in genes:
            _promoter_list.append(promoter_id_list[genes_i])
            _Tu_list.append(TU_id[genes_i])
            _regulation_list.append(regulation_list[genes_i])
    gene_promoter_list.append([gene, _promoter_list, _Tu_list, _regulation_list])

gene_promoter_array = np.array(gene_promoter_list, dtype=object)

def check_list_empty(list):
    if len(list) == 0:
        return True
    else:
        for item in list:
            if len(item) > 0:
                return False
        return True

unregulated_list = []
for i in range(len(gene_promoter_array)):
    _reg_list = gene_promoter_array[i, 3]
    if check_list_empty(_reg_list):
        print(f"{gene_promoter_array[i, 0]} has no regulation")
        unregulated_list.append(i)

unreg_gene_promoter_array = gene_promoter_array[unregulated_list, :]
#%% create a list for unregulated genes
tu_xml_list = os.listdir('./exported_data/Transcription_Unit_xml')
tu_xml_name = [xml_file.strip('.xml') for xml_file in tu_xml_list]
gene_xml_files_dir = '.\exported_data\Genes_info_xml'
gene_xml_list = os.listdir(gene_xml_files_dir)
gene_xml_names = [xml_file.strip('.xml') for xml_file in gene_xml_list]

product_list = []
gene_url = []
gene_name_list = []
promoters_list = []
sigma_factors_list = []
regulators_list = []
binding_list = []
for gene_id in unreg_gene_promoter_array[:, 0]:
    if gene_id in gene_xml_names:
        _xml_data = open(os.path.join(gene_xml_files_dir, gene_id + '.xml'), 'r').read()
    else:
        quar = f'https://websvc.biocyc.org/getxml?ECOLI:{gene_id}&detail=full'
        try:
            _xml_data = s.get(quar).text
            with open(os.path.join(gene_xml_files_dir, gene_id + '.xml'), 'w') as file:
                file.write(_xml_data)
        except SSLError:
            print(f"SSL error in {gene_id}")
            continue
    _gene_soup = BeautifulSoup(_xml_data, 'xml')
    try:
        _product = _gene_soup.find_all('product')[0].find_all()[0].name
    except:
        _product = _gene_soup.find_all('product')[0].text
    _gene_name = _gene_soup.find_all('common-name')[0].text
    _TU_IDs = _gene_soup.find_all('component-of')[0].find_all('Transcription-Unit')
    _TU_IDs = [TU.attrs['frameid'] for TU in _TU_IDs]
    _promoters = []
    _sigma_factors = []
    _regulators = []
    _binding_sites = []
    print(_TU_IDs)
    # find promoters and its sigma factors
    for TU_ID in _TU_IDs:
        # get xml of TU
        if TU_ID in tu_xml_name:
            TU_data = open(os.path.join('./exported_data/Transcription_Unit_xml', TU_ID + '.xml'), 'r').read()
        else:
            quar = f'https://websvc.biocyc.org/getxml?ECOLI:{TU_ID}&detail=full'
            try:
                TU_data = s.get(quar).text
                with open(os.path.join('./exported_data/Transcription_Unit_xml', TU_ID + '.xml'), 'w') as file:
                    file.write(TU_data)
            except SSLError:
                print(f"SSL error in {TU_ID}")
                continue
        TU_soup = BeautifulSoup(TU_data, 'xml')
        # find promoters
        _component = TU_soup.find_all('component')
        for comp in _component:
            if comp.find('Promoter'):
                promoter_id = comp.find('Promoter').attrs['frameid']
                _promoters.append(promoter_id)
                promoters_table = tables_dic['PROMOTERS']
                _sigma_factor = promoters_table[promoters_table['*UNIQUE-ID'] == promoter_id]['SIGMA-FACTOR'].values[0]
                _sigma_factors.append(_sigma_factor)
                promoter_pos = promoters_table[promoters_table['*UNIQUE-ID'] == promoter_id]['POSITION'].values[0]
                transcrbinding_table = tables_dic['Transcription-Factor-Binding-Sites']
                trancrbbinding_loc = transcrbinding_table['CENTER'].values
                distence = np.abs(trancrbbinding_loc - promoter_pos)
                filtered_binding = transcrbinding_table['TF-NAME'][distence <= 50].to_list()
                _binding_sites += filtered_binding
        _binding_sites = list(set(_binding_sites))
        # find regulator
        _regulator = TU_soup.find_all('regulator')
        for reg in _regulator:
            _regulator_id = [reg.attrs['frameid'] for reg in reg.find_all()]
            _regulators.append(_regulator_id)

    _url = f'https://ecocyc.org/gene?orgid=ECOLI&id={gene_id}'
    product_list.append(_product)
    gene_url.append(_url)
    gene_name_list.append(_gene_name)
    promoters_list.append(_promoters)
    sigma_factors_list.append(_sigma_factors)
    regulators_list.append(_regulators)
    binding_list.append(_binding_sites)
unregulated_dataFrame = pd.DataFrame(data=[list(unreg_gene_promoter_array[:, 0]),
                                           gene_name_list,
                                           product_list,
                                           gene_url, promoters_list, sigma_factors_list,
                                           regulators_list, binding_list],
                                     index=['GeneID', 'GeneName' ,'Product', 'URL', 'Promoters', 'SigmaFactors', 'Regulators',
                                            'Transcription-Factor-Binding-Sites']).T
#%%
