


#%%
fime_path = r'D:\python_code\Biocyc_adapter\exported_data\ajax-replicon-genbro-colf'


# read the file\
tables_dic = {}
with open(fime_path, 'r') as file:
    # read by line
    lines = file.readlines()
    for line in lines:
        # remove \n
        line = line.replace('\n', '')
        # skip comment lines 
        if line.startswith('#'):
            continue
        # $ table header
        if line.startswith('$'):
            table_header = line.strip('$').split(' ')[0]
            tables_dic[table_header] = []
            print('Reading table:', table_header)
            continue

        # read the table content
        tables_dic[table_header].append(line.split('\t'))


        
# %%

table = tables_dic[table_name]

table_df_dict = {}

for table_name, table in tables_dic.items():
    table_df = pd.DataFrame(table[1:], columns=table[0])
    print(table_name)
    print(table_df.head())
    print('-------------------')
    table_df_dict[table_name] = table_df

# %%
# save all dataframe to xlsx
import pandas as pd
with pd.ExcelWriter('./exported_data/genome_browser_features.xlsx') as writer:
    for table_name, table_df in table_df_dict.items():
        table_df.to_excel(writer, sheet_name=table_name, index=False)

# %%
