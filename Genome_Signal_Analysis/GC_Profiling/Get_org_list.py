import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 100)
from pybiomart import Server

# plants.ensembl.org for plants
# http://www.ensembl.org for eukaryote

server = Server(host='http://plants.ensembl.org')
server.verbose = True
print(server.list_marts())

mart = server['plants_mart']
# print(type(mart))
df = mart.list_datasets()
print(df)
df['id'] = df.index
df['id'] = df['id'] +1

df = df[['id', 'name', 'display_name']]
df.to_csv('all_plant_list.txt', index=False, header=False)
