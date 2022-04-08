import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 100)
from pybiomart import Server

server = Server(host='http://www.ensembl.org')
server.verbose = True
print(server.list_marts())

mart = server['ENSEMBL_MART_ENSEMBL']
# print(type(mart))
df = mart.list_datasets()
df.to_csv('all_org_list.txt')
