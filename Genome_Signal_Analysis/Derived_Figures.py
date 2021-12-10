import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import matplotlib.pyplot as plt



Element_GC_content = pd.read_csv('GC_element_density.csv')
Element_density = pd.read_csv('../elements_density.csv')
Avg_GC_density = pd.read_csv('Avg_GC_element_density.csv')

Element_GC_content.drop(['Unnamed: 0'], axis=1, inplace=True)

Element_density.drop(['Unnamed: 0'], axis=1, inplace=True)

Avg_GC_density.drop(['Unnamed: 0'], axis=1, inplace=True)


df_GC_Element =  Element_GC_content.merge(Element_density, on=['position', 'variable'], how='left')
df_GC_Element['density'] = df_GC_Element['density_x']*df_GC_Element['density_y']

plt.figure(figsize=(20, 10))
sns.scatterplot(data=df_GC_Element, x='position', y="density", hue='variable')

plt.xlabel("Position w.r.t TSS")
plt.ylabel("Density_GC*Density_Element")
plt.title("")
plt.xticks(rotation=0)
# plt.show()
plt.savefig('Chart4_GC_element.png')


absolute_GC = df_GC_Element.groupby('position').agg({'density':sum}).reset_index()

plt.figure(figsize=(20, 10))
sns.scatterplot(data=absolute_GC, x='position', y="density")

plt.xlabel("Position w.r.t TSS")
plt.ylabel("Sum(Density_GC*Density_Element)")
plt.title("")
plt.xticks(rotation=0)
# plt.show()
plt.savefig('Chart6_Absolute_GC.png')


def multiply_operation(el, density, tab):
    return density*tab.loc[tab['element']==el, 'Avg_GC'].values[0]


Element_density['Homoginized_density_GC'] = Element_density.apply(lambda x: multiply_operation(x['variable'],
                                                                                               x['density'], Avg_GC_density),
                                                                  axis=1)

plt.figure(figsize=(20, 10))
sns.scatterplot(data=Element_density, x='position', y="Homoginized_density_GC", hue='variable')

plt.xlabel("Position w.r.t TSS")
plt.ylabel("AVG_Density_GC*Density_Element")
plt.title("")
plt.xticks(rotation=0)
# plt.show()
plt.savefig('Chart5_avg_GC_element.png')

AVG_absolute_GC = Element_density.groupby('position').agg({'Homoginized_density_GC':sum}).reset_index()

plt.figure(figsize=(20, 10))
sns.scatterplot(data=AVG_absolute_GC, x='position', y="Homoginized_density_GC")

plt.xlabel("Position w.r.t TSS")
plt.ylabel("Sum(Avg_Density_GC*Density_Element)")
plt.title("")
plt.xticks(rotation=0)
# plt.show()
plt.savefig('Chart7_Avg_Absolute_GC.png')
