import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import matplotlib.pyplot as plt

element_GC_density_perbasepose_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\NUC_element_GC_density_perbasepos.csv"
element_density_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\NUC_elements_density.csv"
Avg_GC_per_element_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\NUC_Avg_GC_per_element.csv"

Results_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa"


Element_GC_content = pd.read_csv(element_GC_density_perbasepose_path)
Element_density = pd.read_csv(element_density_path)
Avg_GC_density = pd.read_csv(Avg_GC_per_element_path)

Element_GC_content.drop(['Unnamed: 0'], axis=1, inplace=True)
Element_density.drop(['Unnamed: 0'], axis=1, inplace=True)
Avg_GC_density.drop(['Unnamed: 0'], axis=1, inplace=True)


df_GC_Element = Element_GC_content.merge(Element_density, on=['position', 'variable'], how='left')
df_GC_Element['density'] = df_GC_Element['density_x']*df_GC_Element['density_y']

plt.figure(figsize=(20, 10))
sns.scatterplot(data=df_GC_Element, x='position', y="density", hue='variable')

plt.xlabel("Position w.r.t TSS")
plt.ylabel("% Nucleosomal GC content * % Nucleosomal Occurrence")
plt.title("%GC distribution by Region per base position, with Average across ALL the Transcripts")
plt.xticks(rotation=0)
# plt.show()
plt.savefig(Results_path+'\Charts\/NUC_Chart4_GC_Density_and_element_density.png')


absolute_GC = df_GC_Element.groupby('position').agg({'density':sum}).reset_index()

plt.figure(figsize=(20, 10))
sns.scatterplot(data=absolute_GC, x='position', y="density")

plt.xlabel("Position w.r.t TSS")
plt.ylabel("Sum(% Nucleosomal GC content * % Nucleosomal Occurrence)")
plt.title("Absolute GC content averaged across ALL transcripts")
plt.xticks(rotation=0)
# plt.show()
plt.savefig(Results_path+'\Charts\/NUC_Chart6_Derived_Absolute_GC.png')


def multiply_operation(el, density, tab):
    return density*tab.loc[tab['element']==el, 'Avg_GC'].values[0]


Element_density['Homoginized_density_GC'] = Element_density.apply(lambda x: multiply_operation(x['variable'],
                                                                                               x['density'], Avg_GC_density),
                                                                  axis=1)

plt.figure(figsize=(20, 10))
sns.scatterplot(data=Element_density, x='position', y="Homoginized_density_GC", hue='variable')

plt.xlabel("Position w.r.t TSS")
plt.ylabel("% Nucleosomal GC content * % Nucleosomal Occurrence")
plt.title("Average %GC distribution by Region per base position")
plt.xticks(rotation=0)
# plt.show()
plt.savefig(Results_path+'\Charts\/NUC_Chart5_Avg_GC_Density_and_element_density.png')

AVG_absolute_GC = Element_density.groupby('position').agg({'Homoginized_density_GC':sum}).reset_index()

plt.figure(figsize=(20, 10))
sns.scatterplot(data=AVG_absolute_GC, x='position', y="Homoginized_density_GC")

plt.xlabel("Position w.r.t TSS")
plt.ylabel("Sum(% Nucleosomal GC content * % Nucleosomal Occurrence)")
plt.title("Average %GC distribution per base position")
plt.xticks(rotation=0)
# plt.show()
plt.savefig(Results_path+'\Charts\/NUC_Chart7_Derived_Avg_Absolute_GC.png')
