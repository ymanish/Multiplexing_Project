import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import matplotlib.pyplot as plt

element_GC_density_perbasepose_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\NUC_element_GC_density_perbasepos.csv"
element_density_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\NUC_elements_density.csv"
Avg_GC_per_element_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\NUC_Avg_GC_per_element.csv"

Absolute_GC_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\NUC_absolute_gc_content.csv"

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
Element_density.to_csv(Results_path + '\Files\Temp_Chart_5.csv')

plt.figure(figsize=(20, 10))
sns.scatterplot(data=Element_density, x='position', y="Homoginized_density_GC", hue='variable')

plt.xlabel("Position w.r.t TSS")
plt.ylabel("% Nucleosomal GC content * % Nucleosomal Occurrence")
plt.title("Average %GC distribution by Region per base position")
plt.xticks(rotation=0)
# plt.show()
plt.savefig(Results_path+'\Charts\/NUC_Chart5_Avg_GC_Density_and_element_density.png')

AVG_absolute_GC = Element_density.groupby('position').agg({'Homoginized_density_GC':sum}).reset_index()
AVG_absolute_GC.to_csv(Results_path + '\Files\Temp_Chart_7.csv')

plt.figure(figsize=(20, 10))
sns.scatterplot(data=AVG_absolute_GC, x='position', y="Homoginized_density_GC")

plt.xlabel("Position w.r.t TSS")
plt.ylabel("Sum(% Nucleosomal GC content * % Nucleosomal Occurrence)")
plt.title("Average %GC distribution per base position")
plt.xticks(rotation=0)
# plt.show()
plt.savefig(Results_path+'\Charts\/NUC_Chart7_Derived_Avg_Absolute_GC.png')


def create_intra_regional_charts(DF_element_GC_den, DF_Avg_GC_density ):
    DF_element_GC_den = DF_element_GC_den.rename({'variable': 'element',
                                          'density_x': 'element_GC_density',
                                          'density_y': 'element_density'},
                                         axis=1)

    DF = DF_element_GC_den.merge(DF_Avg_GC_density,
                             on='element',
                             how='left')

    DF['intra_signal'] = DF['element_GC_density'] - DF['Avg_GC']
    DF.to_csv(Results_path+'\Files\Temp_Chart_8.csv')

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=DF, x='position', y="intra_signal", hue='element')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("mean-centered %GC")
    plt.title("Mean centered %GC content distribution by Region per base position, with average across Region specific Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(Results_path + '\Charts\/NUC_Chart8_Centered_GC_Region_density.png')

    # plt.show()

    DF['intra_signal'] = DF['intra_signal'] * DF['element_density']
    DF.to_csv(Results_path+'\Files\Temp_Chart_9.csv')

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=DF, x='position', y="intra_signal", hue='element')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("mean-centered %GC")
    plt.title("Mean centered %GC distribution by Region per base position, with Average across ALL the Transcripts")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(Results_path + '\Charts\/NUC_Chart9_Centered_GC_density_and element_density.png')

    intra_absolute_GC = DF.groupby('position').agg({'intra_signal': sum}).reset_index()
    intra_absolute_GC.to_csv(Results_path+'\Files\Temp_Chart_10.csv')

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=intra_absolute_GC, x='position', y="intra_signal")

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("mean-centered %GC")
    plt.title("Mean Centered Absolute %GC content averaged across ALL transcripts")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(Results_path + '\Charts\/NUC_Chart10_Centered_Derived_Absolute_GC.png')

    DF['inter_signal'] = DF['Avg_GC'] * DF['element_density']
    return DF


DF_intra = create_intra_regional_charts(df_GC_Element, Avg_GC_density)
absolute_GC_file = pd.read_csv(Absolute_GC_path)


absolute_GC_file.drop(['Unnamed: 0'], axis=1, inplace=True)

absolute_GC_file.rename({'0':'Absolute_GC'},
                        axis=1,
                        inplace=True)

absolute_GC_file['position'] = absolute_GC_file.index
absolute_GC_file['position'] = absolute_GC_file['position']-1000

DF_absolute = DF_intra.merge(absolute_GC_file,
                         on='position',
                         how='left')

DF_absolute['Derived_GC'] = DF_absolute['intra_signal'] + DF_absolute['inter_signal']

df_ = DF_absolute.groupby('position').agg({'Derived_GC':sum,
                                   'Absolute_GC':'first' }).reset_index()

df_melt = pd.melt(df_,
                  id_vars=['position'],
                  value_vars=['Derived_GC', 'Absolute_GC'])

plt.figure(figsize=(20, 10))
sns.scatterplot(data=df_melt, x='position', y="value", hue='variable')
plt.xlabel("Position w.r.t TSS")
plt.ylabel("%GC")
plt.title("Comparison of Derived GC and Absolute GC")
plt.xticks(rotation=0)
# plt.show()
plt.savefig(Results_path + '\Charts\/NUC_Chart11_Derived_and_Absolute_GC.png')