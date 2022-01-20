import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Initialise_Script import *


# output_file_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\osativa\output_fasta"
# Results_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa"
# Total_Input_Files = 1
# COL = [str(i) for i in range(2000)]
#
#
# element_GC_density_perbasepose_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\MRNA_element_GC_density_perbasepos.csv"
# element_density_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\MRNA_elements_density.csv"
# Avg_GC_per_element_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\MRNA_Avg_GC_per_element.csv"
#
# Absolute_GC_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa\Files\MRNA_absolute_gc_content.csv"



def density_GC_Region_and_Avg(df_region, region_code, df_seq):
    df_region_matrix = df_region.replace(["U", "E", "D", "I", "F", "f"], region_code)

    print('Region dataframe size: {} and Sequence dataframe size: {}'.format(df_region.shape, df_seq.shape))
    if df_region.shape[0] != df_seq.shape[0]:
        print('Warning, Region and Seq does not match')


    df_seq_region = df_region_matrix.multiply(df_seq)

    df_seq_region = df_seq_region.replace([1, 2, 3, 4], [0, 1, 1, 0])

    region_GC_per_pos = df_seq_region.sum(axis=0)  # Total GC bases in the region per position
    region_bases_per_pos = df_region_matrix.sum(axis=0)  # Total bases in the region per position
    density_per_pos = region_GC_per_pos / region_bases_per_pos

    region_GC = region_GC_per_pos.sum()  # Total GC bases in the region
    region_total = region_bases_per_pos.sum()  # Total bases in the region
    avg_region_GC = region_GC / region_total

    return density_per_pos, avg_region_GC


def density_array(T_df, to_rep, with_rep):
    exon_den_df = T_df.replace(to_rep, with_rep)
    density_list = exon_den_df.sum(axis=0)
    return density_list


def Region_Density(DF_region):
    to_rep = ["U", "E", "D", "I", "F", "f"]

    total_trans = DF_region.shape[0]

    exon_density = density_array(DF_region, to_rep, [0, 1, 0, 0, 0, 0])
    UUTR_density = density_array(DF_region, to_rep, [1, 0, 0, 0, 0, 0])
    DUTR_density = density_array(DF_region, to_rep, [0, 0, 1, 0, 0, 0])

    exon_density = exon_density / total_trans
    UUTR_density = UUTR_density / total_trans
    DUTR_density = DUTR_density / total_trans

    densities = pd.DataFrame([exon_density, UUTR_density, DUTR_density],
                             index=['exon', 'UUTR', 'DUTR'])

    densities_T = densities.T

    densities_T['position'] = densities_T.index
    densities_T_melt = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)

    densities_T_melt['position'] = densities_T_melt['position'].astype(int)
    densities_T_melt['position'] = densities_T_melt['position'] - 1000

    densities_T_melt.to_csv(ELEMENT_DENSITY_FILE)


    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt, x='position', y='density', hue="variable")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%Occurrence")
    plt.title("MRNA: Element Density per Base position")
    plt.xticks(rotation=0)
    plt.savefig(ELEMENT_DENSITY_CHART)

    return None


def Regional_and_Avg_GC(DF_seq, DF_region):
    DF_SEQ_T = DF_seq.replace(["A", "C", "G", "T", "N"], [1, 2, 3, 4, 0])

    exon_code = [0, 1, 0, 0, 0, 0]
    UUTR_code = [1, 0, 0, 0, 0, 0]
    DUTR_code = [0, 0, 1, 0, 0, 0]

    DF_REGION.reset_index(inplace=True, drop=True)
    DF_SEQ_T.reset_index(inplace=True, drop=True)

    density_exon_GC, avg_exon_GC = density_GC_Region_and_Avg(DF_region, exon_code, DF_SEQ_T)
    density_UUTR_GC, avg_UUTR_GC = density_GC_Region_and_Avg(DF_region, UUTR_code, DF_SEQ_T)
    density_DUTR_GC, avg_DUTR_GC = density_GC_Region_and_Avg(DF_region, DUTR_code, DF_SEQ_T)

    densities = pd.DataFrame([density_exon_GC, density_UUTR_GC, density_DUTR_GC],
                             index=['exon', 'UUTR', 'DUTR'])

    densities_T = densities.T
    densities_T['position'] = densities_T.index
    densities_T_melt = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)

    avg_GC_df = pd.DataFrame({'element': ['exon', 'UUTR', 'DUTR'],
                              'Avg_GC': [avg_exon_GC, avg_UUTR_GC, avg_DUTR_GC]
                              })

    avg_GC_df.to_csv(AVG_GC_PER_REGION_FILE)

    densities_T_melt['position'] = densities_T_melt['position'].astype(int)
    densities_T_melt['position'] = densities_T_melt['position'] - 1000
    densities_T_melt.to_csv(GC_PER_REGION_PER_BP_FILE)

    plt.figure(figsize=(10, 10))
    sns.barplot(x="element", y="Avg_GC", data=avg_GC_df)

    plt.xlabel("Region")
    plt.ylabel("%GC content")
    plt.title("Average GC content per Region")
    # plt.xticks(rotation = 90)
    # plt.xticks(np.arange(-1000, 1000, 50), rotation=45)
    # plt.show()
    plt.savefig(AVG_GC_PER_REGION_CHART)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt, x='position', y="density", hue='variable')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content")
    plt.title("%GC content distribution by Region per base position, with average across Region specific Transcripts")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(GC_PER_REGION_PER_BP_CHART)
    return None


def main(n):
    df_seq = pd.DataFrame()
    df_region = pd.DataFrame()

    for seq_record in SeqIO.parse(REGION_FILE_PATH+"/mrna_group_"+str(n)+".fasta", "fasta"):
        header = seq_record.id
        TRANSCRIPT_ID = header[0]
        region_seq = list(seq_record.seq.split('|')[0][:2000])
        mrna_seq = list(seq_record.seq.split('|')[1][:2000])
        # print (mrna_seq)
        temp_seq = pd.DataFrame(mrna_seq, index=COL)
        temp_seq = temp_seq.T
        temp_seq['id'] = TRANSCRIPT_ID
        df_seq = pd.concat([df_seq, temp_seq], axis=0)

        temp_region = pd.DataFrame(region_seq, index=COL)
        temp_region = temp_region.T
        temp_region['id'] = TRANSCRIPT_ID
        df_region = pd.concat([df_region, temp_region], axis=0)


    # df_seq = df_seq.replace(["A", "C", "G", "T", "N"], [0, 1, 1, 0, 0])
    df_seq = df_seq.reset_index(drop=True)
    df_seq = df_seq.drop('id', axis=1)

    df_region = df_region.reset_index(drop=True)
    df_region = df_region.drop('id', axis=1)

    return df_seq, df_region

def multiply_operation(el, density, tab):
    return density*tab.loc[tab['element']==el, 'Avg_GC'].values[0]

def derived_Charts():
    Element_GC_content = pd.read_csv(GC_PER_REGION_PER_BP_FILE)
    Element_density = pd.read_csv(ELEMENT_DENSITY_FILE)
    Avg_GC_density = pd.read_csv(AVG_GC_PER_REGION_FILE)

    Element_GC_content = Element_GC_content.drop(['Unnamed: 0'], axis=1)
    Element_density = Element_density.drop(['Unnamed: 0'], axis=1)
    Avg_GC_density = Avg_GC_density.drop(['Unnamed: 0'], axis=1)

    df_GC_Element = Element_GC_content.merge(Element_density, on=['position', 'variable'], how='left')
    df_GC_Element['density'] = df_GC_Element['density_x'] * df_GC_Element['density_y']

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=df_GC_Element, x='position', y="density", hue='variable')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content * %Occurrence")
    plt.title("%GC distribution by Region per base position, with Average across ALL the Transcripts")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(ABS_GC_DENSITY_BY_REGION_CHART)

    absolute_GC = df_GC_Element.groupby('position').agg({'density': sum}).reset_index()

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=absolute_GC, x='position', y="density")

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Sum(%GC content * %Occurrence)")
    plt.title("Absolute GC content averaged across ALL transcripts")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(DERIVED_ABSOLUTE_GC_CHART)

    Element_density['Homoginized_density_GC'] = Element_density.apply(lambda x: multiply_operation(x['variable'],
                                                                                                   x['density'],
                                                                                                   Avg_GC_density),
                                                                      axis=1)
    Element_density.to_csv(AVG_GC_DENSITY_BY_REGION_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=Element_density, x='position', y="Homoginized_density_GC", hue='variable')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content * %Occurrence")
    plt.title("Average %GC distribution by Region per base position")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(AVG_GC_DENSITY_BY_REGION_CHART)

    AVG_absolute_GC = Element_density.groupby('position').agg({'Homoginized_density_GC': sum}).reset_index()

    AVG_absolute_GC.to_csv(DERIVED_AVG_GC_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=AVG_absolute_GC, x='position', y="Homoginized_density_GC")

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Sum(%GC content * %Occurrence)")
    plt.title("Average %GC distribution per base position")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(DERIVED_AVG_GC_CHART)

    return df_GC_Element, Avg_GC_density



def create_intra_regional_charts(DF_element_GC_den, DF_Avg_GC_density ):
    DF_element_GC_den = DF_element_GC_den.rename({'variable': 'element',
                                          'density_x': 'element_GC_density',
                                          'density_y': 'element_density'},
                                         axis=1)

    DF = DF_element_GC_den.merge(DF_Avg_GC_density,
                             on='element',
                             how='left')

    DF['intra_signal'] = DF['element_GC_density'] - DF['Avg_GC']

    DF.to_csv(CENTERED_GC_PER_REGION_PER_BP_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=DF, x='position', y="intra_signal", hue='element')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("mean-centered %GC")
    plt.title("Mean centered %GC content distribution by Region per base position, with average across Region specific Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(CENTERED_GC_PER_REGION_PER_BP_CHART)

    # plt.show()


    DF['intra_signal'] = DF['intra_signal'] * DF['element_density']
    DF.to_csv(CENTERED_GC_DENSITY_BY_REGION_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=DF, x='position', y="intra_signal", hue='element')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("mean-centered %GC")
    plt.title("Mean centered %GC distribution by Region per base position, with Average across ALL the Transcripts")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(CENTERED_GC_DENSITY_BY_REGION_CHART)

    intra_absolute_GC = DF.groupby('position').agg({'intra_signal': sum}).reset_index()

    intra_absolute_GC.to_csv(DERIVED_CENTERED_GC_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=intra_absolute_GC, x='position', y="intra_signal")

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("mean-centered %GC")
    plt.title("Mean Centered Absolute %GC content averaged across ALL transcripts")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(DERIVED_CENTERED_GC_CHART)

    DF['inter_signal'] = DF['Avg_GC'] * DF['element_density']
    return DF


if __name__ == "__main__":
    start = time.perf_counter()
    DF_SEQ = pd.DataFrame()
    DF_REGION = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            df_temp_seq, df_temp_region = j.result()

            DF_SEQ = pd.concat([DF_SEQ, df_temp_seq],axis=0)
            DF_REGION = pd.concat([DF_REGION, df_temp_region], axis=0)


    DF_SEQ.reset_index(inplace=True, drop=True)
    DF_REGION.reset_index(inplace=True, drop=True)

    DF_SEQ_ = DF_SEQ.replace(["A", "C", "G", "T", "N"], [0, 1, 1, 0, 0])

    absolute_gc = DF_SEQ_.sum(axis=0)/len(DF_SEQ_)
    absolute_gc.to_csv(ABSOLUTE_GC_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(x=range(-1000, 1000), y=absolute_gc)

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content")
    plt.title("GC content per Base position")
    # plt.xticks(rotation = 90)
    plt.xticks(np.arange(-1000, 1000, 50), rotation=45)
    # plt.show()
    plt.savefig(ABSOLUTE_GC_CHART)


    Region_Density(DF_REGION)
    Regional_and_Avg_GC(DF_SEQ, DF_REGION)


    ### Derived Charts
    df_GC_Element, Avg_GC_density =  derived_Charts()
    DF_intra = create_intra_regional_charts(df_GC_Element, Avg_GC_density)

    absolute_GC_file = pd.read_csv(ABSOLUTE_GC_FILE)

    absolute_GC_file.drop(['Unnamed: 0'], axis=1, inplace=True)

    absolute_GC_file.rename({'0': 'Absolute_GC'},
                            axis=1,
                            inplace=True)

    absolute_GC_file['position'] = absolute_GC_file.index
    absolute_GC_file['position'] = absolute_GC_file['position'] - 1000

    DF_absolute = DF_intra.merge(absolute_GC_file,
                                 on='position',
                                 how='left')

    DF_absolute['Derived_GC'] = DF_absolute['intra_signal'] + DF_absolute['inter_signal']

    df_ = DF_absolute.groupby('position').agg({'Derived_GC': sum,
                                               'Absolute_GC': 'first'}).reset_index()

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
    plt.savefig(COMPARISON_DERIVED_ABS_GC)



    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')