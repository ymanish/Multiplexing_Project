import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Initialise_GC_profiling import *


def density_array(T_df, to_rep, with_rep):
    exon_den_df = T_df.replace(to_rep, with_rep)
    density_list = exon_den_df.sum(axis=0)
    return density_list

def spare_matrix(file_number):

    df = pd.DataFrame()

    for seq_record in SeqIO.parse(REGION_FILE_PATH+"\/region_group_"+str(file_number)+".fasta", "fasta"):
        region = seq_record.seq[:2000]
        temp = pd.DataFrame(list(region))
        df = pd.concat([df, temp], axis=1)

    trans_df = df.T
    to_rep = ["U", "E", "D", "I", "F", "f"]

    total_trans = trans_df.shape[0]

    exon_density = density_array(trans_df, to_rep, [0, 1, 0, 0, 0, 0])
    intron_density = density_array(trans_df, to_rep, [0, 0, 0, 1, 0, 0])
    UUTR_density = density_array(trans_df, to_rep, [1, 0, 0, 0, 0, 0])
    DUTR_density = density_array(trans_df, to_rep, [0, 0, 1, 0, 0, 0])
    UF_density = density_array(trans_df, to_rep, [0, 0, 0, 0, 1, 0])
    DF_density = density_array(trans_df, to_rep, [0, 0, 0, 0, 0, 1])

    exon_density['length'] = total_trans
    intron_density['length'] = total_trans
    UUTR_density['length'] = total_trans
    DUTR_density['length'] = total_trans
    UF_density['length'] = total_trans
    DF_density['length'] = total_trans

    return exon_density, intron_density, UUTR_density, DUTR_density, UF_density, DF_density

if __name__ == "__main__":


    start = time.perf_counter()

    exon_df = pd.DataFrame()
    intron_df = pd.DataFrame()
    UUTR_df = pd.DataFrame()
    DUTR_df = pd.DataFrame()
    UF_df = pd.DataFrame()
    DF_df = pd.DataFrame()

    print('Generate the Region  Density........')

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(spare_matrix, file_number=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            exon_temp, intron_temp, UUTR_temp, DUTR_temp, UF_temp, DF_temp = j.result()

            exon_df = pd.concat([exon_df, exon_temp], axis=1)
            intron_df = pd.concat([intron_df, intron_temp], axis=1)
            UUTR_df = pd.concat([UUTR_df, UUTR_temp], axis=1)
            DUTR_df = pd.concat([DUTR_df, DUTR_temp], axis=1)
            UF_df = pd.concat([UF_df, UF_temp], axis=1)
            DF_df = pd.concat([DF_df, DF_temp], axis=1)


    exon_all = exon_df.sum(axis=1)
    intron_all = intron_df.sum(axis=1)
    UUTR_all = UUTR_df.sum(axis=1)
    DUTR_all = DUTR_df.sum(axis=1)
    UF_all = UF_df.sum(axis=1)
    DF_all = DF_df.sum(axis=1)

    exon_all_density = exon_all/exon_all['length']
    intron_all_density = intron_all/intron_all['length']
    UUTR_all_density = UUTR_all/UUTR_all['length']
    DUTR_all_density = DUTR_all/DUTR_all['length']
    UF_all_density = UF_all/UF_all['length']
    DF_all_density = DF_all/DF_all['length']

    exon_all_density.drop('length', axis=0, inplace=True)
    intron_all_density.drop('length', axis=0, inplace=True)
    UUTR_all_density.drop('length', axis=0, inplace=True)
    DUTR_all_density.drop('length', axis=0, inplace=True)
    UF_all_density.drop('length', axis=0, inplace=True)
    DF_all_density.drop('length', axis=0, inplace=True)


    densities = pd.DataFrame([exon_all_density, intron_all_density, UUTR_all_density, DUTR_all_density, UF_all_density, DF_all_density],
                             index=['exon', 'intron', 'UUTR', 'DUTR', 'UF', 'DF'])

    densities_T = densities.T
    densities_T['position'] = densities_T.index
    densities_T_melt = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)
    densities_T_melt['position'] = densities_T_melt['position']-1000
    densities_T_melt.to_csv(ELEMENT_DENSITY_FILE)


    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt, x='position', y='density', hue="variable")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%Occurrence")
    plt.title("Element Density per Base position")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(ELEMENT_DENSITY_CHART)

    print('Generate the Rolling Region Density........')
    for k in ['exon', 'intron', 'UUTR', 'DUTR', 'UF', 'DF']:

        densities_T[k] = densities_T[k].rolling(window=147,
                                                               min_periods=1,
                                                               center=True).mean()

    densities_T_melt_rolling = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt_rolling.rename({'value': 'density'}, inplace=True, axis=1)
    densities_T_melt_rolling['position'] = densities_T_melt_rolling['position'] - 1000
    densities_T_melt_rolling.to_csv(NUCLO_ELEMENT_DENSITY_FILE)


    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt_rolling, x='position', y='density', hue="variable")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%Rolling_Occurrence")
    plt.title("Element Density per Base position")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(NUCLO_ELEMENT_DENSITY_CHART)


    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')