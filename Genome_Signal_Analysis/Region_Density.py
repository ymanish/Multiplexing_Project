import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt

def density_array(T_df, to_rep, with_rep):
    exon_den_df = T_df.replace(to_rep, with_rep)
    density_list = exon_den_df.sum(axis=0)
    return density_list

def spare_matrix(file_number):

    df = pd.DataFrame()

    for seq_record in SeqIO.parse(r"C:\Users\maya620d\Documents\Notebooks\Output_files\region_human_"+str(file_number)+".fasta", "fasta"):
        #     upstream=seq_record.seq[:1000]
        down_region = seq_record.seq[1000:2000]
        temp = pd.DataFrame(list(down_region))
        df = pd.concat([df, temp], axis=1)

    trans_df = df.T
    to_rep = ["U", "E", "D", "I", "F", "f"]

    total_trans = trans_df.shape[0]

    exon_density = density_array(trans_df, to_rep, [0, 1, 0, 0, 0, 0])
    intron_density = density_array(trans_df, to_rep, [0, 0, 0, 1, 1, 1])
    UUTR_density = density_array(trans_df, to_rep, [1, 0, 0, 0, 0, 0])
    DUTR_density = density_array(trans_df, to_rep, [0, 0, 1, 0, 0, 0])

    exon_density['length'] = total_trans
    intron_density['length'] = total_trans
    UUTR_density['length'] = total_trans
    DUTR_density['length'] = total_trans

    return exon_density, intron_density, UUTR_density, DUTR_density

if __name__ == "__main__":


    start = time.perf_counter()

    exon_df = pd.DataFrame()
    intron_df = pd.DataFrame()
    UUTR_df = pd.DataFrame()
    DUTR_df = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(0,2000,1000)
        pool = [executor.submit(spare_matrix, file_number=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            #print(f'Return Value: {j.result()}')
            exon_temp, intron_temp, UUTR_temp, DUTR_temp = j.result()

            exon_df = pd.concat([exon_df, exon_temp], axis=1)
            intron_df = pd.concat([intron_df, intron_temp], axis=1)
            UUTR_df = pd.concat([UUTR_df, UUTR_temp], axis=1)
            DUTR_df = pd.concat([DUTR_df, DUTR_temp], axis=1)

            #f.write(i.result())

    # f.close()
    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')

    # print(exon_df)
    print(exon_df.shape)
    print(exon_df)
    print(exon_df.sum(axis=1))

    exon_all = exon_df.sum(axis=1)
    intron_all = intron_df.sum(axis=1)
    UUTR_all = UUTR_df.sum(axis=1)
    DUTR_all = DUTR_df.sum(axis=1)

    exon_all_density = exon_all/exon_all['length']
    intron_all_density = intron_all/intron_all['length']
    UUTR_all_density = UUTR_all/UUTR_all['length']
    DUTR_all_density = DUTR_all/DUTR_all['length']

    exon_all_density.drop('length', axis=0, inplace=True)
    intron_all_density.drop('length', axis=0, inplace=True)
    UUTR_all_density.drop('length', axis=0, inplace=True)
    DUTR_all_density.drop('length', axis=0, inplace=True)

    x = range(1000)

    densities = pd.DataFrame([exon_all_density, intron_all_density, UUTR_all_density, DUTR_all_density],
                             index=['exon', 'intron', 'UUTR', 'DUTR'])

    densities_T = densities.T
    densities_T['position'] = densities_T.index
    densities_T_melt = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)

    densities_T_melt.to_csv("./elements_density.csv")


    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt, x='position', y='density', hue="variable")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%Occurrence")
    plt.title("")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig('density_chart.png')

