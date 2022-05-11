import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt

#
# Input_files = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\osativa\output_fasta"
# Results_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa"
#
Total_Input_Files = 1
FILE_PATH = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\h_sapiens\ucsc_fasta"
REPEAT_DENSITY_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\repeat.csv"
REPEAT_DENSITY_CHART= r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\repeat.png"

def density_array(T_df, to_rep, with_rep):
    exon_den_df = T_df.replace(to_rep, with_rep)
    density_list = exon_den_df.sum(axis=0)
    return density_list

def spare_matrix(file_number):

    df = pd.DataFrame()

    for seq_record in SeqIO.parse(FILE_PATH+"\/group_"+str(file_number)+".fasta", "fasta"):
        #     upstream=seq_record.seq[:1000]
        # down_region = seq_record.seq[1000:2000]
        region = seq_record.seq[:2000]
        temp = pd.DataFrame(list(region))
        df = pd.concat([df, temp], axis=1)

    trans_df = df.T
    to_rep = ['A', 'T', 'G', 'C', 'N', 'a', 'g', 'c', 't', 'n']
    total_trans = trans_df.shape[0]

    repeats_density = density_array(trans_df, to_rep, [0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    a_densities = density_array(trans_df, to_rep, [0, 0, 0, 0, 0, 1, 0, 0, 0, 0])
    t_densities = density_array(trans_df, to_rep, [0, 0, 0, 0, 0, 0, 0, 0, 1, 0])
    c_densities = density_array(trans_df, to_rep, [0, 0, 0, 0, 0, 0, 0, 1, 0, 0])
    g_densities = density_array(trans_df, to_rep, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0])

    print (repeats_density)
    repeats_density['length'] = total_trans
    a_densities['length'] = total_trans
    t_densities['length'] = total_trans
    c_densities['length'] = total_trans
    g_densities['length'] = total_trans

    return repeats_density, a_densities, t_densities, c_densities, g_densities

if __name__ == "__main__":


    start = time.perf_counter()

    repeat_df = pd.DataFrame()
    a_repeat_df = pd.DataFrame()
    t_repeat_df = pd.DataFrame()
    c_repeat_df = pd.DataFrame()
    g_repeat_df = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(spare_matrix, file_number=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            #print(f'Return Value: {j.result()}')
            repeat_temp, a_repeat_temp, t_repeat_temp, c_repeat_temp, g_repeat_temp = j.result()

            repeat_df = pd.concat([repeat_df, repeat_temp], axis=1)
            a_repeat_df = pd.concat([a_repeat_df, a_repeat_temp], axis=1)
            t_repeat_df = pd.concat([t_repeat_df, t_repeat_temp], axis=1)
            c_repeat_df = pd.concat([c_repeat_df, c_repeat_temp], axis=1)
            g_repeat_df = pd.concat([g_repeat_df, g_repeat_temp], axis=1)


    # print(exon_df)
    # print(exon_df.shape)
    # print(exon_df)
    # print(exon_df.sum(axis=1))

    repeat_all = repeat_df.sum(axis=1)
    a_repeat_all = a_repeat_df.sum(axis=1)
    t_repeat_all = t_repeat_df.sum(axis=1)
    c_repeat_all = c_repeat_df.sum(axis=1)
    g_repeat_all = g_repeat_df.sum(axis=1)

    repeat_all_density = repeat_all/repeat_all['length']
    a_repeat_all_density = a_repeat_all/a_repeat_all['length']
    t_repeat_all_density = t_repeat_all/t_repeat_all['length']
    c_repeat_all_density = c_repeat_all/c_repeat_all['length']
    g_repeat_all_density = g_repeat_all/g_repeat_all['length']


    repeat_all_density.drop('length', axis=0, inplace=True)
    a_repeat_all_density.drop('length', axis=0, inplace=True)
    t_repeat_all_density.drop('length', axis=0, inplace=True)
    c_repeat_all_density.drop('length', axis=0, inplace=True)
    g_repeat_all_density.drop('length', axis=0, inplace=True)

    print (repeat_all_density.to_list())

    densities = pd.DataFrame([repeat_all_density, a_repeat_all_density, t_repeat_all_density, c_repeat_all_density, g_repeat_all_density],
                             index=['R', 'A', 'T', 'C', 'G'])

    densities_T = densities.T
    densities_T['position'] = densities_T.index
    densities_T_melt = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)
    densities_T_melt.to_csv(REPEAT_DENSITY_FILE)
    #
    #
    densities_T_melt['position'] = densities_T_melt['position']-1000


    print (densities_T_melt)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt, x='position', y='density', hue="variable")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%Occurrence")
    plt.title("Repeats Density per Base position")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(REPEAT_DENSITY_CHART)

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')