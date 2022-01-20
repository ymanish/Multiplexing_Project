import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Genome_Signal_Analysis.Initialise_Script import *

# Input_files_region = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\osativa\output_fasta"
# Input_files_seq = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\osativa\input_fasta"
#
# Results_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa"
#
# Total_Input_Files = 2
# Group = 'Plant' # 'Eukaryote'


def main(n):
    df_region = pd.DataFrame()
    df_seq = pd.DataFrame()

    for (seq_record_1, seq_record_2) in zip(
            SeqIO.parse(REGION_FILE_PATH+"\/region_group_"+str(n)+".fasta", "fasta"),
            SeqIO.parse(SEQ_FILE_PATH+"\group_"+str(n)+".fasta", "fasta")):

        region = seq_record_1.seq[:2000]
        temp_region = pd.DataFrame(list(region), index=COL)
        temp_region = temp_region.T
        temp_region['id'] = (seq_record_1.id).split('|')[0]

        df_region = pd.concat([df_region, temp_region], axis=0)

        seq = seq_record_2.seq[:2000]
        temp_seq = pd.DataFrame(list(seq), index=COL)
        temp_seq = temp_seq.T
        if GROUP == 'Plant':
            temp_seq['id'] = (seq_record_2.id).split('|')[1]
        else:
            temp_seq['id'] = (seq_record_2.id).split('|')[2]

        df_seq = pd.concat([df_seq, temp_seq], axis=0)

    df_region = df_region.reset_index(drop=True)
    df_seq = df_seq.reset_index(drop=True)

    df_seq = df_seq.drop('id', axis=1)
    df_region = df_region.drop('id', axis=1)

    return df_seq, df_region


# def density_GC_Region_and_Avg(df_region, region_code, df_seq):
#     df_region_matrix = df_region.replace(["U", "E", "D", "I", "F", "f"], region_code)
#
#     df_seq_region = df_seq.mul(df_region_matrix)
#     df_seq_region = df_seq_region.replace([1, 2, 3, 4], [0, 1, 1, 0])
#
#     region_GC_per_pos = df_seq_region.sum(axis=0)  # Total GC bases in the region per position
#     region_bases_per_pos = df_region_matrix.sum(axis=0)  # Total bases in the region per position
#     density_per_pos = region_GC_per_pos / region_bases_per_pos
#
#     region_GC = region_GC_per_pos.sum()  # Total GC bases in the region
#     region_total = region_bases_per_pos.sum()  # Total bases in the region
#     avg_region_GC = region_GC / region_total
#
#     return density_per_pos, avg_region_GC


def density_GC_Region_and_Avg(df_region, region_code, df_seq):
    df_region_matrix = df_region.replace(["U", "E", "D", "I", "F", "f"], region_code)

    print('Region dataframe size: {} and Sequence dataframe size: {}'.format(df_region.shape, df_seq.shape))
    if df_region.shape[0] != df_seq.shape[0]:
        print('Warning, Region and Seq does not match')


    df_seq_region = df_region_matrix.multiply(df_seq)

    # df_u_d_temp = df_region_matrix.merge(df_seq, on="id", how='left')
    # df_u_d_temp = df_u_d_temp.dropna(axis=0)  ###Dropping any short sequences
    # df_seq_region = pd.DataFrame()
    #
    # for pos in range(1000, 2000):
    #     df_seq_region.loc[:, str(pos)] = df_u_d_temp.loc[:, str(pos) + '_x'] * df_u_d_temp.loc[:, str(pos) + '_y']

    df_seq_region = df_seq_region.replace([1, 2, 3, 4], [0, 1, 1, 0])

    # df_region_matrix = df_region_matrix.drop('id', axis=1)

    region_GC_per_pos = df_seq_region.sum(axis=0)  # Total GC bases in the region per position
    region_bases_per_pos = df_region_matrix.sum(axis=0)  # Total bases in the region per position
    density_per_pos = region_GC_per_pos / region_bases_per_pos

    region_GC = region_GC_per_pos.sum()  # Total GC bases in the region
    region_total = region_bases_per_pos.sum()  # Total bases in the region
    avg_region_GC = region_GC / region_total

    return density_per_pos, avg_region_GC



if __name__ == "__main__":
    start = time.perf_counter()

    df_SEQ = pd.DataFrame()
    df_REGION = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            # print(f'Return Value: {i.result()}')
            temp_df_seq, temp_df_region = j.result()

            print('File Reading: Region dataframe size: {} and Sequence dataframe size: {}'.format(temp_df_region.shape, temp_df_seq.shape))
            if temp_df_region.shape[0] != temp_df_seq.shape[0]:
                print('Warning, Region and Seq does not match during File Reading for file {}'.format(j))

            df_REGION = pd.concat([df_REGION, temp_df_region], axis=0)
            df_SEQ = pd.concat([df_SEQ, temp_df_seq], axis=0)

    df_SEQ_T = df_SEQ.replace(["A", "C", "G", "T", "N"], [1, 2, 3, 4, 0])

    exon_code = [0, 1, 0, 0, 0, 0]
    intron_code = [0, 0, 0, 1, 0, 0]
    UUTR_code = [1, 0, 0, 0, 0, 0]
    DUTR_code = [0, 0, 1, 0, 0, 0]

    df_REGION.reset_index(inplace=True, drop=True)
    df_SEQ_T.reset_index(inplace=True, drop=True)
    #
    # print(df_REGION)
    # print(df_SEQ_T)
    # import sys
    # sys.exit()

    density_exon_GC, avg_exon_GC = density_GC_Region_and_Avg(df_REGION, exon_code, df_SEQ_T)
    density_intron_GC, avg_intron_GC = density_GC_Region_and_Avg(df_REGION, intron_code, df_SEQ_T)
    density_UUTR_GC, avg_UUTR_GC = density_GC_Region_and_Avg(df_REGION, UUTR_code, df_SEQ_T)
    density_DUTR_GC, avg_DUTR_GC = density_GC_Region_and_Avg(df_REGION, DUTR_code, df_SEQ_T)

    densities = pd.DataFrame([density_exon_GC, density_intron_GC, density_UUTR_GC, density_DUTR_GC],
                             index=['exon', 'intron', 'UUTR', 'DUTR'])

    densities_T = densities.T
    densities_T['position'] = densities_T.index
    densities_T_melt = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)

    avg_GC_df = pd.DataFrame({'element': ['exon', 'intron', 'UUTR', 'DUTR'],
                              'Avg_GC': [avg_exon_GC, avg_intron_GC, avg_UUTR_GC, avg_DUTR_GC]
                              })

    avg_GC_df.to_csv(AVG_GC_PER_REGION_FILE)

    densities_T_melt['position'] = densities_T_melt['position'].astype(int)
    densities_T_melt.to_csv(GC_PER_REGION_PER_BP_FILE)

    densities_T_melt['position'] = densities_T_melt['position'] - 1000

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

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')