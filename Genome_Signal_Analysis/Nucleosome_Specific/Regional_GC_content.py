import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt

def sum_columns(s):
    return s.rolling(window=147,  min_periods=1, center=True).sum()

def replace_chars(T_df, to_rep, with_rep):
    exon_den_df = T_df.replace(to_rep, with_rep)
    return exon_den_df


def density_per_Region_per_position(trans_df_region, trans_df_seq_, rep_region):
    trans_df_region_ = replace_chars(trans_df_region, ["U", "E", "D", "I", "F", "f"], rep_region)
    trans_df_region_seq = trans_df_region_.mul(trans_df_seq_)
    trans_df_region_seq_gc = trans_df_region_seq.apply(sum_columns, axis=1)

    zero_col = [str(i) for i in range(853, 926)] + [str(j) for j in range(1927, 2000)]
    trans_df_region_seq_gc[zero_col] = 0
    trans_df_region_seq_gc = trans_df_region_seq_gc.reset_index(drop=True)

    trans_df_region_total = trans_df_region_.apply(sum_columns, axis=1)
    trans_df_region_total = trans_df_region_total.reset_index(drop=True)
    trans_df_region_total[zero_col] = 0

    trans_Region = trans_df_region_seq_gc.sum(axis=0) / trans_df_region_total.sum(axis=0)
    trans_Region = trans_Region.fillna(0)

    Avg_GC_Region = trans_df_region_seq_gc.sum(axis=0).sum() / trans_df_region_total.sum(axis=0).sum()

    return trans_Region, Avg_GC_Region



def main(n):
    df_region = pd.DataFrame()
    df_seq = pd.DataFrame()
    COL = [str(i) for i in range(853, 2000)]

    for (seq_record_1, seq_record_2) in zip(SeqIO.parse(r"C:\Users\maya620d\PycharmProjects\Multiplexing\output_fasta\region_human_"+str(n)+".fasta", "fasta"),
                                            SeqIO.parse(r"C:\Users\maya620d\PycharmProjects\Multiplexing\input_fasta\group_"+str(n)+".fasta", "fasta")):
        #     upstream=seq_record.seq[:1000]
        down_region = seq_record_1.seq[853:2000]
        temp_1 = pd.DataFrame(list(down_region), index=COL)
        df_region = pd.concat([df_region, temp_1], axis=1)

        down_sequence = seq_record_2.seq[853:2000]
        temp_2 = pd.DataFrame(list(down_sequence), index=COL)
        df_seq = pd.concat([df_seq, temp_2], axis=1)

    T_df_region = df_region.T
    T_df_seq = df_seq.T

    return T_df_seq, T_df_region




if __name__ == "__main__":
    start = time.perf_counter()

    # df_U_SEQ = pd.DataFrame()
    df_D_SEQ = pd.DataFrame()
    df_D_REGION = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, 5)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            # print(f'Return Value: {i.result()}')
            temp_df_d_seq, temp_df_d_region = j.result()

            print('File Reading: Region dataframe size: {} and Sequence dataframe size: {}'.format(temp_df_d_region.shape, temp_df_d_seq.shape))
            if temp_df_d_region.shape[0] != temp_df_d_seq.shape[0]:
                print('Warning, Region and Seq does not match during File Reading for file {}'.format(j))

            df_D_REGION = pd.concat([df_D_REGION, temp_df_d_region], axis=0)
            df_D_SEQ = pd.concat([df_D_SEQ, temp_df_d_seq], axis=0)
            # df_U_SEQ = pd.concat([df_U_SEQ, temp_df_u_seq], axis=0)

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')

    with_rep_seq = [0, 1, 1, 0, 0]

    df_D_SEQ_T = replace_chars(df_D_SEQ, ["A", "C", "G", "T", "N"], with_rep_seq)

    exon_code = [0, 1, 0, 0, 0, 0]
    intron_code = [0, 0, 0, 1, 0, 0]
    UUTR_code = [1, 0, 0, 0, 0, 0]
    DUTR_code = [0, 0, 1, 0, 0, 0]

    GC_Exon_density_per_position, Avg_GC_Exon = density_per_Region_per_position(df_D_REGION, df_D_SEQ_T,
                                                                                exon_code)
    GC_Intron_density_per_position, Avg_GC_Intron = density_per_Region_per_position(df_D_REGION, df_D_SEQ_T,
                                                                               intron_code)
    GC_UUTR_density_per_position, Avg_GC_UUTR = density_per_Region_per_position(df_D_REGION, df_D_SEQ_T,
                                                                               UUTR_code)
    GC_DUTR_density_per_position, Avg_GC_DUTR = density_per_Region_per_position(df_D_REGION, df_D_SEQ_T,
                                                                               DUTR_code)

    region_density = pd.DataFrame([GC_Exon_density_per_position, GC_Intron_density_per_position, GC_UUTR_density_per_position, GC_DUTR_density_per_position],
                                  index=['exon', 'intron', 'UUTR', 'DUTR'])

    region_density_T = region_density.T
    region_density_T['position'] = region_density_T.index

    region_density_T_melt = pd.melt(region_density_T, id_vars=['position'])
    region_density_T_melt.rename({'value': 'density'}, inplace=True, axis=1)
    region_density_T_melt['position'] = region_density_T_melt['position'].astype(int)

    region_density_T_melt.to_csv("GC_density_Element_wise.csv")

    region_density_T_melt['position'] = region_density_T_melt['position'] - 1000

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=region_density_T_melt, x='position', y='density', hue="variable")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC_content")
    plt.title("")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig('GC_density_Element_wise_chart.png')

    avg_GC_df = pd.DataFrame({'element': ['exon', 'intron', 'UUTR', 'DUTR'],
                              'Avg_GC': [Avg_GC_Exon, Avg_GC_Intron, Avg_GC_UUTR, Avg_GC_DUTR]
                              })

    avg_GC_df.to_csv('Avg_GC_element_density.csv')

    plt.figure(figsize=(10, 10))
    sns.barplot(x="element", y="Avg_GC", data=avg_GC_df)

    plt.xlabel("Region")
    plt.ylabel("Avg GC content")
    plt.title("")
    # plt.xticks(rotation = 90)
    # plt.xticks(np.arange(-1000, 1000, 50), rotation=45)
    # plt.show()
    plt.savefig('Avg_GC_Region.png')
