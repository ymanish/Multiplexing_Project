import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt

def replace_chars(T_df, to_rep, with_rep):
    exon_den_df = T_df.replace(to_rep, with_rep)
    return exon_den_df

def sum_columns(s):
    return s.rolling(window=147,  min_periods=1, center=True).sum()


def main(n):
    df = pd.DataFrame()
    COL = [str(i) for i in range(853, 2000)]

    for seq_record in SeqIO.parse(
            r"C:\Users\maya620d\PycharmProjects\Multiplexing\output_fasta\region_human_" + str(n) + ".fasta", "fasta"):
        #     upstream=seq_record.seq[:1000]
        down_region = seq_record.seq[853:2000]
        temp = pd.DataFrame(list(down_region), index=COL)
        df = pd.concat([df, temp], axis=1)

    trans_df = df.T
    to_rep = ["U", "E", "D", "I", "F", "f"]

    trans_df_exon = replace_chars(trans_df, to_rep, [0, 1, 0, 0, 0, 0])
    trans_df_intron = replace_chars(trans_df, to_rep, [0, 0, 0, 1, 0, 0])
    trans_df_UUTR = replace_chars(trans_df, to_rep, [1, 0, 0, 0, 0, 0])
    trans_df_DUTR = replace_chars(trans_df, to_rep, [0, 0, 1, 0, 0, 0])

    tf_exon = (trans_df_exon.apply(sum_columns, axis=1)) / 147
    tf_intron = (trans_df_intron.apply(sum_columns, axis=1)) / 147
    tf_UUTR = (trans_df_UUTR.apply(sum_columns, axis=1)) / 147
    tf_DUTR = (trans_df_DUTR.apply(sum_columns, axis=1)) / 147

    zero_col = [str(i) for i in range(853, 926)] + [str(j) for j in range(1927, 2000)]

    tf_exon[zero_col] = 0
    tf_intron[zero_col] = 0

    tf_UUTR[zero_col] = 0
    tf_DUTR[zero_col] = 0

    return tf_exon, tf_intron, tf_UUTR, tf_DUTR

if __name__ == "__main__":

    start = time.perf_counter()

    exon_df = pd.DataFrame()
    intron_df = pd.DataFrame()
    UUTR_df = pd.DataFrame()
    DUTR_df = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, 5)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for i in concurrent.futures.as_completed(pool):
            # print(f'Return Value: {i.result()}')
            exon_temp, intron_temp, UUTR_temp, DUTR_temp = i.result()
            print('File Reading Number {}: Exon_shape: {}, Intron_shape: {}, UUTR_shape: {}, DUTR_shape: {}'.format(i, exon_temp.shape, intron_temp.shape, UUTR_temp.shape, DUTR_temp.shape))

            exon_df = pd.concat([exon_df, exon_temp], axis=0)
            intron_df = pd.concat([intron_df, intron_temp], axis=0)
            UUTR_df = pd.concat([UUTR_df, UUTR_temp], axis=0)
            DUTR_df = pd.concat([DUTR_df, DUTR_temp], axis=0)

    total_trans = len(exon_df)
    print (exon_df.shape)
    Exon_density = exon_df.sum(axis=0) / total_trans
    Intron_density = intron_df.sum(axis=0) / total_trans
    UUTR_density = UUTR_df.sum(axis=0) / total_trans
    DUTR_density = DUTR_df.sum(axis=0) / total_trans
    print (Exon_density)
    region_density = pd.DataFrame([Exon_density, Intron_density, UUTR_density, DUTR_density],
                                  index=['exon', 'intron', 'UUTR', 'DUTR'])

    region_density_T = region_density.T
    region_density_T['position'] = region_density_T.index

    region_density_T_melt = pd.melt(region_density_T, id_vars=['position'])
    region_density_T_melt.rename({'value': 'density'}, inplace=True, axis=1)
    region_density_T_melt['position'] = region_density_T_melt['position'].astype(int)

    region_density_T_melt.to_csv("elements_density.csv")

    region_density_T_melt['position'] = region_density_T_melt['position'] - 1000

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=region_density_T_melt, x='position', y='density', hue="variable")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%Occurrence")
    plt.title("")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig('density_chart.png')

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')

