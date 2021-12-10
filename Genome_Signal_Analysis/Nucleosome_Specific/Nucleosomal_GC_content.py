import pandas as pd
import numpy as np
from Bio import SeqIO
import concurrent.futures
import time
import matplotlib.pyplot as plt
import seaborn as sns

COL= [str(i) for i in range(2000)]
ZERO_COL = [str(i) for i in range(73)] + [str(j) for j in range(1927,2000)]

def sum_columns(s):
    ID = s['id']
    s[COL]=s[COL].rolling(window=147,  min_periods=1, center=True).sum()
    return s

def nucleosomal_GC_file(n):
    df_seq = pd.DataFrame()

    for seq_record in SeqIO.parse(r"C:\Users\maya620d\PycharmProjects\Multiplexing\input_fasta\group_"+str(n)+".fasta", "fasta"):

        header = seq_record.id.split('|')
        print(header[2])
        TRANSCRIPT_ID = header[2]
        sequ = list(seq_record.seq[:2000])
        temp_seq = pd.DataFrame(sequ, index=COL)
        temp_seq = temp_seq.T
        temp_seq['id'] = TRANSCRIPT_ID
        df_seq = pd.concat([df_seq, temp_seq], axis=0)

    df_seq = df_seq.replace(["A", "C", "G", "T", "N"], [0, 1, 1, 0, 0])
    df_seq.reset_index(inplace=True, drop=True)
    # df_seq_ = pd.DataFrame()
    df_seq_ = df_seq.apply(sum_columns, axis=1)
    del df_seq
    df_seq_[COL] = df_seq_[COL] / 147
    df_seq_[ZERO_COL] = 0

    df_seq_.to_csv(r"C:\Users\maya620d\PycharmProjects\Multiplexing\Nucleosome_GC\group_"+str(n)+".csv")

    df_seq_ = df_seq_.drop('id', axis=1)
    return df_seq_.sum(axis=0), len(df_seq_)


if __name__ == "__main__":

    start = time.perf_counter()
    len_df = 0
    df = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, 6)
        pool = [executor.submit(nucleosomal_GC_file, n=i) for i in iter_seq]
        for i in concurrent.futures.as_completed(pool):
            print('Return Value: {}', i)
            df_temp, temp_len_df = i.result()
            len_df = len_df + temp_len_df
            df = pd.concat([df, df_temp], axis=1) ##Output is Series so we have to concatenate the Series horizontally
            print(df.shape)

    df = df.reset_index(drop=True)
    print(df.head(3))
    print(len_df)
    # print(df.sum(axis=1)/len_df)
    absolute_gc = df.sum(axis=1)/len_df

    plt.figure(figsize=(20, 10))
    sns.scatterplot(x=range(-1000, 1000), y=absolute_gc)

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content")
    plt.title("")
    # plt.xticks(rotation = 90)
    plt.xticks(np.arange(-1000, 1000, 50), rotation=45)
    plt.savefig('absolute_gc_content.png')
    # plt.show()

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')