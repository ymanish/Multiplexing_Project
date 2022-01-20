import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Genome_Signal_Analysis.Initialise_Script import *

# Group = 'Plant' # 'Eukaryote'
# Total_Input_Files = 3
#
# COL = [str(i) for i in range(2000)]
#
# Input_seq_file = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\osativa\input_fasta"
# Results_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa"

def main(n):
    df_seq = pd.DataFrame()

    for seq_record in SeqIO.parse(SEQ_FILE_PATH+"/group_"+str(n)+".fasta", "fasta"):
        header = seq_record.id.split('|')
        if GROUP == 'Plant':
            print(header[1])
            TRANSCRIPT_ID = header[1]
        else:
            print(header[2])
            TRANSCRIPT_ID = header[2]

        sequ = list(seq_record.seq[:2000])
        temp_seq = pd.DataFrame(sequ, index=COL)
        temp_seq = temp_seq.T
        temp_seq['id'] = TRANSCRIPT_ID
        df_seq = pd.concat([df_seq, temp_seq], axis=0)

    df_seq = df_seq.replace(["A", "C", "G", "T", "N"], [0, 1, 1, 0, 0])
    df_seq.reset_index(inplace=True, drop=True)
    df_seq = df_seq.drop('id', axis=1)

    return df_seq.sum(axis=0), len(df_seq)

if __name__ == "__main__":
    start = time.perf_counter()
    len_df = 0
    df = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):

            df_temp, temp_len_df = j.result()
            len_df = len_df + temp_len_df
            df = pd.concat([df, df_temp], axis=1)  ##Output is Series so we have to concatenate the Series horizontally
            print(df.shape)

    df = df.reset_index(drop=True)
    print(df.head(3))
    print(len_df)
    absolute_gc = df.sum(axis=1) / len_df
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

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')

