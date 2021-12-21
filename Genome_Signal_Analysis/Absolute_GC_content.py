import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt

Group = 'Plant' # 'Eukaryote'
Total_Input_Files = 3

COL = [str(i) for i in range(2000)]

Input_seq_file = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\osativa\input_fasta"
Results_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa"

def main(n):
    df_seq = pd.DataFrame()

    for seq_record in SeqIO.parse(Input_seq_file+"/group_"+str(n)+".fasta", "fasta"):
        header = seq_record.id.split('|')
        if Group == 'Plant':
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



###Calculate the absolute GC content
def absolute_gc_content(df_u, df_d):
    df_u_seq_temp = df_u.replace(["A", "C", "G", "T", "N"], [0, 1, 1, 0, 0])
    df_d_seq_temp = df_d.replace(["A", "C", "G", "T", "N"], [0, 1, 1, 0, 0])

    d_seq_temp_col = [str(i) for i in range(0, 1000)]
    d_seq_temp_col.append('id_down')

    u_seq_temp_col = [str(i) for i in range(-1000, 0)]
    u_seq_temp_col.append('id_up')

    df_u_seq_temp.columns = u_seq_temp_col
    df_d_seq_temp.columns = d_seq_temp_col

    df_u_d_temp = df_u_seq_temp.merge(df_d_seq_temp, left_on="id_up", right_on="id_down", how='left')

    print(df_u_d_temp.shape)
    df_u_d_temp = df_u_d_temp.dropna(axis=0)  ###Dropping any short sequences
    print(df_u_d_temp.shape)

    df_u_d_temp = df_u_d_temp.drop(['id_up', 'id_down'], axis=1)

    df_u_d_temp = df_u_d_temp.astype(float)
    pos_gc_content = (df_u_d_temp.sum(axis=0)) / len(df_u_d_temp)

    print(pos_gc_content.shape)

    return pos_gc_content, df_u_d_temp.columns

if __name__ == "__main__":
    start = time.perf_counter()

    # df_U_SEQ = pd.DataFrame()
    # df_D_SEQ = pd.DataFrame()
    # df_D_REGION = pd.DataFrame()
    len_df = 0
    df = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            # print(f'Return Value: {i.result()}')
            # temp_df_u_seq, temp_df_d_seq, temp_df_d_region = j.result()
            #
            # print('File Reading: Region dataframe size: {} and Up and Downstream Sequence dataframe size: {} and {}'.format(temp_df_d_region.shape, temp_df_u_seq.shape, temp_df_d_seq.shape))
            # if temp_df_d_region.shape[0] != temp_df_d_seq.shape[0]:
            #     print('Warning, Region and Seq does not match during File Reading for file {}'.format(j))
            #
            # elif temp_df_d_region.shape[0] != temp_df_u_seq.shape[0]:
            #     print('Warning, Region and Seq does not match during File Reading for file {}'.format(j))
            #
            # elif temp_df_d_seq.shape[0] != temp_df_u_seq.shape[0]:
            #     print('Warning, Upstream and Downstream Seq does not match during File Reading for file {}'.format(j))
            #
            # df_D_REGION = pd.concat([df_D_REGION, temp_df_d_region], axis=0)
            # df_D_SEQ = pd.concat([df_D_SEQ, temp_df_d_seq], axis=0)
            # df_U_SEQ = pd.concat([df_U_SEQ, temp_df_u_seq], axis=0)

            df_temp, temp_len_df = j.result()
            len_df = len_df + temp_len_df
            df = pd.concat([df, df_temp], axis=1)  ##Output is Series so we have to concatenate the Series horizontally
            print(df.shape)

    df = df.reset_index(drop=True)
    print(df.head(3))
    print(len_df)
    # print(df.sum(axis=1)/len_df)
    absolute_gc = df.sum(axis=1) / len_df
    absolute_gc.to_csv(Results_path+"\Files\/absolute_gc_content.csv")



    # df_U_SEQ_T = df_U_SEQ.replace(["A", "C", "G", "T", "N"], [1, 2, 3, 4, 0])
    # df_D_SEQ_T = df_D_SEQ.replace(["A", "C", "G", "T", "N"], [1, 2, 3, 4, 0])
    #
    # pos_gc_content, x_axis_col = absolute_gc_content(df_U_SEQ, df_D_SEQ)

    plt.figure(figsize=(20, 10))
    # sns.scatterplot(x=x_axis_col, y=pos_gc_content)
    sns.scatterplot(x=range(-1000, 1000), y=absolute_gc)

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content")
    plt.title("GC content per Base position")
    # plt.xticks(rotation = 90)
    plt.xticks(np.arange(-1000, 1000, 50), rotation=45)
    # plt.show()
    plt.savefig(Results_path+'\Charts\Chart0_absolute_gc_content.png')

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')

