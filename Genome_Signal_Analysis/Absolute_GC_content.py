import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt


def main(n):

    df_u_region = pd.DataFrame()
    df_d_region = pd.DataFrame()

    df_u_seq = pd.DataFrame()
    df_d_seq = pd.DataFrame()

    for (seq_record_1, seq_record_2) in zip(
            SeqIO.parse(r"C:\Users\maya620d\Documents\Notebooks\Output_files\region_human_"+str(n)+".fasta", "fasta"),
            SeqIO.parse(r"C:\Users\maya620d\Documents\Notebooks\Input_files\chunk"+str(n)+".fasta", "fasta")):
        #     print(seq_record_1.id)
        #     print(seq_record_2.id)

        #     upstream_region = seq_record_1.seq[:1000]
        downstream_region = seq_record_1.seq[1000:2000]

        #     temp_u_region = pd.DataFrame(list(upstream_region))
        temp_d_region = pd.DataFrame(list(downstream_region))
        temp_d_region = temp_d_region.T
        #     df_u_region=pd.concat([df_u_region, temp_u_region], axis =1)

        temp_d_region['id'] = (seq_record_1.id).split('|')[0]

        df_d_region = pd.concat([df_d_region, temp_d_region], axis=0)

        upstream_seq = seq_record_2.seq[:1000]
        downstream_seq = seq_record_2.seq[1000:2000]

        temp_u_seq = pd.DataFrame(list(upstream_seq))
        temp_d_seq = pd.DataFrame(list(downstream_seq))

        temp_u_seq = temp_u_seq.T
        temp_d_seq = temp_d_seq.T

        temp_u_seq['id'] = (seq_record_2.id).split('|')[2]
        temp_d_seq['id'] = (seq_record_2.id).split('|')[2]

        df_u_seq = pd.concat([df_u_seq, temp_u_seq], axis=0)
        df_d_seq = pd.concat([df_d_seq, temp_d_seq], axis=0)

    return df_u_seq, df_d_seq, df_d_region


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

    df_U_SEQ = pd.DataFrame()
    df_D_SEQ = pd.DataFrame()
    df_D_REGION = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(0, 97000, 1000)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            # print(f'Return Value: {i.result()}')
            temp_df_u_seq, temp_df_d_seq, temp_df_d_region = j.result()

            print('File Reading: Region dataframe size: {} and Up and Downstream Sequence dataframe size: {} and {}'.format(temp_df_d_region.shape, temp_df_u_seq.shape, temp_df_d_seq.shape))
            if temp_df_d_region.shape[0] != temp_df_d_seq.shape[0]:
                print('Warning, Region and Seq does not match during File Reading for file {}'.format(j))

            elif temp_df_d_region.shape[0] != temp_df_u_seq.shape[0]:
                print('Warning, Region and Seq does not match during File Reading for file {}'.format(j))

            elif temp_df_d_seq.shape[0] != temp_df_u_seq.shape[0]:
                print('Warning, Upstream and Downstream Seq does not match during File Reading for file {}'.format(j))

            df_D_REGION = pd.concat([df_D_REGION, temp_df_d_region], axis=0)
            df_D_SEQ = pd.concat([df_D_SEQ, temp_df_d_seq], axis=0)
            df_U_SEQ = pd.concat([df_U_SEQ, temp_df_u_seq], axis=0)

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')

    df_U_SEQ_T = df_U_SEQ.replace(["A", "C", "G", "T", "N"], [1, 2, 3, 4, 0])
    df_D_SEQ_T = df_D_SEQ.replace(["A", "C", "G", "T", "N"], [1, 2, 3, 4, 0])

    pos_gc_content, x_axis_col = absolute_gc_content(df_U_SEQ, df_D_SEQ)

    plt.figure(figsize=(20, 10))
    # sns.scatterplot(x=x_axis_col, y=pos_gc_content)
    sns.scatterplot(x=range(-1000, 1000), y=pos_gc_content)

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content")
    plt.title("")
    # plt.xticks(rotation = 90)
    plt.xticks(np.arange(-1000, 1000, 50), rotation=45)
    # plt.show()
    plt.savefig('absolute_gc_content.png')



