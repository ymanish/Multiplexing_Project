import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt


Input_files_region = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\osativa\output_fasta"
Input_files_seq = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\osativa\input_fasta"

Results_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa"

Total_Input_Files = 2
Group = 'Plant' # 'Eukaryote'




def main(n):
    dol_COL = [str(i) for i in range(1000, 2000)]
    df_u_region = pd.DataFrame()
    df_d_region = pd.DataFrame()

    df_u_seq = pd.DataFrame()
    df_d_seq = pd.DataFrame()

    for (seq_record_1, seq_record_2) in zip(
            SeqIO.parse(Input_files_region+"\/region_group_"+str(n)+".fasta", "fasta"),
            SeqIO.parse(Input_files_seq+"\group_"+str(n)+".fasta", "fasta")):
        #     print(seq_record_1.id)
        #     print(seq_record_2.id)

        #     upstream_region = seq_record_1.seq[:1000]
        downstream_region = seq_record_1.seq[1000:2000]

        #     temp_u_region = pd.DataFrame(list(upstream_region))
        temp_d_region = pd.DataFrame(list(downstream_region), index=dol_COL)
        temp_d_region = temp_d_region.T
        temp_d_region['id'] = (seq_record_1.id).split('|')[0]

        df_d_region = pd.concat([df_d_region, temp_d_region], axis=0)

        # upstream_seq = seq_record_2.seq[:1000]
        downstream_seq = seq_record_2.seq[1000:2000]

        # temp_u_seq = pd.DataFrame(list(upstream_seq))
        temp_d_seq = pd.DataFrame(list(downstream_seq), index=dol_COL)

        # temp_u_seq = temp_u_seq.T
        temp_d_seq = temp_d_seq.T

        # temp_u_seq['id'] = (seq_record_2.id).split('|')[2]
        if Group == 'Plant':
            temp_d_seq['id'] = (seq_record_2.id).split('|')[1]
        else:
            temp_d_seq['id'] = (seq_record_2.id).split('|')[2]

        # df_u_seq = pd.concat([df_u_seq, temp_u_seq], axis=0)
        df_d_seq = pd.concat([df_d_seq, temp_d_seq], axis=0)

    return df_d_seq, df_d_region


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

    df_u_d_temp = df_region_matrix.merge(df_seq, on="id", how='left')
    df_u_d_temp = df_u_d_temp.dropna(axis=0)  ###Dropping any short sequences
    df_seq_region = pd.DataFrame()

    for pos in range(1000, 2000):
        df_seq_region.loc[:, str(pos)] = df_u_d_temp.loc[:, str(pos) + '_x'] * df_u_d_temp.loc[:, str(pos) + '_y']

    df_seq_region = df_seq_region.replace([1, 2, 3, 4], [0, 1, 1, 0])

    df_region_matrix = df_region_matrix.drop('id', axis=1)

    region_GC_per_pos = df_seq_region.sum(axis=0)  # Total GC bases in the region per position
    region_bases_per_pos = df_region_matrix.sum(axis=0)  # Total bases in the region per position
    density_per_pos = region_GC_per_pos / region_bases_per_pos

    region_GC = region_GC_per_pos.sum()  # Total GC bases in the region
    region_total = region_bases_per_pos.sum()  # Total bases in the region
    avg_region_GC = region_GC / region_total

    return density_per_pos, avg_region_GC



if __name__ == "__main__":
    start = time.perf_counter()

    # df_U_SEQ = pd.DataFrame()
    df_D_SEQ = pd.DataFrame()
    df_D_REGION = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
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


    # df_U_SEQ_T = df_U_SEQ.replace(["A", "C", "G", "T", "N"], [1, 2, 3, 4, 0])
    df_D_SEQ_T = df_D_SEQ.replace(["A", "C", "G", "T", "N"], [1, 2, 3, 4, 0])

    # df_D_REGION.drop('id', axis=1, inplace=True)
    # df_D_SEQ.drop('id', axis=1, inplace=True)

    exon_code = [0, 1, 0, 0, 0, 0]
    intron_code = [0, 0, 0, 1, 0, 0]
    UUTR_code = [1, 0, 0, 0, 0, 0]
    DUTR_code = [0, 0, 1, 0, 0, 0]

    df_D_REGION.reset_index(inplace=True, drop=True)
    df_D_SEQ_T.reset_index(inplace=True, drop=True)
    #
    # print(df_D_REGION.index)
    # print(df_D_SEQ_T.index)
    # import sys
    # sys.exit()

    density_exon_GC, avg_exon_GC = density_GC_Region_and_Avg(df_D_REGION, exon_code, df_D_SEQ_T)
    density_intron_GC, avg_intron_GC = density_GC_Region_and_Avg(df_D_REGION, intron_code, df_D_SEQ_T)
    density_UUTR_GC, avg_UUTR_GC = density_GC_Region_and_Avg(df_D_REGION, UUTR_code, df_D_SEQ_T)
    density_DUTR_GC, avg_DUTR_GC = density_GC_Region_and_Avg(df_D_REGION, DUTR_code, df_D_SEQ_T)

    densities = pd.DataFrame([density_exon_GC, density_intron_GC, density_UUTR_GC, density_DUTR_GC],
                             index=['exon', 'intron', 'UUTR', 'DUTR'])

    densities_T = densities.T
    densities_T['position'] = densities_T.index
    densities_T_melt = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)

    avg_GC_df = pd.DataFrame({'element': ['exon', 'intron', 'UUTR', 'DUTR'],
                              'Avg_GC': [avg_exon_GC, avg_intron_GC, avg_UUTR_GC, avg_DUTR_GC]
                              })

    avg_GC_df.to_csv(Results_path+'\Files\Avg_GC_per_element.csv')

    densities_T_melt['position'] = densities_T_melt['position'].astype(int)
    densities_T_melt['position'] = densities_T_melt['position'] - 1000
    densities_T_melt.to_csv(Results_path+"\Files\element_GC_density_perbasepos.csv")

    # x=range(1000)

    plt.figure(figsize=(10, 10))
    sns.barplot(x="element", y="Avg_GC", data=avg_GC_df)

    plt.xlabel("Region")
    plt.ylabel("%GC content")
    plt.title("Average GC content per Region")
    # plt.xticks(rotation = 90)
    # plt.xticks(np.arange(-1000, 1000, 50), rotation=45)
    # plt.show()
    plt.savefig(Results_path+'\Charts\Chart2_Avg_GC_Region.png')


    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt, x='position', y="density", hue='variable')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content")
    plt.title("%GC content distribution by Region per base position, with average across Region specific Transcripts")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(Results_path+'\Charts\Chart3_GC_Region_density.png')

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')