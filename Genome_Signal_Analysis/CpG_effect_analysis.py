import pandas as pd
import numpy as np
from Bio import SeqIO
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Initialise_Script import *
import re
import random
from Dinucleotide_Analysis import di_nucleotide_count


CpG_Eff_DI_CONTENT_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_effect_di_nucleotide_content.csv"
CpG_Eff_DI_CONTENT_CHART = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_effect_di_nucleotide_content.png"
CpG_Eff_SMOTHENED_DI_CONTENT_CHART = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\Nucleosomal_CpG_effect_di_nucleotide_content.png"

CpG_WITHOUT_Eff_DI_CONTENT_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_without_effect_di_nucleotide_content.csv"
CpG_WITHOUT_Eff_DI_CONTENT_CHART = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_without_effect_di_nucleotide_content.png"
CpG_WITHOUT_Eff_SMOTHENED_DI_CONTENT_CHART = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\Nucleosomal_CpG_without_effect_di_nucleotide_content.png"

def main(n):
    DI_POS_DICT = {}

    for i in range(2000):
        DI_POS_DICT[str(i)] = {'CG': 0, 'CC': 0, 'CT': 0, 'CA': 0, 'GG': 0, 'GC': 0, 'GT': 0,
                               'GA': 0, 'TG': 0, 'TC': 0, 'TT': 0, 'TA': 0, 'AG': 0, 'AC': 0, 'AT': 0, 'AA': 0}
    records = 0
    for seq_record in SeqIO.parse(SEQ_FILE_PATH+"/group_"+str(n)+".fasta", "fasta"):
        header = seq_record.id.split('|')
        if GROUP == 'Plant':
            print(header[1])
            TRANSCRIPT_ID = header[1]
        else:
            print(header[2])
            TRANSCRIPT_ID = header[2]

        SEQ = str(seq_record.seq[:2000])

        original_seq_array = []

        for seq in re.split('(CG)', SEQ):
            if seq not in ['CG', '']:

                previous_di_NN = 'CG'
                new_seq = ''
                for i in range(0, len(seq), 2):
                    if i != len(seq) - 1:
                        # print(i)
                        if previous_di_NN[1] == 'C':
                            #                         print('current_di_NN has 11 possibilities')
                            current_di_NN = random.choice(
                                ['CT', 'CA', 'TG', 'TC', 'TA', 'AG', 'AC', 'AT', 'CC', 'TT', 'AA'])
                            new_seq = new_seq + current_di_NN
                        else:
                            #                     print('current_di_NN has 15 possibilities')
                            current_di_NN = random.choice(
                                ['CT', 'CA', 'TG', 'TC', 'TA', 'AG', 'AC', 'AT', 'CC', 'TT', 'AA', 'GG', 'GA', 'GT',
                                 'GC'])
                            new_seq = new_seq + current_di_NN

                    else:
                        if previous_di_NN[1] == 'C':
                            # print('A, T, or C')
                            current_di_NN = random.choice(['C', 'A', 'T'])
                            new_seq = new_seq + current_di_NN

                        else:
                            # print('A, T, or C')
                            current_di_NN = random.choice(['C', 'A', 'T', 'G'])
                            new_seq = new_seq + current_di_NN

                    previous_di_NN = current_di_NN
                # print(new_seq)
                original_seq_array.append(new_seq)


            else:
                original_seq_array.append(seq)

        random_seq = ''.join(original_seq_array)

        assert(len(re.split('(CG)', SEQ)) == len(re.split('(CG)', random_seq)))

        DI_POS_DICT = di_nucleotide_count(random_seq, DI_POS_DICT)
        records = records+1

    return DI_POS_DICT, records


def main_without_CpG(n):
    DI_POS_DICT = {}

    for i in range(2000):
        DI_POS_DICT[str(i)] = {'CG': 0, 'CC': 0, 'CT': 0, 'CA': 0, 'GG': 0, 'GC': 0, 'GT': 0,
                               'GA': 0, 'TG': 0, 'TC': 0, 'TT': 0, 'TA': 0, 'AG': 0, 'AC': 0, 'AT': 0, 'AA': 0}
    records = 0
    for seq_record in SeqIO.parse(SEQ_FILE_PATH + "/group_" + str(n) + ".fasta", "fasta"):
        header = seq_record.id.split('|')
        if GROUP == 'Plant':
            print(header[1])
            TRANSCRIPT_ID = header[1]
        else:
            print(header[2])
            TRANSCRIPT_ID = header[2]

        SEQ = str(seq_record.seq[:2000])
        SEQ_ = SEQ.replace('CG', 'NN')

        DI_POS_DICT = di_nucleotide_count(SEQ_, DI_POS_DICT)

        records = records + 1

    return DI_POS_DICT, records


if __name__ == "__main__":

    WITH_CpG = True

    for j in range(0, 5):

        CpG_Eff_DI_CONTENT_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_effect_di_nucleotide_content"+str(j)+".csv"
        start = time.perf_counter()
        len_df = 0
        df = pd.DataFrame()
        Total_records = 0
        
        with concurrent.futures.ProcessPoolExecutor() as executor:
            iter_seq = range(1, Total_Input_Files+1)
            if WITH_CpG:
                pool = [executor.submit(main, n=i) for i in iter_seq]
            else:
                pool= [executor.submit(main_without_CpG, n=i) for i in iter_seq]

            for j in concurrent.futures.as_completed(pool):

                # df_temp, temp_len_df = j.result()
                temp_di_count_matrix, record_count = j.result()
                df_temp = pd.DataFrame(temp_di_count_matrix)
                df = pd.concat([df, df_temp], axis=0)
                Total_records = Total_records + record_count
                # print(df.shape)
        df['DI_NUCL'] = df.index
        df = df.reset_index(drop=True)
        di_nuc_content_df = df.groupby('DI_NUCL').sum().reset_index()
        di_nuc_content_df[COL] = di_nuc_content_df[COL] / Total_records
        # print (di_nuc_content_df)
        # print (Total_records)

        di_nuc_content_df_T = pd.melt(di_nuc_content_df,
                                      id_vars=['DI_NUCL'],
                                      value_vars=COL)
        di_nuc_content_df_T.rename(columns={'variable': 'position', 'value': 'di_nuc_fraction'}, inplace=True)
        di_nuc_content_df_T['position'] = di_nuc_content_df_T['position'].astype(int)
        di_nuc_content_df_T['position'] = di_nuc_content_df_T['position'] - 1000

        print(di_nuc_content_df_T)

        if WITH_CpG:

            di_nuc_content_df_T.to_csv(CpG_Eff_DI_CONTENT_FILE)
        else:
            di_nuc_content_df_T.to_csv(CpG_WITHOUT_Eff_DI_CONTENT_FILE)

        plt.figure(figsize=(20, 10))
        palette = sns.color_palette("Paired", 16)
        ax0 = sns.lineplot(x='position', y='di_nuc_fraction',
                           data=di_nuc_content_df_T,
                           hue='DI_NUCL', palette=palette)  # y='shifted_GC' , , y='GC_content'

        ax0.set(xlabel='Position w.r.t TSS', ylabel='Fraction')  # ylabel='%GC content'
        ax0.tick_params(axis='x', labelrotation=0, labelsize=12)
        ax0.tick_params(axis='y', labelsize=12)
        ax0.set_title("DI nucleotide content per Base position", size=12)
        ax0.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=12)

        if WITH_CpG:
            plt.savefig(CpG_Eff_DI_CONTENT_CHART)
        else:
            plt.savefig(CpG_WITHOUT_Eff_DI_CONTENT_CHART)

        #####SMOTHENED DI_NUCLEOTIDE CHART

        di_nuc_content_df_T.sort_values(by=['position'], inplace=True)
        di_nuc_content_df_T.reset_index(inplace=True, drop=True)

        di_nuc_content_df_T['Rolling_DI'] = di_nuc_content_df_T.groupby('DI_NUCL')['di_nuc_fraction'].transform(
            lambda s: s.rolling(147).mean())
        di_nuc_content_df_T['shifted_DI'] = di_nuc_content_df_T.groupby('DI_NUCL')['Rolling_DI'].shift(-74)

        mean_data = di_nuc_content_df_T.groupby('DI_NUCL').agg({'shifted_DI': 'mean'}).reset_index()
        mean_data.rename(columns={'shifted_DI': 'avg_DI'}, inplace=True)

        di_nuc_content_df_T_ = di_nuc_content_df_T.merge(mean_data, on='DI_NUCL', how='left')
        di_nuc_content_df_T_['deviation_avg_DI'] = pd.Series.abs(
            di_nuc_content_df_T_['shifted_DI'] - di_nuc_content_df_T_['avg_DI'])
        mean_di_nuc_content_df_T_ = di_nuc_content_df_T_.groupby('DI_NUCL')['deviation_avg_DI'].sum().reset_index()

        plt.figure(figsize=(20, 10))
        palette = sns.color_palette("Paired", 16)
        ax0 = sns.lineplot(x='position', y='shifted_DI',
                           data=di_nuc_content_df_T,
                           hue='DI_NUCL', palette=palette)  # y='shifted_GC' , , y='GC_content'

        ax0.set(xlabel='Position w.r.t TSS', ylabel='Fraction')  # ylabel='%GC content'
        ax0.tick_params(axis='x', labelrotation=0, labelsize=12)
        ax0.tick_params(axis='y', labelsize=12)
        ax0.set_title("DI nucleotide content per Base position", size=12)
        ax0.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=12)

        if WITH_CpG:
            plt.savefig(CpG_Eff_SMOTHENED_DI_CONTENT_CHART)
        else:
            plt.savefig(CpG_WITHOUT_Eff_SMOTHENED_DI_CONTENT_CHART)

        end = time.perf_counter()
        print(f'Finished in {round(end - start, 2)} second(s)')

