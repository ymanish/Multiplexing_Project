import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Initialise_Script import *
import warnings
warnings.filterwarnings('ignore')

### Data downloaded from the Biomart Query
chr_search_file = DATA_DIR + '/h_sapiens_chr_data.txt'
gtf_data = pd.read_csv(chr_search_file)

### Data Downloaded from the UCSC Genome browser
cpg_search_file = DATA_DIR + '/cpgIslandExt.txt'
cpg_data = pd.read_csv(cpg_search_file, sep='\t', header=None)

CPG_DENSITY_AND_GC_FILE =  r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_GC_and_Density.csv"
CpG_GC_CONTENT_DISTRIBUTION  = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_GC_content.png"
CpG_DENSITY = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_Density.png"
CpG_GC_CONTENT_DISTRIBUTION_REGION_WISE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_GC_content_region_wise.png"
SIGNAL_SIZE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\h_sapiens\CpG_GC_Signal_Size.png"


#### Calculate the singal size
def nucleosomal_rolling(Absolute_DI_DF, category, rolling_column):

    Absolute_DI_DF.sort_values(by=['position'], inplace=True)
    Absolute_DI_DF.reset_index(inplace=True, drop=True)

    Absolute_DI_DF['Rolling_signal'] = Absolute_DI_DF.groupby(category)[rolling_column].transform(
        lambda s: s.rolling(147).mean())
    Absolute_DI_DF['shifted_signal'] = Absolute_DI_DF.groupby(category)['Rolling_signal'].shift(-74)

    mean_data = Absolute_DI_DF.groupby(category).agg({'shifted_signal': 'mean'}).reset_index()
    mean_data.rename(columns={'shifted_signal': 'avg_signal'}, inplace=True)

    Absolute_DI_DF_ = Absolute_DI_DF.merge(mean_data, on=category, how='left')
    Absolute_DI_DF_['deviation_avg_signal'] = pd.Series.abs(Absolute_DI_DF_['shifted_signal'] - Absolute_DI_DF_['avg_signal'])
    mean_Absolute_DI_DF_ = Absolute_DI_DF_.groupby(category)['deviation_avg_signal'].sum().reset_index()

    return Absolute_DI_DF, mean_Absolute_DI_DF_

def get_cpg_seq(gene_start, gene_end, cpg_data_chr, length_seq, strand):
    gene_start_ = gene_start - 1000
    gene_end_ = gene_end + 1000

    cgi_islands_u = cpg_data_chr[(cpg_data_chr[3] > gene_start_) & (cpg_data_chr[2] <= gene_start_)]
    cgi_islands_u.loc[:, 2] = gene_start_
    cgi_islands_m = cpg_data_chr[(cpg_data_chr[3] > gene_start_) & (cpg_data_chr[2] >= gene_start_) &
                                 (cpg_data_chr[3] <= gene_end_) & (cpg_data_chr[2] < gene_end_)]

    cgi_islands_d = cpg_data_chr[(cpg_data_chr[3] > gene_end_) & (cpg_data_chr[2] <= gene_end_)]
    cgi_islands_d.loc[:, 3] = gene_end_


    cgi_islands = pd.concat([cgi_islands_u, cgi_islands_m, cgi_islands_d], axis=0)

    cpg_seq = 'S' * (length_seq)

    if strand == 1:
        for p in cgi_islands[[2, 3]].values.tolist():
            s = p[0] - gene_start_
            e = (p[1] - gene_start_) + 1
            #             print(s,e)
            cpg_seq = ('I' * (e - s)).join([cpg_seq[:s], cpg_seq[e:]])
    else:
        for p in cgi_islands[[2, 3]].values.tolist():
            s = gene_end_ - p[1]
            e = (gene_end_ - p[0]) + 1
            cpg_seq = ('I' * (e - s)).join([cpg_seq[:s], cpg_seq[e:]])

    return cpg_seq


def main(n):
    df_seq = pd.DataFrame()
    df_cpg_seq = pd.DataFrame()

    for seq_record in SeqIO.parse(SEQ_FILE_PATH+"/group_"+str(n)+".fasta", "fasta"):

        header = seq_record.id.split('|')
        STABLE_GENE_ID = header[0]
        STABLE_TRANSCRIPT_ID = header[2]
        cod_strand = int(header[11])
        #     print(STABLE_TRANSCRIPT_ID)
        #     print(cod_strand)

        EXON_STARTS = [int(i) for i in header[8].split(';')]  # array
        EXON_ENDS = [int(i) for i in header[9].split(';')]  # array

        #     TSS = int(header[10])

        EXON_RANK = [int(i) for i in header[14].split(';')]  # array
        EXON_STARTS_ = [i for _, i in sorted(zip(EXON_RANK, EXON_STARTS))]
        EXON_ENDS_ = [i for _, i in sorted(zip(EXON_RANK, EXON_ENDS))]

        try:
            CHR = 'chr' + str(gtf_data.loc[gtf_data['Transcript stable ID'] == STABLE_TRANSCRIPT_ID,
                                           'Chromosome/scaffold name'].values[0])
        except IndexError:
            print('The chromosome information not present for ', STABLE_TRANSCRIPT_ID)
            continue

        CPG_chr_data = cpg_data[cpg_data[1] == CHR]

        if cod_strand == 1:
            GENE_START = EXON_STARTS_[0]
            GENE_STOP = EXON_ENDS_[-1]

        else:
            GENE_START = EXON_STARTS_[-1]
            GENE_STOP = EXON_ENDS_[0]

        cpg_full_seq = get_cpg_seq(GENE_START, GENE_STOP, CPG_chr_data, len(seq_record.seq), cod_strand)

        sequ = list(seq_record.seq[:2000])
        temp_seq_df = pd.Series(sequ, index=COL)
        # temp_seq_df = temp_seq_df.T
        df_seq = pd.concat([df_seq, temp_seq_df], axis=1)

        cpg_seq = list(cpg_full_seq[:2000])
        temp_cpg_seq_df = pd.Series(cpg_seq, index=COL)
        # temp_cpg_seq_df = temp_cpg_seq_df.T
        df_cpg_seq = pd.concat([df_cpg_seq, temp_cpg_seq_df], axis=1)

    df_seq = df_seq.T
    df_cpg_seq = df_cpg_seq.T

    return df_cpg_seq, df_seq, df_seq.shape[0]

def replace_count(df_sam, to_replace, with_replace):
    df_sam = df_sam.replace(to_replace, with_replace)
    return df_sam



if __name__ == "__main__":
    start = time.perf_counter()
    len_df = 0
    DF_CpG = pd.DataFrame()
    DF_BASES = pd.DataFrame()
    DF_CpG_BASES_density = pd.DataFrame()


    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(main, n=i, ) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):

            temp_cpg_df, temp_bases_df, temp_len_df = j.result()
            len_df = len_df + temp_len_df
            DF_CpG = pd.concat([DF_CpG, temp_cpg_df], axis=0)
            DF_BASES = pd.concat([DF_BASES, temp_bases_df], axis=0)


    DF_CpG.reset_index(drop=True, inplace=True)
    DF_BASES.reset_index(drop=True, inplace=True)

    # print(DF_CpG)
    # print(DF_BASES)
    # print ('##################################################')
    DF_CpG_ = replace_count(DF_CpG, ["S", "I"], [0, 1])
    DF_BASES = replace_count(DF_BASES, ["A", "C", "G", "T", "N"], [0, 1, 1, 0, 0])
    DF_CpG_BASES = DF_CpG_.multiply(DF_BASES)

    ###Contribution from non-CpG islands

    DF_non_CpG = replace_count(DF_CpG, ["S", "I"], [1, 0])
    DF_non_CpG_BASES = DF_non_CpG.multiply(DF_BASES)



    total_cpg_count = DF_CpG_.sum(axis=0)
    total_non_cpg_count = DF_non_CpG.sum(axis=0)

    total_bases_count = DF_BASES.sum(axis=0)
    total_cpg_GC_count = DF_CpG_BASES.sum(axis=0)
    total_non_cpg_GC_count = DF_non_CpG_BASES.sum(axis=0)



    print(total_cpg_count)
    print(total_bases_count)
    print(total_cpg_GC_count)
    print(len_df)

    cpg_density = total_cpg_count / len_df
    cpg_density.name = 'CpG_density'

    non_cpg_density = total_non_cpg_count / len_df
    non_cpg_density.name = 'Non_CpG_density'

    base_content = total_bases_count/ len_df
    base_content.name = 'Base_content'

    base_content_CpG_region = total_cpg_GC_count/len_df
    base_content_CpG_region.name = 'CpG_Base_content'

    base_content_non_CpG_region = total_non_cpg_GC_count / len_df
    base_content_non_CpG_region.name = 'Non_CpG_Base_content'

    final_data = pd.concat([cpg_density, non_cpg_density, base_content, base_content_CpG_region, base_content_non_CpG_region],
                           axis=1)
    final_data['position'] = final_data.index

    final_data['Base_content_per_CpG_region'] = final_data['CpG_Base_content']/final_data['CpG_density']
    final_data['Base_content_per_non_CpG_region'] = final_data['Non_CpG_Base_content']/final_data['Non_CpG_density']


    final_data.to_csv(CPG_DENSITY_AND_GC_FILE)
    final_data_melted = pd.melt(final_data, id_vars=['position'], value_vars=['CpG_density',
                                                                              'Non_CpG_density',
                                                                              'Base_content',
                                                                              'CpG_Base_content',
                                                                              'Non_CpG_Base_content',
                                                                              'Base_content_per_CpG_region',
                                                                              'Base_content_per_non_CpG_region'
                                                                              ])
    final_data_melted['position'] = final_data_melted['position'].astype(int)

    plt.figure(figsize=(20, 10))
    # palette = sns.color_palette("Paired",)
    ax0 = sns.lineplot(x='position', y='value',
                       data=final_data_melted[final_data_melted['variable'].isin(['Base_content',
                                                                                  'CpG_Base_content',
                                                                                  'Non_CpG_Base_content'])],
                       hue='variable')

    ax0.set(xlabel='Position w.r.t TSS', ylabel='Fraction')
    ax0.tick_params(axis='x', labelrotation=0, labelsize=12)
    ax0.tick_params(axis='y', labelsize=12)
    ax0.set_title(" GC_content contribution from CpG islands per Base position", size=12)
    ax0.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=12)
    plt.savefig(CpG_GC_CONTENT_DISTRIBUTION)

    plt.figure(figsize=(20, 10))
    # palette = sns.color_palette("Paired",)
    ax01 = sns.lineplot(x='position', y='value',
                       data=final_data_melted[final_data_melted['variable'].isin(['Base_content_per_CpG_region',
                                                                                    'Base_content_per_non_CpG_region'])],
                       hue='variable')

    ax01.set(xlabel='Position w.r.t TSS', ylabel='Fraction')
    ax01.tick_params(axis='x', labelrotation=0, labelsize=12)
    ax01.tick_params(axis='y', labelsize=12)
    ax01.set_title("GC_content contribution from CpG islands per Base position, averaged across region specific", size=12)
    ax01.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=12)
    plt.savefig(CpG_GC_CONTENT_DISTRIBUTION_REGION_WISE)

    plt.figure(figsize=(20, 10))
    # palette = sns.color_palette("Paired",)
    ax1 = sns.lineplot(x='position', y='value',
                       data=final_data_melted[final_data_melted['variable'].isin(['CpG_density'])])

    ax1.set(xlabel='Position w.r.t TSS', ylabel='Fraction')  # ylabel='%GC content'
    ax1.tick_params(axis='x', labelrotation=0, labelsize=12)
    ax1.tick_params(axis='y', labelsize=12)
    ax1.set_title("CpG islands density per Base position", size=12)
    ax1.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=12)
    plt.savefig(CpG_DENSITY)

    #### Signal Size

    NUC_Data, signal_size_data = nucleosomal_rolling(final_data_melted, 'variable', 'value')

    ax11 = sns.barplot(x="variable", y="deviation_avg_signal", data=signal_size_data,
                      color=(0.4, 0.4, 0.4, 0.4), edgecolor='blue')
    ax11.set(xlabel='Data_type', ylabel="Signal Size")
    ax11.set_title("Signal Strength", size=12)
    ax11.tick_params(axis='both', labelrotation=90, labelsize=12)
    plt.savefig(SIGNAL_SIZE)

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')

