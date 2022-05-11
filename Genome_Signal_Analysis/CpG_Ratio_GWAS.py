import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Initialise_Script import *



def Read_FASTA(n):

    df_seq = pd.DataFrame()
    df_region = pd.DataFrame()

    for (seq_record_1, seq_record_2) in zip(
            SeqIO.parse(SEQ_FILE_PATH + "\group_" + str(n) + ".fasta", "fasta"),
            SeqIO.parse(REGION_FILE_PATH + "\/region_group_" + str(n) + ".fasta", "fasta")):
        header = seq_record_2.id.split('|')
        TRANSCRIPT_ID = header[0]

        sequ = str(seq_record_1.seq[:2000])
        sequ = sequ.replace('CG', 'XX')
        temp_seq = pd.Series(list(sequ), index=COL)
        temp_seq['id'] = TRANSCRIPT_ID
        df_seq = pd.concat([df_seq, temp_seq], axis=1)

        region_sequ = seq_record_2.seq[:2000]
        temp_seq_region = pd.Series(list(region_sequ), index=COL)
        temp_seq_region['id'] = TRANSCRIPT_ID
        df_region = pd.concat([df_region, temp_seq_region], axis=1)

    df_seq = df_seq.T
    df_region = df_region.T

    df_seq.reset_index(inplace=True, drop=True)
    df_region.reset_index(inplace=True, drop=True)

    return df_seq, df_region, len(df_seq)


def calculate_Cpg_ratios(df_seq, df_region):

    df_seq = df_seq.drop('id', axis=1)
    df_region = df_region.drop('id', axis=1)
    df_seq_ob = df_seq.replace(["A", "C", "G", "T", "N", 'X'], [0, 0, 0, 0, 0, 1])
    df_seq_exp = df_seq.replace(["A", "C", "G", "T", "N", 'X'], [0, 1, 1, 0, 0, 1])

    ob_data = observed_cpg(df_seq_ob)
    exp_data = expected_cpg(df_seq_exp)
    # print(df_seq)
    # print(df_region)
    to_rep = ["U", "E", "D", "I", "F", "f"]

    df_exon_temp = df_region.replace(to_rep, [0, 1, 0, 0, 0, 0])
    df_intron_temp = df_region.replace(to_rep, [0, 0, 0, 1, 0, 0])
    df_UUTR_temp = df_region.replace(to_rep, [1, 0, 0, 0, 0, 0])
    df_DUTR_temp = df_region.replace(to_rep, [0, 0, 1, 0, 0, 0])

    ob_exon_data = observed_cpg_region(df_seq_ob, df_exon_temp)
    exp_exon_data = expected_cpg_region(df_seq_exp, df_exon_temp)

    ob_intron_data = observed_cpg_region(df_seq_ob, df_intron_temp)
    exp_intron_data = expected_cpg_region(df_seq_exp, df_intron_temp)

    ob_UUTR_data = observed_cpg_region(df_seq_ob, df_UUTR_temp)
    exp_UUTR_data = expected_cpg_region(df_seq_exp, df_UUTR_temp)

    ob_DUTR_data = observed_cpg_region(df_seq_ob, df_DUTR_temp)
    exp_DUTR_data = expected_cpg_region(df_seq_exp, df_DUTR_temp)

    exon_D = df_exon_temp.sum(axis=0)
    intron_D = df_intron_temp.sum(axis=0)
    UUTR_D = df_UUTR_temp.sum(axis=0)
    DUTR_D = df_DUTR_temp.sum(axis=0)

    l = len(df_region)
    # print(ob_exon_data)

    return ob_data, exp_data, ob_exon_data, exp_exon_data, ob_intron_data, exp_intron_data, ob_UUTR_data, \
           exp_UUTR_data, ob_DUTR_data, exp_DUTR_data, exon_D, intron_D, UUTR_D, DUTR_D, l


def sum_columns(s):
    s[COL] = s[COL].rolling(window=147,  min_periods=1, center=True).sum()
    return s

def observed_cpg(df):
    # df = df.replace(["A", "C", "G", "T", "N", 'X'], [0, 0, 0, 0, 0, 1])
    df = df.apply(sum_columns, axis=1)
    df.iloc[:,0:73] = 0
    df.iloc[:,-73:] = 0
    df = df/2
    df = df/147
    return df


def expected_cpg(df):
    # df = df.replace(["A", "C", "G", "T", "N", 'X'], [0, 1, 1, 0, 0, 1])
    df = df.apply(sum_columns, axis=1)
    df = df/147
    df = df/2
    df = df**2
    return df


def ROLLING_SUM(d_df, len_data):
    d_df = d_df.sum(axis=1)
    d_df = d_df/len_data
    d_df = d_df.transform(lambda s: s.rolling(147).mean())
    d_df = d_df.shift(-73)
    return d_df


def observed_cpg_region(df, region_data):
    # df = df.replace(["A", "C", "G", "T", "N", 'X'], [0, 0, 0, 0, 0, 1])
    df_region_seq = df.multiply(region_data)
    df_region_seq = df_region_seq.apply(sum_columns, axis=1)
    df_region_seq.iloc[:,0:73] = 0
    df_region_seq.iloc[:,-73:] = 0
    df_region_seq = df_region_seq/2
    df_region_seq = df_region_seq/147

    return df_region_seq


def expected_cpg_region(df, region_data):
    # df = df.replace(["A", "C", "G", "T", "N", 'X'], [0, 1, 1, 0, 0, 1])
    df_region_seq = df.multiply(region_data)
    df_region_seq = df_region_seq.apply(sum_columns, axis=1)
    df_region_seq = df_region_seq/147
    df_region_seq = df_region_seq/2
    df_region_seq = df_region_seq**2
#     print(df_region_seq.iloc[995:,1000:-73])
#     df_region_seq = df_region_seq/147
    return df_region_seq


def obs_region_cpg_content(obs_exon, obs_intron, obs_UUTR, obs_DUTR):
    exon_cpg_content = obs_exon.mean(axis=0)
    intron_cpg_content = obs_intron.mean(axis=0)
    UUTR_cpg_content = obs_UUTR.mean(axis=0)
    DUTR_cpg_content = obs_DUTR.mean(axis=0)


    obs_cpg_region_densities = pd.DataFrame([exon_cpg_content, intron_cpg_content, UUTR_cpg_content, DUTR_cpg_content],
                            index=['exon', 'intron', 'UUTR', 'DUTR'])

    obs_cpg_region_densities_T = obs_cpg_region_densities.T
    obs_cpg_region_densities_T['position'] = obs_cpg_region_densities_T.index
    obs_cpg_region_densities_T_melt = pd.melt(obs_cpg_region_densities_T, id_vars=['position'])
    obs_cpg_region_densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)

    obs_cpg_region_densities_T_melt['position'] = obs_cpg_region_densities_T_melt['position'].astype(int)
    # obs_cpg_region_densities_T_melt.to_csv(CPG_PER_REGION_PER_BP_FILE)

    obs_cpg_region_densities_T_melt['position'] = obs_cpg_region_densities_T_melt['position'] - 1000
    return obs_cpg_region_densities_T_melt

if __name__ == "__main__":

    start_time = time.perf_counter()

    len_df = 0
    DATA_DF = pd.DataFrame()
    DATA_REGION_DF = pd.DataFrame()

    print('Reading the Sequence Files................')

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files + 1)
        pool = [executor.submit(Read_FASTA, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            df_temp, df_region_temp, temp_len_df = j.result()
            len_df = len_df + temp_len_df
            DATA_DF = pd.concat([DATA_DF, df_temp], axis=0)
            DATA_REGION_DF = pd.concat([DATA_REGION_DF, df_region_temp], axis=0)

    DATA_DF.reset_index(inplace=True, drop=True)
    DATA_REGION_DF.reset_index(inplace=True, drop=True)

    print(DATA_DF.shape)
    print(DATA_REGION_DF.shape)
    print (len_df)

    ocpg = pd.DataFrame()
    ecpg = pd.DataFrame()

    obs_exon_cpg = pd.DataFrame()
    exp_exon_cpg = pd.DataFrame()
    obs_intron_cpg = pd.DataFrame()
    exp_intron_cpg = pd.DataFrame()
    obs_UUTR_cpg = pd.DataFrame()
    exp_UUTR_cpg = pd.DataFrame()
    obs_DUTR_cpg = pd.DataFrame()
    exp_DUTR_cpg = pd.DataFrame()
    trans_df_exon = pd.DataFrame()
    trans_df_intron = pd.DataFrame()
    trans_df_UUTR = pd.DataFrame()
    trans_df_DUTR = pd.DataFrame()

    leng = 0
    with concurrent.futures.ProcessPoolExecutor() as executor:
        start = []
        end = []
        for i in range(0, len_df, 1000):
            start.append(i)
        end = start[1:]
        end.append(len_df)
        pool = [executor.submit(calculate_Cpg_ratios, DATA_DF.iloc[start[j]:end[j],],
                                DATA_REGION_DF.iloc[start[j]:end[j],]) for j in range(len(start))]
        for k in concurrent.futures.as_completed(pool):
            ob_temp, exp_temp, ob_exon_temp, exp_exon_temp, ob_intron_temp, exp_intron_temp, ob_UUTR_temp, \
            exp_UUTR_temp, ob_DUTR_temp, exp_DUTR_temp, exon_D_temp, intron_D_temp, UUTR_D_temp, DUTR_D_temp, len_temp = k.result()

            leng = leng + len_temp

            ocpg = pd.concat([ocpg, ob_temp], axis=0)
            ecpg = pd.concat([ecpg, exp_temp], axis=0)

            obs_exon_cpg = pd.concat([obs_exon_cpg, ob_exon_temp], axis=0)
            exp_exon_cpg = pd.concat([exp_exon_cpg, exp_exon_temp], axis=0)

            obs_intron_cpg = pd.concat([obs_intron_cpg, ob_intron_temp], axis=0)
            exp_intron_cpg = pd.concat([exp_intron_cpg, exp_intron_temp], axis=0)

            obs_UUTR_cpg = pd.concat([obs_UUTR_cpg, ob_UUTR_temp], axis=0)
            exp_UUTR_cpg = pd.concat([exp_UUTR_cpg, exp_UUTR_temp], axis=0)

            obs_DUTR_cpg = pd.concat([obs_DUTR_cpg, ob_DUTR_temp], axis=0)
            exp_DUTR_cpg = pd.concat([exp_DUTR_cpg, exp_DUTR_temp], axis=0)

            trans_df_exon = pd.concat([trans_df_exon, exon_D_temp],
                                      axis=1)
            trans_df_intron = pd.concat([trans_df_intron, intron_D_temp],
                                        axis=1)
            trans_df_UUTR = pd.concat([trans_df_UUTR, UUTR_D_temp],
                                      axis=1)
            trans_df_DUTR = pd.concat([trans_df_DUTR, DUTR_D_temp],
                                      axis=1)

    density_exon = ROLLING_SUM(trans_df_exon, leng)
    density_intron = ROLLING_SUM(trans_df_intron, leng)
    density_UUTR = ROLLING_SUM(trans_df_UUTR, leng)
    density_DUTR = ROLLING_SUM(trans_df_DUTR, leng)

    densities = pd.DataFrame([density_exon, density_intron, density_UUTR, density_DUTR],
                             index=['exon', 'intron', 'UUTR', 'DUTR'])

    densities_T = densities.T
    densities_T['position'] = densities_T.index
    densities_T_melt = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)

    densities_T_melt['position'] = densities_T_melt['position'].astype(int)
    densities_T_melt.to_csv(REGION_DENSITY_FILE)

    densities_T_melt['position'] = densities_T_melt['position'] - 1000
    #
    # plt.figure(figsize=(20, 10))
    # sns.lineplot(data=densities_T_melt, x='position', y='density', hue="variable")
    # plt.xlabel("Position w.r.t TSS")
    # plt.ylabel("%Occurrence")
    # plt.title("Element Density per Base position")
    # plt.xticks(rotation=0)
    # plt.show()
    # # plt.savefig(REGION_DENSITY_CHART)
    #
    #
    # print ('Calculating the Observed Cpg Region Ratio....')

    ocpg.reset_index(inplace=True, drop=True)
    ecpg.reset_index(inplace=True, drop=True)

    obs_exon_cpg.reset_index(inplace=True, drop=True)
    exp_exon_cpg.reset_index(inplace=True, drop=True)
    obs_intron_cpg.reset_index(inplace=True, drop=True)
    exp_intron_cpg.reset_index(inplace=True, drop=True)
    obs_UUTR_cpg.reset_index(inplace=True, drop=True)
    exp_UUTR_cpg.reset_index(inplace=True, drop=True)
    obs_DUTR_cpg.reset_index(inplace=True, drop=True)
    exp_DUTR_cpg.reset_index(inplace=True, drop=True)

    # eo_ratio = ocpg / ecpg
    # eo_ratio_average = eo_ratio.mean(axis=0)
    #
    # eo_exon_ratio = obs_exon_cpg / exp_exon_cpg
    # eo_exon_ratio.fillna(0, inplace=True)
    #
    # eo_intron_ratio = obs_intron_cpg / exp_intron_cpg
    # eo_intron_ratio.fillna(0, inplace=True)
    #
    # eo_UUTR_ratio = obs_UUTR_cpg / exp_UUTR_cpg
    # eo_UUTR_ratio.fillna(0, inplace=True)
    #
    # eo_DUTR_ratio = obs_DUTR_cpg / exp_DUTR_cpg
    # eo_DUTR_ratio.fillna(0, inplace=True)
    # print(eo_exon_ratio.iloc[995:,1000:-73])
    # eo_exon_ratio_avg = eo_exon_ratio.mean(axis=0)
    # eo_intron_ratio_avg = eo_intron_ratio.mean(axis=0)
    # eo_UUTR_ratio_avg = eo_UUTR_ratio.mean(axis=0)
    # eo_DUTR_ratio_avg = eo_DUTR_ratio.mean(axis=0)
    # print(eo_exon_ratio_avg)
    #
    # eo_ratio_average.to_csv(ABSOLUTE_CPG_FILE)
    # plt.figure(figsize=(20, 10))
    # sns.lineplot(x=range(-1000, 1000), y=eo_ratio_average)
    #
    # plt.xlabel("Position w.r.t TSS")
    # plt.ylabel("Normalised CpG")
    # plt.title("Normalised CpG Ratio per Base position")
    # # plt.xticks(rotation = 90)
    # plt.xticks(np.arange(-1000, 1000, 50), rotation=45)
    # plt.show()
    # # plt.savefig(ABSOLUTE_CPG_CHART)
    #
    #
    # eo_ratio_region_densities = pd.DataFrame([eo_exon_ratio_avg, eo_intron_ratio_avg, eo_UUTR_ratio_avg, eo_DUTR_ratio_avg],
    #                         index=['exon', 'intron', 'UUTR', 'DUTR'])
    #
    # eo_ratio_region_densities_T = eo_ratio_region_densities.T
    # eo_ratio_region_densities_T['position'] = eo_ratio_region_densities_T.index
    # eo_ratio_region_densities_T_melt = pd.melt(eo_ratio_region_densities_T, id_vars=['position'])
    # eo_ratio_region_densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)
    #
    # eo_ratio_region_densities_T_melt['position'] = eo_ratio_region_densities_T_melt['position'].astype(int)
    # eo_ratio_region_densities_T_melt.to_csv(CPG_PER_REGION_PER_BP_FILE)
    #
    # eo_ratio_region_densities_T_melt['position'] = eo_ratio_region_densities_T_melt['position'] - 1000
    #
    # plt.figure(figsize=(20, 10))
    # sns.lineplot(data=eo_ratio_region_densities_T_melt, x='position', y="density", hue='variable')
    #
    # plt.xlabel("Position w.r.t TSS")
    # plt.ylabel("Normalised CpG")
    # plt.title("CpG content distribution by Region per base position")
    # plt.xticks(rotation=0)
    # plt.show()
    # # plt.savefig(CPG_PER_REGION_PER_BP_CHART)
    #
    # end_time = time.perf_counter()
    # print(f'Finished in {round(end_time - start_time, 2)} second(s)')




    print('Calculate the Observed CpG fraction across all the transcripts per base position')

    # print(ocpg)
    # print(obs_exon_cpg)

    print(obs_exon_cpg.mean(axis=0))

    obs_cpg_content_df = ocpg.mean(axis=0)

    plt.figure(figsize=(20, 10))
    sns.lineplot(x=range(-1000, 1000), y=obs_cpg_content_df)

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("CpG fraction")
    plt.title("Observed CpG fraction per Base position")
    # plt.xticks(rotation = 90)
    plt.xticks(np.arange(-1000, 1000, 50), rotation=45)
    plt.show()
    # plt.savefig()


    obs_cpg_per_region_per_base = obs_region_cpg_content(obs_exon_cpg, obs_intron_cpg, obs_UUTR_cpg, obs_DUTR_cpg)

    plt.figure(figsize=(20, 10))
    sns.lineplot(data=obs_cpg_per_region_per_base, x='position', y="density", hue='variable')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("CpG content")
    plt.title("Observed CpG content distribution by Region per base position")
    plt.xticks(rotation=0)
    plt.show()
    # plt.savefig()


    print(obs_cpg_per_region_per_base)
    print(densities_T_melt)
    obs_cpg_per_region_per_base.rename(columns={'density':'CpG_fraction'}, inplace=True)
    obs_cpg_per_base_by_region = obs_cpg_per_region_per_base.merge(densities_T_melt,
                                      on=['position', 'variable'],
                                      how='left')
    obs_cpg_per_base_by_region['Cpg_fraction_new'] = obs_cpg_per_base_by_region['CpG_fraction']/obs_cpg_per_base_by_region['density']
    plt.figure(figsize=(20, 10))
    sns.lineplot(data=obs_cpg_per_base_by_region, x='position', y="Cpg_fraction_new", hue='variable')

    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("CpG content")
    plt.title("Observed CpG content distribution by Region per base position divided by Region")
    plt.xticks(rotation=0)
    plt.show()

    avg_cpg_content = obs_cpg_per_base_by_region.groupby('variable')['Cpg_fraction_new'].mean()
    avg_cpg_content = avg_cpg_content.to_frame()
    avg_cpg_content.reset_index(inplace=True)

    plt.figure(figsize=(10, 10))
    sns.barplot(x="variable", y="Cpg_fraction_new", data=avg_cpg_content)

    plt.xlabel("Region")
    plt.ylabel("CpG fraction")
    plt.title("Average CpG fraction per Region")
    plt.show()
    # plt.savefig(AVG_GC_PER_REGION_CHART)