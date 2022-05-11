import copy

import pandas as pd
from Bio import SeqIO
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Initialise_GC_profiling import *


def main(n):
    df_region = pd.DataFrame()
    df_seq = pd.DataFrame()

    for (seq_record_1, seq_record_2) in zip(
            SeqIO.parse(REGION_FILE_PATH+"\/region_group_"+str(n)+".fasta", "fasta"),
            SeqIO.parse(SEQ_FILE_PATH+"\group_"+str(n)+".fasta", "fasta")):

        region = seq_record_1.seq[:2000]
        temp_region = pd.DataFrame(list(region), index=COL)
        temp_region = temp_region.T
        temp_region['id'] = (seq_record_1.id).split('|')[0]

        df_region = pd.concat([df_region, temp_region], axis=0)

        seq = seq_record_2.seq[:2000]
        temp_seq = pd.DataFrame(list(seq), index=COL)
        temp_seq = temp_seq.T
        if GROUP == 'Plant':
            temp_seq['id'] = (seq_record_2.id).split('|')[1]
        else:
            temp_seq['id'] = (seq_record_2.id).split('|')[2]

        df_seq = pd.concat([df_seq, temp_seq], axis=0)

    df_region = df_region.reset_index(drop=True)
    df_seq = df_seq.reset_index(drop=True)

    df_seq = df_seq.drop('id', axis=1)
    df_region = df_region.drop('id', axis=1)

    df_seq_exon, exon_bases = Region_Seq_multipliction(df_region, [0, 1, 0, 0, 0, 0], df_seq)
    df_seq_intron, intron_bases = Region_Seq_multipliction(df_region, [0, 0, 0, 1, 0, 0], df_seq)
    df_seq_UUTR, UUTR_bases = Region_Seq_multipliction(df_region, [1, 0, 0, 0, 0, 0], df_seq)
    df_seq_DUTR, DUTR_bases = Region_Seq_multipliction(df_region, [0, 0, 1, 0, 0, 0], df_seq)
    df_seq_UF, UF_bases = Region_Seq_multipliction(df_region, [0, 0, 0, 0, 1, 0], df_seq)
    df_seq_DF, DF_bases = Region_Seq_multipliction(df_region, [0, 0, 0, 0, 0, 1], df_seq)

    return df_seq, df_region, df_seq_exon, df_seq_intron, df_seq_UUTR, df_seq_DUTR, df_seq_UF, df_seq_DF,\
           exon_bases, intron_bases, UUTR_bases, DUTR_bases, UF_bases, DF_bases, len(df_seq)


def Region_Seq_multipliction(DF_region, region_code, DF_seq):
    df_region_matrix = DF_region.replace(["U", "E", "D", "I", "F", "f"], region_code)

    print('Region dataframe size: {} and Sequence dataframe size: {}'.format(DF_region.shape, DF_seq.shape))
    if DF_region.shape[0] != DF_seq.shape[0]:
        print('Warning, Region and Seq does not match')
    DF_seq = DF_seq.replace(["A", "C", "G", "T", "N", "S"], [0, 1, 1, 0, 0, 1])
    DF_seq = DF_seq.replace('[A-Z]', 0, regex=True)
    df_seq_region = df_region_matrix.multiply(DF_seq)
    total_region_bases = df_region_matrix.sum(axis=0).sum()
    return df_seq_region, total_region_bases


def density_GC_Region_and_Avg(df_region, region_code, df_seq):
    df_region_matrix = df_region.replace(["U", "E", "D", "I", "F", "f"], region_code)

    print('Region dataframe size: {} and Sequence dataframe size: {}'.format(df_region.shape, df_seq.shape))
    if df_region.shape[0] != df_seq.shape[0]:
        print('Warning, Region and Seq does not match')
    df_seq = df_seq.replace(["A", "C", "G", "T", "N", "S"], [0, 1, 1, 0, 0, 1])
    df_seq_region = df_region_matrix.multiply(df_seq)

    region_GC_per_pos = df_seq_region.sum(axis=0)  # Total GC bases in the region per position
    region_bases_per_pos = df_region_matrix.sum(axis=0)  # Total bases in the region per position
    density_per_pos = region_GC_per_pos / region_bases_per_pos

    region_GC = region_GC_per_pos.sum()  # Total GC bases in the region
    region_total = region_bases_per_pos.sum()  # Total bases in the region
    avg_region_GC = region_GC / region_total

    return density_per_pos, avg_region_GC

def multiply_operation(el, density, tab):
    return density * tab.loc[tab['element'] == el, 'Avg_GC'].values[0]


if __name__ == "__main__":
    start = time.perf_counter()

    df_SEQ = pd.DataFrame()
    df_REGION = pd.DataFrame()
    df_EXON = pd.DataFrame()
    df_INTRON = pd.DataFrame()
    df_UUTR = pd.DataFrame()
    df_DUTR = pd.DataFrame()
    df_UF = pd.DataFrame()
    df_DF = pd.DataFrame()

    Total_Exon_bases = 0
    Total_Intron_bases = 0
    Total_UUTR_bases = 0
    Total_DUTR_bases = 0
    Total_UF_bases = 0
    Total_DF_bases = 0
    total_records = 0
    print('Reading the Sequence and Region files........')
    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            # print(f'Return Value: {i.result()}')
            temp_df_seq, temp_df_region, temp_df_exon, temp_df_intron, temp_df_UUTR, temp_df_DUTR, temp_df_UF, temp_df_DF,\
            Temp_Exon_bases, Temp_Intron_bases, Temp_UUTR_bases, Temp_DUTR_bases, Temp_UF_bases, Temp_DF_bases, rec = j.result()
            print('File Reading: Region dataframe size: {} and Sequence dataframe size: {}'.format(temp_df_region.shape, temp_df_seq.shape))
            if temp_df_region.shape[0] != temp_df_seq.shape[0]:
                print('Warning, Region and Seq does not match during File Reading for file {}'.format(j))

            df_REGION = pd.concat([df_REGION, temp_df_region], axis=0)
            df_SEQ = pd.concat([df_SEQ, temp_df_seq], axis=0)

            df_EXON = pd.concat([df_EXON, temp_df_exon], axis=0)
            df_INTRON = pd.concat([df_INTRON, temp_df_intron], axis=0)
            df_UUTR = pd.concat([df_UUTR, temp_df_UUTR], axis=0)
            df_DUTR = pd.concat([df_DUTR, temp_df_DUTR], axis=0)
            df_UF = pd.concat([df_UF, temp_df_UF], axis=0)
            df_DF = pd.concat([df_DF, temp_df_DF], axis=0)

            Total_Exon_bases = Total_Exon_bases + Temp_Exon_bases
            Total_Intron_bases = Total_Intron_bases + Temp_Intron_bases
            Total_UUTR_bases = Total_UUTR_bases + Temp_UUTR_bases
            Total_DUTR_bases = Total_DUTR_bases + Temp_DUTR_bases
            Total_UF_bases = Total_UF_bases + Temp_UF_bases
            Total_DF_bases = Total_DF_bases + Temp_DF_bases

            total_records = total_records + rec

    df_REGION.reset_index(inplace=True, drop=True)
    # df_SEQ_T.reset_index(inplace=True, drop=True)
    df_SEQ.reset_index(inplace=True, drop=True)

    df_EXON.reset_index(inplace=True, drop=True)
    df_INTRON.reset_index(inplace=True, drop=True)
    df_UUTR.reset_index(inplace=True, drop=True)
    df_DUTR.reset_index(inplace=True, drop=True)
    df_UF.reset_index(inplace=True, drop=True)
    df_DF.reset_index(inplace=True, drop=True)

    density_exon_GC = df_EXON.sum(axis=0)/total_records
    density_intron_GC = df_INTRON.sum(axis=0)/total_records
    density_UUTR_GC = df_UUTR.sum(axis=0)/total_records
    density_DUTR_GC = df_DUTR.sum(axis=0)/total_records
    density_UF_GC = df_UF.sum(axis=0)/total_records
    density_DF_GC = df_DF.sum(axis=0)/total_records

    print("INTRA-REGIONAL SIGNAL CHARTS AND FILES >>>>>>>>>>>>>>>>>>>")
    region_density_df = pd.read_csv(ELEMENT_DENSITY_FILE)
    region_density_df.drop(['Unnamed: 0'], axis=1, inplace=True)
    region_density_df.rename(columns={'density' : 'region_density'}, inplace=True)

    print('(INTRA_SIGNAL-NOT_CENTERED) Calculating the Region Specific GC content per base pair average across ALL transcripts')

    densities = pd.DataFrame([density_exon_GC, density_intron_GC, density_UUTR_GC, density_DUTR_GC, density_UF_GC, density_DF_GC],
                             index=['exon', 'intron', 'UUTR', 'DUTR', 'UF', 'DF'])

    densities_T = densities.T
    densities_T['position'] = densities_T.index
    densities_T_melt = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt.rename({'value': 'density'}, inplace=True, axis=1)
    densities_T_melt['position'] = densities_T_melt['position'].astype(int)
    densities_T_melt['position'] = densities_T_melt['position'] - 1000
    densities_T_melt.to_csv(GC_PER_REGION_PER_BP_ALL_TRANS_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt, x='position', y="density", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("GC content")
    plt.title("GC content distribution by Region per base position, with average across ALL Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(GC_PER_REGION_PER_BP_ALL_TRANS_CHART)


    print('(INTRA_SIGNAL-CENTERED) Mean Centered GC content Calculation for Region Specific GC content per base pair average across all transcripts......')

    centered_densities_T_melt = copy.deepcopy(densities_T_melt)
    mean_centered_dict_all = centered_densities_T_melt.groupby('variable')['density'].mean().to_dict()
    print(mean_centered_dict_all)
    for e in mean_centered_dict_all.keys():
        centered_densities_T_melt.loc[centered_densities_T_melt['variable'] == e, 'density'] = \
            centered_densities_T_melt.loc[centered_densities_T_melt['variable'] == e, 'density'] - mean_centered_dict_all[e]

    centered_densities_T_melt.to_csv(CENTERED_GC_DENSITY_BY_REGION_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=centered_densities_T_melt, x='position', y="density", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("GC content")
    plt.title("Mean Centered GC content distribution by Region per base position, with average across ALL Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(CENTERED_GC_DENSITY_BY_REGION_CHART)



    print('(CUMULATIVE-INTRA_SIGNAL-NOT_CENTERED) Calculating the Cumulative GC content from the "Region Specific GC content average across all transcripts"...')

    absolute_GC = densities_T_melt.groupby('position').agg({'density': sum}).reset_index()
    absolute_GC.to_csv(DERIVED_ABSOLUTE_GC_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=absolute_GC, x='position', y="density")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Sum(%GC content * %Occurrence)")
    plt.title("Absolute GC content averaged across ALL transcripts")
    plt.xticks(rotation=0)
    plt.savefig(DERIVED_ABSOLUTE_GC_CHART)

    print('(CUMULATIVE-INTRA_SIGNAL-CENTERED) Calculating the Cumulative MEAN CENTERED GC content from the "Region Specific GC content average across all transcripts"......')

    centered_densities_T_melt = centered_densities_T_melt.groupby('position').agg(
        {'density': sum}).reset_index()
    centered_densities_T_melt.to_csv(DERIVED_CENTERED_GC_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=centered_densities_T_melt, x='position', y="density")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Mean Centered GC content")
    plt.title("Absolute Mean Centered GC content averaged across ALL transcripts")
    plt.xticks(rotation=0)
    plt.savefig(DERIVED_CENTERED_GC_CHART)




    print ('(NORMALISED-INTRA_SIGNAL-NOT_CENTERED) Calculating the GC content per Region per bp, average across region specific transcript........')

    GC_per_region_per_bp = densities_T_melt.merge(region_density_df,
                                                   on=['position', 'variable'],
                                                   how='left')

    GC_per_region_per_bp['density_per_region'] = GC_per_region_per_bp['density']/GC_per_region_per_bp['region_density']
    GC_per_region_per_bp.drop(['density', 'region_density'], axis=1, inplace=True)
    GC_per_region_per_bp.to_csv(GC_PER_REGION_PER_BP_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=GC_per_region_per_bp, x='position', y="density_per_region", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content")
    plt.title("%GC content distribution by Region per base position, with average across region specific Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(GC_PER_REGION_PER_BP_CHART)



    print ('(NORMALISED-INTRA_SIGNAL-CENTERED) Calculating the MEAN CENTERED GC content per Region per bp, average across region specific transcript........')

    GC_per_region_per_bp_centered = copy.deepcopy(GC_per_region_per_bp)
    mean_centered_dict = GC_per_region_per_bp_centered.groupby('variable')['density_per_region'].mean().to_dict()

    for e in mean_centered_dict.keys():
        GC_per_region_per_bp_centered.loc[GC_per_region_per_bp_centered['variable'] == e, 'density_per_region'] = \
        GC_per_region_per_bp_centered.loc[GC_per_region_per_bp_centered['variable'] == e, 'density_per_region'] - mean_centered_dict[e]

    GC_per_region_per_bp_centered.to_csv(CENTERED_GC_PER_REGION_PER_BP_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=GC_per_region_per_bp_centered, x='position', y="density_per_region", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Mean Centered GC content")
    plt.title("Mean Centered GC content distribution by Region per base position, with average across region specific Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(CENTERED_GC_PER_REGION_PER_BP_CHART)



    print("INTER-REGIONAL SIGNAL CHARTS AND FILES >>>>>>>>>>>>>>>>>>>")


    print ('Calculating the Avg GC content in each Region........')

    avg_exon_GC = df_EXON.sum(axis=0).sum()/Total_Exon_bases
    avg_intron_GC = df_INTRON.sum(axis=0).sum()/Total_Intron_bases
    avg_UUTR_GC = df_UUTR.sum(axis=0).sum()/Total_UUTR_bases
    avg_DUTR_GC = df_DUTR.sum(axis=0).sum()/Total_DUTR_bases
    avg_UF_GC = df_UF.sum(axis=0).sum()/Total_UF_bases
    avg_DF_GC = df_DF.sum(axis=0).sum()/Total_DF_bases

    avg_GC_df = pd.DataFrame({'element': ['exon', 'intron', 'UUTR', 'DUTR', 'UF', 'DF'],
                              'Avg_GC': [avg_exon_GC, avg_intron_GC, avg_UUTR_GC, avg_DUTR_GC, avg_UF_GC, avg_DF_GC]
                              })

    avg_GC_df.to_csv(AVG_GC_PER_REGION_FILE)
    plt.figure(figsize=(10, 10))
    sns.barplot(x="element", y="Avg_GC", data=avg_GC_df)
    plt.xlabel("Region")
    plt.ylabel("%GC content")
    plt.title("Average GC content per Region")
    plt.savefig(AVG_GC_PER_REGION_CHART)



    print('Calculating average GC content in each region per base pair......... ')

    region_density_df['avg_GC_region_density'] = region_density_df.apply(lambda x: multiply_operation(x['variable'],
                                                                                                   x['region_density'],
                                                                                                   avg_GC_df),
                                                                      axis=1)
    region_density_df.to_csv(AVG_GC_DENSITY_BY_REGION_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=region_density_df, x='position', y="avg_GC_region_density", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("%GC content")
    plt.title("Average %GC distribution by Region per base position")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(AVG_GC_DENSITY_BY_REGION_CHART)

    print('Calculating the Cumulative GC plot ber base from "average GC content in each region per base pair"......... ')

    AVG_absolute_GC = region_density_df.groupby('position').agg({'avg_GC_region_density': sum}).reset_index()

    AVG_absolute_GC.to_csv(DERIVED_CUM_AVG_GC_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=AVG_absolute_GC, x='position', y="avg_GC_region_density")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Sum(%GC content * %Occurrence)")
    plt.title("Average %GC distribution per base position")
    plt.xticks(rotation=0)
    plt.savefig(DERIVED_CUM_AVG_GC_CHART)


    print("NUCLEOSOMAL INTRA-REGIONAL SIGNAL CHARTS AND FILES >>>>>>>>>>>>>>>>>>>")


    print('The NUCLEOSOMAL Calculation for Region Specific GC content average across all transcripts............')

    for k in ['exon', 'intron', 'UUTR', 'DUTR', 'UF', 'DF']:
        densities_T[k] = densities_T[k].rolling(window=147,
                                                min_periods=1,
                                                center=True).mean()

    densities_T_melt_rolling = pd.melt(densities_T, id_vars=['position'])
    densities_T_melt_rolling.rename({'value': 'density'}, inplace=True, axis=1)
    densities_T_melt_rolling['position'] = densities_T_melt_rolling['position'].astype(int)
    densities_T_melt_rolling['position'] = densities_T_melt_rolling['position'] - 1000
    densities_T_melt_rolling.to_csv(NUCLO_GC_PER_REGION_PER_BP_ALL_FILE)


    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt_rolling, x='position', y="density", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Nucleosomal GC content")
    plt.title(" Nucleosomal GC content distribution by Region per base position, with average across ALL Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(NUCLO_GC_PER_REGION_PER_BP_ALL_CHART)



    print('The NUCLEOSOMAL Calculation for the Cumulative GC content from the "Region Specific GC content average across all transcripts"......')

    nuclo_absolute_GC = densities_T_melt_rolling.groupby('position').agg({'density': sum}).reset_index()
    nuclo_absolute_GC.to_csv(NUCLO_DERIVED_ABSOLUTE_GC_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=nuclo_absolute_GC, x='position', y="density")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Sum(%GC content * %Occurrence)")
    plt.title("Absolute Nucleosomal GC content averaged across ALL transcripts")
    plt.xticks(rotation=0)
    plt.savefig(NUCLO_DERIVED_ABSOLUTE_GC_CHART)



    print ('The NUCLEOSOMAL Calculation for the GC content per Region per bp, with average across region specific transcripts........')

    nuclo_region_density_df = pd.read_csv(NUCLO_ELEMENT_DENSITY_FILE)
    nuclo_region_density_df.drop(['Unnamed: 0'], axis=1, inplace=True)
    nuclo_region_density_df.rename(columns={'density': 'region_density'}, inplace=True)

    nuclo_GC_per_region_per_bp = densities_T_melt_rolling.merge(nuclo_region_density_df,
                                                  on=['position', 'variable'],
                                                  how='left')

    nuclo_GC_per_region_per_bp['density_per_region'] = nuclo_GC_per_region_per_bp['density'] / nuclo_GC_per_region_per_bp[
        'region_density']
    nuclo_GC_per_region_per_bp.drop(['density', 'region_density'], axis=1, inplace=True)
    nuclo_GC_per_region_per_bp.to_csv(NUCLO_GC_PER_REGION_PER_BP_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=nuclo_GC_per_region_per_bp, x='position', y="density_per_region", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Nucleosomal GC content")
    plt.title("Nucleosomal GC content distribution by Region per base position, with average across region specific Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(NUCLO_GC_PER_REGION_PER_BP_CHART)


    print ('Calculating the MEAN CENTERED Nucleosomal GC content per Region per bp, average across region specific transcript........')

    nuclo_GC_per_region_per_bp_centered = copy.deepcopy(nuclo_GC_per_region_per_bp)

    nuclo_avg_dict_1 = nuclo_GC_per_region_per_bp_centered.groupby('variable')['density_per_region'].mean().to_dict()

    for e in nuclo_avg_dict_1.keys():
        nuclo_GC_per_region_per_bp_centered.loc[nuclo_GC_per_region_per_bp_centered['variable'] == e, 'density_per_region'] = \
            nuclo_GC_per_region_per_bp_centered.loc[nuclo_GC_per_region_per_bp_centered['variable'] == e, 'density_per_region'] - \
            nuclo_avg_dict_1[e]

    nuclo_GC_per_region_per_bp_centered.to_csv(NUCLO_CENTERED_GC_PER_REGION_PER_BP_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=nuclo_GC_per_region_per_bp_centered, x='position', y="density_per_region", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Mean Centered nucleosomal GC content")
    plt.title("Mean Centered Nucleosomal GC content distribution by Region per base position, with average across region specific Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(NUCLO_CENTERED_GC_PER_REGION_PER_BP_CHART)



    print('The NUCLEOSOMAL Calculation for Region Specific Mean Centered GC content average across all transcripts............')

    nuclo_mean_centered_GC_per_region_mul_denstiy = nuclo_GC_per_region_per_bp_centered.merge(nuclo_region_density_df,
                                                                                  on=['position', 'variable'],
                                                                                  how='left')

    nuclo_mean_centered_GC_per_region_mul_denstiy['density'] = nuclo_mean_centered_GC_per_region_mul_denstiy['density_per_region'] * \
                                                         nuclo_mean_centered_GC_per_region_mul_denstiy['region_density']

    nuclo_mean_centered_GC_per_region_mul_denstiy.drop(['region_density', 'density_per_region'], axis=1, inplace=True)
    nuclo_mean_centered_GC_per_region_mul_denstiy.to_csv(NUCLO_CENTERED_GC_DENSITY_BY_REGION_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=nuclo_mean_centered_GC_per_region_mul_denstiy, x='position', y="density", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Nucleosomal GC content")
    plt.title("Mean Centered Nucleosomal GC content distribution by Region per base position, with average across ALL Transcripts")
    plt.xticks(rotation=0)
    plt.savefig(NUCLO_CENTERED_GC_DENSITY_BY_REGION_CHART)


    print('The NUCLEOSOMAL Calculation for the Cumulative Mean Centered GC content from the "Region Specific Mean Centered GC content average across all transcripts"......')

    densities_T_melt_centered_rolling_ = nuclo_mean_centered_GC_per_region_mul_denstiy.groupby('position').agg({'density': sum}).reset_index()
    densities_T_melt_centered_rolling_.to_csv(NUCLO_DERIVED_CENTERED_GC_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=densities_T_melt_centered_rolling_, x='position', y="density")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Mean Centered Nucleosomal GC content")
    plt.title("Absolute Mean Centered Nucleosomal GC content averaged across ALL transcripts")
    plt.xticks(rotation=0)
    plt.savefig(NUCLO_DERIVED_CENTERED_GC_CHART)


    print("NUCLEOSOMAL INTER-REGIONAL SIGNAL CHARTS AND FILES >>>>>>>>>>>>>>>>>>>")

    print('The NUCLEOSOMAL Calculation for the average GC content in each region per base pair......... ')

    nuclo_region_density_df['avg_GC_region_density'] = nuclo_region_density_df.apply(lambda x: multiply_operation(x['variable'],
                                                                                                      x[
                                                                                                          'region_density'],
                                                                                                      avg_GC_df),
                                                                         axis=1)
    nuclo_region_density_df.to_csv(NUCLO_AVG_GC_DENSITY_BY_REGION_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=nuclo_region_density_df, x='position', y="avg_GC_region_density", hue='variable')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Nucleosomal GC content")
    plt.title("Average %GC distribution by Region per base position")
    plt.xticks(rotation=0)
    # plt.show()
    plt.savefig(NUCLO_AVG_GC_DENSITY_BY_REGION_CHART)

    print('Calculating the Cumulative GC plot ber base from "average GC content in each region per base pair"......... ')

    NUCLO_AVG_absolute_GC = nuclo_region_density_df.groupby('position').agg({'avg_GC_region_density': sum}).reset_index()

    NUCLO_AVG_absolute_GC.to_csv(NUCLO_DERIVED_CUM_AVG_GC_FILE)

    plt.figure(figsize=(20, 10))
    sns.scatterplot(data=NUCLO_AVG_absolute_GC, x='position', y="avg_GC_region_density")
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Sum(%GC content * %Occurrence)")
    plt.title("Average %GC distribution per base position")
    plt.xticks(rotation=0)
    plt.savefig(NUCLO_DERIVED_CUM_AVG_GC_CHART)

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')