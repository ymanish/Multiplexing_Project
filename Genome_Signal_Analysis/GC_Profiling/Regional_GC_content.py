import copy

import pandas as pd
from Bio import SeqIO
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Initialise_GC_profiling import *
import numpy as np


def region_density_func(region_data):
    region_freq_dict = dict()
    for c in COL:
        # print(region_data[str(c)].value_counts().to_dict())
        region_freq_dict[str(c)] = region_data[str(c)].value_counts().to_dict()

    return pd.DataFrame.from_dict(region_freq_dict)


def region_matrix_func(r_data):
    for c in COL:
        r_data[c] = r_data[c].map(REGION_DICT)
    # r_agg_data = amino_acid_density(r_data)

    return r_data


def make_nan(x, value_col):
    if x['region'] == 'F' and x['position']>=0:
        return np.nan
    elif x['region'] != 'F' and x['position']<0:
        return np.nan
    else: return x[value_col]


def broken_axis_plot(data_, plot_title, x_title, y_title, x_column, y_column, hue_column, save_loc, N):

    colormap_region = {'D': 'blue', 'E': 'red', 'I': 'orange', 'U': 'green', 'F': 'cyan', 'f': 'gray'}
    region_labels = {'D': "3'UTR", 'E': "CDS", 'I': "Intron", 'U': "5'UTR", 'F': "5'end", 'f': "3'end"}
    # element_densities = data_.fillna(0)
    element_densities = data_.dropna(subset=[y_column])


    fig, (ax_1, ax_2) = plt.subplots(1, 2, sharey=True, figsize=(10, 7), gridspec_kw={'width_ratios': [0.1, 0.9]})
    fig.subplots_adjust(wspace=0.05)
    ax_1.spines.right.set_visible(False)
    ax_2.spines.left.set_visible(False)

    sns.scatterplot(data=element_densities, x=x_column, y=y_column, hue=hue_column, fc='none',
                    ec=element_densities[hue_column].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_2)
    handles, labels = ax_2.get_legend_handles_labels()
    new_labels = [region_labels.get(item, item) for item in labels]
    ax_2.tick_params(axis='x', labelsize=8)
    ax_2.tick_params(axis='y', colors='white')
    ax_2.set_xticks(ticks=np.arange(element_densities[x_column].min(), element_densities[x_column].max(), 100))
    ax_2.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    ax_2.set_xlim(-100, 1000)
    ax_2.set(xlabel=None, ylabel=None)
    ax_2.grid(True, linestyle='--', linewidth=0.5, axis='both')

    sns.scatterplot(data=element_densities, x=x_column, y=y_column, hue=hue_column, fc='none',
                    ec=element_densities[hue_column].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_1, legend=False)
    ax_1.tick_params(axis='both', labelsize=8)
    ax_1.set_xticks(ticks=np.arange(element_densities[x_column].min(), element_densities[x_column].max(), 100))
    ax_1.set_xlim(-1000, -850)
    ax_1.set(xlabel=None, ylabel=None)
    ax_1.grid(True, linestyle='--', linewidth=0.5, axis='both')

    fig.supxlabel(x_title, fontsize=9, x=0.5, y=0.05)
    fig.supylabel(y_title, fontsize=9, y=0.5, x=0.08)
    fig.suptitle(plot_title, size=10, fontweight='bold', y=0.91, x=0.5)
    d = .01
    kwargs = dict(transform=ax_1.transAxes, color='k', clip_on=False)
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (-d, +d), **kwargs)  ##bottom left
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (1 - d, 1 + d), **kwargs)  ##top left
    kwargs.update(transform=ax_2.transAxes)
    ax_2.plot((-d + .6 * d, +d - .6 * d), (1 - d, 1 + d), **kwargs)  ##top right
    ax_2.plot((-d + .6 * d, +d - .6 * d), (-d, +d), **kwargs)  ##bottom right

    if N:
        ax_2.axvline(73, linestyle='dotted', color='m')
        plt.savefig(save_loc)
    else:
        plt.savefig(save_loc)

    return None


def get_nucleosome_gc_content(reg_data, reg_seq_data):
    nuc_region_data = pd.DataFrame()
    nuc_seq_region_data = pd.DataFrame()

    for i in [1, 2, 3, 4, 5, 6]:

        reg_seq_data_ = pd.DataFrame(np.where(reg_seq_data == i, 1, 0))
        reg_data_ = pd.DataFrame(np.where(reg_data == i, 1, 0))

        reg_seq_data_ = reg_seq_data_.rolling(window=147, min_periods=1, center=True, axis=1).mean()
        reg_data_ = reg_data_.rolling(window=147, min_periods=1, center=True, axis=1).mean()

        #         reg_seq_data_total = pd.DataFrame(reg_seq_data_.sum(axis=0)/obs, columns=[i])
        #         reg_data_total = pd.DataFrame(reg_data_.sum(axis=0)/obs, columns=[i])

        reg_seq_data_total = pd.DataFrame(reg_seq_data_.sum(axis=0), columns=[i])
        reg_data_total = pd.DataFrame(reg_data_.sum(axis=0), columns=[i])

        #         print(reg_seq_data_total)
        nuc_region_data = pd.concat([nuc_region_data, reg_data_total], axis=1)
        nuc_seq_region_data = pd.concat([nuc_seq_region_data, reg_seq_data_total], axis=1)

    return nuc_region_data, nuc_seq_region_data


def combine_absolute_derived_gc(data, nuclo):
    plt.figure(figsize=(10, 7))
    colormap_signals = {'inter_signal': 'red', 'intra_signal': 'blue', 'Derived GC': 'green', 'Actual GC': 'black'}
    data = data.fillna(0)
    data = data[data['variable'].isin(['inter_signal', 'intra_signal', 'Actual GC'])]
    sns.scatterplot(data=data, x='position', y='value', hue='variable', fc='none',
                    ec=data["variable"].map(colormap_signals),
                    palette=colormap_signals, s=1.5, linewidth=1)
    if SEQ_TYPE == 'mRNA_ATG':
        plt.xlabel("position w.r.t TIS", fontsize=9)
    else:
        plt.xlabel("position w.r.t TSS", fontsize=9)

    if nuclo:
        plt.ylabel("Nucleosomal GC content", fontsize=9)
    else:
        plt.ylabel("GC content", fontsize=9)

    plt.xticks(ticks=np.arange(-1000, 1000, 100))
    plt.tick_params(axis='x', labelrotation=45, labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    if nuclo:
        plt.title("Nucleosomal GC content per base-pair position", size=10, fontweight='bold')
    else:
        plt.title("GC content per base-pair position", size=10, fontweight='bold')

    plt.xlim(-1000, 1000)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=8)

    if nuclo:
        plt.axvline(73, linestyle='dotted', color='m')
        plt.savefig(NUCLO_ABSOLUTE_and_SIGNAL_GC_LINE_CHART)
    else:
        plt.savefig(ABSOLUTE_and_SIGNAL_GC_LINE_CHART)


    return None


def scatter_plot_full_axis(data_, plot_title, x_title, y_title, x_column, y_column, hue_column, save_loc, N):
    colormap_region = {'D': 'blue', 'E': 'red', 'I': 'orange', 'U': 'green', 'F': 'cyan', 'f': 'gray'}
    region_labels = {'D': "3'UTR", 'E': "CDS", 'I': "Intron", 'U': "5'UTR", 'F': "5'end", 'f': "3'end"}

    region_density_melted_ = data_.dropna(subset=[y_column])

    fig, ax2 = plt.subplots(figsize=(10, 7))
    sns.scatterplot(data=region_density_melted_, x=x_column, y=y_column, hue=hue_column, fc='none',
                    ec=region_density_melted_["region"].map(colormap_region),
                    palette=colormap_region, s=1.5, linewidth=1, ax=ax2)
    handles, labels = ax2.get_legend_handles_labels()
    new_labels = [region_labels.get(item, item) for item in labels]

    ax2.set_xlabel(x_title, fontsize=9)
    ax2.set_ylabel(y_title, fontsize=9)
    ax2.set_title(plot_title, size=10,
                  fontweight='bold')
    ax2.set_xticks(np.arange(-1000, 1000, 100))
    ax2.tick_params(axis='both', labelrotation=0, labelsize=8)
    ax2.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    ax2.set_xlim(-1000, 1000, 100)
    ax2.grid(linestyle='--', linewidth=0.5)

    if N:
        ax2.axvline(73, linestyle='dotted', color='m')
        plt.savefig(save_loc)
    else:
        plt.savefig(save_loc)

    return None


def main(n):
    df_region = pd.DataFrame()
    df_seq = pd.DataFrame()

    for (seq_record_1, seq_record_2) in zip(
            SeqIO.parse(REGION_FILE_PATH+"\/region_group_"+str(n)+".fasta", "fasta"),
            SeqIO.parse(SEQ_FILE_PATH+"\group_"+str(n)+".fasta", "fasta")):

        header_region = (seq_record_1.id).split('|')
        if SEQ_TYPE == 'mRNA_ATG':
            region = seq_record_1.seq[int(header_region[1]) - 1000:int(header_region[1]) + 1000]

        else:
            region = seq_record_1.seq[:2000]

        temp_region = pd.DataFrame(list(region), index=COL)
        temp_region = temp_region.T
        temp_region['id'] = header_region[0]
        df_region = pd.concat([df_region, temp_region], axis=0)



        header_seq = (seq_record_2.id).split('|')
        if SEQ_TYPE == 'mRNA_ATG':
            seq = seq_record_2.seq[int(header_seq[1]) - 1000:int(header_seq[1]) + 1000]
        else:
            seq = seq_record_2.seq[:2000]

        temp_seq = pd.DataFrame(list(seq), index=COL)
        temp_seq = temp_seq.T
        if GROUP == 'Plant':
            temp_seq['id'] = header_seq[1]
        else:
            temp_seq['id'] = header_seq[2]
        df_seq = pd.concat([df_seq, temp_seq], axis=0)

    df_region = df_region.reset_index(drop=True)
    df_seq = df_seq.reset_index(drop=True)

    df_seq = df_seq.drop('id', axis=1)
    df_region = df_region.drop('id', axis=1)


    region_density_df = region_density_func(df_region)
    region_matrix = region_matrix_func(df_region)

    seq_matrix = df_seq.replace(["A", "C", "G", "T", "N", "S"], [0, 1, 1, 0, 0, 1])
    seq_matrix = seq_matrix.replace('[A-Z]', 0, regex=True)

    print('Region dataframe size: {} and Sequence dataframe size: {}'.format(region_matrix.shape, seq_matrix.shape))
    if region_matrix.shape[0] != seq_matrix.shape[0]:
        print('Warning, Region and Seq does not match')

    df_seq_region = region_matrix.multiply(seq_matrix)
    df_seq_region_count = region_density_func(df_seq_region)

    NUC_region_count, NUC_df_seq_region_count = get_nucleosome_gc_content(region_matrix, df_seq_region)


    return region_density_df, len(df_seq), df_seq_region_count, NUC_region_count, NUC_df_seq_region_count


# def Region_Seq_multipliction(DF_region, region_code, DF_seq):
#     df_region_matrix = DF_region.replace(["U", "E", "D", "I", "F", "f"], region_code)
#
#     print('Region dataframe size: {} and Sequence dataframe size: {}'.format(DF_region.shape, DF_seq.shape))
#     if DF_region.shape[0] != DF_seq.shape[0]:
#         print('Warning, Region and Seq does not match')
#     DF_seq = DF_seq.replace(["A", "C", "G", "T", "N", "S"], [0, 1, 1, 0, 0, 1])
#     DF_seq = DF_seq.replace('[A-Z]', 0, regex=True)
#     df_seq_region = df_region_matrix.multiply(DF_seq)
#     total_region_bases = df_region_matrix.sum(axis=0).sum()
#     return df_seq_region, total_region_bases
#
#
# def density_GC_Region_and_Avg(df_region, region_code, df_seq):
#     df_region_matrix = df_region.replace(["U", "E", "D", "I", "F", "f"], region_code)
#
#     print('Region dataframe size: {} and Sequence dataframe size: {}'.format(df_region.shape, df_seq.shape))
#     if df_region.shape[0] != df_seq.shape[0]:
#         print('Warning, Region and Seq does not match')
#     df_seq = df_seq.replace(["A", "C", "G", "T", "N", "S"], [0, 1, 1, 0, 0, 1])
#     df_seq_region = df_region_matrix.multiply(df_seq)
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

def multiply_operation(el, density, tab):
    return density * tab.loc[tab['element'] == el, 'Avg_GC'].values[0]


if __name__ == "__main__":

    start = time.perf_counter()

    region_density = pd.DataFrame()
    region_gc_density = pd.DataFrame()
    NUC_region_df = pd.DataFrame()
    NUC_region_gc_count = pd.DataFrame()

    total_records = 0

    print('Reading the Sequence and Region files........')
    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files + 1)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):
            temp_region_density, rec, temp_region_gc, temp_nuc_region_count, temp_nuc_region_gc_count = j.result()
            region_density = pd.concat([region_density, temp_region_density], axis=0)
            total_records = total_records + rec

            region_gc_density = pd.concat([region_gc_density, temp_region_gc], axis=0)
            NUC_region_df = pd.concat([NUC_region_df, temp_nuc_region_count], axis=0)
            NUC_region_gc_count = pd.concat([NUC_region_gc_count, temp_nuc_region_gc_count], axis=0)




    print('CALCULATING THE REGION DENSITY')

    region_density.index.name = 'region'
    region_density_grouped = region_density.groupby('region').sum()
    region_density_grouped.reset_index(inplace=True, drop=False)


    region_density_grouped.loc[:, 'total_region'] = region_density_grouped.loc[:, COL].sum(axis=1)

    total_count_regions = region_density_grouped[['region', 'total_region']]

    region_density_grouped.drop(['total_region'], inplace=True, axis=1)
    region_density_grouped[COL] = region_density_grouped[COL] / total_records
    region_density_melted = pd.melt(region_density_grouped, id_vars=['region'], value_vars=COL)
    region_density_melted.rename(columns={'value': 'region_density'}, inplace=True)
    region_density_melted['variable'] = region_density_melted['variable'].astype(int)
    region_density_melted.rename(columns={'variable': 'position'}, inplace=True)
    region_density_melted['position'] = region_density_melted['position'] - 1000

    region_density_melted.to_csv(ELEMENT_DENSITY_FILE)

    if SEQ_TYPE == 'mRNA_ATG':

        scatter_plot_full_axis(data_=region_density_melted,
                     plot_title="Region density per base-pair position",
                     x_title='position w.r.t TIS',
                     y_title='density',
                     x_column='position',
                     y_column='region_density',
                     hue_column='region',
                     save_loc=ELEMENT_DENSITY_LINE_CHART, N=False)


    else:
        broken_axis_plot(data_=region_density_melted,
                     plot_title="Region density per base-pair position",
                     x_title='position w.r.t TSS',
                     y_title='density',
                     x_column='position',
                     y_column='region_density',
                     hue_column='region',
                     save_loc=ELEMENT_DENSITY_LINE_CHART, N=False)



    region_density_melted.sort_values(by=['position'], inplace=True, ascending=True)
    region_density_melted.reset_index(drop=True, inplace=True)

    nuc_region_data = region_density_melted.groupby('region')['region_density'].rolling(window=147,
                                                                                        min_periods=1,
                                                                                        center=True).mean().reset_index()
    nuc_region_data.rename(columns={'level_1': 'index', 'region_density': 'nuc_region_density'}, inplace=True)

    region_density_melted['index'] = region_density_melted.index
    nuclo_region_density_melted = region_density_melted.merge(nuc_region_data, on=['index', 'region'], how='left')
    nuclo_region_density_melted.drop(columns=['index'], inplace=True)

    nuclo_region_density_melted.to_csv(NUCLO_ELEMENT_DENSITY_FILE)

    if SEQ_TYPE == 'mRNA_ATG':

        scatter_plot_full_axis(data_=nuclo_region_density_melted,
                               plot_title="Region density per base-pair position",
                               x_title='position w.r.t TIS',
                               y_title='density',
                               x_column='position',
                               y_column='nuc_region_density',
                               hue_column='region',
                               save_loc=NUCLO_ELEMENT_DENSITY_LINE_CHART, N=True)


    else:
        broken_axis_plot(data_ = nuclo_region_density_melted,
                         plot_title = "Region density per base-pair position",
                         x_title='position w.r.t TSS',
                         y_title = 'density',
                         x_column = 'position',
                         y_column = 'nuc_region_density',
                         hue_column = 'region',
                         save_loc=NUCLO_ELEMENT_DENSITY_LINE_CHART, N=True)




    print('CALCULATING THE AVG GC CONTENT IN EACH REGION')

    region_gc_density.index.name = 'region'

    region_gc_density_grouped = region_gc_density.groupby('region').sum()
    region_gc_density_grouped.reset_index(inplace=True, drop=False)

    region_gc_density_grouped['region'] = region_gc_density_grouped['region'].map(OPP_REGION_DICT)
    region_gc_density_grouped.loc[:, 'gc_count'] = region_gc_density_grouped.loc[:, COL].sum(axis=1)

    region_gc_density_grouped.dropna(subset=['region'], axis=0, inplace=True)
    avg_gc_region_df = region_gc_density_grouped[['region', 'gc_count']]
    region_gc_density_grouped.drop(['gc_count'], axis=1, inplace=True)

    merged_avg_gc_region = avg_gc_region_df.merge(total_count_regions , on='region', how='left')
    merged_avg_gc_region['avg_gc'] = merged_avg_gc_region['gc_count'] / merged_avg_gc_region['total_region']
    merged_avg_gc_region['avg_gc'] = merged_avg_gc_region['avg_gc'].replace([np.inf, -np.inf], np.nan)

    merged_avg_gc_region.to_csv(AVG_GC_PER_REGION_FILE) #### file location

    merged_avg_gc_region_ = copy.deepcopy(merged_avg_gc_region)

    merged_avg_gc_region_.replace({'D': "3'UTR",
                                   'E': "CDS",
                                   'I': "Intron",
                                   'U': "5'UTR",
                                   'F': "5'end",
                                   'f': "3'end"}, inplace=True)

    fig, ax3 = plt.subplots(figsize=(10, 7))
    sns.barplot(data=merged_avg_gc_region_,
                x="region", y="avg_gc", ax=ax3, color='gray')

    ax3.set_xlabel("Region", fontsize=9)
    ax3.set_ylabel("GC content", fontsize=9)
    ax3.set_title("Average GC content per region", size=10,
                  fontweight='bold')
    ax3.tick_params(axis='both', labelrotation=0, labelsize=8)
    ax3.grid(linestyle='--', linewidth=0.5, axis='y')
    plt.savefig(AVG_GC_REGION_CHART)



    print('INTER-REGION SIGNAL PER REGION, CALCULATING THE AVG GC CONTENT FROM EACH REGION PER BASE PAIR')

    inter_region_signal = region_density_melted.merge(merged_avg_gc_region[['region', 'avg_gc']],
                                               on='region',
                                               how='left')
    inter_region_signal['signal'] = inter_region_signal['region_density'] * inter_region_signal['avg_gc']
    inter_region_signal.to_csv(INTER_REGION_SIGNAL_DIST_FILE)


    if SEQ_TYPE == 'mRNA_ATG':
        scatter_plot_full_axis(data_=inter_region_signal,
                               plot_title="Average GC content distribution per region per base-pair position (Inter-region signal)",
                               x_title='position w.r.t TIS',
                               y_title='GC content',
                               x_column='position',
                               y_column='signal',
                               hue_column='region',
                               save_loc=INTER_SIGNAL_REGION_DIST_CHART, N=False)


    else:
        broken_axis_plot(data_=inter_region_signal,
                         plot_title="Average GC content distribution per region per base-pair position (Inter-region signal)",
                         x_title='position w.r.t TSS',
                         y_title='GC content',
                         x_column='position',
                         y_column='signal',
                         hue_column='region',
                         save_loc=INTER_SIGNAL_REGION_DIST_CHART, N=False)




    print('TOTAL INTER-REGION SIGNAL PER BASE PAIR')

    total_inter_region_signal = inter_region_signal.groupby('position')['signal'].sum()
    total_inter_region_signal = total_inter_region_signal.reset_index()
    total_inter_region_signal.rename(columns={'signal': 'inter_signal'}, inplace=True)
    total_inter_region_signal.to_csv(INTER_REGION_SIGNAL_FILE)




    print('INTRA-SIGNAL PER REGION,  Calculating the Region Specific GC content per base pair average across ALL transcripts')

    region_gc_density_grouped[COL] = region_gc_density_grouped[COL] / total_records
    region_gc_density_grouped_melted = pd.melt(region_gc_density_grouped, id_vars=['region'], value_vars=COL)
    region_gc_density_grouped_melted.rename(columns={'value': 'intra_signal',
                                                     'variable': 'position'}, inplace=True)

    region_gc_density_grouped_melted['position'] = region_gc_density_grouped_melted['position'].astype(int)
    region_gc_density_grouped_melted['position'] = region_gc_density_grouped_melted['position'] - 1000

    region_gc_density_grouped_melted.to_csv(INTRA_REGION_SIGNAL_DIST_FILE) ### Intermediate File before centering the data

    intra_signal_dist_data = region_gc_density_grouped_melted.merge(inter_region_signal[['position', 'region', 'signal']],
                                            on=['position', 'region'],
                                           how='left')

    intra_signal_dist_data['actual_intra_signal'] = intra_signal_dist_data['intra_signal']-intra_signal_dist_data['signal']

    intra_signal_dist_data.to_csv(CENTERED_INTRA_REGION_SIGNAL_DIST_FILE)


    if SEQ_TYPE == 'mRNA_ATG':
        x_axis_title = "position w.r.t TIS"
    else:
        x_axis_title = "position w.r.t TSS"

    scatter_plot_full_axis(data_=intra_signal_dist_data,
                           plot_title="Region wise GC content distribution per base-pair position (Intra-region signal)",
                           x_title=x_axis_title,
                           y_title='GC content',
                           x_column='position',
                           y_column='actual_intra_signal',
                           hue_column='region',
                           save_loc=INTRA_SIGNAL_REGION_DIST_CHART, N=False)
    #
    # colormap_region = {'D': 'blue', 'E': 'red', 'I': 'orange', 'U': 'green', 'F': 'cyan', 'f': 'gray'}
    # region_labels = {'D': "3'UTR", 'E': "CDS", 'I': "Intron", 'U': "5'UTR", 'F': "5'end", 'f': "3'end"}
    #
    # intra_signal_dist_data_ = intra_signal_dist_data.dropna(subset=['actual_intra_signal'])
    #
    # fig, ax2 = plt.subplots(figsize=(10, 7))
    # sns.scatterplot(data=intra_signal_dist_data_, x='position', y="actual_intra_signal", hue='region', fc='none',
    #                 ec=intra_signal_dist_data_["region"].map(colormap_region),
    #                 palette=colormap_region, s=1.5, linewidth=1, ax=ax2)
    # handles, labels = ax2.get_legend_handles_labels()
    # new_labels = [region_labels.get(item, item) for item in labels]
    #
    # ax2.set_xlabel("position w.r.t TSS", fontsize=9)
    # ax2.set_ylabel("GC content", fontsize=9)
    # ax2.set_title("Region wise GC content distribution per base-pair position (Intra-region signal)", size=10, fontweight='bold')
    # ax2.set_xticks(np.arange(-1000, 1000, 100))
    # ax2.tick_params(axis='both', labelrotation=0, labelsize=8)
    # ax2.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    # ax2.set_xlim(-1000, 1000, 100)
    # ax2.grid(linestyle='--', linewidth=0.5)
    #
    # plt.savefig(INTRA_SIGNAL_REGION_DIST_CHART)


    print('NORMALISED INTRA-REGION SIGNAL, CALCULATING THE GC CONTENT PER EACH REGION PER BASE PAIR')

    region_gc_density_grouped_melted_joined_region_den = region_gc_density_grouped_melted.merge(region_density_melted,
                                                                                    on=['region', 'position'],
                                                                                    how='left')
    region_gc_density_grouped_melted_joined_region_den['normalised_intra_signal'] = region_gc_density_grouped_melted_joined_region_den[
                                                                              'intra_signal'] / \
                                                                          region_gc_density_grouped_melted_joined_region_den[
                                                                              'region_density']

    region_gc_density_grouped_melted_joined_region_den['normalised_intra_signal'] = region_gc_density_grouped_melted_joined_region_den['normalised_intra_signal'].replace([np.inf, -np.inf], np.nan)
    region_gc_density_grouped_melted_joined_region_den.to_csv(NORM_INTRA_REGION_SIGNAL_DIST_FILE) ### Intermittent file

    cent_region_gc_density_grouped_melted_joined_region_den = region_gc_density_grouped_melted_joined_region_den[['region',
                                                        'position',
                                                        'normalised_intra_signal']].merge(merged_avg_gc_region[['region', 'avg_gc']],
                                                                                          on='region',
                                                                                          how='left')
    cent_region_gc_density_grouped_melted_joined_region_den['cent_norm_intra_signal'] = cent_region_gc_density_grouped_melted_joined_region_den['normalised_intra_signal'] - cent_region_gc_density_grouped_melted_joined_region_den['avg_gc']

    cent_region_gc_density_grouped_melted_joined_region_den.to_csv(CENTERED_NORM_INTRA_REGION_SIGNAL_DIST_FILE)

    if SEQ_TYPE == 'mRNA_ATG':
        x_axis_title = "position w.r.t TIS"
    else:
        x_axis_title = "position w.r.t TSS"

    scatter_plot_full_axis(data_=cent_region_gc_density_grouped_melted_joined_region_den,
                           plot_title="Region wise GC content distribution per region-specific base-pair position (Normalised Intra-region signal)",
                           x_title=x_axis_title,
                           y_title='GC content',
                           x_column='position',
                           y_column='cent_norm_intra_signal',
                           hue_column='region',
                           save_loc=NORM_INTRA_SIGNAL_REGION_DIST_CHART, N=False)
    #
    # fig, ax1 = plt.subplots(figsize=(10, 7))
    # # cent_region_gc_density_grouped_melted_joined_region_den = cent_region_gc_density_grouped_melted_joined_region_den.fillna(0)
    # cent_region_gc_density_grouped_melted_joined_region_den_ = cent_region_gc_density_grouped_melted_joined_region_den.dropna(subset=['cent_norm_intra_signal'])
    #
    # sns.scatterplot(data=cent_region_gc_density_grouped_melted_joined_region_den_, x='position', y="cent_norm_intra_signal", hue='region',
    #                 fc='none', ec=cent_region_gc_density_grouped_melted_joined_region_den_["region"].map(colormap_region),
    #                 palette=colormap_region, s=1.5, linewidth=1, ax=ax1)
    # handles, labels = ax1.get_legend_handles_labels()
    # new_labels = [region_labels.get(item, item) for item in labels]
    # plt.xlabel("position w.r.t TSS", fontsize=9)
    # plt.ylabel("GC content", fontsize=9)
    # plt.title("Region wise GC content distribution per region-specific base-pair position (Normalised Intra-region signal)", size=10, fontweight='bold')
    # plt.xticks(np.arange(-1000, 1000, 100))
    # plt.tick_params(axis='both', labelrotation=0, labelsize=8)
    # plt.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    # plt.xlim(-1000, 1000, 100)
    # plt.grid(linestyle='--', linewidth=0.5)
    # plt.savefig(NORM_INTRA_SIGNAL_REGION_DIST_CHART)


    print('CALCULATING THE TOTAL INTRA REGIONAL SIGNAL')

    total_intra_region_signal = intra_signal_dist_data.groupby('position')['actual_intra_signal'].sum()
    total_intra_region_signal = total_intra_region_signal.reset_index()
    total_intra_region_signal.rename(columns={'actual_intra_signal': 'intra_signal'}, inplace=True)
    total_intra_region_signal.to_csv(INTRA_REGION_SIGNAL_FILE)


    print('N'
          'U'
          'C'
          'L'
          'E'
          'O'
          'S'
          'O'
          'M'
          'A'
          'L###########################')

    print('CALCUALTING NUCLEOSOME AVG GC CONTENT IN EACH REGION')

    region_code = NUC_region_df.columns
    NUC_region_df['position'] = NUC_region_df.index
    NUC_region_gc_count['position'] = NUC_region_gc_count.index

    NUC_region_df = NUC_region_df.groupby(['position']).sum().reset_index()
    NUC_region_gc_count = NUC_region_gc_count.groupby(['position']).sum().reset_index()

    nuc_total_gc = NUC_region_gc_count.sum(axis=0).reset_index()
    nuc_total_gc.rename(columns={0: 'total_gcs', 'index': 'region'}, inplace=True)

    nuc_total_bases = NUC_region_df.sum(axis=0).reset_index()
    nuc_total_bases.rename(columns={0: 'total_bases', 'index': 'region'}, inplace=True)

    nuc_avg_gc_region = nuc_total_gc.merge(nuc_total_bases, on='region', how='left')
    nuc_avg_gc_region['avg_GC'] = nuc_avg_gc_region['total_gcs'] / nuc_avg_gc_region['total_bases']

    nuc_avg_gc_region['region'] = nuc_avg_gc_region['region'].map(OPP_REGION_DICT)
    nuc_avg_gc_region.dropna(inplace=True)

    nuc_avg_gc_region.to_csv(NUCLO_AVG_GC_PER_REGION_FILE)

    nuc_avg_gc_region_ = copy.deepcopy(nuc_avg_gc_region)

    nuc_avg_gc_region_.replace({'D': "3'UTR",
                                   'E': "CDS",
                                   'I': "Intron",
                                   'U': "5'UTR",
                                   'F': "5'end",
                                   'f': "3'end"}, inplace=True)

    fig, ax3 = plt.subplots(figsize=(10, 7))
    sns.barplot(data=nuc_avg_gc_region_,
                x="region", y="avg_GC", ax=ax3, color='gray')
    ax3.set_xlabel("Region", fontsize=9)
    ax3.set_ylabel("Nucleosomal GC content", fontsize=9)
    ax3.set_title("Average nucleosomal GC content per region", size=10,
                  fontweight='bold')
    ax3.tick_params(axis='both', labelrotation=0, labelsize=8)
    ax3.grid(linestyle='--', linewidth=0.5, axis='y')
    plt.savefig(NUCLO_AVG_GC_REGION_CHART)



    print('NUCLEOSOMAL INTER-REGION SIGNAL PER REGION, CALCULATING THE AVG GC CONTENT FROM EACH REGION PER BASE PAIR')

    nuclo_inter_region_signal = nuclo_region_density_melted.merge(nuc_avg_gc_region[['region', 'avg_GC']], on='region',
                                                                  how='left')

    nuclo_inter_region_signal['inter_signal'] = nuclo_inter_region_signal['nuc_region_density'] * nuclo_inter_region_signal[
        'avg_GC']

    nuclo_inter_region_signal.to_csv(NUCLO_INTER_REGION_SIGNAL_DIST_FILE)

    if SEQ_TYPE == 'mRNA_ATG':
        scatter_plot_full_axis(data_=nuclo_inter_region_signal,
                               plot_title="Average nucleosomal GC content distribution per region per base-pair position (Inter-region signal)",
                               x_title='position w.r.t TIS',
                               y_title='Nucleosomal GC content',
                               x_column='position',
                               y_column='inter_signal',
                               hue_column='region',
                               save_loc=NUCLO_INTER_SIGNAL_REGION_DIST_CHART, N=True)
    else:

        broken_axis_plot(data_=nuclo_inter_region_signal,
                         plot_title="Average nucleosomal GC content distribution per region per base-pair position (Inter-region signal)",
                         x_title='position w.r.t TSS',
                         y_title='Nucleosomal GC content',
                         x_column='position',
                         y_column='inter_signal',
                         hue_column='region',
                         save_loc=NUCLO_INTER_SIGNAL_REGION_DIST_CHART, N=True)



    print('TOTAL NUCLEOSOMAL INTER-REGION SIGNAL PER BASE PAIR')

    total_nuclo_inter_region_signal = nuclo_inter_region_signal.groupby('position')['inter_signal'].sum()
    total_nuclo_inter_region_signal = total_nuclo_inter_region_signal.reset_index()
    total_nuclo_inter_region_signal.to_csv(NUCLO_INTER_REGION_SIGNAL_FILE)




    print('NUCLEOSOMAL INTRA-SIGNAL PER REGION, Calculating the Region Specific '
          'GC content per base pair average across ALL transcripts')

    NUC_intra_region_gc_count = copy.deepcopy(NUC_region_gc_count)
    NUC_intra_region_gc_count[region_code] = NUC_intra_region_gc_count[region_code] / total_records

    NUC_intra_region_gc_count_melted = pd.melt(NUC_intra_region_gc_count, id_vars='position', value_vars=region_code)
    NUC_intra_region_gc_count_melted['variable'] = NUC_intra_region_gc_count_melted['variable'].map(OPP_REGION_DICT)
    NUC_intra_region_gc_count_melted.rename(columns={'variable': 'region', 'value': 'intra_signal'}, inplace=True)

    NUC_intra_region_gc_count_melted['position'] = NUC_intra_region_gc_count_melted['position'] - 1000

    NUC_intra_region_gc_count_melted.to_csv(
        NUCLO_INTRA_REGION_SIGNAL_DIST_FILE)  ### Intermediate File before centering the data

    cent_NUC_intra_region_gc_count_melted = NUC_intra_region_gc_count_melted.merge(
        nuclo_inter_region_signal[['position', 'region', 'inter_signal']],
        on=['position', 'region'],
        how='left')

    cent_NUC_intra_region_gc_count_melted['actual_intra_signal'] = cent_NUC_intra_region_gc_count_melted['intra_signal'] - cent_NUC_intra_region_gc_count_melted[
        'inter_signal']


    cent_NUC_intra_region_gc_count_melted.to_csv(NUCLO_CENTERED_INTRA_REGION_SIGNAL_DIST_FILE)


    if SEQ_TYPE == 'mRNA_ATG':
        x_axis_title = "position w.r.t TIS"
    else:
        x_axis_title = "position w.r.t TSS"

    scatter_plot_full_axis(data_=cent_NUC_intra_region_gc_count_melted,
                       plot_title="Region wise nucleosomal GC content distribution per base-pair position (Intra-region signal)",
                       x_title=x_axis_title,
                       y_title="Nucleosomal GC content",
                       x_column='position',
                       y_column='actual_intra_signal',
                       hue_column='region',
                       save_loc=NUCLO_INTRA_SIGNAL_REGION_DIST_CHART, N=True)




    # colormap_region = {'D': 'blue', 'E': 'red', 'I': 'orange', 'U': 'green', 'F': 'cyan', 'f': 'gray'}
    # region_labels = {'D': "3'UTR", 'E': "CDS", 'I': "Intron", 'U': "5'UTR", 'F': "5'end", 'f': "3'end"}
    #
    # cent_NUC_intra_region_gc_count_melted_ = cent_NUC_intra_region_gc_count_melted.dropna(subset=['actual_intra_signal'])
    #
    # fig, ax2 = plt.subplots(figsize=(10, 7))
    # sns.scatterplot(data=cent_NUC_intra_region_gc_count_melted_, x='position', y="actual_intra_signal", hue='region',
    #                 fc='none',
    #                 ec=cent_NUC_intra_region_gc_count_melted_["region"].map(colormap_region),
    #                 palette=colormap_region, s=1.5, linewidth=1, ax=ax2)
    # handles, labels = ax2.get_legend_handles_labels()
    # new_labels = [region_labels.get(item, item) for item in labels]
    # print(handles, labels, new_labels)
    # ax2.set_xlabel("position w.r.t TSS", fontsize=9)
    # ax2.set_ylabel("Nucleosomal GC content", fontsize=9)
    # ax2.set_title("Region wise nucleosomal GC content distribution per base-pair position (Intra-region signal)",
    #               size=10,
    #               fontweight='bold')
    # ax2.set_xticks(np.arange(-1000, 1000, 100))
    # ax2.tick_params(axis='both', labelrotation=0, labelsize=8)
    # ax2.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    # ax2.set_xlim(-1000, 1000, 100)
    # ax2.grid(linestyle='--', linewidth=0.5)
    #
    # plt.axvline(73, linestyle='dotted', color='m')
    # plt.savefig(NUCLO_INTRA_SIGNAL_REGION_DIST_CHART)

    print('NORMALISED NUCLEOSOMAL INTRA-REGION SIGNAL, CALCULATING THE GC CONTENT PER EACH REGION PER BASE PAIR')

    norm_NUC_intra_region_gc_count_melted = NUC_intra_region_gc_count_melted.merge(nuclo_region_density_melted,
                                                                                   on=['position', 'region'],
                                                                                   how='left')
    norm_NUC_intra_region_gc_count_melted['norm_intra_signal'] = norm_NUC_intra_region_gc_count_melted['intra_signal'] / \
                                                                 norm_NUC_intra_region_gc_count_melted[
                                                                     'nuc_region_density']
    norm_NUC_intra_region_gc_count_melted['norm_intra_signal'] = norm_NUC_intra_region_gc_count_melted[
        'norm_intra_signal'].replace([np.inf, -np.inf], np.nan)

    norm_NUC_intra_region_gc_count_melted.to_csv(NUCLO_NORM_INTRA_REGION_SIGNAL_DIST_FILE)  ### Intermittent file

    cent_norm_NUC_intra_region_gc_count_melted = norm_NUC_intra_region_gc_count_melted[
        ['region',
         'position',
         'norm_intra_signal']].merge(nuc_avg_gc_region[['region', 'avg_GC']],
                                           on='region',
                                           how='left')
    cent_norm_NUC_intra_region_gc_count_melted['cent_norm_intra_signal'] = \
    cent_norm_NUC_intra_region_gc_count_melted['norm_intra_signal'] - \
    cent_norm_NUC_intra_region_gc_count_melted['avg_GC']

    cent_norm_NUC_intra_region_gc_count_melted.to_csv(NUCLO_CENTERED_NORM_INTRA_REGION_SIGNAL_DIST_FILE)

    if SEQ_TYPE == 'mRNA_ATG':
        x_axis_title = "position w.r.t TIS"
    else:
        x_axis_title = "position w.r.t TSS"

    scatter_plot_full_axis(data_=cent_norm_NUC_intra_region_gc_count_melted,
                           plot_title="Region wise nucleosomal GC content distribution per region-specific base-pair position (Normalised Intra-region signal)",
                           x_title=x_axis_title,
                           y_title="Nucleosomal GC content",
                           x_column='position',
                           y_column='cent_norm_intra_signal',
                           hue_column='region',
                           save_loc=NUCLO_NORM_INTRA_SIGNAL_REGION_DIST_CHART, N=True)
    # fig, ax1 = plt.subplots(figsize=(10, 7))
    # cent_norm_NUC_intra_region_gc_count_melted_ = cent_norm_NUC_intra_region_gc_count_melted.dropna(subset=['cent_norm_intra_signal'])
    #
    # sns.scatterplot(data=cent_norm_NUC_intra_region_gc_count_melted_, x='position', y="cent_norm_intra_signal",
    #                 hue='region',
    #                 fc='none', ec=cent_norm_NUC_intra_region_gc_count_melted_["region"].map(colormap_region),
    #                 palette=colormap_region, s=1.5, linewidth=1, ax=ax1)
    # handles, labels = ax1.get_legend_handles_labels()
    # new_labels = [region_labels.get(item, item) for item in labels]
    # plt.xlabel("position w.r.t TSS", fontsize=9)
    # plt.ylabel("Nucleosomal GC content", fontsize=9)
    # plt.title("Region wise nucleosomal GC content distribution per region-specific base-pair position (Normalised Intra-region signal)",
    #     size=10, fontweight='bold')
    # plt.xticks(np.arange(-1000, 1000, 100))
    # plt.tick_params(axis='both', labelrotation=0, labelsize=8)
    # plt.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    # plt.xlim(-1000, 1000, 100)
    # plt.grid(linestyle='--', linewidth=0.5)
    #
    # plt.axvline(73, linestyle='dotted', color='m')
    # plt.savefig(NUCLO_NORM_INTRA_SIGNAL_REGION_DIST_CHART)



    print('CALCULATING THE TOTAL INTRA REGIONAL SIGNAL')

    total_nuclo_intra_region_signal = cent_NUC_intra_region_gc_count_melted.groupby('position')[
        'actual_intra_signal'].sum()
    total_nuclo_intra_region_signal = total_nuclo_intra_region_signal.reset_index()
    total_nuclo_intra_region_signal.rename(columns={'actual_intra_signal': 'intra_signal'}, inplace=True)
    total_nuclo_intra_region_signal.to_csv(NUCLO_INTRA_REGION_SIGNAL_FILE)



    print('Combine the Absolute GC, Intra Signal and Inter Signal')
    ##### For actual values
    DERIVED_CUM_GC_PER_BP = total_inter_region_signal.merge(total_intra_region_signal, on='position', how='left')
    DERIVED_CUM_GC_PER_BP['Derived GC'] = DERIVED_CUM_GC_PER_BP['intra_signal'] + DERIVED_CUM_GC_PER_BP['inter_signal']
    DERIVED_CUM_GC_PER_BP_MELT = pd.melt(DERIVED_CUM_GC_PER_BP, id_vars=['position'])

    absolute_gc = pd.read_csv(ABSOLUTE_GC_FILE)
    absolute_gc.drop(['Unnamed: 0'], axis=1, inplace=True)
    absolute_gc['variable'] = 'Actual GC'
    absolute_gc.rename(columns={'0': 'value'}, inplace=True)
    absolute_gc['position'] = absolute_gc.index
    absolute_gc['position'] = absolute_gc['position'] - 1000
    Signal_DF = pd.concat([DERIVED_CUM_GC_PER_BP_MELT, absolute_gc], axis=0)
    Signal_DF.reset_index(inplace=True, drop=True)

    combine_absolute_derived_gc(data = Signal_DF, nuclo=False)

    #####For Nucleosomal gc values

    NUCLO_DERIVED_CUM_GC_PER_BP = total_nuclo_inter_region_signal.merge(total_nuclo_intra_region_signal, on='position', how='left')
    NUCLO_DERIVED_CUM_GC_PER_BP['Derived GC'] = NUCLO_DERIVED_CUM_GC_PER_BP['intra_signal'] + NUCLO_DERIVED_CUM_GC_PER_BP['inter_signal']
    NUCLO_DERIVED_CUM_GC_PER_BP_MELT = pd.melt(NUCLO_DERIVED_CUM_GC_PER_BP, id_vars=['position'])

    nuclo_absolute_gc = pd.read_csv(NUCLO_ABSOLUTE_GC_FILE)
    nuclo_absolute_gc.drop(['Unnamed: 0'], axis=1, inplace=True)
    nuclo_absolute_gc['variable'] = 'Actual GC'
    nuclo_absolute_gc.rename(columns={'0': 'value'}, inplace=True)
    nuclo_absolute_gc['position'] = nuclo_absolute_gc.index
    nuclo_absolute_gc['position'] = nuclo_absolute_gc['position'] - 1000
    NUCLO_Signal_DF = pd.concat([NUCLO_DERIVED_CUM_GC_PER_BP_MELT, nuclo_absolute_gc], axis=0)
    NUCLO_Signal_DF.reset_index(inplace=True, drop=True)

    combine_absolute_derived_gc(data = NUCLO_Signal_DF, nuclo=True)

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')