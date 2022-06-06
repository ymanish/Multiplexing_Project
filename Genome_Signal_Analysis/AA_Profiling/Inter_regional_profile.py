import pandas as pd
import numpy as np
from Bio import SeqIO
import copy
import seaborn as sns
import concurrent.futures
import time
import matplotlib.pyplot as plt
from Initialise_AA_profiling import *


def amino_acid_density(region_data):
    aa_freq_dict = dict()
    for c in COL:
        # print(region_data[str(c)].value_counts().to_dict())
        aa_freq_dict[str(c)] = region_data[str(c)].value_counts().to_dict()

    return pd.DataFrame.from_dict(aa_freq_dict)
def aa_matrix(r_data):
    for c in COL:
        r_data[c] = r_data[c].map(AA_DICT)
    # r_agg_data = amino_acid_density(r_data)

    return r_data



def make_nan(x, value_col):
    if x['amino_acid'] == 'UF' and x['position']>=0:
        return np.nan
    elif x['amino_acid'] != 'UF' and x['position']<0:
        return np.nan
    else: return x[value_col]

def get_nucleosome_gc_content(reg_data, reg_seq_data):
    nuc_region_data = pd.DataFrame()
    nuc_seq_region_data = pd.DataFrame()

    for i in list(OPP_AA_DICT.keys()):
        reg_seq_data_ = pd.DataFrame(np.where(reg_seq_data == i, 1, 0))
        reg_data_ = pd.DataFrame(np.where(reg_data == i, 1, 0))

        reg_seq_data_ = reg_seq_data_.rolling(window=147, min_periods=1, center=True, axis=1).mean()
        reg_data_ = reg_data_.rolling(window=147, min_periods=1, center=True, axis=1).mean()

        reg_seq_data_total = pd.DataFrame(reg_seq_data_.sum(axis=0), columns=[i])
        reg_data_total = pd.DataFrame(reg_data_.sum(axis=0), columns=[i])

        nuc_region_data = pd.concat([nuc_region_data, reg_data_total], axis=1)
        nuc_seq_region_data = pd.concat([nuc_seq_region_data, reg_seq_data_total], axis=1)

    return nuc_region_data, nuc_seq_region_data

def main(n):
    df_region = pd.DataFrame()
    df_seq = pd.DataFrame()

    if SEQ_TYPE == 'Amino_Acid':
        INPUT_FILE_AA = AA_FILE_PATH
        INPUT_FILE_SEQ = pre_mrna_SEQ_FILE_PATH

    else:
        INPUT_FILE_AA = mrna_AA_FILE_PATH
        INPUT_FILE_SEQ = mrna_SEQ_FILE_PATH

    for (seq_record_1) in SeqIO.parse(INPUT_FILE_AA + "\/region_group_" + str(n) + ".fasta", "fasta"):
        region_header = seq_record_1.id.split('|')
        if SEQ_TYPE == 'Amino_Acid_mRNA_ATG':
            region = seq_record_1.seq[int(region_header[1]) - 1000:int(region_header[1]) + 1000]

        else:
            region = seq_record_1.seq[:2000]

        temp_region = pd.DataFrame(list(region), index=COL)
        temp_region = temp_region.T
        df_region = pd.concat([df_region, temp_region], axis=0)
    df_region = df_region.reset_index(drop=True)

    for (seq_record_2) in SeqIO.parse(INPUT_FILE_SEQ + "\/group_" + str(n) + ".fasta", "fasta"):
        seq_header = seq_record_2.id.split('|')
        if SEQ_TYPE == 'Amino_Acid_mRNA_ATG':
            seq = seq_record_2.seq[int(seq_header[1]) - 1000:int(seq_header[1]) + 1000]

        else:
            seq = seq_record_2.seq[:2000]

        temp_seq = pd.DataFrame(list(seq), index=COL)
        temp_seq = temp_seq.T
        df_seq = pd.concat([df_seq, temp_seq], axis=0)
    df_seq = df_seq.reset_index(drop=True)
    records = len(df_seq)
    aa_count = amino_acid_density(df_region)


    region_matrix = aa_matrix(df_region)
    # print(region_matrix)

    seq_matrix = df_seq.replace(["A", "C", "G", "T", "N", "S"], [0, 1, 1, 0, 0, 1])
    seq_matrix = seq_matrix.replace('[A-Z]', 0, regex=True)

    print('Region dataframe size: {} and Sequence dataframe size: {}'.format(region_matrix.shape, seq_matrix.shape))
    if region_matrix.shape[0] != seq_matrix.shape[0]:
        print('Warning, Region and Seq does not match')

    df_seq_region = region_matrix.multiply(seq_matrix)
    df_seq_region_count = amino_acid_density(df_seq_region)
    # print(df_seq_region_count)
    # total_region_bases = df_region_matrix.sum(axis=0).sum()
    NUC_region_count, NUC_df_seq_region_count = get_nucleosome_gc_content(region_matrix, df_seq_region)


    return aa_count, records, df_seq_region_count, NUC_region_count, NUC_df_seq_region_count



def grid_plotting(data, plot_title, x_title, y_title, x_column, y_column, save_loc):
    data = data.dropna(subset=[y_column])
    fig, axes = plt.subplots(9, 3, sharey=True, sharex=True, figsize=(20, 30))
    fig.subplots_adjust(wspace=0.1, hspace=0.25)
    aa_ind = 0
    for ax in axes.reshape(-1):
        data_ = data[data['amino_acid'] == str(AA_include[aa_ind])]
        sns.scatterplot(data=data_,
                        x=x_column, y=y_column, fc='none',
                        s=4, linewidth=1, ax=ax, color='red', edgecolor="black")
        # handles, labels = ax.get_legend_handles_labels()
        # new_labels = [region_labels.get(item, item) for item in labels]
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        if len(data_) > 0:
            ax.set_xticks(ticks=np.arange(data_[x_column].min(), data_[x_column].max(), 100))
        # ax.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
        ax.set_xlim(0, 1000, 100)
        ax.set(xlabel=None, ylabel=None)
        ax.grid(True, linestyle='--', linewidth=0.5, axis='both')
        ax.set_title('AMINO ACID = ' + str(AA_include[aa_ind]), size=9)
        aa_ind = aa_ind + 1
    plt.suptitle(plot_title, y=0.92)
    plt.setp(axes[-1, :], xlabel=x_title)
    plt.setp(axes[:, 0], ylabel=y_title)
    # plt.show()
    plt.savefig(save_loc)
    return None


def broken_axis_plot(data_, plot_title, x_title, y_title, x_column, y_column, hue_column, save_loc, N):
    colormap_region = {'non-CDS_signal': 'blue', 'CDS_signal': 'red'}
    region_labels = {'non-CDS_signal': "non-CDS", 'CDS_signal': "CDS"}
    data_ = data_.dropna(subset=[y_column])
    fig, (ax_1, ax_2) = plt.subplots(1, 2, sharey=True, figsize=(10, 7), gridspec_kw={'width_ratios': [0.1, 0.9]})
    fig.subplots_adjust(wspace=0.05)
    ax_1.spines.right.set_visible(False)
    ax_2.spines.left.set_visible(False)

    sns.scatterplot(data=data_, x=x_column, y=y_column, hue=hue_column, fc='none',
                    ec=data_[hue_column].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_2)
    handles, labels = ax_2.get_legend_handles_labels()
    new_labels = [region_labels.get(item, item) for item in labels]
    ax_2.tick_params(axis='x', labelsize=8)
    ax_2.tick_params(axis='y', colors='white')
    ax_2.set_xticks(ticks=np.arange(data_[x_column].min(), data_[x_column].max(), 100))
    ax_2.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    ax_2.set_xlim(-100, 1000)
    ax_2.set(xlabel=None, ylabel=None)
    ax_2.grid(True, linestyle='--', linewidth=0.5, axis='both')

    sns.scatterplot(data=data_, x=x_column, y=y_column, hue=hue_column, fc='none',
                    ec=data_[hue_column].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_1, legend=False)
    ax_1.tick_params(axis='both', labelsize=8)
    ax_1.set_xticks(ticks=np.arange(data_[x_column].min(), data_[x_column].max(), 100))
    ax_1.set_xlim(-1000, -850)
    ax_1.set(xlabel=None, ylabel=None)
    ax_1.grid(True, linestyle='--', linewidth=0.5, axis='both')

    fig.supxlabel(x_title, fontsize=9, x=0.5, y=0.05)
    fig.supylabel(y_title, fontsize=9, y=0.5, x=0.04)
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



if __name__ == "__main__":
    start = time.perf_counter()

    aa_density = pd.DataFrame()
    aa_gc_density = pd.DataFrame()
    total_records = 0
    NUC_aa_df = pd.DataFrame()
    NUC_aa_gc_count = pd.DataFrame()


    print('Reading the Sequence and Region files........')
    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):

            temp_aa_density, rec, temp_aa_gc,  temp_nuc_aa_count, temp_nuc_aa_gc_count = j.result()
            aa_density = pd.concat([aa_density, temp_aa_density], axis=0)
            total_records = total_records + rec

            aa_gc_density = pd.concat([aa_gc_density, temp_aa_gc], axis=0)

            NUC_aa_df = pd.concat([NUC_aa_df, temp_nuc_aa_count], axis=0)
            NUC_aa_gc_count = pd.concat([NUC_aa_gc_count, temp_nuc_aa_gc_count], axis=0)

    #
    # NUC_aa_df.to_csv('NUC_aa_df.csv')
    # NUC_aa_gc_count.to_csv('NUC_aa_gc_count.csv')


    print('CALCULATING THE AMINO ACID DENSTIY')

    aa_density.index.name = 'amino_acid'
    aa_density_grouped = aa_density.groupby('amino_acid').sum()
    aa_density_grouped.reset_index(inplace=True, drop=False)

    aa_density_grouped['amino_acid'] = aa_density_grouped['amino_acid'].replace(['z', 'i', 'u', 'd', 'f', 'x'],
                                                                                ['UF', 'Intron', 'UUTR', 'DUTR', 'DF',
                                                                                 'Stop'])
    aa_density_grouped.loc[: ,'total_aa'] = aa_density_grouped.loc[:,COL].sum(axis=1)
    total_aa_acids = aa_density_grouped[['amino_acid', 'total_aa']]


    aa_density_grouped.drop(['total_aa'], inplace=True, axis=1)
    aa_density_grouped[COL] = aa_density_grouped[COL]/total_records
    aa_density_melted = pd.melt(aa_density_grouped, id_vars=['amino_acid'], value_vars=COL)
    aa_density_melted.rename(columns={'value':'aa_density'}, inplace=True)

    aa_density_melted['variable'] = aa_density_melted['variable'].astype(int)
    aa_density_melted.rename(columns={'variable': 'position'}, inplace=True)
    aa_density_melted['position'] = aa_density_melted['position']-1000
    aa_density_melted.to_csv(AMINO_ACID_DENSITY_FILE)

    # ['A', 'C', 'D', 'E', 'F', 'G', 'H',
    # 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
    # 'R', 'S', 'T', 'U', 'V', 'W',
    # 'Y', 'X', 'B', 'Z', 'J', 'O',
    # 'x', 'u', 'i', 'd', 'f', 'z']

    AA_include = ['A', 'C', 'D', 'E', 'F', 'G', 'H',
                'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                'R', 'S', 'T', 'U', 'V', 'W',
                'Y', 'X', 'B', 'Z', 'J', 'O', 'Stop']

    if SEQ_TYPE == 'Amino_Acid_mRNA_ATG':
        x_axis_title = 'position w.r.t TIS'
    else:
        x_axis_title = 'position w.r.t TSS'

    grid_plotting(data = aa_density_melted,
                  plot_title ="Amino acid density" ,
                  x_title=x_axis_title,
                  y_title='density',
                  x_column='position',
                  y_column='aa_density',
                  save_loc=AMINO_ACID_DENSITY_CHART)

    aa_density_melted.sort_values(by=['position'], inplace=True, ascending=True)
    aa_density_melted.reset_index(drop=True, inplace=True)

    nuc_aa_data = aa_density_melted.groupby('amino_acid')['aa_density'].rolling(window=147,
                                                                                        min_periods=1,
                                                                                        center=True).mean().reset_index()
    nuc_aa_data.rename(columns={'level_1': 'index', 'aa_density': 'nuc_aa_density'}, inplace=True)

    aa_density_melted['index'] = aa_density_melted.index
    nuclo_aa_density_melted = aa_density_melted.merge(nuc_aa_data, on=['index', 'amino_acid'], how='left')
    nuclo_aa_density_melted.drop(columns=['index'], inplace=True)

    nuclo_aa_density_melted.to_csv(NUCLO_AMINO_ACID_DENSITY_FILE)

    grid_plotting( data = nuclo_aa_density_melted,
                    plot_title = "Amino acid density",
                    x_title = x_axis_title,
                    y_title = 'density',
                    x_column = 'position',
                    y_column = 'nuc_aa_density',
                    save_loc = NUCLO_AMINO_ACID_DENSITY_CHART)


    print('CALCULATING THE AVG GC CONTENT IN EACH AMINO ACID')

    aa_gc_density.index.name = 'amino_acid'

    aa_gc_density_grouped = aa_gc_density.groupby('amino_acid').sum()
    aa_gc_density_grouped.reset_index(inplace=True, drop=False)

    aa_gc_density_grouped['amino_acid'] = aa_gc_density_grouped['amino_acid'].map(OPP_AA_DICT)
    aa_gc_density_grouped.loc[:, 'gc_count'] = aa_gc_density_grouped.loc[:, COL].sum(axis=1)
    aa_gc_density_grouped['amino_acid'] = aa_gc_density_grouped['amino_acid'].replace(['z', 'i', 'u', 'd', 'f', 'x'],
                                                                               ['UF', 'Intron', 'UUTR', 'DUTR', 'DF',
                                                                                'Stop'])

    aa_gc_density_grouped.dropna(subset=['amino_acid'], axis=0, inplace=True)
    avg_gc_aa_df = aa_gc_density_grouped[['amino_acid', 'gc_count']]
    aa_gc_density_grouped.drop(['gc_count'], axis=1, inplace=True)

    merged_avg_gc_aa = avg_gc_aa_df.merge(total_aa_acids, on='amino_acid', how='left')
    merged_avg_gc_aa['avg_gc'] = merged_avg_gc_aa['gc_count']/merged_avg_gc_aa['total_aa']
    merged_avg_gc_aa['avg_gc'] = merged_avg_gc_aa['avg_gc'].replace([np.inf, -np.inf], np.nan)
    merged_avg_gc_aa.to_csv(AVG_GC_AA_FILE)

    fig, ax3 = plt.subplots(figsize=(10, 7))
    sns.barplot(data=merged_avg_gc_aa[merged_avg_gc_aa['amino_acid'].isin(AA_include)],
                x="amino_acid", y="avg_gc", ax=ax3, color='gray')
    plt.xlabel("amino acid", fontsize=9)
    plt.ylabel("GC content", fontsize=9)
    plt.title("Average GC content per amino acid", size=10, fontweight='bold')
    plt.tick_params(axis='both', labelrotation=0, labelsize=8)
    plt.grid(linestyle='--', linewidth=0.5, axis='y')
    plt.savefig(AVG_GC_AA_CHART)


    print('CALCULATING THE GC CONTENT FROM THE CHOICE OF AMINO ACID FOR EACH AMINO ACID')

    signal_choice_aa = aa_density_melted.merge(merged_avg_gc_aa[['amino_acid', 'avg_gc']],
                                               on='amino_acid',
                                               how='left')
    signal_choice_aa['signal'] = signal_choice_aa['aa_density']*signal_choice_aa['avg_gc']

    signal_choice_aa.to_csv(AA_CHOICE_SIGNAL_FILE)

    grid_plotting(data=signal_choice_aa,
                  plot_title="Average GC content per bp position for each amino acid",
                  x_title=x_axis_title,
                  y_title='GC content',
                  x_column='position',
                  y_column='signal',
                  save_loc=AA_CHOICE_SIGNAL_CHART)


    print('CALCULATING THE TOTAL SIGNAL FROM CHOICE OF AMINO ACID')

    cds_inter_signal = signal_choice_aa[signal_choice_aa['amino_acid'].isin(AA_include)]
    cds_inter_signal_grouped = cds_inter_signal.groupby('position')['signal'].sum()
    cds_inter_signal_grouped = cds_inter_signal_grouped.reset_index()
    cds_inter_signal_grouped.rename(columns={'signal': 'CDS_signal'}, inplace=True)

    not_cds_inter_signal = signal_choice_aa[signal_choice_aa['amino_acid'].isin(['UF', 'Intron', 'UUTR', 'DUTR', 'DF'])]
    not_cds_inter_signal_grouped = not_cds_inter_signal.groupby('position')['signal'].sum()
    not_cds_inter_signal_grouped = not_cds_inter_signal_grouped.reset_index()
    not_cds_inter_signal_grouped.rename(columns={'signal':'non-CDS_signal'}, inplace=True)
    total_inter_signal = not_cds_inter_signal_grouped.merge(cds_inter_signal_grouped, on='position', how='left')
    total_inter_signal = pd.melt(total_inter_signal,
                                 id_vars=['position'],
                                 value_vars=['non-CDS_signal', 'CDS_signal'])

    total_inter_signal.to_csv(TOTAL_AA_CHOICE_SIGNAL_FILE)

    broken_axis_plot(data_=total_inter_signal,
                     plot_title="Total signal from choice of amino acid",
                     x_title=x_axis_title,
                     y_title='GC content',
                     x_column='position',
                     y_column='value',
                     hue_column='variable',
                     save_loc=TOTAL_AA_CHOICE_SIGNAL_CHART, N=False)

    ####### Subtract the inter region signal from the choice of amino acid signal to get the position dependence effect.

    pre_mrna_inter_signal_region_dist = pd.read_csv(CDS_INTER_REGION_FILE)
    pre_mrna_inter_signal_region_dist.drop(columns=['Unnamed: 0', 'index', 'region_density', 'avg_gc'], inplace=True)
    cds_inter_signal_subtracted = cds_inter_signal_grouped.merge(pre_mrna_inter_signal_region_dist[pre_mrna_inter_signal_region_dist['region']=='E'],
                                                                       on=['position'],
                                                                       how='left')
    cds_inter_signal_subtracted['pd_choice_aa_signal'] = cds_inter_signal_subtracted['CDS_signal'] - cds_inter_signal_subtracted['signal']

    cds_inter_signal_subtracted.to_csv(CORRECTED_TOTAL_AA_CHOICE_SIGNAL_FILE)

    # broken_axis_plot(data_=cds_inter_signal_subtracted,
    #                  plot_title="Total signal from choice of amino acid",
    #                  x_title='position w.r.t TSS',
    #                  y_title='GC content',
    #                  x_column='position',
    #                  y_column='pd_choice_aa_signal',
    #                  hue_column='variable',
    #                  save_loc=CORRECTED_TOTAL_AA_CHOICE_SIGNAL_CHART)

    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=cds_inter_signal_subtracted[cds_inter_signal_subtracted['position'] >= -100], x='position', y='pd_choice_aa_signal',
                    fc='none', s=1.5, linewidth=1, edgecolor="black")


    plt.xlabel(x_axis_title, fontsize=9)
    plt.ylabel("GC content", fontsize=9)
    # plt.xticks(ticks=np.arange(-1000, 1000, 100))
    plt.tick_params(axis='x', labelrotation=45, labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    plt.title("Total signal from choice of amino acid-Adjusted", size=10, fontweight='bold')
    plt.xlim(-100, 1000, 100)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=8)
    plt.savefig(CORRECTED_TOTAL_AA_CHOICE_SIGNAL_CHART)


    print('CALCULATING THE GC CONTENT FROM THE SYNONYMOUS CODONS FOR EACH AMINO ACID')

    aa_gc_density_grouped[COL] = aa_gc_density_grouped[COL]/total_records
    aa_gc_density_grouped_melted = pd.melt(aa_gc_density_grouped, id_vars=['amino_acid'], value_vars=COL)
    aa_gc_density_grouped_melted.rename(columns={'value':'syn_signal'}, inplace=True)
    aa_gc_density_grouped_melted.rename(columns={'variable': 'position'}, inplace=True)

    aa_gc_density_grouped_melted['position'] = aa_gc_density_grouped_melted['position'].astype(int)
    aa_gc_density_grouped_melted['position'] = aa_gc_density_grouped_melted['position'] - 1000

    aa_gc_density_grouped_melted.to_csv(SYN_CODON_SIGNAL_FILE)

    # grid_plotting(data=aa_gc_density_grouped_melted,
    #               plot_title="GC content per bp position for each amino acid due to synonymous codon",
    #               x_title='position w.r.t TSS',
    #               y_title='GC content',
    #               x_column='position',
    #               y_column='syn_signal',
    #               save_loc=SYN_CODON_SIGNAL_CHART)
    #


    print('CALCULATING THE MEAN CENTERED SIGNAL FROM THE SYNONYMOUS CODONS FOR EACH AMINO ACID')
    cent_aa_gc_density_grouped_melted = aa_gc_density_grouped_melted[['position',
                                                                      'amino_acid',
                                                                      'syn_signal']].merge(signal_choice_aa[['position',
                                                                                                             'amino_acid',
                                                                                                             'signal']],
                                                                                           on=['position',
                                                                                               'amino_acid'],
                                                                                           how='left')
    cent_aa_gc_density_grouped_melted['corrected_syn_signal'] = cent_aa_gc_density_grouped_melted['syn_signal'] - \
                                                                cent_aa_gc_density_grouped_melted['signal']
    # aa_gc_density_grouped_melted['syn_signal'] = aa_gc_density_grouped_melted.apply(make_nan, value_col='syn_signal', axis=1)
    #
    # syn_mean_data = aa_gc_density_grouped_melted.groupby(['amino_acid'])['syn_signal'].mean()
    # syn_mean_data = syn_mean_data.reset_index()
    #
    # syn_mean_data.rename(columns={'syn_signal': 'mean_syn_signal'}, inplace=True)
    #
    # cent_aa_gc_density_grouped_melted = aa_gc_density_grouped_melted.merge(syn_mean_data, on='amino_acid', how='left')
    # cent_aa_gc_density_grouped_melted['cent_syn_signal'] = cent_aa_gc_density_grouped_melted['syn_signal'] - cent_aa_gc_density_grouped_melted['mean_syn_signal']
    #
    cent_aa_gc_density_grouped_melted.to_csv(CEN_SYN_CODON_SIGNAL_FILE)

    grid_plotting(data=cent_aa_gc_density_grouped_melted,
                  plot_title="Mean Centered GC content per bp position for each amino acid due to synonymous codon",
                  x_title=x_axis_title,
                  y_title='GC content',
                  x_column='position',
                  y_column='corrected_syn_signal',
                  save_loc=CEN_SYN_CODON_SIGNAL_CHART)



    print('CALCULATING THE NORMALISED GC CONTENT FROM THE CHOICE OF SYNONYMOUS CODON')


    aa_gc_density_grouped_melted_joined_aa_den = aa_gc_density_grouped_melted.merge(aa_density_melted,
                                                                                     on=['amino_acid', 'position'],
                                                                                     how='left')
    aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal'] = aa_gc_density_grouped_melted_joined_aa_den['syn_signal']/\
                                                                          aa_gc_density_grouped_melted_joined_aa_den['aa_density']

    aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal'] = aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal'].replace([np.inf, -np.inf], np.nan)
    aa_gc_density_grouped_melted_joined_aa_den.to_csv(NORMALISED_SYN_CODON_SIGNAL_FILE)

    # grid_plotting(data=aa_gc_density_grouped_melted_joined_aa_den,
    #               plot_title="Normalised GC content per bp position for each amino acid due to synonymous codon",
    #               x_title='position w.r.t TSS',
    #               y_title='GC content',
    #               x_column='position',
    #               y_column='normalised_syn_signal',
    #               save_loc=NORMALISED_SYN_CODON_SIGNAL_CHART)


    print('CALCULATING THE MEAN CENTERED NORMALISED GC CONTENT FROM THE CHOICE OF SYNONYMOUS CODON')

    cent_aa_gc_density_grouped_melted_joined_aa_den = aa_gc_density_grouped_melted_joined_aa_den[['position',
                                                                    'amino_acid',
                                                                    'normalised_syn_signal']].merge(merged_avg_gc_aa[['amino_acid', 'avg_gc']],
                                                                                                    on=['amino_acid'],
                                                                                how='left')
    cent_aa_gc_density_grouped_melted_joined_aa_den['cent_norm_intra_signal'] = cent_aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal']- cent_aa_gc_density_grouped_melted_joined_aa_den['avg_gc']
    cent_aa_gc_density_grouped_melted_joined_aa_den.to_csv(CENT_NORMALISED_SYN_CODON_SIGNAL_FILE)

    grid_plotting(data=cent_aa_gc_density_grouped_melted_joined_aa_den,
                  plot_title= "Mean centered normalised GC content per bp position for each amino acid due to synonymous codon",
                  x_title=x_axis_title,
                  y_title='GC content',
                  x_column='position',
                  y_column='cent_norm_intra_signal',
                  save_loc=CENT_NORMALISED_SYN_CODON_SIGNAL_CHART)



    print('CALCULATING THE TOTAL INTRA REGIONAL SIGNAL (Signal from synonymous codons)')


    cds_intra_signal = cent_aa_gc_density_grouped_melted[cent_aa_gc_density_grouped_melted['amino_acid'].isin(AA_include)]
    cds_intra_signal_grouped = cds_intra_signal.groupby('position')['corrected_syn_signal'].sum()
    cds_intra_signal_grouped = cds_intra_signal_grouped.reset_index()
    cds_intra_signal_grouped.rename(columns={'corrected_syn_signal': 'CDS_signal'}, inplace=True)

    not_cds_intra_signal = cent_aa_gc_density_grouped_melted[cent_aa_gc_density_grouped_melted['amino_acid'].isin(['UF', 'Intron', 'UUTR', 'DUTR', 'DF'])]
    not_cds_intra_signal_grouped = not_cds_intra_signal.groupby('position')['corrected_syn_signal'].sum()
    not_cds_intra_signal_grouped = not_cds_intra_signal_grouped.reset_index()
    not_cds_intra_signal_grouped.rename(columns={'corrected_syn_signal': 'non-CDS_signal'}, inplace=True)

    total_intra_signal = not_cds_intra_signal_grouped.merge(cds_intra_signal_grouped, on='position', how='left')
    total_intra_signal = pd.melt(total_intra_signal, id_vars=['position'], value_vars=['non-CDS_signal', 'CDS_signal'])
    total_intra_signal = total_intra_signal.fillna(0)

    total_intra_signal.to_csv(TOTAL_SYN_CODON_SIGNAL_FILE)

    broken_axis_plot(data_=total_intra_signal,
                     plot_title="Total signal from choice of synonymous codons",
                     x_title=x_axis_title,
                     y_title='GC content',
                     x_column='position',
                     y_column='value',
                     hue_column='variable',
                     save_loc=TOTAL_SYN_CODON_SIGNAL_CHART, N=False)


    print('PLOT THE CHOICE OF AA SIGNAL AND SYNONYMOUS CODON SIGNAL')

    aa_signal = cds_inter_signal_subtracted[['position', 'pd_choice_aa_signal']].merge(cds_intra_signal_grouped,
                                      on='position',
                                      how='left')
    aa_signal.rename(columns={'CDS_signal': 'syn_signal'}, inplace=True)

    aa_signal['CDS_intra_signal'] = aa_signal['syn_signal'] + aa_signal['pd_choice_aa_signal']

    aa_signal = pd.melt(aa_signal, id_vars=['position'],
                        value_vars=['pd_choice_aa_signal', 'syn_signal', 'CDS_intra_signal'])

    colormap_region = {'pd_choice_aa_signal': 'blue', 'syn_signal': 'red', 'CDS_intra_signal': 'black'}
    region_labels = {'pd_choice_aa_signal': "aa choice signal", 'syn_signal': "syn codon signal",
                     'CDS_intra_signal': 'CDS Intra-region signal'}

    fig, (ax_1, ax_2) = plt.subplots(1, 2, sharey=True, figsize=(10, 7), gridspec_kw={'width_ratios': [0.1, 0.9]})
    fig.subplots_adjust(wspace=0.05)
    ax_1.spines.right.set_visible(False)
    ax_2.spines.left.set_visible(False)

    sns.scatterplot(data=aa_signal, x='position', y='value', hue='variable', fc='none',
                    ec=aa_signal['variable'].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_2)
    handles, labels = ax_2.get_legend_handles_labels()
    new_labels = [region_labels.get(item, item) for item in labels]
    ax_2.tick_params(axis='x', labelsize=8)
    ax_2.tick_params(axis='y', colors='white')
    ax_2.set_xticks(ticks=np.arange(aa_signal['position'].min(), aa_signal['position'].max(), 100))
    ax_2.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    ax_2.set_xlim(-100, 1000)
    ax_2.set(xlabel=None, ylabel=None)
    ax_2.grid(True, linestyle='--', linewidth=0.5, axis='both')

    sns.scatterplot(data=aa_signal, x='position', y='value', hue='variable', fc='none',
                    ec=aa_signal['variable'].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_1, legend=False)
    ax_1.tick_params(axis='both', labelsize=8)
    ax_1.set_xticks(ticks=np.arange(aa_signal['position'].min(), aa_signal['position'].max(), 100))
    ax_1.set_xlim(-1000, -850)
    ax_1.set(xlabel=None, ylabel=None)
    ax_1.grid(True, linestyle='--', linewidth=0.5, axis='both')

    fig.supxlabel(x_axis_title, fontsize=9, x=0.5, y=0.05)
    fig.supylabel('GC content', fontsize=9, y=0.5, x=0.04)
    fig.suptitle('Position dependent signal from choice of amino acid and synonymous codon',
                 size=10, fontweight='bold', y=0.91, x=0.5)
    d = .01
    kwargs = dict(transform=ax_1.transAxes, color='k', clip_on=False)
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (-d, +d), **kwargs)  ##bottom left
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (1 - d, 1 + d), **kwargs)  ##top left
    kwargs.update(transform=ax_2.transAxes)
    ax_2.plot((-d + .6 * d, +d - .6 * d), (1 - d, 1 + d), **kwargs)  ##top right
    ax_2.plot((-d + .6 * d, +d - .6 * d), (-d, +d), **kwargs)  ##bottom right
    plt.savefig(AA_CHOICE_AND_SYN_SIGNAL)



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

    print('CALCUALTING NUCLEOSOME AVG GC CONTENT FOR EACH AMINO ACID')

    aa_code = NUC_aa_df.columns
    NUC_aa_df['position'] = NUC_aa_df.index
    NUC_aa_gc_count['position'] = NUC_aa_gc_count.index

    NUC_aa_df = NUC_aa_df.groupby(['position']).sum().reset_index()
    NUC_aa_gc_count = NUC_aa_gc_count.groupby(['position']).sum().reset_index()

    nuc_total_gc = NUC_aa_gc_count.sum(axis=0).reset_index()
    nuc_total_gc.rename(columns={0: 'total_gcs', 'index': 'amino_acid'}, inplace=True)

    nuc_total_bases = NUC_aa_df.sum(axis=0).reset_index()
    nuc_total_bases.rename(columns={0: 'total_bases', 'index': 'amino_acid'}, inplace=True)

    nuc_avg_gc_aa = nuc_total_gc.merge(nuc_total_bases, on='amino_acid', how='left')
    nuc_avg_gc_aa['avg_GC'] = nuc_avg_gc_aa['total_gcs'] / nuc_avg_gc_aa['total_bases']

    # nuc_avg_gc_aa = nuc_avg_gc_aa[nuc_avg_gc_aa['amino_acid'] != 'position']
    # nuc_avg_gc_aa['amino_acid'] = nuc_avg_gc_aa['amino_acid'].astype(int)

    nuc_avg_gc_aa['amino_acid'] = nuc_avg_gc_aa['amino_acid'].map(OPP_AA_DICT)
    nuc_avg_gc_aa['amino_acid'] = nuc_avg_gc_aa['amino_acid'].replace(['z', 'i', 'u', 'd', 'f', 'x'],
                                                                        ['UF', 'Intron', 'UUTR', 'DUTR', 'DF', 'Stop'])
    nuc_avg_gc_aa.dropna(inplace=True)

    nuc_avg_gc_aa.to_csv(NUCLO_AVG_GC_AA_FILE)

    fig, ax3 = plt.subplots(figsize=(10, 7))
    sns.barplot(data=nuc_avg_gc_aa[nuc_avg_gc_aa['amino_acid'].isin(AA_include)],
                x="amino_acid", y="avg_GC", ax=ax3, color='gray')
    plt.xlabel("amino acid", fontsize=9)
    plt.ylabel("Nucleosomal GC content", fontsize=9)
    plt.title("Average nucleosomal GC content per amino acid", size=10, fontweight='bold')
    plt.tick_params(axis='both', labelrotation=0, labelsize=8)
    plt.grid(linestyle='--', linewidth=0.5, axis='y')
    plt.savefig(NUCLO_AVG_GC_AA_CHART)

    print('CALCULATING THE NUCLEOSOMAL GC CONTENT FROM THE CHOICE OF AMINO ACID FOR EACH AMINO ACID')

    nuclo_signal_choice_aa = nuclo_aa_density_melted.merge(nuc_avg_gc_aa[['amino_acid', 'avg_GC']],
                                                           on='amino_acid',
                                                           how='left')
    nuclo_signal_choice_aa['signal'] = nuclo_signal_choice_aa['nuc_aa_density'] * nuclo_signal_choice_aa['avg_GC']

    nuclo_signal_choice_aa.to_csv(NUCLO_AA_CHOICE_SIGNAL_FILE)

    grid_plotting(data=nuclo_signal_choice_aa,
                  plot_title="Average nucleosomal GC content per bp position for each amino acid",
                  x_title=x_axis_title,
                  y_title='Nucleosomal GC content',
                  x_column='position',
                  y_column='signal',
                  save_loc=NUCLO_AA_CHOICE_SIGNAL_CHART)



    print('CALCULATING THE TOTAL NUCLEOSOMAL SIGNAL FROM CHOICE OF AMINO ACID')

    nuclo_cds_inter_signal = nuclo_signal_choice_aa[nuclo_signal_choice_aa['amino_acid'].isin(AA_include)]
    nuclo_cds_inter_signal_grouped = nuclo_cds_inter_signal.groupby('position')['signal'].sum()
    nuclo_cds_inter_signal_grouped = nuclo_cds_inter_signal_grouped.reset_index()
    nuclo_cds_inter_signal_grouped.rename(columns={'signal': 'CDS_signal'}, inplace=True)

    nuclo_not_cds_inter_signal = nuclo_signal_choice_aa[
        nuclo_signal_choice_aa['amino_acid'].isin(['UF', 'Intron', 'UUTR', 'DUTR', 'DF'])]
    nuclo_not_cds_inter_signal_grouped = nuclo_not_cds_inter_signal.groupby('position')['signal'].sum()
    nuclo_not_cds_inter_signal_grouped = nuclo_not_cds_inter_signal_grouped.reset_index()
    nuclo_not_cds_inter_signal_grouped.rename(columns={'signal': 'non-CDS_signal'}, inplace=True)
    nuclo_total_inter_signal = nuclo_not_cds_inter_signal_grouped.merge(nuclo_cds_inter_signal_grouped, on='position',
                                                                        how='left')
    nuclo_total_inter_signal = pd.melt(nuclo_total_inter_signal,
                                       id_vars=['position'],
                                       value_vars=['non-CDS_signal', 'CDS_signal'])

    nuclo_total_inter_signal.to_csv(NUCLO_TOTAL_AA_CHOICE_SIGNAL_FILE)

    broken_axis_plot(data_=nuclo_total_inter_signal,
                     plot_title="Total signal from choice of amino acid",
                     x_title=x_axis_title,
                     y_title='Nucleosomal GC content',
                     x_column='position',
                     y_column='value',
                     hue_column='variable',
                     save_loc=NUCLO_TOTAL_AA_CHOICE_SIGNAL_CHART, N=True)

    ####### Subtract the inter region signal from the choice of amino acid signal to get the position dependence effect.

    pre_mrna_NUCLO_inter_signal_region_dist = pd.read_csv(NUCLO_CDS_INTER_REGION_FILE)
    pre_mrna_NUCLO_inter_signal_region_dist.drop(columns=['Unnamed: 0', 'region_density', 'avg_GC'], inplace=True)
    nuclo_cds_inter_signal_subtracted = nuclo_cds_inter_signal_grouped.merge(
        pre_mrna_NUCLO_inter_signal_region_dist[pre_mrna_NUCLO_inter_signal_region_dist['region'] == 'E'],
        on=['position'],
        how='left')
    nuclo_cds_inter_signal_subtracted['pd_choice_aa_signal'] = nuclo_cds_inter_signal_subtracted['CDS_signal'] - \
                                                               nuclo_cds_inter_signal_subtracted['inter_signal']

    nuclo_cds_inter_signal_subtracted.to_csv(NUCLO_CORRECTED_TOTAL_AA_CHOICE_SIGNAL_FILE)


    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=nuclo_cds_inter_signal_subtracted[nuclo_cds_inter_signal_subtracted['position'] >= -100],
                    x='position', y='pd_choice_aa_signal',
                    fc='none', s=1.5, linewidth=1, edgecolor="black")

    plt.xlabel(x_axis_title, fontsize=9)
    plt.ylabel("Nucleosomal GC content", fontsize=9)
    # plt.xticks(ticks=np.arange(-1000, 1000, 100))
    plt.tick_params(axis='x', labelrotation=45, labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    plt.title("Total signal from choice of amino acid-Adjusted", size=10, fontweight='bold')
    plt.xlim(-100, 1000, 100)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=8)
    plt.axvline(73, linestyle='dotted', color='m')
    plt.savefig(NUCLO_CORRECTED_TOTAL_AA_CHOICE_SIGNAL_CHART)



    print('CALCULATING THE GC CONTENT FROM THE SYNONYMOUS CODONS FOR EACH AMINO ACID')

    NUC_aa_gc_count_melted = pd.melt(NUC_aa_gc_count, id_vars='position')
    NUC_aa_gc_count_melted.rename(columns={'variable':'amino_acid'}, inplace=True)
    NUC_aa_gc_count_melted['amino_acid'] = NUC_aa_gc_count_melted['amino_acid'].astype(int)
    NUC_aa_gc_count_melted['amino_acid'] = NUC_aa_gc_count_melted['amino_acid'].map(OPP_AA_DICT)
    NUC_aa_gc_count_melted['amino_acid'] = NUC_aa_gc_count_melted['amino_acid'].replace(['z', 'i', 'u', 'd', 'f', 'x'],
                                                                        ['UF', 'Intron', 'UUTR', 'DUTR', 'DF', 'Stop'])

    NUC_aa_gc_count_melted['value'] = NUC_aa_gc_count_melted['value']/total_records

    NUC_aa_gc_count_melted.rename(columns={'value':'syn_signal'}, inplace=True)
    NUC_aa_gc_count_melted['position'] = NUC_aa_gc_count_melted['position'].astype(int)
    NUC_aa_gc_count_melted['position'] = NUC_aa_gc_count_melted['position'] - 1000

    NUC_aa_gc_count_melted.to_csv(NUCLO_SYN_CODON_SIGNAL_FILE)



    print('CALCULATING THE MEAN CENTERED SIGNAL FROM THE SYNONYMOUS CODONS FOR EACH AMINO ACID')


    # NUC_aa_gc_count_melted['syn_signal'] = NUC_aa_gc_count_melted.apply(make_nan, value_col='syn_signal', axis=1)
    #
    # nuc_syn_mean_data = NUC_aa_gc_count_melted.groupby(['amino_acid'])['syn_signal'].mean()
    # nuc_syn_mean_data = nuc_syn_mean_data.reset_index()
    #
    # nuc_syn_mean_data.rename(columns={'syn_signal': 'mean_syn_signal'}, inplace=True)
    #
    # nuc_cent_aa_gc_density_grouped_melted = NUC_aa_gc_count_melted.merge(nuc_syn_mean_data, on='amino_acid', how='left')
    # nuc_cent_aa_gc_density_grouped_melted['cent_syn_signal'] = nuc_cent_aa_gc_density_grouped_melted['syn_signal'] - \
    #                                                            nuc_cent_aa_gc_density_grouped_melted['mean_syn_signal']

    nuc_cent_aa_gc_density_grouped_melted = NUC_aa_gc_count_melted[['position',
                                                                      'amino_acid',
                                                                      'syn_signal']].merge(nuclo_signal_choice_aa[['position',
                                                                                                             'amino_acid',
                                                                                                             'signal']],
                                                                                           on=['position',
                                                                                               'amino_acid'],
                                                                                           how='left')
    nuc_cent_aa_gc_density_grouped_melted['corrected_syn_signal'] = nuc_cent_aa_gc_density_grouped_melted['syn_signal'] - \
                                                                nuc_cent_aa_gc_density_grouped_melted['signal']




    nuc_cent_aa_gc_density_grouped_melted.to_csv(NUCLO_CEN_SYN_CODON_SIGNAL_FILE)

    grid_plotting(data=nuc_cent_aa_gc_density_grouped_melted,
                  plot_title="Mean Centered GC content per bp position for each amino acid due to synonymous codon",
                  x_title=x_axis_title,
                  y_title='Nucleosomal GC content',
                  x_column='position',
                  y_column='corrected_syn_signal',
                  save_loc=NUCLO_CEN_SYN_CODON_SIGNAL_CHART)




    print('CALCULATING THE NORMALISED GC CONTENT FROM THE CHOICE OF SYNONYMOUS CODON')


    NUC_aa_gc_count_melted_joined_nuc_aa_den = NUC_aa_gc_count_melted.merge(nuclo_aa_density_melted,
                                                                                     on=['amino_acid', 'position'],
                                                                                     how='left')
    NUC_aa_gc_count_melted_joined_nuc_aa_den['normalised_syn_signal'] = NUC_aa_gc_count_melted_joined_nuc_aa_den['syn_signal']/\
                                                                          NUC_aa_gc_count_melted_joined_nuc_aa_den['nuc_aa_density']

    NUC_aa_gc_count_melted_joined_nuc_aa_den['normalised_syn_signal'] = NUC_aa_gc_count_melted_joined_nuc_aa_den['normalised_syn_signal'].replace([np.inf, -np.inf], np.nan)

    NUC_aa_gc_count_melted_joined_nuc_aa_den.to_csv(NUCLO_NORMALISED_SYN_CODON_SIGNAL_FILE)



    print('CALCULATING THE MEAN CENTERED NORMALISED GC CONTENT FROM THE CHOICE OF SYNONYMOUS CODON')

    # NUC_aa_gc_count_melted_joined_nuc_aa_den['normalised_syn_signal'] = NUC_aa_gc_count_melted_joined_nuc_aa_den.apply(
    #     make_nan,
    #     value_col='normalised_syn_signal',
    #     axis=1)
    # nuc_norm_intra_mean_data = NUC_aa_gc_count_melted_joined_nuc_aa_den.groupby(['amino_acid'])[
    #     'normalised_syn_signal'].mean()
    # nuc_norm_intra_mean_data = nuc_norm_intra_mean_data.reset_index()
    # nuc_norm_intra_mean_data.rename(columns={'normalised_syn_signal': 'mean_intra_signal'}, inplace=True)
    #
    # cent_NUC_aa_gc_count_melted_joined_nuc_aa_den = NUC_aa_gc_count_melted_joined_nuc_aa_den.merge(
    #     nuc_norm_intra_mean_data, on='amino_acid',
    #     how='left')
    #
    # cent_NUC_aa_gc_count_melted_joined_nuc_aa_den['cent_norm_intra_signal'] = \
    #     cent_NUC_aa_gc_count_melted_joined_nuc_aa_den['normalised_syn_signal'] - \
    #     cent_NUC_aa_gc_count_melted_joined_nuc_aa_den[
    #         'mean_intra_signal']

    cent_NUC_aa_gc_count_melted_joined_nuc_aa_den = NUC_aa_gc_count_melted_joined_nuc_aa_den[['position',
                                                                                                  'amino_acid',
                                                                                                  'normalised_syn_signal']].merge(
        nuc_avg_gc_aa[['amino_acid', 'avg_GC']],
        on=['amino_acid'],
        how='left')
    cent_NUC_aa_gc_count_melted_joined_nuc_aa_den['cent_norm_intra_signal'] = \
    cent_NUC_aa_gc_count_melted_joined_nuc_aa_den['normalised_syn_signal'] - \
    cent_NUC_aa_gc_count_melted_joined_nuc_aa_den['avg_GC']


    cent_NUC_aa_gc_count_melted_joined_nuc_aa_den.to_csv(NUCLO_CENT_NORMALISED_SYN_CODON_SIGNAL_FILE)

    grid_plotting(data=cent_NUC_aa_gc_count_melted_joined_nuc_aa_den,
                  plot_title="Mean centered normalised GC content per bp position for each amino acid due to synonymous codon",
                  x_title=x_axis_title,
                  y_title='Nucleosomal GC content',
                  x_column='position',
                  y_column='cent_norm_intra_signal',
                  save_loc=NUCLO_CENT_NORMALISED_SYN_CODON_SIGNAL_CHART)



    print('CALCULATING THE TOTAL INTRA REGIONAL SIGNAL (Signal from synonymous codons)')

    NUC_cds_intra_signal = nuc_cent_aa_gc_density_grouped_melted[
        nuc_cent_aa_gc_density_grouped_melted['amino_acid'].isin(AA_include)]
    NUC_cds_intra_signal_grouped = NUC_cds_intra_signal.groupby('position')['corrected_syn_signal'].sum()
    NUC_cds_intra_signal_grouped = NUC_cds_intra_signal_grouped.reset_index()
    NUC_cds_intra_signal_grouped.rename(columns={'corrected_syn_signal': 'CDS_signal'}, inplace=True)

    NUC_not_cds_intra_signal = nuc_cent_aa_gc_density_grouped_melted[
        nuc_cent_aa_gc_density_grouped_melted['amino_acid'].isin(['UF', 'Intron', 'UUTR', 'DUTR', 'DF'])]
    NUC_not_cds_intra_signal_grouped = NUC_not_cds_intra_signal.groupby('position')['corrected_syn_signal'].sum()
    NUC_not_cds_intra_signal_grouped = NUC_not_cds_intra_signal_grouped.reset_index()
    NUC_not_cds_intra_signal_grouped.rename(columns={'corrected_syn_signal': 'non-CDS_signal'}, inplace=True)

    NUC_total_intra_signal = NUC_not_cds_intra_signal_grouped.merge(NUC_cds_intra_signal_grouped, on='position', how='left')
    NUC_total_intra_signal = pd.melt(NUC_total_intra_signal, id_vars=['position'], value_vars=['non-CDS_signal', 'CDS_signal'])
    NUC_total_intra_signal = NUC_total_intra_signal.fillna(0)

    NUC_total_intra_signal.to_csv(NUCLO_TOTAL_SYN_CODON_SIGNAL_FILE)

    broken_axis_plot(data_=NUC_total_intra_signal,
                     plot_title="Total signal from choice of synonymous codons",
                     x_title=x_axis_title,
                     y_title='Nucleosomal GC content',
                     x_column='position',
                     y_column='value',
                     hue_column='variable',
                     save_loc=NUCLO_TOTAL_SYN_CODON_SIGNAL_CHART, N=True)




    print('PLOT THE CHOICE OF AA SIGNAL AND SYNONYMOUS CODON SIGNAL')

    NUCLO_aa_signal = nuclo_cds_inter_signal_subtracted[['position', 'pd_choice_aa_signal']].merge(NUC_cds_intra_signal_grouped,
                                                                                       on='position',
                                                                                       how='left')
    NUCLO_aa_signal.rename(columns={'CDS_signal': 'syn_signal'}, inplace=True)

    NUCLO_aa_signal['CDS_intra_signal'] = NUCLO_aa_signal['pd_choice_aa_signal'] + NUCLO_aa_signal['syn_signal']

    NUCLO_aa_signal = pd.melt(NUCLO_aa_signal, id_vars=['position'],
                        value_vars=['pd_choice_aa_signal', 'syn_signal', 'CDS_intra_signal'])

    colormap_region = {'pd_choice_aa_signal': 'blue', 'syn_signal': 'red', 'CDS_intra_signal':'black'}
    region_labels = {'pd_choice_aa_signal': "aa choice signal", 'syn_signal': "syn codon signal", 'CDS_intra_signal':'CDS Intra-region signal'}

    fig, (ax_1, ax_2) = plt.subplots(1, 2, sharey=True, figsize=(10, 7), gridspec_kw={'width_ratios': [0.1, 0.9]})
    fig.subplots_adjust(wspace=0.05)
    ax_1.spines.right.set_visible(False)
    ax_2.spines.left.set_visible(False)

    sns.scatterplot(data=NUCLO_aa_signal, x='position', y='value', hue='variable', fc='none',
                    ec=NUCLO_aa_signal['variable'].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_2)
    handles, labels = ax_2.get_legend_handles_labels()
    new_labels = [region_labels.get(item, item) for item in labels]
    ax_2.tick_params(axis='x', labelsize=8)
    ax_2.tick_params(axis='y', colors='white')
    ax_2.set_xticks(ticks=np.arange(NUCLO_aa_signal['position'].min(), NUCLO_aa_signal['position'].max(), 100))
    ax_2.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    ax_2.set_xlim(-100, 1000)
    ax_2.set(xlabel=None, ylabel=None)
    ax_2.grid(True, linestyle='--', linewidth=0.5, axis='both')

    sns.scatterplot(data=NUCLO_aa_signal, x='position', y='value', hue='variable', fc='none',
                    ec=NUCLO_aa_signal['variable'].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_1, legend=False)
    ax_1.tick_params(axis='both', labelsize=8)
    ax_1.set_xticks(ticks=np.arange(NUCLO_aa_signal['position'].min(), NUCLO_aa_signal['position'].max(), 100))
    ax_1.set_xlim(-1000, -850)
    ax_1.set(xlabel=None, ylabel=None)
    ax_1.grid(True, linestyle='--', linewidth=0.5, axis='both')

    fig.supxlabel(x_axis_title, fontsize=9, x=0.5, y=0.05)
    fig.supylabel('Nucleosomal GC content', fontsize=9, y=0.5, x=0.04)
    fig.suptitle('Position dependent signal from choice of amino acid and synonymous codon',
                 size=10, fontweight='bold', y=0.91, x=0.5)
    d = .01
    kwargs = dict(transform=ax_1.transAxes, color='k', clip_on=False)
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (-d, +d), **kwargs)  ##bottom left
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (1 - d, 1 + d), **kwargs)  ##top left
    kwargs.update(transform=ax_2.transAxes)
    ax_2.plot((-d + .6 * d, +d - .6 * d), (1 - d, 1 + d), **kwargs)  ##top right
    ax_2.plot((-d + .6 * d, +d - .6 * d), (-d, +d), **kwargs)  ##bottom right
    plt.savefig(NUCLO_AA_CHOICE_AND_SYN_SIGNAL)
