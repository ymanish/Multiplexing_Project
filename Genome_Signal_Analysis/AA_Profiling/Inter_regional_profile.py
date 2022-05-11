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



def make_nan(x):
    if x['amino_acid'] == 'UF' and x['position']>=0:
        return np.nan
    elif x['amino_acid'] != 'UF' and x['position']<0:
        return np.nan
    else: return x['syn_signal']


def main(n):
    df_region = pd.DataFrame()
    df_seq = pd.DataFrame()
    for (seq_record_1) in SeqIO.parse(AA_FILE_PATH + "\/region_group_" + str(n) + ".fasta", "fasta"):
        region_header = seq_record_1.id
        region = seq_record_1.seq[:2000]

        temp_region = pd.DataFrame(list(region), index=COL)
        temp_region = temp_region.T
        df_region = pd.concat([df_region, temp_region], axis=0)
    df_region = df_region.reset_index(drop=True)

    for (seq_record_2) in SeqIO.parse(SEQ_FILE_PATH + "\/group_" + str(n) + ".fasta", "fasta"):
        seq_header = seq_record_2.id
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


    return aa_count, records, df_seq_region_count



def grid_plotting(data, plot_title, x_title, y_title, x_column, y_column, save_loc):
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
    print(save_loc)
    plt.savefig(save_loc)
    return None


def broken_axis_plot(data_, plot_title, x_title, y_title, x_column, y_column, hue_column, save_loc):
    colormap_region = {'non-CDS_signal': 'blue', 'CDS_signal': 'red'}
    region_labels = {'non-CDS_signal': "non-CDS", 'CDS_signal': "CDS"}

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
    fig.supylabel(y_title, fontsize=9, y=0.5, x=0.08)
    fig.suptitle(plot_title, size=10, fontweight='bold', y=0.91, x=0.5)
    d = .01
    kwargs = dict(transform=ax_1.transAxes, color='k', clip_on=False)
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (-d, +d), **kwargs)  ##bottom left
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (1 - d, 1 + d), **kwargs)  ##top left
    kwargs.update(transform=ax_2.transAxes)
    ax_2.plot((-d + .6 * d, +d - .6 * d), (1 - d, 1 + d), **kwargs)  ##top right
    ax_2.plot((-d + .6 * d, +d - .6 * d), (-d, +d), **kwargs)  ##bottom right

    plt.savefig(save_loc)

    return None



if __name__ == "__main__":
    start = time.perf_counter()

    aa_density = pd.DataFrame()
    aa_gc_density = pd.DataFrame()
    total_records = 0

    print('Reading the Sequence and Region files........')
    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for j in concurrent.futures.as_completed(pool):

            temp_aa_density, rec, temp_aa_gc = j.result()
            aa_density = pd.concat([aa_density, temp_aa_density], axis=0)
            total_records = total_records + rec

            aa_gc_density = pd.concat([aa_gc_density, temp_aa_gc], axis=0)

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


    grid_plotting(data = aa_density_melted,
                  plot_title ="Amino acid density" ,
                  x_title='position w.r.t TSS',
                  y_title='density',
                  x_column='position',
                  y_column='aa_density',
                  save_loc=AMINO_ACID_DENSITY_CHART)


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

    ### CALCULATE THE INTERREGIONAL SIGNAL FROM THE AMINO ACID (Signal from choice of amino acid)#########

    signal_choice_aa = aa_density_melted.merge(merged_avg_gc_aa[['amino_acid', 'avg_gc']],
                                               on='amino_acid',
                                               how='left')
    signal_choice_aa['signal'] = signal_choice_aa['aa_density']*signal_choice_aa['avg_gc']

    signal_choice_aa.to_csv(AA_CHOICE_SIGNAL_FILE)

    grid_plotting(data=signal_choice_aa,
                  plot_title="Average GC content per bp position for each amino acid",
                  x_title='position w.r.t TSS',
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
    total_inter_signal = pd.melt(total_inter_signal,id_vars=['position'], value_vars=['non-CDS_signal', 'CDS_signal'])
    # total_inter_signal = total_inter_signal.fillna(0)

    total_inter_signal.to_csv(TOTAL_INTER_SIGNAL_FILE)


    broken_axis_plot(data_=total_inter_signal,
                     plot_title="Total signal from choice of amino acid",
                     x_title='position w.r.t TSS',
                     y_title='GC content',
                     x_column='position',
                     y_column='value',
                     hue_column='variable',
                     save_loc=TOTAL_INTER_SIGNAL_CHART)




    print('CALCULATING THE GC CONTENT FROM THE SYNONYMOUS CODONS FOR EACH AMINO ACID')

    aa_gc_density_grouped[COL] = aa_gc_density_grouped[COL]/total_records
    aa_gc_density_grouped_melted = pd.melt(aa_gc_density_grouped, id_vars=['amino_acid'], value_vars=COL)
    aa_gc_density_grouped_melted.rename(columns={'value':'syn_signal'}, inplace=True)
    aa_gc_density_grouped_melted.rename(columns={'variable': 'position'}, inplace=True)

    aa_gc_density_grouped_melted['position'] = aa_gc_density_grouped_melted['position'].astype(int)
    aa_gc_density_grouped_melted['position'] = aa_gc_density_grouped_melted['position'] - 1000

    aa_gc_density_grouped_melted.to_csv(SYN_CODON_SIGNAL_FILE)

    grid_plotting(data=aa_gc_density_grouped_melted,
                  plot_title="GC content per bp position for each amino acid due to synonymous codon",
                  x_title='position w.r.t TSS',
                  y_title='GC content',
                  x_column='position',
                  y_column='syn_signal',
                  save_loc=SYN_CODON_SIGNAL_CHART)



    print('CALCULATING THE MEAN CENTERED SIGNAL FROM THE SYNONYMOUS CODONS FOR EACH AMINO ACID')

    aa_gc_density_grouped_melted['syn_signal'] = aa_gc_density_grouped_melted.apply(make_nan, axis=1)

    syn_mean_data = aa_gc_density_grouped_melted.groupby(['amino_acid'])['syn_signal'].mean()
    syn_mean_data = syn_mean_data.reset_index()

    syn_mean_data.rename(columns={'syn_signal': 'mean_syn_signal'}, inplace=True)

    cent_aa_gc_density_grouped_melted = aa_gc_density_grouped_melted.merge(syn_mean_data, on='amino_acid', how='left')
    cent_aa_gc_density_grouped_melted['cent_syn_signal'] = cent_aa_gc_density_grouped_melted['syn_signal'] - cent_aa_gc_density_grouped_melted['mean_syn_signal']

    cent_aa_gc_density_grouped_melted.to_csv(CEN_SYN_CODON_SIGNAL_FILE)

    grid_plotting(data=cent_aa_gc_density_grouped_melted,
                  plot_title="Mean Centered GC content per bp position for each amino acid due to synonymous codon",
                  x_title='position w.r.t TSS',
                  y_title='GC content',
                  x_column='position',
                  y_column='cent_syn_signal',
                  save_loc=CEN_SYN_CODON_SIGNAL_CHART)



    print('CALCULATING THE NORMALISED GC CONTENT FROM THE CHOICE OF SYNONYMOUS CODON')


    aa_gc_density_grouped_melted_joined_aa_den = aa_gc_density_grouped_melted.merge(aa_density_melted,
                                                                                     on=['amino_acid', 'position'],
                                                                                     how='left')
    aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal'] = aa_gc_density_grouped_melted_joined_aa_den['syn_signal']/\
                                                                          aa_gc_density_grouped_melted_joined_aa_den['aa_density']

    aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal'].fillna(0, inplace=True)
    aa_gc_density_grouped_melted_joined_aa_den.to_csv(NORMALISED_SYN_CODON_SIGNAL_FILE)

    grid_plotting(data=aa_gc_density_grouped_melted_joined_aa_den,
                  plot_title="Normalised GC content per bp position for each amino acid due to synonymous codon",
                  x_title='position w.r.t TSS',
                  y_title='GC content',
                  x_column='position',
                  y_column='normalised_syn_signal',
                  save_loc=NORMALISED_SYN_CODON_SIGNAL_CHART)


    print('CALCULATING THE MEAN CENTERED NORMALISED GC CONTENT FROM THE CHOICE OF SYNONYMOUS CODON')

    cent_aa_gc_density_grouped_melted_joined_aa_den = cent_aa_gc_density_grouped_melted.merge(aa_density_melted,
                                                                                    on=['amino_acid', 'position'],
                                                                                    how='left')
    cent_aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal'] = cent_aa_gc_density_grouped_melted_joined_aa_den[
                                                                              'cent_syn_signal'] / \
                                                                          cent_aa_gc_density_grouped_melted_joined_aa_den[
                                                                              'aa_density']

    cent_aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal'] = cent_aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal'].replace([np.inf, -np.inf], np.nan)

    cent_aa_gc_density_grouped_melted_joined_aa_den['normalised_syn_signal'].fillna(0, inplace=True)
    cent_aa_gc_density_grouped_melted_joined_aa_den.to_csv(CENT_NORMALISED_SYN_CODON_SIGNAL_FILE)

    grid_plotting(data=cent_aa_gc_density_grouped_melted_joined_aa_den,
                  plot_title= "Mean centered normalised GC content per bp position for each amino acid due to synonymous codon",
                  x_title='position w.r.t TSS',
                  y_title='GC content',
                  x_column='position',
                  y_column='normalised_syn_signal',
                  save_loc=CENT_NORMALISED_SYN_CODON_SIGNAL_CHART)



    print('CALCULATING THE TOTAL INTRA REGIONAL SIGNAL (Signal from synonymous codons)')


    ##### CALCULATE THE TOTAL INTERREGIONAL SIGNAL (Total CDS signal)
    cds_intra_signal = cent_aa_gc_density_grouped_melted[cent_aa_gc_density_grouped_melted['amino_acid'].isin(AA_include)]
    cds_intra_signal_grouped = cds_intra_signal.groupby('position')['cent_syn_signal'].sum()
    cds_intra_signal_grouped = cds_intra_signal_grouped.reset_index()
    cds_intra_signal_grouped.rename(columns={'cent_syn_signal': 'CDS_signal'}, inplace=True)

    not_cds_intra_signal = cent_aa_gc_density_grouped_melted[cent_aa_gc_density_grouped_melted['amino_acid'].isin(['UF', 'Intron', 'UUTR', 'DUTR', 'DF'])]
    not_cds_intra_signal_grouped = not_cds_intra_signal.groupby('position')['cent_syn_signal'].sum()
    not_cds_intra_signal_grouped = not_cds_intra_signal_grouped.reset_index()
    not_cds_intra_signal_grouped.rename(columns={'cent_syn_signal': 'non-CDS_signal'}, inplace=True)

    total_intra_signal = not_cds_intra_signal_grouped.merge(cds_intra_signal_grouped, on='position', how='left')
    total_intra_signal = pd.melt(total_intra_signal, id_vars=['position'], value_vars=['non-CDS_signal', 'CDS_signal'])
    total_intra_signal = total_intra_signal.fillna(0)

    total_intra_signal.to_csv(TOTAL_INTRA_SIGNAL_FILE)

    broken_axis_plot(data_=total_intra_signal,
                     plot_title="Total signal from choice of synonymous codons",
                     x_title='position w.r.t TSS',
                     y_title='GC content',
                     x_column='position',
                     y_column='value',
                     hue_column='variable',
                     save_loc=TOTAL_INTRA_SIGNAL_CHART)