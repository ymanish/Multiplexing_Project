import pandas as pd
import numpy as np
import seaborn as sns
import time
import os
import matplotlib.pyplot as plt
from Initialise_GC_profiling import *


def main(File1, File2, File3, File4, File5, File6, File7, File8, NUC):
    print ('Reading Nucleosomal GC content Files and Generate the charts')


    absolute_gc = pd.read_csv(File1)
    absolute_gc.drop(['Unnamed: 0'], axis=1, inplace=True)

    element_densities = pd.read_csv(File2)
    element_densities.drop(['Unnamed: 0'], axis=1, inplace=True)

    print('Reading INTRA REGIONAL SIGNAL Files......')

    INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP = pd.read_csv(File3)
    INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP.drop(['Unnamed: 0'], axis=1, inplace=True)

    INTRA_GC_CONTENT_PER_REGION_PER_BP = pd.read_csv(File4)
    INTRA_GC_CONTENT_PER_REGION_PER_BP.drop(['Unnamed: 0'], axis=1, inplace=True)

    INTRA_CUM_GC_PER_BP = pd.read_csv(File5)
    INTRA_CUM_GC_PER_BP.drop(['Unnamed: 0'], axis=1, inplace=True)

    print('Reading INTER REGIONAL SIGNAL Files......')

    INTER_AVG_GC_EACH_REGION = pd.read_csv(File6)
    INTER_AVG_GC_EACH_REGION.drop(['Unnamed: 0'], axis=1, inplace=True)

    INTER_AVG_GC_EACH_REGION_PER_BP = pd.read_csv(File7)
    INTER_AVG_GC_EACH_REGION_PER_BP.drop(['Unnamed: 0'], axis=1, inplace=True)

    INTER_CUM_GC_PER_BP = pd.read_csv(File8)
    INTER_CUM_GC_PER_BP.drop(['Unnamed: 0'], axis=1, inplace=True)

    ######ELEMENT DENSITY CHART######
    # hueorder = ['DUTR', 'exon', 'intron', 'UUTR']
    colormap_region = {'DUTR': 'blue', 'exon': 'red', 'intron': 'orange', 'UUTR': 'green', 'UF':'cyan', 'DF':'gray'}
    region_labels = {'DUTR': "3'UTR", 'exon': "CDS", 'intron': "Intron", 'UUTR': "5'UTR", 'UF':"5'end", 'DF':"3'end"}
    element_densities =  element_densities.fillna(0)
    # element_densities['position'] = element_densities['position']/1000
    # fig, ax = plt.subplots(figsize=(10, 7))
    fig, (ax_1, ax_2) = plt.subplots(1, 2, sharey=True, figsize=(10, 7), gridspec_kw={'width_ratios': [0.1, 0.9]})
    fig.subplots_adjust(wspace=0.05)
    ax_1.spines.right.set_visible(False)
    ax_2.spines.left.set_visible(False)

    sns.scatterplot(data=element_densities, x='position', y='density', hue="variable", fc='none',
                    ec=element_densities["variable"].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_2)
    handles, labels = ax_2.get_legend_handles_labels()
    new_labels = [region_labels.get(item, item) for item in labels]
    ax_2.tick_params(axis='x', labelsize=8)
    ax_2.tick_params(axis='y', colors='white')
    ax_2.set_xticks(ticks=np.arange(element_densities['position'].min(), element_densities['position'].max(), 100))
    ax_2.legend(handles= handles, loc='upper right', fontsize=8, labels=new_labels)
    ax_2.set_xlim(-100, 1000)
    ax_2.set(xlabel=None, ylabel=None)
    ax_2.grid(True, linestyle='--', linewidth=0.5, axis='both')

    sns.scatterplot(data=element_densities, x='position', y='density', hue="variable", fc='none',
                    ec=element_densities["variable"].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_1, legend=False)
    ax_1.tick_params(axis='both', labelsize=8)
    ax_1.set_xticks(ticks=np.arange(element_densities['position'].min(), element_densities['position'].max(), 100))
    ax_1.set_xlim(-1000, -850)
    ax_1.set(xlabel=None, ylabel=None)
    ax_1.grid(True,linestyle='--', linewidth=0.5, axis='both')

    fig.supxlabel('position w.r.t TSS',  fontsize=9, x=0.5, y=0.05)
    fig.supylabel('density',  fontsize=9, y=0.5, x=0.08)
    fig.suptitle("Region density per base-pair position", size=10, fontweight='bold', y=0.91, x=0.5)
    d = .01
    kwargs = dict(transform=ax_1.transAxes, color='k', clip_on=False)
    ax_1.plot((1 - d - 3*d, 1 + d + 3*d), (-d, +d), **kwargs) ##bottom left
    ax_1.plot((1 - d - 3*d, 1 + d + 3*d), (1 - d, 1 + d), **kwargs) ##top left
    kwargs.update(transform=ax_2.transAxes)
    ax_2.plot((-d+.6*d, +d-.6*d), (1 - d, 1 + d), **kwargs) ##top right
    ax_2.plot((-d+.6*d, +d-.6*d), (-d, +d), **kwargs) ##bottom right
    ####Old part of the code
    # plt.xlabel("Position w.r.t TSS (in kb)", fontweight='bold')
    # plt.ylabel("Density", fontweight='bold')
    # plt.tick_params(axis='both', labelsize=10)
    # plt.xticks(ticks=np.arange(element_densities['position'].min(), element_densities['position'].max(), 100))
    # plt.title("Region density per base-pair position", size=10, fontweight='bold')
    # plt.legend(handles= handles, loc='upper right', fontsize=8, labels=new_labels, title='Region')
    # plt.ylim(element_densities['density'].min()-element_densities['density'].min()/10, element_densities['density'].max()+element_densities['density'].max()/10)
    # plt.xlim(-100, 1000)
    # plt.grid(linestyle='--', linewidth=0.5)
    print(ELEMENT_DENSITY_LINE_CHART)
    if NUC:
        ax_2.axvline(73, linestyle='dotted', color='m')
        plt.savefig(NUCLO_ELEMENT_DENSITY_LINE_CHART)
    else:
        plt.savefig(ELEMENT_DENSITY_LINE_CHART)




    #####ABSOLUTE GC AND SIGNAL CONTRIBUTION CHART##########
    INTER_CUM_GC_PER_BP.rename(columns={'avg_GC_region_density': 'Inter GC'}, inplace=True)
    INTRA_CUM_GC_PER_BP.rename(columns={'density': 'Intra GC'}, inplace=True)

    DERIVED_CUM_GC_PER_BP = INTER_CUM_GC_PER_BP.merge(INTRA_CUM_GC_PER_BP, on='position', how='left')
    DERIVED_CUM_GC_PER_BP['Derived GC'] = DERIVED_CUM_GC_PER_BP['Intra GC'] + DERIVED_CUM_GC_PER_BP['Inter GC']
    DERIVED_CUM_GC_PER_BP_MELT = pd.melt(DERIVED_CUM_GC_PER_BP, id_vars=['position'])

    absolute_gc['variable'] = 'Actual GC'
    absolute_gc.rename(columns={'0':'value'}, inplace=True)
    absolute_gc['position'] = absolute_gc.index
    absolute_gc['position'] = absolute_gc['position'] - 1000
    Signal_DF = pd.concat([DERIVED_CUM_GC_PER_BP_MELT, absolute_gc], axis=0)
    Signal_DF.reset_index(inplace=True, drop=True)


    plt.figure(figsize=(10, 7))
    colormap_signals = {'Inter GC': 'red', 'Intra GC': 'blue', 'Derived GC': 'green', 'Actual GC':'black'}
    Signal_DF = Signal_DF.fillna(0)
    Signal_DF = Signal_DF[Signal_DF['variable'].isin(['Inter GC', 'Intra GC', 'Actual GC'])]
    sns.scatterplot(data=Signal_DF, x='position', y='value', hue='variable', fc='none',
                    ec=Signal_DF["variable"].map(colormap_signals),
                    palette=colormap_signals, s=1.5, linewidth=1)
    plt.xlabel("position w.r.t TSS", fontsize=9)
    plt.ylabel("GC content", fontsize=9)
    plt.xticks(ticks=np.arange(-1000, 1000, 100))
    plt.tick_params(axis='x', labelrotation=45, labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    plt.title("GC content per base-pair position", size=10, fontweight='bold')
    plt.xlim(-1000, 1000)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=8)

    if NUC:
        plt.axvline(73, linestyle='dotted', color='m')
        plt.savefig(NUCLO_ABSOLUTE_GC_LINE_CHART)
    else:
        plt.savefig(ABSOLUTE_GC_LINE_CHART)


    #####NORMALISED INTRA REGIONAL SIGNAL##########

    fig, ax1 = plt.subplots(figsize=(10, 7))
    y_intra_min_lim = INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP.loc[INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP['position'] >=-1000, 'density_per_region'].min()
    y_intra_max_lim = INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP.loc[INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP['position'] >=-1000, 'density_per_region'].max()

    INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP = INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP.fillna(0)

    sns.scatterplot(data=INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP, x='position', y="density_per_region", hue='variable',
                    fc='none', ec=INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP["variable"].map(colormap_region),
                    palette=colormap_region, s=1.5, linewidth=1, ax=ax1)
    handles, labels = ax1.get_legend_handles_labels()
    new_labels = [region_labels.get(item, item) for item in labels]
    plt.xlabel("position w.r.t TSS", fontsize=9)
    plt.ylabel("GC content", fontsize=9)
    plt.title("Region wise GC content distribution per region-specific base-pair position", size=10, fontweight='bold')
    plt.xticks(np.arange(-1000, 1000, 100))
    plt.tick_params(axis='both', labelrotation=0, labelsize=8)
    plt.legend(handles= handles, loc='upper right', fontsize=8, labels=new_labels)
    # plt.ylim(y_intra_min_lim-(y_intra_min_lim/10), y_intra_max_lim+(y_intra_max_lim/10), 0.1)
    plt.xlim(-1000, 1000, 100)
    plt.grid(linestyle='--', linewidth=0.5)

    if NUC:
        plt.axvline(73, linestyle='dotted', color='m')
        plt.savefig(NUCLO_INTRA_SIGNAL_REGION_DIST_PERBP_CHART)
    else:
        plt.savefig(INTRA_SIGNAL_REGION_DIST_PERBP_CHART)

    #####INTRA REGIONAL SIGNAL##########

    INTRA_GC_CONTENT_PER_REGION_PER_BP = INTRA_GC_CONTENT_PER_REGION_PER_BP.fillna(0)
    y_intra_min_lim_1 = INTRA_GC_CONTENT_PER_REGION_PER_BP.loc[
        INTRA_GC_CONTENT_PER_REGION_PER_BP['position'] >= -1000, 'density'].min()
    y_intra_max_lim_1 = INTRA_GC_CONTENT_PER_REGION_PER_BP.loc[
        INTRA_GC_CONTENT_PER_REGION_PER_BP['position'] >= -1000, 'density'].max()

    fig, ax2 = plt.subplots(figsize=(10, 7))
    sns.scatterplot(data=INTRA_GC_CONTENT_PER_REGION_PER_BP, x='position', y="density", hue='variable', fc='none',
                 ec=INTRA_GC_CONTENT_PER_REGION_PER_BP["variable"].map(colormap_region),
                    palette=colormap_region, s=1.5, linewidth=1, ax=ax2)
    handles, labels = ax2.get_legend_handles_labels()
    new_labels = [region_labels.get(item, item) for item in labels]
    # ax2.set(xlabel='Position w.r.t TSS', ylabel="GC content")
    ax2.set_xlabel("position w.r.t TSS", fontsize=9)
    ax2.set_ylabel("GC content", fontsize=9)
    ax2.set_title("Region wise GC content distribution per base-pair position", size=10, fontweight='bold')
    ax2.set_xticks(np.arange(-1000, 1000, 100))
    ax2.tick_params(axis='both', labelrotation=0, labelsize=8)
    ax2.legend(handles= handles, loc='upper right', fontsize=8, labels=new_labels)
    # ax2.set_ylim(y_intra_min_lim_1 - (y_intra_min_lim_1 / 10), y_intra_max_lim_1 + (y_intra_max_lim_1 / 10))
    ax2.set_xlim(-1000, 1000, 100)
    ax2.grid(linestyle='--', linewidth=0.5)

    if NUC:
        plt.axvline(73, linestyle='dotted', color='m')
        plt.savefig(NUCLO_INTRA_SIGNAL_REGION_DIST_CHART)
    else:
        plt.savefig(INTRA_SIGNAL_REGION_DIST_CHART)



    ######AVG GC CONTENT PER REGION#########
    # y_inter_min_lim = INTER_AVG_GC_EACH_REGION['Avg_GC'].min()
    # y_inter_max_lim = INTER_AVG_GC_EACH_REGION['Avg_GC'].max()
    fig, ax3 = plt.subplots(figsize=(10, 7))
    colormap_avg_region = {"3'UTR": 'blue', "CDS": 'red', "Intron": 'orange', "5'UTR": 'green', "5'end":'cyan', "3'end":'gray'}
    INTER_AVG_GC_EACH_REGION = INTER_AVG_GC_EACH_REGION.fillna(0)
    INTER_AVG_GC_EACH_REGION.replace(region_labels, inplace=True)
    sns.barplot(data=INTER_AVG_GC_EACH_REGION, x="element", y="Avg_GC", palette=colormap_avg_region, ax=ax3)
    plt.xlabel("region", fontsize=9)
    plt.ylabel("GC content", fontsize=9)
    plt.title("Average GC content per Region", size=10, fontweight='bold')
    plt.tick_params(axis='both', labelrotation=0, labelsize=8)
    plt.grid(linestyle='--', linewidth=0.5, axis='y')
    if NUC:
        plt.savefig(NUCLO_AVG_GC_REGION_CHART)
    else:
        plt.savefig(AVG_GC_REGION_CHART)

    #####INTER REGIONAL SIGNAL##########
    INTER_AVG_GC_EACH_REGION_PER_BP = INTER_AVG_GC_EACH_REGION_PER_BP.fillna(0)

    fig, (ax_1, ax_2) = plt.subplots(1, 2, sharey=True, figsize=(10, 7), gridspec_kw={'width_ratios': [0.1, 0.9]})
    fig.subplots_adjust(wspace=0.05)
    ax_1.spines.right.set_visible(False)
    ax_2.spines.left.set_visible(False)
    sns.scatterplot(data=INTER_AVG_GC_EACH_REGION_PER_BP, x='position', y="avg_GC_region_density", hue='variable',
                    fc='none',
                    ec=INTER_AVG_GC_EACH_REGION_PER_BP["variable"].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_2)

    handles, labels = ax_2.get_legend_handles_labels()
    new_labels = [region_labels.get(item, item) for item in labels]
    ax_2.tick_params(axis='x', labelsize=8)
    ax_2.tick_params(axis='y', colors='white')
    ax_2.set_xticks(np.arange(-1000, 1000, 100))
    ax_2.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels)
    ax_2.set_xlim(-100, 1000)
    ax_2.set(xlabel=None, ylabel=None)
    ax_2.grid(True, linestyle='--', linewidth=0.5, axis='both')

    sns.scatterplot(data=INTER_AVG_GC_EACH_REGION_PER_BP, x='position', y="avg_GC_region_density", hue='variable',
                    fc='none',
                    ec=INTER_AVG_GC_EACH_REGION_PER_BP["variable"].map(colormap_region),
                    palette=colormap_region, s=1, linewidth=0.5, ax=ax_1, legend=False)
    ax_1.tick_params(axis='both', labelsize=8)
    ax_1.set_xticks(np.arange(-1000, 1000, 100))
    ax_1.set_xlim(-1000, -850)
    ax_1.set(xlabel=None, ylabel=None)
    ax_1.grid(True, linestyle='--', linewidth=0.5, axis='both')

    fig.supxlabel('position w.r.t TSS', fontsize=9, x=0.5, y=0.05)
    fig.supylabel('GC content', fontsize=9, y=0.5, x=0.08)
    fig.suptitle("Average GC content distribution per region per base-pair position", size=10, fontweight='bold', y=0.91, x=0.5)
    d = .01
    kwargs = dict(transform=ax_1.transAxes, color='k', clip_on=False)
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (-d, +d), **kwargs)  ##bottom left
    ax_1.plot((1 - d - 3 * d, 1 + d + 3 * d), (1 - d, 1 + d), **kwargs)  ##top left
    kwargs.update(transform=ax_2.transAxes)
    ax_2.plot((-d + .6 * d, +d - .6 * d), (1 - d, 1 + d), **kwargs)  ##top right
    ax_2.plot((-d + .6 * d, +d - .6 * d), (-d, +d), **kwargs)  ##bottom right

    if NUC:
        ax_2.axvline(73, linestyle='dotted', color='m')
        plt.savefig(NUCLO_INTER_SIGNAL_REGION_DIST_CHART)
    else:
        plt.savefig(INTER_SIGNAL_REGION_DIST_CHART)


    ######Old part of the code
    # INTER_AVG_GC_EACH_REGION_PER_BP = INTER_AVG_GC_EACH_REGION_PER_BP.fillna(0)
    # fig, ax4 = plt.subplots(figsize=(10, 7))
    # sns.scatterplot(data=INTER_AVG_GC_EACH_REGION_PER_BP, x='position', y="avg_GC_region_density", hue='variable',
    #                 fc='none',
    #                 ec=INTER_AVG_GC_EACH_REGION_PER_BP["variable"].map(colormap_region),
    #                 palette=colormap_region, s=1, linewidth=0.5, ax=ax4)
    # handles, labels = ax4.get_legend_handles_labels()
    # new_labels = [region_labels.get(item, item) for item in labels]
    # ax4.set_xlabel("Position w.r.t TSS", fontweight='bold')
    # ax4.set_ylabel("GC content", fontweight='bold')
    # ax4.set_title("Average GC content distribution per region per base-pair position", size=10, fontweight='bold')
    # ax4.set_xticks(np.arange(-1000, 1000, 100))
    # ax4.tick_params(axis='both', labelrotation=0, labelsize=10)
    # ax4.legend(handles=handles, loc='upper right', fontsize=8, labels=new_labels, title='Region')
    # ax4.set_xlim(-100, 1000, 100)
    # # ax4.set_ylim(Signal_DF.loc[Signal_DF['variable'] == 'Inter GC', 'value'].min() - (Signal_DF.loc[Signal_DF['variable'] == 'Inter GC', 'value'].min() / 10),
    # #              Signal_DF.loc[Signal_DF['variable'] == 'Inter GC', 'value'].max() + (Signal_DF.loc[Signal_DF['variable'] == 'Inter GC', 'value'].max() / 10))
    # ax4.grid(linestyle='--', linewidth=0.5)
    #
    # if NUC:
    #     plt.axvline(73, linestyle='dotted', color='m')
    #     plt.savefig(NUCLO_INTER_SIGNAL_REGION_DIST_CHART)
    # else:
    #     plt.savefig(INTER_SIGNAL_REGION_DIST_CHART)
    return None


main(File1 = ABSOLUTE_GC_FILE,
    File2 = ELEMENT_DENSITY_FILE,
    File3 = CENTERED_GC_PER_REGION_PER_BP_FILE,
    File4 = CENTERED_GC_DENSITY_BY_REGION_FILE,
    File5 = DERIVED_CENTERED_GC_FILE,
    File6 = AVG_GC_PER_REGION_FILE,
    File7 = AVG_GC_DENSITY_BY_REGION_FILE,
    File8 = DERIVED_CUM_AVG_GC_FILE, NUC=False)


main(File1 = NUCLO_ABSOLUTE_GC_FILE,
    File2 = NUCLO_ELEMENT_DENSITY_FILE,
    File3 = NUCLO_CENTERED_GC_PER_REGION_PER_BP_FILE,
    File4 = NUCLO_CENTERED_GC_DENSITY_BY_REGION_FILE,
    File5 = NUCLO_DERIVED_CENTERED_GC_FILE,
    File6 = AVG_GC_PER_REGION_FILE,
    File7 = NUCLO_AVG_GC_DENSITY_BY_REGION_FILE,
    File8 = NUCLO_DERIVED_CUM_AVG_GC_FILE, NUC=True)

