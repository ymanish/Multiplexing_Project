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



    hueorder = ['DUTR', 'exon', 'intron', 'UUTR']
    # fig, (ax1, ax2) = plt.subplots(2, figsize=(16, 9))

    plt.figure(figsize=(16, 9))
    sns.lineplot(x=range(-1000, 1000), y=absolute_gc['0'].values)
    # plt.set(xlabel='Position w.r.t TSS', ylabel='GC content')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("GC content")
    plt.xticks(ticks=np.arange(-1000, 1000, 100))
    plt.tick_params(axis='x', labelrotation=0, labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    plt.title("GC content per base position", size=10)
    # plt.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)
    if NUC:
        plt.savefig(NUCLO_ABSOLUTE_GC_LINE_CHART)
    else:
        plt.savefig(ABSOLUTE_GC_LINE_CHART)


    plt.figure(figsize=(16, 9))
    sns.lineplot(data=element_densities, x='position', y='density', hue="variable", hue_order=hueorder)
    # plt.set(xlabel='Position w.r.t TSS', ylabel='%Occurrence')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("Density")
    plt.tick_params(axis='both', labelsize=8)
    plt.xticks(ticks=np.arange(element_densities['position'].min(), element_densities['position'].max(), 100))
    plt.title("Element Density per base position", size=10)
    plt.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)
    if NUC:
        plt.savefig(NUCLO_ELEMENT_DENSITY_LINE_CHART)
    else:
        plt.savefig(ELEMENT_DENSITY_LINE_CHART)


    INTER_CUM_GC_PER_BP.rename(columns={'avg_GC_region_density' : 'INTER_GC'}, inplace=True)
    INTRA_CUM_GC_PER_BP.rename(columns={'density': 'INTRA_GC'}, inplace=True)

    DERIVED_CUM_GC_PER_BP = INTER_CUM_GC_PER_BP.merge(INTRA_CUM_GC_PER_BP, on='position', how='left')
    DERIVED_CUM_GC_PER_BP['Total GC'] = DERIVED_CUM_GC_PER_BP['INTRA_GC'] + DERIVED_CUM_GC_PER_BP['INTER_GC']
    DERIVED_CUM_GC_PER_BP_MELT = pd.melt(DERIVED_CUM_GC_PER_BP, id_vars=['position'])


    plt.figure(figsize=(16, 9))
    sns.lineplot(data=DERIVED_CUM_GC_PER_BP_MELT, x='position', y='value', hue="variable")
    # plt.set(xlabel='Position w.r.t TSS', ylabel='%Occurrence')
    plt.xlabel("Position w.r.t TSS")
    plt.ylabel("GC content")
    plt.tick_params(axis='both', labelsize=8)
    plt.xticks(ticks=np.arange(DERIVED_CUM_GC_PER_BP_MELT['position'].min(),
                               DERIVED_CUM_GC_PER_BP_MELT['position'].max()+100, 100))
    plt.title("Intra and Inter Regional GC signal per base position", size=10)
    plt.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)
    plt.ylim(DERIVED_CUM_GC_PER_BP_MELT['value'].min(), 0.8, 0.1)
    if NUC:
        plt.savefig(NUCLO_DERIVED_CUM_GC_LINE_CHART)
    else:
        plt.savefig(DERIVED_CUM_GC_LINE_CHART)



    ######Intra_Regional Signals
    fig, ax = plt.subplots(3, 1, figsize=(16, 9))

    sns.lineplot(data=INTRA_GC_CONTENT_PER_REGION_SPECIFIC_BP, x='position', y="density_per_region", hue='variable', hue_order=hueorder, ax=ax[0])
    ax[0].set(xlabel='Position w.r.t TSS', ylabel="GC content")
    ax[0].set_title("GC content distribution per region specific base pair", size=10)
    ax[0].tick_params(axis='both', labelrotation=0, labelsize=8)
    ax[0].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

    sns.lineplot(data=INTRA_GC_CONTENT_PER_REGION_PER_BP, x='position', y="density", hue='variable', hue_order=hueorder, ax=ax[1])
    ax[1].set(xlabel='Position w.r.t TSS', ylabel="GC content")
    ax[1].set_title("GC content distribution per region per base pair", size=10)
    ax[1].tick_params(axis='both', labelrotation=0, labelsize=8)
    ax[1].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

    sns.lineplot(data=INTRA_CUM_GC_PER_BP, x='position', y="INTRA_GC",  ax=ax[2])
    ax[2].set(xlabel='Position w.r.t TSS', ylabel="GC content")
    ax[2].set_title("GC signal per base pair from cumulative Intra Regional signal", size=10)
    ax[2].tick_params(axis='both', labelrotation=0, labelsize=8)
    ax[2].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

    fig.tight_layout()
    if NUC:
        plt.savefig(NUCLO_INTRA_REGION_CHARTS)
    else:
        plt.savefig(INTRA_REGION_CHARTS)


    ######Inter_Regional Signals
    fig, ax = plt.subplots(3, 1, figsize=(16, 9))

    sns.barplot(x="element", y="Avg_GC", data=INTER_AVG_GC_EACH_REGION, hue_order=hueorder, ax=ax[0])
    ax[0].set(xlabel='Region', ylabel="GC content")
    ax[0].set_title("Average GC content per Region", size=10)
    ax[0].tick_params(axis='both', labelrotation=0, labelsize=8)
    ax[0].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

    sns.lineplot(data=INTER_AVG_GC_EACH_REGION_PER_BP, x='position', y="avg_GC_region_density", hue='variable', hue_order=hueorder, ax=ax[1])
    ax[1].set(xlabel='Position w.r.t TSS', ylabel="GC content")
    ax[1].set_title("Average GC content distribution per region per base pair", size=10)
    ax[1].tick_params(axis='both', labelrotation=0, labelsize=8)
    ax[1].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

    sns.lineplot(data=INTER_CUM_GC_PER_BP, x='position', y="INTER_GC",  ax=ax[2])
    ax[2].set(xlabel='Position w.r.t TSS', ylabel="GC content")
    ax[2].set_title("GC signal per base pair from cumulative Inter Regional signal", size=10)
    ax[2].tick_params(axis='both', labelrotation=0, labelsize=8)
    ax[2].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

    fig.tight_layout()
    if NUC:
        plt.savefig(NUCLO_INTER_REGION_CHARTS)
    else:
        plt.savefig(INTER_REGION_CHARTS)

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