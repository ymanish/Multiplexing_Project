import pandas as pd
import numpy as np
import seaborn as sns
import time
import os

import matplotlib.pyplot as plt
from Genome_Signal_Analysis.Initialise_Script import *

# Results_path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\osativa"

absolute_gc = pd.read_csv(ABSOLUTE_GC_FILE)
absolute_gc.drop(['Unnamed: 0'], axis=1, inplace=True)

element_densities = pd.read_csv(ELEMENT_DENSITY_FILE)
element_densities.drop(['Unnamed: 0'], axis=1, inplace=True)

Intra_Regional_Chart8 = pd.read_csv(CENTERED_GC_PER_REGION_PER_BP_FILE)
Intra_Regional_Chart8.drop(['Unnamed: 0'], axis=1, inplace=True)

Intra_Regional_Chart9 = pd.read_csv(CENTERED_GC_DENSITY_BY_REGION_FILE)
Intra_Regional_Chart9.drop(['Unnamed: 0'], axis=1, inplace=True)

Intra_Regional_Chart10 = pd.read_csv(DERIVED_CENTERED_GC_FILE)
Intra_Regional_Chart10.drop(['Unnamed: 0'], axis=1, inplace=True)

Inter_Regional_Chart2 = pd.read_csv(AVG_GC_PER_REGION_FILE)
Inter_Regional_Chart2.drop(['Unnamed: 0'], axis=1, inplace=True)
print (Inter_Regional_Chart2)

Inter_Regional_Chart5 = pd.read_csv(AVG_GC_DENSITY_BY_REGION_FILE)
Inter_Regional_Chart5.drop(['Unnamed: 0'], axis=1, inplace=True)

Inter_Regional_Chart7 = pd.read_csv(DERIVED_AVG_GC_FILE)
Inter_Regional_Chart7.drop(['Unnamed: 0'], axis=1, inplace=True)


hueorder = ['DUTR', 'exon', 'intron', 'UUTR']
fig, (ax1, ax2) = plt.subplots(2, figsize=(16, 9))

sns.lineplot(x=range(-1000, 1000), y=absolute_gc['0'].values, ax=ax1)
ax1.set(xlabel='Position w.r.t TSS', ylabel='%GC content')
ax1.set_xticks(ticks=np.arange(-1000, 1000, 200))
ax1.tick_params(axis='x', labelrotation=0, labelsize=8)
ax1.tick_params(axis='y', labelsize=8)
ax1.set_title("GC content per Base position", size=10)
ax1.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

sns.lineplot(data=element_densities, x='position', y='density', hue="variable", hue_order=hueorder, ax=ax2)
ax2.set(xlabel='Position w.r.t TSS', ylabel='%Occurrence')
ax2.tick_params(axis='both', labelsize=8)
ax2.set_title("Element Density per Base position", size=10)
ax2.legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

fig.tight_layout()
plt.savefig(ABSOLUTE_GC_AND_REGION_DENSITY_CHARTS)


######Intra_Regional Signals
fig, ax = plt.subplots(3, 2, figsize=(20, 10))

sns.lineplot(data=Intra_Regional_Chart8, x='position', y="intra_signal", hue='element', hue_order=hueorder, ax=ax[0][0])
ax[0][0].set(xlabel='Position w.r.t TSS', ylabel="mean-centered %GC")
ax[0][0].set_title("Mean centered %GC content distribution by Region per base position, with average across Region specific Transcripts", size=10)
ax[0][0].tick_params(axis='both', labelrotation=0, labelsize=8)
ax[0][0].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

sns.lineplot(data=Intra_Regional_Chart9, x='position', y="intra_signal", hue='element', hue_order=hueorder, ax=ax[1][0])
ax[1][0].set(xlabel='Position w.r.t TSS', ylabel="mean-centered %GC")
ax[1][0].set_title("Mean centered %GC distribution by Region per base position, with Average across ALL the Transcripts", size=10)
ax[1][0].tick_params(axis='both', labelrotation=0, labelsize=8)
ax[1][0].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

sns.lineplot(data=Intra_Regional_Chart10, x='position', y="intra_signal",  ax=ax[2][0])
ax[2][0].set(xlabel='Position w.r.t TSS', ylabel="mean-centered %GC")
ax[2][0].set_title("Mean Centered Absolute %GC content averaged across ALL transcripts", size=10)
ax[2][0].tick_params(axis='both', labelrotation=0, labelsize=8)
ax[2][0].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

sns.barplot(x="element", y="Avg_GC", data=Inter_Regional_Chart2, hue_order=hueorder, ax=ax[0][1])
ax[0][1].set(xlabel='Region', ylabel="%GC content")
ax[0][1].set_title("Average GC content per Region", size=10)
ax[0][1].tick_params(axis='both', labelrotation=0, labelsize=8)

# ax[0][1].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

sns.lineplot(data=Inter_Regional_Chart5, x='position', y="Homoginized_density_GC", hue='variable', hue_order=hueorder, ax=ax[1][1])
ax[1][1].set(xlabel='Position w.r.t TSS', ylabel="%GC content * %Occurrence")
ax[1][1].set_title("Average %GC distribution by Region per base position", size=10)
ax[1][1].tick_params(axis='both', labelrotation=0, labelsize=8)
ax[1][1].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)

sns.lineplot(data=Inter_Regional_Chart7, x='position', y="Homoginized_density_GC", ax=ax[2][1])
ax[2][1].set(xlabel='Position w.r.t TSS', ylabel="Sum(%GC content * %Occurrence)")
ax[2][1].set_title("Average %GC distribution per base position", size=10)
ax[2][1].tick_params(axis='both', labelrotation=0, labelsize=8)
ax[2][1].legend(handletextpad=0, loc='upper right', markerscale=0.5, fontsize=7)


fig.tight_layout()
plt.savefig(INTRA_AND_INTER_REGION_CHARTS)

# plt.show()

# os.remove(CENTERED_GC_PER_REGION_PER_BP_FILE)
# os.remove(CENTERED_GC_DENSITY_BY_REGION_FILE)
# os.remove(DERIVED_CENTERED_GC_FILE)
#
# os.remove(AVG_GC_DENSITY_BY_REGION_FILE)
# os.remove(DERIVED_AVG_GC_FILE)


