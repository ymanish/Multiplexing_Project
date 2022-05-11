import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 100)
import numpy as np
import scipy.stats

GROUP = 'Yeast' #'Vertebrate' or 'Plant' or 'Yeast'

path = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" + GROUP + '\/'
org_data = pd.read_csv('all_org_list.txt', header=None)
# org_data = pd.read_csv('sample_org_list.txt', header=None)

# not_org_data = pd.read_csv('run_again_org.txt', header=None)
# print(not_org_data)
print(org_data)
SEQ_TYPE = 'pre_mRNA'

intra_region_per_bp = 'Centered_GC_Region_density.csv'
intra_region = 'Centered_GC_density_and element_density.csv'
inter_region = 'Avg_GC_per_element.csv'


intra_signal_perbp_df = pd.DataFrame()
intra_signal_df = pd.DataFrame()
inter_signal_df = pd.DataFrame()

exclude_list = []
# exclude_list = not_org_data[1].tolist()
exclude_list.append('pleo_gene_ensembl')
exclude_list.append('celegans_gene_ensembl')
exclude_list.append('scerevisiae_gene_ensembl')
exclude_list.append('dmelanogaster_gene_ensembl')

def signal_size(DF, colu, org):

    mean_ber_bp = DF.groupby('variable')[colu].mean().reset_index()
    mean_ber_bp = mean_ber_bp.rename(columns={colu: 'mean'})
    DF = DF.merge(mean_ber_bp, on=['variable'], how='left')
    DF['signal'] = DF[colu] - DF['mean']
    DF['signal'] = DF['signal'].abs()
    signal = DF.groupby('variable')['signal'].sum().reset_index()
    signal['organism'] = org
    return signal

for i in org_data[1]:
    if i not in exclude_list:
        result_path = path + i + r'\/' + SEQ_TYPE + r'\Files\/'
        DF_INTRA_SIGNAL = pd.read_csv(result_path+intra_region)
        DF_INTRA_PER_BP_SIGNAL = pd.read_csv(result_path+intra_region_per_bp)
        DF_INTER_SIGNAL = pd.read_csv(result_path+inter_region)
        print(DF_INTER_SIGNAL)
        DF_INTER_SIGNAL['organism'] = i

        signal_perbp_temp = signal_size(DF_INTRA_PER_BP_SIGNAL[DF_INTRA_PER_BP_SIGNAL['position']>=0], 'density_per_region', i)
        signal_temp = signal_size(DF_INTRA_SIGNAL[DF_INTRA_SIGNAL['position']>=0], 'density', i)
        # inter_signal_temp = signal_size(DF_INTER_SIGNAL[DF_INTER_SIGNAL['position']>=0], 'avg_GC_region_density')

        intra_signal_perbp_df = pd.concat([intra_signal_perbp_df, signal_perbp_temp], axis=0)
        intra_signal_df = pd.concat([intra_signal_df, signal_temp], axis=0)
        inter_signal_df = pd.concat([inter_signal_df, DF_INTER_SIGNAL], axis=0)


intra_signal_perbp_df.reset_index(inplace=True, drop=True)
intra_signal_df.reset_index(inplace=True, drop=True)
inter_signal_df.reset_index(inplace=True, drop=True)

intra_signal_perbp_df_pivot = intra_signal_perbp_df.pivot(index='organism', columns='variable', values='signal')
intra_signal_perbp_df_pivot.reset_index(inplace=True)
# print(intra_signal_perbp_df_pivot)

coeff_1 = scipy.stats.pearsonr(intra_signal_perbp_df_pivot['intron'].values, intra_signal_perbp_df_pivot['exon'].values)

plt.figure(figsize=(20, 10))
# sns.scatterplot(data=intra_signal_perbp_df_pivot, x='intron', y="exon")
sns.regplot(data=intra_signal_perbp_df_pivot, x='intron', y="exon", x_estimator=np.mean, scatter=True, fit_reg=True)
plt.xlabel("Exon Intra-regional signal size per bp")
plt.ylabel("Intron Intra-regional signal size per bp")
plt.title("Correlation between Exon and Intron Intraregional signals")
plt.xticks(rotation=0)
plt.grid()
plt.text(intra_signal_perbp_df_pivot['intron'].values.min()+1, intra_signal_perbp_df_pivot['exon'].values.max()-1, 'Pearson Corr: '+str(round(coeff_1[0],2)), horizontalalignment='left', size='medium', color='black', weight='semibold')
# plt.show()
plt.savefig(path+'Exon_Intron_Intra_perbp.png')


intra_signal_df_pivot = intra_signal_df.pivot(index='organism', columns='variable', values='signal')
intra_signal_df_pivot.reset_index(inplace=True)
# print(intra_signal_df_pivot)
coeff_2 = scipy.stats.pearsonr(intra_signal_df_pivot['intron'].values, intra_signal_df_pivot['exon'].values)

plt.figure(figsize=(20, 10))
# sns.scatterplot(data=intra_signal_df_pivot, x='intron', y="exon")
sns.regplot(data=intra_signal_df_pivot, x='intron', y="exon", x_estimator=np.mean, scatter=True, fit_reg=True)

plt.xlabel("Exon Intra-regional signal size")
plt.ylabel("Intron Intra-regional signal size")
plt.title("Correlation between Exon and Intron Intraregional signals")
plt.xticks(rotation=0)
plt.grid()
plt.text(intra_signal_df_pivot['intron'].values.min()+1, intra_signal_df_pivot['exon'].values.max()-1,  'Pearson Corr: '+str(round(coeff_2[0],2)), horizontalalignment='left', size='medium', color='black', weight='semibold')
# plt.show()
plt.savefig(path+'Exon_Intron_Intra.png')


inter_signal_df_pivot = inter_signal_df.pivot(index='organism', columns='element', values='Avg_GC')
inter_signal_df_pivot.reset_index(inplace=True)

inter_signal_df_pivot = inter_signal_df_pivot[inter_signal_df_pivot['exon']>0]
inter_signal_df_pivot = inter_signal_df_pivot[inter_signal_df_pivot['intron']>0]

coeff_3 = scipy.stats.pearsonr(inter_signal_df_pivot['intron'].values, inter_signal_df_pivot['exon'].values)
plt.figure(figsize=(20, 10))
# sns.scatterplot(data=inter_signal_df_pivot, x='intron', y="exon")
sns.regplot(data=inter_signal_df_pivot, x='intron', y="exon", x_estimator=np.mean, scatter=True, fit_reg=True)

plt.xlabel("Exon Inter-regional signal size")
plt.ylabel("Intron Inter-regional signal size")
plt.title("Correlation between Exon and Intron Inter Regional signals")
plt.xticks(rotation=0)
plt.grid()
plt.text(inter_signal_df_pivot['intron'].values.min()+0.01, inter_signal_df_pivot['exon'].values.max()-0.01,  'Pearson Corr: '+str(round(coeff_3[0],2)), horizontalalignment='left', size='medium', color='black', weight='semibold')
# plt.show()
plt.savefig(path+'Exon_Intron_Inter.png')
