"""
IMPORT DEPENDENCIES
"""
import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
matplotlib.rcParams['font.sans-serif'] = 'Arial'
import seaborn as sns
%matplotlib inline

from xpressplot import r_fpkm, count_table

"""
FUNCTIONS
"""
def make_figure4(
    data,
    tcga_data,
    file_list,
    title):

    fig, axes = plt.subplots(
        nrows = 3,
        ncols = 4,
        figsize = (20, 15),
        subplot_kw = {
            'facecolor':'none'},
        sharex=True, sharey=True) # Create shared axis for cleanliness
    plt.subplots_adjust(
        bottom = 0.1)
    plt.yticks([0,1,2,3,4,5,6]) # Limit axis labels to ints
    plt.xticks([0,1,2,3,4,5,6])

    file_number = 0

    for x in file_list:

        # Get data as array-like for samples being compared
        data_c1 = data.copy()
        data_c2 = tcga_data.copy()

        sample_a = data_c1[x].values.tolist()
        sample_a = [x + 1 for x in sample_a]
        sample_a = np.array(sample_a).astype(np.float)
        sample_a = np.ndarray.tolist(sample_a)

        sample_b = data_c2[x].values.tolist()
        sample_b = [x + 1 for x in sample_b]
        sample_b = np.array(sample_b).astype(np.float)
        sample_b = np.ndarray.tolist(sample_b)

        # Run Spearman R linreg for non-normal data
        rho, p_value = stats.spearmanr(sample_a, sample_b)

        # Determine subplot location
        if file_number in [0,1,2,3]:
            ax_x = file_number % 4
            ax_y = 0
        elif file_number in [4,5,6,7]:
            ax_x = file_number - 4
            ax_y = 1
        elif file_number in [8,9,10,11]:
            ax_x = file_number - 8
            ax_y = 2
        elif file_number in [12,13,14,15]:
            ax_x = file_number - 12
            ax_y = 3
        else:
            print('oops')

        # Format p value
        if p_value.astype('float') < 0.001:
            p_val = '< 0.001'
        else:
            p_val = round(p_value.astype('float'), 4).astype('str')

        # Plot data
        axes[ax_y, ax_x].scatter(np.log10(sample_a), np.log10(sample_b), s=1,c='black')
        axes[ax_y, ax_x].set_title('R = ' + round(rho.astype('float'), 2).astype('str') + '\nP ' + p_val, y=0.1, x=0.9, fontsize=16) # Format titles
        axes[ax_y, ax_x].axhline(0, ls='-', color='black', xmin=0.05, xmax=1) # Create axis lines
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymin=0.05, ymax=1)
        file_number += 1 # Plot counter
        print(round(rho.astype('float'), 2).astype('str'))

    # Create shared row/column titles
    cols = file_list[:4]
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=24)
        ax.xaxis.set_label_position('top')

    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(counts)', fontsize=16)

    cols = file_list[4:8]
    for ax, col in zip(axes[1], cols):
        ax.set_xlabel(col, fontsize=24)
        ax.xaxis.set_label_position('top')

    cols = file_list[8:11]
    for ax, col in zip(axes[2], cols):
        ax.set_xlabel(col, fontsize=24)
        ax.xaxis.set_label_position('top')

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/tcga_data/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')


def make_figure4S(
    data,
    tcga_data,
    x,
    title,
    file_number,
    axes,
    name):

    # Get data as array-like for samples being compared
    data_c1 = data.copy()
    data_c2 = tcga_data.copy()

    sample_a = data_c1[x].values.tolist()
    sample_a = [x + 1 for x in sample_a]
    sample_a = np.array(sample_a).astype(np.float)
    sample_a = np.ndarray.tolist(sample_a)

    sample_b = data_c2[x].values.tolist()
    sample_b = [x + 1 for x in sample_b]
    sample_b = np.array(sample_b).astype(np.float)
    sample_b = np.ndarray.tolist(sample_b)

    # Run Spearman R linreg for non-normal data
    rho, p_value = stats.spearmanr(sample_a, sample_b)

    # Determine subplot location
    if file_number in [0,1,2,3]:
        ax_x = file_number % 4
        ax_y = 0
    elif file_number in [4,5,6,7]:
        ax_x = file_number - 4
        ax_y = 1
    elif file_number in [8,9,10,11]:
        ax_x = file_number - 8
        ax_y = 2
    elif file_number in [12,13,14,15]:
        ax_x = file_number - 12
        ax_y = 3
    else:
        print('oops')

    # Format p value
    if p_value.astype('float') < 0.001:
        p_val = '< 0.001'
    else:
        p_val = round(p_value.astype('float'), 4).astype('str')

    # Plot data
    axes[ax_y, ax_x].scatter(np.log10(sample_a), np.log10(sample_b), s=1,c='black')
    axes[ax_y, ax_x].set_title('R = ' + round(rho.astype('float'), 2).astype('str') + '\nP ' + p_val, y=0.1, x=0.9, fontsize=16) # Format titles
    axes[ax_y, ax_x].axhline(0, ls='-', color='black', xmin=0.05, xmax=1) # Create axis lines
    axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymin=0.05, ymax=1)

    print(round(rho.astype('float'), 2).astype('str'))


    axes[ax_y, ax_x].set_xlabel(name, fontsize=12)
    axes[ax_y, ax_x].xaxis.set_label_position('top')

    axes[ax_y, ax_x].set_ylabel('log$_1$$_0$(counts)', fontsize=10)

    return axes


"""
IMPORT METADATA
"""
meta = '~/Desktop/xpressyourself_manuscript/tcga_data/gdc_sample_sheet.2019-06-21.tsv'
meta = pd.read_csv(meta, sep='\t')

metadata = pd.DataFrame()
metadata['sample_id'] = meta['Sample ID']
metadata['count_path'] = meta['File ID']
metadata['count_name'] = meta['File Name']

"""
IMPORT XPRESSPIPE GENERATED DATA
"""
x_path = '/Users/jordan/Desktop/xpressyourself_manuscript/tcga_data/xpresspipe/'

x_samples = [
    'G17189.TCGA-06-0132-01A-02R-1849-01.2.bam__Aligned.sort.tsv',
    'G17190.TCGA-06-0174-01A-01R-1849-01.2.bam__Aligned.sort.tsv',
    'G17193.TCGA-06-0743-01A-01R-1849-01.2.bam__Aligned.sort.tsv',
    'G17195.TCGA-06-0138-01A-02R-1849-01.2.bam__Aligned.sort.tsv',
    'G17197.TCGA-06-0211-01B-01R-1849-01.2.bam__Aligned.sort.tsv',
    'G17199.TCGA-06-0744-01A-01R-1849-01.2.bam__Aligned.sort.tsv',
    'G17202.TCGA-06-0184-01A-01R-1849-01.2.bam__Aligned.sort.tsv',
    'G17203.TCGA-06-0211-02A-02R-2005-01.2.bam__Aligned.sort.tsv',
    'G17204.TCGA-08-0386-01A-01R-1849-01.2.bam__Aligned.sort.tsv',
    'G17205.TCGA-06-0745-01A-01R-1849-01.2.bam__Aligned.sort.tsv',
    'G17206.TCGA-06-0125-02A-11R-2005-01.2.bam__Aligned.sort.tsv',
    'G17207.TCGA-06-0156-01A-03R-1849-01.2.bam__Aligned.sort.tsv',
    ]

x_path_samples = [str(x_path) + x for x in x_samples]
xpresspipe = count_table(x_path_samples)

di = {}
for x in x_samples:
    name = x.split('.')[0]
    id = x.split('.')[1][:-12]
    di[name] = id

xpresspipe.columns = xpresspipe.columns.to_series().map(di)
xpresspipe.shape

"""
IMPORT TCGA GENERATED DATA
"""
directory2 = '/Users/jordan/Desktop/xpressyourself_manuscript/tcga_data/tcga_counts/'
file_list2 = []
di2 = {}

for index, row in metadata.iterrows():

    file_list2.append(str(directory2) + row[2][:-3])
    di2[row[2].split('.')[0]] = row[0]

tcga = count_table(file_list2)
tcga.columns = tcga.columns.to_series().map(di2)
tcga.index = tcga.index.str.split('.').str[0]
tcga = tcga.loc[:,~tcga.columns.duplicated(keep=False)]
tcga.shape

"""
GET COMMON SAMPLES AND GENES
"""
xpresspipe_genes = xpresspipe.index.tolist()
tcga_genes = tcga.index.tolist()
len(xpresspipe_genes)
len(tcga_genes)

xpresspipe_cols = xpresspipe.columns.tolist()
tcga_cols = tcga.columns.tolist()
len(xpresspipe_cols)
len(tcga_cols)

common_genes = list(set(xpresspipe_genes).intersection(tcga_genes))
len(common_genes)

xpresspipe_common = xpresspipe.reindex(index = common_genes)
tcga_common = tcga.reindex(index = common_genes)

common_cols = list(set(xpresspipe_cols).intersection(tcga_cols))
len(common_cols)

# Missing for some reason
set(xpresspipe_cols) - set(common_cols)

xpresspipe_common = xpresspipe_common.reindex(columns = common_cols)
tcga_common = tcga_common.reindex(columns = common_cols)

xpresspipe_common.shape
tcga_common.shape

xpresspipe_common = xpresspipe_common.reindex(sorted(xpresspipe_common.columns), axis=1)
tcga_common = tcga_common.reindex(sorted(tcga_common.columns), axis=1)

"""
PLOT FIGURE 4
"""
sample_list = ['TCGA-06-0125-02A',
             'TCGA-06-0132-01A',
             'TCGA-06-0138-01A',
             'TCGA-06-0174-01A',
             'TCGA-06-0184-01A',
             'TCGA-06-0211-01B',
             'TCGA-06-0211-02A',
             'TCGA-06-0743-01A',
             'TCGA-06-0744-01A',
             'TCGA-06-0745-01A',
             'TCGA-08-0386-01A']

make_figure4(
    xpresspipe_common,
    tcga_common,
    sample_list,
    'xpresspipe_vs_tcga_counts_more.png')




"""
Plot a sample using v79 Ensemble build (gencode v22)
"""
# Sample quantified to v22 gencode
v22_path = '/Users/jordan/Desktop/xpressyourself_manuscript/tcga_data/count_gencode22/'
v22_counts = pd.read_csv(
    str(v22_path) + 'G17189.TCGA-06-0132-01A_justsorted_gencodev22.tsv',
    sep='\t',
    header=None,
    index_col=0)

del v22_counts.index.name
v22_counts.columns = ['TCGA-06-0132-01A']
v22_counts.index = v22_counts.index.str.split('.').str[0]

v22_counts = v22_counts.iloc[:-5]
v22_counts.shape

v22_counts_genes = v22_counts.index.tolist()
tcga_v22_genes = tcga[['TCGA-06-0132-01A']].index.tolist()
len(v22_counts_genes)
len(tcga_v22_genes)

common_v22_genes = list(set(v22_counts_genes).intersection(tcga_v22_genes))
len(common_v22_genes)

v22_counts_common = v22_counts.reindex(index = common_v22_genes)
tcga_v22_common = tcga[['TCGA-06-0132-01A']].reindex(index = common_v22_genes)

v22_counts_common.shape
tcga_v22_common.shape

v22_counts_common = v22_counts_common.reindex(sorted(v22_counts_common.columns), axis=1)
tcga_v22_common = tcga_v22_common.reindex(sorted(tcga_v22_common.columns), axis=1)


make_figure4(
    v22_counts_common,
    tcga_v22_common,
    ['TCGA-06-0132-01A'],
    'xpresspipe_vs_tcga_counts_v22_comparison.png')



comp_dir = '/Users/jordan/Desktop/xpressyourself_manuscript/tcga_data/tcga_comp/'
comp_samples = [
    'v79_looseTrim_longestGTF_multimappers_count_table',
    'v79_looseTrim_normalGTF_multimappers_count_table',
    'v79_strictTrim_longestGTF_multimappers_count_table',
    'v79_strictTrim_normalGTF_multimappers_count_table',
    'v96_looseTrim_normalGTF_multimappers_count_table',
    'v96_strictTrim_longestGTF_multimappers_count_table',
    'v96_strictTrim_normalGTF_multimappers_count_table',
    'v79_looseTrim_normalGTF_uniqueOnly_count_table',
    'htseq_oldtest_count_table'
]

comp_files = []
for x in comp_samples:
    comp_files.append([str(x)[:-12], str(comp_dir) + str(x) + '.tsv'])

fig, axes = plt.subplots(
    nrows = 4,
    ncols = 4,
    figsize = (20, 20),
    subplot_kw = {
        'facecolor':'none'},
    sharex=True, sharey=True) # Create shared axis for cleanliness
plt.subplots_adjust(
    bottom = 0.1)
plt.yticks([0,1,2,3,4,5,6]) # Limit axis labels to ints
plt.xticks([0,1,2,3,4,5,6])
file_number = 0

for x in comp_files:

    print(x[0])

    counts = pd.read_csv(
        str(x[1]),
        sep='\t',
        header=None,
        index_col=0)

    del counts.index.name
    counts.columns = [str(x[0])]
    counts.index = counts.index.str.split('.').str[0]
    counts = counts.dropna()
    counts = counts.iloc[:-5]
    counts.columns = ['TCGA-06-0132-01A']

    counts_genes = counts.index.tolist()
    tcga_v22_genes = tcga[['TCGA-06-0132-01A']].index.tolist()
    len(counts_genes)
    len(tcga_v22_genes)

    common_v22_genes = list(set(counts_genes).intersection(tcga_v22_genes))
    len(common_v22_genes)

    counts_common = counts.reindex(index = common_v22_genes)
    tcga_v22_common = tcga[['TCGA-06-0132-01A']].reindex(index = common_v22_genes)

    counts_common.shape
    tcga_v22_common.shape

    counts_common = counts_common.reindex(sorted(counts_common.columns), axis=1)
    tcga_v22_common = tcga_v22_common.reindex(sorted(tcga_v22_common.columns), axis=1)

    axes = make_figure4S(
        counts_common,
        tcga_v22_common,
        'TCGA-06-0132-01A',
        str(x[0]) + '.png',
        file_number,
        axes,
        x[0])

    file_number += 1
fig.savefig(
    '/Users/jordan/Desktop/xpressyourself_manuscript/tcga_data/tcga_comparisons_table.png',
    dpi = 600,
    bbox_inches = 'tight')




import plotly_express as px
best_file = '/Users/jordan/Desktop/xpressyourself_manuscript/tcga_data/tcga_comp/v79_looseTrim_normalGTF_multimappers_count_table.tsv'


worst_file = '/Users/jordan/Desktop/xpressyourself_manuscript/tcga_data/tcga_comp/v96_strictTrim_longestGTF_multimappers_count_table.tsv'



best_name = 'v79_looseTrim_normalGTF_multimappers'
best_counts = pd.read_csv(
    str(best_file),
    sep='\t',
    header=None,
    index_col=0)

del best_counts.index.name
best_counts.columns = [best_name]
best_counts.index = best_counts.index.str.split('.').str[0]
best_counts = best_counts.dropna()
best_counts = best_counts.iloc[:-5]

counts_genes = best_counts.index.tolist()
tcga_v22_genes = tcga[['TCGA-06-0132-01A']].index.tolist()
len(counts_genes)
len(tcga_v22_genes)

common_v22_genes = list(set(counts_genes).intersection(tcga_v22_genes))
len(common_v22_genes)

counts_common = best_counts.reindex(index = common_v22_genes)
tcga_v22_common = tcga[['TCGA-06-0132-01A']].reindex(index = common_v22_genes)

counts_common.shape
tcga_v22_common.shape

counts_common = counts_common.reindex(sorted(counts_common.columns), axis=1)
tcga_v22_common = tcga_v22_common.reindex(sorted(tcga_v22_common.columns), axis=1)

merged_best = pd.concat([counts_common, tcga_v22_common], axis=1, sort=False)
merged_best.head()
merged_best.columns.tolist()[0]

from xpressplot import convert_names
merged_best = convert_names(
    merged_best,
    '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.96.gtf')
merged_best['genes'] = merged_best.index


# Import reference GTF
gtf = pd.read_csv(str('/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.96.gtf'),sep='\t',comment='#', low_memory=False, header=None)

# Parse out old and new names from GTF
gtf_genes = gtf.loc[gtf[2] == 'gene']
gtf_genes['id'] = gtf[8].str.split(';').str[0]
gtf_genes['original'] = gtf[8].str.split(';').str[2]
gtf_genes['new'] = gtf[8].str.split(';').str[4]
gtf_genes['id'] = gtf_genes['id'].map(lambda x: x.lstrip(str('gene_id \"')).rstrip('\"').rstrip(' '))
gtf_genes['original'] = gtf_genes['original'].map(lambda x: x.lstrip(str('gene_name \"')).rstrip('\"').rstrip(' '))
gtf_genes['new'] = gtf_genes['new'].map(lambda x: x.lstrip(str('gene_biotype \"')).rstrip('\"').rstrip(' '))
gtf_genes = gtf_genes[['id','original','new']].copy()

gtf_genes.head()

gene_dict = pd.Series(gtf_genes['new'].values,index=gtf_genes['original']).to_dict()
merged_best['color'] = merged_best.index.to_series().map(gene_dict)

merged_best.head()

merged_best_c = merged_best.copy()
merged_best_c = merged_best_c.reset_index()
del merged_best_c['index']
merged_best_c = merged_best_c.fillna('other')

merged_best_c.loc[merged_best_c['color'].str.contains('coding'), 'color'] = 'protein_coding'
merged_best_c.loc[merged_best_c['color'].str.contains('pseudogene'), 'color'] = 'pseudogene'
merged_best_c.loc[(~merged_best_c['color'].str.contains('coding')) & (~merged_best_c['color'].str.contains('pseudogene')), 'color'] = 'other'
merged_best_c.index = merged_best_c['genes']
del merged_best_c.index.name

merged_best_c.head()

set(merged_best_c.color.tolist())


labels = {
    'TCGA-06-0132-01A':'TCGA processed',
    best_name: 'xpresspipe: GTFv96, PHRED>28, longestGTF, allow multimappers'
}



sc = px.scatter(
    merged_best_c,
    x=merged_best_c.columns.tolist()[0],
    y=merged_best_c.columns.tolist()[1],
    color='color',
    hover_name='genes',
    log_x=True,
    log_y=True,
    opacity=0.4,
    width=1400,
    height=1000,
    labels=labels,
    title='TCGA-06-0132-01A')

import plotly
import plotly.offline as py
py.offline.plot(sc, filename='/Users/jordan/Desktop/tcga_comparison_v96_stricttrim_longestGTF_multimappers.html')
