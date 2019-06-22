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
