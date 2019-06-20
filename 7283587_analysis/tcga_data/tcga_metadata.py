

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
        nrows = 2,
        ncols = 4,
        figsize = (20, 10),
        subplot_kw = {
            'facecolor':'none'},
        sharex=True, sharey=True) # Create shared axis for cleanliness
    plt.subplots_adjust(
        bottom = 0.1)
    plt.yticks([0,1,2,3,4,5]) # Limit axis labels to ints
    plt.xticks([0,1,2,3,4,5])

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
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymin=0.05, ymax=0.88)
        file_number += 1 # Plot counter
        print(round(rho.astype('float'), 2).astype('str'))

    # Create shared row/column titles
    count_label = ['log$_1$$_0$(counts)','log$_1$$_0$(counts)','log$_1$$_0$(counts)','log$_1$$_0$(counts)']

    cols = file_list[:4]
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=24)
        ax.xaxis.set_label_position('top')

    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(FPKM)', fontsize=16)

    cols = file_list[4:8]
    for ax, col in zip(axes[1], cols):
        ax.set_xlabel(col, fontsize=24)
        ax.xaxis.set_label_position('top')

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/tcga_data/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

"""
IMPORT METADATA
"""
meta = '~/Desktop/xpressyourself_manuscript/7283587_analysis/tcga_data/gdc_sample_sheet.2019-06-19.tsv'
meta = pd.read_csv(meta, sep='\t')

metadata = pd.DataFrame()
metadata['sample_id'] = meta['Sample ID']
metadata['fpkm_path'] = meta['File ID']
metadata['fpkm_name'] = meta['File Name']

"""
IMPORT XPRESSPIPE GENERATED DATA
"""
directory = '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/tcga_data/xpresspipe/'
file_list = []
di = {}

for file in os.listdir(directory):
    if file.endswith('sort.tsv') and os.path.isfile(str(directory) + str(file)) == True:
        file_list.append(str(directory) + str(file))
        di[file.split('.')[0]] = file.split('.')[1][:-12]

xpresspipe = count_table(file_list)
xpresspipe.columns = xpresspipe.columns.to_series().map(di)
xpresspipe = xpresspipe.reindex(sorted(xpresspipe.columns), axis=1)

gtf = '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.96.gtf'
xpresspipe_f = r_fpkm(xpresspipe, gtf)

"""
IMPORT TCGA GENERATED DATA
"""
directory2 = '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/tcga_data/tcga_fpkms/'
file_list2 = []
di2 = {}

for index, row in metadata.iterrows():

    file_list2.append(str(directory2) + str(row[1]) + '/' + row[2][:-3])
    di2[row[2][:-12]] = row[0]

tcga_fpkm = count_table(file_list2)
tcga_fpkm.columns = tcga_fpkm.columns.to_series().map(di2)
tcga_fpkm.index = tcga_fpkm.index.str.split('.').str[0]
tcga_fpkm = tcga_fpkm.loc[:,~tcga_fpkm.columns.duplicated(keep=False)]


"""
GET COMMON SAMPLES AND GENES
"""
xpresspipe_genes = xpresspipe_f.index.tolist()
tcga_fpkm_genes = tcga_fpkm.index.tolist()
len(xpresspipe_genes)
len(tcga_fpkm_genes)

xpresspipe_cols = xpresspipe.columns.tolist()
tcga_fpkm_cols = tcga_fpkm.columns.tolist()
len(xpresspipe_cols)
len(tcga_fpkm_cols)

common_genes = list(set(xpresspipe_genes).intersection(tcga_fpkm_genes))
len(common_genes)

xpresspipe_common = xpresspipe_f.reindex(index = common_genes)
tcga_fpkm_common = tcga_fpkm.reindex(index = common_genes)

common_cols = list(set(xpresspipe_cols).intersection(tcga_fpkm_cols))
len(common_cols)

xpresspipe_common = xpresspipe_common.reindex(columns = common_cols)
tcga_fpkm_common = tcga_fpkm_common.reindex(columns = common_cols)

xpresspipe_common.shape
tcga_fpkm_common.shape

xpresspipe_common = xpresspipe_common.reindex(sorted(xpresspipe_common.columns), axis=1)
tcga_fpkm_common = tcga_fpkm_common.reindex(sorted(tcga_fpkm_common.columns), axis=1)

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
             'TCGA-06-0743-01A',]
             #'TCGA-06-0744-01A',
             #'TCGA-06-0745-01A',
             #'TCGA-08-0386-01A']

make_figure4(
    xpresspipe_common,
    tcga_fpkm_common,
    sample_list,
    'xpresspipe_vs_tcga.png')
