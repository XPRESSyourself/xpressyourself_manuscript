"""
XPRESSpipe
An alignment and analysis pipeline for RNAseq data
alias: xpresspipe

Copyright (C) 2019  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import print_function

"""Import dependencies
"""
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
matplotlib.rcParams['font.sans-serif'] = 'Arial'
from matplotlib_venn import venn2, venn3
import seaborn as sns
from xpressplot import rpm

"""Set globals
"""
ingolia_tm = [
    'FAM50B','CHST2','LNP1','AX747756','ACRC',
    'ATF5','ANKRD23','ATF4','C19orf48','TEX14',
    'SLC35A4','DDIT3','RAB35','AK293147','PTP4A1',
    'ANKRD36C','NPIP','ANKMY1','UCP2','NPIPL1',
    'PNRC2','C7orf31','PPP1R15A','IFRD1','BCL2L11',
    'OBSCN','SLC24A1','NPIPL1','NPIPL1','ZNF22',
    'NAIP','NPIP','DNAH1','KIAA1377','FMO4',
    'CCDC150','NPIP','AK310228','ZSCAN5A','FAM13B',
    'NPIPL1','ZNF79','HOXA6','NPIPL2','SLMO1',
    'ARMCX4','AK310228','LOC100288332','SAT1','FLJ00285',
    'AK097143','AX748249','AGAP4','FLJ00322','GOLGA8A',
    'GOLGA8B','ZNF432','FSBP','CPS1','MATK',
    'EPB41L4B','HCN2','GLRB','KIAA1841','MYO5B',
    'EHBP1L1','MAPK8IP1','SH3RF3','DET1','PTPN21',
    'VANGL2','MMP15','MAP3K10','PAQR5','GAL',
    'RPP25','FLJ44955','DQ576756']

ingolia_tmisrib = [
    'AK310228','AK310228','NPIP','AK302511','NPIP',
    'LOC100288332','STARD9','PKD1','ANKRD36']

ingolia_isrib = [
    'AK310228','AK310228','AK302511','NPIP']

# Initialize lists of Ingolia highlights
isr = [
    'ATF4',
    'ATF5',
    'PPP1R15A',
    'DDIT3']

uorf_targets = [
    'SLC35A4',
    'PTP4A1',
    'UCP2',
    'C7orf31',
    'BCL2L11',
    'PNRC2',
    'SAT1']

"""Prep data vectors
"""
def prep_data(
        data,
        y_name,
        x_name):

    ribo_values = data[[y_name]].values.tolist()
    ribo_values = np.array(ribo_values).astype(np.float)
    ribo_values = np.ndarray.tolist(ribo_values)

    rna_values = data[[x_name]].values.tolist()
    rna_values = np.array(rna_values).astype(np.float)
    rna_values = np.ndarray.tolist(rna_values)

    return ribo_values, rna_values

"""Plot optional data
"""
def plot_optional(
        ax,
        ribo_values,
        rna_values,
        data,
        data_subset,
        color,
        alpha=1,
        regulation_list=None,
        significance_list=None,
        label=False,
        optional_genes_1=[],
        optional_genes_2=[]):

    # Plot first pass data
    ax.scatter(
        rna_values,
        ribo_values,
        s=80,
        c=color,
        alpha=alpha)

    # Plot significance coloring
    for index, row in data.iterrows():
        if index in regulation_list:
            if index in significance_list:
                ax.scatter(
                    row[1],
                    row[0],
                    s=20,
                    c='black',
                    alpha=1)

            else:
                ax.scatter(
                    row[1],
                    row[0],
                    s=20,
                    c='gray',
                    alpha=1)

        else:
            pass

    # Add labels
    if label == True:

        for index, row in data_subset.iterrows():

            if index in optional_genes_1 + optional_genes_2:

                for x in optional_genes_1:

                    if index == x:
                        ax.text(
                            row[1] + 0.1,
                            row[0] - 0.07,
                            str(index),
                            horizontalalignment='left',
                            size='medium',
                            color=color,
                            weight='semibold')

                for y in optional_genes_2:

                    if index == y:
                        ax.text(
                            row[1] - 0.05,
                            row[0] - 0.2,
                            str(index),
                            horizontalalignment='right',
                            size='medium',
                            color=color,
                            weight='semibold')

            else:
                ax.text(
                    row[1] - 0.1,
                    row[0] - 0.07,
                    str(index),
                    horizontalalignment='right',
                    size='medium',
                    color=color,
                    weight='semibold')

    return ax

"""Plot ribo-seq data on X/Y figure
"""
def rp_plot(
        data,
        y_fc,
        x_fc,
        y_p,
        x_p,
        y_name,
        x_name,
        te_data,
        te_p,
        title,
        list_up,
        list_down,
        list_down_custom,
        genes=[],
        label_up=True,
        label_down=True,
        label_down_custom=True):

    # Get relevant data
    data_fig = data[[y_fc, x_fc, y_p, x_p]].copy()

    # Remove infinite and missing values
    data_fig[[y_fc, x_fc]] = data_fig[[y_fc, x_fc]].replace([np.inf, -np.inf], np.nan)
    data_fig = data_fig.dropna(subset=[y_fc, x_fc])


    # Initialize plotting space
    fig, ax = plt.subplots(1,1, figsize=(7,7))
    plt.grid(False)
    ax.axhline(0, ls='-', color='black') # Create central axis line
    ax.axvline(0, ls='-', color='black')
    rect = patches.Rectangle((-1,-1),2,2,linewidth=1.5,edgecolor='gray',facecolor='none') # Create significance zone
    ax.add_patch(rect)
    ax.set_ylabel('Ribo-Seq log$_2$(' + str(y_name) + '/' + str(x_name) + ')', fontsize=16) # Set axis labels
    ax.set_xlabel('mRNA-Seq log$_2$(' + str(y_name) + '/' + str(x_name) + ')', fontsize=16)
    ax.set_xlim(-3.5,3.5) # Set axis limits
    ax.set_ylim(-3.5,3.5)
    x = [-3,-2,-1,0,1,2,3] # Set axis spacing
    ticks = ["-3","-2","-1","0","1","2","3"]
    ax.set_facecolor("#FFFFFF") # Set background color

    # Set X and Y axis scales
    ax.set_xticks(x)
    ax.set_xticklabels(ticks, fontsize=16)
    ax.set_xticklabels([""]*len(x), minor=True)

    ax.set_yticks(x)
    ax.set_yticklabels(ticks, fontsize=16)
    ax.set_yticklabels([""]*len(x), minor=True)

    # Prep data for plotting
    ribo_all, rna_all = prep_data(
        data=data_fig,
        y_name=y_fc,
        x_name=x_fc)

    # Prep significant hits for plotting
    sig_list = te_data.loc[te_data[te_p] < 0.1].index.tolist()
    data_sig = data_fig.reindex(sig_list)
    ribo_sig, rna_sig = prep_data(
        data=data_sig,
        y_name=y_fc,
        x_name=x_fc)

    # Prep option gene sets
    if list_up != None:
        data_fig_up = data_fig.reindex(list_up)
        ribo_up, rna_up = prep_data(
            data=data_fig_up,
            y_name=y_fc,
            x_name=x_fc)

    if list_down != None:
        data_fig_down = data_fig.reindex(list_down)
        ribo_down, rna_down = prep_data(
            data=data_fig_down,
            y_name=y_fc,
            x_name=x_fc)

    if list_down_custom != None:
        data_fig_down_custom = data_fig.reindex(list_down_custom)
        ribo_down_custom, rna_down_custom = prep_data(
            data=data_fig_down_custom,
            y_name=y_fc,
            x_name=x_fc)

    #Plot data
    ax.scatter(rna_all, ribo_all, s=2.5,c='gray',alpha=0.5)
    ax.scatter(rna_sig, ribo_sig, s=5,c='black',alpha=1)

    if list_down != None:

        ax = plot_optional(
            ax=ax,
            ribo_values=ribo_down,
            rna_values=rna_down,
            data=data_fig,
            data_subset=data_fig_down,
            color='#1b9e77',
            regulation_list=list_down,
            significance_list=sig_list,
            label=label_down,
            optional_genes_1=[
                'PNRC2',
                'SAT1']
            )

    if list_down_custom != None:

        ax = plot_optional(
            ax=ax,
            ribo_values=ribo_down_custom,
            rna_values=rna_down_custom,
            data=data_fig,
            data_subset=data_fig_down_custom,
            color='#d95f02',
            regulation_list=list_down_custom,
            significance_list=sig_list,
            label=label_down_custom,
            optional_genes_1=[
                'RPS15A',
                'HIST2H3D',
                'RPL27'],
            optional_genes_2=[
                'HSPA8']
            )

    if list_up != None:

        ax = plot_optional(
            ax=ax,
            ribo_values=ribo_up,
            rna_values=rna_up,
            data=data_fig,
            data_subset=data_fig_up,
            color='#7570b3',
            regulation_list=list_up,
            significance_list=sig_list,
            label=label_up,
            optional_genes_1=[
                'ATF4',
                'DDIT3',
                'PPP1R15A']
            )

    plt.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/' + str(title),
        dpi=3600,
        bbox_inches='tight')

    plt.close()

"""Import DE data
"""
def parse_de_data(
        directory,
        check_file):

    check_data = pd.read_csv(
        str(directory) + str(check_file),
        sep='\t',
        index_col=0)
    check_list = check_data.index.tolist()

    # Import DESeq2 TE data
    tm_data = pd.read_csv(str(directory) + 'tm_counts_diffx.tsv', sep='\t')
    tmisrib_data = pd.read_csv(str(directory) + 'tmisrib_counts_diffx.tsv', sep='\t')
    isrib_data = pd.read_csv(str(directory) + 'isrib_counts_diffx.tsv', sep='\t')

    ### Clean up this data
    tm_data = tm_data.drop(['baseMean','lfcSE','stat'],axis=1)
    tm_data.columns = ['tm_log2FC','tm_p','tm_padj']

    tmisrib_data = tmisrib_data.drop(['baseMean','lfcSE','stat'],axis=1)
    tmisrib_data.columns = ['tmisrib_log2FC','tmisrib_p','tmisrib_padj']

    isrib_data = isrib_data.drop(['baseMean','lfcSE','stat'],axis=1)
    isrib_data.columns = ['isrib_log2FC','isrib_p','isrib_padj']

    merged_data = pd.concat([tm_data, tmisrib_data, isrib_data], axis=1, sort=False)

    # Import DESeq2 RNA and RPF data
    tm_ribo_data = pd.read_csv(str(directory) + 'tm_ribo_counts_diffx.tsv', sep='\t')
    tm_rna_data = pd.read_csv(str(directory) + 'tm_rna_counts_diffx.tsv', sep='\t')

    tmisrib_ribo_data = pd.read_csv(str(directory) + 'tmisrib_ribo_counts_diffx.tsv', sep='\t')
    tmisrib_rna_data = pd.read_csv(str(directory) + 'tmisrib_rna_counts_diffx.tsv', sep='\t')

    isrib_ribo_data = pd.read_csv(str(directory) + 'isrib_ribo_counts_diffx.tsv', sep='\t')
    isrib_rna_data = pd.read_csv(str(directory) + 'isrib_rna_counts_diffx.tsv', sep='\t')

    ### Clean up this data
    tm_ribo_data = tm_ribo_data.drop(['baseMean','lfcSE','stat','pvalue'],axis=1)
    tm_ribo_data.columns = ['tm_ribo_log2FC','tm_ribo_padj']
    tm_rna_data = tm_rna_data.drop(['baseMean','lfcSE','stat','pvalue'],axis=1)
    tm_rna_data.columns = ['tm_rna_log2FC','tm_rna_padj']

    tmisrib_ribo_data = tmisrib_ribo_data.drop(['baseMean','lfcSE','stat','pvalue'],axis=1)
    tmisrib_ribo_data.columns = ['tmisrib_ribo_log2FC','tmisrib_ribo_padj']
    tmisrib_rna_data = tmisrib_rna_data.drop(['baseMean','lfcSE','stat','pvalue'],axis=1)
    tmisrib_rna_data.columns = ['tmisrib_rna_log2FC','tmisrib_rna_padj']

    isrib_ribo_data = isrib_ribo_data.drop(['baseMean','lfcSE','stat','pvalue'],axis=1)
    isrib_ribo_data.columns = ['isrib_ribo_log2FC','isrib_ribo_padj']
    isrib_rna_data = isrib_rna_data.drop(['baseMean','lfcSE','stat','pvalue'],axis=1)
    isrib_rna_data.columns = ['isrib_rna_log2FC','isrib_rna_padj']

    merged_data_split = pd.concat([tm_ribo_data, tm_rna_data, tmisrib_ribo_data, tmisrib_rna_data, isrib_ribo_data, isrib_rna_data], axis=1, sort=False)

    return merged_data, merged_data_split, check_data

"""TE sanity check
"""
def te_sanity(
        data,
        check_data,
        gene_name):

    check_data_rpm = rpm(check_data)

    print(gene_name)

    print('DESeq2 estimate: ', data.loc[gene_name]['tm_log2FC'])

    # Get gene info for all samples
    gene = check_data_rpm.loc[gene_name][['ribo_tm_a', 'ribo_tm_b', 'ribo_untr_a', 'ribo_untr_b','tm_a_hek','tm_b_hek','untr_a_hek','untr_b_hek']]

    # Calculate Tm TE by summing replicates and calculating RPF / RNA
    gene['tm_te'] = (gene['ribo_tm_a'] + gene['ribo_tm_b']) / (gene['tm_a_hek'] + gene['tm_b_hek'])
    gene['untr_te'] = (gene['ribo_untr_a'] + gene['ribo_untr_b']) / (gene['untr_a_hek'] + gene['untr_b_hek'])

    # Calculate FC of TEs for Treatment / not-treated and take the log2 value of that outpu
    gene['fold_change'] = np.log2(gene['tm_te'] / gene['untr_te'])
    print('Sanity estimate: ', gene['fold_change'])

"""Plot TE for genes of interest
"""
def te_plot(
        data,
        list_up,
        list_down,
        directory,
        plot_all=False,
        dpi=3600):

    # Init plotting space
    fig, ax = plt.subplots()
    plt.yticks([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
    ax.tick_params(axis='x', which='minor',length=0)

    ax.set_facecolor('white')
    ax.grid(color='grey', axis='y')

    # Plot all genes' TE
    if plot_all == True:
        data_de_plot.T.plot.line(
            legend=False,
            color='lightgrey',
            linewidth=0.5,
            ax=ax)

    # Set grid lines
    ax.axvline(0.01, ls='-', color='black')
    ax.axvline(1, ls='-', color='black')
    ax.axvline(2, ls='-', color='black')
    ax.axvline(2.99, ls='-', color='black')
    ax.axhline(0, ls='-', color='black')
    ax.axhline(1, ls='--', color='black')
    ax.axhline(-1, ls='--', color='black')

    # Plot genes of interest
    data_de_plot.loc[list_up].T.plot.line(
        legend=False,
        color='#7570b3',
        linewidth=1,
        ax=ax)
    data_de_plot.loc[list_down].T.plot.line(
        legend=False,
        color='#d95f02',
        linewidth=1,
        ax=ax)
    ax.set_ylabel(u'log$_2$(ΔTE)')

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color='gray', lw=2, label='All'),
                       Line2D([0], [0], color='#7570b3', lw=2, label='ISR'),
                       Line2D([0], [0], color='#d95f02', lw=2, label='Other')]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.savefig(
        str(directory) + 'te_analysis_down_strict.png',
        dpi=dpi,
        bbox_inches='tight')

    plt.close()

"""Plot Venn diagrams
"""
def venn_plot(
        directory,
        figure_title,
        output_name,
        gene_list_A,
        gene_list_B,
        label_A,
        label_B,
        color_A='#d95f02',
        color_B='#7570b3',
        dpi=3600):

    plt.title(
        figure_title,
        fontsize = 28)

    c = venn2(
        [set(gene_list_A), set(gene_list_B)],
        set_labels = (label_A, label_B))
    c.get_patch_by_id('10').set_color(color_A)
    c.get_patch_by_id('01').set_color(color_B)
    c.get_patch_by_id('10').set_edgecolor('black')
    c.get_patch_by_id('01').set_edgecolor('black')

    for text in c.set_labels:

        text.set_fontsize(22)

    try:
        for text in c.subset_labels:

            text.set_fontsize(18)

    except:
        pass

    plt.savefig(
        directory + output_name,
        dpi=dpi,
        bbox_inches='tight')

    plt.close()


"""Import XPRESSpipe counts data for any checking
"""
dir = '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/xpresspipe_data_deseq2/'

merged_data_XP, merged_data_split_XP, check_data_XP = parse_de_data(
        directory=dir,
        check_file='ISRIBxpresspipe_thresholded_counts.tsv')

"""
Some of the genes in plotting cluster on a single y-axis coordinate
This is due to insufficient RPF counts after RNA thresholding of the data
"""
merged_data_split_XP[(merged_data_split_XP['tm_ribo_log2FC'] > 2.7) & (merged_data_split_XP['tm_ribo_log2FC'] < 2.75)][['tm_ribo_log2FC', 'tm_ribo_padj']]
""""""

# Check fold changes of targets for comparison
print('ATF4')
print(2**merged_data_split_XP.loc['ATF4']['tm_ribo_log2FC'])
print('ATF5')
print(2**merged_data_split_XP.loc['ATF5']['tm_ribo_log2FC'])
print('PPP1R15A')
print(2**merged_data_split_XP.loc['PPP1R15A']['tm_ribo_log2FC'])
print('DDIT3')
print(2**merged_data_split_XP.loc['DDIT3']['tm_ribo_log2FC'])

te_sanity(
    data=merged_data_XP[['tm_log2FC', 'tm_padj']],
    check_data=check_data_XP,
    gene_name='ATF4')
te_sanity(
    data=merged_data_XP[['tm_log2FC', 'tm_padj']],
    check_data=check_data_XP,
    gene_name='DDIT3')
te_sanity(
    data=merged_data_XP[['tm_log2FC', 'tm_padj']],
    check_data=check_data_XP,
    gene_name='UCP2')
te_sanity(
    data=merged_data_XP[['tm_log2FC', 'tm_padj']],
    check_data=check_data_XP,
    gene_name='POMGNT1')
te_sanity(
    data=merged_data_XP[['tm_log2FC', 'tm_padj']],
    check_data=check_data_XP,
    gene_name='NDUFA11')

# TE analysis
data_de_plot = merged_data_XP.copy()

data_de_plot['Untreated'] = 0
data_de_plot = data_de_plot[['Untreated','tm_log2FC', 'tmisrib_log2FC','isrib_log2FC']]
data_de_plot.columns = ['Untreated','Tm', 'Tm + ISRIB','ISRIB']

down_strict_1_XP = merged_data_XP.loc[
    (merged_data_XP['tm_log2FC'] - merged_data_XP['tmisrib_log2FC'] <= -1) \
    & (merged_data_XP['tm_padj'] <= 0.001) \
    & (merged_data_XP['tmisrib_padj'] <= 0.001) \
    & ((merged_data_split_XP['tm_ribo_padj'] <= 0.1))]

down_strict_2_XP = merged_data_XP.loc[
    (merged_data_XP['tm_log2FC'] - merged_data_XP['tmisrib_log2FC'] <= -1) \
    & (merged_data_split_XP['tm_ribo_log2FC'] <= -0.9) \
    & (merged_data_split_XP['tm_ribo_padj'] <= 0.1)]

down_strict_XP = down_strict_1_XP.index.tolist() + down_strict_2_XP.index.tolist()

# Sanity check count data
print(check_data_XP.loc[down_strict_XP].min().min())

merged_data_XP.loc[down_strict_XP][['tm_log2FC', 'tm_padj', 'tmisrib_log2FC', 'tmisrib_padj']]

2**.51

# Considering alpha-debt of dataset decay (https://doi.org/10.1101/801696)
# Number of PubMed citations for ISRIB paper == 199, but in random sample of 20 (~10%), only 2 used for some sort of re-analysis. Will assume a 1/10 use-rate.
n = 199 * (2/20)
alpha = 0.05
from math import ceil
for x in range(0, ceil(n)):
    alpha = alpha / 2
alpha


# Make ribo-seq vs rna-seq plots with hightlights
# Tm
rp_plot(
    merged_data_split_XP,
    'tm_ribo_log2FC',
    'tm_rna_log2FC',
    'tm_ribo_padj',
    'tm_rna_padj',
    'Tm',
    'Untr',
    merged_data_XP[['tm_log2FC','tm_padj']],
    'tm_padj',
    'Tm_vs_Untr_deseq.png',
    isr,
    uorf_targets,
    down_strict_XP)

# Tm + ISRIB
rp_plot(
    merged_data_split_XP,
    'tmisrib_ribo_log2FC',
    'tmisrib_rna_log2FC',
    'tmisrib_ribo_padj',
    'tmisrib_rna_padj',
    'Tm + ISRIB',
    'Untr',
    merged_data_XP[['tmisrib_log2FC','tmisrib_padj']],
    'tmisrib_padj',
    'TmISRIB_vs_Untr_deseq.png',
    isr,
    uorf_targets,
    down_strict_XP,
    label_up=False,
    label_down=False,
    label_down_custom=False)

# ISRIB
rp_plot(
    merged_data_split_XP,
    'isrib_ribo_log2FC',
    'isrib_rna_log2FC',
    'isrib_ribo_padj',
    'isrib_rna_padj',
    'ISRIB',
    'Untr',
    merged_data_XP[['isrib_log2FC','isrib_padj']],
    'isrib_padj',
    'ISRIB_vs_Untr_deseq.png',
    isr,
    uorf_targets,
    down_strict_XP,
    label_up=False,
    label_down=False,
    label_down_custom=False)


# Plot TEs
te_plot(
    data=data_de_plot,
    list_up=isr,
    list_down=down_strict_XP,
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    plot_all=True)

"""Comparing outputs
"""

"""XPRESSpipe vs Original w/ DESeq1
"""
# Tm
tm_fc_XP = merged_data_split_XP.loc[(merged_data_split_XP['tm_ribo_log2FC'] >= 1) | (merged_data_split_XP['tm_ribo_log2FC'] <= -1)].index.tolist()
tm_te_XP = merged_data_XP.loc[merged_data_XP['tm_padj'] < 0.1].index.tolist()
tm_common_XP = list(set(tm_fc_XP).intersection(tm_te_XP))

venn_plot(
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    figure_title='Tm',
    gene_list_A=tm_common_XP,
    gene_list_B=ingolia_tm,
    label_A='This paper',
    label_B='Original paper',
    color_A='#d95f02',
    color_B='#1b9e77',
    output_name='tm_deseq2XP_venn.png',
    dpi=3600)

# Tm + ISRIB
tmisrib_fc_XP = merged_data_split_XP.loc[(merged_data_split_XP['tmisrib_ribo_log2FC'] >= 1) | (merged_data_split_XP['tmisrib_ribo_log2FC'] <= -1)].index.tolist()
tmisrib_te_XP = merged_data_XP.loc[merged_data_XP['tmisrib_padj'] < 0.1].index.tolist()
tmisrib_common_XP = list(set(tmisrib_fc_XP).intersection(tmisrib_te_XP))

venn_plot(
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    figure_title='Tm + ISRIB',
    gene_list_A=tmisrib_common_XP,
    gene_list_B=ingolia_tmisrib,
    label_A='This paper',
    label_B='Original paper',
    color_A='#d95f02',
    color_B='#1b9e77',
    output_name='tm_isrib_deseq2XP_venn.png',
    dpi=3600)

# ISRIB
isrib_fc_XP = merged_data_split_XP.loc[(merged_data_split_XP['isrib_ribo_log2FC'] >= 1) | (merged_data_split_XP['isrib_ribo_log2FC'] <= -1)].index.tolist()
isrib_te_XP = merged_data_XP.loc[merged_data_XP['isrib_padj'] < 0.1].index.tolist()
isrib_common_XP = list(set(isrib_fc_XP).intersection(isrib_te_XP))

venn_plot(
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    figure_title='ISRIB',
    gene_list_A=isrib_common_XP,
    gene_list_B=ingolia_isrib,
    label_A='This paper',
    label_B='Original paper',
    color_A='#d95f02',
    color_B='#1b9e77',
    output_name='isrib_deseq2XP_venn.png',
    dpi=3600)

"""Original data DESeq2 vs Original data DESeq1
"""
dir_og2 = '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/original_data_deseq2/'

merged_data_OG2, merged_data_split_OG2, check_data_OG2 = parse_de_data(
    directory=dir_og2,
    check_file='ingolia_counts_deduplicated.txt')

# Tm
tm_fc_og2 = merged_data_split_OG2.loc[(merged_data_split_OG2['tm_ribo_log2FC'] >= 1) | (merged_data_split_OG2['tm_ribo_log2FC'] <= -1)].index.tolist()
tm_te_og2 = merged_data_OG2.loc[merged_data_OG2['tm_padj'] < 0.1].index.tolist()
tm_common_og2 = list(set(tm_fc_og2).intersection(tm_te_og2))

venn_plot(
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    figure_title='Tm',
    gene_list_A=tm_common_og2,
    gene_list_B=ingolia_tm,
    label_A='DESeq2',
    label_B='DESeq1',
    color_A='#7570b3',
    color_B='#1b9e77',
    output_name='tm_deseq2OG_venn.png',
    dpi=3600)

# Tm + ISRIB
tmisrib_fc_og2 = merged_data_split_OG2.loc[(merged_data_split_OG2['tmisrib_ribo_log2FC'] >= 1) | (merged_data_split_OG2['tmisrib_ribo_log2FC'] <= -1)].index.tolist()
tmisrib_te_og2 = merged_data_OG2.loc[merged_data_OG2['tmisrib_padj'] < 0.1].index.tolist()
tmisrib_common_og2 = list(set(tmisrib_fc_og2).intersection(tmisrib_te_og2))

venn_plot(
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    figure_title='Tm + ISRIB',
    gene_list_A=tmisrib_common_og2,
    gene_list_B=ingolia_tmisrib,
    label_A='DESeq2',
    label_B='DESeq1',
    color_A='#7570b3',
    color_B='#1b9e77',
    output_name='tm_isrib_deseq2OG_venn.png',
    dpi=3600)

# ISRIB
isrib_fc_og2 = merged_data_split_OG2.loc[(merged_data_split_OG2['isrib_ribo_log2FC'] >= 1) | (merged_data_split_OG2['isrib_ribo_log2FC'] <= -1)].index.tolist()
isrib_te_og2 = merged_data_OG2.loc[merged_data_OG2['isrib_padj'] < 0.1].index.tolist()
isrib_common_og2 = list(set(isrib_fc_og2).intersection(isrib_te_og2))

venn_plot(
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    figure_title='ISRIB',
    gene_list_A=isrib_common_og2,
    gene_list_B=ingolia_isrib,
    label_A='DESeq2',
    label_B='DESeq1',
    color_A='#7570b3',
    color_B='#1b9e77',
    output_name='isrib_deseq2OG_venn.png',
    dpi=3600)

"""XPRESSpipe DESeq2 vs Original w/ DESeq2
"""
# Tm
venn_plot(
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    figure_title='Tm',
    gene_list_A=tm_common_XP,
    gene_list_B=tm_common_og2,
    label_A='This paper',
    label_B='Original paper',
    color_A='#d95f02',
    color_B='#7570b3',
    output_name='tm_deseq2both_venn.png',
    dpi=3600)

# Tm + ISRIB
venn_plot(
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    figure_title='Tm + ISRIB',
    gene_list_A=tmisrib_common_XP,
    gene_list_B=tmisrib_common_og2,
    label_A='This paper',
    label_B='Original paper',
    color_A='#d95f02',
    color_B='#7570b3',
    output_name='tm_isrib_deseq2both_venn.png',
    dpi=3600)

# ISRIB
venn_plot(
    directory='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/',
    figure_title='ISRIB',
    gene_list_A=isrib_common_XP,
    gene_list_B=isrib_common_og2,
    label_A='This paper',
    label_B='Original paper',
    color_A='#d95f02',
    color_B='#7570b3',
    output_name='isrib_deseq2both_venn.png',
    dpi=3600)
