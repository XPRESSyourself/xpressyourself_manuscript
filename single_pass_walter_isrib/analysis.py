import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
matplotlib.rcParams['font.sans-serif'] = 'Arial'
import seaborn as sns
from xpresstools import convert_names, rpm, check_samples, threshold
from scipy import stats
from statsmodels.stats.multitest import multipletests
%matplotlib inline

"""
READ IN DATA
"""
# Read in single-pass XPRESSpipe processed read data
data = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/isrib_riboprof_counts_table.tsv',
    sep = '\t',
    index_col = 0)

data = convert_names(
    data,
    '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.95.gtf')

data.columns = data.columns.str.replace('_1_dedupRemoved', '')

data.head()

# Combine lanes
sra_info = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/GSE65778_table.txt',
    sep = '\t')

sra_info.head()

name_dict = pd.Series(sra_info.ingolia_name.values,index=sra_info.Run).to_dict()
data_c = data.copy()
data_c.columns = data_c.columns.to_series().map(name_dict)
data_sum = data_c.groupby(data_c.columns, axis=1).sum()
data_sum.head()

# Read in Ingolia raw counts from Elife supplement table
ingolia = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/elife_ingolia_raw_reads_isrib.txt',
    sep = '\t',
    index_col = 0)

ingolia.head()

# Combine duplicate genes (assuming mappings of different isoforms)
data_dups = data_sum.groupby(level=0).sum()
data_dups_number = data_dups.groupby(level=0).filter(lambda x: len(x) > 1)
len(data_dups_number)

ingolia_dups = ingolia.groupby(level=0).sum()
ingolia_dups_number = ingolia_dups.groupby(level=0).filter(lambda x: len(x) > 1)
len(ingolia_dups_number)

data_threshold = data_dups[data_dups.min(axis=1) > 25]
ingolia_threshold = ingolia_dups[ingolia_dups.min(axis=1) > 25]

data_threshold10 = data_dups[data_dups.min(axis=1) > 10]
ingolia_threshold10 = ingolia_dups[ingolia_dups.min(axis=1) > 10]

# Make dataframes for each dataset for matching genes
jordan_genes = data_threshold.index.tolist()
ingolia_genes = ingolia_threshold.index.tolist()
len(jordan_genes)
len(ingolia_genes)

common_genes = list(set(jordan_genes).intersection(ingolia_genes))
len(common_genes)

data_common = data_threshold.reindex(index = common_genes)
ingolia_common = ingolia_threshold.reindex(index = common_genes)

# Check all genes are common
jordan_commongenes = data_common.index.tolist()
ingolia_commongenes = ingolia_common.index.tolist()
common_all = list(set(jordan_commongenes).intersection(ingolia_commongenes))

if len(common_all) == len(common_genes):
    print('Data looks good')
else:
    print('Looks like the datasets are still not equal')

# Normalize datasets
data_rpm = rpm(data_threshold)
ingolia_rpm = rpm(ingolia_threshold)

data_rpm10 = rpm(data_threshold10)
ingolia_rpm10 = rpm(ingolia_threshold10)

"""
SUPPLEMENTAL FIGURE 2A
"""
# Run correlations between sample alignments
fig, axes = plt.subplots(
    nrows = 4,
    ncols = 4,
    figsize = (20, 20),
    subplot_kw = {
        'facecolor':'none'},
    sharex=True, sharey=True) # Create shared axis for cleanliness
plt.subplots_adjust(
    bottom = 0.1)
plt.yticks([1,2,3,4,5]) # Limit axis labels to ints
plt.xticks([1,2,3,4])

file_number = 0
file_list = [
    'ribo_untr_a',
    'ribo_untr_b',
    'untr_a_hek',
    'untr_b_hek',
    'ribo_tm_a',
    'ribo_tm_b',
    'tm_a_hek',
    'tm_b_hek',
    'ribo_tmisrib_a',
    'ribo_tmisrib_b',
    'tmisrib_a_hek',
    'tmisrib_b_hek',
    'ribo_isrib_a',
    'ribo_isrib_b',
    'isrib_a_hek',
    'isrib_b_hek'] # Designate sample order

for x in file_list:

    # Get data as array-like for samples being compared
    data_c1 = data_common.copy()
    data_c2 = ingolia_common.copy()

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
    axes[ax_y, ax_x].axhline(1, ls='-', color='black') # Create axis lines
    axes[ax_y, ax_x].axvline(1, ls='-', color='black', ymax=0.88)
    file_number += 1 # Plot counter

# Create shared row/column titles
cols = ['RPF RepA','RPF RepB','mRNA RepA','mRNA RepB']
for ax, col in zip(axes[0], cols):
    ax.set_xlabel(col, fontsize=24)
    ax.xaxis.set_label_position('top')

rows = ['Untreated','Tm','Tm + ISRIB','ISRIB']
for ax, row in zip(axes[:,0], rows):
    ax.set_ylabel(row, fontsize=24)

fig.savefig(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/plots/figures/external_correlations_summary.png',
    dpi = 600,
    bbox_inches = 'tight')

"""
SUPPLEMENTAL FIGURE 2B
"""
# Run correlations between sample alignments
fig, axes = plt.subplots(
    nrows = 2,
    ncols = 4,
    figsize = (20, 10),
    subplot_kw = {
        'facecolor':'none'},
    sharex=True, sharey=True) # Create shared axis for cleanliness
plt.subplots_adjust(
    bottom = 0.1)
plt.yticks([1,2,3,4,5]) # Limit axis labels to ints
plt.xticks([1,2,3,4,5])

x = 0
file_number = 0
file_list = [
    'ribo_untr_a',
    'ribo_untr_b',
    'untr_a_hek',
    'untr_b_hek',
    'ribo_tm_a',
    'ribo_tm_b',
    'tm_a_hek',
    'tm_b_hek',
    'ribo_tmisrib_a',
    'ribo_tmisrib_b',
    'tmisrib_a_hek',
    'tmisrib_b_hek',
    'ribo_isrib_a',
    'ribo_isrib_b',
    'isrib_a_hek',
    'isrib_b_hek'] # Designate sample order

for y in range(int(len(file_list)/2)):

    # Get data as array-like for samples being compared
    data_c1 = data_threshold.copy()

    sample_a = data_c1[file_list[x]].values.tolist()
    sample_a = [x + 1 for x in sample_a]
    sample_a = np.array(sample_a).astype(np.float)
    sample_a = np.ndarray.tolist(sample_a)

    sample_b = data_c1[file_list[x+1]].values.tolist()
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
    else:
        print('oops')

    # Format p value
    if p_value.astype('float') < 0.001:
        p_val = '< 0.001'
    else:
        p_val = round(p_value.astype('float'), 4).astype('str')


    if rho >= 0.99:
        rh = '> 0.99'
    else:
        rh = round(rho.astype('float'), 2).astype('float').astype('str')

    # Plot data
    axes[ax_y, ax_x].scatter(np.log10(sample_a), np.log10(sample_b), s=1,c='black')
    axes[ax_y, ax_x].set_title('R ' + rh + '\nP ' + p_val, y=0.1, x=0.9, fontsize=16) # Format titles
    axes[ax_y, ax_x].axhline(1, ls='-', color='black') # Create axis lines
    axes[ax_y, ax_x].axvline(1, ls='-', color='black', ymax=0.88)
    file_number += 1 # Plot counter
    x += 2

# Create shared row/column titles
cols = ['RPF Untreated','mRNA Untreated','RPF Tm','mRNA Tm']
for ax, col in zip(axes[0], cols):
    ax.set_xlabel(col, fontsize=24)
    ax.xaxis.set_label_position('top')

cols = ['RPF Tm + ISRIB','mRNA Tm + ISRIB','RPF ISRIB','mRNA ISRIB']
for ax, col in zip(axes[1], cols):
    ax.set_xlabel(col, fontsize=24)
    ax.xaxis.set_label_position('top')

fig.savefig(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/plots/figures/internal_correlations_summary.png',
    dpi = 600,
    bbox_inches = 'tight')

"""
SUPPLEMENTAL FIGURE 2C
My alignments appear to have less outliers, indicating a better read processing overall
"""
# Get reference line
ref_line = np.log10(data_rpm + 1).mean().mean()
fig, ax = plt.subplots(1,1, figsize=(12.5,5))
plt.grid(False)
names = [
    'ribo_untr_a',
    'ribo_untr_b',
    'untr_a_hek',
    'untr_b_hek',
    'ribo_tm_a',
    'ribo_tm_b',
    'tm_a_hek',
    'tm_b_hek',
    'ribo_tmisrib_a',
    'ribo_tmisrib_b',
    'tmisrib_a_hek',
    'tmisrib_b_hek',
    'ribo_isrib_a',
    'ribo_isrib_b',
    'isrib_a_hek',
    'isrib_b_hek']
ax = sns.violinplot(data=np.log10(data_rpm + 1), order = names)
plt.xticks(rotation=45)
ax.set_ylabel('Normalized Counts (RPM)')
ax.set_xlabel('Samples')
ax.axhline(ref_line, ls='--', color='grey')
plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/plots/figures/xpresspipe_run_distributions.png', dpi=1800, bbox_inches='tight')
plt.close

# Check Ingolia
check_samples(np.log10(ingolia_rpm + 1))
sns.violinplot(data=np.log10(ingolia_rpm + 1))

"""
FIGURE 2A
My alignments have less "Super-translators", which is actually better and expected
Not sure why the mRNA levels seem to be more constained on the positive side
"""
# Make Walter ISRIB paper figures for ribosome profiling data
def rp_plot(
    data,
    conditionB_mrna,
    conditionB_ribo,
    conditionA_mrna,
    conditionA_ribo,
    conditionB_name,
    conditionA_name,
    save_fig,
    genes=[]): #File name

    data_c = data.copy()

    # Get samples of interest
    samples = (
        conditionB_mrna
        + conditionB_ribo
        + conditionA_mrna
        + conditionA_ribo)

    data_fig = data_c[samples]

    # Calculate fold changes
    data_fig['conditionB_mrna'] = data_fig[conditionB_mrna[0]] + data_fig[conditionB_mrna[1]]
    data_fig['conditionB_ribo'] = data_fig[conditionB_ribo[0]] + data_fig[conditionB_ribo[1]]
    data_fig['conditionA_mrna'] = data_fig[conditionA_mrna[0]] + data_fig[conditionA_mrna[1]]
    data_fig['conditionA_ribo'] = data_fig[conditionA_ribo[0]] + data_fig[conditionA_ribo[1]]

    data_fig['rna_fc'] = data_fig['conditionB_mrna'] / data_fig['conditionA_mrna']
    data_fig['ribo_fc'] = data_fig['conditionB_ribo'] / data_fig['conditionA_ribo']

    # Calculate p-values
    drop_index = []

    data_fig['p_value'] = '' # Initialize column to store p-values

    # Calculate p-value using One-way ANOVA
    for index, row_data in data_fig.iterrows():
        condB_mrna = [data_fig[conditionB_mrna[0]][index], data_fig[conditionB_mrna[1]][index]]
        condB_ribo = [data_fig[conditionB_ribo[0]][index], data_fig[conditionB_ribo[1]][index]]
        condA_mrna = [data_fig[conditionA_mrna[0]][index], data_fig[conditionA_mrna[1]][index]]
        condA_ribo = [data_fig[conditionA_ribo[0]][index], data_fig[conditionA_ribo[1]][index]]

        # Append p_value to dataframe
        try:
            statistic, p_value = stats.f_oneway(condB_mrna, condB_ribo, condA_mrna, condA_ribo)
            data_fig.loc[index,'p_value'] = p_value
        except:
            drop_index.append(index)
            print(index)

    data_fig = data_fig.drop(labels=drop_index, axis=0) # Drop bad indices

    # Calculate adjusted p via Ben-Hoch FDR
    data_fig['p_adj'] = multipletests(
        data_fig['p_value'].values,
        method = 'fdr_bh')[1]

    # Print interesting ribo fc results
    print('Translation targets')
    print(data_fig[(data_fig['p_adj'] < 0.01) & (data_fig['ribo_fc'] > 2)].index.tolist())
    print(data_fig[(data_fig['p_adj'] < 0.01) & (data_fig['ribo_fc'] < 0.5)].index.tolist())

    print('Transcription targets')
    print(data_fig[(data_fig['p_adj'] < 0.01) & (data_fig['rna_fc'] > 2)].index.tolist())
    print(data_fig[(data_fig['p_adj'] < 0.01) & (data_fig['rna_fc'] < 0.5)].index.tolist())

    # Initialize figure
    fig, ax = plt.subplots(1,1, figsize=(7,7))
    plt.grid(False)
    ax.axhline(1, ls='-', color='black') # Create central axis line
    ax.axvline(1, ls='-', color='black')
    rect = patches.Rectangle((.5,.5),1.5,1.5,linewidth=1.5,edgecolor='gray',facecolor='none') # Create significance zone
    ax.add_patch(rect)
    ax.set_ylabel('Ribo-Seq (' + str(conditionB_name) + '/' + str(conditionA_name) + ')', fontsize=16) # Set axis labels
    ax.set_xlabel('mRNA-Seq (' + str(conditionB_name) + '/' + str(conditionA_name) + ')', fontsize=16)
    ax.set_xlim(.25,4.5) # Set axis limits
    ax.set_ylim(.25,4.5)
    x = [.250,.500,1,2,4] # Set axis spacing
    ticks = ["1/4","1/2","1","2","4"]
    ax.set_facecolor("#FFFFFF") # Set background color

    # Set X and Y axis scales
    ax.set_xscale('log', basex=2)
    ax.set_xticks(x)
    ax.set_xticklabels(ticks, fontsize=16)
    ax.set_xticklabels([""]*len(x), minor=True)
    ax.set_yscale('log', basey=2)
    ax.set_yticks(x)
    ax.set_yticklabels(ticks, fontsize=16)
    ax.set_yticklabels([""]*len(x), minor=True)

    # Prep data for plotting
    ribo_all = data_fig[['ribo_fc']].sum(axis=1).values.tolist()
    ribo_all = np.array(ribo_all).astype(np.float)
    ribo_all = np.ndarray.tolist(ribo_all)

    rna_all = data_fig[['rna_fc']].sum(axis=1).values.tolist()
    rna_all = np.array(rna_all).astype(np.float)
    rna_all = np.ndarray.tolist(rna_all)

    # Prep significant hits for plotting
    data_sig = data_fig[data_fig['p_adj'] < 0.01]

    ribo_sig = data_sig[['ribo_fc']].sum(axis=1).values.tolist()
    ribo_sig = np.array(ribo_sig).astype(np.float)
    ribo_sig = np.ndarray.tolist(ribo_sig)

    rna_sig = data_sig[['rna_fc']].sum(axis=1).values.tolist()
    rna_sig = np.array(rna_sig).astype(np.float)
    rna_sig = np.ndarray.tolist(rna_sig)

    isr = ['ATF5','ATF4']
    data_fig_isr = data_fig.loc[isr]
    ribo_isr = data_fig_isr[['ribo_fc']].sum(axis=1).values.tolist()
    ribo_isr = np.array(ribo_isr).astype(np.float)
    ribo_isr = np.ndarray.tolist(ribo_isr)

    rna_isr = data_fig_isr[['rna_fc']].sum(axis=1).values.tolist()
    rna_isr = np.array(rna_isr).astype(np.float)
    rna_isr = np.ndarray.tolist(rna_isr)

    other = ['CCNG2', 'ENPP5', 'FGFRL1', 'GALNS', 'KIAA1841', 'MYO5B', 'NGEF', 'NOTCH3', 'POMGNT1', 'SMURF2']
    data_fig_other = data_fig.loc[other]
    ribo_other = data_fig_other[['ribo_fc']].sum(axis=1).values.tolist()
    ribo_other = np.array(ribo_other).astype(np.float)
    ribo_other = np.ndarray.tolist(ribo_other)

    rna_other = data_fig_other[['rna_fc']].sum(axis=1).values.tolist()
    rna_other = np.array(rna_other).astype(np.float)
    rna_other = np.ndarray.tolist(rna_other)

    # Label points of interest
    for index, row in data_fig_isr.iterrows():
        ax.text(row[12] + 0.05, row[13] + 0.05, str(index), horizontalalignment='left', size='large', color='#be00be', weight='semibold')

    #Plot data
    ax.scatter(rna_all, ribo_all, s=5,c='gray')
    ax.scatter(rna_sig, ribo_sig, s=5,c='black',alpha=0.8)
    ax.scatter(rna_isr, ribo_isr, s=40,c='#be00be')
    ax.scatter(rna_other, ribo_other, s=40,c='green')

    plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/plots/elife_figures/' + str(save_fig), dpi=1800, bbox_inches='tight')
    plt.close()

    return data_fig

# Tm treatment vs Untreated
untr_mrna = ['untr_a_hek', 'untr_b_hek']
untr_ribo = ['ribo_untr_a', 'ribo_untr_b']
tm_mrna = ['tm_a_hek', 'tm_b_hek']
tm_ribo = ['ribo_tm_a', 'ribo_tm_b']
tmisrib_mrna = ['tmisrib_a_hek', 'tmisrib_b_hek']
tmisrib_ribo = ['ribo_tmisrib_a', 'ribo_tmisrib_b']
isrib_mrna = ['isrib_a_hek', 'isrib_b_hek']
isrib_ribo = ['ribo_isrib_a', 'ribo_isrib_b']

os.system('mkdir plots/elife_figures')

"""Tm vs Untr"""
tm_data = rp_plot(data_rpm, tm_mrna, tm_ribo, untr_mrna, untr_ribo, 'Tm', 'Untr', 'Tm_vs_Untr_jordan.png')

#tm_data10 = rp_plot(data_rpm10, tm_mrna, tm_ribo, untr_mrna, untr_ribo, 'Tm', 'Untr', 'Tm_vs_Untr_threshold10_jordan.png')
#tm_data10.loc[['ATF5','ATF4','DDIT3','PPP1R15A']]

tm_up = ['ATF4', 'ATF5', 'ISY1', 'PMAIP1', 'PNRC2', 'TMSB10', 'ZNF155']
# ISY1 = mRNA splicing
# PMAIP1 = diseases include Ischemic Neuropathy, promotes mitochondrial membrane changes and efflux of apoptogenic proteins from the mitochondria
# PNRC2 = Plays a role in glucocorticoid receptor-mediated mRNA degradation by interacting with the glucocorticoid receptor NR3C1 in a ligand-dependent manner when it is bound to the 5' UTR of target mRNAs and recruiting the RNA helicase UPF1 and the mRNA-decapping enzyme DCP1A, leading to RNA decay; may explain why we see a general slight decrease in mRNA levels in Tm conditions
# ZNF155 = May be involved in transcriptional regulation

tm_down = ['AVEN', 'BMP6', 'CABLES2', 'CCNG2', 'DDN', 'EHBP1L1', 'ENPP5', 'FGFRL1', 'GALNS', 'HSPA2', 'KIAA1841', 'LIMD1', 'LOXL2', 'LY75-CD302', 'MYB', 'MYO5B', 'NGEF', 'NOTCH3', 'POMGNT1', 'RASGRF2', 'SMURF2']

tm_down_neuro_strict = ['CCNG2', 'ENPP5', 'FGFRL1', 'GALNS', 'KIAA1841', 'MYO5B', 'NGEF', 'NOTCH3', 'POMGNT1', 'SMURF2']
tm_down_neuro = ['CCNG2', 'ENPP5', 'FGFRL1', 'GALNS', 'KIAA1841', 'MYO5B', 'NGEF', 'NOTCH3','POMGNT1','SMURF2']
tm_down_transc = ['BMP6', 'DDN', 'LIMD1', 'LOXL2', 'MYB', 'SMURF2']
# AVEN = associated with Schizoid Personality Disorder
# ENPP5 = Studies in rat suggest the encoded protein may play a role in neuronal cell communications
# FGFRL1 = behavior/neurological phenotype
# CCNG2 = overexpressed in brain, tightly regulated in cell cycle
# KIAA1841 = ID'd in screen of proteins expressed in brain
# MYO5B = ID'd in screen for large proteins expressed in brain, utilize energy from ATP hydrolysis to generate mechanical force, required for the recycling of transferrin receptor back to the plasma membrane through an endocytotic recycling compartment in nonpolarized, but not polarized, cells.
# NGEF = axon guidance regulating ephrin-induced growth cone collapse and dendritic spine morphogenesis
# ***NOTCH3 = intercellular signalling pathway that plays a key role in neural development, mutations in NOTCH3 have been identified as the underlying cause of cerebral autosomal dominant arteriopathy with subcortical infarcts and leukoencephalopathy
# POMGNT1 = Mutations in this gene may be associated with muscle-eye-brain disease and several congenital muscular dystrophies
# RASGRF2 = coordinating the signaling of distinct mitogen-activated protein kinase pathways, synaptic plasticity by contributing to the induction of long term potentiation

# BMP6 = Ligands of this family bind various TGF-beta receptors leading to recruitment and activation of SMAD family transcription factors that regulate gene expression
# DDN = Gene Ontology (GO) annotations related to this gene include RNA polymerase II proximal promoter sequence-specific DNA binding and transcription factor activity; over-expressed in brain
# LIMD1 = Gene Ontology (GO) annotations related to this gene include transcription corepressor activity.
# LOXL2 = Gene Ontology (GO) annotations related to this gene include chromatin binding and electron transfer activity; Also mediates deamination of methylated TAF10, a member of the transcription factor IID (TFIID) complex, which induces release of TAF10 from promoters, leading to inhibition of TFIID-dependent transcription; SNAI1 recruits LOXL2 to pericentromeric regions to oxidize histone H3 and repress transcription which leads to release of heterochromatin component CBX5/HP1A, enabling chromatin reorganization and acquisition of mesenchymal traits; a susceptibility gene to intracranial aneurysms
# MYB = contains three HTH DNA-binding domains that functions as a transcription regulator; considered an oncogene, related to Notch signaling,
# SMURF2 = Interacts with SMAD proteins; functions in regulation of neuronal and planar cell polarity, induction of senescence, and tumor suppression

# HSPA2 = Heat Shock Protein Family A (Hsp70) Member 2;  including protection of the proteome from stress, folding and transport of newly synthesized polypeptides, activation of proteolysis of misfolded proteins and the formation and dissociation of protein complexes. Plays a pivotal role in the protein quality control system, ensuring the correct folding of proteins, the re-folding of misfolded proteins and controlling the targeting of proteins for subsequent degradation; normally over-expressed in brain
# GALNS = Involved in breaking down and recycling molecules in lysosome

"""Tm + ISRIB vs Untr"""
tmisrib_data = rp_plot(data_rpm, tmisrib_mrna, tmisrib_ribo, untr_mrna, untr_ribo, 'Tm+ISRIB', 'Untr', 'TmISRIB_vs_Untr_jordan.png')
tmisrib_up = ['CHD5', 'STARD9']
# CHD5 = Neuron-specific protein that may function in chromatin remodeling and gene transcription. This gene is a potential tumor suppressor gene that may play a role in the development of neuroblastoma. Diseases associated with CHD5 include Neuroblastoma. Overexpressed in brain

tmisrib_down = ['CHAC1']
# CHAC1 = Shown to promote neuronal differentiation by deglycination of the Notch receptor, which prevents receptor maturation and inhibits Notch signaling. This protein may also play a role in the unfolded protein response, and in regulation of glutathione levels and oxidative balance in the cell. Elevated expression of this gene may indicate increased risk of cancer recurrence among breast and ovarian cancer patients.

"""ISRIB vs Untr"""
isrib_data = rp_plot(data_rpm, isrib_mrna, isrib_ribo, untr_mrna, untr_ribo, 'ISRIB', 'Untr', 'ISRIB_vs_Untr_jordan.png')
isrib_up = []
isrib_down = ['AC008758.1', 'CA13', 'EHBP1L1', 'ENPP5', 'ZBTB37']
# ***CA13 = related pathways are Metabolism and Nitrogen metabolism
# EHBP1L1 = May act as Rab effector protein and play a role in vesicle trafficking; several Rab related proteins were regulated in Tm
# ENPP5 = DOWNREGULATED IN Tm;
# ZBTB37 = May be involved in transcriptional regulation

# Ingolia plots
rp_plot(ingolia_rpm, tm_mrna, tm_ribo, untr_mrna, untr_ribo, 'Tm', 'Untr', 'Tm_vs_Untr_ingolia.png')
tm_up_ingolia = ['AGAP4', 'AK097143', 'ATF4', 'BCL2L11', 'CPS1', 'DDIT3', 'FAM13B', 'FAM198A', 'GOLGA8A', 'GOLGA8B', 'HOXB2', 'NGRN', 'OAZ1', 'PNRC2', 'PTP4A1', 'SAT1', 'UCP2', 'ZNF432', 'ZSCAN5A']
tm_down_ingolia = ['KIAA1841', 'MYO5B']

rp_plot(ingolia_rpm, tmisrib_mrna, tmisrib_ribo, untr_mrna, untr_ribo, 'Tm+ISRIB', 'Untr', 'TmISRIB_vs_Untr_ingolia.png')
tmisrib_up_ingolia = ['AGAP4', 'AGAP5', 'AGAP6', 'AGAP7', 'AGAP8', 'AGAP9', 'AK097143', 'AK304826', 'GOLGA8A', 'GOLGA8B', 'LOC100132247', 'LRP1B', 'MAPK8IP3', 'NPIPL3', 'PKD1']
tmisrib_down_ingolia = ['ZNF449']

rp_plot(ingolia_rpm, isrib_mrna, isrib_ribo, untr_mrna, untr_ribo, 'ISRIB', 'Untr', 'ISRIB_vs_Untr_ingolia.png')
isrib_up_ingolia = []
isrib_down_ingolia = []

"""
FIGURE 2B
"""
# Track neuronal related genes across conditions to see if they are significantly changed between conditions
tm_data_set = tm_data.reindex(labels=tm_down_neuro_strict)
tmisrib_data_set = tmisrib_data.reindex(labels=tm_down_neuro_strict)
isrib_data_set = isrib_data.reindex(labels=tm_down_neuro_strict)
tm_data_set['Tm'] = np.log2(tm_data_set['ribo_fc'] / tm_data_set['rna_fc'])
tmisrib_data_set['Tm + ISRIB'] = np.log2(tmisrib_data_set['ribo_fc'] / tmisrib_data_set['rna_fc'])
isrib_data_set['ISRIB'] = np.log2(isrib_data_set['ribo_fc'] / isrib_data_set['rna_fc'])

dataset = pd.concat([tm_data_set[['Tm']], tmisrib_data_set[['Tm + ISRIB']], isrib_data_set[['ISRIB']]], axis=1, sort=False)
dataset['diff'] = tmisrib_data_set['Tm + ISRIB'] - tm_data_set['Tm']
2**dataset['diff'].min()
2**dataset['diff'].max()

dataset

ax = dataset.T.plot.line(figsize=(5,7.5), grid = False)
ax.axvline(1, ls='-', color='white')
ax.set_ylabel('log$_2$TE')
plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/plots/elife_figures/translation_efficiencies_interest.png', dpi=1800, bbox_inches='tight')


"""
FIGURE 2C
"""
# Analyze bigwig files for genes of interest to make sure they are reliable hits
