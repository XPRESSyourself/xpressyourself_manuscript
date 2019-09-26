import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
matplotlib.rcParams['font.sans-serif'] = 'Arial'
import seaborn as sns
from xpressplot import rpm

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

    data_fig = data[[y_fc, x_fc, y_p, x_p]].copy()

    data_fig[[y_fc, x_fc]] = data_fig[[y_fc, x_fc]].replace([np.inf, -np.inf], np.nan)
    data_fig = data_fig.dropna(subset=[y_fc, x_fc])

    # Initialize figure
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
    #ax.set_xscale('log', basex=2)
    ax.set_xticks(x)
    ax.set_xticklabels(ticks, fontsize=16)
    ax.set_xticklabels([""]*len(x), minor=True)
    #ax.set_yscale('log', basey=2)
    ax.set_yticks(x)
    ax.set_yticklabels(ticks, fontsize=16)
    ax.set_yticklabels([""]*len(x), minor=True)

    # Prep data for plotting
    ribo_all = data_fig[[y_fc]].values.tolist()
    ribo_all = np.array(ribo_all).astype(np.float)
    ribo_all = np.ndarray.tolist(ribo_all)

    rna_all = data_fig[[x_fc]].values.tolist()
    rna_all = np.array(rna_all).astype(np.float)
    rna_all = np.ndarray.tolist(rna_all)

    # Prep significant hits for plotting
    sig_list = te_data.loc[te_data[te_p] < 0.1].index.tolist()
    data_sig = data_fig.reindex(sig_list)

    ribo_sig = data_sig[[y_fc]].values.tolist()
    ribo_sig = np.array(ribo_sig).astype(np.float)
    ribo_sig = np.ndarray.tolist(ribo_sig)

    rna_sig = data_sig[[x_fc]].values.tolist()
    rna_sig = np.array(rna_sig).astype(np.float)
    rna_sig = np.ndarray.tolist(rna_sig)

    if list_up != None:
        data_fig_up = data_fig.reindex(list_up)
        ribo_up = data_fig_up[[y_fc]].values.tolist()
        ribo_up = np.array(ribo_up).astype(np.float)
        ribo_up = np.ndarray.tolist(ribo_up)

        rna_up = data_fig_up[[x_fc]].values.tolist()
        rna_up = np.array(rna_up).astype(np.float)
        rna_up = np.ndarray.tolist(rna_up)

    if list_down != None:
        data_fig_down = data_fig.reindex(list_down)
        ribo_down = data_fig_down[[y_fc]].values.tolist()
        ribo_down = np.array(ribo_down).astype(np.float)
        ribo_down = np.ndarray.tolist(ribo_down)

        rna_down = data_fig_down[[x_fc]].values.tolist()
        rna_down = np.array(rna_down).astype(np.float)
        rna_down = np.ndarray.tolist(rna_down)

    if list_down_custom != None:
        data_fig_down_custom = data_fig.reindex(list_down_custom)
        ribo_down_custom = data_fig_down_custom[[y_fc]].values.tolist()
        ribo_down_custom = np.array(ribo_down_custom).astype(np.float)
        ribo_down_custom = np.ndarray.tolist(ribo_down_custom)

        rna_down_custom = data_fig_down_custom[[x_fc]].values.tolist()
        rna_down_custom = np.array(rna_down_custom).astype(np.float)
        rna_down_custom = np.ndarray.tolist(rna_down_custom)

    #Plot data
    ax.scatter(rna_all, ribo_all, s=2.5,c='gray',alpha=0.5)
    ax.scatter(rna_sig, ribo_sig, s=5,c='black',alpha=1)

    if list_down != None:
        ax.scatter(rna_down, ribo_down, s=80,c='#1b9e77',alpha=1)

        for index, row in data_fig.iterrows():
            if index in list_down:
                if index in sig_list:
                    ax.scatter(row[1], row[0], s=20,c='black',alpha=1)
                else:
                    ax.scatter(row[1], row[0], s=20,c='gray',alpha=1)
            else:
                pass

        if label_down == True:
            for index, row in data_fig_down.iterrows():
                if index == 'PNRC2' or index == 'SAT1':
                    ax.text(row[1] + 0.1, row[0] - 0.07, str(index), horizontalalignment='left', size='medium', color='#117a5b', weight='semibold')
                else:
                    ax.text(row[1] - 0.1, row[0] - 0.07, str(index), horizontalalignment='right', size='medium', color='#117a5b', weight='semibold')

    if list_down_custom != None:
        ax.scatter(rna_down_custom, ribo_down_custom, s=80,c='#d95f02',alpha=1)

        for index, row in data_fig.iterrows():
            if index in list_down_custom:
                if index in sig_list:
                    ax.scatter(row[1], row[0], s=20,c='black',alpha=1)
                else:
                    ax.scatter(row[1], row[0], s=20,c='gray',alpha=1)
            else:
                pass

        if label_down_custom == True:
            for index, row in data_fig_down_custom.iterrows():
                if index == 'RPS15A' or index == 'HIST2H3D' or index == 'RPL27':
                    ax.text(row[1] + 0.1, row[0] - 0.07, str(index), horizontalalignment='left', size='medium', color='#c25502', weight='semibold')
                elif index == 'HSPA8':
                    ax.text(row[1] - 0.05, row[0] - 0.2, str(index), horizontalalignment='right', size='medium', color='#c25502', weight='semibold')
                else:
                    ax.text(row[1] - 0.1, row[0] - 0.07, str(index), horizontalalignment='right', size='medium', color='#c25502', weight='semibold')

    if list_up != None:
        ax.scatter(rna_up, ribo_up, s=80,c='#7570b3',alpha=1)

        for index, row in data_fig.iterrows():
            if index in list_up:
                if index in sig_list:
                    ax.scatter(row[1], row[0], s=20,c='black',alpha=1)
                else:
                    ax.scatter(row[1], row[0], s=20,c='gray',alpha=1)
            else:
                pass

        if label_up == True:
            for index, row in data_fig_up.iterrows():
                if index == 'DDIT3':
                    ax.text(row[1] - 0.1, row[0] - 0.07, str(index), horizontalalignment='right', size='medium', color='#544db3', weight='semibold')
                else:
                    ax.text(row[1] + 0.1, row[0] - 0.07, str(index), horizontalalignment='left', size='medium', color='#544db3', weight='semibold')

    plt.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/' + str(title),
        dpi=3600,
        bbox_inches='tight')
    #plt.show()
    plt.close()

# Import threshold counts data for any checking
dir = '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/'

check_data = pd.read_csv(str(dir) + 'ISRIBxpresspipe_thresholded_counts.tsv', sep='\t', index_col=0)
check_list = check_data.index.tolist()

# Import DESeq2 TE data
tm_data = pd.read_csv(str(dir) + 'tm_ISRIBxpresspipe_processed_counts_diffx.tsv', sep='\t')
tmisrib_data = pd.read_csv(str(dir) + 'tmisrib_ISRIBxpresspipe_processed_counts_diffx.tsv', sep='\t')
isrib_data = pd.read_csv(str(dir) + 'isrib_ISRIBxpresspipe_processed_counts_diffx.tsv', sep='\t')

### Clean up this data
tm_data = tm_data.drop(['baseMean','lfcSE','stat'],axis=1)
tm_data.columns = ['tm_log2FC','tm_p','tm_padj']

tmisrib_data = tmisrib_data.drop(['baseMean','lfcSE','stat'],axis=1)
tmisrib_data.columns = ['tmisrib_log2FC','tmisrib_p','tmisrib_padj']

isrib_data = isrib_data.drop(['baseMean','lfcSE','stat'],axis=1)
isrib_data.columns = ['isrib_log2FC','isrib_p','isrib_padj']

merged_data = pd.concat([tm_data, tmisrib_data, isrib_data], axis=1, sort=False)

# Import DESeq2 RNA and RPF data
tm_ribo_data = pd.read_csv(str(dir) + 'tm_ribo_ISRIBxpresspipe_processed_counts_diffx.tsv', sep='\t')
tm_rna_data = pd.read_csv(str(dir) + 'tm_rna_ISRIBxpresspipe_processed_counts_diffx.tsv', sep='\t')

tmisrib_ribo_data = pd.read_csv(str(dir) + 'tmisrib_ribo_ISRIBxpresspipe_processed_counts_diffx.tsv', sep='\t')
tmisrib_rna_data = pd.read_csv(str(dir) + 'tmisrib_rna_ISRIBxpresspipe_processed_counts_diffx.tsv', sep='\t')

isrib_ribo_data = pd.read_csv(str(dir) + 'isrib_ribo_ISRIBxpresspipe_processed_counts_diffx.tsv', sep='\t')
isrib_rna_data = pd.read_csv(str(dir) + 'isrib_rna_ISRIBxpresspipe_processed_counts_diffx.tsv', sep='\t')

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

# Initialize lists of Ingolia highlights
isr = ['ATF4','ATF5','PPP1R15A','DDIT3']
uorf_targets = ['SLC35A4','PTP4A1','UCP2','C7orf31','BCL2L11','PNRC2','SAT1']

"""
Some of the genes in plotting cluster on a single y-axis coordinate
This is due to insufficient RPF counts after RNA thresholding of the data
"""
tm_ribo_data[(tm_ribo_data['tm_ribo_log2FC'] > 2.7) & (tm_ribo_data['tm_ribo_log2FC'] < 2.75)]
check_data.loc['SLC10A5'][['ribo_tm_a', 'ribo_tm_b', 'ribo_untr_a', 'ribo_untr_b']]
""""""

# Check fold changes of targets for comparison
print('ATF4')
print(2**merged_data_split.loc['ATF4']['tm_ribo_log2FC'])
print('ATF5')
print(2**merged_data_split.loc['ATF5']['tm_ribo_log2FC'])
print('PPP1R15A')
print(2**merged_data_split.loc['PPP1R15A']['tm_ribo_log2FC'])
print('DDIT3')
print(2**merged_data_split.loc['DDIT3']['tm_ribo_log2FC'])

"""
TE sanity check
"""
def te_sanity(data, check_data, gene_name):

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

te_sanity(tm_data, check_data, 'ATF4')
te_sanity(tm_data, check_data, 'DDIT3')
te_sanity(tm_data, check_data, 'UCP2')
te_sanity(tm_data, check_data, 'ID2')

te_sanity(tm_data, check_data, 'POMGNT1')

# NDUFA11 pattern really interesting
te_sanity(tm_data, check_data, 'NDUFA11')

# TE analysis
data_de_plot = merged_data.copy()

data_de_plot['Untreated'] = 0
data_de_plot = data_de_plot[['Untreated','tm_log2FC', 'tmisrib_log2FC','isrib_log2FC']]
data_de_plot.columns = ['Untreated','Tm', 'Tm + ISRIB','ISRIB']

down_strict_1 = merged_data.loc[
    (merged_data['tm_log2FC'] - merged_data['tmisrib_log2FC'] <= -1) \
    & (merged_data['tm_padj'] <= 0.001) \
    & (merged_data['tmisrib_padj'] <= 0.001) \
    & ((merged_data_split['tm_ribo_padj'] <= 0.1))]

down_strict_2 = merged_data.loc[
    (merged_data['tm_log2FC'] - merged_data['tmisrib_log2FC'] <= -0.58) \
    & (merged_data_split['tm_ribo_log2FC'] <= -0.9) \
    & (merged_data_split['tm_ribo_padj'] <= 0.1)]


merged_data.loc[set(down_strict_L)]



merged_data_split.loc[set(down_strict_L)]


# All hits have at least 50 counts per RNA
# Thresholds were designed to maximize
down_strict_L = down_strict_1.index.tolist() + down_strict_2.index.tolist()

# Make ribo-seq vs rna-seq plots with hightlights
rp_plot(
    merged_data_split,
    'tm_ribo_log2FC',
    'tm_rna_log2FC',
    'tm_ribo_padj',
    'tm_rna_padj',
    'Tm',
    'Untr',
    tm_data,
    'tm_padj',
    'Tm_vs_Untr_deseq.png',
    isr,
    uorf_targets,
    down_strict_L)

rp_plot(
    merged_data_split,
    'tmisrib_ribo_log2FC',
    'tmisrib_rna_log2FC',
    'tmisrib_ribo_padj',
    'tmisrib_rna_padj',
    'Tm + ISRIB',
    'Untr',
    tmisrib_data,
    'tmisrib_padj',
    'TmISRIB_vs_Untr_deseq.png',
    isr,
    uorf_targets,
    down_strict_L,
    label_up=False,
    label_down=False,
    label_down_custom=False)

rp_plot(
    merged_data_split,
    'isrib_ribo_log2FC',
    'isrib_rna_log2FC',
    'isrib_ribo_padj',
    'isrib_rna_padj',
    'ISRIB',
    'Untr',
    isrib_data,
    'isrib_padj',
    'ISRIB_vs_Untr_deseq.png',
    isr,
    uorf_targets,
    down_strict_L,
    label_up=False,
    label_down=False,
    label_down_custom=False)


# Plot TEs
fig, ax = plt.subplots()
plt.yticks([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
ax.tick_params(axis='x', which='minor',length=0)

ax.set_facecolor('white')
ax.grid(color='grey', axis='y')

data_de_plot.T.plot.line(legend=False, color='lightgrey', linewidth=0.5, ax=ax)

ax.axvline(0.01, ls='-', color='black')
ax.axvline(1, ls='-', color='black')
ax.axvline(2, ls='-', color='black')
ax.axvline(2.99, ls='-', color='black')
ax.axhline(0, ls='-', color='black')
ax.axhline(1, ls='--', color='black')
ax.axhline(-1, ls='--', color='black')

data_de_plot.loc[isr].T.plot.line(legend=False, color='#7570b3', linewidth=1, ax=ax)
data_de_plot.loc[down_strict_L].T.plot.line(legend=False, color='#d95f02', linewidth=1, ax=ax)
ax.set_ylabel(u'log$_2$(ΔTE)')

from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], color='gray', lw=2, label='All'),
                   Line2D([0], [0], color='#7570b3', lw=2, label='ISR'),
                   Line2D([0], [0], color='#d95f02', lw=2, label='Other')]

ax.legend(handles=legend_elements, loc='upper right')
plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/te_analysis_down_strict.png', dpi=3600, bbox_inches='tight')

plt.close()

"""
Comparing outputs
"""
from matplotlib_venn import venn2, venn3

# Tm
tm_fc = merged_data_split.loc[(merged_data_split['tm_ribo_log2FC'] >= 1) | (merged_data_split['tm_ribo_log2FC'] <= -1)].index.tolist()
tm_te = merged_data.loc[merged_data['tm_padj'] < 0.1].index.tolist()
tm_common = list(set(tm_fc).intersection(tm_te))
ingolia_tm = [
    'FAM50B',
    'CHST2',
    'LNP1',
    'AX747756',
    'ACRC',
    'ATF5',
    'ANKRD23',
    'ATF4',
    'C19orf48',
    'TEX14',
    'SLC35A4',
    'DDIT3',
    'RAB35',
    'AK293147',
    'PTP4A1',
    'ANKRD36C',
    'NPIP',
    'ANKMY1',
    'UCP2',
    'NPIPL1',
    'PNRC2',
    'C7orf31',
    'PPP1R15A',
    'IFRD1',
    'BCL2L11',
    'OBSCN',
    'SLC24A1',
    'NPIPL1',
    'NPIPL1',
    'ZNF22',
    'NAIP',
    'NPIP',
    'DNAH1',
    'KIAA1377',
    'FMO4',
    'CCDC150',
    'NPIP',
    'AK310228',
    'ZSCAN5A',
    'FAM13B',
    'NPIPL1',
    'ZNF79',
    'HOXA6',
    'NPIPL2',
    'SLMO1',
    'ARMCX4',
    'AK310228',
    'LOC100288332',
    'SAT1',
    'FLJ00285',
    'AK097143',
    'AX748249',
    'AGAP4',
    'FLJ00322',
    'GOLGA8A',
    'GOLGA8B',
    'ZNF432',
    'FSBP',
    'CPS1',
    'MATK',
    'EPB41L4B',
    'HCN2',
    'GLRB',
    'KIAA1841',
    'MYO5B',
    'EHBP1L1',
    'MAPK8IP1',
    'SH3RF3',
    'DET1',
    'PTPN21',
    'VANGL2',
    'MMP15',
    'MAP3K10',
    'PAQR5',
    'GAL',
    'RPP25',
    'FLJ44955',
    'DQ576756']
plt.title("Tm", fontsize = 28)
c = venn2(
    [set(tm_common), set(ingolia_tm)],
    set_labels = ('This paper', 'Original paper'))
c.get_patch_by_id('10').set_color('#d95f02')
c.get_patch_by_id('01').set_color('#1b9e77')
c.get_patch_by_id('10').set_edgecolor('black')
c.get_patch_by_id('01').set_edgecolor('black')
for text in c.set_labels:
    text.set_fontsize(22)
try:
    for text in c.subset_labels:
        text.set_fontsize(18)
except:
    pass
plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/tm_venn.png', dpi=3600, bbox_inches='tight')


# Tm + ISRIB
tmisrib_fc = merged_data_split.loc[(merged_data_split['tmisrib_ribo_log2FC'] >= 1) | (merged_data_split['tmisrib_ribo_log2FC'] <= -1)].index.tolist()
tmisrib_te = merged_data.loc[merged_data['tmisrib_padj'] < 0.1].index.tolist()
tmisrib_common = list(set(tmisrib_fc).intersection(tmisrib_te))
ingolia_tmisrib = [
    'AK310228',
    'AK310228',
    'NPIP',
    'AK302511',
    'NPIP',
    'LOC100288332',
    'STARD9',
    'PKD1',
    'ANKRD36']
plt.title("Tm + ISRIB", fontsize = 28)
c = venn2([set(tmisrib_common), set(ingolia_tmisrib)], set_labels = ('This paper', 'Original paper'))
c.get_patch_by_id('10').set_color('#d95f02')
c.get_patch_by_id('01').set_color('#1b9e77')
c.get_patch_by_id('10').set_edgecolor('black')
c.get_patch_by_id('01').set_edgecolor('black')
for text in c.set_labels:
    text.set_fontsize(22)
try:
    for text in c.subset_labels:
        text.set_fontsize(18)
except:
    pass
plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/tm_isrib_venn.png', dpi=3600, bbox_inches='tight')

# ISRIB
isrib_fc = merged_data_split.loc[(merged_data_split['isrib_ribo_log2FC'] >= 1) | (merged_data_split['isrib_ribo_log2FC'] <= -1)].index.tolist()
isrib_te = merged_data.loc[merged_data['isrib_padj'] < 0.1].index.tolist()
isrib_common = list(set(isrib_fc).intersection(isrib_te))
ingolia_isrib = [
    'AK310228',
    'AK310228',
    'AK302511',
    'NPIP']
plt.title("ISRIB", fontsize = 28)
c = venn2([set(isrib_common), set(ingolia_isrib)], set_labels = ('This paper', 'Original paper'))
c.get_patch_by_id('10').set_color('#d95f02')
c.get_patch_by_id('01').set_color('#1b9e77')
c.get_patch_by_id('10').set_edgecolor('black')
c.get_patch_by_id('01').set_edgecolor('black')
for text in c.set_labels:
    text.set_fontsize(22)
try:
    for text in c.subset_labels:
        text.set_fontsize(18)
except:
    pass
plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/plots/isrib_venn.png', dpi=3600, bbox_inches='tight')
