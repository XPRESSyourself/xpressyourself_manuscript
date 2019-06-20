import os
import sys
import pandas as pd
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
matplotlib.rcParams['font.sans-serif'] = 'Arial'
import seaborn as sns
from xpressplot import convert_names, rpm, check_samples, threshold
from scipy import stats
from statsmodels.stats.multitest import multipletests
%matplotlib inline
import random
import ast
from urllib.request import Request, urlopen

"""
FUNCTIONS
"""
# Collapse technical replicates and rename samples
def name_map(
    data,
    map):

    name_dict = pd.Series(map.ingolia_name.values,index=map.Run).to_dict()
    data_c = data.copy()
    data_c.columns = data_c.columns.to_series().map(name_dict)
    data_sum = data_c.groupby(data_c.columns, axis=1).sum()

    return data_sum

# Make figure 2A
def make_figure2A(
    data,
    ingolia_data,
    title):

    fig, axes = plt.subplots(
        nrows = 4,
        ncols = 4,
        figsize = (20, 20),
        subplot_kw = {
            'facecolor':'none'},
        sharex=True, sharey=True) # Create shared axis for cleanliness
    plt.subplots_adjust(
        bottom = 0.1)
    plt.yticks([0,1,2,3,4,5]) # Limit axis labels to ints
    plt.xticks([0,1,2,3,4,5])

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
        data_c1 = data.copy()
        data_c2 = ingolia_data.copy()

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

    cols = ['RPF RepA','RPF RepB','mRNA RepA','mRNA RepB']
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=24)
        ax.xaxis.set_label_position('top')

    for ax in axes[3]:
        ax.set_xlabel('log$_1$$_0$(counts)', fontsize=16)
        ax.xaxis.set_label_position('bottom')

    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(counts)', fontsize=16)

    rows = ['Untreated','Tm','Tm + ISRIB','ISRIB']
    for ax, rows in zip(axes[:,3], rows):
        ax.set_ylabel(rows, fontsize=24, labelpad=40, rotation=270)
        ax.yaxis.set_label_position('right')

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

# Make Figure 2B
def make_figure2B(
    data,
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
        data_c1 = data.copy()

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
        axes[ax_y, ax_x].axhline(0, ls='-', color='black', xmin=0.05, xmax=1) # Create axis lines
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymin=0.05, ymax=1)
        file_number += 1 # Plot counter
        x += 2
        print(rh)

    # Create shared row/column titles

    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(counts)', fontsize=16)

    cols = ['RPF Untreated','mRNA Untreated','RPF Tm','mRNA Tm']
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=24)
        ax.xaxis.set_label_position('top')

    cols = ['RPF Tm + ISRIB','mRNA Tm + ISRIB','RPF ISRIB','mRNA ISRIB']
    for ax, col in zip(axes[1], cols):
        ax.set_xlabel(col, fontsize=24)
        ax.xaxis.set_label_position('top')

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

# Make Walter ISRIB paper figures for ribosome profiling data
def rp_plot(
    data,
    conditionB_mrna,
    conditionB_ribo,
    conditionA_mrna,
    conditionA_ribo,
    conditionB_name,
    conditionA_name,
    title,
    list_up,
    list_down,
    genes=[],
    label_up=True,
    label_down=False): #File name

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
    ax.set_xlim(.125,8.5) # Set axis limits
    ax.set_ylim(.125,8.5)
    x = [.125,.250,.500,1,2,4,8] # Set axis spacing
    ticks = ["1/8","1/4","1/2","1","2","4","8"]
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
    data_sig = data_fig[data_fig['p_adj'] < 0.05]

    ribo_sig = data_sig[['ribo_fc']].sum(axis=1).values.tolist()
    ribo_sig = np.array(ribo_sig).astype(np.float)
    ribo_sig = np.ndarray.tolist(ribo_sig)

    rna_sig = data_sig[['rna_fc']].sum(axis=1).values.tolist()
    rna_sig = np.array(rna_sig).astype(np.float)
    rna_sig = np.ndarray.tolist(rna_sig)

    data_fig_up = data_fig.loc[list_up]
    ribo_up = data_fig_up[['ribo_fc']].sum(axis=1).values.tolist()
    ribo_up = np.array(ribo_up).astype(np.float)
    ribo_up = np.ndarray.tolist(ribo_up)

    rna_up = data_fig_up[['rna_fc']].sum(axis=1).values.tolist()
    rna_up = np.array(rna_up).astype(np.float)
    rna_up = np.ndarray.tolist(rna_up)

    data_fig_down = data_fig.loc[list_down]
    ribo_down = data_fig_down[['ribo_fc']].sum(axis=1).values.tolist()
    ribo_down = np.array(ribo_down).astype(np.float)
    ribo_down = np.ndarray.tolist(ribo_down)

    rna_down = data_fig_down[['rna_fc']].sum(axis=1).values.tolist()
    rna_down = np.array(rna_down).astype(np.float)
    rna_down = np.ndarray.tolist(rna_down)

    # Label points of interest
    if label_up == True:
        for index, row in data_fig_up.iterrows():
            ax.text(row[12] + 0.05, row[13] + 0.05, str(index), horizontalalignment='left', size='large', color='#be00be', weight='semibold')

    if label_down == True:
        for index, row in data_fig_down.iterrows():
            ax.text(row[12] + 0.05, row[13] + 0.05, str(index), horizontalalignment='left', size='large', color='green', weight='semibold')

    #Plot data
    ax.scatter(rna_all, ribo_all, s=5,c='gray')
    ax.scatter(rna_sig, ribo_sig, s=5,c='black',alpha=0.8)
    ax.scatter(rna_up, ribo_up, s=40,c='#be00be')
    ax.scatter(rna_down, ribo_down, s=40,c='green')

    plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/plots/' + str(title), dpi=1800, bbox_inches='tight')
    plt.close()

    return data_fig


# Read in gene dictionary
def make_ref(ref_file, gene_list):

    df = pd.read_csv(
        str(ref_file) ,
        sep='\t',
        index_col=2)
    del df.index.name
    df = df[['NCBI gene ID', 'UniProt accession']]
    df['NCBI gene ID'] = df['NCBI gene ID'].astype(pd.Int32Dtype())

    df['index'] = df.index
    df_keep = df.drop_duplicates(subset='index', keep='first')
    df_genes = df_keep.reindex(gene_list)
    df_genes = df_genes.drop('index', axis=1)

    return df_genes

def run_scraper(ref):
    hot_terms_1 = [
        'brain',
        'neuron',
        'central nervous system',
        'synaptic',
        'sympathetic',
        'adrenergic',
        'cerebellum',
        'cerebral',
        'dendritic',
        'astrocyte',
        'excitatory',
        'substantia nigra',
        'red nucleus',
        'hippocampus',
        'spinal cord',
        'motor neuron',
        'cerebral cortex',
        'forebrain',
        'nerve',
        'neuro',
        'neurological',
        'postsynaptic',
        'synaptic',
        'neurotransmitter',
        'neurotoxic',
        'neuron migration',
        'neuron projection',
        'neurodegenerative',
        'psychomotor',
        'psych',
        'noradrenergic',
        'sympathetic',
        'synaptic',
        'epilepsy',
        'alzheimer',
        'encephalopathy',
        'schizophrenia',
        'encephalomyopathic',
        'parkinson',
        'retardation',
        'autistic',
        'autism',
    ]

    hot_terms_2 = [
        'calcium',
        'ca2+',
        'ca(2+)',
        'sodium',
        'nacl',
        'potassium',
        'k+',
        'k(+)',
        'channels',
        'chemotactic',
        'receptor',
        'signaling',
        'signal transduction',
        'glutamate',
        'behavior'
    ]

    penalty = [
        'associated',
        'association',
        'associate'
    ]

    # Score gene based on NCBI gene summary
    for gene_name in ref.index.tolist():

        gene_id = ref.at[gene_name, 'NCBI gene ID']
        url = 'https://www.ncbi.nlm.nih.gov/gene/' + str(gene_id)

        try:
            req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            webpage = urlopen(req).read()
            page = webpage.decode('utf-8').splitlines()

            score_1 = None
            score_2 = None
            score_3 = None
            for line in range(len(page) + 1):
                if '</html>' in page[line]:
                    break

                if '<dt>Summary</dt>' in page[line]:
                    score_1 = 0
                    score_2 = 0
                    score_3 = 0
                    for x in hot_terms_1:
                        if x.lower() in page[line+1].lower():
                            score_1 += 2

                    for y in hot_terms_2:
                        if y.lower() in page[line+1].lower():
                            score_2 += 1

                    for z in penalty:
                        if z.lower() in page[line+1].lower():
                            score_3 -= 0.5

                    # If only hit the penalty, revert back to zero
                    if score < 0:
                        score = 0
                    ref.at[gene_name, 'NCBI_SCORE_1'] = score_1
                    ref.at[gene_name, 'NCBI_SCORE_2'] = score_2
                    ref.at[gene_name, 'NCBI_SCORE_PENALTY'] = score_3
                    break

                if 'var tissues_data =' in page[line]:
                    exp_data = ast.literal_eval(page[line][27:-1])
                    exp_df = pd.DataFrame.from_dict(exp_data,orient='index')
                    exp_len = len(exp_df.index.tolist())
                    exp_loc = pd.DataFrame.from_dict(exp_data,orient='index').sort_values('full_rpkm', ascending=False).index.get_loc('brain')
                    exp_rank = (exp_len - exp_loc) / exp_len
                    ref.at[gene_name, 'NCBI_EXPRESSION_RANK'] = exp_rank

        except:
            ref.at[gene_name, 'NCBI_SCORE_1'] = None
            ref.at[gene_name, 'NCBI_SCORE_2'] = None
            ref.at[gene_name, 'NCBI_SCORE_PENALTY'] = None
            ref.at[gene_name, 'NCBI_EXPRESSION_RANK'] = None

    # Score gene based on UNIPROT gene summary
    for gene_name in ref.index.tolist():

        uniprot_id = ref.at[gene_name, 'UniProt accession']
        url = 'https://www.uniprot.org/uniprot/' + str(uniprot_id)

        try:
            req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            webpage = urlopen(req).read()
            page = webpage.decode('utf-8').splitlines()

            score_1 = None
            score_2 = None
            score_3 = None
            for line in range(len(page) + 1):
                if '</html>' in page[line]:
                    break

                if 'This section provides any useful information about the protein, mostly biological knowledge' in page[line] or 'This subsection of the ‘Pathology and Biotech’ section provides information' in page[line]:
                    score_1 = 0
                    score_2 = 0
                    score_3 = 0
                    for x in hot_terms_1:
                        if x.lower() in page[line].lower():
                            score_1 += 2

                    for y in hot_terms_2:
                        if y.lower() in page[line].lower():
                            score_2 += 1

                    for z in penalty:
                        if z.lower() in page[line].lower():
                            score_3 -= 0.5

                    # If only hit the penalty, revert back to zero
                    if score_3 < 0:
                        score_3 = 0
                    ref.at[gene_name, 'UNIPROT_SCORE1'] = score_1
                    ref.at[gene_name, 'UNIPROT_SCORE2'] = score_2
                    ref.at[gene_name, 'UNIPROT_SCORE_PENALTY'] = score_3
                    break

                exp_score = 0
                if 'This subsection of the ‘Expression’ section provides information on the expression of the gene product' in page[line]:
                    for x in hot_terms_1:
                        if x.lower() in page[line].lower():
                            exp_score += 1

                ref.at[gene_name, 'UNIPROT_EXPRESSION'] = exp_score

        except:
            ref.at[gene_name, 'UNIPROT_SCORE1'] = None
            ref.at[gene_name, 'UNIPROT_SCORE2'] = None
            ref.at[gene_name, 'UNIPROT_SCORE_PENALTY'] = None
            ref.at[gene_name, 'UNIPROT_EXPRESSION'] = None

    return ref.loc[(ref['NCBI_SCORE_1'] > 0) | (ref['UNIPROT_SCORE1'] > 0) | (ref['NCBI_EXPRESSION_RANK'] > 0.75)].index.tolist(), len(ref.loc[(ref['NCBI_SCORE_1'] > 0) | (ref['UNIPROT_SCORE1'] > 0) | (ref['NCBI_EXPRESSION_RANK'] > 0.75)].index.tolist())



"""
READ IN DATA
"""
# Read in single-pass XPRESSpipe processed read data quantified with HTSeq using non-de-deuplicated alignments
data = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/isrib_riboseq_xpresspipe_count_table.tsv',
    sep = '\t',
    index_col = 0)

data.head()

# Accessed via:
# $ curl -O ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz
# $ gzip -d Homo_sapiens.GRCh38.96.gtf.gz
data = convert_names(
    data,
    '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.96.gtf')
data.columns = data.columns.str.replace('_1_Aligned', '')
data.shape

# Combine lanes
sra_info = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/GSE65778_table.txt',
    sep = '\t')
data = name_map(data, sra_info)
data.head()

# Clean up data
data = data.groupby(level=0).sum() # Combine duplicate named genes
data_threshold = data[data[['untr_a_hek', 'untr_b_hek', 'tm_a_hek', 'tm_b_hek', 'tmisrib_a_hek', 'tmisrib_b_hek', 'isrib_a_hek', 'isrib_b_hek']].min(axis=1) > 25] # Apply threshold to data
data_rpm = rpm(data_threshold)

"""
READ IN TOPHAT SAMPLES
"""
data_tophat = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/tophat_isrib_counts_table.tsv',
    sep = '\t',
    index_col = 0)

data_tophat.head()

data_tophat = convert_names(
    data_tophat,
    '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.96.gtf')
data_tophat.columns = data_tophat.columns.str.replace('_1_Aligned', '')
data_tophat.shape

# Combine lanes
sra_info = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/GSE65778_table.txt',
    sep = '\t')
data_tophat = name_map(data_tophat, sra_info)
data_tophat.head()


"""
READ IN INGOLIA DATA
"""
# Read in Ingolia raw counts from Elife supplement table
ingolia = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/ingolia_counts_table.txt',
    sep = '\t',
    index_col = 0)
ingolia = ingolia.drop('size', axis = 1)

# Clean up data
ingolia = ingolia.groupby(level=0).sum() # Combine duplicate named genes
ingolia_threshold = ingolia[ingolia[['untr_a_hek', 'untr_b_hek', 'tm_a_hek', 'tm_b_hek', 'tmisrib_a_hek', 'tmisrib_b_hek', 'isrib_a_hek', 'isrib_b_hek']].min(axis=1) > 25]
ingolia_threshold = ingolia[ingolia.min(axis=1) > 25] # Apply threshold to data
ingolia_rpm = rpm(ingolia_threshold)


"""
GET COMMON GENE SET BETWEEN DATASETS
"""
# Clean NaNs
data_rpm.shape
ingolia_rpm.shape

data_rpm = data_rpm.dropna()
ingolia_rpm = ingolia_rpm.dropna()

data_rpm.shape
ingolia_rpm.shape


# Make dataframes for each dataset for matching genes
data_genes = data.index.tolist()
ingolia_genes = ingolia.index.tolist()
len(data_genes)
len(ingolia_genes)

common_genes = list(set(data_genes).intersection(ingolia_genes))
len(common_genes)

data_common = data.reindex(index = common_genes)
ingolia_common = ingolia.reindex(index = common_genes)


"""
FIGURE 2A
"""
# Run correlations between sample alignments
make_figure2A(
    data_common,
    ingolia_common,
    'external_correlations_summary_htseq.png')

"""
FIGURE 2B
"""
# Run correlations between sample alignments
make_figure2B(
    data,
    'internal_correlations_summary_htseq.png')

"""
FIGURE 3A
"""
# Tm treatment vs Untreated
untr_mrna = ['untr_a_hek', 'untr_b_hek']
untr_ribo = ['ribo_untr_a', 'ribo_untr_b']
tm_mrna = ['tm_a_hek', 'tm_b_hek']
tm_ribo = ['ribo_tm_a', 'ribo_tm_b']
tmisrib_mrna = ['tmisrib_a_hek', 'tmisrib_b_hek']
tmisrib_ribo = ['ribo_tmisrib_a', 'ribo_tmisrib_b']
isrib_mrna = ['isrib_a_hek', 'isrib_b_hek']
isrib_ribo = ['ribo_isrib_a', 'ribo_isrib_b']

"""Tm vs Untr"""
isr = ['ATF5','ATF4','PPP1R15A','DDIT3']
down_intersection = []
tm_data = rp_plot(data_rpm, tm_mrna, tm_ribo, untr_mrna, untr_ribo, 'Tm', 'Untr', 'Tm_vs_Untr_jordan.png', isr, down_intersection)

tm_data.loc['ATF4']['ribo_fc']
tm_data.loc['ATF5']['ribo_fc']
tm_data.loc['PPP1R15A']['ribo_fc']
tm_data.loc['DDIT3']['ribo_fc']

"""Tm + ISRIB vs Untr"""
tmisrib_data = rp_plot(data_rpm, tmisrib_mrna, tmisrib_ribo, untr_mrna, untr_ribo, 'Tm+ISRIB', 'Untr', 'TmISRIB_vs_Untr_jordan.png', isr, down_intersection)

"""ISRIB vs Untr"""
isrib_data = rp_plot(data_rpm, isrib_mrna, isrib_ribo, untr_mrna, untr_ribo, 'ISRIB', 'Untr', 'ISRIB_vs_Untr_jordan.png', isr, down_intersection)

"""
FIGURE 3B
"""



"""
TE Analysis
"""
order = [
    'ribo_untr_a',
    'untr_a_hek',
    'ribo_untr_b',
    'untr_b_hek',
    'ribo_tm_a',
    'tm_a_hek',
    'ribo_tm_b',
    'tm_b_hek',
    'ribo_tmisrib_a',
    'tmisrib_a_hek',
    'ribo_tmisrib_b',
    'tmisrib_b_hek',
    'ribo_isrib_a',
    'isrib_a_hek',
    'ribo_isrib_b',
    'isrib_b_hek']

data_rpm_c = data_rpm.copy()
data_ordered = data_rpm_c[order]

data_te = pd.DataFrame()
data_te['untr_a'] = data_ordered['ribo_untr_a'] / data_ordered['untr_a_hek']
data_te['untr_b'] = data_ordered['ribo_untr_b'] / data_ordered['untr_b_hek']
data_te['tm_a'] = data_ordered['ribo_tm_a'] / data_ordered['tm_a_hek']
data_te['tm_b'] = data_ordered['ribo_tm_b'] / data_ordered['tm_b_hek']
data_te['tmisrib_a'] = data_ordered['ribo_tmisrib_a'] / data_ordered['tmisrib_a_hek']
data_te['tmisrib_b'] = data_ordered['ribo_tmisrib_b'] / data_ordered['tmisrib_b_hek']
data_te['isrib_a'] = data_ordered['ribo_isrib_a'] / data_ordered['isrib_a_hek']
data_te['isrib_b'] = data_ordered['ribo_isrib_b'] / data_ordered['isrib_b_hek']
data_te += 0.01
data_te = np.log2(data_te)

check_samples(data_te)

data_te_normed = pd.DataFrame()
data_te_normed['untr'] = (data_te['untr_a'] + data_te['untr_a']) / 2
data_te_normed['tm'] = (data_te['tm_a'] + data_te['tm_b']) / 2
data_te_normed['tmisrib'] = (data_te['tmisrib_a'] + data_te['tmisrib_b']) / 2
data_te_normed['isrib'] = (data_te['isrib_a'] + data_te['isrib_b']) / 2

data_te_zero = data_te_normed.sub(data_te_normed['untr'], axis=0)

check_samples(data_te_zero)

down_list = data_te_zero.loc[(data_te_zero['tm'] <= -1) & (data_te_zero['tmisrib'] >= 0)].index.tolist()
down_list_neuro = [
 'ADAP1',
 'ADRA2C',
 'GPR183',
 'HTRA1',
 'KIFC2',
 'RGPD8',
 'SLC1A1']

data_te_zero_plot = data_te_zero.copy()
data_te_zero_plot.columns = ['Untreated','Tm','Tm + ISRIB','ISRIB']

fig, ax = plt.subplots()
plt.yticks([-4,-3,-2,-1,0,1,2,3,4,5,6])
ax.axvline(0.02, ls='-', color='black')
ax.axvline(1, ls='-', color='black')
ax.axvline(2, ls='-', color='black')
ax.axvline(2.99, ls='-', color='black')

ax.set_facecolor('white')
ax.grid(color='grey', axis='y')
data_te_zero_plot.T.plot.line(legend=False, color='lightgrey', ax=ax)
data_te_zero_plot.loc[isr].T.plot.line(legend=False, color='#be00be', ax=ax)
data_te_zero_plot.loc[down_list].T.plot.line(legend=False, color='#98fb98', ax=ax)
data_te_zero_plot.loc[down_list_neuro].T.plot.line(legend=False, color='#004600', ax=ax)
ax.set_ylabel('log$_2$(TE)')

ax.axhline(0, ls='-', color='black')
ax.axhline(1, ls='--', color='black')
ax.axhline(-1, ls='--', color='black')

from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], color='gray', lw=2, label='All'),
                   Line2D([0], [0], color='#be00be', lw=2, label='ISR'),
                   Line2D([0], [0], color='green', lw=2, label='Neuro')]

ax.legend(handles=legend_elements, loc='upper right')
plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/plots/te_analysis.png', dpi=1800, bbox_inches='tight')
plt.close()

"""
down_list ->
['AC097634.4',
 'ADAP1',
 'ADRA2C',
 'AL109811.3',
 'AL358472.7',
 'CMYA5',
 'CXorf40B',
 'DCDC2B',
 'ECHDC2',
 'EFCAB13',
 'GNB3',
 'GPR183',
 'HTRA1',
 'KIFC2',
 'LRP5L',
 'MYO5B', --> shows up in Ingolia hits
 'RGPD8',
 'RPS27',
 'SCART1',
 'SLC1A1',
 'SLC9A3',
 'TLL1',
 'TMEM110-MUSTN1',
 'WNK4']

down_list_neuro ->
 ['ADAP1',
 'ADRA2C',
 'GPR183',
 'HTRA1',
 'KIFC2',
 'RGPD8',
 'SLC1A1']
"""

ref = make_ref('/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/human_biomart.txt', down_list)
ref, ref_len = run_scraper(ref)
print(ref_len)
ref

del val_te

# Get random subset same size as down_list and compare annotations
val_random = []
for x in range(25):

    gene_list = random.sample(data_te_zero.index.tolist(), k=len(down_list))
    ref = make_ref('/Users/jordan/Desktop/xpressyourself_manuscript/7283587_analysis/human_biomart.txt', gene_list)
    val_te = run_scraper(ref)
    val_random.append(val_te)


# From 25 iterations performed 19 Jun 2019
[
(['CAMK2D', 'THAP9'], 2),
(['ZDHHC21', 'NCAPH2', 'PLCXD2', 'RNF219'], 4),
(['ZNF786', 'PITX1', 'C5orf51'], 3),
(['AMN1', 'SCNN1D', 'PWWP2A', 'CEP97', 'UBE2O', 'HEY1'], 6),
(['BCAS3', 'DNM2', 'LRRC47', 'SERPINF1', 'DAXX', 'AP4S1'], 6),
(['MGARP', 'MB21D2', 'UNC119', 'C17orf75', 'RADIL', 'ZNF461', 'GALNT13'], 7),
(['GAPVD1', 'CLCC1', 'HEATR5B', 'ZNF711', 'TEX2', 'FAXC', 'KDM1A', 'HS6ST1'], 8),
(['FAM229B', 'ANKIB1', 'BTBD10', 'PLAT', 'PRX'], 5),
(['ZDHHC17', 'FAM102A', 'JAKMIP1', 'FOXK1', 'MMD', 'PPIG', 'STK24'], 7),
(['ICA1L', 'ZNF594', 'ADAL', 'INA', 'CEBPZOS', 'RADX', 'AKT1S1'], 7),
(['LRRC27', 'MTMR2', 'MCUR1', 'SEMA3A', 'EPHA7'], 5),
(['CRIM1', 'MZT1', 'UBALD1'], 3),
(['UNC13B', 'RAD54L2', 'LIN7A', 'B4GALT6', 'GALNT13', 'ZBTB6', 'RABEPK', 'SLC17A5', 'AMZ2', 'CEBPZOS', 'FAM169A', 'GSK3B'], 12),
(['COPG2', 'ZNF461'], 2),
(['PCYOX1L', 'RNF115', 'EIF1', 'SYK', 'SLC2A13', 'CDK17'], 6),
(['HEATR5B', 'WASHC3', 'STRN4', 'GNAQ', 'EPB41L3'], 5),
(['RNASEK', 'PLCXD2', 'RRN3', 'PRX', 'AMZ2'], 5),
(['DLL1', 'FAM199X', 'TRAPPC6B', 'C9orf40', 'QKI', 'ADGRL3'], 6),
(['PSAP', 'UBE2O', 'SMARCD1'], 3),
(['PI4KB', 'HEY1', 'NMU', 'OCIAD1', 'IWS1'], 5),
(['HMGCS1', 'MICOS13', 'GLUD1', 'SYNGR1', 'HDX', 'LRFN4'], 6),
(['ANKRD6', 'ZNF614', 'SH3RF3', 'RBM48', 'FAM172A'], 5),
(['KCNQ2', 'NOLC1', 'ZNF614', 'MDK', 'HCFC1R1', 'PAIP2'], 6),
(['ZNF341', 'QSOX2', 'GPR155', 'NDNF'], 4),
(['KHDRBS1', 'SFR1', 'CC2D1A'], 3)]

sum([2,4,3,6,6,7,8,5,7,7,5,3,12,2,6,5,5,6,3,5,6,5,6,4,3]) / 25
