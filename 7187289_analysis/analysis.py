import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
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
    plt.yticks([1,2,3,4,5]) # Limit axis labels to ints
    plt.xticks([0,1,2,3,4])

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
        axes[ax_y, ax_x].axhline(1, ls='-', color='black') # Create axis lines
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymax=0.88)
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
        '/Users/jordan/Desktop/xpressyourself_manuscript/7187289_analysis/plots/' + str(title),
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
        axes[ax_y, ax_x].axhline(0, ls='-', color='black') # Create axis lines
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymax=0.88)
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
        '/Users/jordan/Desktop/xpressyourself_manuscript/7187289_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

# Make distribution
def make_figure2C(
    data,
    title):

    ref_line = np.log10(data + 1).mean().mean()
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
    ax = sns.violinplot(data=np.log10(data + 1), order = names)
    plt.xticks(rotation=45)
    ax.set_ylabel('Normalized Counts (RPM)')
    ax.set_xlabel('Samples')
    ax.axhline(ref_line, ls='--', color='grey')
    plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/7187289_analysis/plots/' + str(title), dpi=1800, bbox_inches='tight')
    plt.close

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

    plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/7187289_analysis/plots/' + str(title), dpi=1800, bbox_inches='tight')
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
        ''

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




    return len(ref.loc[(ref['NCBI_SCORE_1'] > 0) | (ref['UNIPROT_SCORE1'] > 0) | (ref['NCBI_EXPRESSION_RANK'] > 0.75)].index.tolist())



"""
READ IN DATA
"""
# Read in single-pass XPRESSpipe processed read data quantified with HTSeq using non-de-deuplicated alignments
data = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/7187289_analysis/isrib_riboprof_count_table.tsv',
    sep = '\t',
    index_col = 0)

data.head()

data = convert_names(
    data,
    '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.95.gtf')
data.columns = data.columns.str.replace('_1_Aligned', '')
data.shape

# Combine lanes
sra_info = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/GSE65778_table.txt',
    sep = '\t')
data = name_map(data, sra_info)
data.head()

# Clean up data
data = data.groupby(level=0).sum() # Combine duplicate named genes
data_threshold = data[data[['untr_a_hek', 'untr_b_hek', 'tm_a_hek', 'tm_b_hek', 'tmisrib_a_hek', 'tmisrib_b_hek', 'isrib_a_hek', 'isrib_b_hek']].min(axis=1) > 25] # Apply threshold to data
data_rpm = rpm(data_threshold)


"""
READ IN DE-DUPLICATED DATA
"""
# Read in single-pass XPRESSpipe processed read data quantified with HTSeq using de-deuplicated alignments
"""data_deduped = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/isrib_riboprof_count_table_7155666_deduplicated_htseq.tsv',
    sep = '\t',
    index_col = 0)
data_deduped = convert_names(
    data_deduped,
    '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.95.gtf')
data_deduped.columns = data_deduped.columns.str.replace('_1_dedupRemoved', '')
data_deduped.shape

# Combine lanes
sra_info = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/GSE65778_table.txt',
    sep = '\t')
data_deduped = name_map(data_deduped, sra_info)
data_deduped.head()

# Clean up data
data_deduped = data_deduped.groupby(level=0).sum() # Combine duplicate named genes
data_deduped_threshold = data_deduped[data_deduped.min(axis=1) > 25] # Apply threshold to data
data_deduped_rpm = rpm(data_deduped_threshold)"""


"""
READ IN INGOLIA DATA
"""
# Read in Ingolia raw counts from Elife supplement table
ingolia = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/ingolia_counts_table.txt',
    sep = '\t',
    index_col = 0)
ingolia = ingolia.drop('size', axis = 1)

# Clean up data
ingolia = ingolia.groupby(level=0).sum() # Combine duplicate named genes
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
data_genes = data_threshold.index.tolist()
# data_deduped_genes = data_deduped_threshold.index.tolist()
ingolia_genes = ingolia_threshold.index.tolist()
len(data_genes)
# len(data_deduped_genes)
len(ingolia_genes)

common_genes = list(set(data_genes).intersection(ingolia_genes))
len(common_genes)
# deduped_common_genes = list(set(data_deduped_genes).intersection(ingolia_genes))
# len(deduped_common_genes)

data_common = data_threshold.reindex(index = common_genes)
ingolia_common = ingolia_threshold.reindex(index = common_genes)

# data_deduped_common = data_deduped_threshold.reindex(index = deduped_common_genes)
# ingolia_deduped_common = ingolia_threshold.reindex(index = deduped_common_genes)

"""
FIGURE 2A
"""
# Run correlations between sample alignments
make_figure2A(
    data_common,
    ingolia_common,
    'external_correlations_summary_htseq_non-deduplicated.png')

#make_figure2A(
#    data_deduped_common,
#    ingolia_deduped_common,
#    'external_correlations_summary_htseq_deduplicated.png')

"""
FIGURE 2B
"""
# Run correlations between sample alignments
make_figure2B(
    data,
    'internal_correlations_summary.png')

# make_figure2B(
#    data_deduped,
#    'internal_correlations_summary_deduped.png')

"""
FIGURE 2C
My alignments appear to have less outliers, indicating a better read processing overall
"""
make_figure2C(
    data_rpm,
    'xpresspipe_run_distributions.png')

# make_figure2C(
#    data_deduped_rpm,
#    'xpresspipe_run_distributions_deduped.png')

make_figure2C(
    ingolia_rpm,
    'ingolia_run_distributions.png')

"""
FIGURE 3A
My alignments have less "Super-translators", which is actually better and expected
Not sure why the mRNA levels seem to be more constained on the positive side
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
down_intersection = ['MYO5B', 'NOTCH3', 'NGEF', 'POMGNT1', 'ENPP5', 'CCNG2', 'SMURF2', 'ARHGEF6', 'FGFRL1']

up_ingolia = ['ATF4','ATF5']
down_ingolia = ['MYO5B']

tm_data = rp_plot(data_rpm, tm_mrna, tm_ribo, untr_mrna, untr_ribo, 'Tm', 'Untr', 'Tm_vs_Untr_jordan.png', isr, down_intersection)

# tm_data_deduped = rp_plot(data_deduped_rpm, tm_mrna, tm_ribo, untr_mrna, untr_ribo, 'Tm', 'Untr', 'Tm_vs_Untr_jordan_deduped.png', isr, down_intersection)
tm_data.loc['ATF4']['ribo_fc']
tm_data.loc['ATF5']['ribo_fc']
tm_data_down = tm_data.loc[(tm_data['ribo_fc'] < 0.6) & (tm_data['p_adj'] < 0.05)].index.tolist()

# tm_data_deduped_down = tm_data_deduped.loc[(tm_data_deduped['ribo_fc'] < 0.5) & (tm_data_deduped['p_adj'] < 0.05)].index.tolist()
# dup_dedup_down_list = list(set(tm_data_down) & set(tm_data_deduped_down))
# tm_data.loc[dup_dedup_down_list]
# tm_data_deduped.loc[dup_dedup_down_list]


"""Tm + ISRIB vs Untr"""
tmisrib_data = rp_plot(data_rpm, tmisrib_mrna, tmisrib_ribo, untr_mrna, untr_ribo, 'Tm+ISRIB', 'Untr', 'TmISRIB_vs_Untr_jordan.png', isr, down_intersection)
# tmisrib_data_deduped = rp_plot(data_deduped_rpm, tmisrib_mrna, tmisrib_ribo, untr_mrna, untr_ribo, 'Tm+ISRIB', 'Untr', 'TmISRIB_vs_Untr_jordan_deduped.png', isr, down_intersection)

"""ISRIB vs Untr"""
isrib_data = rp_plot(data_rpm, isrib_mrna, isrib_ribo, untr_mrna, untr_ribo, 'ISRIB', 'Untr', 'ISRIB_vs_Untr_jordan.png', isr, down_intersection)
# isrib_data_deduped = rp_plot(data_deduped_rpm, isrib_mrna, isrib_ribo, untr_mrna, untr_ribo, 'ISRIB', 'Untr', 'ISRIB_vs_Untr_jordan_deduped.png', isr, down_intersection)


# Ingolia plots
ingolia_rpm.shape
ingolia_rpm = ingolia_rpm.dropna()
ingolia_rpm.shape

ingolia_tm = rp_plot(ingolia_rpm, tm_mrna, tm_ribo, untr_mrna, untr_ribo, 'Tm', 'Untr', 'Tm_vs_Untr_ingolia.png', up_ingolia, down_ingolia)
ingolia_tm.loc['ATF4']['ribo_fc']
ingolia_tm.loc['ATF5']['ribo_fc']

rp_plot(ingolia_rpm, tmisrib_mrna, tmisrib_ribo, untr_mrna, untr_ribo, 'Tm+ISRIB', 'Untr', 'TmISRIB_vs_Untr_ingolia.png', isr, down_intersection)

rp_plot(ingolia_rpm, isrib_mrna, isrib_ribo, untr_mrna, untr_ribo, 'ISRIB', 'Untr', 'ISRIB_vs_Untr_ingolia.png', isr, down_intersection)


"""
FIGURE 3B
"""
# Track neuronal related genes across conditions to see if they are significantly changed between conditions
tm_data_set = tm_data.reindex(labels=down_intersection)
tmisrib_data_set = tmisrib_data.reindex(labels=down_intersection)
isrib_data_set = isrib_data.reindex(labels=down_intersection)

tm_data_set['untr_te'] =  (tm_data_set['conditionA_ribo']) / (tm_data_set['conditionA_mrna'])
tm_data_set['tm_te'] = (tm_data_set['conditionB_ribo']) / (tm_data_set['conditionB_mrna'])
tm_data_set['tmisrib_te'] = (tmisrib_data_set['conditionB_ribo']) / (tmisrib_data_set['conditionB_mrna'])
tm_data_set['isrib_te'] = (isrib_data_set['conditionB_ribo']) / (isrib_data_set['conditionB_mrna'])
tm_data_set['Tm'] = np.log2(tm_data_set['tm_te'] / tm_data_set['untr_te'])
tm_data_set['Tm + ISRIB'] = np.log2(tm_data_set['tmisrib_te'] / tm_data_set['untr_te'])
tm_data_set['ISRIB'] = np.log2(tm_data_set['isrib_te'] / tm_data_set['untr_te'])

dataset = pd.concat([tm_data_set[['Tm']], tm_data_set[['Tm + ISRIB']], tm_data_set[['ISRIB']]], axis=1, sort=False)
dataset['Tm + ISRIB'] - dataset['Tm']


ax = dataset.T.plot.line(figsize=(5,7.5), grid = False)
ax.axvline(1, ls='-', color='white')
ax.axhline(0, ls='--', color='red')
ax.set_ylabel('log$_2$$\Delta$TE (vs Untreated)')
plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/plots/elife_figures/translation_efficiencies_interest_down_intersection.png', dpi=1800, bbox_inches='tight')


tm_data_set = tm_data.reindex(labels=isr)
tmisrib_data_set = tmisrib_data.reindex(labels=isr)
isrib_data_set = isrib_data.reindex(labels=isr)

tm_data_set['untr_te'] =  (tm_data_set['conditionA_ribo']) / (tm_data_set['conditionA_mrna'])
tm_data_set['tm_te'] = (tm_data_set['conditionB_ribo']) / (tm_data_set['conditionB_mrna'])
tm_data_set['tmisrib_te'] = (tmisrib_data_set['conditionB_ribo']) / (tmisrib_data_set['conditionB_mrna'])
tm_data_set['isrib_te'] = (isrib_data_set['conditionB_ribo']) / (isrib_data_set['conditionB_mrna'])
tm_data_set['Tm'] = np.log2(tm_data_set['tm_te'] / tm_data_set['untr_te'])
tm_data_set['Tm + ISRIB'] = np.log2(tm_data_set['tmisrib_te'] / tm_data_set['untr_te'])
tm_data_set['ISRIB'] = np.log2(tm_data_set['isrib_te'] / tm_data_set['untr_te'])

dataset = pd.concat([tm_data_set[['Tm']], tm_data_set[['Tm + ISRIB']], tm_data_set[['ISRIB']]], axis=1, sort=False)
dataset['Tm + ISRIB'] - dataset['Tm']

dataset


ax = dataset.T.plot.line(figsize=(5,7.5), grid = False)
ax.axvline(1, ls='-', color='white')
ax.axhline(0, ls='--', color='red')
ax.set_ylabel('log$_2$$\Delta$TE (vs Untreated)')
plt.savefig('/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/plots/elife_figures/translation_efficiencies_interest_isr.png', dpi=1800, bbox_inches='tight')


"""
FIGURE 3C
"""
# Analyze bigwig files for genes of interest to make sure they are reliable hits



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
data_te.mean()

data_te.head()

check_samples(data_te)

data_te_normed = pd.DataFrame()
data_te_normed['untr'] = (data_te['untr_a'] + data_te['untr_a']) / 2
data_te_normed['tm'] = (data_te['tm_a'] + data_te['tm_b']) / 2
data_te_normed['tmisrib'] = (data_te['tmisrib_a'] + data_te['tmisrib_b']) / 2
data_te_normed['isrib'] = (data_te['isrib_a'] + data_te['isrib_b']) / 2

data_te_zero = data_te_normed.sub(data_te_normed['untr'], axis=0)

data_te_zero.head()

check_samples(data_te_zero)

down_list = data_te_zero.loc[(data_te_zero['tm'] <= -1) & (data_te_zero['tmisrib'] >= 0)].index.tolist()
up_list = data_te_zero.loc[(data_te_zero['tm'] >= 1) & (data_te_zero['tmisrib'] >= 0.6)].index.tolist()
neuro_diff = []
down_up_up = data_te_zero.loc[(data_te_zero['tm'] <= -0.9) & (data_te_zero['tmisrib'] >= -0.3) & (data_te_zero['isrib'] >= 0)].index.tolist()
up_up = data_te_zero.loc[(data_te_zero['tm'] <= 0.1) & (data_te_zero['tmisrib'] >= 0.6) & (data_te_zero['isrib'] >= 0.6)].index.tolist()

fig, ax = plt.subplots()
#data_te_zero.T.plot.line(legend=False, color='lightgrey', ax=ax)
data_te_zero.loc[isr].T.plot.line(legend=False, color='red', ax=ax)
data_te_zero.loc[down_list].T.plot.line(legend=False, color='green', ax=ax)












ref = make_ref('/Users/jordan/Desktop/xpressyourself_manuscript/7187289_analysis/human_biomart.txt', down_list)
val_te = run_scraper(ref)
print(val_te)
del val_te

# Get random subset same size as down_list and compare annotations
val_random = []
for x in range(10):

    gene_list = random.sample(data_te_zero.index.tolist(), k=len(down_list))
    ref = make_ref('/Users/jordan/Desktop/xpressyourself_manuscript/7187289_analysis/human_biomart.txt', gene_list)
    val_te = run_scraper(ref)
    val_random.append(val_te)

print(val_random)

sum(val_random) / len(val_random)
