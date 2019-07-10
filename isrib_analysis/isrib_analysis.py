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
import plotly
import plotly.offline as py
import plotly_express as px

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
    title,
    interactive=False):

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

        if interactive == True:
            merged_best = pd.concat([data[[str(x)]], ingolia_data[[str(x)]]], axis=1, sort=False)
            merged_best.columns = ['xpresspipe', 'ingolia']
            merged_best['genes'] = merged_best.index
            sc = px.scatter(
                merged_best,
                x='xpresspipe',
                y='ingolia',
                hover_name='genes',
                log_x=True,
                log_y=True,
                opacity=0.4,
                width=1400,
                height=1000,
                title=str(x))

            py.offline.plot(sc, filename='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(title)[:-4] + '.html')

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
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

def make_figure2S(
    data,
    ingolia_data,
    title):

    fig, axes = plt.subplots(
        nrows = 2,
        ncols = 2,
        figsize = (10, 10),
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
        'untr_a_hek'] # Designate sample order

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
        ax_y = 0
        if file_number == 0:
            ax_x = 0
        elif file_number == 1:
            ax_x = 1
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

    cols = ['Untr RPF RepA','Untr mRNA RepA']
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=24)
        ax.xaxis.set_label_position('top')

    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(counts)', fontsize=16)

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(title),
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
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

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

def interactive_scatter(sample, data, ingolia_data, data1='xpresspipe', data2='ingolia'):

    merged_best = pd.concat([data[[sample]], ingolia_data[[sample]]], axis=1, sort=False)
    merged_best.columns = [data1,data2]
    merged_best['genes'] = merged_best.index.tolist()

    merged_best['color'] = 'others'
    merged_best.loc[merged_best['genes'] == 'EIF3C', 'color'] = 'EIF3C'
    merged_best.loc[merged_best['genes'] == 'TUBA1A', 'color'] = 'TUBA1A'
    merged_best.loc[merged_best['genes'] == 'NOMO2', 'color'] = 'NOMO2'
    merged_best.loc[merged_best['genes'] == 'RPL39', 'color'] = 'RPL39'
    merged_best.loc[merged_best['genes'] == 'KCNQ2', 'color'] = 'KCNQ2'

    sc = px.scatter(
        merged_best,
        x=merged_best.columns.tolist()[0],
        y=merged_best.columns.tolist()[1],
        hover_name='genes',
        color='color',
        color_discrete_sequence = ['lightgrey','red','blue','green','black','orange'],
        category_orders = {'others':0,'EIF3C':1,'TUBA1A':2,'NOMO2':3,'RPL39':4,'KCNQ2':5},
        log_x=True,
        log_y=True,
        opacity=1,
        width=1400,
        height=1000,
        title=sample)

    py.offline.plot(sc, filename='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(sample) + '.html')


# Read in single-pass XPRESSpipe processed read data quantified with HTSeq using non-de-deuplicated alignments
def get_data(file, sample_suffix='_1_Aligned'):
    data = pd.read_csv(
        file,
        sep = '\t',
        index_col = 0)
    data.index = data.index.str.split('.').str[0]

    # Accessed via:
    # $ curl -O ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz
    # $ gzip -d Homo_sapiens.GRCh38.96.gtf.gz
    data = convert_names(
        data,
        '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.96.gtf')
    data.columns = data.columns.str.replace(str(sample_suffix), '')
    data.shape

    # Combine lanes
    sra_info = pd.read_csv(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/GSE65778_table.txt',
        sep = '\t')
    data = name_map(data, sra_info)
    data = data.groupby(level=0).sum() # Combine duplicate named genes

    return data

file = '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_riboseq_xpresspipe_count_table.tsv'
data = get_data(file)


"""
READ IN DATA
"""
# Clean up data
data_threshold = data.loc[data[['untr_a_hek', 'untr_b_hek', 'tm_a_hek', 'tm_b_hek', 'tmisrib_a_hek', 'tmisrib_b_hek', 'isrib_a_hek', 'isrib_b_hek']].min(axis=1) >= 10] # Apply threshold to data

data_rpm = rpm(data_threshold)

# Export for DESeq2
data_threshold.to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/isribxpresspipe_thresholded_counts.tsv',sep='\t')

untr_mrna = ['untr_a_hek', 'untr_b_hek']
untr_ribo = ['ribo_untr_a', 'ribo_untr_b']
tm_mrna = ['tm_a_hek', 'tm_b_hek']
tm_ribo = ['ribo_tm_a', 'ribo_tm_b']
tmisrib_mrna = ['tmisrib_a_hek', 'tmisrib_b_hek']
tmisrib_ribo = ['ribo_tmisrib_a', 'ribo_tmisrib_b']
isrib_mrna = ['isrib_a_hek', 'isrib_b_hek']
isrib_ribo = ['ribo_isrib_a', 'ribo_isrib_b']

data_threshold[untr_ribo + tm_ribo].to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/tm_ribo_ISRIBxpresspipe_processed_counts.tsv',sep='\t')
data_threshold[untr_ribo + tmisrib_ribo].to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/tmisrib_ribo_ISRIBxpresspipe_processed_counts.tsv',sep='\t')
data_threshold[untr_ribo + isrib_ribo].to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/isrib_ribo_ISRIBxpresspipe_processed_counts.tsv',sep='\t')

data_threshold[untr_mrna + tm_mrna].to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/tm_rna_ISRIBxpresspipe_processed_counts.tsv',sep='\t')
data_threshold[untr_mrna + tmisrib_mrna].to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/tmisrib_rna_ISRIBxpresspipe_processed_counts.tsv',sep='\t')
data_threshold[untr_mrna + isrib_mrna].to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/isrib_rna_ISRIBxpresspipe_processed_counts.tsv',sep='\t')

data_threshold[untr_ribo + tm_ribo + untr_mrna + tm_mrna].to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/tm_ISRIBxpresspipe_processed_counts.tsv',sep='\t')
data_threshold[untr_ribo + tmisrib_ribo + untr_mrna + tmisrib_mrna].to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/tmisrib_ISRIBxpresspipe_processed_counts.tsv',sep='\t')
data_threshold[untr_ribo + isrib_ribo + untr_mrna + isrib_mrna].to_csv('/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/isrib_ISRIBxpresspipe_processed_counts.tsv',sep='\t')

"""
READ IN TOPHAT SAMPLES
"""
tophat_file = '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/tophat_isrib_counts_count_table.tsv'
data_tophat = pd.read_csv(
    tophat_file,
    sep = '\t',
    index_col = 0)

data_tophat = convert_names(
    data_tophat,
    '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.96.gtf')

data_tophat.columns = data_tophat.columns.str.replace('_accepted_hits', '')
data_tophat.shape

# Combine lanes
sra_info = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/GSE65778_table.txt',
    sep = '\t')
data_tophat = name_map(data_tophat, sra_info)
data_tophat = data_tophat.groupby(level=0).sum()
data_tophat_threshold = data_tophat[data_tophat[['untr_a_hek']].min(axis=1) > 25] # Apply threshold to data


"""
READ IN INGOLIA DATA
"""
# Read in Ingolia raw counts from Elife supplement table
ingolia = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/ingolia_counts_table.txt',
    sep = '\t',
    index_col = 0)
ingolia = ingolia.drop('size', axis = 1)

# Clean up data
ingolia = ingolia.groupby(level=0).sum() # Combine duplicate named genes
ingolia_threshold = ingolia[ingolia[['untr_a_hek', 'untr_b_hek', 'tm_a_hek', 'tm_b_hek', 'tmisrib_a_hek', 'tmisrib_b_hek', 'isrib_a_hek', 'isrib_b_hek']].min(axis=1) > 50]
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


# Tophat
data_tophat_genes = data_tophat.index.tolist()
ingolia_genes_thcomp = ingolia.index.tolist()
len(data_tophat_genes)
len(ingolia_genes_thcomp)

common_genes_th = list(set(data_tophat_genes).intersection(ingolia_genes_thcomp))
len(common_genes_th)

data_tophat_common = data_tophat.reindex(index = common_genes_th)
ingolia_common_thcomp = ingolia.reindex(index = common_genes_th)



"""
FIGURE 2A
"""
# Run correlations between sample alignments
def cycle_fig2a(file, ingolia):

    data = get_data(file, sample_suffix='__Aligned')

    data_genes = data.index.tolist()
    ingolia_genes = ingolia.index.tolist()

    common_genes = list(set(data_genes).intersection(ingolia_genes))

    data_common = data.reindex(index = common_genes)
    ingolia_common = ingolia.reindex(index = common_genes)

    make_figure2A(
        data_common,
        ingolia_common,
        str(file.split('/')[-1][:-4]) + '_external_correlations_summary_htseq.png')


file_list = [
    'isrib_comp_v76_truncated_count_table.tsv',
    'isrib_comp_v76_normal_count_table.tsv',
    'isrib_comp_v76_longest_truncated_count_table.tsv',
    'isrib_comp_v96_truncated_count_table.tsv',
    'isrib_comp_v96_normal_count_table.tsv',
    'isrib_comp_v96_longest_truncated_count_table.tsv',
    'isrib_comp_v72_truncated_count_table.tsv',
    'isrib_comp_v72_normal_count_table.tsv',
    'isrib_comp_v72_longest_truncated_count_table.tsv'
]

file_list = ['/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_comp_test/' + str(f) for f in file_list]

for file in file_list:

    cycle_fig2a(file, ingolia)


make_figure2A(
    data_common,
    ingolia_common,
    'external_correlations_summary_htseq.png')

make_figure2S(
    data_tophat_common,
    ingolia_common_thcomp,
    'external_correlations_summary_tophat.png')


interactive_scatter('ribo_untr_a', data_common, ingolia_common)
interactive_scatter('ribo_tm_a', data_common, ingolia_common)
interactive_scatter('untr_a_hek', data_common, ingolia_common)
interactive_scatter('tm_a_hek', data_common, ingolia_common)



"""
FIGURE 2B
"""
# Run correlations between sample alignments
make_figure2B(
    data,
    'internal_correlations_summary_htseq.png')
