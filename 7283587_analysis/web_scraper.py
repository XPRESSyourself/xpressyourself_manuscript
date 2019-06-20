# Gene dictionary accessed from https://biomart.genenames.org
import os
import sys
import ast
import pandas as pd
import numpy as np
from urllib.request import Request, urlopen

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/Desktop/xpressyourself_manuscript/7187289_analysis/'
+ 'human_biomart.txt'

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
