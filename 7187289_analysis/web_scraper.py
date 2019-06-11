


# Gene dictionary accessed from https://biomart.genenames.org

import os
import sys
import pandas as pd

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
__path__ = '/Users/jordan/Desktop/xpressyourself_manuscript/7187289_analysis/'

# Read in gene dictionary
df = pd.read_csv(
    str(__path__) + 'human_biomart.txt',
    sep='\t',
    index_col=2)
del df.index.name
df = df[['NCBI gene ID', 'UniProt accession']]
df['NCBI gene ID'] = df['NCBI gene ID'].astype(pd.Int32Dtype())

df.head()

# Set gene list 
gene_list = [
    'AC097634.4',
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
    'MYO5B',
    'RAD51',
    'RGPD8',
    'RPS27',
    'SCART1',
    'SLC1A1',
    'SLC9A3',
    'TLL1',
    'TMEM110-MUSTN1',
    'WNK4']


gene_id = '29974'
url = "https://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&dopt=full_report&list_uids=" + str(gene_id)

from urllib.request import Request, urlopen

req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})

webpage = urlopen(req).read()







go to 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + str(gene_name)

    go to 'href="http://www.ncbi.nlm.nih.gov/' URL
    go to 'href="http://www.uniprot.org/' URL

get UNIPROT
    read in gene dictionary
    search genes and scrape html
    'http://www.uniprot.org/uniprot/' + str(uni_id) + '#function'










    parse sentance by sentence
        look for items in hit_list
            if 'association' or 'associated' in sentance
                (1)
            else:
                (2)



    if 'behavior/neurological phenotype' or 'nervous system phenotype' in file

    go to 'mRNA expression in embryonic tissues'
        hit_list = brain, neural, cerebellum, cerebral, neuron
        (1) cehck for 'mRNA expression in embryonic tissues' for hit_list

        (2) check for 'This gene is overexpressed' for hit_list

        (1) check for 'Gene Ontology (GO) annotations related' for hit_list
