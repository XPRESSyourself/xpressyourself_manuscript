import pandas as pd

data = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/original/ingolia_counts.txt',
    sep='\t',
    index_col=0)

data.shape
data.head()
data['name'] = data.index
data = data.groupby(['name']).agg('sum')
del data.index.name

data.shape
data.head()


data.to_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/original/ingolia_counts_deduplicated.txt',
    sep='\t')
    
