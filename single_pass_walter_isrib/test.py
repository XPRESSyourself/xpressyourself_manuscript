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

"""
READ IN DATA
"""
# Read in single-pass XPRESSpipe processed read data
data = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/single_pass_walter_isrib/isrib_merged_counts_htseq_count_table.tsv',
    sep = '\t',
    index_col = 0)

data.head()


data_named = convert_names(
    data,
    '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.95.gtf')


data2  = data_named.copy()

data_named.shape

data2.loc[~(data2==0).all(axis=1)].shape
data2[data2.min(axis=1) > 3].shape
thresh = data2[data2.min(axis=1) > 5]
thresh_rpm = rpm(thresh)

plt.close
ref_line = np.log10(thresh_rpm + 1).mean().mean()
fig, ax = plt.subplots(1,1, figsize=(12.5,5))
plt.grid(False)
ax = sns.violinplot(data=np.log10(thresh_rpm + 1))
plt.xticks(rotation=45)
ax.set_ylabel('Normalized Counts (RPM)')
ax.set_xlabel('Samples')
ax.axhline(ref_line, ls='--', color='grey')
plt.show()
