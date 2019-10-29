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
import os
import sys
import pandas as pd
import numpy as np
import matplotlib
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


"""Collapse technical replicates and rename samples
"""
def name_map(
        data,
        map):

    name_dict = pd.Series(map.ingolia_name.values,index=map.Run).to_dict()
    data_c = data.copy()
    data_c.columns = data_c.columns.to_series().map(name_dict)
    data_sum = data_c.groupby(data_c.columns, axis=1).sum()

    return data_sum

"""Make supplemental figure 4
"""
def make_paralogs_analysis(
    data,
    original_data,
    title):

    # Init plotting space
    fig, axes = plt.subplots(
        nrows = 2,
        ncols = 2,
        figsize = (20, 20),
        subplot_kw = {
            'facecolor':'none'},
        sharey=True,
        sharex=True) # Create shared axis for cleanliness
    plt.subplots_adjust(
        bottom = 0.1)
    plt.yticks([0,1,2,3,4,5]) # Limit axis labels to ints
    plt.xticks([0,1,2,3,4,5])

    x = 0
    file_number = 0
    file_list = [
        'ribo_untr_a',
        'ribo_tm_a',
        'untr_a_hek',
        'tm_a_hek'] # Designate sample order

    for y in range(4):

        # Get data as array-like for samples being compared
        data_c1 = data.copy()
        data_c2 = original_data.copy()

        sample_a = data_c1[file_list[file_number]].values.tolist()
        sample_a = [x + 1 for x in sample_a]
        sample_a = np.array(sample_a).astype(np.float)
        sample_a = np.ndarray.tolist(sample_a)

        sample_b = data_c2[file_list[file_number]].values.tolist()
        sample_b = [x + 1 for x in sample_b]
        sample_b = np.array(sample_b).astype(np.float)
        sample_b = np.ndarray.tolist(sample_b)

        # Determine subplot location
        if file_number in [0,1]:
            ax_x = file_number % 2
            ax_y = 0
        elif file_number in [2,3]:
            ax_x = file_number % 2
            ax_y = 1
        else:
            print('oops')

        # Plot data
        axes[ax_y, ax_x].tick_params(axis='both', labelsize=32)
        axes[ax_y, ax_x].scatter(np.log10(sample_a), np.log10(sample_b), s=5,c='grey')

        genes = ['EIF3C','NOMO2','KCNQ2','RPL39','TUBA1A']
        colors = ['red','blue','green','purple','orange']
        for x in range(5):

            s1 = (data_c1 + 1).at[genes[x], file_list[file_number]]
            s2 = (data_c2 + 1).at[genes[x], file_list[file_number]]

            axes[ax_y, ax_x].scatter(np.log10(s1), np.log10(s2), s=100,c=colors[x])

        axes[ax_y, ax_x].axhline(0, ls='-', color='black', xmin=0.0457, xmax=1) # Create axis lines
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymin=0.0457, ymax=1)
        file_number += 1 # Plot counter


    # Create shared row/column titles
    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(counts)', fontsize=32)

    cols = ['Untreated','Tm']
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=40)
        ax.xaxis.set_label_position('top')

    rows = ['Ribo-seq','RNA-seq']
    for ax, rows in zip(axes[:,1], rows):
        ax.set_ylabel(rows, fontsize=40, labelpad=40, rotation=270)
        ax.yaxis.set_label_position('right')

    for ax in axes[1]:
        ax.set_xlabel('log$_1$$_0$(counts)', fontsize=32)
        ax.xaxis.set_label_position('bottom')

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

# Make figure 2B
def make_external_correlations(
        data,
        original_data,
        title):

    fig, axes = plt.subplots(
        nrows = 2,
        ncols = 2,
        figsize = (20, 20),
        subplot_kw = {
            'facecolor':'none'},
        sharey=True) # Create shared axis for cleanliness
    plt.subplots_adjust(
        bottom = 0.1)
    plt.yticks([0,1,2,3,4,5]) # Limit axis labels to ints
    plt.xticks([0,1,2,3,4,5])
    plt.axis('equal')

    x = 0
    file_number = 0
    file_list = [
        'ribo_tm_a',
        'isrib_b_hek'] # Designate sample order

    for y in range(2):

        # Get data as array-like for samples being compared
        data_c1 = data.copy()
        data_c2 = original_data.copy()

        sample_a = data_c1[file_list[file_number]].values.tolist()
        sample_a = [x + 1 for x in sample_a]
        sample_a = np.array(sample_a).astype(np.float)
        sample_a = np.ndarray.tolist(sample_a)

        sample_b = data_c2[file_list[file_number]].values.tolist()
        sample_b = [x + 1 for x in sample_b]
        sample_b = np.array(sample_b).astype(np.float)
        sample_b = np.ndarray.tolist(sample_b)

        # Run Spearman R linreg for non-normal data
        rho, p_value = stats.spearmanr(sample_a, sample_b)

        # Determine subplot location
        if file_number in [0,1]:
            ax_x = file_number % 2
            ax_y = 0
        else:
            print('oops')

        # Format p value
        if p_value < 0.001:
            p_val = '< 0.001'
        else:
            p_val = '= {:.3f}'.format(round(p_value, 3))

        rho = '{:.3f}'.format(round(rho, 3))

        # Plot data
        axes[ax_y, ax_x].tick_params(axis='both', labelsize=32)
        axes[ax_y, ax_x].scatter(np.log10(sample_a), np.log10(sample_b), s=1,c='black')
        axes[ax_y, ax_x].set_title(r"$\rho$" + ' = ' + str(rho) + '\nP ' + p_val, y=0.1, x=0.9, fontsize=32) # Format titles
        axes[ax_y, ax_x].axhline(0, ls='-', color='black', xmin=0.0457, xmax=1) # Create axis lines
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymin=0.0457, ymax=1)
        file_number += 1 # Plot counter
        x += 2
        print(rho)

    # Create shared row/column titles
    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(counts)', fontsize=32)

    cols = ['RPF Tm RepA','mRNA ISRIB RepB']
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=40)
        ax.xaxis.set_label_position('top')
    cols = ['log$_1$$_0$(counts)','log$_1$$_0$(counts)']
    for ax, col in zip(axes[1], cols):
        ax.set_xlabel(col, fontsize=32, labelpad=30)
        ax.xaxis.set_label_position('top')

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

def make_external_correlations_supplement(
        data,
        original_data,
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
            merged_best = pd.concat([data[[str(x)]], original_data[[str(x)]]], axis=1, sort=False)
            merged_best.columns = ['xpresspipe', 'original']
            merged_best['genes'] = merged_best.index
            sc = px.scatter(
                merged_best,
                x='xpresspipe',
                y='original',
                hover_name='genes',
                log_x=True,
                log_y=True,
                opacity=0.4,
                width=1400,
                height=1000,
                title=str(x))

            py.offline.plot(sc, filename='/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/interactive/interactive' + str(x) + '.html')

        # Get data as array-like for samples being compared
        data_c1 = data.copy()
        data_c2 = original_data.copy()

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
        if p_value < 0.001:
            p_val = '< 0.001'
        else:
            p_val = '= {:.3f}'.format(round(p_value, 3))

        rho = '{:.3f}'.format(round(rho, 3))

        # Plot data
        axes[ax_y, ax_x].tick_params(axis='both', labelsize=32)
        axes[ax_y, ax_x].scatter(np.log10(sample_a), np.log10(sample_b), s=1,c='black')
        axes[ax_y, ax_x].set_title(r"$\rho$" + ' = ' + str(rho) + '\nP ' + p_val, y=0.1, x=0.8, fontsize=20) # Format titles
        axes[ax_y, ax_x].axhline(0, ls='-', color='black', xmin=0.05, xmax=1) # Create axis lines
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymin=0.05, ymax=0.88)
        file_number += 1 # Plot counter
        print(rho)

    # Create shared row/column titles
    count_label = ['log$_1$$_0$(counts)','log$_1$$_0$(counts)','log$_1$$_0$(counts)','log$_1$$_0$(counts)']

    cols = ['RPF RepA','RPF RepB','mRNA RepA','mRNA RepB']
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=40)
        ax.xaxis.set_label_position('top')

    for ax in axes[3]:
        ax.set_xlabel('log$_1$$_0$(counts)', fontsize=32)
        ax.xaxis.set_label_position('bottom')

    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(counts)', fontsize=32)

    rows = ['Untreated','Tm','Tm + ISRIB','ISRIB']
    for ax, rows in zip(axes[:,3], rows):
        ax.set_ylabel(rows, fontsize=40, labelpad=40, rotation=270)
        ax.yaxis.set_label_position('right')

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')


def make_spearman_pearson_correlations(
        data,
        original_data,
        title):

    data_t = np.log10(data + 1)
    original_data_t = np.log10(original_data + 1)

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

    # Run Spearman on raw counts
    spear_list = []
    for x in file_list:

        # Get data as array-like for samples being compared
        data_c1 = data.copy()
        data_c2 = original_data.copy()

        sample_a = data_c1[x].values.tolist()
        sample_a = [x + 1 for x in sample_a]
        sample_a = np.array(sample_a).astype(np.float)
        sample_a = np.ndarray.tolist(sample_a)

        sample_b = data_c2[x].values.tolist()
        sample_b = [x + 1 for x in sample_b]
        sample_b = np.array(sample_b).astype(np.float)
        sample_b = np.ndarray.tolist(sample_b)

        # Run Spearman R linreg for non-normal data
        rho, spear_p_value = stats.spearmanr(sample_a, sample_b)
        spear_list.append(rho)

    # Run Pearson on log transformed data
    pear_list = []
    for x in file_list:

        # Get data as array-like for samples being compared
        data_c1 = data_t.copy()
        data_c2 = original_data_t.copy()

        sample_a = data_c1[x].values.tolist()
        sample_a = [x + 1 for x in sample_a]
        sample_a = np.array(sample_a).astype(np.float)
        sample_a = np.ndarray.tolist(sample_a)

        sample_b = data_c2[x].values.tolist()
        sample_b = [x + 1 for x in sample_b]
        sample_b = np.array(sample_b).astype(np.float)
        sample_b = np.ndarray.tolist(sample_b)

        r, pear_p_value = stats.pearsonr(sample_a, sample_b)
        pear_list.append(r)

    spear_list_reps = []
    file_number = 0
    for x in range(int(len(file_list)/2)):

        # Get data as array-like for samples being compared
        data_c1 = data.copy()

        sample_a = data_c1[file_list[file_number]].values.tolist()
        sample_a = [x + 1 for x in sample_a]
        sample_a = np.array(sample_a).astype(np.float)
        sample_a = np.ndarray.tolist(sample_a)

        sample_b = data_c1[file_list[file_number+1]].values.tolist()
        sample_b = [x + 1 for x in sample_b]
        sample_b = np.array(sample_b).astype(np.float)
        sample_b = np.ndarray.tolist(sample_b)

        file_number += 2

        # Run Spearman R linreg for non-normal data
        rho, spear_p_value = stats.spearmanr(sample_a, sample_b)
        spear_list_reps.append(rho)


    pear_list_reps = []
    file_number = 0
    for x in range(int(len(file_list)/2)):

        # Get data as array-like for samples being compared
        data_c1 = data_t.copy()

        sample_a = data_c1[file_list[file_number]].values.tolist()
        sample_a = [x + 1 for x in sample_a]
        sample_a = np.array(sample_a).astype(np.float)
        sample_a = np.ndarray.tolist(sample_a)

        sample_b = data_c1[file_list[file_number+1]].values.tolist()
        sample_b = [x + 1 for x in sample_b]
        sample_b = np.array(sample_b).astype(np.float)
        sample_b = np.ndarray.tolist(sample_b)

        file_number += 2

        r, pear_p_value = stats.pearsonr(sample_a, sample_b)
        pear_list_reps.append(r)

    data = [
        spear_list_reps,
        pear_list_reps,
        spear_list,
        pear_list]

    names = [
        'Spearman Replicates',
        'Pearson Replicates',
        'Spearman Method Comparison',
        'Pearson Method Comparison']

    counter = 0
    df = pd.DataFrame(columns = ['Value', 'Type'])

    for x in range(4):
        l = data[x]
        n = names[x]

        for y in l:
            df.loc[counter] = [float(y), str(n)]
            counter += 1

    df1 = df.loc[df['Type'].str.contains('Replicates')]
    df1['Type'] = df1['Type'].str.split(' ').str[0]
    df1.columns = ['Coefficient Value', 'Replicates']
    ax1 = sns.violinplot(x="Coefficient Value", y="Replicates", data=df1)
    plt.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/replicates_' + str(title),
        dpi = 1800,
        bbox_inches = 'tight')
    plt.close()

    df2 = df.loc[df['Type'].str.contains('Method')]
    df2['Type'] = df2['Type'].str.split(' ').str[0]
    df2.columns = ['Coefficient Value', 'Method Comparison']
    ax2 = sns.violinplot(x="Coefficient Value", y="Method Comparison", data=df2)
    plt.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/method_' + str(title),
        dpi = 1800,
        bbox_inches = 'tight')
    plt.close()

"""Make Figure 2A
"""
def make_internal_correlations(
        data,
        title):

    fig, axes = plt.subplots(
        nrows = 2,
        ncols = 2,
        figsize = (20, 20),
        subplot_kw = {
            'facecolor':'none'}) # Create shared axis for cleanliness
    plt.subplots_adjust(
        bottom = 0.1)
    plt.yticks([0,1,2,3,4,5]) # Limit axis labels to ints
    plt.xticks([0,1,2,3,4,5])

    x = 0
    file_number = 0
    file_list = [
        'ribo_tm_a',
        'ribo_tm_b',
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
        if file_number in [0,1]:
            ax_x = file_number % 2
            ax_y = 0
        else:
            print('oops')

        # Format p value
        if p_value < 0.001:
            p_val = '< 0.001'
        else:
            p_val = '= {:.3f}'.format(round(p_value, 3))

        rho = '{:.3f}'.format(round(rho, 3))

        # Plot data
        axes[ax_y, ax_x].tick_params(axis='both', labelsize=32)
        axes[ax_y, ax_x].scatter(np.log10(sample_a), np.log10(sample_b), s=1,c='black')
        axes[ax_y, ax_x].set_title(r"$\rho$" + ' = ' + str(rho) + '\nP ' + p_val, y=0.1, x=0.9, fontsize=32) # Format titles
        axes[ax_y, ax_x].axhline(0, ls='-', color='black', xmin=0.0457, xmax=1) # Create axis lines
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymin=0.0457, ymax=1)
        file_number += 1 # Plot counter
        x += 2
        print(rho)

    # Create shared row/column titles

    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(counts)', fontsize=32)

    cols = ['RPF Tm','mRNA ISRIB']
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=40)
        ax.xaxis.set_label_position('top')
    cols = ['log$_1$$_0$(counts)','log$_1$$_0$(counts)']
    for ax, col in zip(axes[1], cols):
        ax.set_xlabel(col, fontsize=32, labelpad=30)
        ax.xaxis.set_label_position('top')

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

def make_internal_correlations_supplement(
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
        if p_value < 0.001:
            p_val = '< 0.001'
        else:
            p_val = '= {:.3f}'.format(round(p_value, 3))

        rho = '{:.3f}'.format(round(rho, 3))

        # Plot data
        axes[ax_y, ax_x].tick_params(axis='both', labelsize=32)
        axes[ax_y, ax_x].scatter(np.log10(sample_a), np.log10(sample_b), s=1,c='black')
        axes[ax_y, ax_x].set_title(r"$\rho$" + ' = ' + str(rho) + '\nP ' + p_val, y=0.1, x=0.8, fontsize=20) # Format titles
        axes[ax_y, ax_x].axhline(0, ls='-', color='black', xmin=0.05, xmax=1) # Create axis lines
        axes[ax_y, ax_x].axvline(0, ls='-', color='black', ymin=0.05, ymax=1)
        file_number += 1 # Plot counter
        x += 2
        print(rho)

    # Create shared row/column titles

    for ax in axes[:,0]:
        ax.set_ylabel('log$_1$$_0$(counts)', fontsize=32)

    cols = ['RPF Untreated','mRNA Untreated','RPF Tm','mRNA Tm']
    for ax, col in zip(axes[0], cols):
        ax.set_xlabel(col, fontsize=32)
        ax.xaxis.set_label_position('top')

    cols = ['RPF Tm + ISRIB','mRNA Tm + ISRIB','RPF ISRIB','mRNA ISRIB']
    for ax, col in zip(axes[1], cols):
        ax.set_xlabel(col, fontsize=32)
        ax.xaxis.set_label_position('top')

    fig.savefig(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/plots/' + str(title),
        dpi = 600,
        bbox_inches = 'tight')

# Read in gene dictionary
def make_ref(
        ref_file,
        gene_list):

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

"""Read in single-pass XPRESSpipe processed read data quantified with HTSeq using non-de-deuplicated alignments
"""
def get_data(
        file,
        sample_suffix='_1_Aligned'):

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
        '/Users/jordan/Desktop/reference/Homo_sapiens.GRCh38.98.gtf')
    data.columns = data.columns.str.replace(str(sample_suffix), '')
    data.shape

    # Combine lanes
    sra_info = pd.read_csv(
        '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/GSE65778_table.txt',
        sep = '\t')
    data = name_map(data, sra_info)
    data = data.groupby(level=0).sum() # Combine duplicate named genes

    return data

"""Run correlations between sample alignments
"""
def cycle_external_correlations_supplement(
        file,
        original):

    try:
        data = get_data(file, sample_suffix='_1_Aligned')
    except:
        data = get_data(file, sample_suffix='__Aligned')

    data_threshold = data.loc[data[['untr_a_hek', 'untr_b_hek', 'tm_a_hek', 'tm_b_hek', 'tmisrib_a_hek', 'tmisrib_b_hek', 'isrib_a_hek', 'isrib_b_hek']].min(axis=1) >= 10] # Apply threshold to data
    data_rpm = rpm(data_threshold)
    data_genes = data_rpm.index.tolist()

    original_rpm = rpm(original)
    original_genes = original_rpm.index.tolist()

    common_genes = list(set(data_genes).intersection(original_genes))

    data_common = data_rpm.reindex(index = common_genes)
    original_common = original_rpm.reindex(index = common_genes)

    make_external_correlations_supplement(
        data_common,
        original_common,
        str(file.split('/')[-1][:-4]) + '_external_correlations_summary_htseq_all.png')

"""Read in data
"""
# Input file is protein coding only and truncated, not parsed for longest transcript only
file = '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_comp_test/isrib_comp_v98_truncated_count_table.tsv'
data = get_data(file, sample_suffix='_1_Aligned')

# Clean up data
data_threshold_deseq = data.loc[data[['untr_a_hek', 'untr_b_hek', 'tm_a_hek', 'tm_b_hek', 'tmisrib_a_hek', 'tmisrib_b_hek', 'isrib_a_hek', 'isrib_b_hek']].min(axis=1) >= 25] # Apply threshold to data

data_threshold = data.loc[data[['untr_a_hek', 'untr_b_hek', 'tm_a_hek', 'tm_b_hek', 'tmisrib_a_hek', 'tmisrib_b_hek', 'isrib_a_hek', 'isrib_b_hek']].min(axis=1) >= 10] # Apply threshold to data

data_rpm = rpm(data_threshold)

"""Export for DESeq2
"""
deseq_directory = '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_de/xpresspipe_data_deseq2/'

data_threshold_deseq.to_csv(
    deseq_directory + 'ISRIBxpresspipe_thresholded_counts.tsv',
    sep='\t')

untr_mrna = ['untr_a_hek', 'untr_b_hek']
untr_ribo = ['ribo_untr_a', 'ribo_untr_b']
tm_mrna = ['tm_a_hek', 'tm_b_hek']
tm_ribo = ['ribo_tm_a', 'ribo_tm_b']
tmisrib_mrna = ['tmisrib_a_hek', 'tmisrib_b_hek']
tmisrib_ribo = ['ribo_tmisrib_a', 'ribo_tmisrib_b']
isrib_mrna = ['isrib_a_hek', 'isrib_b_hek']
isrib_ribo = ['ribo_isrib_a', 'ribo_isrib_b']

data_threshold_deseq[untr_ribo + tm_ribo].to_csv(
    deseq_directory + 'tm_ribo_counts.tsv',
    sep='\t')
data_threshold_deseq[untr_ribo + tmisrib_ribo].to_csv(
    deseq_directory + 'tmisrib_ribo_counts.tsv',
    sep='\t')
data_threshold_deseq[untr_ribo + isrib_ribo].to_csv(
    deseq_directory + 'isrib_ribo_counts.tsv',
    sep='\t')

data_threshold_deseq[untr_mrna + tm_mrna].to_csv(
    deseq_directory + 'tm_rna_counts.tsv',
    sep='\t')
data_threshold_deseq[untr_mrna + tmisrib_mrna].to_csv(
    deseq_directory + 'tmisrib_rna_counts.tsv',
    sep='\t')
data_threshold_deseq[untr_mrna + isrib_mrna].to_csv(
    deseq_directory + 'isrib_rna_counts.tsv',
    sep='\t')

data_threshold_deseq[untr_ribo + tm_ribo + untr_mrna + tm_mrna].to_csv(
    deseq_directory + 'tm_counts.tsv',
    sep='\t')
data_threshold_deseq[untr_ribo + tmisrib_ribo + untr_mrna + tmisrib_mrna].to_csv(
    deseq_directory + 'tmisrib_counts.tsv',
    sep='\t')
data_threshold_deseq[untr_ribo + isrib_ribo + untr_mrna + isrib_mrna].to_csv(
    deseq_directory + 'isrib_counts.tsv',
    sep='\t')

"""Read in original data
"""
# Read in original raw counts from Elife supplement table
original = pd.read_csv(
    '/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/ingolia_counts_table.txt',
    sep = '\t',
    index_col = 0)
original = original.drop('size', axis = 1)

# Clean up data
original = original.groupby(level=0).sum() # Combine duplicate named genes
original_rpm = rpm(original)

"""Get common genes between datasets
"""
# Clean NaNs
data_rpm.shape
original_rpm.shape

data_rpm = data_rpm.dropna()
original_rpm = original_rpm.dropna()

data_rpm.shape
original_rpm.shape

# Make dataframes for each dataset for matching genes
data_genes = data_rpm.index.tolist()
original_genes = original_rpm.index.tolist()
len(data_genes)
len(original_genes)

common_genes = list(set(data_genes).intersection(original_genes))
len(common_genes)

data_common = data_rpm.reindex(index = common_genes)
original_common = original_rpm.reindex(index = common_genes)

"""Over-representation analysis
TopHat2 is .48% more likely to align a read that should be counted as ambiguous
https://www.nature.com/articles/nmeth.4106, supplementary table 5
"""
tophat_overcount_rate = '0.48%'

#RPF Untr A -> started with ~135,000,000 reads
rpf_start_reads = 135000000 #approx

rpf_test = pd.concat([data_common[['ribo_untr_a']], original_common[['ribo_untr_a']]], axis=1, sort=False)
rpf_test.columns = ['xpresspipe','original']
rpf_test = rpf_test[(rpf_test['xpresspipe'] >= 10) & (rpf_test['original'] >= 10)]

rpf_test.shape
rpf_test.sum()
rpf_log = np.log10(rpf_test+1)

plt.scatter(rpf_log['xpresspipe'].values, rpf_log['original'], s=1)
x, y = [0, 4.9], [-0.05, 5.05]
plt.plot(x, y, c='black')

m = (y[1] - y[0]) / (x[1] - x[0])
b = y[0]

rpf_above = rpf_log.loc[rpf_log['original'] > ((rpf_log['xpresspipe'] * m) + b)]
rpf_above.shape
rpf_above['diff'] = 10**rpf_above['original'] - 10**rpf_above['xpresspipe']
overcounts = rpf_above['diff'].sum()
overcount_rate = overcounts / rpf_start_reads
print('RPF_Untr_A')
print('# Reads overcounted: ' + str(int(round(overcounts))))
print('Overcount rate: ' + str(overcount_rate * 100) + '%')
print('TopHat2 vs STAR overcount ambiguous rate increase: ' + tophat_overcount_rate)

#RNA Untr A
rna_start_reads = 48000000

rna_test = pd.concat([data_common[['untr_a_hek']], original_common[['untr_a_hek']]], axis=1, sort=False)
rna_test.columns = ['xpresspipe','original']
rna_test.shape
rna_test.sum()

rna_log = np.log10(rna_test+1)

plt.scatter(rna_log['xpresspipe'].values, rna_log['original'], s=1)
x, y = [0, 4.9], [0.18, 5.03]
plt.plot(x, y, c='black')

m = (y[1] - y[0]) / (x[1] - x[0])
b = y[0]

rna_above = rna_log.loc[rna_log['original'] > ((rna_log['xpresspipe'] * m) + b)]
print(rna_above.shape)
rna_above['diff'] = 10**rna_above['original'] - 10**rna_above['xpresspipe']
overcounts = rna_above['diff'].sum()

rna_above['diff'].index.tolist()


plt.scatter(rna_log['xpresspipe'].loc[rna_above['diff'].index.tolist()].values, rna_log['original'].loc[rna_above['diff'].index.tolist()], s=5)
x, y = [0, 4.9], [0.18, 5.03]
plt.plot(x, y, c='red')

len(rna_above['diff'].index.tolist())

overcount_rate = overcounts / rna_start_reads
print('RNA_Untr_A')
print('# Reads overcounted: ' + str(int(round(overcounts))))
print('Overcount rate: ' + str(overcount_rate * 100) + '%')
print('TopHat2 vs STAR overcount ambiguous rate increase: ' + tophat_overcount_rate)

# Some further analysis
len(rpf_above['diff'].index.tolist())
len(rna_above['diff'].index.tolist())

common_above = list(set(rpf_above['diff'].index.tolist()).intersection(rna_above['diff'].index.tolist()))
len(common_above)

from urllib.request import Request, urlopen

target = "important paralog"
counter = 0
gene_number = 0
for gene in common_above:

    gene_number += 1

    url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + gene
    req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    webpage = urlopen(req).read()
    page = webpage.decode('utf-8').splitlines()

    if any(target in s for s in page) == True:
        counter += 1
        print(gene_number)

print(counter)
print(len(common_above))

"""Internal correlations
"""
# Run correlations between sample alignments
make_internal_correlations(
    data_rpm,
    'internal_correlations_summary_htseq_protein_truncated_v98.png')

make_internal_correlations_supplement(
    data_rpm,
    'internal_correlations_summary_htseq_protein_truncated_v98_all.png')

make_internal_correlations_supplement(
    original_rpm,
    'internal_correlations_summary_ingolia.png')

"""External correlations
"""
file_list = [
    'isrib_comp_v98_truncated_count_table.tsv',
    'isrib_comp_v98_normal_count_table.tsv',
    'isrib_comp_v98_longest_truncated_count_table.tsv',
    'isrib_comp_v72_normal_count_table.tsv',
    'isrib_comp_v72_longest_truncated_count_table.tsv',
    'isrib_comp_v72_truncated_count_table.tsv']

file_list = ['/Users/jordan/Desktop/xpressyourself_manuscript/isrib_analysis/isrib_comp_test/' + str(f) for f in file_list]

for file in file_list:

    cycle_external_correlations_supplement(file, original)

make_external_correlations(
    data_common,
    original_common,
    'external_correlations_summary_htseq_protein_truncated_v98.png')

make_external_correlations_supplement(
    data_common,
    original_common,
    'external_correlations_summary_htseq_protein_truncated_v98_all.png')

make_external_correlations_supplement(
    data_common,
    original_common,
    '',
    interactive=True)

make_spearman_pearson_correlations(
    data_common,
    original_common,
    'spearman_pearson_correlations_htseq_protein_truncated_v98.png')

make_paralogs_analysis(
    data_common,
    original_common,
    'highlight_paralogs_comparison_isrib.png')
