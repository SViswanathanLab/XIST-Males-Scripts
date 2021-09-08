#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 16:37:46 2021

@author: ananthansadagopan
"""

import pandas as pd
import collections
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from statannot import add_stat_annotation
import statistics
import collections
import matplotlib

overall_size=12 #font
top_font_size=7.5
plt.rcParams["font.family"] = "Arial"

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(10,5)})
sns.set(font_scale=0.8)
plt.rcParams.update({'font.size': 14})

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/Titan_Absolute_CN_ratio_data.txt", sep="\t", header=None)

df.columns = ['Sample', 'Chromosome', 'Unrounded_Ratio', 'Rounded_Ratio']

temp_vals = df['Unrounded_Ratio'].tolist()

new_vals = []

for a in temp_vals:
    new_vals.append(math.log(a, 2))
    
df['Unrounded_Ratio'] = new_vals

df = df[df['Sample']!="TCGA-HC-7740"] #XIST+ normal, XIST- tumor
df = df[df['Sample']!="TCGA-19-4065"] #XIST+ tumor is -02, CN calculated for -01

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df = df[~(df['Sample'].str.contains(a))]
    
    
uq_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

df_chr = df['Chromosome'].tolist()
df_ratio = df['Unrounded_Ratio'].tolist()

master_df = pd.DataFrame()

for a in uq_chr:
    df_temp = df[df['Chromosome']==a]
    temp_list = df_temp['Unrounded_Ratio'].tolist()
    master_df[a] = temp_list

fig = plt.figure()
    
gs = gridspec.GridSpec(1,1)

ax_stripplot = fig.add_subplot(gs[0])
plt.tick_params(bottom='on', left='on')
    
sns.stripplot(data=master_df, color="black", s=2, ax=ax_stripplot, jitter=0.2, edgecolor="none", alpha=0.9, zorder=2, linewidth=0.15)
sns.boxplot(data=master_df, ax=ax_stripplot, fliersize=10, showfliers=False, color="white", linewidth=0.5)

"""
for i,artist in enumerate(ax_stripplot.artists):
    # Set the linecolor on the artist to the facecolor, and set the facecolor to None
    col = artist.get_facecolor()
    artist.set_edgecolor(col)
    #artist.set_facecolor()
"""

ax_stripplot.set(xlabel='', ylabel='log$_2$(CN$_T$/CN$_A$)')

ax_stripplot.grid(False)
ax_stripplot.set_facecolor("white")

#ax_stripplot.axhline(2, c="black")

ax_stripplot.spines['bottom'].set_color('0')
ax_stripplot.spines['left'].set_color('0')

#ax_stripplot.spines['bottom'].set_lw()

#ax_stripplot.set_ylim([0, 2.5])
ax_stripplot.set_ylim([-1.5, 1.3])

fig.tight_layout()

dpi_set = 72

plt.tick_params(bottom='on', left='on')

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/log2_TITAN_ABSOLUTE_copy_number_ratio_jitter.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/log2_TITAN_ABSOLUTE_copy_number_ratio_jitter.pdf", dpi=dpi_set, bbox_inches = 'tight')


    