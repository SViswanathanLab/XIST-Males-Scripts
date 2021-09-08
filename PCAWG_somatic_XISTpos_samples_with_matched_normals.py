#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 12:00:30 2021

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

plt.rcParams["font.family"] = "Arial"

hfont = {'fontname':'Arial'}
sns.set(font_scale=0.7)

plt.rcParams['axes.labelsize'] = 3

figsize = (3,3) # the size of the figure - changes the shape of the squares in the comut
dpi_set = 72 # change the output resolution

x_padding = 0.01 # the x distance between patches in comut
y_padding = 0.01 # the y distance between patches in comut
tri_padding = 0.03 # the distance between triangles in comut



df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/XIST_expression_summary.xlsx")

df_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/project_code_sp", sep="\t")

ids_ref = df_ref['icgc_specimen_id'].tolist()
proj_code = df_ref['dcc_project_code'].tolist()

id_dict = dict(zip(ids_ref, proj_code))

df_barcode = df['Barcode'].tolist()

ids_to_remove = []

for a in df_barcode:
    if id_dict[a][-2:] == "US":
        ids_to_remove.append(a)
        
df = df[~(df['Barcode'].isin(ids_to_remove))]

valid_ids = ['DO45197', 'DO46877']

df_XIST = df

fig, axs = plt.subplots(1, 2, sharex=True, tight_layout=True, figsize=figsize)

q=0
for a in valid_ids:

    df_temp = df_XIST[df_XIST['Patient']==a]
     
    sample_name = df_temp['Barcode'].tolist()
    sample_exp = df_temp['XIST_FPKM_UQ'].tolist()
    
    sample_dict = dict(zip(sample_name, sample_exp))
            
    for w in sample_name:
        if "SP99117" not in w and "SP102901" not in w:
            normal_val = math.log(sample_dict[w]+1, 2)
        else:
            tumor_val = math.log(sample_dict[w]+1, 2)
    
    ind = np.arange(2)
    width = 0.5
    axs[q].bar(ind, [normal_val, tumor_val], width, color='b')

    axs[q].set_ylim([0, 8.1])
    axs[q].grid(False)
    axs[q].set_facecolor("white")
    axs[q].spines['bottom'].set_color('0')
    axs[q].spines['left'].set_color('0')
    axs[q].tick_params(bottom='on', left='on', which="both")
    
    axs[q].xaxis.set_tick_params(width=0.2, size=2, pad=2)
    axs[q].yaxis.set_tick_params(width=0.2, size=2, pad=2)
    axs[q].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")

    axs[q].xaxis.labelpad = 1
    axs[q].yaxis.labelpad = 8.5
    
    axs[q].text(0.2, 8, a,
     fontsize=9, verticalalignment='top', color="black")

    axs[q].set_xticks([0, 1])
    axs[q].set_xticklabels(['N', 'T'])

    if q>0:
        axs[q].set_yticklabels([])

    print(q)

    q=q+1
    
axs[0].set_ylabel("XIST Expression (log$_2$(FPKM-UQ+1))", size=9)

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PCAWG_somatic_XIST_exp_examples.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PCAWG_somatic_XIST_exp_examples.pdf", dpi=dpi_set, bbox_inches = 'tight')

