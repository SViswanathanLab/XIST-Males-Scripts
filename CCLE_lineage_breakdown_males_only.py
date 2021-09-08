#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 19:06:35 2021

@author: ananthansadagopan
"""

import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
import statistics
import seaborn as sns
import collections
import matplotlib.patches as mpatches
from matplotlib.ticker import FormatStrFormatter


def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.linewidth'] = 0.1
sns.set(font_scale=0.4)

figsize = (6,20)
dpi_set = 72

temp_font = 14

df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/CCLE_XIST_summary.xlsx")

val = df['XIST Expression (TPM)'].tolist()

XIST_log2 = []

for a in val:
    XIST_log2.append(math.log(a+1, 2))

df['XIST_log2TPM_plus_1'] = XIST_log2

dataset_unp = df['cell line'].tolist()

dataset = []

for a in dataset_unp:
    full_end = split_advanced(a, "_", 1)[1]
    segments = full_end.split("_")
    new_info = ""
    for w in segments:
        if w != "AND":
            new_info = new_info + w.lower().capitalize() + " "
        else:
            new_info = new_info + w.lower() + "\n"
    dataset.append(new_info)

df['Classification'] = dataset

df2 = df[df['Gender']=="FEMALE"]
df = df[df['Gender']=="MALE"]

df_class = df['Classification'].tolist()
df2_class = df2['Classification'].tolist()

all_male_lineages = list(set(df_class))
all_female_lineages = list(set(df2_class))

shared_lineages = list(set(all_male_lineages) & set(all_female_lineages))
uq_male_lineages = np.setdiff1d(all_male_lineages,all_female_lineages)
uq_female_lineages = np.setdiff1d(all_female_lineages,all_male_lineages)

N_axs = len(shared_lineages) + len(uq_male_lineages)

bins_list = np.arange(0,10.01,0.1)

fig, axs = plt.subplots(N_axs, 1, sharex=True, tight_layout=True, figsize=figsize)

a=0
while a<len(shared_lineages):
    
    df_temp = df[df['Classification']==shared_lineages[a]]
    
    male_temp_unp = df_temp['XIST Expression (TPM)'].tolist()
    
    male_temp = []
    
    for b in male_temp_unp:
        if b>=0:
            male_temp.append(math.log(b+1, 2))
        else:
            male_temp.append(0)

    axs[a].hist(male_temp, density=False, color="#B7ADCF", bins=bins_list, linewidth=0.5, edgecolor="none")

    axs[a].set_ylabel("M-" + str(shared_lineages[int(a)]) + " (" + str(len(male_temp)) + ")", rotation=0, labelpad=30, fontsize=temp_font)
        
    axs[a].grid(False)
    axs[a].set_facecolor("white")
    axs[a].spines['bottom'].set_color('0')
    axs[a].spines['left'].set_color('0')
    axs[a].tick_params(bottom='on', left='on')
 
    #plt.setp(axs[a].yaxis.get_majorticklabels(), ha="right", va="bottom")
    #plt.setp(axs[a+1].yaxis.get_majorticklabels(), ha="right", va="bottom")
    
    axs[a].tick_params(axis='x', labelsize=4)
    #axs[a].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #axs[a+1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    
    axs[a].spines['bottom'].set_lw(0.1)
    axs[a].spines['left'].set_lw(0.1)
    
    axs[a].xaxis.set_tick_params(width=0.1, size=1, pad=3)    
    axs[a].yaxis.set_tick_params(width=0.1, size=1, pad=2)
    
    axs[a].xaxis.labelpad = 3
    axs[a].yaxis.labelpad = 120

    axs[a].set_ylim([0.01, 30])
    
    axs[a].set_yticks([5])
    
    #axs[a].set_yscale('log')
    
    #axs[a+1].set_yscale('log')
    
    #axs[a].set_xticks([1, 10])
    #axs[a+1].set_xticks([1, 10])
    
    """
    for tick in axs[a].yaxis.get_major_ticks():
        tick.label.set_fontsize(1)     
        
    for tick in axs[a+1].yaxis.get_major_ticks():
        tick.label.set_fontsize(1)     
    """
    
    a=a+1
    
    print(a)


while (a-len(shared_lineages))<len(uq_male_lineages):
    
    df_temp = df[df['Classification']==uq_male_lineages[a-len(shared_lineages)]]

    male_temp_unp = df_temp['XIST Expression (TPM)'].tolist()
    
    male_temp = []
    
    for b in male_temp_unp:
        if b>=0:
            male_temp.append(math.log(b+1, 2))
        else:
            male_temp.append(0)
    
    axs[a].hist(male_temp, density=False, bins=bins_list, color="#B7ADCF", linewidth=0.5, edgecolor="none")
    axs[a].set_ylabel("M-" + str(uq_male_lineages[a-len(shared_lineages)]) + " (" + str(len(male_temp)) + ")", rotation=0, labelpad=30, fontsize=temp_font)
        
    axs[a].grid(False)
    axs[a].set_facecolor("white")
    axs[a].spines['bottom'].set_color('0')
    axs[a].spines['left'].set_color('0')
    axs[a].tick_params(bottom='on', left='on')
    
    axs[a].spines['bottom'].set_lw(0.1)
    axs[a].spines['left'].set_lw(0.1)
    
    axs[a].xaxis.set_tick_params(width=0.1, size=1, pad=1)    
    axs[a].yaxis.set_tick_params(width=0.1, size=1, pad=1)
    
    axs[a].xaxis.labelpad = 3
    axs[a].yaxis.labelpad = 120

    axs[a].set_ylim([0.01, 30])
    
    axs[a].set_yticks([5])
    
    """
    for tick in axs[a].yaxis.get_major_ticks():
        tick.label.set_fontsize(1)     
    """
    
    #axs[a].set_yscale('log')

    #axs[a].set_xticks([1, 10])

    #plt.setp(axs[a].yaxis.get_majorticklabels(), ha="right", va="bottom")

    #axs[a].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    a=a+1
    
    print(a)

    
plt.xlabel("XIST Expression (log$_2$(TPM+1))", fontsize=temp_font)

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

plt.xlim(0, 10)

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_XIST_exp_by_lineage_CCLE_male_and_female_histogram.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_XIST_exp_by_lineage_CCLE_male_and_female_histogram.pdf", dpi=dpi_set, bbox_inches = 'tight')


