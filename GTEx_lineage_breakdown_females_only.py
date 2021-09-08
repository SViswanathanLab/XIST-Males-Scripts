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
from matplotlib.ticker import FormatStrFormatter

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.linewidth'] = 0.3
sns.set(font_scale=0.4)

figsize = (6,20)
dpi_set = 72 # change the output resolution

temp_font = 14

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/GTEx/male_and_female_XIST_expression_GTEx_summary.csv")

orig_cols = df.columns.tolist()

df = df[df['Classification']!="UNKNOWN"]
df = df[df['Gender']!="UNKNOWN"]

curr_barcode = df['Barcode'].tolist()

trunc_barcode = []

for a in curr_barcode:
    trunc_barcode.append(split_advanced(a, "-", 3)[0])

dup_items = list(set([item for item, count in collections.Counter(trunc_barcode).items() if count > 1]))

new_rows = []

for a in dup_items:
    temp_df = df[df['Barcode'].str.contains(a)]
    vals = temp_df['XIST_Expression'].tolist()
    
    if len(vals)<2:
        print("ERROR")
    
    temp_average = sum(vals)/len(vals)
    
    temp_gender = temp_df['Gender'].tolist()[0]
    temp_classification = temp_df['Classification'].tolist()[0]
    
    new_rows.append([str(a)+"-AVERAGE", temp_average, temp_classification, temp_gender]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)
df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Barcode'].str.contains(a))]

df = pd.concat([df, df_from_new_rows], ignore_index=True)

df2 = df[df['Gender']=="FEMALE"]
df = df[df['Gender']=="MALE"]

df_class = df['Classification'].tolist()
df2_class = df2['Classification'].tolist()

male_exp_unp = df['XIST_Expression'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

female_exp_unp = df2['XIST_Expression'].tolist()

female_exp = []

for a in female_exp_unp:
    if a>=0:
        female_exp.append(math.log(a+1, 2))
    else:
        female_exp.append(0)
        
all_male_lineages = list(set(df_class))
all_female_lineages = list(set(df2_class))

shared_lineages = list(set(all_male_lineages) & set(all_female_lineages))
uq_male_lineages = np.setdiff1d(all_male_lineages,all_female_lineages)
uq_female_lineages = np.setdiff1d(all_female_lineages,all_male_lineages)

N_axs = len(shared_lineages) + len(uq_female_lineages)

bins_list = np.arange(0,10.01,0.1)

fig, axs = plt.subplots(N_axs, 1, sharex=True, tight_layout=True, figsize=figsize)

a=0
while a<len(shared_lineages):
    
    df_temp = df[df['Classification']==shared_lineages[int(a)]]
    df2_temp = df2[df2['Classification']==shared_lineages[int(a)]]
    
    male_temp_unp = df_temp['XIST_Expression'].tolist()
    
    male_temp = []
    
    for b in male_temp_unp:
        if b>=0:
            male_temp.append(math.log(b+1, 2))
        else:
            male_temp.append(0)

    female_temp_unp = df2_temp['XIST_Expression'].tolist()
    
    female_temp = []
    
    for b in female_temp_unp:
        if b>=0:
            female_temp.append(math.log(b+1, 2))
        else:
            female_temp.append(0)
    
    axs[a].hist(female_temp, density=False, bins=bins_list, color="#92DCE5", linewidth=0.5, edgecolor="none")

    axs[a].set_ylabel("F-" + str(shared_lineages[int(a)]) + " (" + str(len(female_temp)) + ")", rotation=0, labelpad=30, fontsize=temp_font)
            
    axs[a].grid(False)
    axs[a].set_facecolor("white")
    axs[a].spines['bottom'].set_color('0')
    axs[a].spines['left'].set_color('0')
    axs[a].tick_params(bottom='on', left='on')

    axs[a].tick_params(axis='x', labelsize=4)
    
    #axs[a].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #axs[a+1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    
    axs[a].spines['bottom'].set_lw(0.1)
    axs[a].spines['left'].set_lw(0.1)
    
    axs[a].xaxis.set_tick_params(width=0.1, size=1, pad=1)    
    axs[a].yaxis.set_tick_params(width=0.1, size=1, pad=1)
    
    axs[a].xaxis.labelpad = 3
    axs[a].yaxis.labelpad = 82

    #axs[a].set_ylim([0.01, 20])

    axs[a].set_yticks([5.00])
    
    if str(shared_lineages[int(a)]) == "Bladder" or str(shared_lineages[int(a)]) == "Kidney":
            print("YES")
            axs[a].set_yticks([1])
    
    #axs[a].set_yscale('log')
    
    
    a#xs[a+1].set_ylim([0.01, 500])
    
    #axs[a+1].set_yticks([100])
    
    #axs[a+1].set_yscale('log')
    
    #axs[a].set_xticks([1, 10])
    #axs[a+1].set_xticks([1, 10])
    
    #axs[a].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    #axs[a+1].yaxis.set_major_formatter(FormatStrFormatter('%03d'))
    
    a=a+1
    
    print(a)


while (a-len(shared_lineages))<len(uq_female_lineages):
    
    df2_temp = df2[df2['Classification']==uq_female_lineages[a-len(shared_lineages)]]

    female_temp_unp = df2_temp['XIST_Expression'].tolist()
    
    female_temp = []
    
    for b in female_temp_unp:
        if b>=0:
            female_temp.append(math.log(b+1, 2))
        else:
            female_temp.append(0)
    
    axs[a].hist(female_temp, density=False, bins=bins_list, color="#92DCE5", linewidth=0.5, edgecolor="none")
    axs[a].set_ylabel("F-" + str(uq_female_lineages[a-len(shared_lineages)]) + " (" + str(len(female_temp)) + ")", rotation=0, labelpad=30, fontsize=temp_font)
        
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
    axs[a].yaxis.labelpad = 82

    #axs[a].set_ylim([0.01, 20])
    
    axs[a].set_yticks([2.00])
    
    #axs[a].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    
    a=a+1
    
    print(a)

    
plt.xlabel("XIST Expression (log$_2$(TPM+1))", fontsize=temp_font)

plt.xlim(0, 10)

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/FEMALES_XIST_exp_by_lineage_website_GTEx_male_and_female_histogram.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/FEMALES_XIST_exp_by_lineage_website_GTEx_male_and_female_histogram.pdf", dpi=dpi_set, bbox_inches = 'tight')




