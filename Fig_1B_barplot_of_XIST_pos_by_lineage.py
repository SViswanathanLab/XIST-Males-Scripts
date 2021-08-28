#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 19:06:35 2021

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

overall_size=2 #font
top_font_size=overall_size
plt.rcParams["font.family"] = "Arial"

sns.set_style(style="white")

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
figsize = (4,2.7)
sns.set(font_scale=0.4)
plt.rcParams.update({'font.size': 2})

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

df_class_temp = df['Classification'].tolist()

orig_cols = df.columns.tolist()

df_barcode = df['Barcode'].tolist()

end_barcode = []

for a in df_barcode:
    end_barcode.append(int(split_advanced(a, "-", 3)[1]))

df['End_Val'] = end_barcode

df = df[df['End_Val']<10]
df = df[df['Classification'] != "UNKNOWN"]
df = df[df['Gender'] != "UNKNOWN"]

curr_barcode = df['Barcode'].tolist()

trunc_barcode = []

for a in curr_barcode:
    trunc_barcode.append(split_advanced(a, "-", 3)[0])

dup_items = list(set([item for item, count in collections.Counter(trunc_barcode).items() if count > 1]))

new_rows = []

for a in dup_items:
    temp_df = df[df['Barcode'].str.contains(a)]
    vals = temp_df['XIST_TPM'].tolist()
    
    if len(vals)<2:
        print("ERROR")
    
    temp_average = sum(vals)/len(vals)
    temp_log2 = math.log(temp_average, 2)
    
    temp_gender = temp_df['Gender'].tolist()[0]
    temp_classification = temp_df['Classification'].tolist()[0]
    temp_secondary_classification = temp_df['Secondary_Class'].tolist()[0]
    
    new_rows.append([str(a)+"-09", temp_log2, temp_average, temp_classification, temp_gender, temp_secondary_classification, float("nan")]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)
df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Barcode'].str.contains(a))]

df = pd.concat([df, df_from_new_rows], ignore_index=True)

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

for a in invalid_ids:
    df = df[~(df['Barcode'].str.contains(a))]

df = df[df['Gender']=="MALE"]

temp_classes = df['Classification'].tolist()

secondary_classes = df['Secondary_Class'].tolist()

new_class_list = []

a=0
while a<len(temp_classes):
    if temp_classes[a] != "TGCT":
        new_class_list.append(temp_classes[a])
    else:
        if secondary_classes[a] == secondary_classes[a]:
            new_class_list.append(secondary_classes[a])
        else:
            new_class_list.append("FAKE")            
    a=a+1

df['Classification'] = new_class_list

df = df[df['Classification']!="FAKE"]

male_exp_unp = df['XIST_TPM'].tolist()

dataset = df['Classification'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

uq_datasets = list(set(dataset))
n_max = []

new_df_rows = []

temp_pcts = []

all_lengths = []

for w in uq_datasets:
    temp_list = []
    a=0
    while a<len(dataset):
        if dataset[a] == w:
            temp_list.append(male_exp[a])
        a=a+1
    above_cutoff = [x for x in temp_list if x>2]
    
    temp_length = len(temp_list)
    
    all_lengths.append(temp_length)
    
    len_above_cutoff = len(above_cutoff)
    
    median_temp = statistics.median(temp_list)
    
    try:
        median_above = statistics.median(above_cutoff)
    except:
        median_above = 0
    
    pct_val = len_above_cutoff/temp_length
    
    temp_pcts.append(pct_val)
    
    new_cutoff_val = median_temp/1000000 + pct_val
    
    n_max.append(new_cutoff_val)
    new_df_rows.append(temp_list)

pct_dict = dict(zip(uq_datasets,temp_pcts))
lengths_dict = dict(zip(uq_datasets, all_lengths))
    
master_df_unp = pd.DataFrame(new_df_rows).T
master_df_unp.columns = uq_datasets

sorted_columns = [x for _, x in sorted(zip(n_max, uq_datasets))][::-1]

master_df = pd.DataFrame()

for a in sorted_columns:
    master_df[a] = master_df_unp[a].tolist()

fig = plt.figure(figsize=figsize)

gs = gridspec.GridSpec(1,1)

ax_stripplot = fig.add_subplot(gs[0])
plt.tick_params(bottom='on', left='on')

rotate_true = 0

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    
    if not ax:
        ax = plt.gca()

    im = ax.imshow(data, **kwargs)

    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))

    ax.set_xticklabels(col_labels,size=overall_size)
    ax.set_yticklabels(row_labels, size=overall_size)

    ax.tick_params(top=False, bottom=False,
                   labeltop=False, labelbottom=False)

    if rotate_true == 0:
        plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
                 rotation_mode="anchor")
        ax.tick_params(top=False)

    else:   
        plt.setp(ax.get_xticklabels(), rotation=45, ha="left",
                 rotation_mode="anchor")

    for edge, spine in ax.spines.items():
        spine.set_color("none")

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="white", linestyle='-', linewidth=0)
    ax.grid(False)
    ax.set_facecolor("white")
        
    return im

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color="black")
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

temp_cols = master_df.columns.tolist()

new_pcts = []
new_lengths = []

for a in temp_cols:
    new_pcts.append(pct_dict[a]*100)
    new_lengths.append(lengths_dict[a])

pct_df = pd.DataFrame([new_pcts])
pct_df.index = ['XIST+ (%)']

new_names = []

a=0
while a<len(temp_cols):
    
    if temp_cols[a] != "TGCT-S" and temp_cols[a] != "TGCT-NS":
        new_names.append(temp_cols[a] + " (" + str(new_lengths[a]) + ")")
    elif temp_cols[a] == "TGCT-S":
        new_names.append("SGCT (" + str(new_lengths[a]) + ")")
    elif temp_cols[a] == "TGCT-NS":
        new_names.append("NSGCT (" + str(new_lengths[a]) + ")")
    a=a+1


temp_vals = 0

a=0
while a<len(new_pcts):
    temp_vals = temp_vals + new_pcts[a]/100*new_lengths[a]
    a=a+1

pct_df.columns = new_names

col_dict={100:"#EBEBEB", 0:"#EBEBEB"}

cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

new_vals = []
for a in new_names:
    new_vals.append(100)

new_df = pd.DataFrame(new_vals).T
temp_plot = sns.barplot(data=pct_df, color="#008BB5", lw=0)

ax_stripplot.set(xlabel='')
ax_stripplot.set_ylabel('Percentage of XIST+ Male Cancer Samples', labelpad=50, fontsize=7)

ax_stripplot.set_yticklabels([0, 20, 40, 60, 80, 100], fontsize=7)


ax_stripplot.grid(False)
ax_stripplot.set_facecolor("white")

ax_stripplot.spines['bottom'].set_color('0')
ax_stripplot.spines['left'].set_color('0')

ax_stripplot.spines['bottom'].set_lw(0.1)

ax_stripplot.spines['bottom'].set_lw(0.2)
ax_stripplot.spines['left'].set_lw(0.2)

XISTpos_patch = mpatches.Patch(color='#9D2727', label='XIST+')
XISTneg_patch = mpatches.Patch(color='#008BB5', label='XIST-')

ax_stripplot.xaxis.set_tick_params(width=0.2, size=2, pad=1)
ax_stripplot.yaxis.set_tick_params(width=0.2, size=2, pad=2)

ax_stripplot.xaxis.labelpad = 1
ax_stripplot.yaxis.labelpad = 4

plt.xticks(rotation = 90, ha="center")

ax_stripplot.set_ylim([0, 100.55])

fig.tight_layout()

dpi_set = 72

plt.tick_params(bottom='on', left='on')

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/pct_barplot_TCGA_Xena_male_XIST_above_3_TPM_by_cancer_type.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/pct_barplot_TCGA_Xena_male_XIST_above_3_TPM_by_cancer_type.pdf", dpi=dpi_set, bbox_inches = 'tight')




