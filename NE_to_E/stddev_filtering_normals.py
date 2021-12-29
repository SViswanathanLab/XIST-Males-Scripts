#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 25 19:06:23 2021

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
import scipy
import scipy.stats as stats
from scipy.stats import skewtest

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.linewidth'] = 0.3

figsize = (4.11,3.52) # the size of the figure - changes the shape of the squares in the comut
dpi_set = 72 # change the output resolution
sns.set(font_scale=0.6)

name = "/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/skewness_male_normals_normalized_NE_expression_TPM_geq200"

testing_val = "Male"
gene_type = "Non-Escapee"
cutoff = 0.125 #0.125 for Escapee, 0.25 for Non-Escape




df_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

df = df_ref[df_ref['Classification']!="UNKNOWN"]

df_barcode = df['Barcode'].tolist()

end_barcode = []

for a in df_barcode:
    end_barcode.append(int(split_advanced(a, "-", 3)[1]))

df['End_Val'] = end_barcode

df = df[df['End_Val']==11]

del df['End_Val']

temp_class = []

orig_cols = df.columns.tolist()

curr_barcode = df['Barcode'].tolist()

"""
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
    temp_second_class = temp_df['Secondary_Class'].tolist()[0]
    
    new_rows.append([str(a)+"-09", temp_log2, temp_average, temp_classification, temp_gender, temp_second_class, "UNKNOWN"]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)

df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Barcode'].str.contains(a))]

df_ref = pd.concat([df, df_from_new_rows], ignore_index=True)
"""

df_ref = df

bar = df_ref['Barcode'].tolist()
XIST_TPM = df_ref['XIST_TPM'].tolist()
class_ref = df_ref['Classification'].tolist()
sex_ref = df_ref['Gender'].tolist()
second_class_ref = df_ref['Secondary_Class'].tolist()

TPM_dict = dict(zip(bar, XIST_TPM))
class_dict = dict(zip(bar, class_ref))
sex_dict = dict(zip(bar, sex_ref))
second_class_dict = dict(zip(bar, second_class_ref))

uq_classes = list(set(class_ref))


df = pd.read_csv(name + ".csv")
cols = df.columns.tolist()

old_cols = df.iloc[:, 0].tolist()

del cols[0]

stats = []

valid_classes = []

output_df = pd.DataFrame()

for a in uq_classes:
    valid_ids = []
    for b in cols:
        if class_dict[b] == a:
            valid_ids.append(b)   
    temp_df = df[valid_ids]
    
    if temp_df.empty:
        continue
    
    valid_classes.append(a)
    
    vals = temp_df.values.tolist()

    b=0
    used_vals = []
        
    while b<len(vals):
        temp_list = vals[b]
        del temp_list[0]
        temp_list = [x for x in temp_list if x == x]
        sem_test = scipy.stats.sem(temp_list)
        if sem_test > cutoff:
            used_vals.append(b)
        b=b+1
    
    temp_df.loc[(np.array(used_vals)), :] = np.nan
    
    output_df = pd.concat([output_df, temp_df], axis=1)
    
output_df.index = old_cols

output_df.to_csv(name + "_stddev_filtered_" + str(cutoff) + ".csv")

df = output_df

cols = df.columns.tolist()

male_valid = ['BLCA', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'PRAD']

female_valid = ['UVM', 'OV', 'SKCM', 'ACC', 'KICH', 'KIRC', 'KIRP', 'UCS', 'HNSC', 'GBM', 'COAD', 'LUSC', 'LUAD', 'SARC', 'STAD', 'LAML', 'BRCA', 'UCEC', 'CESC']

valid_ids = []

if testing_val == "Male":
    for a in cols:
        if class_dict[a] in male_valid:
            valid_ids.append(a)
else:
    for a in cols:
        if class_dict[a] in female_valid:
            valid_ids.append(a)        

cols = list(set(cols) & set(valid_ids))

stats = []

for a in cols:

    curr_vals = [x for x in df[a].tolist() if x == x]
    
    try:
        c = scipy.stats.skew(curr_vals)
        stats.append(c)
    except ValueError:
        pass
    
    
fig, ax = plt.subplots(1, 1, sharex=True, tight_layout=True, figsize=figsize)
    
bins_list = np.arange(-10,10,0.1)

ax.hist(stats, density=False, bins=bins_list, color="black", linewidth=0.5, edgecolor="none")

ax.grid(False)
ax.set_facecolor("white")
ax.spines['bottom'].set_color('0')
ax.spines['left'].set_color('0')

ax.tick_params(bottom='on', left='on', which="both")

ax.set_xlabel("Fisher-Pearson Coefficient of Skewness for\nNormalized " + gene_type + " Expression in " + testing_val + " Normals")
ax.set_ylabel("Number of Samples")

fig.savefig(name + "_NORMALS_stddev_filtered_" + str(cutoff) + ".pdf", dpi=dpi_set, bbox_inches = 'tight')





