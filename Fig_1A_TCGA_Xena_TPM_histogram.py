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
from matplotlib.ticker import StrMethodFormatter, NullFormatter, ScalarFormatter, FormatStrFormatter
import matplotlib

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.linewidth'] = 0.3

figsize = (4.11,3.52)
dpi_set = 72
sns.set(font_scale=0.6)

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

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
    
    new_rows.append([str(a)+"-09", temp_log2, temp_average, temp_classification, temp_gender, temp_classification, "UNKNOWN"]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)

df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Barcode'].str.contains(a))]

df = pd.concat([df, df_from_new_rows], ignore_index=True)

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

#Not Biological Males: 'TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM'
#Unclassified: 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862'

for a in invalid_ids:
    df = df[~(df['Barcode'].str.contains(a))]

df2 = df[df['Gender']=="FEMALE"]
df = df[df['Gender']=="MALE"]

df_class = df['Classification'].tolist()

df_second_class = df['Secondary_Class'].tolist()

df2_class = df['Classification'].tolist()

male_exp_unp = df['XIST_TPM'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

TGCT_exp = []
PanCan_exp = []
a=0
while a<len(male_exp_unp):
    if df_class[a] == "TGCT":
        if male_exp_unp[a] >= 0:
            TGCT_exp.append(math.log(male_exp_unp[a]+1, 2))
        else:
            TGCT_exp.append(0)
    else:
        if male_exp_unp[a] >= 0:
            PanCan_exp.append(math.log(male_exp_unp[a]+1, 2))
        else:
            PanCan_exp.append(0)        
    a=a+1

female_exp_unp = df2['XIST_TPM'].tolist()

female_exp = []

for a in female_exp_unp:
    if a>=0:
        female_exp.append(math.log(a+1, 2))
    else:
        female_exp.append(0)

n_bins = 55

fig, axs = plt.subplots(4, 1, sharex=True, tight_layout=True, figsize=figsize)

bins_list = np.arange(0,12.51,0.1)

axs[3].hist(male_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")

axs[1].hist(female_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")

tgct_yes = 0
tgct_no = 0

other_yes_male = 0
other_no_male = 0

a=0
while a<len(male_exp):
    if df_class[a] == "TGCT":
        if male_exp[a] >= 2:
            tgct_yes = tgct_yes + 1
        else:
            tgct_no = tgct_no + 1
    else:
        if male_exp[a] >= 2:
            other_yes_male = other_yes_male + 1
        else:
            other_no_male = other_no_male + 1
    a=a+1

male_XISTpos_exp = []

a=0
while a<len(male_exp):
    if df_second_class[a] == "TGCT-S" and male_exp[a] >= 2:
        male_XISTpos_exp.append(male_exp[a])
    a=a+1

other_yes_female = 0
other_no_female = 0

a=0
while a<len(female_exp):
    if female_exp[a] >= 2:
        other_yes_female = other_yes_female + 1
    else:
        other_no_female = other_no_female + 1
    a=a+1

axs[3].axvline(2, c="black", ls="--", lw=0.2)
axs[1].axvline(2, c="black", ls="--", lw=0.2)

axs[3].grid(False)
axs[3].set_facecolor("white")
axs[3].spines['bottom'].set_color('0')
axs[3].spines['left'].set_color('0')

axs[3].tick_params(bottom='on', left='on', which="both")


axs[1].grid(False)
axs[1].set_facecolor("white")
axs[1].spines['bottom'].set_color('0')
axs[1].spines['left'].set_color('0')

axs[1].tick_params(bottom='on', left='on', which="both")





#Repeat for normals

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

orig_cols = df.columns.tolist()

df_barcode = df['Barcode'].tolist()

end_barcode = []

for a in df_barcode:
    end_barcode.append(int(split_advanced(a, "-", 3)[1]))

df['End_Val'] = end_barcode

df = df[df['End_Val']==11] #only 11s are present for normal samples
df = df[df['Classification'] != "UNKNOWN"]
df = df[df['Gender'] != "UNKNOWN"]

#No duplicates in normals

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

for a in invalid_ids:
    df = df[~(df['Barcode'].str.contains(a))]

df2 = df[df['Gender']=="FEMALE"]
df = df[df['Gender']=="MALE"]

df_class = df['Classification'].tolist()
df2_class = df['Classification'].tolist()

male_exp_unp = df['XIST_TPM'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

female_exp_unp = df2['XIST_TPM'].tolist()

female_exp = []

for a in female_exp_unp:
    if a>=0:
        female_exp.append(math.log(a+1, 2))
    else:
        female_exp.append(0)

axs[2].hist(male_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")

axs[0].hist(female_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")

tgct_yes = 0
tgct_no = 0

other_yes_male = 0
other_no_male = 0

a=0
while a<len(male_exp):
    if df_class[a] == "TGCT":
        if male_exp[a] >= 2:
            tgct_yes = tgct_yes + 1
        else:
            tgct_no = tgct_no + 1
    else:
        if male_exp[a] >= 2:
            other_yes_male = other_yes_male + 1
        else:
            other_no_male = other_no_male + 1
    a=a+1

other_yes_female = 0
other_no_female = 0

a=0
while a<len(female_exp):
    if female_exp[a] >= 2:
        other_yes_female = other_yes_female + 1
    else:
        other_no_female = other_no_female + 1
    a=a+1

plt.xlim(0, 12.5)

axs[2].axvline(2, c="black", ls="--", lw=0.2)
axs[0].axvline(2, c="black", ls="--", lw=0.2)

axs[2].grid(False)
axs[2].set_facecolor("white")
axs[2].spines['bottom'].set_color('0')
axs[2].spines['left'].set_color('0')

axs[2].tick_params(bottom='on', left='on', which="both")

axs[0].grid(False)
axs[0].set_facecolor("white")
axs[0].spines['bottom'].set_color('0')
axs[0].spines['left'].set_color('0')

axs[0].tick_params(bottom='on', left='on', which="both")

axs[0].set_yscale('log')
axs[2].set_yscale('log')
axs[1].set_yscale('log')
axs[3].set_yscale('log')

axs[3].set_ylabel("Number of\nSamples")
axs[1].set_ylabel("Number of\nSamples")
axs[2].set_ylabel("Number of\nSamples")
axs[3].set_xlabel("XIST Expression (log$_2$(TPM+1))")
axs[0].set_ylabel("Number of\nSamples")

axs[0].spines['bottom'].set_lw(0.2)
axs[0].spines['left'].set_lw(0.2)
axs[1].spines['bottom'].set_lw(0.2)
axs[1].spines['left'].set_lw(0.2)
axs[2].spines['bottom'].set_lw(0.2)
axs[2].spines['left'].set_lw(0.2)
axs[3].spines['bottom'].set_lw(0.2)
axs[3].spines['left'].set_lw(0.2)

axs[0].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[0].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[0].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")
axs[0].set_ylim([0, 22])

axs[0].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20], minor=True)

axs[2].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[2].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[2].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")
axs[2].set_ylim([0, 350])

axs[2].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300], minor=True)

axs[1].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[1].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[1].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")
axs[1].set_ylim([0, 250])

axs[1].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200], minor=True)

axs[3].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[3].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[3].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")
axs[3].set_ylim([0, 3400])

axs[3].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900, 2000, 3000], minor=True)

axs[0].xaxis.labelpad = 1
axs[0].yaxis.labelpad = 8.5
axs[1].xaxis.labelpad = 1
axs[1].yaxis.labelpad = 8.5
axs[2].xaxis.labelpad = 1
axs[2].yaxis.labelpad = 8.5
axs[3].xaxis.labelpad = 1
axs[3].yaxis.labelpad = 1

axs[0].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))

axs[1].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))

axs[2].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))

axs[3].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/OP5_hist_log_normal_and_cancer_XIST_exp_Xena_TCGA_male_and_female.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/OP5_hist_log_normal_and_cancer_XIST_exp_Xena_TCGA_male_and_female.pdf", dpi=dpi_set, bbox_inches = 'tight')




