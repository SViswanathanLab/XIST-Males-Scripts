#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 07:58:26 2021

@author: ananthansadagopan
"""

import pandas as pd 
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors
import matplotlib.patches as mpatches
import collections
import scipy
from matplotlib.ticker import StrMethodFormatter, NullFormatter, ScalarFormatter, FormatStrFormatter
from matplotlib.ticker import PercentFormatter

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(3.61,4.41)})
sns.set(font_scale=0.6)
dpi_set = 72

bins_list = np.arange(0,3.01,0.03)

final_name = "XISTpos_females"

if final_name == "XISTpos_males" or final_name == "XISTneg_males":
    df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_male_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2_revised.csv")    
elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
    df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_female_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2_revised.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

for a in invalid_ids:
    df = df[~(df['Truncated_Barcodes'].str.contains(a))]


df['total_RNA_counts'] = df['RNA_alt_counts'] + df['RNA_ref_counts']
df['total_DNA_counts'] = df['DNA_alt_counts'] + df['DNA_ref_counts']
df['DNA_VAF_purity_ratio'] = df['DNA_VAF']/df['Purity']

cutoff = 20

df = df[df['total_RNA_counts']>=cutoff]
df = df[df['total_DNA_counts']>=cutoff]

df = df[df['DNA_VAF_purity_ratio']>=0.4]

"""
df = df[df['RNA_alt_counts']>=cutoff]
df = df[df['RNA_ref_counts']>=cutoff]
df = df[df['DNA_alt_counts']>=cutoff]
df = df[df['DNA_ref_counts']>=cutoff]
"""

if final_name == "XISTpos_males" or final_name == "XISTpos_females":
    df = df[df['XIST_positive_1_is_yes']==1]
elif final_name == "XISTneg_females" or final_name == "XISTneg_males":
    df = df[df['XIST_positive_1_is_yes']==0]

df = df[df['Variant_Classification'].isin(['Missense_Mutation', 'Silent'])]
#df = df[df['Variant_Classification'].isin(['Silent'])]

xs = df['DNA_VAF'].tolist()
ys = df['RNA_VAF'].tolist()

if len(xs) != len(ys):
    print("ERROR")

ratio_list = []

b=0
while b<len(xs):
    try:
        ratio_list.append(ys[b]/xs[b])
    except ZeroDivisionError:
        b=b+1
        continue
    b=b+1

n_bins = 100

fig, axs = plt.subplots(4, 1, sharex=True, tight_layout=True)

"""
x_unp, y_unp = sns.distplot(ratio_list, hist=False, kde=True, ax=axs[0],
             bins=n_bins, color = 'black',
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

print(len(ratio_list))

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 4:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1

axs[0].fill_between(x, 0,  y, facecolor='black', alpha=0.2)
"""

axs[0].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="black", linewidth=0.5, edgecolor="none")


XISTpos_female_ratio = ratio_list

if final_name == "XISTneg_males":
    n_bins = 100    
elif final_name == "XISTpos_females":
    n_bins = 150    
elif final_name == "XISTpos_males":
    n_bins = 60
elif final_name == "XISTneg_females":
    n_bins = 60

#ax.hist(ratio_list, density=False, bins=n_bins, color="black", linewidth=0.5, edgecolor="none")

plt.xlim([0, 3])

axs[0].grid(False)
axs[0].set_facecolor("white")
axs[0].spines['bottom'].set_color('0')
axs[0].spines['left'].set_color('0')

axs[0].tick_params(bottom='on', left='on')

plt.xlabel("RNA VAF / DNA VAF")

plt.ylabel("Number of Variants")









final_name = "XISTneg_females"

if final_name == "XISTpos_males" or final_name == "XISTneg_males":
    df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_male_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2_revised.csv")
elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
    df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_female_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2_revised.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

for a in invalid_ids:
    df = df[~(df['Truncated_Barcodes'].str.contains(a))]

df['total_RNA_counts'] = df['RNA_alt_counts'] + df['RNA_ref_counts']
df['total_DNA_counts'] = df['DNA_alt_counts'] + df['DNA_ref_counts']
df['DNA_VAF_purity_ratio'] = df['DNA_VAF']/df['Purity']

cutoff = 20

df = df[df['total_RNA_counts']>=cutoff]
df = df[df['total_DNA_counts']>=cutoff]

df = df[df['DNA_VAF_purity_ratio']>=0.4]

"""
df = df[df['RNA_alt_counts']>=cutoff]
df = df[df['RNA_ref_counts']>=cutoff]
df = df[df['DNA_alt_counts']>=cutoff]
df = df[df['DNA_ref_counts']>=cutoff]
"""

if final_name == "XISTpos_males" or final_name == "XISTpos_females":
    df = df[df['XIST_positive_1_is_yes']==1]
elif final_name == "XISTneg_females" or final_name == "XISTneg_males":
    df = df[df['XIST_positive_1_is_yes']==0]

df = df[df['Variant_Classification'].isin(['Missense_Mutation', 'Silent'])]
#df = df[df['Variant_Classification'].isin(['Silent'])]

xs = df['DNA_VAF'].tolist()
ys = df['RNA_VAF'].tolist()

if len(xs) != len(ys):
    print("ERROR")

ratio_list = []

b=0
while b<len(xs):
    try:
        ratio_list.append(ys[b]/xs[b])
    except ZeroDivisionError:
        b=b+1
        continue
    b=b+1

n_bins = 100

"""
x_unp, y_unp = sns.distplot(ratio_list, hist=False, kde=True, ax=axs[1],
             bins=n_bins, color = 'orange',
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()
    
print(len(ratio_list))

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 4:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[1].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)
"""
axs[1].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="orange", linewidth=0.5, edgecolor="none")

XISTneg_female_ratio = ratio_list

if final_name == "XISTneg_males":
    n_bins = 100    
elif final_name == "XISTpos_females":
    n_bins = 150    
elif final_name == "XISTpos_males":
    n_bins = 60
elif final_name == "XISTneg_females":
    n_bins = 60

#ax.hist(ratio_list, density=False, bins=n_bins, color="black", linewidth=0.5, edgecolor="none")

plt.xlim([0, 3])

axs[1].grid(False)
axs[1].set_facecolor("white")
axs[1].spines['bottom'].set_color('0')
axs[1].spines['left'].set_color('0')

axs[1].tick_params(bottom='on', left='on')

plt.xlabel("RNA VAF / DNA VAF")

plt.ylabel("Number of Variants")









plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(5,5)})
sns.set(font_scale=1)
dpi_set = 72

final_name = "XISTpos_males"

if final_name == "XISTpos_males" or final_name == "XISTneg_males":
    df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_male_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2_revised.csv")
elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
    df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_female_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2_revised.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

for a in invalid_ids:
    df = df[~(df['Truncated_Barcodes'].str.contains(a))]

df['total_RNA_counts'] = df['RNA_alt_counts'] + df['RNA_ref_counts']
df['total_DNA_counts'] = df['DNA_alt_counts'] + df['DNA_ref_counts']
df['DNA_VAF_purity_ratio'] = df['DNA_VAF']/df['Purity']

cutoff = 20

df = df[df['total_RNA_counts']>=cutoff]
df = df[df['total_DNA_counts']>=cutoff]

df = df[df['DNA_VAF_purity_ratio']>=0.4]

"""
df = df[df['RNA_alt_counts']>=cutoff]
df = df[df['RNA_ref_counts']>=cutoff]
df = df[df['DNA_alt_counts']>=cutoff]
df = df[df['DNA_ref_counts']>=cutoff]
"""

if final_name == "XISTpos_males" or final_name == "XISTpos_females":
    df = df[df['XIST_positive_1_is_yes']==1]
elif final_name == "XISTneg_females" or final_name == "XISTneg_males":
    df = df[df['XIST_positive_1_is_yes']==0]

df = df[df['Variant_Classification'].isin(['Missense_Mutation', 'Silent'])]
#df = df[df['Variant_Classification'].isin(['Silent'])]

xs = df['DNA_VAF'].tolist()
ys = df['RNA_VAF'].tolist()

if len(xs) != len(ys):
    print("ERROR")

ratio_list = []

b=0
while b<len(xs):
    try:
        ratio_list.append(ys[b]/xs[b])
    except ZeroDivisionError:
        b=b+1
        continue
    b=b+1

n_bins = 100
"""
x_unp, y_unp = sns.distplot(ratio_list, hist=False, kde=True, ax=axs[2],
             bins=n_bins, color = 'red',
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

print(len(ratio_list))

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 4:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1

axs[2].fill_between(x, 0,  y, facecolor='red', alpha=0.2)
"""
axs[2].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="red", linewidth=0.5, edgecolor="none")

XISTpos_male_ratio = ratio_list

if final_name == "XISTneg_males":
    n_bins = 100    
elif final_name == "XISTpos_females":
    n_bins = 150    
elif final_name == "XISTpos_males":
    n_bins = 60
elif final_name == "XISTneg_females":
    n_bins = 60

#ax.hist(ratio_list, density=False, bins=n_bins, color="black", linewidth=0.5, edgecolor="none")

plt.xlim([0, 3])

axs[2].grid(False)
axs[2].set_facecolor("white")
axs[2].spines['bottom'].set_color('0')
axs[2].spines['left'].set_color('0')

axs[2].tick_params(bottom='on', left='on')

plt.xlabel("RNA VAF / DNA VAF")

plt.ylabel("Number of Variants")










plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(5,5)})
sns.set(font_scale=1)
dpi_set = 72

final_name = "XISTneg_males"

if final_name == "XISTpos_males" or final_name == "XISTneg_males":
    df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_male_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2_revised.csv")
elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
    df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_female_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2_revised.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

for a in invalid_ids:
    df = df[~(df['Truncated_Barcodes'].str.contains(a))]

df['total_RNA_counts'] = df['RNA_alt_counts'] + df['RNA_ref_counts']
df['total_DNA_counts'] = df['DNA_alt_counts'] + df['DNA_ref_counts']
df['DNA_VAF_purity_ratio'] = df['DNA_VAF']/df['Purity']

cutoff = 20

df = df[df['total_RNA_counts']>=cutoff]
df = df[df['total_DNA_counts']>=cutoff]

df = df[df['DNA_VAF_purity_ratio']>=0.4]

"""
df = df[df['RNA_alt_counts']>=cutoff]
df = df[df['RNA_ref_counts']>=cutoff]
df = df[df['DNA_alt_counts']>=cutoff]
df = df[df['DNA_ref_counts']>=cutoff]
"""

if final_name == "XISTpos_males" or final_name == "XISTpos_females":
    df = df[df['XIST_positive_1_is_yes']==1]
elif final_name == "XISTneg_females" or final_name == "XISTneg_males":
    df = df[df['XIST_positive_1_is_yes']==0]

df = df[df['Variant_Classification'].isin(['Missense_Mutation', 'Silent'])]
#df = df[df['Variant_Classification'].isin(['Silent'])]

xs = df['DNA_VAF'].tolist()
ys = df['RNA_VAF'].tolist()

if len(xs) != len(ys):
    print("ERROR")

ratio_list = []

b=0
while b<len(xs):
    try:
        ratio_list.append(ys[b]/xs[b])
    except ZeroDivisionError:
        b=b+1
        continue
    b=b+1

n_bins = 100
"""
x_unp, y_unp = sns.distplot(ratio_list, hist=False, kde=True, ax=axs[3],
             bins=n_bins, color = 'blue',
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

print(len(ratio_list))

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 4:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1

axs[3].fill_between(x, 0,  y, facecolor='blue', alpha=0.2)
"""
axs[3].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="blue", linewidth=0.5, edgecolor="none")

XISTneg_male_ratio = ratio_list

if final_name == "XISTneg_males":
    n_bins = 100    
elif final_name == "XISTpos_females":
    n_bins = 150    
elif final_name == "XISTpos_males":
    n_bins = 60
elif final_name == "XISTneg_females":
    n_bins = 60

#ax.hist(ratio_list, density=False, bins=n_bins, color="black", linewidth=0.5, edgecolor="none")

plt.xlim([0, 3])

axs[3].grid(False)
axs[3].set_facecolor("white")
axs[3].spines['bottom'].set_color('0')
axs[3].spines['left'].set_color('0')

axs[3].tick_params(bottom='on', left='on')

plt.xlabel("Ratio of VAF for chrX Somatic Variants (RNA-Seq / WES)")

xist_pos_female_patch = mpatches.Patch(color='black', label='XIST+ Females')
xist_neg_female_patch = mpatches.Patch(color='orange', label='XIST- Females')
xist_pos_male_patch = mpatches.Patch(color='red', label='XIST+ Males')
xist_neg_male_patch = mpatches.Patch(color='blue', label='XIST- Males')


"""
leg = ax.legend(handles=[xist_pos_female_patch, xist_neg_female_patch], facecolor='white', bbox_to_anchor=(1, 1), fontsize=7)
    leg = ax.legend(handles=[xist_pos_male_patch, xist_neg_male_patch], facecolor='white', bbox_to_anchor=(1, 1), fontsize=7)

frame = leg.get_frame()
frame.set_edgecolor("black")
frame.set_linewidth(0.5)
"""

q=0
while q<4:    
    axs[q].xaxis.labelpad = 4
    axs[q].yaxis.labelpad = 4
    axs[q].spines['bottom'].set_lw(0.2)
    axs[q].spines['left'].set_lw(0.2)
    axs[q].xaxis.set_tick_params(width=0.2, size=3, pad=2)
    axs[q].yaxis.set_tick_params(width=0.2, size=3, pad=2)
    axs[q].set_ylabel("Percentage\nof Variants")
    
    axs[q].set_yticks([0, 0.2, 0.4])    
    axs[q].set_xticks([0, 1, 2, 3])
    #axs[q].yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
    axs[q].yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    
    """
    if q == 1:
        axs[q].yaxis.labelpad = 12
    if q == 2:
        axs[q].yaxis.labelpad = 8.2
    """        
    #plt.ylim([0, 1.02])
    q=q+1
    
    
    
    

print("FEMALE - KS TEST")

print("Autosome distribution")

a, b = scipy.stats.ks_2samp(XISTpos_female_ratio, XISTneg_female_ratio)

print(a)
print(b)
    
    
print("MALE - KS TEST")

print("Autosome distribution")

a, b = scipy.stats.ks_2samp(XISTpos_male_ratio, XISTneg_male_ratio)

print(a)
print(b)
    

"""

counter=0

for a in ratio_list:
    if a<0.5:
        counter=counter+1

leq_point_5_ratio = str("{:.1f}".format(counter/len(ratio_list)*100) + "% of Variants with RNA VAF / DNA VAF < 0.5")

ax.text(0.07, 1, leq_point_5_ratio, transform=ax.transAxes,
     fontsize=12, verticalalignment='top', color="black")
"""

plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_and_FEMALES_VAF_density_plots_test.png", dpi=dpi_set)
plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_and_FEMALES_VAF_density_plots_test.pdf", dpi=dpi_set)
