#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 13:12:54 2022

@author: ananthansadagopan
"""

#skewness_female_tumors_normalized_E_expression_TPM_geq0.csv
#skewness_male_tumors_normalized_NE_expression_TPM_geq0.csv
#skewness_female_tumors_normalized_E_expression_TPM_geq0.csv
#skewness_male_tumors_normalized_NE_expression_TPM_geq0.csv

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

"""
NOT centered at 0 bc of log transform, take average expression in XIST-, average(a - mean) across
all a is entered at 0; however, average(log2(a)-log2(mean)) is not
"""
#Plotting stdev vs. location might be the move?

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.linewidth'] = 0.3

sns.set(rc={'figure.figsize':(6.4,5)})
dpi_set = 72 # change the output resolution
sns.set(font_scale=1)

gene_type = "ESCAPE"
plot_type = "MALE"

if gene_type == "NONESCAPE":
    m_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/skewness_male_tumors_normalized_NE_expression_TPM_geq0.csv")
    f_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/skewness_female_tumors_normalized_NE_expression_TPM_geq0.csv")
elif gene_type == "ESCAPE":
    m_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/skewness_male_tumors_normalized_E_expression_TPM_geq0.csv")
    f_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/skewness_female_tumors_normalized_E_expression_TPM_geq0.csv")

m_gene_list = m_df['Unnamed: 0'].tolist()
f_gene_list = f_df['Unnamed: 0'].tolist()

if m_gene_list != f_gene_list:
    print("ERROR gene lists not equivalent")

m_cols = m_df.columns.tolist()
f_cols = f_df.columns.tolist()

df_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

#df_ref = df_ref[df_ref['Classification'].isin(["OV"])]

df_h = df_ref[df_ref['XIST_TPM']>=3]
df_l = df_ref[df_ref['XIST_TPM']<3]

df_h_samples = df_h['Barcode'].tolist()
df_l_samples = df_l['Barcode'].tolist()

ref1 = list(set(df_h_samples) & set(m_cols))
ref2 = list(set(df_l_samples) & set(m_cols))
ref3 = list(set(df_h_samples) & set(f_cols))
ref4 = list(set(df_l_samples) & set(f_cols))

hm_df = m_df[ref1]
lm_df = m_df[ref2]
hf_df = f_df[ref3]
lf_df = f_df[ref4]


numeric_cols = [col for col in hm_df if hm_df[col].dtype.kind != 'O']
hm_df[numeric_cols] = hm_df[numeric_cols].apply(lambda x: 2**x)

numeric_cols = [col for col in lm_df if lm_df[col].dtype.kind != 'O']
lm_df[numeric_cols] = lm_df[numeric_cols].apply(lambda x: 2**x)

numeric_cols = [col for col in hf_df if hf_df[col].dtype.kind != 'O']
hf_df[numeric_cols] = hf_df[numeric_cols].apply(lambda x: 2**x)

numeric_cols = [col for col in lf_df if lf_df[col].dtype.kind != 'O']
lf_df[numeric_cols] = lf_df[numeric_cols].apply(lambda x: 2**x)




df_RNA_temp = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/gene_id_name_chr_chrX_biomart.csv")

location = df_RNA_temp['Gene start (bp hg19)'].tolist()
gene_symbol = df_RNA_temp['Gene name'].tolist()

loc_dict = dict(zip(gene_symbol,location))

loc_vals = []

for a in m_gene_list:
    loc_vals.append(loc_dict[a])

"""
hm_df_mean = list(hm_df.median(axis = 1))
lm_df_mean = list(lm_df.median(axis = 1))
hf_df_mean = list(hf_df.median(axis = 1))
lf_df_mean = list(lf_df.median(axis = 1))
"""

hm_df_mean = list(hm_df.mean(axis = 1))
lm_df_mean = list(lm_df.mean(axis = 1))
hf_df_mean = list(hf_df.mean(axis = 1))
lf_df_mean = list(lf_df.mean(axis = 1))




fig = plt.figure()

gene_list2 = []

q=0
while q<len(lm_df_mean):
    if hm_df_mean[q] <= 2 and hm_df_mean[q] >= 0.5:
        gene_list2.append(m_gene_list[q])
    q=q+1


males_invalid = np.setdiff1d(m_gene_list, gene_list2)
print(len(males_invalid))
print(len(m_gene_list))

gene_list3 = []

q=0
while q<len(hf_df_mean):
    if lf_df_mean[q] <= 2 and lf_df_mean[q] >= 0.5:
        gene_list3.append(m_gene_list[q])
    q=q+1

females_invalid = np.setdiff1d(f_gene_list, gene_list3)
print(len(females_invalid))
print(len(f_gene_list))

valid_genes_df = pd.DataFrame([gene_list2, gene_list3]).T
valid_genes_df.columns = ['XISTneg_male_ref', 'XISTpos_female_ref']

if gene_type == "NONESCAPE":
    valid_genes_df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/ne_valid_genes_0.5_to_2_cutoff.csv", index=False)
elif gene_type == "ESCAPE":
    valid_genes_df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/e_valid_genes_0.5_to_2_cutoff.csvv", index=False)


if plot_type == "FEMALE":
    plt.scatter(loc_vals,hf_df_mean, c="black", alpha=1, s=1)
    plt.scatter(loc_vals,lf_df_mean, c="red", alpha=1, s=1, zorder=10)
elif plot_type == "MALE":
    plt.scatter(loc_vals,hm_df_mean, c="red", alpha=1, s=1)
    plt.scatter(loc_vals,lm_df_mean, c="black", alpha=1, s=1, zorder=10)
    
ax = plt.gca()
#ax.set_ylim([0, 2])


"""
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
ax.set_aspect('equal')

ax.set_xlim([80, 805000])
ax.set_ylim([80, 805000])

ax.set_yscale('log')
ax.set_xscale('log')
"""

#ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
#ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

#plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


patch_list = []


if plot_type == "FEMALE":
    patch_list.append(mpatches.Patch(color='red', label="XIST-\nFemale"))
    patch_list.append(mpatches.Patch(color='black', label="XIST+\nFemale"))
elif plot_type == "MALE":
    patch_list.append(mpatches.Patch(color='red', label="XIST+\nMale"))
    patch_list.append(mpatches.Patch(color='black', label="XIST-\nMale"))

ax.legend(handles=patch_list, facecolor='white', loc="upper right", fontsize=10, bbox_to_anchor=[1.25, 0.97])

ax.set_facecolor("white")
ax.grid(False)
ax.tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
ax.tick_params(axis='y', which='major', left=True, labelsize=10, size=4)

ax.spines['bottom'].set_color('0')
ax.spines['left'].set_color('0')

ax.set_yscale('log', basey=2)

plt.xlabel("hg38 Genomic Location on chrX")
    
plt.ylabel("Average Normalized Gene Expression")

plt.axhline(0.5, c="black", ls="--", lw=0.4)
plt.axhline(2, c="black", ls="--", lw=0.4)

fig.tight_layout()

dpi_set = 72
plt.tick_params(bottom='on', left='on')

if gene_type == "NONESCAPE" and plot_type == "FEMALE":
    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/normalized_ne_expression_vs_chrX_location_female.pdf", dpi=dpi_set, bbox_inches = 'tight')
elif gene_type == "ESCAPE" and plot_type == "FEMALE":
    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/normalized_e_expression_vs_chrX_location_female.pdf", dpi=dpi_set, bbox_inches = 'tight')
elif gene_type == "NONESCAPE" and plot_type == "MALE":
    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/normalized_ne_expression_vs_chrX_location_male.pdf", dpi=dpi_set, bbox_inches = 'tight')
elif gene_type == "ESCAPE" and plot_type == "MALE":
    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/normalized_e_expression_vs_chrX_location_male.pdf", dpi=dpi_set, bbox_inches = 'tight')






