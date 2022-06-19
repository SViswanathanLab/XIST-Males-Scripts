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
from statsmodels.stats.proportion import proportions_ztest
import matplotlib.gridspec as gridspec
import os
import random

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(3.61*4.41/2.396,4.41)})
sns.set(font_scale=0.6)
dpi_set = 72

DNA_CCF_cutoff = 1

#gene_class = "Non-Escapee"
gene_class = "Escapee"
cn_to_analyze = 2

n_cutoff = 0.15
n_iterations = 100000

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/merged_TITAN_best_solutions.txt", sep="\t")
df = df[df['Sample']!="Sample"]

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

for a in invalid_ids:
    df = df[~(df['Sample'].str.contains(a))]

df['Start'] = df['Start'].astype("float64")
df['End'] = df['End'].astype("float64")
df['Copy_Number'] = df['Copy_Number'].astype("float64")

df_seg = df
df_seg = df_seg[df_seg['Chromosome']=="X"]

cn_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/TITAN_final_solutions_plus_annotations.csv")

purity_dict = dict(zip(cn_df['Sub_ID'],cn_df['TITAN_purity']))

cn2_df = cn_df[cn_df['Good?']!="No"]
cn2_df = cn_df[cn_df['Rounded_chrX_CN']!="Unknown"]
cn2_df['Rounded_chrX_CN'] = cn2_df['Rounded_chrX_CN'].astype("int")
if cn_to_analyze == 2:
    cn2_df = cn2_df[cn2_df['Rounded_chrX_CN']==2]
elif cn_to_analyze == 1:
    cn2_df = cn2_df[cn2_df['Rounded_chrX_CN']==1]
unp_samples = cn2_df['ID'].tolist()

unp_samples = cn2_df['ID'].tolist()

cn_select_samples = []

for a in unp_samples:
    cn_select_samples.append(a.split("_")[0])

PAR_genes = ['PLCXD1', 'GTPBP6', 'PPP2R3B', 'SHOX', 'CRLF2', 'CSF2RA', 'IL3RA', 'SLC25A6', 'ASMTL', 'P2RY8', 'CXYorf3', 'AKAP17A', 'ASMT', 'DHRSXY', 'ZBED1', 'CD99', 'XG']
df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")
PAR_genes = [x for x in df_gene_class['Escaping_with_PARs'].tolist() if x == x]
non_escapee = [x for x in df_gene_class['Inactivated'].tolist() if x == x]

if gene_class == "Escapee":
    genes_of_interest = PAR_genes
elif gene_class == "Non-Escapee":
    genes_of_interest = non_escapee

bins_list = np.arange(0,3.02,0.02)

final_name = "XISTpos_females"
fname = "/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/out_" + final_name + "_" + gene_class + "_chrX_CN" + str(cn_to_analyze) + ".csv"

if os.path.isfile(fname):
    df = pd.read_csv(fname)
else:
    if final_name == "XISTpos_males" or final_name == "XISTneg_males":
        df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_chrX_muts_REVISED_males_unaveraged_TPM_geq2.csv")    
    elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
        df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_chrX_muts_REVISED_females_unaveraged_TPM_geq2.csv")
    invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

    df = df[(df['Hugo_Symbol'].isin(genes_of_interest))]   
    df.drop_duplicates(subset=['Hugo_Symbol','Start_Position','Tumor_Sample_Barcode'], keep="first", inplace=True)

    for a in invalid_ids:
        df = df[~(df['Truncated_Barcodes'].str.contains(a))]

    secondary_barcodes = []

    for a in df['Truncated_Barcodes'].tolist():
        secondary_barcodes.append(split_advanced(a, "-", 3)[0])

    df['Secondary_Barcodes'] = secondary_barcodes

    df = df[df['Truncated_Barcodes'].isin(cn_select_samples)]


    df['total_RNA_counts'] = df['RNA_alt_counts'] + df['RNA_ref_counts']
    df['total_DNA_counts'] = df['DNA_alt_counts'] + df['DNA_ref_counts']
    df['DNA_VAF_purity_ratio'] = df['DNA_VAF']/df['Purity']

    cutoff = 20

    df = df[df['total_RNA_counts']>=cutoff]
    df = df[df['total_DNA_counts']>=cutoff]

    #df = df[df['DNA_VAF_purity_ratio']>=0.4]

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

    ids = df['Truncated_Barcodes'].tolist()
    start = df['Start_Position'].tolist()

    cn_temp = []

    print(len(ids))

    new_purity = []

    a=0
    while a<len(ids):
        temp_id = split_advanced(ids[a], "-", 3)[0]
        new_purity.append(purity_dict[temp_id])
        temp_df = df_seg[df_seg['Sample'].str.contains(temp_id)]
        start_vals = temp_df['Start'].tolist()
        cn_vals = temp_df['Copy_Number'].tolist()
        q=0
        while (q-1)<len(start_vals):
            if q == len(start_vals):
                cn_temp.append(cn_vals[q-1])
                break
            elif start[a] > start_vals[q]:
                q=q+1
            else:
                cn_temp.append(cn_vals[q-1])
                break
        a=a+1
        if a % 100 == 0:
            print(a)

    df['Purity'] = new_purity
    df['CN_temp'] = cn_temp

    #CCF calculations: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5538405/

    DNA_VAF = df['DNA_VAF'].tolist()
    RNA_VAF = df['RNA_VAF'].tolist()
    purity = df['Purity'].tolist()

    DNA_CCF = []
    RNA_CCF = []
    DNA_ui = []
    RNA_ui = []
    DNA_multiplicity = []
    RNA_multiplicity = []

    if final_name == "XISTpos_males" or final_name == "XISTneg_males":
        normal_cn = 1
    elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
        normal_cn = 2

    a=0
    while a<len(cn_temp):
        temp_DNA_ui = DNA_VAF[a]*1/purity[a]*(purity[a]*cn_temp[a]+(1-purity[a])*normal_cn)
        temp_RNA_ui = RNA_VAF[a]*1/purity[a]*(purity[a]*cn_temp[a]+(1-purity[a])*normal_cn)

        DNA_ui.append(temp_DNA_ui)
        RNA_ui.append(temp_RNA_ui)    

        if temp_DNA_ui >= 1:
            DNA_multiplicity.append(temp_DNA_ui)
            DNA_CCF.append(1)
        else:
            DNA_multiplicity.append(1)
            DNA_CCF.append(temp_DNA_ui)

        if temp_RNA_ui >= 1:
            RNA_multiplicity.append(temp_RNA_ui)
            RNA_CCF.append(1)
        else:
            RNA_multiplicity.append(1)
            RNA_CCF.append(temp_RNA_ui)

        a=a+1

    df['DNA_ui'] = DNA_ui
    df['RNA_ui'] = RNA_ui
    df['DNA_multiplicity'] = DNA_multiplicity
    df['RNA_multiplicity'] = RNA_multiplicity
    df['DNA_CCF'] = DNA_CCF
    df['RNA_CCF'] = RNA_CCF

    df = df[df['DNA_CCF']>=DNA_CCF_cutoff]

    df.to_csv(fname, index=False)

#df = df[df['Variant_Classification'].isin(['Silent'])]

df = df[df['Hugo_Symbol']!="MXRA5"]

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

fig, axs = plt.subplots(4, 1)

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

#axs[0].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="black", linewidth=0.5, edgecolor="none")
#axs[1].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="black", linewidth=0.5, edgecolor="none")

sns.histplot(x=ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), kde=False, bins=bins_list, ax=axs[0], color="black", linewidth=0.5, edgecolor="none", line_kws={'lw':0.4}, kde_kws={'bw_adjust':0.23})
#sns.histplot(x=ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), kde=False, bins=bins_list, ax=axs[1], color="black", linewidth=0.5, edgecolor="none", line_kws={'lw':0.4}, kde_kws={'bw_adjust':0.23})

#sns.distplot(..., hist_kws={'weights': your weights array}, ...)

# zoom-in / limit the view to different portions of the data

# hide the spines between ax and ax2
#axs[0].spines['bottom'].set_visible(False)
#axs[1].spines['top'].set_visible(False)
#axs[0].xaxis.tick_top()
#axs[0].tick_params(labeltop=False)  # don't put tick labels at the top
#axs[1].xaxis.tick_bottom()

#plt.xlim([0, 3])

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

plt.xlabel("RNA VAF / DNA VAF")

plt.ylabel("Number of Variants")









final_name = "XISTneg_females"
fname = "/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/out_" + final_name + "_" + gene_class + "_chrX_CN" + str(cn_to_analyze) + ".csv"

if os.path.isfile(fname):
    df = pd.read_csv(fname)
else:
    if final_name == "XISTpos_males" or final_name == "XISTneg_males":
        df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_chrX_muts_REVISED_males_unaveraged_TPM_geq2.csv")    
    elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
        df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_chrX_muts_REVISED_females_unaveraged_TPM_geq2.csv")
        
    invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

    df = df[(df['Hugo_Symbol'].isin(genes_of_interest))]
    df.drop_duplicates(subset=['Hugo_Symbol','Start_Position','Tumor_Sample_Barcode'], keep="first", inplace=True)

    for a in invalid_ids:
        df = df[~(df['Truncated_Barcodes'].str.contains(a))]

    secondary_barcodes = []

    for a in df['Truncated_Barcodes'].tolist():
        secondary_barcodes.append(split_advanced(a, "-", 3)[0])

    df['Secondary_Barcodes'] = secondary_barcodes

    df = df[df['Truncated_Barcodes'].isin(cn_select_samples)]

    df['total_RNA_counts'] = df['RNA_alt_counts'] + df['RNA_ref_counts']
    df['total_DNA_counts'] = df['DNA_alt_counts'] + df['DNA_ref_counts']
    df['DNA_VAF_purity_ratio'] = df['DNA_VAF']/df['Purity']

    cutoff = 20

    df = df[df['total_RNA_counts']>=cutoff]
    df = df[df['total_DNA_counts']>=cutoff]

    #df = df[df['DNA_VAF_purity_ratio']>=0.4]

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

    ids = df['Truncated_Barcodes'].tolist()
    start = df['Start_Position'].tolist()

    cn_temp = []

    new_purity = []

    a=0
    while a<len(ids):
        temp_id = split_advanced(ids[a], "-", 3)[0]
        new_purity.append(purity_dict[temp_id])
        temp_df = df_seg[df_seg['Sample'].str.contains(temp_id)]
        start_vals = temp_df['Start'].tolist()
        cn_vals = temp_df['Copy_Number'].tolist()
        q=0
        while (q-1)<len(start_vals):
            if q == len(start_vals):
                cn_temp.append(cn_vals[q-1])
                break
            elif start[a] > start_vals[q]:
                q=q+1
            else:
                cn_temp.append(cn_vals[q-1])
                break
        a=a+1
        if a % 100 == 0:
            print(a)

    df['Purity'] = new_purity
    df['CN_temp'] = cn_temp

    #CCF calculations: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5538405/

    DNA_VAF = df['DNA_VAF'].tolist()
    RNA_VAF = df['RNA_VAF'].tolist()
    purity = df['Purity'].tolist()

    DNA_CCF = []
    RNA_CCF = []
    DNA_ui = []
    RNA_ui = []
    DNA_multiplicity = []
    RNA_multiplicity = []

    if final_name == "XISTpos_males" or final_name == "XISTneg_males":
        normal_cn = 1
    elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
        normal_cn = 2

    a=0
    while a<len(cn_temp):
        temp_DNA_ui = DNA_VAF[a]*1/purity[a]*(purity[a]*cn_temp[a]+(1-purity[a])*normal_cn)
        temp_RNA_ui = RNA_VAF[a]*1/purity[a]*(purity[a]*cn_temp[a]+(1-purity[a])*normal_cn)

        DNA_ui.append(temp_DNA_ui)
        RNA_ui.append(temp_RNA_ui)    

        if temp_DNA_ui >= 1:
            DNA_multiplicity.append(temp_DNA_ui)
            DNA_CCF.append(1)
        else:
            DNA_multiplicity.append(1)
            DNA_CCF.append(temp_DNA_ui)

        if temp_RNA_ui >= 1:
            RNA_multiplicity.append(temp_RNA_ui)
            RNA_CCF.append(1)
        else:
            RNA_multiplicity.append(1)
            RNA_CCF.append(temp_RNA_ui)

        a=a+1

    df['DNA_ui'] = DNA_ui
    df['RNA_ui'] = RNA_ui
    df['DNA_multiplicity'] = DNA_multiplicity
    df['RNA_multiplicity'] = RNA_multiplicity
    df['DNA_CCF'] = DNA_CCF
    df['RNA_CCF'] = RNA_CCF

    df = df[df['DNA_CCF']>=DNA_CCF_cutoff]

    df.to_csv(fname, index=False)

#df = df[df['Variant_Classification'].isin(['Silent'])]
df = df[df['Hugo_Symbol']!="MXRA5"]

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
#axs[2].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="orange", linewidth=0.5, edgecolor="none")
#axs[3].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="orange", linewidth=0.5, edgecolor="none")

sns.histplot(x=ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), kde=False, bins=bins_list, ax=axs[1], color="orange", linewidth=0.5, edgecolor="none", line_kws={'lw':0.4}, kde_kws={'bw_adjust':0.15})
#sns.histplot(x=ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), kde=False, bins=bins_list, ax=axs[3], color="orange", linewidth=0.5, edgecolor="none", line_kws={'lw':0.4}, kde_kws={'bw_adjust':0.15})

XISTneg_female_ratio = ratio_list

if final_name == "XISTneg_males":
    n_bins = 100    
elif final_name == "XISTpos_females":
    n_bins = 150    
elif final_name == "XISTpos_males":
    n_bins = 60
elif final_name == "XISTneg_females":
    n_bins = 60

"""
#ax.hist(ratio_list, density=False, bins=n_bins, color="black", linewidth=0.5, edgecolor="none")
axs[2].spines['bottom'].set_visible(False)
axs[3].spines['top'].set_visible(False)
axs[2].xaxis.tick_top()
axs[2].tick_params(labeltop=False)  # don't put tick labels at the top
axs[3].xaxis.tick_bottom()


plt.xlim([0, 3])
"""

plt.xlabel("RNA VAF / DNA VAF")

plt.ylabel("Number of Variants")









plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(5,5)})
sns.set(font_scale=1)
dpi_set = 72

final_name = "XISTpos_males"
fname = "/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/out_" + final_name + "_" + gene_class + "_chrX_CN" + str(cn_to_analyze) + ".csv"

if os.path.isfile(fname):
    df = pd.read_csv(fname)
else:
    if final_name == "XISTpos_males" or final_name == "XISTneg_males":
        df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_chrX_muts_REVISED_males_unaveraged_TPM_geq2.csv")    
    elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
        df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_chrX_muts_REVISED_females_unaveraged_TPM_geq2.csv")
    
    df = df[(df['Hugo_Symbol'].isin(genes_of_interest))] 
    df.drop_duplicates(subset=['Hugo_Symbol','Start_Position','Tumor_Sample_Barcode'], keep="first", inplace=True)
    
    invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

    for a in invalid_ids:
        df = df[~(df['Truncated_Barcodes'].str.contains(a))]

    secondary_barcodes = []

    for a in df['Truncated_Barcodes'].tolist():
        secondary_barcodes.append(split_advanced(a, "-", 3)[0])

    df['Secondary_Barcodes'] = secondary_barcodes

    df = df[df['Truncated_Barcodes'].isin(cn_select_samples)]

    df['total_RNA_counts'] = df['RNA_alt_counts'] + df['RNA_ref_counts']
    df['total_DNA_counts'] = df['DNA_alt_counts'] + df['DNA_ref_counts']
    df['DNA_VAF_purity_ratio'] = df['DNA_VAF']/df['Purity']

    cutoff = 20

    df = df[df['total_RNA_counts']>=cutoff]
    df = df[df['total_DNA_counts']>=cutoff]

    #df = df[df['DNA_VAF_purity_ratio']>=0.4]

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

    ids = df['Truncated_Barcodes'].tolist()
    start = df['Start_Position'].tolist()

    cn_temp = []

    print(len(ids))

    new_purity = []

    a=0
    while a<len(ids):
        temp_id = split_advanced(ids[a], "-", 3)[0]
        new_purity.append(purity_dict[temp_id])
        temp_df = df_seg[df_seg['Sample'].str.contains(temp_id)]
        start_vals = temp_df['Start'].tolist()
        cn_vals = temp_df['Copy_Number'].tolist()
        q=0
        while (q-1)<len(start_vals):
            if q == len(start_vals):
                cn_temp.append(cn_vals[q-1])
                break
            elif start[a] > start_vals[q]:
                q=q+1
            else:
                cn_temp.append(cn_vals[q-1])
                break
        a=a+1
        if a % 100 == 0:
            print(a)

    df['Purity'] = new_purity
    df['CN_temp'] = cn_temp


    #CCF calculations: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5538405/

    DNA_VAF = df['DNA_VAF'].tolist()
    RNA_VAF = df['RNA_VAF'].tolist()
    purity = df['Purity'].tolist()

    DNA_CCF = []
    RNA_CCF = []
    DNA_ui = []
    RNA_ui = []
    DNA_multiplicity = []
    RNA_multiplicity = []

    if final_name == "XISTpos_males" or final_name == "XISTneg_males":
        normal_cn = 1
    elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
        normal_cn = 2

    a=0
    while a<len(cn_temp):
        temp_DNA_ui = DNA_VAF[a]*1/purity[a]*(purity[a]*cn_temp[a]+(1-purity[a])*normal_cn)
        temp_RNA_ui = RNA_VAF[a]*1/purity[a]*(purity[a]*cn_temp[a]+(1-purity[a])*normal_cn)

        DNA_ui.append(temp_DNA_ui)
        RNA_ui.append(temp_RNA_ui)    

        if temp_DNA_ui >= 1:
            DNA_multiplicity.append(temp_DNA_ui)
            DNA_CCF.append(1)
        else:
            DNA_multiplicity.append(1)
            DNA_CCF.append(temp_DNA_ui)

        if temp_RNA_ui >= 1:
            RNA_multiplicity.append(temp_RNA_ui)
            RNA_CCF.append(1)
        else:
            RNA_multiplicity.append(1)
            RNA_CCF.append(temp_RNA_ui)

        a=a+1

    df['DNA_ui'] = DNA_ui
    df['RNA_ui'] = RNA_ui
    df['DNA_multiplicity'] = DNA_multiplicity
    df['RNA_multiplicity'] = RNA_multiplicity
    df['DNA_CCF'] = DNA_CCF
    df['RNA_CCF'] = RNA_CCF

    df = df[df['DNA_CCF']>=DNA_CCF_cutoff]

    df.to_csv(fname, index=False)

#df = df[df['Variant_Classification'].isin(['Silent'])]
df = df[df['Hugo_Symbol']!="MXRA5"]

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
#axs[4].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="red", linewidth=0.5, edgecolor="none")
#axs[5].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="red", linewidth=0.5, edgecolor="none")

sns.histplot(x=ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), kde=False, bins=bins_list, ax=axs[2], color="red", linewidth=0.5, edgecolor="none", line_kws={'lw':0.4}, kde_kws={'bw_adjust':0.15})
#sns.histplot(x=ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), kde=False, bins=bins_list, ax=axs[5], color="red", linewidth=0.5, edgecolor="none", line_kws={'lw':0.4}, kde_kws={'bw_adjust':0.15})

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

"""
axs[4].spines['bottom'].set_visible(False)
axs[5].spines['top'].set_visible(False)
axs[4].xaxis.tick_top()
axs[4].tick_params(labeltop=False)  # don't put tick labels at the top
axs[5].xaxis.tick_bottom()

plt.xlim([0, 3])
"""

plt.xlabel("RNA VAF / DNA VAF")

plt.ylabel("Number of Variants")










plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(5,5)})
sns.set(font_scale=1)
dpi_set = 72

final_name = "XISTneg_males"
fname = "/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/out_" + final_name + "_" + gene_class + "_chrX_CN" + str(cn_to_analyze) + ".csv"

if os.path.isfile(fname):
    df = pd.read_csv(fname)
else:
    if final_name == "XISTpos_males" or final_name == "XISTneg_males":
        df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_chrX_muts_REVISED_males_unaveraged_TPM_geq2.csv")    
    elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
        df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/IMRAN_output_chrX_muts_REVISED_females_unaveraged_TPM_geq2.csv")
        
    invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

    df = df[(df['Hugo_Symbol'].isin(genes_of_interest))] 
    df.drop_duplicates(subset=['Hugo_Symbol','Start_Position','Tumor_Sample_Barcode'], keep="first", inplace=True)

    for a in invalid_ids:
        df = df[~(df['Truncated_Barcodes'].str.contains(a))]

    secondary_barcodes = []

    for a in df['Truncated_Barcodes'].tolist():
        secondary_barcodes.append(split_advanced(a, "-", 3)[0])

    df['Secondary_Barcodes'] = secondary_barcodes

    df = df[df['Truncated_Barcodes'].isin(cn_select_samples)]


    df['total_RNA_counts'] = df['RNA_alt_counts'] + df['RNA_ref_counts']
    df['total_DNA_counts'] = df['DNA_alt_counts'] + df['DNA_ref_counts']
    df['DNA_VAF_purity_ratio'] = df['DNA_VAF']/df['Purity']

    cutoff = 20

    df = df[df['total_RNA_counts']>=cutoff]
    df = df[df['total_DNA_counts']>=cutoff]

    #df = df[df['DNA_VAF_purity_ratio']>=0.4]

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

    ids = df['Truncated_Barcodes'].tolist()
    start = df['Start_Position'].tolist()

    cn_temp = []

    print(len(ids))

    new_purity = []

    a=0
    while a<len(ids):
        temp_id = split_advanced(ids[a], "-", 3)[0]
        new_purity.append(purity_dict[temp_id])
        temp_df = df_seg[df_seg['Sample'].str.contains(temp_id)]
        start_vals = temp_df['Start'].tolist()
        cn_vals = temp_df['Copy_Number'].tolist()
        q=0
        while (q-1)<len(start_vals):
            if q == len(start_vals):
                cn_temp.append(cn_vals[q-1])
                break
            elif start[a] > start_vals[q]:
                q=q+1
            else:
                cn_temp.append(cn_vals[q-1])
                break
        a=a+1
        if a % 100 == 0:
            print(a)

    df['Purity'] = new_purity
    df['CN_temp'] = cn_temp

    #CCF calculations: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5538405/

    DNA_VAF = df['DNA_VAF'].tolist()
    RNA_VAF = df['RNA_VAF'].tolist()
    purity = df['Purity'].tolist()

    DNA_CCF = []
    RNA_CCF = []
    DNA_ui = []
    RNA_ui = []
    DNA_multiplicity = []
    RNA_multiplicity = []

    if final_name == "XISTpos_males" or final_name == "XISTneg_males":
        normal_cn = 1
    elif final_name == "XISTpos_females" or final_name == "XISTneg_females":
        normal_cn = 2

    a=0
    while a<len(cn_temp):
        temp_DNA_ui = DNA_VAF[a]*1/purity[a]*(purity[a]*cn_temp[a]+(1-purity[a])*normal_cn)
        temp_RNA_ui = RNA_VAF[a]*1/purity[a]*(purity[a]*cn_temp[a]+(1-purity[a])*normal_cn)

        DNA_ui.append(temp_DNA_ui)
        RNA_ui.append(temp_RNA_ui)    

        if temp_DNA_ui >= 1:
            DNA_multiplicity.append(temp_DNA_ui)
            DNA_CCF.append(1)
        else:
            DNA_multiplicity.append(1)
            DNA_CCF.append(temp_DNA_ui)

        if temp_RNA_ui >= 1:
            RNA_multiplicity.append(temp_RNA_ui)
            RNA_CCF.append(1)
        else:
            RNA_multiplicity.append(1)
            RNA_CCF.append(temp_RNA_ui)

        a=a+1

    df['DNA_ui'] = DNA_ui
    df['RNA_ui'] = RNA_ui
    df['DNA_multiplicity'] = DNA_multiplicity
    df['RNA_multiplicity'] = RNA_multiplicity
    df['DNA_CCF'] = DNA_CCF
    df['RNA_CCF'] = RNA_CCF

    df = df[df['DNA_CCF']>=DNA_CCF_cutoff]

    df.to_csv(fname, index=False)

#df = df[df['Variant_Classification'].isin(['Silent'])]
df = df[df['Hugo_Symbol']!="MXRA5"]

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
#axs[6].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="blue", linewidth=0.5, edgecolor="none")
#axs[7].hist(ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), density=False, bins=bins_list, color="blue", linewidth=0.5, edgecolor="none")

sns.histplot(x=ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), kde=False, bins=bins_list, ax=axs[3], color="blue", linewidth=0.5, edgecolor="none", line_kws={'lw':0.4}, kde_kws={'bw_adjust':0.23})
#sns.histplot(x=ratio_list, weights=np.ones(len(ratio_list)) / len(ratio_list), kde=False, bins=bins_list, ax=axs[7], color="blue", linewidth=0.5, edgecolor="none", line_kws={'lw':0.4}, kde_kws={'bw_adjust':0.23})

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
"""
axs[6].spines['bottom'].set_visible(False)
axs[7].spines['top'].set_visible(False)
axs[6].xaxis.tick_top()
axs[6].tick_params(labeltop=False)  # don't put tick labels at the top
axs[7].xaxis.tick_bottom()

plt.xlim([0, 3])
"""

plt.xlabel(r'VAF$_{\rmRNA–Seq/WES}$ for chrX Coding Somatic ' + gene_class + ' Variants', fontsize=9)

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
    
    axs[q].set_xlim([0,3.02])
    axs[q].set_ylim([0,0.41])



    """
    if q == 1 or q == 5:
        axs[q].set_yticks([0, 0.05])     
    elif q == 3:
        axs[q].set_yticks([0, 0.08])             
    elif q == 7:
        axs[q].set_yticks([0, 0.05])  
    elif q == 0 or q == 2:
        axs[q].set_yticks([0.18, 0.30])  
    else:
        axs[q].set_yticks([0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3, 0.33, 0.36, 0.39])    
    """

    axs[q].set_xticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
    #axs[q].yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
    axs[q].yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        
    axs[q].axvline(n_cutoff, lw=0.4, c="black", ls="--")
    
    axs[q].tick_params(left='on')
    axs[q].tick_params(bottom='on')
    axs[q].set_ylabel("Percentage\nof Variants", fontsize=9, labelpad=10)
    if q != 3:
        axs[q].set_xticks([])
    
    axs[q].grid(False)
    axs[q].set_facecolor("white")
    axs[q].spines['bottom'].set_color('0')
    axs[q].spines['left'].set_color('0')

    
    """
    if q == 1:
        axs[q].yaxis.labelpad = 12
    if q == 2:
        axs[q].yaxis.labelpad = 8.2
    """        
    #plt.ylim([0, 1.02])
    q=q+1

counter1 = 0

for a in XISTneg_male_ratio:
    if a < n_cutoff:
        counter1 = counter1 + 1

counter2 = 0

for a in XISTpos_male_ratio:
    if a < n_cutoff:
        counter2 = counter2 + 1

counter3 = 0

for a in XISTneg_female_ratio:
    if a < n_cutoff:
        counter3 = counter3 + 1

counter4 = 0

for a in XISTpos_female_ratio:
    if a < n_cutoff:
        counter4 = counter4 + 1

# can we assume anything from our sample
# our samples - 82% are good in one, and ~79% are good in the other
# note - the samples do not need to be the same size
sample_success_a, sample_size_a = (counter3, len(XISTneg_female_ratio))
sample_success_b, sample_size_b = (counter4, len(XISTpos_female_ratio))

# check our sample against Ho for Ha != Ho
successes = np.array([sample_success_a, sample_success_b])
samples = np.array([sample_size_a, sample_size_b])
# note, no need for a Ho value here - it's derived from the other parameters
stat, p_value = proportions_ztest(count=successes, nobs=samples,  alternative='two-sided')
# report
print('FEMALE - z_stat: %0.3f, p_value: %0.3f' % (stat, p_value))
print(successes)
print(samples)

group1 = XISTneg_female_ratio
group2 = XISTpos_female_ratio
counter_success = counter3
N_variants_to_sample = len(group1)

prop_val = counter_success/len(group1)

counter_above = 0

a=0
while a<n_iterations:
    new_list = random.sample(group2, N_variants_to_sample)
    prop_to_compare = len([x for x in new_list if x < n_cutoff])/N_variants_to_sample
    if prop_val >= prop_to_compare:
        counter_above = counter_above + 1
    if a % 10000 == 0:
        print(a)
    a=a+1

print("Resampled P-val: " + str(counter_above/n_iterations))

# can we assume anything from our sample
# our samples - 82% are good in one, and ~79% are good in the other
# note - the samples do not need to be the same size
sample_success_a, sample_size_a = (counter1, len(XISTneg_male_ratio))
sample_success_b, sample_size_b = (counter2, len(XISTpos_male_ratio))
# check our sample against Ho for Ha != Ho
successes = np.array([sample_success_a, sample_success_b])
samples = np.array([sample_size_a, sample_size_b])

# note, no need for a Ho value here - it's derived from the other parameters
stat, p_value = proportions_ztest(count=successes, nobs=samples,  alternative='two-sided')
# report
print('MALE - z_stat: %0.3f, p_value: %0.3f' % (stat, p_value))
print(successes)
print(samples)

group1 = XISTpos_male_ratio
group2 = XISTneg_male_ratio
counter_success = counter2
N_variants_to_sample = len(group1)

prop_val = counter_success/len(group1)

counter_above = 0

a=0
while a<n_iterations:
    new_list = random.sample(group2, N_variants_to_sample)
    prop_to_compare = len([x for x in new_list if x < n_cutoff])/N_variants_to_sample
    if prop_val <= prop_to_compare:
        counter_above = counter_above + 1
    if a % 10000 == 0:
        print(a)
    a=a+1

print("Resampled P-val: " + str(counter_above/n_iterations))

"""
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
    

print(len(XISTneg_male_ratio))
print(len(XISTpos_male_ratio))
print(len(XISTneg_female_ratio))
print(len(XISTpos_female_ratio))

"""
"""

counter=0

for a in ratio_list:
    if a<0.5:
        counter=counter+1

leq_point_5_ratio = str("{:.1f}".format(counter/len(ratio_list)*100) + "% of Variants with RNA VAF / DNA VAF < 0.5")

ax.text(0.07, 1, leq_point_5_ratio, transform=ax.transAxes,
     fontsize=12, verticalalignment='top', color="black")
"""

plt.subplots_adjust(wspace=0, hspace=0.3)

if cn_to_analyze == 2:
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_and_FEMALES_CCF_density_plots_test_cn2_" + gene_class + ".png", dpi=dpi_set)
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_and_FEMALES_CCF_density_plots_test_cn2_" + gene_class + ".pdf", dpi=dpi_set)
elif cn_to_analyze == 1:
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_and_FEMALES_CCF_density_plots_test_cn1_" + gene_class + ".png", dpi=dpi_set)
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_and_FEMALES_CCF_density_plots_test_cn1_" + gene_class + ".pdf", dpi=dpi_set)

