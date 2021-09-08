#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 13:34:25 2021

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

cutoff_value = 0

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(3.3,8.5*9/7)})
sns.set(font_scale=0.6)
plt.rcParams.update({'font.size': 14})
plt.rcParams['axes.linewidth'] = 0.3
dpi_set = 72 # change the output resolution

full_class_list = ['PanCan']

df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_averaged_TCGA_Xena_TPM.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]

    
"""

if temp_class != 'TGCT-NS':

    df_XIST = df_XIST[df_XIST['Classification']==temp_class]
    
else:
    
    df_XIST = df_XIST[df_XIST['Secondary_Class']==temp_class]  
    
"""

#df_XIST = df_XIST[df_XIST['Secondary_Class']!=temp_class] 

    #df_XIST = df_XIST[df_XIST['Classification'].isin(['LUSC', 'LUAD', 'SKCM', 'LIHC', 'ESCA', 'STAD', ''])] #Pan-lung, SKCM, LIHC, and Pan-GI (top 4)

df_high = df_XIST[df_XIST['XIST_TPM']>=3]
df_low = df_XIST[df_XIST['XIST_TPM']<3]

h_samples = ['TCGA-2G-AAFZ-01'] #TCGA-2G-AAFZ-01

l_samples = ['TCGA-2G-AAL5-01']

df_high = h_samples

df_low = l_samples
        
df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.processed.txt", sep="\t")

uq_samples = df2['array'].tolist()

chrX_cn = df2['purity'].tolist()

cn_dict = dict(zip(uq_samples, chrX_cn))

pur_high = []

for a in df_high:
    try:
        pur_high.append(cn_dict[a])
    except KeyError:
        print("KEY ERROR " + str(a))
        continue
        
pur_low = []

for a in df_low:
    try:
        pur_low.append(cn_dict[a])
    except KeyError:
        print("KEY ERROR " + str(a))
        continue
    
print(statistics.median(pur_high))
print(statistics.median(pur_low))

print("READING methylation")

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_only_averaged_jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.txt", sep="\t")

print("READ methylation")

df_annot = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_annotation_ref.csv")

#df_annot = df_annot[df_annot['Feature_Type'] == "Island"] #Restrict to CpG islands

comp_element_ref = df_annot['Composite Element REF'].tolist()
chromosome = df_annot['Chromosome'].tolist()
start = df_annot['Start'].tolist()
end = df_annot['End'].tolist()
gene_list = df_annot['Gene_Symbol'].tolist()
cgi_list = df_annot['CGI_Coordinate'].tolist()
feature_list = df_annot['Feature_Type'].tolist()

chr_dict = dict(zip(comp_element_ref, chromosome))
start_dict = dict(zip(comp_element_ref, start))
end_dict = dict(zip(comp_element_ref, end))
gene_dict = dict(zip(comp_element_ref, gene_list))
cgi_dict = dict(zip(comp_element_ref, cgi_list))
feature_dict = dict(zip(comp_element_ref, feature_list))

df_RNA = df_RNA[df_RNA['Composite_Element_REF'].isin(comp_element_ref)]

orig_comp_elements = df_RNA['Composite_Element_REF'].tolist()

chr_new = []
start_new = []
end_new = []
gene_new = []
cgi_new = []
feature_new = []

for a in orig_comp_elements:
    chr_new.append(chr_dict[a])
    start_new.append(start_dict[a])
    end_new.append(end_dict[a])
    gene_new.append(gene_dict[a])
    cgi_new.append(cgi_dict[a])
    feature_new.append(feature_dict[a])

print("inserting columns")

df_RNA.insert(loc=0, column='End', value=end_new)
df_RNA.insert(loc=0, column='Start', value=start_new)
df_RNA.insert(loc=0, column='Chromosome', value=chr_new)
df_RNA.insert(loc=0, column='Gene', value=gene_new)
df_RNA.insert(loc=0, column='CGI_Coordinate', value=cgi_new)
df_RNA.insert(loc=0, column='Feature_Type', value=feature_new)

#df_RNA = df_RNA[df_RNA['Feature_Type']=="Island"]

df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))
df_low = list(set(df_low) & set(df_RNA_temp_cols))

cpg_promoter_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_non_escaping_CpG_promoter_islands.csv")

islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()

#df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_X = df_RNA[df_RNA['Composite_Element_REF'].isin(islands_of_interest)]

df_autosome = df_RNA[~(df_RNA['Chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, temp_sample):
    
    sample_val = []
    
    gene_ids = df_inp['Composite_Element_REF'].tolist()
    
    temp_df = df_inp[temp_sample]
    
    temp_vals = temp_df.values.tolist()
    
    temp_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        if temp_vals[a] >= cutoff_value:
            
            temp_total_TPM = temp_total_TPM + temp_vals[a]
            
            sample_val.append(temp_vals[a])
        
        a=a+1
        
    median_val = sample_val

    return median_val, temp_total_TPM

h_median_list = []
h_total_TPM_list = []

for a in df_high:
    temp_median, total_TPM = DGE(df_X, a)
    h_median_list.append(temp_median)
    h_total_TPM_list.append(total_TPM)
    
l_median_list = []
l_total_TPM_list = []

for a in df_low:
    temp_median, total_TPM = DGE(df_X, a)
    l_median_list.append(temp_median)
    l_total_TPM_list.append(total_TPM)

h_median_list_auto = []
h_total_TPM_list_auto = []

for a in df_high:
    temp_median, total_TPM = DGE(df_autosome, a)
    h_median_list_auto.append(temp_median)
    h_total_TPM_list_auto.append(total_TPM)

l_median_list_auto = []
l_total_TPM_list_auto = []

for a in df_low:
    temp_median, total_TPM = DGE(df_autosome, a)
    l_median_list_auto.append(temp_median)
    l_total_TPM_list_auto.append(total_TPM)  
    
X_high = h_median_list[0]
X_low = l_median_list[0]
other_high = h_median_list_auto[0]
other_low = l_median_list_auto[0]

X_high_med = statistics.median(X_high)

X_low_med = statistics.median(X_low)

n_bins = 100

fig, axs = plt.subplots(9, 1, sharex=True, tight_layout=True)

x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[4],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[4].fill_between(x, 0,  y, facecolor='#DACAFF', alpha=0.8)




x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[4],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[4].fill_between(x, 0,  y, facecolor='#008BB5', alpha=0.2)





x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[3],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[3].fill_between(x, 0,  y, facecolor='#DACAFF', alpha=0.8)



x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[3],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
     
    
axs[3].fill_between(x, 0,  y, facecolor='#008BB5', alpha=0.2)





other_high_med = statistics.median(other_high)

other_low_med = statistics.median(other_low)
high_ratio = "chrX Median β = " + str("{:.2f}".format(X_high_med))
auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_high_med))

"""
axs[4].text(0.7, 1.05, high_ratio, transform=axs[4].transAxes,
     fontsize=6, verticalalignment='top')

axs[4].text(0.62, 0.91, auto_ratio, transform=axs[4].transAxes,
     fontsize=6, verticalalignment='top')
"""

low_ratio = "chrX Median β = " + str("{:.2f}".format(X_low_med))
auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_low_med))

"""
axs[3].text(0.62, 0.91, auto_ratio, transform=axs[3].transAxes,
     fontsize=6, verticalalignment='top')

axs[3].text(0.7, 1.05, low_ratio, transform=axs[3].transAxes,
     fontsize=6, verticalalignment='top')
"""

axs[4].axvline(other_high_med, c="#C4ABFF", ls="--")
axs[3].axvline(other_low_med, c="#C4ABFF", ls="--")

axs[4].axvline(X_high_med, c="#008BB5", ls="--")
axs[3].axvline(X_low_med, c="#008BB5", ls="--")

axs[4].set_ylabel("Density")
axs[3].set_ylabel("Density")

#plt.xlim(0, 12.5)

axs[4].grid(False)
axs[4].set_facecolor("white")
axs[4].spines['bottom'].set_color('0')
axs[4].spines['left'].set_color('0')

axs[4].tick_params(bottom='on', left='on')

axs[3].grid(False)
axs[3].set_facecolor("white")
axs[3].spines['bottom'].set_color('0')
axs[3].spines['left'].set_color('0')

axs[5].tick_params(bottom='on', left='on')
    




#LUNG


full_class_list = ['LUSC'] #SKCM, LIHC, Pan-lung, Pan-GI (top 4)

h_samples = ['TCGA-56-7822-01']

l_samples = ['TCGA-77-6843-01'] #TCGA-85-8277-01, TCGA-77-A5G3-01

df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.processed.txt", sep="\t")

uq_samples = df2['array'].tolist()

chrX_cn = df2['purity'].tolist()

cn_dict = dict(zip(uq_samples, chrX_cn))

pur_high = []

#h_samples = df_high['Barcode'].tolist()

#l_samples = df_low['Barcode'].tolist()

df_high = h_samples

df_low = l_samples


cpg_promoter_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_non_escaping_CpG_promoter_islands.csv")

islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()

#df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_X = df_RNA[df_RNA['Composite_Element_REF'].isin(islands_of_interest)]

df_autosome = df_RNA[~(df_RNA['Chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, temp_sample):
    
    sample_val = []
    
    gene_ids = df_inp['Composite_Element_REF'].tolist()
    
    temp_df = df_inp[temp_sample]
    
    temp_vals = temp_df.values.tolist()
    
    temp_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        if temp_vals[a] >= cutoff_value:
            
            temp_total_TPM = temp_total_TPM + temp_vals[a]
            
            sample_val.append(temp_vals[a])
        
        a=a+1
        
    median_val = sample_val

    return median_val, temp_total_TPM

h_median_list = []
h_total_TPM_list = []

for a in df_high:
    temp_median, total_TPM = DGE(df_X, a)
    h_median_list.append(temp_median)
    h_total_TPM_list.append(total_TPM)
    
l_median_list = []
l_total_TPM_list = []

for a in df_low:
    temp_median, total_TPM = DGE(df_X, a)
    l_median_list.append(temp_median)
    l_total_TPM_list.append(total_TPM)

h_median_list_auto = []
h_total_TPM_list_auto = []

for a in df_high:
    temp_median, total_TPM = DGE(df_autosome, a)
    h_median_list_auto.append(temp_median)
    h_total_TPM_list_auto.append(total_TPM)

l_median_list_auto = []
l_total_TPM_list_auto = []

for a in df_low:
    temp_median, total_TPM = DGE(df_autosome, a)
    l_median_list_auto.append(temp_median)
    l_total_TPM_list_auto.append(total_TPM)  
    
X_high = h_median_list[0]
X_low = l_median_list[0]
other_high = h_median_list_auto[0]
other_low = l_median_list_auto[0]

X_high_med = statistics.median(X_high)

X_low_med = statistics.median(X_low)

n_bins = 100

x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[6],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[6].fill_between(x, 0,  y, facecolor='#DACAFF', alpha=0.8)






x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[6],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[6].fill_between(x, 0,  y, facecolor='#008BB5', alpha=0.2)






x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
        
axs[5].fill_between(x, 0,  y, facecolor='#DACAFF', alpha=0.8)




x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[1].get_data()


x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[5].fill_between(x, 0,  y, facecolor='#008BB5', alpha=0.2)




other_high_med = statistics.median(other_high)

other_low_med = statistics.median(other_low)

high_ratio = "chrX Median β = " + str("{:.2f}".format(X_high_med))

auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_high_med))

"""
axs[6].text(0.62, 0.91, auto_ratio, transform=axs[6].transAxes,
     fontsize=6, verticalalignment='top')

axs[6].text(0.7, 1.05, high_ratio, transform=axs[6].transAxes,
     fontsize=6, verticalalignment='top')
"""

low_ratio = "chrX Median β = " + str("{:.2f}".format(X_low_med))

auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_low_med))

"""
axs[5].text(0.62, 0.91, auto_ratio, transform=axs[5].transAxes,
     fontsize=6, verticalalignment='top')

axs[5].text(0.7, 1.05, low_ratio, transform=axs[5].transAxes,
     fontsize=6, verticalalignment='top')
"""

axs[6].axvline(other_high_med, c="#C4ABFF", ls="--")
axs[5].axvline(other_low_med, c="#C4ABFF", ls="--")

axs[6].axvline(X_high_med, c="#008BB5", ls="--")
axs[5].axvline(X_low_med, c="#008BB5", ls="--")

axs[6].set_ylabel("Density")
axs[6].set_xlabel("β at Probe")
axs[5].set_ylabel("Density")

#plt.xlim(0, 12.5)

axs[6].grid(False)
axs[6].set_facecolor("white")
axs[6].spines['bottom'].set_color('0')
axs[6].spines['left'].set_color('0')

axs[6].tick_params(bottom='on', left='on')

axs[5].grid(False)
axs[5].set_facecolor("white")
axs[5].spines['bottom'].set_color('0')
axs[5].spines['left'].set_color('0')
axs[5].tick_params(bottom='on', left='on')






#SGCT

temp_class = 'TGCT-S'

df_high = ['TCGA-XY-A8S3-01']

df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.processed.txt", sep="\t")

uq_samples = df2['array'].tolist()

chrX_cn = df2['purity'].tolist()

cn_dict = dict(zip(uq_samples, chrX_cn))

pur_high = []

for a in df_high:
    try:
        pur_high.append(cn_dict[a])
    except KeyError:
        continue
        
pur_low = []

for a in df_low:
    try:
        pur_low.append(cn_dict[a])
    except KeyError:
        continue
    
df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))

cpg_promoter_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_non_escaping_CpG_promoter_islands.csv")

islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()

#df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_X = df_RNA[df_RNA['Composite_Element_REF'].isin(islands_of_interest)]

df_autosome = df_RNA[~(df_RNA['Chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, temp_sample):
    
    sample_val = []
    
    gene_ids = df_inp['Composite_Element_REF'].tolist()
    
    temp_df = df_inp[temp_sample]
    
    temp_vals = temp_df.values.tolist()
    
    temp_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        if temp_vals[a] >= cutoff_value:
            
            temp_total_TPM = temp_total_TPM + temp_vals[a]
            
            sample_val.append(temp_vals[a])
        
        a=a+1
        
    median_val = sample_val

    return median_val, temp_total_TPM

h_median_list = []
h_total_TPM_list = []

for a in df_high:
    temp_median, total_TPM = DGE(df_X, a)
    h_median_list.append(temp_median)
    h_total_TPM_list.append(total_TPM)

h_median_list_auto = []
h_total_TPM_list_auto = []

for a in df_high:
    temp_median, total_TPM = DGE(df_autosome, a)
    h_median_list_auto.append(temp_median)
    h_total_TPM_list_auto.append(total_TPM)  
    
X_high = h_median_list[0]
other_high = h_median_list_auto[0]

X_high_med = statistics.median(X_high)

n_bins = 100

x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[2],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()


x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[2].fill_between(x, 0,  y, facecolor='#DACAFF', alpha=0.8)




x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[2],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[2].fill_between(x, 0,  y, facecolor='#008BB5', alpha=0.2)



other_high_med = statistics.median(other_high)

high_ratio = "chrX Median β = " + str("{:.2f}".format(X_high_med))

auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_high_med))

"""
axs[2].text(0.62, 0.91, auto_ratio, transform=axs[2].transAxes,
     fontsize=6, verticalalignment='top')

axs[2].text(0.7, 1.05, high_ratio, transform=axs[2].transAxes,
     fontsize=6, verticalalignment='top')
"""

axs[2].axvline(other_high_med, c="#C4ABFF", ls="--")
axs[2].axvline(X_high_med, c="#008BB5", ls="--")

axs[2].set_ylabel("Density")

#plt.xlim(0, 12.5)

axs[2].grid(False)
axs[2].set_facecolor("white")
axs[2].spines['bottom'].set_color('0')
axs[2].spines['left'].set_color('0')

axs[2].tick_params(bottom='on', left='on')







"""
#SGCT

temp_class = 'TGCT-NS'

df_high = ['TCGA-SN-A6IS-01']

df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.processed.txt", sep="\t")

uq_samples = df2['array'].tolist()

chrX_cn = df2['purity'].tolist()

cn_dict = dict(zip(uq_samples, chrX_cn))

pur_high = []

for a in df_high:
    try:
        pur_high.append(cn_dict[a])
    except KeyError:
        continue
        
pur_low = []

for a in df_low:
    try:
        pur_low.append(cn_dict[a])
    except KeyError:
        continue
    
df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))

cpg_promoter_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_non_escaping_CpG_promoter_islands.csv")

islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()

#df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_X = df_RNA[df_RNA['Composite_Element_REF'].isin(islands_of_interest)]

df_autosome = df_RNA[~(df_RNA['Chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, temp_sample):
    
    sample_val = []
    
    gene_ids = df_inp['Composite_Element_REF'].tolist()
    
    temp_df = df_inp[temp_sample]
    
    temp_vals = temp_df.values.tolist()
    
    temp_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        if temp_vals[a] >= cutoff_value:
            
            temp_total_TPM = temp_total_TPM + temp_vals[a]
            
            sample_val.append(temp_vals[a])
        
        a=a+1
        
    median_val = sample_val

    return median_val, temp_total_TPM

h_median_list = []
h_total_TPM_list = []

for a in df_high:
    temp_median, total_TPM = DGE(df_X, a)
    h_median_list.append(temp_median)
    h_total_TPM_list.append(total_TPM)

h_median_list_auto = []
h_total_TPM_list_auto = []

for a in df_high:
    temp_median, total_TPM = DGE(df_autosome, a)
    h_median_list_auto.append(temp_median)
    h_total_TPM_list_auto.append(total_TPM)  
    
X_high = h_median_list[0]
other_high = h_median_list_auto[0]

X_high_med = statistics.median(X_high)

n_bins = 100

sns.distplot(other_high, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

sns.distplot(X_high, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

other_high_med = statistics.median(other_high)

high_ratio = "chrX Median β = " + str("{:.2f}".format(X_high_med))

auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_high_med))

axs[5].text(0.623, 0.8, auto_ratio, transform=axs[5].transAxes,
     fontsize=6, verticalalignment='top')

axs[5].text(0.7, 0.98, high_ratio, transform=axs[5].transAxes,
     fontsize=6, verticalalignment='top')

axs[5].axvline(other_high_med, c="#DACAFF", ls="--")
axs[5].axvline(X_high_med, c="#008BB5", ls="--")

axs[5].set_ylabel("Density")

#plt.xlim(0, 12.5)

axs[5].grid(False)
axs[5].set_facecolor("white")
axs[5].spines['bottom'].set_color('0')
axs[5].spines['left'].set_color('0')

axs[5].tick_params(bottom='on', left='on')





temp_class = 'LUSC'

df_high = ['TCGA-37-4130-01']

df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.processed.txt", sep="\t")

uq_samples = df2['array'].tolist()

chrX_cn = df2['purity'].tolist()

cn_dict = dict(zip(uq_samples, chrX_cn))

pur_high = []

for a in df_high:
    try:
        pur_high.append(cn_dict[a])
    except KeyError:
        continue
        
pur_low = []

for a in df_low:
    try:
        pur_low.append(cn_dict[a])
    except KeyError:
        continue
    
df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))

cpg_promoter_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_non_escaping_CpG_promoter_islands.csv")

islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()

#df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_X = df_RNA[df_RNA['Composite_Element_REF'].isin(islands_of_interest)]

df_autosome = df_RNA[~(df_RNA['Chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, temp_sample):
    
    sample_val = []
    
    gene_ids = df_inp['Composite_Element_REF'].tolist()
    
    temp_df = df_inp[temp_sample]
    
    temp_vals = temp_df.values.tolist()
    
    temp_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        if temp_vals[a] >= cutoff_value:
            
            temp_total_TPM = temp_total_TPM + temp_vals[a]
            
            sample_val.append(temp_vals[a])
        
        a=a+1
        
    median_val = sample_val

    return median_val, temp_total_TPM

h_median_list = []
h_total_TPM_list = []

for a in df_high:
    temp_median, total_TPM = DGE(df_X, a)
    h_median_list.append(temp_median)
    h_total_TPM_list.append(total_TPM)

h_median_list_auto = []
h_total_TPM_list_auto = []

for a in df_high:
    temp_median, total_TPM = DGE(df_autosome, a)
    h_median_list_auto.append(temp_median)
    h_total_TPM_list_auto.append(total_TPM)  
    
X_high = h_median_list[0]
other_high = h_median_list_auto[0]

X_high_med = statistics.median(X_high)

n_bins = 100

sns.distplot(other_high, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

sns.distplot(X_high, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

other_high_med = statistics.median(other_high)

high_ratio = "chrX Median β = " + str("{:.2f}".format(X_high_med))

auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_high_med))

axs[8].text(0.623, 0.8, auto_ratio, transform=axs[8].transAxes,
     fontsize=6, verticalalignment='top')

axs[8].text(0.7, 0.98, high_ratio, transform=axs[8].transAxes,
     fontsize=6, verticalalignment='top')

axs[8].axvline(other_high_med, c="#DACAFF", ls="--")
axs[8].axvline(X_high_med, c="#008BB5", ls="--")

axs[8].set_ylabel("Density")

#plt.xlim(0, 12.5)

axs[8].grid(False)
axs[8].set_facecolor("white")
axs[8].spines['bottom'].set_color('0')
axs[8].spines['left'].set_color('0')

axs[8].tick_params(bottom='on', left='on')

"""




#Normals (Cancer-Adjacent)

df_low = ['TCGA-22-5482-11']

print("READING methylation")

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_normals_only_averaged_jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.txt", sep="\t")

print("READ methylation")

df_annot = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_annotation_ref.csv")

#df_annot = df_annot[df_annot['Feature_Type'] == "Island"] #Restrict to CpG islands

comp_element_ref = df_annot['Composite Element REF'].tolist()
chromosome = df_annot['Chromosome'].tolist()
start = df_annot['Start'].tolist()
end = df_annot['End'].tolist()
gene_list = df_annot['Gene_Symbol'].tolist()
cgi_list = df_annot['CGI_Coordinate'].tolist()
feature_list = df_annot['Feature_Type'].tolist()

chr_dict = dict(zip(comp_element_ref, chromosome))
start_dict = dict(zip(comp_element_ref, start))
end_dict = dict(zip(comp_element_ref, end))
gene_dict = dict(zip(comp_element_ref, gene_list))
cgi_dict = dict(zip(comp_element_ref, cgi_list))
feature_dict = dict(zip(comp_element_ref, feature_list))

df_RNA = df_RNA[df_RNA['Composite_Element_REF'].isin(comp_element_ref)]

orig_comp_elements = df_RNA['Composite_Element_REF'].tolist()

chr_new = []
start_new = []
end_new = []
gene_new = []
cgi_new = []
feature_new = []

for a in orig_comp_elements:
    chr_new.append(chr_dict[a])
    start_new.append(start_dict[a])
    end_new.append(end_dict[a])
    gene_new.append(gene_dict[a])
    cgi_new.append(cgi_dict[a])
    feature_new.append(feature_dict[a])

print("inserting columns")

df_RNA.insert(loc=0, column='End', value=end_new)
df_RNA.insert(loc=0, column='Start', value=start_new)
df_RNA.insert(loc=0, column='Chromosome', value=chr_new)
df_RNA.insert(loc=0, column='Gene', value=gene_new)
df_RNA.insert(loc=0, column='CGI_Coordinate', value=cgi_new)
df_RNA.insert(loc=0, column='Feature_Type', value=feature_new)

#df_RNA = df_RNA[df_RNA['Feature_Type']=="Island"]

cpg_promoter_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_non_escaping_CpG_promoter_islands.csv")

islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()

#df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_X = df_RNA[df_RNA['Composite_Element_REF'].isin(islands_of_interest)]

df_autosome = df_RNA[~(df_RNA['Chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, temp_sample):
    
    sample_val = []
    
    gene_ids = df_inp['Composite_Element_REF'].tolist()
    
    temp_df = df_inp[temp_sample]
    
    temp_vals = temp_df.values.tolist()
    
    temp_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        if temp_vals[a] >= cutoff_value:
            
            temp_total_TPM = temp_total_TPM + temp_vals[a]
            
            sample_val.append(temp_vals[a])
        
        a=a+1
        
    median_val = sample_val

    return median_val, temp_total_TPM
    
l_median_list = []
l_total_TPM_list = []

for a in df_low:
    temp_median, total_TPM = DGE(df_X, a)
    l_median_list.append(temp_median)
    l_total_TPM_list.append(total_TPM)

l_median_list_auto = []
l_total_TPM_list_auto = []

for a in df_low:
    temp_median, total_TPM = DGE(df_autosome, a)
    l_median_list_auto.append(temp_median)
    l_total_TPM_list_auto.append(total_TPM)  
    
X_low = l_median_list[0]
other_low = l_median_list_auto[0]

X_low_med = statistics.median(X_low)

n_bins = 100

x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[1],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[1].fill_between(x, 0,  y, facecolor='#DACAFF', alpha=0.8)






x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[1],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[1].fill_between(x, 0,  y, facecolor='#008BB5', alpha=0.2)





male_autosome_dist = other_low
male_chrX_dist = X_low

other_low_med = statistics.median(other_low)

low_ratio = "chrX Median β = " + str("{:.2f}".format(X_low_med))

auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_low_med))

"""
axs[1].text(0.62, 0.91, auto_ratio, transform=axs[1].transAxes,
     fontsize=6, verticalalignment='top')

axs[1].text(0.7, 1.05, low_ratio, transform=axs[1].transAxes,
     fontsize=6, verticalalignment='top')
"""

axs[1].axvline(other_low_med, c="#C4ABFF", ls="--")
axs[1].axvline(X_low_med, c="#008BB5", ls="--")
axs[1].set_ylabel("Density")

axs[1].grid(False)
axs[1].set_facecolor("white")
axs[1].spines['bottom'].set_color('0')
axs[1].spines['left'].set_color('0')

axs[1].tick_params(bottom='on', left='on')







#Normals (Cancer-Adjacent) female

df_low = ['TCGA-E9-A1NA-11']

print("READING methylation")

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_normals_only_averaged_jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.txt", sep="\t")

print("READ methylation")

df_annot = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_annotation_ref.csv")

#df_annot = df_annot[df_annot['Feature_Type'] == "Island"] #Restrict to CpG islands

comp_element_ref = df_annot['Composite Element REF'].tolist()
chromosome = df_annot['Chromosome'].tolist()
start = df_annot['Start'].tolist()
end = df_annot['End'].tolist()
gene_list = df_annot['Gene_Symbol'].tolist()
cgi_list = df_annot['CGI_Coordinate'].tolist()
feature_list = df_annot['Feature_Type'].tolist()

chr_dict = dict(zip(comp_element_ref, chromosome))
start_dict = dict(zip(comp_element_ref, start))
end_dict = dict(zip(comp_element_ref, end))
gene_dict = dict(zip(comp_element_ref, gene_list))
cgi_dict = dict(zip(comp_element_ref, cgi_list))
feature_dict = dict(zip(comp_element_ref, feature_list))

df_RNA = df_RNA[df_RNA['Composite_Element_REF'].isin(comp_element_ref)]

orig_comp_elements = df_RNA['Composite_Element_REF'].tolist()

chr_new = []
start_new = []
end_new = []
gene_new = []
cgi_new = []
feature_new = []

for a in orig_comp_elements:
    chr_new.append(chr_dict[a])
    start_new.append(start_dict[a])
    end_new.append(end_dict[a])
    gene_new.append(gene_dict[a])
    cgi_new.append(cgi_dict[a])
    feature_new.append(feature_dict[a])

print("inserting columns")

df_RNA.insert(loc=0, column='End', value=end_new)
df_RNA.insert(loc=0, column='Start', value=start_new)
df_RNA.insert(loc=0, column='Chromosome', value=chr_new)
df_RNA.insert(loc=0, column='Gene', value=gene_new)
df_RNA.insert(loc=0, column='CGI_Coordinate', value=cgi_new)
df_RNA.insert(loc=0, column='Feature_Type', value=feature_new)

#df_RNA = df_RNA[df_RNA['Feature_Type']=="Island"]

cpg_promoter_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_non_escaping_CpG_promoter_islands.csv")

islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()

#df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_X = df_RNA[df_RNA['Composite_Element_REF'].isin(islands_of_interest)]

df_autosome = df_RNA[~(df_RNA['Chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, temp_sample):
    
    sample_val = []
    
    gene_ids = df_inp['Composite_Element_REF'].tolist()
    
    temp_df = df_inp[temp_sample]
    
    temp_vals = temp_df.values.tolist()
    
    temp_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        if temp_vals[a] >= cutoff_value:
            
            temp_total_TPM = temp_total_TPM + temp_vals[a]
            
            sample_val.append(temp_vals[a])
        
        a=a+1
        
    median_val = sample_val

    return median_val, temp_total_TPM
    
l_median_list = []
l_total_TPM_list = []

for a in df_low:
    temp_median, total_TPM = DGE(df_X, a)
    l_median_list.append(temp_median)
    l_total_TPM_list.append(total_TPM)

l_median_list_auto = []
l_total_TPM_list_auto = []

for a in df_low:
    temp_median, total_TPM = DGE(df_autosome, a)
    l_median_list_auto.append(temp_median)
    l_total_TPM_list_auto.append(total_TPM)  
    
X_low = l_median_list[0]
other_low = l_median_list_auto[0]

X_low_med = statistics.median(X_low)

n_bins = 100

x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[0],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[0].fill_between(x, 0,  y, facecolor='#DACAFF', alpha=0.8)






x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[0],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[0].fill_between(x, 0,  y, facecolor='#008BB5', alpha=0.2)






male_autosome_dist = other_low
male_chrX_dist = X_low

other_low_med = statistics.median(other_low)

low_ratio = "chrX Median β = " + str("{:.2f}".format(X_low_med))

auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_low_med))

"""
axs[0].text(0.62, 0.91, auto_ratio, transform=axs[0].transAxes,
     fontsize=6, verticalalignment='top')

axs[0].text(0.7, 1.05, low_ratio, transform=axs[0].transAxes,
     fontsize=6, verticalalignment='top')
"""

axs[0].axvline(other_low_med, c="#C4ABFF", ls="--")
axs[0].axvline(X_low_med, c="#008BB5", ls="--")
axs[0].set_ylabel("Density")

axs[0].grid(False)
axs[0].set_facecolor("white")
axs[0].spines['bottom'].set_color('0')
axs[0].spines['left'].set_color('0')

axs[0].tick_params(bottom='on', left='on')

axs[3].tick_params(bottom='on', left='on')


chrX_patch = mpatches.Patch(color='#008BB5', label='chrX Non-Escapee\nPromoter CGIs')
autosome_patch = mpatches.Patch(color='#DACAFF', label='Autosomes')

leg = axs[0].legend(handles=[chrX_patch, autosome_patch], facecolor='white', fontsize=6, bbox_to_anchor=(0.59, 1))

frame = leg.get_frame()
frame.set_edgecolor("black")
frame.set_linewidth(0.5)

axs[5].tick_params(bottom='on', left='on')








df_high = df_XIST[df_XIST['XIST_TPM']>=3]
df_low = df_XIST[df_XIST['XIST_TPM']<3]

h_samples = ['TCGA-23-2078-01']

l_samples = ['TCGA-29-1702-01'] #TCGA-XE-AANI-01

df_high = h_samples

df_low = l_samples

print("READING methylation")

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_tumors_only_averaged_jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.txt", sep="\t")

print("READ methylation")

df_annot = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_annotation_ref.csv")

#df_annot = df_annot[df_annot['Feature_Type'] == "Island"] #Restrict to CpG islands

comp_element_ref = df_annot['Composite Element REF'].tolist()
chromosome = df_annot['Chromosome'].tolist()
start = df_annot['Start'].tolist()
end = df_annot['End'].tolist()
gene_list = df_annot['Gene_Symbol'].tolist()
cgi_list = df_annot['CGI_Coordinate'].tolist()
feature_list = df_annot['Feature_Type'].tolist()

chr_dict = dict(zip(comp_element_ref, chromosome))
start_dict = dict(zip(comp_element_ref, start))
end_dict = dict(zip(comp_element_ref, end))
gene_dict = dict(zip(comp_element_ref, gene_list))
cgi_dict = dict(zip(comp_element_ref, cgi_list))
feature_dict = dict(zip(comp_element_ref, feature_list))

df_RNA = df_RNA[df_RNA['Composite_Element_REF'].isin(comp_element_ref)]

orig_comp_elements = df_RNA['Composite_Element_REF'].tolist()

chr_new = []
start_new = []
end_new = []
gene_new = []
cgi_new = []
feature_new = []

for a in orig_comp_elements:
    chr_new.append(chr_dict[a])
    start_new.append(start_dict[a])
    end_new.append(end_dict[a])
    gene_new.append(gene_dict[a])
    cgi_new.append(cgi_dict[a])
    feature_new.append(feature_dict[a])

print("inserting columns")

df_RNA.insert(loc=0, column='End', value=end_new)
df_RNA.insert(loc=0, column='Start', value=start_new)
df_RNA.insert(loc=0, column='Chromosome', value=chr_new)
df_RNA.insert(loc=0, column='Gene', value=gene_new)
df_RNA.insert(loc=0, column='CGI_Coordinate', value=cgi_new)
df_RNA.insert(loc=0, column='Feature_Type', value=feature_new)

#df_RNA = df_RNA[df_RNA['Feature_Type']=="Island"]

df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))
df_low = list(set(df_low) & set(df_RNA_temp_cols))

cpg_promoter_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_non_escaping_CpG_promoter_islands.csv")

islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()

#df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_X = df_RNA[df_RNA['Composite_Element_REF'].isin(islands_of_interest)]

df_autosome = df_RNA[~(df_RNA['Chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, temp_sample):
    
    sample_val = []
    
    gene_ids = df_inp['Composite_Element_REF'].tolist()
    
    temp_df = df_inp[temp_sample]
    
    temp_vals = temp_df.values.tolist()
    
    temp_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        if temp_vals[a] >= cutoff_value:
            
            temp_total_TPM = temp_total_TPM + temp_vals[a]
            
            sample_val.append(temp_vals[a])
        
        a=a+1
        
    median_val = sample_val

    return median_val, temp_total_TPM

h_median_list = []
h_total_TPM_list = []

for a in df_high:
    temp_median, total_TPM = DGE(df_X, a)
    h_median_list.append(temp_median)
    h_total_TPM_list.append(total_TPM)
    
l_median_list = []
l_total_TPM_list = []

for a in df_low:
    temp_median, total_TPM = DGE(df_X, a)
    l_median_list.append(temp_median)
    l_total_TPM_list.append(total_TPM)

h_median_list_auto = []
h_total_TPM_list_auto = []

for a in df_high:
    temp_median, total_TPM = DGE(df_autosome, a)
    h_median_list_auto.append(temp_median)
    h_total_TPM_list_auto.append(total_TPM)

l_median_list_auto = []
l_total_TPM_list_auto = []

for a in df_low:
    temp_median, total_TPM = DGE(df_autosome, a)
    l_median_list_auto.append(temp_median)
    l_total_TPM_list_auto.append(total_TPM)  
    
X_high = h_median_list[0]
X_low = l_median_list[0]
other_high = h_median_list_auto[0]
other_low = l_median_list_auto[0]

X_high_med = statistics.median(X_high)

X_low_med = statistics.median(X_low)

n_bins = 100







#Check high vs. low





x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[8].fill_between(x, 0,  y, facecolor='#DACAFF', alpha=0.8)




x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[8].fill_between(x, 0,  y, facecolor='#008BB5', alpha=0.2)





x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[7],
             bins=n_bins, color = '#DACAFF', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[7].fill_between(x, 0,  y, facecolor='#DACAFF', alpha=0.8)



x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[7],
             bins=n_bins, color = '#008BB5', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
     
    
axs[7].fill_between(x, 0,  y, facecolor='#008BB5', alpha=0.2)





other_high_med = statistics.median(other_high)

other_low_med = statistics.median(other_low)
high_ratio = "chrX Median β = " + str("{:.2f}".format(X_high_med))
auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_high_med))

low_ratio = "chrX Median β = " + str("{:.2f}".format(X_low_med))
auto_ratio = "Autosome Median β = " + str("{:.2f}".format(other_low_med))

axs[8].axvline(other_high_med, c="#C4ABFF", ls="--")
axs[7].axvline(other_low_med, c="#C4ABFF", ls="--")

axs[8].axvline(X_high_med, c="#008BB5", ls="--")
axs[7].axvline(X_low_med, c="#008BB5", ls="--")

axs[7].set_ylabel("Density")
axs[8].set_ylabel("Density")

#plt.xlim(0, 12.5)

axs[8].grid(False)
axs[8].set_facecolor("white")
axs[8].spines['bottom'].set_color('0')
axs[8].spines['left'].set_color('0')

axs[8].tick_params(bottom='on', left='on')

axs[7].grid(False)
axs[7].set_facecolor("white")
axs[7].spines['bottom'].set_color('0')
axs[7].spines['left'].set_color('0')














a=0
while a<9:
    axs[a].xaxis.set_tick_params(width=0.4, size=3, pad=2)
    axs[a].yaxis.set_tick_params(width=0.4, size=3, pad=2)
    a=a+1

axs[0].set_ylim([0, 11])

q=0
while q<9:
    axs[q].set_ylim([0, 15.5])
    q=q+1

plt.xlim(0, 1)

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/LUSC_and_TGCT_fig_1C_nonescaping_chrX_methylation_histograms_3_3_XIST_cutoff.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/LUSC_and_TGCT_fig_1C_nonescaping_chrX_methylation_histograms_3_3_XIST_cutoff.pdf", dpi=dpi_set, bbox_inches = 'tight')
    

