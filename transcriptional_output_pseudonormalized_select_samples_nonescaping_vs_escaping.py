#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 16:20:37 2021

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


def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

main_cutoff = 0

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(3.4,8.5*9/7)})
sns.set(font_scale=0.6)
plt.rcParams.update({'font.size': 10})
plt.rcParams['axes.linewidth'] = 0.3
dpi_set = 72 # change the output resolution

#full_class_list = ['LUSC', 'LUAD', 'SKCM', 'ESCA', 'KIRC', 'STAD', 'SKCM']

#ploidy_list = [3, 3, 3, 2, 2, 3, 2]





#CONVERT RATIO TO AUC RATIO





#NSGCTs:    

temp_class = 'TGCT-NS'

df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_averaged_TCGA_Xena_TPM.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]

if temp_class != 'TGCT-NS':

    df_XIST = df_XIST[df_XIST['Classification']==temp_class]
    
else:
    
    df_XIST = df_XIST[df_XIST['Secondary_Class']==temp_class]        

df_high = df_XIST[df_XIST['XIST_TPM']>=3]
df_low = df_XIST[df_XIST['XIST_TPM']<3]

h_samples = ['TCGA-2G-AAFZ-01']

l_samples = ['TCGA-2G-AAL5-01'] #TCGA-XE-AANI-01

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
        continue
        
pur_low = []

for a in df_low:
    try:
        pur_low.append(cn_dict[a])
    except KeyError:
        continue
    
print(statistics.median(pur_high))
print(statistics.median(pur_low))   

print("READING 1st RNA")

df_RNA_temp = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/gene_id_name_and_chr_all_biomart.csv")

print("READ 1st RNA")

gene_id = df_RNA_temp['Gene stable ID'].tolist()
chromosome = df_RNA_temp['Chromosome/scaffold name'].tolist()
gene_symbol = df_RNA_temp['Gene name'].tolist()

chr_dict = dict(zip(gene_id, chromosome))
symbol_dict = dict(zip(gene_id, gene_symbol))

print("READING 2nd RNA")

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_only_averaged_tcga_RSEM_gene_tpm.txt", sep="\t")

print("READ 2nd RNA")
    
gene_ids_temp_unp = df_RNA['sample'].tolist()

gene_ids_temp = []

for a in gene_ids_temp_unp:
    gene_ids_temp.append(a.split(".")[0])

new_chr = []
new_genes = []

for a in gene_ids_temp:
    try:
        new_chr.append(chr_dict[a])
    except KeyError:
        new_chr.append("UNKNOWN")
    try:
        new_genes.append(symbol_dict[a])
    except KeyError:
        new_genes.append("UNKNOWN")
    
df_RNA['chromosome'] = new_chr
df_RNA['gene_symbol'] = new_genes

df_RNA = df_RNA[df_RNA['chromosome']!="UNKNOWN"]
df_RNA = df_RNA[df_RNA['gene_symbol']!="UNKNOWN"]

print(df_RNA)

df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))
df_low = list(set(df_low) & set(df_RNA_temp_cols))

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

#df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

df_autosome = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

    
df_e_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/all_redone_e_pseudo_reference_XISTneg_males.csv")
df_ne_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/all_redone_ne_pseudo_reference_XISTneg_males.csv")

def generate_vals(df_inp, vals_ref, sample):
        
        gene_ids = df_inp['gene_symbol'].tolist()
        
        print(gene_ids)
        print(sample)
        print(df_inp)
        print(df_inp.columns.tolist())
        
        temp_vals = df_inp[sample[0]].tolist()
                                
        avg_ratio = []
        
        a=0
        while a<len(gene_ids):
            avg_ratio.append(math.log(temp_vals[a]/vals_ref[a], 2))
            a=a+1
        
        return avg_ratio

df_e_vals = df_e_ref['TGCT-NS'].tolist()
df_ne_vals = df_ne_ref['TGCT-NS'].tolist()
    
X_high = generate_vals(df_X, df_ne_vals, df_high)
other_high = generate_vals(df_autosome, df_e_vals, df_high)

X_low = generate_vals(df_X, df_ne_vals, df_low)
other_low = generate_vals(df_autosome, df_e_vals, df_low)

n_bins = np.arange(0, 15.01, 0.05)

fig, axs = plt.subplots(9, 1, sharex=True, tight_layout=True)

max_val = 4

plt.xlim(-1*max_val, max_val)

x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[4],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[4].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

X_high_med = np.trapz(y,x)




x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[4],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[4].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_high_med = np.trapz(y,x)





x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[3],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[3].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_low_med = np.trapz(y,x)



x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[3],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
X_low_med = np.trapz(y,x)
 
    
    
axs[3].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

#high_ratio = "ΔAUC$_{NE-E}$ = " + str("{:.3f}".format(X_high_med - other_high_med))
high_ratio = "Normalized to:\nXIST- NSGCT"
low_ratio = "Normalized to:\nXIST- NSGCT"

axs[4].text(0.02, 0.92, high_ratio, transform=axs[4].transAxes,
     fontsize=6, verticalalignment='top')


axs[3].text(0.02, 0.92, low_ratio, transform=axs[3].transAxes,
     fontsize=6, verticalalignment='top')

"""
axs[3].axvline(float(math.log(X_high_med+1, 2)), c="orange", ls="--")
axs[3].axvline(float(math.log(other_high_med+1, 2)), c="purple", ls="--")

axs[5].axvline(float(math.log(X_low_med+1, 2)), c="orange", ls="--")
axs[5].axvline(float(math.log(other_low_med+1, 2)), c="purple", ls="--")
"""

axs[4].axvline(statistics.median(other_high), c="purple", ls="--")
axs[3].axvline(statistics.median(other_low), c="purple", ls="--")

axs[4].axvline(statistics.median(X_high), c="orange", ls="--")
axs[3].axvline(statistics.median(X_low), c="orange", ls="--")

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

axs[3].tick_params(bottom='on', left='on')


print("MALE NSGCT - KS TEST")

print("Autosome distribution")

a, b = scipy.stats.ks_2samp(other_high, other_low)

print(a)
print(b)

print("chrX distribution")

a, b = scipy.stats.ks_2samp(X_high, X_low)

print(a)
print(b)

print("MALE NSGCT - CRAMER VON MISES TEST")

print("Autosome distribution")

a = scipy.stats.cramervonmises_2samp(other_high, other_low)

print(a)

print("chrX distribution")

a = scipy.stats.cramervonmises_2samp(X_high, X_low)

print(a)

print("MALE NSGCT - EPPS SINGLETON TEST")

print("Autosome distribution")

a, b = scipy.stats.epps_singleton_2samp(other_high, other_low)

print(a)
print(b)

print("chrX distribution")

a, b = scipy.stats.epps_singleton_2samp(X_high, X_low)

print(a)
print(b)

print("MALE NSGCT - ANDERSON DARLING TEST")

print("Autosome distribution")

a, b, c = scipy.stats.anderson_ksamp(np.array([other_high, other_low]))

print(a)
print(b)
print(c)

print("chrX distribution")

a, b, c = scipy.stats.anderson_ksamp(np.array([X_high, X_low]))

print(a)
print(b)
print(c)





mean_other_high = sum(other_high)/len(other_high)
std_other_high = np.std(other_high)
mean_other_low = sum(other_low)/len(other_low)
std_other_low = np.std(other_low)

mean_X_high = sum(X_high)/len(X_high)
std_X_high = np.std(X_high)
mean_X_low = sum(X_low)/len(X_low)
std_X_low = np.std(X_low)

diff_X = mean_X_high - mean_X_low

print(diff_X)

diff_auto = mean_other_high - mean_other_low

print(diff_auto)

std_X = math.sqrt(std_X_high**2/len(X_high) + std_X_low**2/len(X_low))

print(std_X)

std_auto = math.sqrt(std_other_high**2/len(other_high) + std_other_low**2/len(other_low))

print(std_auto)

z_X = diff_X/std_X

print(z_X)

z_auto = diff_auto/std_auto

print(z_auto)

p_X = scipy.stats.norm.sf(abs(z_X))*2

print("chrX:")
print(p_X)

p_auto = scipy.stats.norm.sf(abs(z_auto))*2

print("Autosome:")
print(p_auto)



full_class_list = ['LUSC'] #SKCM, LIHC, Pan-lung, Pan-GI (top 4)

#full_class_list = ['LUAD', 'LUSC'] #SKCM, LIHC, Pan-lung, Pan-GI (top 4)

df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_averaged_TCGA_Xena_TPM.csv")

for a in invalid_ids:
    df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]

df_XIST = df_XIST[df_XIST['Classification'].isin(full_class_list)]

df_high = df_XIST[df_XIST['XIST_TPM']>=3]

df_low = df_XIST[df_XIST['XIST_TPM']<3]

"""

df_high_lineages = df_high['Classification'].tolist()

n_SKCM = 0
n_LIHC = 0
n_LUSC = 0
n_LUAD = 0
n_STAD = 0
n_ESCA = 0

for a in df_high_lineages:
    if a == "SKCM":
        n_SKCM = n_SKCM+1
    elif a == "LIHC":
        n_LIHC = n_LIHC+1
    elif a == "LUSC":
        n_LUSC = n_LUSC+1
    elif a == "LUAD":
        n_LUAD = n_LUAD+1
    elif a == "STAD":
        n_STAD = n_STAD+1
    elif a == "ESCA":
        n_ESCA = n_ESCA+1

ratio = 20
n_SKCM = n_SKCM*ratio
n_LIHC = n_LIHC*ratio
n_LUSC = n_LUSC*ratio
n_LUAD = n_LUAD*ratio
n_STAD = n_STAD*ratio
n_ESCA = n_ESCA*ratio

df_low = df_low.sort_values(by=['XIST_TPM'], ascending=True)

df_low_all_samples = df_low['Barcode'].tolist()
df_low_all_lineages = df_low['Classification'].tolist()

dn_SKCM = 0
dn_LIHC = 0
dn_LUSC = 0
dn_LUAD = 0
dn_STAD = 0
dn_ESCA = 0

df_low_subset = []

a=0
while a<len(df_low_all_samples):
    if df_low_all_lineages[a] == "SKCM" and dn_SKCM<n_SKCM:
        df_low_subset.append(df_low_all_samples[a])
        dn_SKCM = dn_SKCM+1
    elif df_low_all_lineages[a] == "LIHC" and dn_LIHC<n_LIHC:
        df_low_subset.append(df_low_all_samples[a])
        dn_LIHC = dn_LIHC+1
    elif df_low_all_lineages[a] == "LUSC" and dn_LUSC<n_LUSC:
        df_low_subset.append(df_low_all_samples[a])
        dn_LUSC = dn_LUSC+1
    elif df_low_all_lineages[a] == "LUAD" and dn_LUAD<n_LUAD:
        df_low_subset.append(df_low_all_samples[a])
        dn_LUAD = dn_LUAD+1
    elif df_low_all_lineages[a] == "STAD" and dn_STAD<n_STAD:
        df_low_subset.append(df_low_all_samples[a])
        dn_STAD = dn_STAD+1
    elif df_low_all_lineages[a] == "ESCA" and dn_ESCA<n_ESCA:
        df_low_subset.append(df_low_all_samples[a])
        dn_ESCA = dn_ESCA+1
    a=a+1

"""

h_samples = ['TCGA-56-7822-01']

l_samples = ['TCGA-77-6843-01'] #TCGA-94-A5I4-01;       TCGA-77-A5G3-01, TCGA-85-8277-01, TCGA-18-3409-01


df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.processed.txt", sep="\t")

uq_samples = df2['array'].tolist()

chrX_cn = df2['purity'].tolist()

cn_dict = dict(zip(uq_samples, chrX_cn))

pur_high = []

#h_samples = df_high['Barcode'].tolist()

#l_samples = df_low['Barcode'].tolist()

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
        continue
        
pur_low = []

for a in df_low:
    try:
        pur_low.append(cn_dict[a])
    except KeyError:
        continue
    
print("hi")
    
print(statistics.median(pur_high))
print(statistics.median(pur_low))

df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))
df_low = list(set(df_low) & set(df_RNA_temp_cols))

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

#df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

df_autosome = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

df_e_vals = df_e_ref['LUSC'].tolist()
df_ne_vals = df_ne_ref['LUSC'].tolist()
    
X_high = generate_vals(df_X, df_ne_vals, df_high)
other_high = generate_vals(df_autosome, df_e_vals, df_high)

X_low = generate_vals(df_X, df_ne_vals, df_low)
other_low = generate_vals(df_autosome, df_e_vals, df_low)

x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[6],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[6].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_high_med = np.trapz(y,x)






x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[6],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[6].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

X_high_med = np.trapz(y,x)






x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
        
axs[5].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_low_med = np.trapz(y,x)




x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[5].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

X_low_med = np.trapz(y,x)


high_ratio = "Normalized to:\nXIST- M-LUSC"
low_ratio = "Normalized to:\nXIST- M-LUSC"

axs[6].text(0.02, 0.92, high_ratio, transform=axs[6].transAxes,
     fontsize=6, verticalalignment='top')

axs[5].text(0.02, 0.92, low_ratio, transform=axs[5].transAxes,
     fontsize=6, verticalalignment='top')
"""
axs[6].axvline(float(math.log(X_high_med+1, 2)), c="orange", ls="--")
axs[6].axvline(float(math.log(other_high_med+1, 2)), c="purple", ls="--")

axs[8].axvline(float(math.log(X_low_med+1, 2)), c="orange", ls="--")
axs[8].axvline(float(math.log(other_low_med+1, 2)), c="purple", ls="--")


axs[4].axvline(other_high_med, c="blue", ls="--")
axs[5].axvline(other_low_med, c="blue", ls="--")

axs[4].axvline(X_high_med, c="red", ls="--")
axs[5].axvline(X_low_med, c="red", ls="--")
"""

axs[6].axvline(statistics.median(other_high), c="purple", ls="--")
axs[5].axvline(statistics.median(other_low), c="purple", ls="--")

axs[6].axvline(statistics.median(X_high), c="orange", ls="--")
axs[5].axvline(statistics.median(X_low), c="orange", ls="--")

axs[6].set_ylabel("Density")
axs[6].set_xlabel("log$_2$(Normalized mRNA Expression)")
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


print("MALE PANCAN - KS TEST")

print("Autosome distribution")

a, b = scipy.stats.ks_2samp(other_high, other_low)

print(a)
print(b)

print("chrX distribution")

a, b = scipy.stats.ks_2samp(X_high, X_low)

print(a)
print(b)

print("MALE NSGCT - CRAMER VON MISES TEST")

print("Autosome distribution")

a = scipy.stats.cramervonmises_2samp(other_high, other_low)

print(a)

print("chrX distribution")

a = scipy.stats.cramervonmises_2samp(X_high, X_low)

print(a)

print("MALE NSGCT - EPPS SINGLETON TEST")

print("Autosome distribution")

a, b = scipy.stats.epps_singleton_2samp(other_high, other_low)

print(a)
print(b)

print("chrX distribution")

a, b = scipy.stats.epps_singleton_2samp(X_high, X_low)

print(a)
print(b)

print("MALE NSGCT - ANDERSON DARLING TEST")

print("Autosome distribution")

a, b, c = scipy.stats.anderson_ksamp(np.array([other_high, other_low]))

print(a)
print(b)
print(c)

print("chrX distribution")

a, b, c = scipy.stats.anderson_ksamp(np.array([X_high, X_low]))

print(a)
print(b)
print(c)


mean_other_high = sum(other_high)/len(other_high)
std_other_high = np.std(other_high)
mean_other_low = sum(other_low)/len(other_low)
std_other_low = np.std(other_low)

mean_X_high = sum(X_high)/len(X_high)
std_X_high = np.std(X_high)
mean_X_low = sum(X_low)/len(X_low)
std_X_low = np.std(X_low)

diff_X = mean_X_high - mean_X_low

print(diff_X)

diff_auto = mean_other_high - mean_other_low

print(diff_auto)

std_X = math.sqrt(std_X_high**2/len(X_high) + std_X_low**2/len(X_low))

print(std_X)

std_auto = math.sqrt(std_other_high**2/len(other_high) + std_other_low**2/len(other_low))

print(std_auto)

z_X = diff_X/std_X

print(z_X)

z_auto = diff_auto/std_auto

print(z_auto)

p_X = scipy.stats.norm.sf(abs(z_X))*2

print("chrX:")
print(p_X)

p_auto = scipy.stats.norm.sf(abs(z_auto))*2

print("Autosome:")
print(p_auto)



#SGCT

temp_class = 'TGCT-S'

df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_averaged_TCGA_Xena_TPM.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]

if temp_class != 'TGCT-S':

    df_XIST = df_XIST[df_XIST['Classification']==temp_class]
    
else:
    
    df_XIST = df_XIST[df_XIST['Secondary_Class']==temp_class]        

df_high = df_XIST[df_XIST['XIST_TPM']>=3]

h_samples = df_high['Barcode'].tolist()

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

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

#df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

df_autosome = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

df_e_vals = df_e_ref['TGCT-NS'].tolist()
df_ne_vals = df_ne_ref['TGCT-NS'].tolist()
    
X_high = generate_vals(df_X, df_ne_vals, df_high)
other_high = generate_vals(df_autosome, df_e_vals, df_high)

x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[2],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[2].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_high_med = np.trapz(y,x)




x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[2],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[2].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

X_high_med = np.trapz(y,x)




high_ratio = "Normalized to:\nXIST- NSGCT"

axs[2].text(0.02, 0.92, high_ratio, transform=axs[2].transAxes,
     fontsize=6, verticalalignment='top')
"""
axs[2].axvline(float(math.log(X_high_med+1, 2)), c="orange", ls="--")
axs[2].axvline(float(math.log(other_high_med+1, 2)), c="purple", ls="--")


axs[1].axvline(other_high_med, c="blue", ls="--")
axs[1].axvline(X_high_med, c="red", ls="--")
"""

axs[2].axvline(statistics.median(other_high), c="purple", ls="--")
axs[2].axvline(statistics.median(X_high), c="orange", ls="--")

axs[2].set_ylabel("Density")

#plt.xlim(0, 12.5)

axs[2].grid(False)
axs[2].set_facecolor("white")
axs[2].spines['bottom'].set_color('0')
axs[2].spines['left'].set_color('0')

axs[2].tick_params(bottom='on', left='on')












"""
temp_class = 'TGCT-NS'

df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_averaged_TCGA_Xena_TPM.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]

if temp_class != 'TGCT-NS':

    df_XIST = df_XIST[df_XIST['Classification']==temp_class]
    
else:
    
    df_XIST = df_XIST[df_XIST['Secondary_Class']==temp_class]        

df_high = df_XIST[df_XIST['XIST_TPM']>=3]

h_samples = df_high['Barcode'].tolist()

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

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

#df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

df_autosome = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

df_e_vals = df_e_ref['TGCT-NS'].tolist()
df_ne_vals = df_ne_ref['TGCT-NS'].tolist()
    
X_high = generate_vals(df_X, df_ne_vals, df_high)
other_high = generate_vals(df_autosome, df_e_vals, df_high)

x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 2:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[5].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_high_med = np.trapz(y,x)




x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 2:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[5].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

X_high_med = np.trapz(y,x)



high_ratio = "ΔAUC$_{NE-E}$ = " + str("{:.3f}".format(X_high_med - other_high_med))

"""
#axs[4].axvline(float(math.log(X_high_med+1, 2)), c="orange", ls="--")
#axs[4].axvline(float(math.log(other_high_med+1, 2)), c="purple", ls="--")
"""

axs[5].text(0.73, 0.98, high_ratio, transform=axs[5].transAxes,
     fontsize=6, verticalalignment='top')

"""
#axs[1].axvline(other_high_med, c="blue", ls="--")
#axs[1].axvline(X_high_med, c="red", ls="--")
"""

axs[5].set_ylabel("Density")


axs[5].grid(False)
axs[5].set_facecolor("white")
axs[5].spines['bottom'].set_color('0')
axs[5].spines['left'].set_color('0')

axs[5].tick_params(bottom='on', left='on')







temp_class = 'LUSC'

df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_averaged_TCGA_Xena_TPM.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]

if temp_class != 'TGCT-S':

    df_XIST = df_XIST[df_XIST['Classification']==temp_class]
    
else:
    
    df_XIST = df_XIST[df_XIST['Secondary_Class']==temp_class]        

df_high = df_XIST[df_XIST['XIST_TPM']>=3]

h_samples = df_high['Barcode'].tolist()

df_high = ['TCGA-37-4130-01'] #BAD - TCGA-18-3410-01, TCGA-60-2703-01

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

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

#df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

df_autosome = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

df_e_vals = df_e_ref['LUSC'].tolist()
df_ne_vals = df_ne_ref['LUSC'].tolist()
    
X_high = generate_vals(df_X, df_ne_vals, df_high)
other_high = generate_vals(df_autosome, df_e_vals, df_high)

x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 2:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[8].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_high_med = np.trapz(y,x)





x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 2:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[8].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

X_high_med = np.trapz(y,x)



high_ratio = "ΔAUC$_{NE-E}$ = " + str("{:.3f}".format(X_high_med - other_high_med))

"""
#axs[7].axvline(float(math.log(X_high_med+1, 2)), c="orange", ls="--")
#axs[7].axvline(float(math.log(other_high_med+1, 2)), c="purple", ls="--")
"""

axs[8].text(0.73, 0.98, high_ratio, transform=axs[8].transAxes,
     fontsize=6, verticalalignment='top')

"""
#axs[1].axvline(other_high_med, c="blue", ls="--")
#axs[1].axvline(X_high_med, c="red", ls="--")
"""

axs[8].set_ylabel("Density")


axs[8].grid(False)
axs[8].set_facecolor("white")
axs[8].spines['bottom'].set_color('0')
axs[8].spines['left'].set_color('0')

axs[8].tick_params(bottom='on', left='on')
"""





df_RNA_symbols = df_RNA['gene_symbol'].tolist()
df_RNA_gene_ids = df_RNA['sample'].tolist()
df_RNA_chromosomes = df_RNA['chromosome'].tolist()

symbol_dict = dict(zip(df_RNA_gene_ids, df_RNA_symbols))
chr_dict = dict(zip(df_RNA_gene_ids, df_RNA_chromosomes))

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_normals_only_averaged_tcga_RSEM_gene_tpm.txt", sep="\t")

for a in invalid_ids:
    df_RNA = df_RNA[~(df_RNA['sample'].str.contains(a))]

df_RNA_temp_cols = df_RNA.columns.tolist()

del df_RNA_temp_cols[0]

df_low = ['TCGA-22-5482-11'] #TCGA-22-5482-11 > TCGA-22-5489-11, TCGA-56-8083-11

df_RNA_old_symbols = df_RNA['sample'].tolist()

df_RNA_new_symbols = []
df_RNA_new_chr = []

for a in df_RNA_old_symbols:
    try:
        df_RNA_new_symbols.append(symbol_dict[a])
        df_RNA_new_chr.append(chr_dict[a])
    except KeyError:
        df_RNA_new_symbols.append("UNKNOWN")
        df_RNA_new_chr.append("UNKNOWN")        

df_RNA['gene_symbol'] = df_RNA_new_symbols
df_RNA['chromosome'] = df_RNA_new_chr

df_RNA['gene_id'] = df_RNA['sample']
del df_RNA['sample']

df_RNA = df_RNA[df_RNA['gene_symbol']!="UNKNOWN"]
df_RNA = df_RNA[df_RNA['chromosome']!="UNKNOWN"]

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

#df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

df_autosome = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

df_e_vals = df_e_ref['XIST_neg_male_normals'].tolist()
df_ne_vals = df_ne_ref['XIST_neg_male_normals'].tolist()

X_low = generate_vals(df_X, df_ne_vals, df_low)
other_low = generate_vals(df_autosome, df_e_vals, df_low)

x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[1],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[1].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_low_med = np.trapz(y,x)






x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[1],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[1].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

X_low_med = np.trapz(y,x)







low_ratio = " Normalized to:\nXIST- M-Normal"

"""
axs[1].axvline(float(math.log(X_low_med+1, 2)), c="orange", ls="--")
axs[1].axvline(float(math.log(other_low_med+1, 2)), c="purple", ls="--")
"""

axs[1].text(0.02, 0.92, low_ratio, transform=axs[1].transAxes,
     fontsize=6, verticalalignment='top')

"""
axs[0].axvline(other_low_med, c="blue", ls="--")
axs[0].axvline(X_low_med, c="red", ls="--")
"""

axs[1].axvline(statistics.median(other_low), c="purple", ls="--")
axs[1].axvline(statistics.median(X_low), c="orange", ls="--")

axs[1].set_ylabel("Density")

axs[1].grid(False)
axs[1].set_facecolor("white")
axs[1].spines['bottom'].set_color('0')
axs[1].spines['left'].set_color('0')

axs[1].tick_params(bottom='on', left='on')











df_RNA_symbols = df_RNA['gene_symbol'].tolist()
df_RNA_gene_ids = df_RNA['gene_id'].tolist()
df_RNA_chromosomes = df_RNA['chromosome'].tolist()

symbol_dict = dict(zip(df_RNA_gene_ids, df_RNA_symbols))
chr_dict = dict(zip(df_RNA_gene_ids, df_RNA_chromosomes))

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_normals_only_averaged_tcga_RSEM_gene_tpm.txt", sep="\t")

for a in invalid_ids:
    df_RNA = df_RNA[~(df_RNA['sample'].str.contains(a))]

df_RNA_temp_cols = df_RNA.columns.tolist()

del df_RNA_temp_cols[0]

df_low = ['TCGA-E9-A1NA-11'] #TCGA-BH-A209-11 > Good TCGA-BH-A18M-11, - TCGA-BH-A0DO-11; TCGA-55-6981-11 > TCGA-44-5645-11

df_RNA_old_symbols = df_RNA['sample'].tolist()

df_RNA_new_symbols = []
df_RNA_new_chr = []

for a in df_RNA_old_symbols:
    try:
        df_RNA_new_symbols.append(symbol_dict[a])
        df_RNA_new_chr.append(chr_dict[a])
    except KeyError:
        df_RNA_new_symbols.append("UNKNOWN")
        df_RNA_new_chr.append("UNKNOWN")        

df_RNA['gene_symbol'] = df_RNA_new_symbols
df_RNA['chromosome'] = df_RNA_new_chr

df_RNA['gene_id'] = df_RNA['sample']
del df_RNA['sample']

df_RNA = df_RNA[df_RNA['gene_symbol']!="UNKNOWN"]
df_RNA = df_RNA[df_RNA['chromosome']!="UNKNOWN"]

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

#df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

df_autosome = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

df_e_vals = df_e_ref['XIST_pos_female_normals'].tolist()
df_ne_vals = df_ne_ref['XIST_pos_female_normals'].tolist()

X_low = generate_vals(df_X, df_ne_vals, df_low)
other_low = generate_vals(df_autosome, df_e_vals, df_low)

x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[0],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[0].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_low_med = np.trapz(y,x)





x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[0],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[0].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

X_low_med = np.trapz(y,x)




low_ratio = " Normalized to:\nXIST+ F-Normal"

"""
axs[0].axvline(float(math.log(X_low_med+1, 2)), c="orange", ls="--")
axs[0].axvline(float(math.log(other_low_med+1, 2)), c="purple", ls="--")
"""

axs[0].text(0.02, 0.92, low_ratio, transform=axs[0].transAxes,
     fontsize=6, verticalalignment='top')

"""
axs[0].axvline(other_low_med, c="blue", ls="--")
axs[0].axvline(X_low_med, c="red", ls="--")
"""

axs[0].axvline(statistics.median(other_low), c="purple", ls="--")
axs[0].axvline(statistics.median(X_low), c="orange", ls="--")

axs[0].set_ylabel("Density")

axs[0].grid(False)
axs[0].set_facecolor("white")
axs[0].spines['bottom'].set_color('0')
axs[0].spines['left'].set_color('0')

axs[0].tick_params(bottom='on', left='on')

nonescape_patch = mpatches.Patch(color='orange', label='chrX Non-Escapees')
escape_patch = mpatches.Patch(color='purple', label='chrX Escapees')

leg = axs[0].legend(handles=[nonescape_patch, escape_patch], facecolor='white', fontsize=6, bbox_to_anchor=(0.615, 0.945))

frame = leg.get_frame()
frame.set_edgecolor("black")
frame.set_linewidth(0.5)

axs[3].tick_params(bottom='on', left='on')









temp_class = 'OV'

df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_tumors_averaged_TCGA_Xena_TPM.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]

if temp_class != 'TGCT-NS':

    df_XIST = df_XIST[df_XIST['Classification']==temp_class]
    
else:
    
    df_XIST = df_XIST[df_XIST['Secondary_Class']==temp_class]        

df_high = df_XIST[df_XIST['XIST_TPM']>=3]
df_low = df_XIST[df_XIST['XIST_TPM']<3]

h_samples = ['TCGA-23-2078-01']

l_samples = ['TCGA-29-1702-01'] #TCGA-XE-AANI-01

df_high = h_samples

df_low = l_samples

print("READING 1st RNA")

df_RNA_temp = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/gene_id_name_and_chr_all_biomart.csv")

print("READ 1st RNA")

gene_id = df_RNA_temp['Gene stable ID'].tolist()
chromosome = df_RNA_temp['Chromosome/scaffold name'].tolist()
gene_symbol = df_RNA_temp['Gene name'].tolist()

chr_dict = dict(zip(gene_id, chromosome))
symbol_dict = dict(zip(gene_id, gene_symbol))

print("READING 2nd RNA")

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_tumors_only_averaged_tcga_RSEM_gene_tpm.txt", sep="\t")

print("READ 2nd RNA")
    
gene_ids_temp_unp = df_RNA['sample'].tolist()

gene_ids_temp = []

for a in gene_ids_temp_unp:
    gene_ids_temp.append(a.split(".")[0])

new_chr = []
new_genes = []

for a in gene_ids_temp:
    try:
        new_chr.append(chr_dict[a])
    except KeyError:
        new_chr.append("UNKNOWN")
    try:
        new_genes.append(symbol_dict[a])
    except KeyError:
        new_genes.append("UNKNOWN")
    
df_RNA['chromosome'] = new_chr
df_RNA['gene_symbol'] = new_genes

df_RNA = df_RNA[df_RNA['chromosome']!="UNKNOWN"]
df_RNA = df_RNA[df_RNA['gene_symbol']!="UNKNOWN"]

print(df_RNA)

df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))
df_low = list(set(df_low) & set(df_RNA_temp_cols))

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

#df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

df_autosome = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

    
df_e_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/e_pseudo_reference_XISTpos_females.csv")
df_ne_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/ne_pseudo_reference_XISTpos_females.csv")

def generate_vals(df_inp, vals_ref, sample):
        
        gene_ids = df_inp['gene_symbol'].tolist()
        
        print(gene_ids)
        print(sample)
        print(df_inp)
        print(df_inp.columns.tolist())
        
        temp_vals = df_inp[sample[0]].tolist()
                                
        avg_ratio = []
        
        a=0
        while a<len(gene_ids):
            avg_ratio.append(math.log(temp_vals[a]/vals_ref[a], 2))
            a=a+1
        
        return avg_ratio

df_e_vals = df_e_ref['OV'].tolist()
df_ne_vals = df_ne_ref['OV'].tolist()
    
X_high = generate_vals(df_X, df_ne_vals, df_high)
other_high = generate_vals(df_autosome, df_e_vals, df_high)

X_low = generate_vals(df_X, df_ne_vals, df_low)
other_low = generate_vals(df_autosome, df_e_vals, df_low)

n_bins = np.arange(0, 15.01, 0.05)

max_val = 4

plt.xlim(-1*max_val, max_val)

x_unp, y_unp = sns.distplot(X_high, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[8].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

X_high_med = np.trapz(y,x)




x_unp, y_unp = sns.distplot(other_high, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[8].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_high_med = np.trapz(y,x)





x_unp, y_unp = sns.distplot(other_low, hist=False, kde=True, ax=axs[7],
             bins=n_bins, color = 'purple', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[0].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
axs[7].fill_between(x, 0,  y, facecolor='purple', alpha=0.2)

other_low_med = np.trapz(y,x)



x_unp, y_unp = sns.distplot(X_low, hist=False, kde=True, ax=axs[7],
             bins=n_bins, color = 'orange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5, 'clip': (-1*max_val, max_val)}).get_lines()[1].get_data()

x = []
y = []
a=0
while a<len(x_unp):
    if x_unp[a] < 5:
        x.append(x_unp[a])
        y.append(y_unp[a])
    a=a+1
    
X_low_med = np.trapz(y,x)
 
    
    
axs[7].fill_between(x, 0,  y, facecolor='orange', alpha=0.2)

#high_ratio = "ΔAUC$_{NE-E}$ = " + str("{:.3f}".format(X_high_med - other_high_med))
high_ratio = "Normalized to:\nXIST+ OV"
low_ratio = "Normalized to:\nXIST+ OV"

axs[8].text(0.02, 0.92, high_ratio, transform=axs[8].transAxes,
     fontsize=6, verticalalignment='top')


axs[7].text(0.02, 0.92, low_ratio, transform=axs[7].transAxes,
     fontsize=6, verticalalignment='top')

"""
axs[3].axvline(float(math.log(X_high_med+1, 2)), c="orange", ls="--")
axs[3].axvline(float(math.log(other_high_med+1, 2)), c="purple", ls="--")

axs[5].axvline(float(math.log(X_low_med+1, 2)), c="orange", ls="--")
axs[5].axvline(float(math.log(other_low_med+1, 2)), c="purple", ls="--")
"""

axs[7].axvline(statistics.median(other_high), c="purple", ls="--")
axs[8].axvline(statistics.median(other_low), c="purple", ls="--")

axs[7].axvline(statistics.median(X_high), c="orange", ls="--")
axs[8].axvline(statistics.median(X_low), c="orange", ls="--")

axs[8].set_ylabel("Density")
axs[7].set_ylabel("Density")

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

axs[7].tick_params(bottom='on', left='on')






a=0
while a<9:
    axs[a].xaxis.set_tick_params(width=0.4, size=3, pad=2)
    axs[a].yaxis.set_tick_params(width=0.4, size=3, pad=2)
    a=a+1



"""
#Female normals

df_RNA_symbols = df_RNA['gene_symbol'].tolist()
df_RNA_gene_ids = df_RNA['gene_id'].tolist()
df_RNA_chromosomes = df_RNA['chromosome'].tolist()

symbol_dict = dict(zip(df_RNA_gene_ids, df_RNA_symbols))
chr_dict = dict(zip(df_RNA_gene_ids, df_RNA_chromosomes))

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_normals_only_averaged_tcga_RSEM_gene_tpm.txt", sep="\t")

for a in invalid_ids:
    df_RNA = df_RNA[~(df_RNA['sample'].str.contains(a))]

df_RNA_temp_cols = df_RNA.columns.tolist()

del df_RNA_temp_cols[0]

df_low = df_RNA_temp_cols

df_RNA_old_symbols = df_RNA['sample'].tolist()

df_RNA_new_symbols = []
df_RNA_new_chr = []

for a in df_RNA_old_symbols:
    try:
        df_RNA_new_symbols.append(symbol_dict[a])
        df_RNA_new_chr.append(chr_dict[a])
    except KeyError:
        df_RNA_new_symbols.append("UNKNOWN")
        df_RNA_new_chr.append("UNKNOWN")        

df_RNA['gene_symbol'] = df_RNA_new_symbols
df_RNA['chromosome'] = df_RNA_new_chr

df_RNA['gene_id'] = df_RNA['sample']
del df_RNA['sample']

df_RNA = df_RNA[df_RNA['gene_symbol']!="UNKNOWN"]
df_RNA = df_RNA[df_RNA['chromosome']!="UNKNOWN"]

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, low_class):
    low_list = []
    
    gene_ids = df_inp['gene_id'].tolist()
    
    low_df = df_inp[low_class]
    
    low_vals = low_df.values.tolist()
    
    low_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        curr_low = [x for x in low_vals[a] if x==x]

        avg_low = sum(curr_low)/len(curr_low)      
    
        if avg_low >= main_cutoff:     
            
            low_total_TPM = low_total_TPM + avg_low
            
            low_list.append(float(math.log(avg_low+1, 2)))
        
        a=a+1
        
    return low_list, low_total_TPM

X_low,X_total_low  = DGE(df_X, df_low)

X_low_med = statistics.median(X_low)

other_low,other_total_low = DGE(df_autosome, df_low)

sns.distplot(other_low, hist=False, kde=True, ax=axs[1],
             bins=n_bins, color = 'blue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

sns.distplot(X_low, hist=False, kde=True, ax=axs[1],
             bins=n_bins, color = 'red', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

other_low_med = statistics.median(other_low)

low_ratio = "Ratio = " + str("{:.2f}".format((2**X_low_med-1)/(2**other_low_med-1)))

female_autosome_dist = other_low
female_chrX_dist = X_low

axs[1].text(0.85, 0.98, low_ratio, transform=axs[1].transAxes,
     fontsize=18, verticalalignment='top')

axs[1].axvline(other_low_med, c="blue", ls="--")
axs[1].axvline(X_low_med, c="red", ls="--")
axs[1].set_ylabel("Density (XIST+ F-Normals)")

axs[1].grid(False)
axs[1].set_facecolor("white")
axs[1].spines['bottom'].set_color('0')
axs[1].spines['left'].set_color('0')

axs[1].tick_params(bottom='on', left='on')




print("NORMALS - KS TEST")

print("Autosome distribution")

a, b = scipy.stats.ks_2samp(male_autosome_dist, female_autosome_dist)

print(a)
print(b)

print("chrX distribution")

a, b = scipy.stats.ks_2samp(male_chrX_dist, female_chrX_dist)

print(a)
print(b)


axs[1].text(0.85, 0.98, "Autosome: p = " + str(), transform=axs[1].transAxes,
     fontsize=18, verticalalignment='top')



#Female PanCan

full_class_list = ['UVM', 'OV', 'SKCM', 'ACC'] #UVM, OV, SKCM, ACC

df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_tumors_averaged_TCGA_Xena_TPM.csv")

for a in invalid_ids:
    df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]

df_XIST = df_XIST[df_XIST['Classification'].isin(full_class_list)]

df_high = df_XIST[df_XIST['XIST_TPM']>=3]
df_low = df_XIST[df_XIST['XIST_TPM']<3]

h_samples = df_high['Barcode'].tolist()

l_samples = df_low['Barcode'].tolist()

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

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/averaged_female_all_chr_segment_mean_normalized_expression_rsem_tpm.txt", sep="\t")

df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))
df_low = list(set(df_low) & set(df_RNA_temp_cols))

df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

def DGE(df_inp, high_class, low_class):
    high_list = []
    low_list = []
    
    gene_ids = df_inp['gene_id'].tolist()
    
    high_df = df_inp[high_class]
    low_df = df_inp[low_class]
    
    
    
    high_vals = high_df.values.tolist()
    low_vals = low_df.values.tolist()
    
    high_total_TPM = 0
    low_total_TPM = 0
    
    a=0
    while a<len(gene_ids):
        
        curr_high = [x for x in high_vals[a] if x==x]
        curr_low = [x for x in low_vals[a] if x==x]
        
        total_list = curr_high + curr_low
        
        avg_total = sum(total_list)/len(total_list)
        
        avg_high = sum(curr_high)/len(curr_high)
        avg_low = sum(curr_low)/len(curr_low)      
    
        if avg_high >= main_cutoff and avg_low >= main_cutoff:     
            
            high_total_TPM = high_total_TPM + avg_high
            low_total_TPM = low_total_TPM + avg_low
            
            high_list.append(float(math.log(avg_high+1, 2)))
            low_list.append(float(math.log(avg_low+1, 2)))
        
        a=a+1
        
    return high_list, low_list, high_total_TPM, low_total_TPM

print('\nChrX')

X_high,X_low,X_total_high,X_total_low  = DGE(df_X, df_high, df_low)

print("Ratio of XIST+ and XIST- total TPM: " + str(X_total_high/X_total_low))

X_high_med = statistics.median(X_high)

X_low_med = statistics.median(X_low)

print("XIST+ Median: " + str(X_high_med) + "; equal to: " + str(2**X_high_med-1))
print("XIST- Median: " + str(X_low_med) + "; equal to: " + str(2**X_low_med-1))
print("Delta Median: " + str(X_low_med-X_high_med) + "; equal to: " + str(2**(X_low_med-X_high_med)))

print('\nAutosomes')

other_high,other_low,other_total_high,other_total_low = DGE(df_autosome, df_high, df_low)

print("Ratio of XIST+ and XIST- total TPM: " + str(other_total_high/other_total_low))

sns.distplot(other_high, hist=False, kde=True, ax=axs[7],
             bins=n_bins, color = 'blue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

sns.distplot(X_high, hist=False, kde=True, ax=axs[7],
             bins=n_bins, color = 'red', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

sns.distplot(other_low, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = 'blue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

sns.distplot(X_low, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = 'red', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

other_high_med = statistics.median(other_high)

other_low_med = statistics.median(other_low)

print("XIST+ Median: " + str(other_high_med) + "; equal to: " + str(2**other_high_med-1))
print("XIST- Median: " + str(other_low_med) + "; equal to: " + str(2**other_low_med-1))
print("Delta Median: " + str(other_low_med-other_high_med) + "; equal to: " + str(2**(other_low_med-other_high_med)))

high_ratio = "Ratio = " + str("{:.2f}".format((2**X_high_med-1)/(2**other_high_med-1)))

axs[7].text(0.85, 0.98, high_ratio, transform=axs[7].transAxes,
     fontsize=18, verticalalignment='top')

low_ratio = "Ratio = " + str("{:.2f}".format((2**X_low_med-1)/(2**other_low_med-1)))

axs[8].text(0.85, 0.98, low_ratio, transform=axs[8].transAxes,
     fontsize=18, verticalalignment='top')

axs[7].axvline(other_high_med, c="blue", ls="--")
axs[8].axvline(other_low_med, c="blue", ls="--")

axs[7].axvline(X_high_med, c="red", ls="--")
axs[8].axvline(X_low_med, c="red", ls="--")

axs[7].set_ylabel("Density (XIST+ F-PanCan)")
axs[8].set_xlabel("Average mRNA Expression (log$_2$(TPM+1))")
axs[8].set_ylabel("Density (XIST- F-PanCan)")

#plt.xlim(0, 12.5)

axs[7].grid(False)
axs[7].set_facecolor("white")
axs[7].spines['bottom'].set_color('0')
axs[7].spines['left'].set_color('0')

axs[7].tick_params(bottom='on', left='on')

axs[8].grid(False)
axs[8].set_facecolor("white")
axs[8].spines['bottom'].set_color('0')
axs[8].spines['left'].set_color('0')
axs[8].tick_params(bottom='on', left='on')
"""

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/LUSC_and_TGCT_fig_1C_nonescaping_chrX_transcriptional_output_histograms_3_3_XIST_cutoff.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/LUSC_and_TGCT_fig_1C_nonescaping_chrX_transcriptional_output_histograms_3_3_XIST_cutoff.pdf", dpi=dpi_set, bbox_inches = 'tight')
    