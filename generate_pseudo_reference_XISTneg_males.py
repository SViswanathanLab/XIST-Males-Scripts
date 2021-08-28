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

cutoff_value = 0

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(10,10)})
sns.set(font_scale=1)
plt.rcParams.update({'font.size': 14})
plt.rcParams['axes.linewidth'] = 0.3
dpi_set = 72 # change the output resolution

ne_values_list = []
ne_ids_list = []

e_values_list = []
e_ids_list = []

autosome_values_list = []
autosome_ids_list = []
    
df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]
    
df_barcode = df_XIST['Barcode'].tolist()

end_barcode = []

for a in df_barcode:
    end_barcode.append(int(split_advanced(a, "-", 3)[1]))

df_XIST['End_Val'] = end_barcode

df_XIST = df_XIST[df_XIST['End_Val']<10]
df_XIST = df_XIST[df_XIST['Classification'] != "UNKNOWN"]
df_XIST = df_XIST[df_XIST['Gender'] == "MALE"]
    
df_second = df_XIST['Secondary_Class'].tolist()
df_primary = df_XIST['Classification'].tolist()

df_new = []
a=0
while a<len(df_primary):
    if df_primary[a] == "TGCT":
        if df_second[a] == df_second[a]:
            df_new.append(df_second[a])
        else:
            df_new.append("TGCT-U")     
    else:
        df_new.append(df_primary[a])      
    a=a+1
    
df_XIST['Combined'] = df_new

#df_XIST = df_XIST[df_XIST['XIST_TPM']>=3]

full_class_list = list(set(df_XIST['Combined'].tolist()))

#full_class_list.remove("TGCT-S")

print("READING 1st RNA")

df_RNA_temp = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/gene_id_name_and_chr_all_biomart.csv")

print("READ 1st RNA")

gene_id = df_RNA_temp['Gene stable ID'].tolist()
chromosome = df_RNA_temp['Chromosome/scaffold name'].tolist()
gene_symbol = df_RNA_temp['Gene name'].tolist()

chr_dict = dict(zip(gene_id, chromosome))
symbol_dict = dict(zip(gene_id, gene_symbol))

print("READING 2nd RNA")

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/tcga_RSEM_gene_tpm", sep="\t")

print("READ 2nd RNA")

numeric_cols = [col for col in df_RNA if df_RNA[col].dtype.kind != 'O']

print(df_RNA)

df_RNA[numeric_cols] = df_RNA[numeric_cols].apply(lambda x: 2**x)

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

www=0
while www<len(full_class_list):
    
    df_XIST = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")
    
    invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']
    
    for a in invalid_ids:
        df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]
        
    df_barcode = df_XIST['Barcode'].tolist()
    
    end_barcode = []
    
    for a in df_barcode:
        end_barcode.append(int(split_advanced(a, "-", 3)[1]))
    
    df_XIST['End_Val'] = end_barcode
    
    df_XIST = df_XIST[df_XIST['End_Val']<10]
    df_XIST = df_XIST[df_XIST['Classification'] != "UNKNOWN"]
    df_XIST = df_XIST[df_XIST['Gender'] == "MALE"]
        
    df_second = df_XIST['Secondary_Class'].tolist()
    df_primary = df_XIST['Classification'].tolist()
    
    df_new = []
    a=0
    while a<len(df_primary):
        if df_primary[a] == "TGCT":
            if df_second[a] == df_second[a]:
                df_new.append(df_second[a])
            else:
                df_new.append("TGCT-U")     
        else:
            df_new.append(df_primary[a])      
        a=a+1
        
    df_XIST['Combined'] = df_new

    temp_class = full_class_list[www]
    
    df_XIST = df_XIST[df_XIST['Combined']==temp_class] 
    
    if full_class_list[www] != "TGCT-S":

        df_low = df_XIST[df_XIST['XIST_TPM']<3]

    else:

        df_low = df_XIST[df_XIST['XIST_TPM']>=3]
    
    l_samples = df_low['Barcode'].tolist()
        
    df_low = l_samples
 
    
    df_low = list(set(df_low) & set(df_RNA_temp_cols))
    
    df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")
    
    inactivated_genes = df_gene_class['Inactivated'].tolist()
    escaping_genes = df_gene_class['Escaping'].tolist()
    
    all_chrX_valid_genes = inactivated_genes + escaping_genes
        
    df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]
        
    df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

    df_escape = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

    def average_vals(df_inp, cols):
        
        gene_ids = df_inp['sample'].tolist()
        
        temp_df = df_inp[cols]
                        
        temp_vals = temp_df.values.tolist()
        
        avg_vals = []
        
        a=0
        while a<len(gene_ids):
            avg_vals.append(sum(temp_vals[a])/len(temp_vals[a]))
            a=a+1
        
        return gene_ids, avg_vals
    
    ne_ids, ne_avg = average_vals(df_X, df_low)
    autosome_ids, autosome_avg = average_vals(df_autosome, df_low)
    e_ids, e_avg = average_vals(df_escape, df_low)
    
    ne_values_list.append(ne_avg)
    ne_ids_list.append(ne_ids)
    
    e_values_list.append(e_avg)
    e_ids_list.append(e_ids)
    
    autosome_values_list.append(autosome_avg)
    autosome_ids_list.append(autosome_ids)
    
    www=www+1
    print(www)
    
ne_df = pd.DataFrame(ne_values_list).T
ne_df.columns = full_class_list
ne_df.index = ne_ids_list[0]

autosome_df = pd.DataFrame(autosome_values_list).T
autosome_df.columns = full_class_list
autosome_df.index = autosome_ids_list[0]

e_df = pd.DataFrame(e_values_list).T
e_df.columns = full_class_list
e_df.index = e_ids_list[0]

a=0
while a<len(ne_ids_list)-1:
    if ne_ids_list[a] != ne_ids_list[a+1]:
        print("ERROR")
    if autosome_ids_list[a] != autosome_ids_list[a+1]:
        print("ERROR")
    if e_ids_list[a] != e_ids_list[a+1]:
        print("ERROR")
    a=a+1

ne_df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_averaging_ne_pseudo_reference_XISTneg_males.csv")
autosome_df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_averaging_autosome_pseudo_reference_XISTneg_males.csv")
e_df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_averaging_e_pseudo_reference_XISTneg_males.csv")
