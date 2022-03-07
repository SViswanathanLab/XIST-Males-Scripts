#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 22:25:37 2021

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

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

df3 = pd.read_csv("../VAF_Analysis/Other_Input/TCGA_mastercalls.abs_tables_JSedit.fixed.processed.txt", sep="\t")

uq_samples = df3['array'].tolist()

chrX_cn = df3['purity'].tolist()

cn_dict = dict(zip(uq_samples, chrX_cn))

df_ploidy = pd.read_csv("../VAF_Analysis/Other_Input/median_autosome_copies_TCGA_median_ploidy.csv")

ploidy = df_ploidy['median_of_rounded_avg_modal_copies_per_chr'].tolist()

sample_test = df_ploidy['sample'].tolist()

ploidy_dict = dict(zip(sample_test, ploidy))

df2 = pd.read_csv("../VAF_Analysis/Other_Input/chrX_muts_REVISED_8_28_females_unaveraged.txt", sep="\t")

df2 = df2[df2['Variant_Classification']!="Intron"]

#df2 = df2[df2['Variant_Classification'].isin(['Missense_Mutation', 'Silent', "3'UTR", "5'UTR", "Splice_Site"])]

valid_samples = list(set(df2['Truncated_Barcodes'].tolist()))

df_gene_class = pd.read_excel("../VAF_Analysis/Other_Input/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()



df = pd.read_csv("../VAF_Analysis/Other_Input/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df = df[~(df['Barcode'].str.contains(a))]

df_barcode = df['Barcode'].tolist()

end_barcode = []

for a in df_barcode:
    end_barcode.append(int(split_advanced(a, "-", 3)[1]))

df['End_Val'] = end_barcode

df = df[df['End_Val']<10]
df = df[df['Classification'] != "UNKNOWN"]
df = df[df['Gender'] != "UNKNOWN"]
df = df[df['Gender'] == "FEMALE"]

df_barcodes = df['Barcode'].tolist()
df_XIST = df['XIST_TPM'].tolist()
df_class = df['Classification'].tolist()
df_second_class = df['Secondary_Class'].tolist()

XIST_dict = dict(zip(df_barcodes, df_XIST))
classification_dict = dict(zip(df_barcodes, df_class))
second_classification_dict = dict(zip(df_barcodes, df_second_class))

print("READING 1st RNA")
    
df_RNA_temp = pd.read_csv("../VAF_Analysis/Other_Input/gene_id_name_and_chr_all_biomart.csv")

print("READ 1st RNA")

gene_id = df_RNA_temp['Gene stable ID'].tolist()
chromosome = df_RNA_temp['Chromosome/scaffold name'].tolist()
gene_symbol = df_RNA_temp['Gene name'].tolist()

chr_dict = dict(zip(gene_id, chromosome))
symbol_dict = dict(zip(gene_id, gene_symbol))

print("READING 2nd RNA")

df_RNA = pd.read_csv("../VAF_Analysis/Other_Input/female_tumors_only_averaged_tcga_RSEM_gene_tpm.txt", sep="\t")

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

df_RNA_cols = df_RNA.columns.tolist()

valid_samples = list(set(df_RNA_cols) & set(valid_samples))

df_X = df_RNA

valid_samples.append('gene_symbol')

df_X = df_X[valid_samples]

print(df_X)

df_X_genes = list(set(df_X['gene_symbol'].tolist()))

#df2 = df2[df2['Hugo_Symbol'].isin(df_X_genes)]
df2 = df2[df2['Truncated_Barcodes'].isin(valid_samples)]

print(df2)

all_hugos = df2['Hugo_Symbol'].tolist()
all_barcodes = df2['Truncated_Barcodes'].tolist()

TPM_exp = []

print(len(all_hugos))

b=0
while b<len(all_hugos):
    df_temp = df_X[df_X['gene_symbol']==all_hugos[b]]
        
    val = df_temp[all_barcodes[b]].tolist()
    
    if len(val) > 1:
        print("ERROR")
    
    try:
        TPM_exp.append(val[0])
    except IndexError:
        TPM_exp.append(0)        
    
    b=b+1
    
    if b % 1000 == 0:
    
        print(str(b))

df2['Gene_TPM_Xena'] = TPM_exp

print(df2)

df2 = df2[df2['Gene_TPM_Xena']>=2]

print(df2)

df2_barcodes = df2['Truncated_Barcodes'].tolist()

pur_vals = []

for a in df2_barcodes:
    try:
        pur_vals.append(cn_dict[a])
    except KeyError:
        pur_vals.append(float("nan"))        

df2['Purity'] = pur_vals

ploidy_vals = []

for a in df2_barcodes:
    try:
        ploidy_vals.append(ploidy_dict[a])
    except KeyError:
        ploidy_vals.append(float("nan"))        

df2['Ploidy'] = ploidy_vals

XIST_exp_vals = []
class_vals = []
second_class_vals = []

XIST_pos = []

for a in df2_barcodes:
    XIST_exp_vals.append(XIST_dict[a])
    
    if XIST_dict[a] >= 3:
        XIST_pos.append(1)
    else:
        XIST_pos.append(0)
    
    class_vals.append(classification_dict[a])
    second_class_vals.append(second_classification_dict[a])

df2['XIST_TPM'] = XIST_exp_vals
df2['XIST_positive_1_is_yes'] = XIST_pos
df2['Classification'] = class_vals
df2['Secondary_Class'] = second_class_vals

predicted_X = []

for a in ploidy_vals:
    if a <= 2:
        predicted_X.append(1)
    elif a <= 4:
        predicted_X.append(2)
    elif a <= 6:
        predicted_X.append(3)
    elif a <= 8:
        predicted_X.append(4)        
    else:
        predicted_X.append(float("nan"))      

#df2['Predicted_X_copies'] = predicted_X

df2_ref = df2['t_ref_count'].tolist()
df2_alt = df2['t_alt_count'].tolist()

df2_VAF = []

a=0
while a<len(df2_ref):
    df2_VAF.append(df2_alt[a]/(df2_alt[a] + df2_ref[a]))
    a=a+1

df2['VAF'] = df2_VAF

#corrected VAF: corrected_VAF = (VAF)*(1+purity*(chrX ploidy-1))/(purity*(chrX ploidy))

corrected_VAF = []

a=0
while a<len(pur_vals):
    new_VAF = df2_VAF[a]*(1+pur_vals[a]*(predicted_X[a]-1))/(pur_vals[a]*predicted_X[a])
    corrected_VAF.append(new_VAF)
    a=a+1

#df2['Corrected_VAF'] = corrected_VAF

#df2 = df2[df2['Purity']>=0.5]

print(len(list(set(df2['Truncated_Barcodes'].tolist()))))

df2.to_csv("../VAF_Analysis/Output_Files/chrX_muts_REVISED_females_unaveraged_TPM_geq2.csv", index=False)

#ploidy_dict
