#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 21:38:32 2021

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

df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/mc3.v0.2.8.PUBLIC.maf", sep="\t")

w = df2['Variant_Classification'].tolist()
print(list(set(w)))

df2 = df2[df2['Chromosome'].isin([23, "23", "X"])]

#df2 = df2[df2['Variant_Classification'].isin(['Missense_Mutation', 'Frame_Shift_Ins', 'Frame_Shift_Del', 'Translation_Start_Site', 'Nonsense_Mutation', 'Splice_Site', 'Nonstop_Mutation', 'In_Frame_Del', 'In_Frame_Ins', 'Silent', 'Targeted_Region', 'RNA'])]

df_barcodes = df2['Tumor_Sample_Barcode'].tolist()

trunc_barcodes = []

for a in df_barcodes:
    trunc_barcodes.append(split_advanced(a, "-", 4)[0][:-1])
    
df2['Truncated_Barcodes'] = trunc_barcodes

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

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
df = df[df['Gender'] == "MALE"]

male_valid = ['BLCA', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'PRAD']

female_valid = ['UVM', 'OV', 'SKCM', 'ACC', 'KICH', 'KIRC', 'KIRP', 'UCS', 'HNSC', 'GBM', 'COAD', 'LUSC', 'LUAD', 'SARC', 'STAD', 'LAML', 'BRCA', 'UCEC', 'CESC']

df = df[df['Classification'].isin(male_valid)]

id_list = df['Barcode'].tolist()

df2 = df2[df2['Truncated_Barcodes'].isin(id_list)]

df2.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_muts_REVISED_8_28_males_unaveraged.txt", sep="\t", index=False)

