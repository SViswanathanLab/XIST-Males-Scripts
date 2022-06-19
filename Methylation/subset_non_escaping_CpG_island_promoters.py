#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 21:49:23 2021

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

df_annot = pd.read_csv("../Methylation/Other_Input/methylation_annotation_ref.csv")

df_annot = df_annot[df_annot['Chromosome']=="chrX"]
df_annot = df_annot[df_annot['Feature_Type']=="Island"]

df_gene_class = pd.read_csv("../Methylation/Other_Input/25381334_escaping_and_nonescaping_TSSs_methylation.csv")
    
df_inactivated = df_gene_class[df_gene_class['450k_DNAm_status']=="subject to XCI in 27 tissues"]
#df_inactivated = df_gene_class[df_gene_class['450k_DNAm_status']=="escape from XCI in 27 tissues"]

df_escaping = df_gene_class[~(df_gene_class['450k_DNAm_status']=="subject to XCI in 27 tissues")]
#df_escaping = df_gene_class[~(df_gene_class['450k_DNAm_status']=="escape from XCI in 27 tissues")]

inactivated_genes = df_inactivated['TSS_name'].tolist()

escaping_genes = df_escaping['TSS_name'].tolist()

compound_list = df_annot['Gene_Symbol'].tolist()
ref_list = df_annot['Composite Element REF'].tolist()
position_toTSS = df_annot['Position_to_TSS'].tolist()

islands_of_interest = []

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3202277/ fig. 3 - +/- 125

w=0
while w<len(compound_list):
    curr_list = (compound_list[w]).split(";")
    curr_pos = (position_toTSS[w]).split(";")
    abs_pos = []
    
    for b in curr_pos:
        if b != ".":
            abs_pos.append(abs(float(b)))
    
    try:
        min_pos = min(abs_pos)
    except ValueError:
        w=w+1
        continue
    
    uq_items = list(set(curr_list))
        
    escaping_genes_counter = 0
    inactivated_genes_counter = 0
    
    for q in uq_items:
        if q in escaping_genes:
            escaping_genes_counter = escaping_genes_counter + 1
        elif q in inactivated_genes:
            inactivated_genes_counter = inactivated_genes_counter + 1
    
    if inactivated_genes_counter >= 1 and escaping_genes_counter == 0 and min_pos <= 125:
        islands_of_interest.append(ref_list[w])
    
    w=w+1

df_out = df_annot[df_annot['Composite Element REF'].isin(islands_of_interest)]
df_out.to_csv("../Methylation/Output_Files/chrX_non_escaping_CpG_promoter_islands.csv", index=False)
