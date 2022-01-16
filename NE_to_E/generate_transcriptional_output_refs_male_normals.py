#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 25 21:12:46 2021

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
from scipy.stats import skewtest

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.linewidth'] = 0.3

figsize = (4.11,3.52) # the size of the figure - changes the shape of the squares in the comut
dpi_set = 72 # change the output resolution
sns.set(font_scale=0.6)

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/skewness_male_normals_normalized_E_expression_TPM_geq0.csv")
cols = df.columns.tolist()

del cols[0]

median_vals = []

for a in cols:
    curr_vals = [x for x in df[a].tolist() if x == x]
    median_vals.append(statistics.median(curr_vals))
    
E_dict = dict(zip(cols, median_vals))


df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/skewness_male_normals_normalized_NE_expression_TPM_geq0.csv")
cols = df.columns.tolist()

del cols[0]

median_vals = []

for a in cols:
    curr_vals = [x for x in df[a].tolist() if x == x]
    median_vals.append(statistics.median(curr_vals))

NE_dict = dict(zip(cols, median_vals))

df_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

df = df_ref[df_ref['Classification']!="UNKNOWN"]

df_barcode = df['Barcode'].tolist()

end_barcode = []

for a in df_barcode:
    end_barcode.append(int(split_advanced(a, "-", 3)[1]))

df['End_Val'] = end_barcode

df = df[df['End_Val']==11]

del df['End_Val']

temp_class = []

orig_cols = df.columns.tolist()

curr_barcode = df['Barcode'].tolist()

trunc_barcode = []

for a in curr_barcode:
    trunc_barcode.append(split_advanced(a, "-", 3)[0])

"""
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
    temp_second_class = temp_df['Secondary_Class'].tolist()[0]
    
    new_rows.append([str(a)+"-09", temp_log2, temp_average, temp_classification, temp_gender, temp_second_class, "UNKNOWN"]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)

df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Barcode'].str.contains(a))]

df_ref = pd.concat([df, df_from_new_rows], ignore_index=True)
"""

bar = df_ref['Barcode'].tolist()
XIST_TPM = df_ref['XIST_TPM'].tolist()
class_ref = df_ref['Classification'].tolist()
sex_ref = df_ref['Gender'].tolist()
second_class_ref = df_ref['Secondary_Class'].tolist()

TPM_dict = dict(zip(bar, XIST_TPM))
class_dict = dict(zip(bar, class_ref))
sex_dict = dict(zip(bar, sex_ref))
second_class_dict = dict(zip(bar, second_class_ref))

NE_median_gene_level_ratio = []
E_median_gene_level_ratio = []
ratio_ne_e = []
XIST_TPM = []
Classification = []
Secondary_Class = []

for a in cols:
    NE_median_gene_level_ratio.append(2**NE_dict[a])
    E_median_gene_level_ratio.append(2**E_dict[a])
    ratio_ne_e.append(2**NE_dict[a]/2**E_dict[a])
    XIST_TPM.append(TPM_dict[a])
    Classification.append(class_dict[a])
    Secondary_Class.append(second_class_dict[a])
    
df = pd.DataFrame([cols, NE_median_gene_level_ratio, E_median_gene_level_ratio, ratio_ne_e, XIST_TPM, Classification, Secondary_Class]).T

df.columns = ["sample", "NE_median_gene_level_ratio", "E_median_gene_level_ratio", "ratio_ne_e", "XIST_TPM", "Classification", "Secondary_Class"]

trunc_barcode = []

for a in df["sample"].tolist():
    trunc_barcode.append(split_advanced(a, "-", 3)[0])
    
dup_items = list(set([item for item, count in collections.Counter(trunc_barcode).items() if count > 1]))

if dup_items:
    
    new_rows = []
    
    for a in dup_items:
        temp_df = df[df['sample'].str.contains(a)]
        vals = temp_df['NE_median_gene_level_ratio'].tolist()
        
        if len(vals)<2:
            print("ERROR")
        
        vals2 = temp_df['E_median_gene_level_ratio'].tolist()
        vals3 = temp_df['XIST_TPM'].tolist()
        
        temp_vals_average = sum(vals)/len(vals)
        temp_vals2_average = sum(vals2)/len(vals2)
        temp_vals3_average = sum(vals3)/len(vals3)
        
        temp_vals_vals2_ratio = temp_vals_average/temp_vals2_average
        
        temp_second_class = temp_df['Secondary_Class'].tolist()[0]
        temp_classification = temp_df['Classification'].tolist()[0]
        
        new_rows.append([str(a)+"-09", temp_vals_average, temp_vals2_average, temp_vals_vals2_ratio, temp_vals3_average, temp_classification, temp_second_class])
    
    df_from_new_rows = pd.DataFrame(new_rows)
    
    df_from_new_rows.columns = ["sample", "NE_median_gene_level_ratio", "E_median_gene_level_ratio", "ratio_ne_e", "XIST_TPM", "Classification", "Secondary_Class"]
    
    for a in dup_items:
        df = df[~(df['Barcode'].str.contains(a))]
    
    df = pd.concat([df, df_from_new_rows], ignore_index=True)

df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/normals_only_MALES_resubset_09averaged_no_PAR_no_segment_mean_normalization_median_ne_e_pseudo_normalized_sample_level_ratio.csv", index=False)


