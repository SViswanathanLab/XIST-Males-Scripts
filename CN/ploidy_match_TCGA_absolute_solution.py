#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 12:45:22 2022

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
import requests, sys
import glob,os
import scipy.stats
import statistics
from os.path import exists
from itertools import islice

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(9,5)})
sns.set(font_scale=1)
plt.rcParams.update({'font.size': 14})

ploidy_vals_to_test = [2, 3, 4, 5]

path = "../CN/Other_Input/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"

df = pd.read_csv(path, sep="\t")
df = df[df['ploidy'].notna()]

sample = df['array'].tolist()
ploidy = df['ploidy'].tolist()

ploidy_dict = dict(zip(sample, ploidy))

"""
new_ploidy = []

for a in ploidy:
    rounded_ploidy = round(a)
    if rounded_ploidy > 5:
        new_ploidy.append(5)
    elif rounded_ploidy < 2:
        new_ploidy.append(2)        
    else:
        new_ploidy.append(rounded_ploidy)

ploidy_dict = dict(zip(sample, new_ploidy))
"""

path2 = "../CN/TCGA_VAF_all_samples_ref_for_Titan_baitset_annotated.csv"

df2 = pd.read_csv(path2)
df2 = df2[df2['Baitset']!="Unknown"]

print(len(df2.index.tolist()))

id_val = df2['IDs'].tolist()
sub_ids = df2['Sub_ID'].tolist()
baitset = df2['Baitset'].tolist()

df_ref = pd.read_csv("../CN/Other_Input/table_s1_v4.csv")

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

for a in invalid_ids:
    df_ref = df_ref[~(df_ref['Barcode'].str.contains(a))]

df_males = df_ref[df_ref['Sex_Custom'] == "MALE"]
df_females = df_ref[df_ref['Sex_Custom'] == "FEMALE"]

male_barcodes = df_males['Barcode'].tolist()
female_barcodes = df_females['Barcode'].tolist()

sex_vals = []

for a in id_val:
    if a in male_barcodes:
        sex_vals.append("MALE")
    elif a in female_barcodes:
         sex_vals.append("FEMALE")
    else:
         sex_vals.append("UNKNOWN")
     
ploidy_val = []   
     
for a in id_val:
    try:
        ploidy_val.append(ploidy_dict[a])
    except:
        ploidy_val.append("Unknown, Used TITAN optimal solution")        


output_dir_list = []

male_path = "TCGA_TITAN_r2_male/"
female_path = "TCGA_TITAN_r2/"

agilent1_path = "Agilent_1_incorrect_baitset/"
agilent2_path = "Agilent_2/"
hg18_path = "hg18_baitset/"
SQEZ2_path = "SQEZ2_baitset/"
SQEZ3_path = "SQEZ3_baitset/"
SureSelect_path = "SureSelect_baitset/"
VCRome_path = "HGSC_baitset/"

path = male_path + agilent1_path + "titan/hmm/optimalClusterSolution/"
agilent1_male_ids = [f.split("_Tumor")[0] for f in os.listdir(path) if f.endswith('.segs.txt')]

path = female_path + agilent1_path + "titan/hmm/optimalClusterSolution/"
agilent1_female_ids = [f.split("_Tumor")[0] for f in os.listdir(path) if f.endswith('.segs.txt')]

path = male_path + agilent2_path + "titan/hmm/optimalClusterSolution/"
agilent2_male_ids = [f.split("_Tumor")[0] for f in os.listdir(path) if f.endswith('.segs.txt')]

path = female_path + agilent2_path + "titan/hmm/optimalClusterSolution/"
agilent2_female_ids = [f.split("_Tumor")[0] for f in os.listdir(path) if f.endswith('.segs.txt')]

optimal_files = []
avg_chrX_CN = []
rounded_avg_chrX_CN = []

a=0
while a<len(id_val):
    
    valid_paths = []
    
    if baitset[a] == "Agilent Custom":
        if sub_ids[a] in agilent1_male_ids or sub_ids[a] in agilent1_female_ids:   
            baitset_path = agilent1_path
        elif sub_ids[a] in agilent2_male_ids or sub_ids[a] in agilent2_female_ids: 
            baitset_path = agilent2_path            
    elif baitset[a] == "Nimblegen hg18":
        baitset_path = hg18_path
    elif baitset[a] == "Nimblegen HGSC":
        baitset_path = VCRome_path
    elif baitset[a] == "Nimblegen.SQEZ2":
        baitset_path = SQEZ2_path
    elif baitset[a] == "Nimblegen.SQEZ3":
        baitset_path = SQEZ3_path
    elif baitset[a] == "Sureselect.38":
        baitset_path = SureSelect_path
    else:
        optimal_files.append("No Baitset Found")                   
        a=a+1
        continue
    
    if sex_vals[a] == "MALE":
        base_path = male_path
    elif sex_vals[a] == "FEMALE":
        base_path = female_path
    else:
        optimal_files.append("No Sex Annotated")    
        a=a+1
        continue
    
    ploidy_vals_new = []

    ploidy_ref = ploidy_val[a]
    if ploidy_ref == "Unknown, Used TITAN optimal solution":
        total_path_file_clus1 = base_path + baitset_path + "titan/hmm/optimalClusterSolution/" + sub_ids[a] + "_Tumor_cluster1.titan.ichor.seg.txt" 
        total_path_file_clus2 = base_path + baitset_path + "titan/hmm/optimalClusterSolution/" + sub_ids[a] + "_Tumor_cluster2.titan.ichor.seg.txt" 
        
        if exists(total_path_file_clus1):
            optimal_files.append(total_path_file_clus1)
        elif exists(total_path_file_clus1):
            optimal_files.append(total_path_file_clus2)   
        else:
            optimal_files.append("TITAN not run on Sample (unknown ABSOLUTE ploidy)")                   

    else:
        for q in ploidy_vals_to_test:
            total_path_clus1 = base_path + baitset_path + "titan/hmm/titanCNA_ploidy" + str(q) + "/" + sub_ids[a] + "_Tumor_cluster1.params.txt"
            total_path_clus2 = base_path + baitset_path + "titan/hmm/titanCNA_ploidy" + str(q) + "/" + sub_ids[a] + "_Tumor_cluster2.params.txt"
            
            valid_paths.append(total_path_clus1)
            valid_paths.append(total_path_clus2)
        
        existing_path = []
        
        for b in valid_paths:
            if exists(b):
                existing_path.append(b)
                with open(b) as fin:
                    for line in islice(fin, 1, 2):
                        ploidy_vals_new.append(float(line.split(":	")[1]))
            
        diff_list = []
        
        if ploidy_vals_new:
            for c in ploidy_vals_new:
                diff_list.append(abs(c-ploidy_ref))

            min_value = min(diff_list)
            min_index = diff_list.index(min_value)
                        
            temp_path = existing_path[min_index].split(".params.txt")[0]
            optimal_files.append(temp_path + ".titan.ichor.seg.txt")
        
        else:
            optimal_files.append("TITAN not run on Sample (known ABSOLUTE ploidy)")
            
    if exists(optimal_files[a]):
        combined_df = pd.read_csv(optimal_files[a], sep="\t")
        combined_df = combined_df[combined_df['Chromosome']=="X"]
        combined_df['Seg_Length'] = combined_df['End']-combined_df['Start']
        combined_df['Seg_Product'] = combined_df['Seg_Length']*combined_df['Corrected_Copy_Number']
        sum_length = sum(combined_df['Seg_Length'].tolist())
        pdt_length = sum(combined_df['Seg_Product'].tolist())
        try:
            val = pdt_length/sum_length
            avg_chrX_CN.append(val)
            rounded_avg_chrX_CN.append(round(val))
        except ZeroDivisionError:
            avg_chrX_CN.append("Unknown")
            rounded_avg_chrX_CN.append("Unknown")     
            
    else:
        avg_chrX_CN.append("Unknown")
        rounded_avg_chrX_CN.append("Unknown")
    
    if a % 100 == 1:
        
        df_out = pd.DataFrame([id_val, sub_ids, optimal_files, avg_chrX_CN, rounded_avg_chrX_CN, ploidy_val, baitset, sex_vals]).T
        cols = ['ID', 'Sub_ID', 'Best_Solution', 'Unrounded_chrX_CN', 'Rounded_chrX_CN', 'ABSOLUTE_ploidy', 'Baitset', 'Sex']
        df_out.columns = cols
        
        df_out.to_csv("../CN/Output_Files/TITAN_final_solutions.csv", index=False)
    
    a=a+1
    print(a)

df_out = pd.DataFrame([id_val, sub_ids, optimal_files, avg_chrX_CN, rounded_avg_chrX_CN, ploidy_val, baitset, sex_vals]).T
cols = ['ID', 'Sub_ID', 'Best_Solution', 'Unrounded_chrX_CN', 'Rounded_chrX_CN', 'ABSOLUTE_ploidy', 'Baitset', 'Sex']
df_out.columns = cols

df_out.to_csv("../CN/Output_Files/TITAN_final_solutions.csv", index=False)



