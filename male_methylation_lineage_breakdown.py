#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 11:39:53 2021

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

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

main_interest = "MALES"
geneset_of_interest = "AUTOSOMAL"


color_to_use = "#B2D3A8"
#color_to_use = "#EDE5A6"
#color_to_use = "#FCB1A6"
#color_to_use = "#C5D8D1"



dpi_set = 72
sns.set(rc={'figure.figsize':(15,5)})
sns.set(font_scale=1)

if geneset_of_interest == "AUTOSOMAL":
    sns.set(rc={'figure.figsize':(15,5)})


if main_interest == "MALES":

    entire_dataset = [['BLCA', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'PRAD']]

elif main_interest == "FEMALES":

    entire_dataset = [['UVM', 'OV', 'SKCM', 'ACC', 'KICH', 'KIRC', 'KIRP', 'UCS', 'HNSC', 'GBM', 'COAD', 'LUSC', 'LUAD', 'SARC', 'STAD', 'LAML', 'BRCA', 'UCEC', 'CESC']]
 
#entire_dataset = [['BLCA'], ['ESCA'], ['HNSC'], ['KICH'], ['KIRC'], ['LAML'], ['LGG'], ['LIHC'], ['LUAD'], ['LUSC'], ['PAAD'], ['PCPG'], ['PRAD'], ['READ'], ['SARC'], ['SKCM'], ['STAD'], ['THCA']]
                  
aaaaa=0

final_df = pd.DataFrame()

fig, ax = plt.subplots(1, 1, sharey=True, tight_layout=True)

for all_classes in entire_dataset:
                
    combined_XISTpos_methylation_vals = []
    combined_XISTneg_methylation_vals = []
    combined_XISTpos_females = []
    combined_XISTneg_females = []
    
    NSGCT_XISTpos_vals = []
    NSGCT_XISTneg_vals = []    
        
    aaa=0
    
    for www in all_classes:
        
        class_of_interest = www
        
        if class_of_interest != "Cancer-Adjacent Normal":
            
            main_cutoff = 1
            
            plt.rcParams["font.family"] = "Arial"
            sns.set(rc={'figure.figsize':(5,8)})
            sns.set(font_scale=1)
            plt.rcParams.update({'font.size': 14})
            plt.rcParams['axes.linewidth'] = 0.3
            dpi_set = 72 # change the output resolution
            
            df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/NEW_male_tumors_sample_level_methylation_XISTpos_samples.csv")
            df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/NEW_male_tumors_sample_level_methylation_XISTneg_samples.csv")
            df3 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/NEW_female_tumors_sample_level_methylation_XISTpos_samples.csv")
            df4 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/NEW_female_tumors_sample_level_methylation_XISTneg_samples.csv")

            df4 = df4[df4['sample']!="TCGA-25-1870-01"]
  
            invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

            for a in invalid_ids:
                df = df[~(df['sample'].str.contains(a))]
                df2 = df2[~(df2['sample'].str.contains(a))]
                df3 = df3[~(df3['sample'].str.contains(a))]
                df4 = df4[~(df4['sample'].str.contains(a))]          
  
            if www != "TGCT-NS":
                df = df[df['Classification']==www]
                df2 = df2[df2['Classification']==www]
                df3 = df3[df3['Classification']==www]
                df4 = df4[df4['Classification']==www]
            else:
                df = df[df['Secondary_Class']==www]
                df2 = df2[df2['Secondary_Class']==www]                
                df3 = df3[df3['Secondary_Class']==www]
                df4 = df4[df4['Secondary_Class']==www]
                
            orig_df_ids = df['sample'].tolist()
            orig_df2_ids = df2['sample'].tolist()
            orig_df3_ids = df3['sample'].tolist()
            orig_df4_ids = df4['sample'].tolist()
                        
            concat_df = pd.concat([df, df2, df3, df4])
            
            if geneset_of_interest == "AUTOSOMAL":
            
                all_methylation = concat_df['autosome_mean_beta'].tolist()
            
            elif geneset_of_interest == "NONESCAPE":
            
                all_methylation = concat_df['chrX_non_escaping_promoter_CpG_islands_mean_beta'].tolist()
         
                        
            unprocessed_values = all_methylation
            
            mean = sum(unprocessed_values) / len(unprocessed_values) 
            variance = sum([((x - mean) ** 2) for x in unprocessed_values]) / (len(unprocessed_values)-1) 
            standard_dev = variance**0.5
            
            ys_old = []
            for t in unprocessed_values:
                #z_score = (t-mean)/standard_dev
                #z_score = math.log((t/mean), 2)
                z_score = t
                ys_old.append(z_score)
                
            ys_methylation = ys_old
            
            concat_df['z_scored_methylation'] = ys_methylation
            
            df_with_z = concat_df[concat_df['sample'].isin(orig_df_ids)]
            df2_with_z = concat_df[concat_df['sample'].isin(orig_df2_ids)]
            
            if www != "TGCT-NS":            
                df3_with_z = concat_df[concat_df['sample'].isin(orig_df3_ids)]
                df4_with_z = concat_df[concat_df['sample'].isin(orig_df4_ids)]
            
            vals1 = df_with_z['z_scored_methylation'].tolist()
            vals2 = df2_with_z['z_scored_methylation'].tolist()
        
            if www != "TGCT-NS":
                vals3 = df3_with_z['z_scored_methylation'].tolist()
                vals4 = df4_with_z['z_scored_methylation'].tolist()
                
            vals1_new = []
            vals2_new = []
            vals3_new = []
            vals4_new = []
                
            q=0
            while q<2000:
                
                if q<len(vals1):
                    vals1_new.append(vals1[q])
                else:
                    vals1_new.append(float("nan"))

                if q<len(vals2):
                    vals2_new.append(vals2[q])
                else:
                    vals2_new.append(float("nan"))
                    
                if q<len(vals3):
                    vals3_new.append(vals3[q])
                else:
                    vals3_new.append(float("nan"))

                if q<len(vals4):
                    vals4_new.append(vals4[q])
                else:
                    vals4_new.append(float("nan"))
                    
                q=q+1
                    
            if main_interest == "MALES":    
    
                final_df['M-' + str(class_of_interest) + '-'] = vals2_new
                final_df['M-' + str(class_of_interest) + '+'] = vals1_new
                
            elif main_interest == "FEMALES":    
                
                final_df['F-' + str(class_of_interest) + '+'] = vals3_new
                final_df['F-' + str(class_of_interest) + '-'] = vals4_new
            
            #ax.set_ylabel("log$_2$(Lineage-Mean Normalized chrX Non-Escape : Escape Transcriptional Output Ratio)")
        
        aaa=aaa+1

    new_df = pd.DataFrame()

    temp_final_df_cols = final_df.columns.tolist()
    
    col_num = []
    diff_list = []

    j=0
    while j<len(temp_final_df_cols):
        col1 = final_df[temp_final_df_cols[j]].tolist()
        col2 = final_df[temp_final_df_cols[j+1]].tolist()
        
        col1 = [x for x in col1 if x==x]
        col2 = [x for x in col2 if x==x]
        
        diff = statistics.median(col2) - statistics.median(col1)
        
        if main_interest == "MALES":
            diff = -1*diff
        
        col_num.append(j)
        diff_list.append(diff)
        j=j+2
        
    col_dict = dict(zip(diff_list, col_num))
    
    diff_list = sorted(diff_list)    
    
    j=0
    while j<len(diff_list):
        value_1 = col_dict[diff_list[j]]
        value_2 = value_1 + 1   
        
        name_1 = temp_final_df_cols[value_1]
        name_2 = temp_final_df_cols[value_2]
        
        new_df[name_1] = final_df[name_1]
        new_df[name_2] = final_df[name_2]
        
        j=j+1
    
    final_df = new_df
    
    temp_labels = []
    
    for a in final_df.columns.tolist():
        temp_list = final_df[a].tolist()
        proc_list = [x for x in temp_list if x == x]
        temp_labels.append(a + " (" + str(len(proc_list)) + ")")
    
    final_df.columns = temp_labels

    if main_interest == "MALES":
        temp_palette = []
        for a in temp_labels:
            if "+" in a:
                temp_palette.append("#C17C74")
            else:
                temp_palette.append("#669BBC")
    else:
        temp_palette = []
        for a in temp_labels:
            if "+" in a:
                temp_palette.append("#8AA29E")
            else:
                temp_palette.append("#F6CA83")   
    
    ax.set_xticklabels(final_df.columns.tolist(), rotation=90)

    sns.boxplot(data=final_df, ax=ax, showfliers=False, palette=temp_palette)

    sns.stripplot(data=final_df, color="black", s=1.5, ax=ax, jitter=0.2, edgecolor="none", alpha=0.9, zorder=1000, linewidth=0.15)
    
    ax.grid(False)
    ax.set_facecolor("white")
    ax.spines['bottom'].set_color('0')
    ax.spines['left'].set_color('0')
    
    ax.tick_params(bottom='on', left='on')
    
    if geneset_of_interest == "NONESCAPE":
    
        ax.set_ylabel("Median chrX Methylation at Non-Escaping\nPromoter CpG Islands in Sample (β)")

    elif geneset_of_interest == "AUTOSOMAL":

        ax.set_ylabel("Median Autosome Methylation in Sample (β)")
        
    if geneset_of_interest == "NONESCAPE":
        a=0
        while a<len(final_df.columns.tolist()):
            add_stat_annotation(ax, data=final_df,
                            box_pairs=[(final_df.columns.tolist()[a], final_df.columns.tolist()[a+1])],
                            test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.2, line_offset_to_box=1/20+0.2, line_offset=0.3, fontsize=8)
            a=a+2
    else:
        a=0
        while a<len(final_df.columns.tolist()):
            add_stat_annotation(ax, data=final_df,
                            box_pairs=[(final_df.columns.tolist()[a], final_df.columns.tolist()[a+1])],
                            test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.2, line_offset_to_box=1/20+0.2, line_offset=0, fontsize=8)
            a=a+2
    #ax.set_ylim([-2, 4])
    
    aaaaa = aaaaa + 1

ax.set_xticklabels(final_df.columns.tolist(), rotation=90)

plt.ylim([0, 0.7])
    
if geneset_of_interest == "AUTOSOMAL" and main_interest == "FEMALES":

    plt.ylim([0, 0.39])

if geneset_of_interest == "AUTOSOMAL" and main_interest == "MALES":

    plt.ylim([0, 0.39])

if geneset_of_interest == "NONESCAPE" and main_interest == "MALES":

    plt.ylim([0, 0.55])
#plt.ylim([0, 6])

plt.tight_layout()

if geneset_of_interest == "NONESCAPE":
    
    if main_interest == "MALES":
        
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_nonescape_lineage_breakdown_of_methylation_all_categories.png", dpi=dpi_set)
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_nonescape_lineage_breakdown_of_methylation_all_categories.pdf", dpi=dpi_set)
        
    elif main_interest == "FEMALES":
      
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/FEMALES_nonescape_lineage_breakdown_of_methylation_all_categories.png", dpi=dpi_set)
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/FEMALES_nonescape_lineage_breakdown_of_methylation_all_categories.pdf", dpi=dpi_set)

elif geneset_of_interest == "AUTOSOMAL":
        
    if main_interest == "MALES":
        
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_autosomal_lineage_breakdown_of_methylation_all_categories.png", dpi=dpi_set)
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/MALES_autosomal_lineage_breakdown_of_methylation_all_categories.pdf", dpi=dpi_set)
        
    elif main_interest == "FEMALES":
      
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/FEMALES_autosomal_lineage_breakdown_of_methylation_all_categories.png", dpi=dpi_set)
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/FEMALES_autosomal_lineage_breakdown_of_methylation_all_categories.pdf", dpi=dpi_set)
        