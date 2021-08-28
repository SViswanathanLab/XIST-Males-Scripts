#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 11:39:53 2021

@author: ananthansadagopan
"""


import pandas as pd
import collections
import math
import numpy as npi
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

plt.rcParams["font.family"] = "Arial"

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

geneset_of_interest = "ESCAPE" #ESCAPE

dpi_set = 72
#sns.set(rc={'figure.figsize':(4.21,2.37)})
#sns.set(font_scale=0.66)

sns.set(rc={'figure.figsize':(4.29,2.396)})
sns.set(font_scale=0.5)

male_valid = ['BLCA', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'PRAD']

female_valid = ['UVM', 'OV', 'SKCM', 'ACC', 'KICH', 'KIRC', 'KIRP', 'UCS', 'HNSC', 'GBM', 'COAD', 'LUSC', 'LUAD', 'SARC', 'STAD', 'LAML', 'BRCA', 'UCEC', 'CESC']

"""
male_valid = list(set(male_valid) & set(female_valid))
female_valid = male_valid
"""

all_valid = list(set(male_valid + female_valid)) 

entire_dataset = ['TGCT-NS']

for a in all_valid:
    entire_dataset.append(a)

entire_dataset = [entire_dataset]
      
aaaaa=0

fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=True)

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
        #curr_color = new_color[aaa]
        
        if class_of_interest != "Cancer-Adjacent Normal":
            
            main_cutoff = 1
            
            plt.rcParams["font.family"] = "Arial"
            sns.set(rc={'figure.figsize':(5,8)})
            sns.set(font_scale=1)
            plt.rcParams.update({'font.size': 14})
            plt.rcParams['axes.linewidth'] = 0.3
            dpi_set = 72 
            
            
            df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_XISTpos_samples.csv")
            df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_XISTneg_samples.csv")
            df3 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/FEMALES_no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_XISTpos_samples.csv")
            df4 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/FEMALES_no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_XISTneg_samples.csv")

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
            
            if geneset_of_interest == "ESCAPE":
            
                all_methylation = concat_df['ratio_ne_e'].tolist()
                
            elif geneset_of_interest == "AUTOSOMAL":
                
                all_methylation = concat_df['ratio_ne_a'].tolist()                
                        
            unprocessed_values = all_methylation
            
            unprocessed_values = [math.log(x,2) for x in unprocessed_values]
            
            mean = sum(unprocessed_values) / len(unprocessed_values) 
            variance = sum([((x - mean) ** 2) for x in unprocessed_values]) / (len(unprocessed_values)-1) 
            standard_dev = variance**0.5
            
            ys_old = []
            for t in unprocessed_values:
                #z_score = (t-mean)/standard_dev
                #z_score = math.log((t/mean), 2)
                #z_score = math.log(t, 2)
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
            
            if www != "TGCT-NS":
                
                if www not in male_valid:
                    master_df = pd.DataFrame([[float("nan")], [float("nan")], [float("nan")], [float("nan")], vals3, vals4]).T
                    master_df.columns = ['XIST- Male\nNSGCT', 'XIST+ Male\nNSGCT', 'XIST- Male\nPanCan', 'XIST+ Male\nPanCan', 'XIST+ Female\nPanCan', 'XIST- Female\nPanCan']
                elif www not in female_valid:
                    master_df = pd.DataFrame([[float("nan")], [float("nan")], vals2, vals1, [float("nan")], [float("nan")]]).T
                    master_df.columns = ['XIST- Male\nNSGCT', 'XIST+ Male\nNSGCT', 'XIST- Male\nPanCan', 'XIST+ Male\nPanCan', 'XIST+ Female\nPanCan', 'XIST- Female\nPanCan']
                else:              
                    master_df = pd.DataFrame([[float("nan")], [float("nan")], vals2, vals1, vals3, vals4]).T
                    master_df.columns = ['XIST- Male\nNSGCT', 'XIST+ Male\nNSGCT', 'XIST- Male\nPanCan', 'XIST+ Male\nPanCan', 'XIST+ Female\nPanCan', 'XIST- Female\nPanCan']
            
            else:
                master_df = pd.DataFrame([vals2, vals1, [float("nan")], [float("nan")], [float("nan")], [float("nan")]]).T
                master_df.columns = ['XIST- Male\nNSGCT', 'XIST+ Male\nNSGCT', 'XIST- Male\nPanCan', 'XIST+ Male\nPanCan', 'XIST+ Female\nPanCan', 'XIST- Female\nPanCan']
                
            if www == "TGCT-NS":
                for abc in vals1:
                    NSGCT_XISTpos_vals.append(abc)
                for abc in vals2:             
                    NSGCT_XISTneg_vals.append(abc)
                    
            if www != "TGCT-NS":
            
                if www not in male_valid:
                    for abc in vals3:
                        combined_XISTpos_females.append(abc)
                    for abc in vals4:
                        combined_XISTneg_females.append(abc)
                elif www not in female_valid:
                    for abc in vals1:
                        combined_XISTpos_methylation_vals.append(abc)
                    for abc in vals2:
                        combined_XISTneg_methylation_vals.append(abc)
                else:              
                    for abc in vals1:
                        combined_XISTpos_methylation_vals.append(abc)
                    for abc in vals2:
                        combined_XISTneg_methylation_vals.append(abc)
                    for abc in vals3:
                        combined_XISTpos_females.append(abc)
                    for abc in vals4:
                        combined_XISTneg_females.append(abc)

            curr_color="black"
            
            df1 = master_df[['XIST- Male\nNSGCT', 'XIST+ Male\nNSGCT']]
            df2 = master_df[['XIST- Male\nPanCan', 'XIST+ Male\nPanCan']]
            df3 = master_df[['XIST+ Female\nPanCan', 'XIST- Female\nPanCan']]
            
            sns.stripplot(data=df1, color=curr_color, s=1, ax=axs[0], jitter=0.2, edgecolor="black", alpha=0.9, zorder=2, linewidth=0.15)
            sns.stripplot(data=df2, color=curr_color, s=1, ax=axs[1], jitter=0.2, edgecolor="black", alpha=0.9, zorder=2, linewidth=0.15)
            sns.stripplot(data=df3, color=curr_color, s=1, ax=axs[2], jitter=0.2, edgecolor="black", alpha=0.9, zorder=2, linewidth=0.15)
            
            q=0
            while q<3:
            
                axs[q].grid(False)
                axs[q].set_facecolor("white")
                axs[q].spines['bottom'].set_color('0')
                axs[q].spines['left'].set_color('0')
                axs[q].tick_params(bottom='on', left='on')
                q=q+1
                    
        aaa=aaa+1

    final_df = pd.DataFrame([NSGCT_XISTneg_vals, NSGCT_XISTpos_vals, combined_XISTneg_methylation_vals, combined_XISTpos_methylation_vals, combined_XISTpos_females, combined_XISTneg_females]).T

    final_df.columns = ['XIST-\n(' + str(len(NSGCT_XISTneg_vals)) + ")", 'XIST+\n(' + str(len(NSGCT_XISTpos_vals)) + ")", 'XIST-\n(' + str(len(combined_XISTneg_methylation_vals)) + ")", 'XIST+\n(' + str(len(combined_XISTpos_methylation_vals)) + ")", 'XIST+\n(' + str(len(combined_XISTpos_females)) + ")", 'XIST-\n(' + str(len(combined_XISTneg_females)) + ")"]
    
    palette1 = ['#BEB2C8', '#A8C686']
    palette2 = ['#669BBC', '#C17C74']
    palette3 = ['#8AA29E', '#F6CA83']

    df1 = final_df[['XIST-\n(' + str(len(NSGCT_XISTneg_vals)) + ")", 'XIST+\n(' + str(len(NSGCT_XISTpos_vals)) + ")"]]
    df2 = final_df[['XIST-\n(' + str(len(combined_XISTneg_methylation_vals)) + ")", 'XIST+\n(' + str(len(combined_XISTpos_methylation_vals)) + ")"]]
    df3 = final_df[['XIST+\n(' + str(len(combined_XISTpos_females)) + ")", 'XIST-\n(' + str(len(combined_XISTneg_females)) + ")"]]
    
    sns.boxplot(data=df1, ax=axs[0], showfliers=False, color="white", linewidth=0.5, palette=palette1)
    sns.boxplot(data=df2, ax=axs[1], showfliers=False, color="white", linewidth=0.5, palette=palette2)
    sns.boxplot(data=df3, ax=axs[2], showfliers=False, color="white", linewidth=0.5, palette=palette3)
    
    q=0
    while q<3:
        for i,artist in enumerate(axs[q].artists):
            col = "black"
            artist.set_edgecolor(col)
            j=0
            while j<6:
                line = axs[q].lines[j]
                line.set_color(col)
                line.set_mfc(col)
                line.set_mec(col)
                j=j+1
        q=q+1
    
    if geneset_of_interest == "ESCAPE":
    
        axs[0].set_ylabel("Normalized Non-Escapee to\nEscapee Expression (log$_2$Ratio)", fontsize=7)        

    elif geneset_of_interest == "AUTOSOMAL":

        axs[0].set_ylabel("Normalized Non-Escapee to\nAutosomal Gene Expression (log$_2$Ratio)", fontsize=7)

    add_stat_annotation(axs[0], data=df1,
                    box_pairs=[('XIST-\n(' + str(len(NSGCT_XISTneg_vals)) + ")", 'XIST+\n(' + str(len(NSGCT_XISTpos_vals)) + ")")],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.01, fontsize=7)
    add_stat_annotation(axs[1], data=df2,
                            box_pairs=[('XIST-\n(' + str(len(combined_XISTneg_methylation_vals)) + ")", 'XIST+\n(' + str(len(combined_XISTpos_methylation_vals)) + ")")],
                            test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.01, fontsize=7)
    add_stat_annotation(axs[2], data=df3,
                    box_pairs=[('XIST+\n(' + str(len(combined_XISTpos_females)) + ")", 'XIST-\n(' + str(len(combined_XISTneg_females)) + ")")],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.01, fontsize=7)
    
    aaaaa = aaaaa + 1

q=0
while q<3:
    for line in axs[q].get_lines():
        line.set_color('black')
        line.set_linewidth(0.5)
    q=q+1
  
q=0
while q<3:
    axs[q].xaxis.labelpad = 1
    axs[q].yaxis.labelpad = 1
    axs[q].spines['bottom'].set_lw(0.5)
    axs[q].spines['left'].set_lw(0.5)
    axs[q].xaxis.set_tick_params(width=0.5, size=3, pad=2)
    axs[q].yaxis.set_tick_params(width=0.5, size=3, pad=2)
    q=q+1

plt.setp(axs[0].get_yticklabels(), fontsize=7)

plt.tight_layout()

if geneset_of_interest == "ESCAPE":
    
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/gene_level_nonescaping_v_escaping_no_segment_normalization_male_and_female_all_categories.png", dpi=dpi_set)
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/gene_level_nonescaping_v_escaping_no_segment_normalization_male_and_female_all_categories.pdf", dpi=dpi_set)
    
elif geneset_of_interest == "AUTOSOMAL":

    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/gene_level_nonescaping_v_autosomal_no_segment_normalization_male_and_female_all_categories.png", dpi=dpi_set)
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/gene_level_nonescaping_v_autosomal_no_segment_normalization_male_and_female_all_categories.pdf", dpi=dpi_set)
    
             