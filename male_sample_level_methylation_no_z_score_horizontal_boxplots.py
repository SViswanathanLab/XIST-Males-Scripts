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
from matplotlib.ticker import StrMethodFormatter, NullFormatter, ScalarFormatter, FormatStrFormatter
import matplotlib

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

interest = "X"

main_size = 7
dpi_set = 72
sns.set(rc={'figure.figsize':(3.1849,5.46)})
sns.set(font_scale=0.45)

if interest == "AUTOSOME":
    sns.set(rc={'figure.figsize':(4.1,9)})    
    sns.set(font_scale=0.8)
    main_size = 9

male_valid = ['BLCA', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'PRAD']

female_valid = ['UVM', 'OV', 'SKCM', 'ACC', 'KICH', 'KIRC', 'KIRP', 'UCS', 'HNSC', 'GBM', 'COAD', 'LUSC', 'LUAD', 'SARC', 'STAD', 'LAML', 'BRCA', 'UCEC', 'CESC']

all_valid = list(set(male_valid + female_valid)) 

entire_dataset = ['TGCT-NS', 'TGCT-S']

for a in all_valid:
    entire_dataset.append(a)

entire_dataset = [entire_dataset]

#Removed PRAD because we are looking for overlapping in male / female

#entire_dataset = [['BLCA'], ['ESCA'], ['HNSC'], ['KICH'], ['KIRC'], ['LAML'], ['LGG'], ['LIHC'], ['LUAD'], ['LUSC'], ['PAAD'], ['PCPG'], ['PRAD'], ['READ'], ['SARC'], ['SKCM'], ['STAD'], ['THCA']]
                  
aaaaa=0

fig, ax = plt.subplots(1, 1, sharex=True, tight_layout=True)

for all_classes in entire_dataset:
                
    combined_XISTpos_methylation_vals = []
    combined_XISTneg_methylation_vals = []
    combined_XISTpos_females = []
    combined_XISTneg_females = []
    
    NSGCT_XISTpos_vals = []
    NSGCT_XISTneg_vals = []   
    
    SGCT_XISTpos_vals = []
    
    """
    df3 = pd.read_excel('../Modified_Input/color_list.xlsx')
    type_color = df3['Type'].tolist()
    color_color = df3['Color'].tolist()
    
    color_dict = dict(zip(type_color, color_color))
    
    new_color = []
    
    for a in all_classes:
        new_color.append(color_dict[a])
    """
    
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
            
            if www != "TGCT-NS" and www != "TGCT-S":
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
            
            if interest == "AUTOSOME":
            
                all_methylation = concat_df['autosome_mean_beta'].tolist()
                
            elif interest == "X":
                
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
            
            if www != "TGCT-NS" and www != "TGCT-S":            
                df3_with_z = concat_df[concat_df['sample'].isin(orig_df3_ids)]
                df4_with_z = concat_df[concat_df['sample'].isin(orig_df4_ids)]
            
            vals1 = df_with_z['z_scored_methylation'].tolist()
            vals2 = df2_with_z['z_scored_methylation'].tolist()
        
            if www != "TGCT-NS" and www != "TGCT-S":
                vals3 = df3_with_z['z_scored_methylation'].tolist()
                vals4 = df4_with_z['z_scored_methylation'].tolist()
            
            
            """
            if www != "TGCT-NS":
                
                if www not in male_valid:
                    master_df = pd.DataFrame([[float("nan")], [float("nan")], [float("nan")], [float("nan")], vals3, vals4]).T
                    master_df.columns = ['XIST+ Male\nNSGCT', 'XIST- Male\nNSGCT', 'XIST+ Male\nPanCan', 'XIST- Male\nPanCan', 'XIST+ Female\nPanCan', 'XIST- Female\nPanCan']
                elif www not in female_valid:
                    master_df = pd.DataFrame([[float("nan")], [float("nan")], vals2, vals1, [float("nan")], [float("nan")]]).T
                    master_df.columns = ['XIST+ Male\nNSGCT', 'XIST- Male\nNSGCT', 'XIST+ Male\nPanCan', 'XIST- Male\nPanCan', 'XIST+ Female\nPanCan', 'XIST- Female\nPanCan']
                else:              
                    master_df = pd.DataFrame([[float("nan")], [float("nan")], vals2, vals1, vals3, vals4]).T
                    master_df.columns = ['XIST+ Male\nNSGCT', 'XIST- Male\nNSGCT', 'XIST+ Male\nPanCan', 'XIST- Male\nPanCan', 'XIST+ Female\nPanCan', 'XIST- Female\nPanCan']
            
            else:
                master_df = pd.DataFrame([vals2, vals1, [float("nan")], [float("nan")], [float("nan")], [float("nan")]]).T
                master_df.columns = ['XIST+ Male\nNSGCT', 'XIST- Male\nNSGCT', 'XIST+ Male\nPanCan', 'XIST- Male\nPanCan', 'XIST+ Female\nPanCan', 'XIST- Female\nPanCan']
            """
            
            if www == "TGCT-NS":
                for abc in vals1:
                    NSGCT_XISTpos_vals.append(abc)
                for abc in vals2:             
                    NSGCT_XISTneg_vals.append(abc)

            if www == "TGCT-S":
                for abc in vals1:
                    SGCT_XISTpos_vals.append(abc)
                    
            if www != "TGCT-NS" and www != "TGCT-S":
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
            
            #sns.stripplot(data=master_df, color=curr_color, s=3, ax=ax, jitter=0.2, edgecolor="black", alpha=0.9, zorder=2, linewidth=0.15)
            
            #ax.set_ylabel("log$_2$(Lineage-Mean Normalized chrX Non-Escape : Escape Transcriptional Output Ratio)")
        
        aaa=aaa+1
            
    #final_df = pd.DataFrame([NSGCT_XISTpos_vals, NSGCT_XISTneg_vals, combined_XISTpos_methylation_vals, combined_XISTneg_methylation_vals, combined_XISTpos_females, combined_XISTneg_females]).T

    final_df = pd.DataFrame([SGCT_XISTpos_vals, NSGCT_XISTneg_vals, NSGCT_XISTpos_vals, combined_XISTneg_methylation_vals, combined_XISTpos_methylation_vals, combined_XISTpos_females, combined_XISTneg_females]).T
        
    final_df.columns = ['XIST+ SGCT (' + str(len(SGCT_XISTpos_vals)) + ")", 'XIST- NSGCT (' + str(len(NSGCT_XISTneg_vals)) + ")", 'XIST+ NSGCT (' + str(len(NSGCT_XISTpos_vals)) + ")", 'XIST- Male\nPan-Cancer (' + str(len(combined_XISTneg_methylation_vals)) + ")", 'XIST+ Male\nPan-Cancer (' + str(len(combined_XISTpos_methylation_vals)) + ")", 'XIST+ Female\nPan-Cancer (' + str(len(combined_XISTpos_females)) + ")", 'XIST- Female\nPan-Cancer (' + str(len(combined_XISTneg_females)) + ")"]
    #palette = ['#8AA29E', '#F6CA83', '#BEB2C8', '#669BBC', '#C17C74', '#A8C686']

    palette1 = ['#BEB2C8', '#A8C686']
    palette2 = ['#669BBC', '#C17C74']
    palette3 = ['#8AA29E', '#F6CA83']

    palette_tot = palette1 + palette2 + palette3

    sns.stripplot(data=final_df, orient="h", color=curr_color, s=1, ax=ax, jitter=0.2, edgecolor="black", alpha=0.9, zorder=2, linewidth=0.15)
        
    sns.boxplot(data=final_df, orient="h", ax=ax, showfliers=False, color="white", linewidth=0.5, palette=palette_tot)
    
    q=0
    while q<1:
        for i,artist in enumerate(ax.artists):
            # Set the linecolor on the artist to the facecolor, and set the facecolor to None
            col = "black"
            artist.set_edgecolor(col)
            #artist.set_facecolor('None')
            j=0
            while j<6:
                line = ax.lines[j]
                line.set_color(col)
                line.set_mfc(col)
                line.set_mec(col)
                j=j+1
        q=q+1
    
    """
    ax.set_xticks([1])
    
    if len(all_classes) > 1:
    
        ax.set_xticklabels(["PanCan"], fontsize='x-large')
        
    elif all_classes[0] == "TGCT-NS":
        
        ax.set_xticklabels(["NSGCT"], fontsize='x-large')      
        
    elif all_classes[0] == "TGCT-S":
        
        ax.set_xticklabels(["SGCT"], fontsize='x-large') 
        
    elif all_classes[0] == "Cancer-Adjacent Normal":
        
        ax.set_xticklabels(["Cancer-Adjacent Normal"], fontsize='x-large')   
        #all_valid_classes
        patch_list = []
        for n in all_valid_classes:
            test_patch = mpatches.Patch(color=color_dict[n], label=n)
            patch_list.append(test_patch)
    
        plt.legend(handles=patch_list, bbox_to_anchor=(1.04, 0.72), fontsize=10, facecolor='white', framealpha=1)
    """
    
    
    if interest == "AUTOSOME":
    
        ax.set_xlabel("Median Autosome Methylation in Sample (β)", fontsize=main_size)
        
    elif interest == "X":
        
        ax.set_xlabel("Median chrX Methylation at NE Promoter CGIs in Sample (β)", fontsize=main_size)      
      
    """
    add_stat_annotation(ax, data=final_df,
                    box_pairs=[('XIST+', 'XIST-')],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.01, fontsize=7)
    add_stat_annotation(ax, data=final_df,
                            box_pairs=[('XIST+', 'XIST-')],
                            test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.01, fontsize=7)
    add_stat_annotation(ax, data=final_df,
                    box_pairs=[('XIST+', 'XIST-')],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.01, fontsize=7)
    """
    
    """

    add_stat_annotation(ax, data=final_df,
                    box_pairs=[('XIST- Male\nNSGCT', 'XIST+ Male\nNSGCT')],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.2, fontsize=7.5)    
    add_stat_annotation(ax, data=final_df,
                    box_pairs=[('XIST- Female\nPanCan', 'XIST+ Female\nPanCan')],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.2, fontsize=7.5)
    add_stat_annotation(ax, data=final_df,
                    box_pairs=[('XIST- Male\nPanCan', 'XIST+ Male\nPanCan')],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.1, fontsize=7.5)
    """
    """
    #ax.set_ylim([-2, 4])
    
    aaaaa = aaaaa + 1

q=0
while q<3:
    for line in axs[q].get_lines():
        line.set_color('black')
        line.set_linewidth(0.5)
    q=q+1

"""

"""
if geneset_of_interest == "ESCAPE":
    ax.set_yticks([-4, -2, 0, 2, 4, 6, 8])
elif geneset_of_interest == "AUTOSOMAL":
    ax.set_yticks([-6, -4, -2, 0, 2, 4])
"""

ax.grid(False)
ax.set_facecolor("white")
ax.spines['bottom'].set_color('0')
ax.spines['left'].set_color('0')
ax.tick_params(bottom='on', left='on', which="both")
  
q=0
while q<1:
    ax.xaxis.labelpad = 1
    ax.yaxis.labelpad = 1
    ax.spines['bottom'].set_lw(0.5)
    ax.spines['left'].set_lw(0.5)
    ax.xaxis.set_tick_params(width=0.5, size=3, pad=2)
    ax.yaxis.set_tick_params(width=0.5, size=3, pad=2)
    q=q+1

plt.setp(ax.get_yticklabels(), fontsize=main_size)
#plt.ylim([-4.2, 9.8])
#plt.ylim([0, 6])

plt.tight_layout()

cutoff = max(combined_XISTneg_methylation_vals)+0.000000000001

if interest == "X":

    ax.axvline([cutoff], c="black", lw=0.2, ls="--")
  

final_df.columns = ['XIST+ SGCT', 'XIST- NSGCT', 'XIST+ NSGCT', 'XIST- Male Pan-Cancer', 'XIST+ Male Pan-Cancer', 'XIST+ Female Pan-Cancer', 'XIST- Female Pan-Cancer']

XIST_NSGCT_vals = final_df['XIST+ NSGCT'].tolist()
XIST_PanCan_vals = final_df['XIST+ Male Pan-Cancer'].tolist()
XIST_f_vals = final_df['XIST+ Female Pan-Cancer'].tolist()


XIST_NSGCT_vals = [x for x in XIST_NSGCT_vals if x==x]
XIST_PanCan_vals = [x for x in XIST_PanCan_vals if x==x]
XIST_f_vals =  [x for x in XIST_f_vals if x==x]


counter = 0

for a in XIST_NSGCT_vals:
    if a>=cutoff:
        counter = counter +1
XISTpos_NSGCT_count = counter
XISTpos_NSGCT_tot = len(XIST_NSGCT_vals)

counter = 0

for a in XIST_PanCan_vals:
    if a>=cutoff:
        counter = counter +1
XISTpos_PanCan_count = counter
XISTpos_PanCan_tot = len(XIST_PanCan_vals)


counter = 0

for a in XIST_f_vals:
    if a>=cutoff:
        counter = counter +1
XISTpos_f_count = counter
XISTpos_f_tot = len(XIST_f_vals)


val1 = XIST_NSGCT_vals
val2 = XIST_PanCan_vals
val3 = XIST_f_vals





XIST_NSGCT_vals = final_df['XIST- NSGCT'].tolist()
XIST_PanCan_vals = final_df['XIST- Male Pan-Cancer'].tolist()
XIST_neg_f_vals = final_df['XIST- Female Pan-Cancer'].tolist()

XIST_NSGCT_vals = [x for x in XIST_NSGCT_vals if x==x]
XIST_PanCan_vals = [x for x in XIST_PanCan_vals if x==x]
XIST_neg_f_vals = [x for x in XIST_neg_f_vals if x==x]



u, p = scipy.stats.mannwhitneyu(val1, XIST_NSGCT_vals)

print(u)
print(p)

u, p = scipy.stats.mannwhitneyu(val2, XIST_PanCan_vals)

print(u)
print(p)

u, p = scipy.stats.mannwhitneyu(val3, XIST_neg_f_vals)

print(u)
print(p)



counter = 0

for a in XIST_NSGCT_vals:
    if a>=cutoff:
        counter = counter +1
XISTneg_NSGCT_count = counter
XISTneg_NSGCT_tot = len(XIST_NSGCT_vals)

counter = 0

for a in XIST_PanCan_vals:
    if a>=cutoff:
        counter = counter +1
XISTneg_PanCan_count = counter
XISTneg_PanCan_tot = len(XIST_PanCan_vals)

counter = 0

for a in XIST_neg_f_vals:
    if a>=cutoff:
        counter = counter +1
XIST_neg_f_vals_count = counter
XIST_neg_f_vals_tot = len(XIST_neg_f_vals)





fisher_test_list_NSGCT = [[XISTpos_NSGCT_count, XISTpos_NSGCT_tot-XISTpos_NSGCT_count], [XISTneg_NSGCT_count, XISTneg_NSGCT_tot-XISTneg_NSGCT_count]] 

o, p = scipy.stats.fisher_exact(fisher_test_list_NSGCT, alternative='two-sided')

print(o)
print(p)

fisher_test_list_PanCan = [[XISTpos_PanCan_count, XISTpos_PanCan_tot-XISTpos_PanCan_count], [XISTneg_PanCan_count, XISTneg_PanCan_tot-XISTneg_PanCan_count]] 

o, p = scipy.stats.fisher_exact(fisher_test_list_PanCan, alternative='two-sided')

print(o)
print(p)




fisher_test_list_PanCan = [[XISTpos_f_count, XISTpos_f_tot-XISTpos_f_count], [XIST_neg_f_vals_count, XIST_neg_f_vals_tot-XIST_neg_f_vals_count]] 

print(fisher_test_list_PanCan)

o, p = scipy.stats.fisher_exact(fisher_test_list_PanCan, alternative='two-sided')

print(o)
print(p)




if interest == "X":

    plt.xlim([0, 0.62])

if interest == "AUTOSOME":
    
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/autosomal_methylation_male_and_female_all_categories_horizontal_boxplots.png", dpi=dpi_set)
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/autosomal_methylation_male_and_female_all_categories_horizontal_boxplots.pdf", dpi=dpi_set)
    
elif interest == "X":
    
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/chrX_methylation_male_and_female_all_categories_horizontal_boxplots.png", dpi=dpi_set)
    plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/chrX_methylation_male_and_female_all_categories_horizontal_boxplots.pdf", dpi=dpi_set)
        
         