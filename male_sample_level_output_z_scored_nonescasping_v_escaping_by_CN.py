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
import sys




"""





Manually computes p-values for PanCan due to recursion error






"""



def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

geneset_of_interest = "ESCAPE"

dpi_set = 72
sns.set(rc={'figure.figsize':(3.95,3.30)})
sns.set(font_scale=0.55)

male_valid = ['BLCA', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'PRAD']

female_valid = ['UVM', 'OV', 'SKCM', 'ACC', 'KICH', 'KIRC', 'KIRP', 'UCS', 'HNSC', 'GBM', 'COAD', 'LUSC', 'LUAD', 'SARC', 'STAD', 'LAML', 'BRCA', 'UCEC', 'CESC']

all_valid = list(set(male_valid + female_valid)) 


entire_dataset = [['TGCT-NS']]
#entire_dataset = [male_valid]
#entire_dataset = [['SKCM']]




#chrX CN df generation

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

figsize = (10,7) # the size of the figure - changes the shape of the squares in the comut
dpi = 72 # change the output resolution

x_padding = 0.01 # the x distance between patches in comut
y_padding = 0.01 # the y distance between patches in comut
tri_padding = 0.03 # the distance between triangles in comut

maf_path = '/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/comut_modified_MAF.csv' # change this to the path of your MAF, e.g. 'mutation.maf'

mutation_data = pd.read_csv(maf_path)

for a in invalid_ids:

    mutation_data = mutation_data[~(mutation_data['Tumor_Sample_Barcode'].str.contains(a))]

df_color = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/depmap/cell_lines/Modified_Input/color_list.xlsx")

mut_mapping = dict(zip(df_color['Type'].tolist(), df_color['Color'].tolist()))

type_mapping = {'Tumor':'violet', 'Normal':'khaki'}

tn_df = mutation_data[['Tumor_Sample_Barcode', 'Type_col', 'Tumor_Normal']]

tn_df['sample'] = tn_df['Tumor_Sample_Barcode']
tn_df['category'] = tn_df['Type_col']
tn_df['value'] = tn_df['Tumor_Normal']

del tn_df['Tumor_Sample_Barcode'] 
del tn_df['Type_col']
del tn_df['Tumor_Normal']

mutation_df = mutation_data[['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification']]
mutation_df['sample'] = mutation_df['Tumor_Sample_Barcode']
mutation_df['category'] = mutation_df['Hugo_Symbol']
mutation_df['value'] = mutation_df['Variant_Classification']

del mutation_df['Tumor_Sample_Barcode'] 
del mutation_df['Hugo_Symbol']
del mutation_df['Variant_Classification']

df_TGCT = pd.read_csv('/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/tgct_tcga_pan_can_atlas_2018/data_clinical_sample.txt', sep="\t")

df_TGCT = df_TGCT[df_TGCT['CANCER_TYPE']=="Non-Seminomatous Germ Cell Tumor"]

all_TGCT_ids = df_TGCT['PATIENT_ID'].tolist()
oncotree = df_TGCT['ONCOTREE_CODE'].tolist()

oncotree_dict = dict(zip(all_TGCT_ids,oncotree))

all_barcodes = mutation_df['sample'].tolist()
trunc_barcodes = []

for a in all_barcodes:
    trunc_barcodes.append(split_advanced(a, "-", 3)[0])

class_list = mutation_df['value'].tolist()

subclass_vals = []
subclass_cat = []

a=0
while a<len(trunc_barcodes):
    if class_list[a] == "NSGCT":
        try:
            subclass_vals.append(oncotree_dict[trunc_barcodes[a]])
        except KeyError:
            subclass_vals.append("UNKNOWN")
    else:
        subclass_vals.append("N/A")
    a=a+1
    subclass_cat.append("NSGCT Subtype")

subclass_df = pd.DataFrame()
subclass_df['sample'] = all_barcodes
subclass_df['category'] = subclass_cat
subclass_df['value'] = subclass_vals

df_cn = pd.read_csv('/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/Titan_CN_Xist_pos.txt', sep="\t")

df_cn = df_cn[df_cn['Sample']!="TCGA-HC-7740"] #XIST+ normal, XIST- tumor
df_cn = df_cn[df_cn['Sample']!="TCGA-19-4065"] #XIST+ tumor is -02, CN calculated for -01
#df_cn = df_cn[df_cn['Sample']!="TCGA-2G-AAGY"] #this is an XIST- sample, but is annotated as such, therefore, no problem

df_cn2 = pd.read_csv('/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/Titan_CN_Xist_neg_nsgct.txt', sep="\t", header=None)

df_cn2.columns = ['Sample', 'Unrounded_chrX_CN', 'Rounded_chrX_CN']

df_cn2 = df_cn2[df_cn2['Sample']!="TCGA-2G-AAGY"] #-05 extension, XIST- -01, solved later

df_cn = pd.concat([df_cn, df_cn2])

df_cn_sample = df_cn['Sample'].tolist()
cn_val = df_cn['Rounded_chrX_CN'].tolist()

cn_dict = dict(zip(df_cn_sample,cn_val))

tumor_normal = mutation_data['Tumor_Normal'].tolist()

new_cn = []
new_cat = []

a=0
while a<len(trunc_barcodes):
    if tumor_normal[a] == "Normal":
        new_cn.append(1)
    elif int(all_barcodes[a][-2:]) == 5:
        new_cn.append(2) #this happens to be true 
    elif trunc_barcodes[a] == "TCGA-19-4065":
        new_cn.append(2)
    elif trunc_barcodes[a] == "TCGA-HC-7740":
        new_cn.append(3)        
    else:
        try:
            new_cn.append(int(cn_dict[trunc_barcodes[a]]))
        except KeyError:
            new_cn.append("Unknown")
            print(trunc_barcodes[a])
    a=a+1
    new_cat.append("chrX Copy Number")

cn_df = pd.DataFrame()
cn_df['sample'] = all_barcodes
cn_df['category'] = new_cat
cn_df['value'] = new_cn
















      
aaaaa=0

fig, ax = plt.subplots(1, 1, sharey=True, tight_layout=True)

for all_classes in entire_dataset:
                
    combined_cn1 = []   
    combined_cn2 = []   
    combined_cn3 = []   
    combined_cn4 = []   
    combined_xistneg = []   
    
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
            
            if geneset_of_interest == "ESCAPE":
            
                #df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/PAR_removed_mean_TPM_nonescape_v_escape_no_segment_mean_normalization_total_TPM_sample_level_ratio_XISTpos_samples.csv")
                #df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/PAR_removed_mean_TPM_nonescape_v_escape_no_segment_mean_normalization_total_TPM_sample_level_ratio_XISTneg_samples.csv")
                #df3 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/PAR_removed_mean_TPM_nonescaping_vs_escape_FEMALES_no_segment_mean_normalization_total_TPM_sample_level_ratio_XISTpos_samples.csv")
                #df4 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/PAR_removed_mean_TPM_nonescaping_vs_escape_FEMALES_no_segment_mean_normalization_total_TPM_sample_level_ratio_XISTneg_samples.csv")
                       
                
                df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_XISTpos_samples.csv")
                df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_XISTneg_samples.csv")

            elif geneset_of_interest == "AUTOSOMAL":
                
                #df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/mean_TPM_nonescape_v_autosome_no_segment_mean_normalization_total_TPM_sample_level_ratio_XISTpos_samples.csv")
                #df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/mean_TPM_nonescape_v_autosome_no_segment_mean_normalization_total_TPM_sample_level_ratio_XISTneg_samples.csv")
                #df3 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/mean_TPM_nonescaping_vs_autosome_FEMALES_no_segment_mean_normalization_total_TPM_sample_level_ratio_XISTpos_samples.csv")
                #df4 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/mean_TPM_nonescaping_vs_autosome_FEMALES_no_segment_mean_normalization_total_TPM_sample_level_ratio_XISTneg_samples.csv")
                  
                df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_XISTpos_samples.csv")
                df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_XISTneg_samples.csv")

            invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

            for a in invalid_ids:
                df = df[~(df['sample'].str.contains(a))]
                df2 = df2[~(df2['sample'].str.contains(a))]
                
            #df4 = df4[df4['sample']!="TCGA-25-1870-01"]
            
            if www != "TGCT-NS":
                df = df[df['Classification']==www]
                df2 = df2[df2['Classification']==www]
                #df3 = df3[df3['Classification']==www]
                #df4 = df4[df4['Classification']==www]
            else:
                df = df[df['Secondary_Class']==www]
                df2 = df2[df2['Secondary_Class']==www]                
                #df3 = df3[df3['Secondary_Class']==www]
                #df4 = df4[df4['Secondary_Class']==www]
                
            orig_df_ids = df['sample'].tolist()
            orig_df2_ids = df2['sample'].tolist()
            
            cn1_df = cn_df[cn_df['value']==1]
            cn2_df = cn_df[cn_df['value']==2]
            cn3_df = cn_df[cn_df['value']==3]
            cn4_df = cn_df[cn_df['value']==4]
            
            cn1_ids = list(set(cn1_df['sample'].tolist()) & set(orig_df_ids))
            cn2_ids = list(set(cn2_df['sample'].tolist()) & set(orig_df_ids))
            cn3_ids = list(set(cn3_df['sample'].tolist()) & set(orig_df_ids))
            cn4_ids = list(set(cn4_df['sample'].tolist()) & set(orig_df_ids))
            
            #orig_df3_ids = df3['sample'].tolist()
            #orig_df4_ids = df4['sample'].tolist()
                        
            concat_df = pd.concat([df, df2])
            
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
                z_score = t
                ys_old.append(z_score)
                
            ys_methylation = ys_old
            
            concat_df['z_scored_methylation'] = ys_methylation
            
            cn1_df_with_z = concat_df[concat_df['sample'].isin(cn1_ids)]
            cn2_df_with_z = concat_df[concat_df['sample'].isin(cn2_ids)]
            cn3_df_with_z = concat_df[concat_df['sample'].isin(cn3_ids)]
            cn4_df_with_z = concat_df[concat_df['sample'].isin(cn4_ids)]
            df2_with_z = concat_df[concat_df['sample'].isin(orig_df2_ids)]
            
            
            cn1_vals = cn1_df_with_z['z_scored_methylation'].tolist()
            cn2_vals = cn2_df_with_z['z_scored_methylation'].tolist()
            cn3_vals = cn3_df_with_z['z_scored_methylation'].tolist()
            cn4_vals = cn4_df_with_z['z_scored_methylation'].tolist()
            xistneg_vals = df2_with_z['z_scored_methylation'].tolist()
            
            master_df = pd.DataFrame([xistneg_vals, cn1_vals, cn2_vals, cn3_vals, cn4_vals]).T
            
            if www != "TGCT-NS":
            
                master_df.columns = ['XIST- Male\nPanCan', 'XIST+ Male\nPanCan\n(chrX CN=1)', 'XIST+ Male\nPanCan\n(chrX CN=2)', 'XIST+ Male\nPanCan\n(chrX CN=3)', 'XIST+ Male\nPanCan\n(chrX CN=4)']

            if www == "TGCT-NS":

                master_df.columns = ['XIST- Male\nPanCan', 'XIST+ Male\nNSGCT\n(chrX CN=1)', 'XIST+ Male\nNSGCT\n(chrX CN=2)', 'XIST+ Male\nNSGCT\n(chrX CN=3)', 'XIST+ Male\nNSGCT\n(chrX CN=4)']
            
            for a in cn1_vals:
                combined_cn1.append(a)
            for a in cn2_vals:
                combined_cn2.append(a)
            for a in cn3_vals:
                combined_cn3.append(a)
            for a in cn4_vals:
                combined_cn4.append(a)
            for a in xistneg_vals:
                combined_xistneg.append(a)
            
            curr_color="black"
            
            sns.stripplot(data=master_df, color=curr_color, s=1.5, ax=ax, jitter=0.2, edgecolor="black", alpha=0.9, zorder=2, linewidth=0.15)
            
            ax.grid(False)
            ax.set_facecolor("white")
            ax.spines['bottom'].set_color('0')
            ax.spines['left'].set_color('0')
            
            ax.tick_params(bottom='on', left='on')
            
            #ax.set_ylabel("log$_2$(Lineage-Mean Normalized chrX Non-Escape : Escape Transcriptional Output Ratio)")
        
        aaa=aaa+1
            
    #final_df = pd.DataFrame([NSGCT_XISTpos_vals, NSGCT_XISTneg_vals, combined_XISTpos_methylation_vals, combined_XISTneg_methylation_vals, combined_XISTpos_females, combined_XISTneg_females]).T
        
    #final_df.columns = ['XIST+ Male\nNSGCT', 'XIST- Male\nNSGCT', 'XIST+ Male\nPanCan', 'XIST- Male\nPanCan', 'XIST+ Female\nPanCan', 'XIST- Female\nPanCan']
    
    
    
    final_df = pd.DataFrame([combined_xistneg, combined_cn1, combined_cn2, combined_cn3, combined_cn4]).T
    
    if www != "TGCT-NS":
    
        final_df.columns = ['XIST-\nM-Pan-Cancer' + "\n(" + str(len(combined_xistneg)) + ")", 'XIST+\nM-Pan-Cancer\nchrX CN=1' + "\n(" + str(len(combined_cn1)) + ")", 'XIST+\nM-Pan-Cancer\n(chrX CN=2)' + "\n(" + str(len(combined_cn2)) + ")", 'XIST+\nM-Pan-Cancer\n(chrX CN=3)' + "\n(" + str(len(combined_cn3)) + ")", 'XIST+\nM-Pan-Cancer\n(chrX CN=4)' + "\n(" + str(len(combined_cn4)) + ")"]

    if www == "TGCT-NS":

        final_df.columns = ['XIST-\nNSGCT' + "\n(" + str(len(combined_xistneg)) + ")", 'XIST+\nNSGCT\nchrX CN=1' + "\n(" + str(len(combined_cn1)) + ")", 'XIST+\nNSGCT\nchrX CN=2' + "\n(" + str(len(combined_cn2)) + ")", 'XIST+\nNSGCT\nchrX CN=3' + "\n(" + str(len(combined_cn3)) + ")", 'XIST+\nNSGCT\nchrX CN=4' + "\n(" + str(len(combined_cn4)) + ")"]
 
    final_df.to_csv("/Users/ananthansadagopan/Downloads/temp_output.csv", index=False)   
 
    #palette = ['#8AA29E', '#F6CA83', '#BEB2C8', '#669BBC', '#C17C74', '#A8C686']
    
    palette = ['#BFB5AF', '#ECE2D0', '#D5B9B2', '#A26769', '#582C4D']
    
    sns.boxplot(data=final_df, ax=ax, showfliers=False, color="white", linewidth=0.5, palette=palette)
    
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
    
    if geneset_of_interest == "ESCAPE":
    
        ax.set_ylabel("Normalized Non-Escapee to\nEscapee Expression (log$_2$Ratio)")

    elif geneset_of_interest == "AUTOSOMAL":

        ax.set_ylabel("ERROR")


    
    if www != "TGCT-NS":
        
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[(final_df.columns.tolist()[1], final_df.columns.tolist()[4])],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.136, fontsize=7.5)
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[(final_df.columns.tolist()[0], final_df.columns.tolist()[4])],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.196, fontsize=7.5)

        """
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[('XIST+\nM-PanCan\n(chrX CN=3)', 'XIST-\nM-PanCan')],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.10, fontsize=6)
        add_stat_annotation(ax, data=final_df,
                            box_pairs=[('XIST+\nM-PanCan\n(chrX CN=4)', 'XIST-\nM-PanCan')],
                            test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.05, fontsize=6)
        """
        
    if www == "TGCT-NS":    

        """
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[('XIST+\nM-NSGCT\n(chrX CN=1)', 'XIST-\nM-NSGCT')],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.17, fontsize=7.5)                        
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[('XIST+\nM-NSGCT\n(chrX CN=2)', 'XIST-\nM-NSGCT')],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.11, fontsize=7.5)
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[('XIST+\nM-NSGCT\n(chrX CN=3)', 'XIST-\nM-NSGCT')],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.055, fontsize=7.5)
        """
        
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[(final_df.columns.tolist()[4], final_df.columns.tolist()[0])],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0, fontsize=7.5)
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[(final_df.columns.tolist()[3], final_df.columns.tolist()[0])],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.08, fontsize=7.5)
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[(final_df.columns.tolist()[2], final_df.columns.tolist()[0])],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.16, fontsize=7.5)
        add_stat_annotation(ax, data=final_df,
                        box_pairs=[(final_df.columns.tolist()[1], final_df.columns.tolist()[0])],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.24, fontsize=7.5)
    #ax.set_ylim([-2, 4])
    
    aaaaa = aaaaa + 1

for line in ax.get_lines():
    line.set_color('black')
    line.set_linewidth(0.5)

ax.xaxis.labelpad = 1
ax.yaxis.labelpad = 1
ax.spines['bottom'].set_lw(0.5)
ax.spines['left'].set_lw(0.5)
ax.xaxis.set_tick_params(width=0.5, size=3, pad=2)
ax.yaxis.set_tick_params(width=0.5, size=3, pad=2)


#plt.ylim([-4.2, 9.8])
#plt.ylim([0, 6])




"""


https://www.statskingdom.com/170median_mann_whitney.html


"""


if www != "TGCT-NS":

    plt.ylim([-1, 2])    

    x1, x2 = 2, 4   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    y, h, col = 1.7, 0.06, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.5, c="black")
    plt.text((x1+x2)*.5, y+h, "6.16e-01", ha='center', va='bottom', color="black", fontsize=7.5)

    x1, x2 = 3, 4   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    y, h, col = 1.55, 0.06, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.5, c="black")
    plt.text((x1+x2)*.5, y+h, "2.89e-01", ha='center', va='bottom', color="black", fontsize=7.5)


plt.tight_layout()

if geneset_of_interest == "ESCAPE":
    
    if www != "TGCT-NS":    
    
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PanCan_nonescaping_v_escaping_no_segment_normalization_male_and_female_all_categories.png", dpi=dpi_set)
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PanCan_nonescaping_v_escaping_no_segment_normalization_male_and_female_all_categories.pdf", dpi=dpi_set)
    
    if www == "TGCT-NS":    
        
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/NSGCT_nonescaping_v_escaping_no_segment_normalization_male_and_female_all_categories.png", dpi=dpi_set)
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/NSGCT_nonescaping_v_escaping_no_segment_normalization_male_and_female_all_categories.pdf", dpi=dpi_set)
        
        
elif geneset_of_interest == "AUTOSOMAL":

    if www != "TGCT-NS":     
    
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PanCan_nonescaping_v_autosomal_no_segment_normalization_male_and_female_all_categories.png", dpi=dpi_set)
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PanCan_nonescaping_v_autosomal_no_segment_normalization_male_and_female_all_categories.pdf", dpi=dpi_set)

    if www == "TGCT-NS":     
    
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/NSGCT_nonescaping_v_autosomal_no_segment_normalization_male_and_female_all_categories.png", dpi=dpi_set)
        plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/NSGCT_nonescaping_v_autosomal_no_segment_normalization_male_and_female_all_categories.pdf", dpi=dpi_set)
                
                 