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
import numpy as np

plt.rcParams["font.family"] = "Arial"

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

geneset_of_interest = "ESCAPE"

log_scale = True

dpi_set = 72
sns.set(rc={'figure.figsize':(7,6)})
sns.set(font_scale=1.6)

male_valid = ['BLCA', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'PRAD']

female_valid = ['UVM', 'OV', 'SKCM', 'ACC', 'KICH', 'KIRC', 'KIRP', 'UCS', 'HNSC', 'GBM', 'COAD', 'LUSC', 'LUAD', 'SARC', 'STAD', 'LAML', 'BRCA', 'UCEC', 'CESC']

print(list(set(male_valid) & set(female_valid)))

NSGCT_valid = ['TGCT-NS']

#SGCT_valid = ['TGCT-S']


#list_of_list = [male_valid, female_valid, NSGCT_valid]
list_of_list = [male_valid, NSGCT_valid]








#CN dict generation

#replace ABSOLUTE purity with TITAN purity


#MISSING CN data:
#TCGA-B0-4846, TCGA-CJ-4642, TCGA-CZ-4862, TCGA-B0-4696 - KIRC
#TCGA-CV-7428 - HNSC
#TCGA-AB-2872 - LAML missing - justified (no data confirmed)


def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

invalid_ids = invalid_ids + ['TCGA-HC-7740-01', 'TCGA-HC-7740-11', 'TCGA-GN-A4U8-06', 'TCGA-GN-A4U8-11']

plt.rcParams["font.family"] = "Arial"

hfont = {'fontname':'Arial'}
sns.set(font_scale=1)

plt.rcParams['axes.labelsize'] = 10

figsize = (10,6) # the size of the figure - changes the shape of the squares in the comut
dpi = 72 # change the output resolution

x_padding = 0.01 # the x distance between patches in comut
y_padding = 0.01 # the y distance between patches in comut
tri_padding = 0.03 # the distance between triangles in comut

maf_path = '/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/PanCan_full_MAF.csv' # change this to the path of your MAF, e.g. 'mutation.maf'

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


df4 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/Titan_CN_Xist_neg_pan_can.txt", sep="\t", header=None)

df4.columns = ['Sample', 'Unrounded_chrX_CN', 'Rounded_chrX_CN']

df4 = df4[df4['Sample']!="TCGA-AB-2940"] #XIST+ tumor is -03, CN calculated for -01 -- this is a guess, must confirm

cn_bar = df4['Sample'].tolist()
cn_vals = df4['Rounded_chrX_CN'].tolist()

cn_dict_rev = dict(zip(cn_bar, cn_vals))







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
            try:
                new_cn.append(cn_dict_rev[trunc_barcodes[a]])
            except KeyError:
                new_cn.append("Unknown")
            print(trunc_barcodes[a])
    a=a+1
    new_cat.append("chrX Copy Number")

cn_df = pd.DataFrame()
cn_df['sample'] = all_barcodes
cn_df['category'] = new_cat
cn_df['value'] = new_cn


temp_ids_cn1 = cn_df['sample'].tolist()
temp_ids_cn2 = cn_df['value'].tolist()











#NSGCT CN


#replace ABSOLUTE purity with TITAN purity


#MISSING CN data:
#TCGA-B0-4846, TCGA-CJ-4642, TCGA-CZ-4862, TCGA-B0-4696 - KIRC
#TCGA-CV-7428 - HNSC
#TCGA-AB-2872 - LAML missing - justified (no data confirmed)


def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

plt.rcParams["font.family"] = "Arial"

hfont = {'fontname':'Arial'}
sns.set(font_scale=1)

plt.rcParams['axes.labelsize'] = 10

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
            new_cn.append(cn_dict[trunc_barcodes[a]])
        except KeyError:
            new_cn.append("Unknown")
            print(trunc_barcodes[a])
    a=a+1
    new_cat.append("chrX Copy Number")

cn_df = pd.DataFrame()
cn_df['sample'] = all_barcodes
cn_df['category'] = new_cat
cn_df['value'] = new_cn

cn_dict = dict(zip(temp_ids_cn1+all_barcodes, temp_ids_cn2+new_cn))





















for ww in list_of_list:
    
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
        
    if ww == ['TGCT-NS']:
        df = df[df['Secondary_Class'].isin(ww)]
        df2 = df2[df2['Secondary_Class'].isin(ww)]
    elif 'PRAD' in ww:
        df = df[df['Classification'].isin(ww)]
        df2 = df2[df2['Classification'].isin(ww)]
    else:
        df = df3[df3['Classification'].isin(ww)]
        df2 = df4[df4['Classification'].isin(ww)]        


    df_male_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_averaged_TCGA_Xena_TPM.csv")
    df_female_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_tumors_averaged_TCGA_Xena_TPM.csv")

    XIST_exp = df_male_ref['XIST_TPM'].tolist() + df_female_ref['XIST_TPM'].tolist()
    
    
    if log_scale:
        XIST_exp = [math.log(x+1, 2) for x in XIST_exp]
    
    barcode = df_male_ref['Barcode'].tolist() + df_female_ref['Barcode'].tolist()

    XIST_dict = dict(zip(barcode, XIST_exp))

    df_final_exp = []
    df2_final_exp = []
    
    for a in df['sample'].tolist():
        df_final_exp.append(XIST_dict[a])

    for a in df2['sample'].tolist():
        df2_final_exp.append(XIST_dict[a])

    test = ['df1', 'df2']
        
    for q in test:
        if q == 'df1' and geneset_of_interest == "ESCAPE":
            xs = df_final_exp
            ys = df['ratio_ne_e'].tolist()
            sample_list = df['sample'].tolist()
        if q == 'df1' and geneset_of_interest == "AUTOSOMAL":
            xs = df_final_exp
            ys = df['ratio_ne_a'].tolist()  
            sample_list = df['sample'].tolist()
        if q == 'df2' and geneset_of_interest == "ESCAPE":
            xs = df2_final_exp
            ys = df2['ratio_ne_e'].tolist()
            sample_list = df2['sample'].tolist()
        if q == 'df2' and geneset_of_interest == "AUTOSOMAL":
            xs = df2_final_exp
            ys = df2['ratio_ne_a'].tolist()  
            sample_list = df2['sample'].tolist()
            
        new_cn = []
        
        for a in sample_list:
            try:
                new_cn.append(int(cn_dict[a]))
            except KeyError:
                new_cn.append("Unknown")
                    
        cn_mapping = {1: "#84A9C0", 2: "#FFD6AF", 3: '#5F6062', 4:"#C0DA74", "Unknown": "gainsboro"}

        final_colors = []

        for a in new_cn:
            final_colors.append(cn_mapping[a])
            
            
        temp1, temp2 = scipy.stats.spearmanr(xs, ys)
        print(temp1)
        print(temp2)
            
        plt.scatter(xs,ys, c=final_colors, alpha=1, s=30, edgecolor="black", lw=0.2)
        zs = np.polyfit(xs, ys, 1)
        p = np.poly1d(zs)
        #plt.plot(xs,p(xs),"r--")
        
        plt.plot(xs,ys, linestyle='None')
        zs = np.polyfit(xs, ys, 1)
        y_hat = np.poly1d(zs)(xs)
        plt.plot(xs, y_hat, "r-", lw=1)
        
        ax = plt.gca()
        
        ax.set_facecolor("white")
        ax.grid(False)
        ax.tick_params(axis='x', which='major', bottom=True, labelsize=11, size=3)
        ax.tick_params(axis='y', which='major', left=True, labelsize=11, size=3)
        
        ax.spines['bottom'].set_color('0')
        ax.spines['left'].set_color('0')
        
        r = np.corrcoef(xs, ys)
        
        text = "ρ" + " (" + str(len(xs)-2) + ")" + " = " + "{:.3f}".format(temp1)
        text2 = "p = " + "{:.2e}".format(temp2)

        plt.gca().text(0.72, 0.98, text, transform=plt.gca().transAxes,
             fontsize=15, verticalalignment='top')        

        plt.gca().text(0.76, 0.92, text2, transform=plt.gca().transAxes,
             fontsize=15, verticalalignment='top')
        
        plt.xlabel("XIST Expression (TPM)", size=15)

        if log_scale:
            plt.xlabel("XIST Expression (log$_2$(TPM+1))", size=15) 
        
        if geneset_of_interest == "ESCAPE":
        
            plt.ylabel("Normalized Non-Escapee to\nEscapee Expression (log$_2$Ratio)", size=15)

        if geneset_of_interest == "AUTOSOMAL":

            plt.ylabel("Normalized Non-Escapee to\nAutosomal Gene Expression (log$_2$Ratio)", size=15)

        for line in ax.get_lines():
            line.set_color('black')
            line.set_linewidth(0.5)
          

        ax.xaxis.labelpad = 5
        ax.yaxis.labelpad = 5
        ax.spines['bottom'].set_lw(0.5)
        ax.spines['left'].set_lw(0.5)
        ax.xaxis.set_tick_params(width=0.5, size=3, pad=2)
        ax.yaxis.set_tick_params(width=0.5, size=3, pad=2)
        
        cn1_patch = mpatches.Patch(color=cn_mapping[1], label='1')
        cn2_patch = mpatches.Patch(color=cn_mapping[2], label='2')
        cn3_patch = mpatches.Patch(color=cn_mapping[3], label='3')
        cn4_patch = mpatches.Patch(color=cn_mapping[4], label='4')
        cnunk_patch = mpatches.Patch(color=cn_mapping["Unknown"], label='Unknown')
        
        leg = ax.legend(handles=[cn1_patch, cn2_patch, cn3_patch, cn4_patch, cnunk_patch], facecolor='white', title="chrX CN", fontsize=10, loc=(1.05, 0.75))
        
        frame = leg.get_frame()
        frame.set_edgecolor("black")
        frame.set_linewidth(0.2)

        #plt.ylim([-4.2, 9.8])
        #plt.ylim([0, 6])
        
        plt.tight_layout()
        
        if ww == ['TGCT-NS']:
            modifier = 'NSGCT_'
        elif 'PRAD' in ww:
            modifier = 'M_PanCan_'
        else:
            modifier = 'F_PanCan_'   
            
        if log_scale:
            modifier = modifier + "log2_"
                    
        if q == 'df1' and geneset_of_interest == "ESCAPE":
            plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XISTpos_" + modifier + "scatter_XIST_expression_v_NE_to_E_ratio.png", dpi=dpi_set)
            plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XISTpos_" + modifier + "scatter_XIST_expression_v_NE_to_E_ratio.pdf", dpi=dpi_set)
        if q == 'df1' and geneset_of_interest == "AUTOSOMAL":
            plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XISTpos_" + modifier + "scatter_XIST_expression_v_NE_to_A_ratio.png", dpi=dpi_set)
            plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XISTpos_" + modifier + "scatter_XIST_expression_v_NE_to_A_ratio.pdf", dpi=dpi_set)
        if q == 'df2' and geneset_of_interest == "ESCAPE":
            pass
            #plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XISTneg_" +  modifier + "scatter_XIST_expression_v_NE_to_E_ratio.png", dpi=dpi_set)
            #plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XISTneg_" + modifier + "scatter_XIST_expression_v_NE_to_E_ratio.pdf", dpi=dpi_set)
        if q == 'df2' and geneset_of_interest == "AUTOSOMAL":
            pass
            #plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XISTneg_" + modifier + "scatter_XIST_expression_v_NE_to_A_ratio.png", dpi=dpi_set)
            #plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XISTneg_" + modifier + "scatter_XIST_expression_v_NE_to_A_ratio.pdf", dpi=dpi_set)
    
        plt.clf()

             