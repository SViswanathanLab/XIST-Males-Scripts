#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 17:38:22 2021

@author: ananthansadagopan
"""

from comut import comut
from comut import fileparsers
import palettable
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib
import math


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

cn_df.to_csv("/Users/ananthansadagopan/Downloads/temp.csv", index=False)
















































df_transcriptional_output = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_averaging_median_MALES_no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_all_samples.csv")


df_ne_a = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/no_averaging_median_MALES_no_PAR_no_segment_mean_normalization_median_ne_e_autosome_pseudo_normalized_sample_level_ratio_all_samples.csv")

nea_dict = dict(zip(df_ne_a['sample'].tolist(), df_ne_a['ratio_ne_a'].tolist()))

output_dict = dict(zip(df_transcriptional_output['sample'].tolist(), df_transcriptional_output['ratio_ne_e'].tolist()))



df_methylation = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/new_for_comut_sample_level_methylation_all_samples.csv")

methylation_dict = dict(zip(df_methylation['sample'].tolist(), df_methylation['chrX_non_escaping_promoter_CpG_islands_mean_beta'].tolist()))

auto_dict = dict(zip(df_methylation['sample'].tolist(), df_methylation['autosome_mean_beta'].tolist()))




all_barcodes = mutation_df['sample'].tolist()

methylation_vals = []
output_vals = []
nea_vals = []
auto_vals = []

methylation_cat = []
output_cat = []
nea_cat = []
auto_cat = []

for a in all_barcodes:
    
    try:
        methylation_vals.append(methylation_dict[a])
    except KeyError:
        methylation_vals.append("UNKNOWN")

    try:
        output_vals.append(math.log(output_dict[a], 2))
    except KeyError:
        output_vals.append("UNKNOWN")

    try:
        nea_vals.append(math.log(nea_dict[a], 2))
    except KeyError:
        nea_vals.append("UNKNOWN")

    try:
        auto_vals.append(auto_dict[a])
    except KeyError:
        auto_vals.append("UNKNOWN")
    
    methylation_cat.append("chrX NE Promoter\nCGI Methylation")
    output_cat.append("log$_2$(NE:E Ratio)")
    nea_cat.append("log$_2$(NE:A Ratio)")
    auto_cat.append("Autosome Methylation")

nea_df = pd.DataFrame()
nea_df['sample'] = all_barcodes
nea_df['category'] = nea_cat
nea_df['value'] = nea_vals

auto_df = pd.DataFrame()
auto_df['sample'] = all_barcodes
auto_df['category'] = auto_cat
auto_df['value'] = auto_vals

methylation_df = pd.DataFrame()
output_df = pd.DataFrame()
methylation_df['sample'] = all_barcodes
output_df['sample'] = all_barcodes

methylation_df['category'] = methylation_cat
output_df['category'] = output_cat

methylation_df['value'] = methylation_vals
output_df['value'] = output_vals

df_new = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.processed.txt", sep="\t")

uq_samples = df_new['array'].tolist()

chrX_cn = df_new['purity'].tolist()

purity_dict = dict(zip(uq_samples, chrX_cn))


purity_vals = []

for a in all_barcodes:
    try:
        purity_vals.append(purity_dict[a])
    except KeyError:
        purity_vals.append("UNKNOWN")
        
purity_df = pd.DataFrame()

purity_df['sample'] = all_barcodes

purity_cat_val = []

for a in all_barcodes:
    purity_cat_val.append("Purity")

purity_df['category'] = purity_cat_val

purity_df['value'] = purity_vals

exp_df = mutation_data[['Tumor_Sample_Barcode', 'XIST_log2TPM_plus_1']]
exp_df['sample'] = exp_df['Tumor_Sample_Barcode']
del exp_df['Tumor_Sample_Barcode'] 
exp_df['value'] = exp_df['XIST_log2TPM_plus_1']
del exp_df['XIST_log2TPM_plus_1']

exp_category = []

secondary_val = []

second_category = []

a=0
while a<len(exp_df['value'].tolist()):
    #exp_category.append("log$_2$(XIST TPM+1)")
    exp_category.append("XIST Expression")
    second_category.append("XIST Status")
    if exp_df['value'].tolist()[a] >= 2:
        secondary_val.append("Positive")
    else:
        secondary_val.append("Negative")
    a=a+1

exp_df['category'] = exp_category

second_df = pd.DataFrame()

second_df['sample'] = exp_df['sample']
second_df['value'] = secondary_val
second_df['category'] = second_category

second_df = second_df[['sample', 'category', 'value']]

mutation_data['Tumor_Sample_Barcode'] = mutation_data['Tumor_Sample_Barcode']

example_comut = comut.CoMut()


# add data to the CoMut object
boolean_mapping = {'Positive': "red", 'Negative': "blue"}

subclass_mapping = {'MGCT': "#D00000", 'EMBCA': "#3A86FF", 'TYST': "#1C3144", 'TT': "#FFC857", 'Unknown': "#1D7874", 'N/A': "gainsboro"}

cn_mapping = {1: "#84A9C0", 2: "#FFD6AF", 3: '#5F6062', 4:"#C0DA74", "Unknown": "gainsboro"}

indicator_df = mutation_data[['Tumor_Sample_Barcode', 'Group']]
indicator_df['sample'] = indicator_df['Tumor_Sample_Barcode']
indicator_df['group'] = indicator_df['Group']

del indicator_df['Group']
del indicator_df['Tumor_Sample_Barcode']

indicator_kwargs = {'color': 'black', 'marker': 'o', 'linewidth': 0.4, 'markersize': 0.2}

example_comut.add_sample_indicators(indicator_df, name = 'Same Patient', plot_kwargs = indicator_kwargs)

temp = [x for x in output_df['value'].tolist() if x!="UNKNOWN"]

print(min(temp))

print(max(temp))


#example_comut.add_continuous_data(auto_df, name = 'Autosome Methylation', mapping = 'Greys', value_range = (0.03, 0.11))

example_comut.add_continuous_data(methylation_df, name = 'chrX NE Promoter\nCGI Methylation', mapping = 'Blues', value_range = (0, 0.41))

#example_comut.add_continuous_data(nea_df, name = 'NE:A Ratio', mapping = 'Greens', value_range = (-1.01, 2))

#example_comut.add_continuous_data(output_df, name = 'NE:E Ratio', mapping = 'Reds', value_range = (-0.7,2.1))

example_comut.add_continuous_data(purity_df, name = 'Purity', mapping = 'Purples', value_range = (0, 1))

example_comut.add_categorical_data(cn_df, name='chrX Copy Number', mapping=cn_mapping)

example_comut.add_categorical_data(second_df, name='XIST Status', mapping=boolean_mapping)

example_comut.add_continuous_data(exp_df, name = 'log$_2$(XIST TPM+1)', mapping = 'Oranges', value_range = (0, 10))

example_comut.add_categorical_data(tn_df, name = 'Tumor/Normal', mapping = type_mapping)

#example_comut.add_categorical_data(subclass_df, name = 'NSGCT Subtype', mapping = subclass_mapping)

example_comut.add_categorical_data(mutation_df, name = 'Cancer Type', mapping = mut_mapping)


second_df = second_df.reset_index()
del second_df['index']
exp_df = exp_df.reset_index()
del exp_df['index']
tn_df = tn_df.reset_index()
del tn_df['index']
mutation_df = mutation_df.reset_index()
del mutation_df['index']

example_comut.plot_comut(x_padding = x_padding, y_padding = y_padding, hspace = 0.05, wspace=0.01, tri_padding = tri_padding, figsize = figsize, widths = [5, 0.5])

legend = example_comut.add_unified_legend(bbox_to_anchor=(1.02,1), ncol=2, fontsize=7)

full_output_df = pd.concat([auto_df, methylation_df, nea_df, output_df, purity_df, cn_df, second_df, exp_df, tn_df, subclass_df, mutation_df], axis=1)

temp_cols = full_output_df['sample']

temp_cols.columns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

tempref_val = temp_cols[0].tolist()

a=0
while a<11:
    if temp_cols[a].tolist() != tempref_val:
        print(str(a) + " ERROR")
    tempref_val = temp_cols[a].tolist()
    a=a+1

full_output_df.to_csv("/Users/ananthansadagopan/Downloads/output_summary.csv", index=False)









font_size = 6

fig = example_comut.figure

# color bars must be added manually based on figure coordinates - [left, bottom, width, height]
purity_ax = fig.add_axes([0.93, 0.4, 0.06, 0.01])

# purity ranges from 0 to 1
norm = matplotlib.colors.Normalize(vmin=0, vmax=1)

# create the colorbar with colormap used when the continuous data was added (purp_7)
purity_colorbar = example_comut.figure.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap="Purples"),
                                                 cax=purity_ax, orientation='horizontal')


purity_colorbar.ax.tick_params(labelsize=font_size)

# remove tick marks and move tick labels slightly inwards. Also remove black border
purity_colorbar.ax.tick_params(size=2)
purity_colorbar.set_ticks([0,0.5,1])
purity_colorbar.set_ticklabels([0,0.5,1])
purity_colorbar.outline.set_visible(False)

# set title of colorbar to line up with other legend elements
purity_colorbar.set_label('Purity', fontsize = font_size)

# color bars must be added manually based on figure coordinates - [left, bottom, width, height]
exp_ax = fig.add_axes([0.85, 0.4, 0.06, 0.01])

# purity ranges from 0 to 1
norm = matplotlib.colors.Normalize(vmin=0, vmax=10)

# create the colorbar with colormap used when the continuous data was added (purp_7)
exp_colorbar = example_comut.figure.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap="Oranges"),
                                                 cax=exp_ax, orientation='horizontal')


exp_colorbar.ax.tick_params(labelsize=font_size)

# remove tick marks and move tick labels slightly inwards. Also remove black border
exp_colorbar.ax.tick_params(size=2)
exp_colorbar.set_ticks([0,5,10])
exp_colorbar.set_ticklabels([0,5,10])
exp_colorbar.outline.set_visible(False)

# set title of colorbar to line up with other legend elements
exp_colorbar.set_label('XIST Expression\n(log$_2$(TPM+1))', fontsize = font_size)


"""
# color bars must be added manually based on figure coordinates - [left, bottom, width, height]
nee_ax = fig.add_axes([0.85, 0.29, 0.06, 0.01])

# purity ranges from 0 to 1
norm = matplotlib.colors.Normalize(vmin=-0.7, vmax=2.1)

# create the colorbar with colormap used when the continuous data was added (purp_7)
nee_colorbar = example_comut.figure.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap="Reds"),
                                                 cax=nee_ax, orientation='horizontal')

font_size = 6

nee_colorbar.ax.tick_params(labelsize=font_size)

# remove tick marks and move tick labels slightly inwards. Also remove black border
nee_colorbar.ax.tick_params(size=2)
nee_colorbar.set_ticks([-0.7,0.7,2.1])
nee_colorbar.set_ticklabels([-0.7,0.7,2.1])
nee_colorbar.outline.set_visible(False)

# set title of colorbar to line up with other legend elements
nee_colorbar.set_label('log$_2$(NE:E Ratio)', fontsize = font_size)



# color bars must be added manually based on figure coordinates - [left, bottom, width, height]
nea_ax = fig.add_axes([0.93, 0.29, 0.06, 0.01])

# purity ranges from 0 to 1
norm = matplotlib.colors.Normalize(vmin=-1.01, vmax=2)

# create the colorbar with colormap used when the continuous data was added (purp_7)
nea_colorbar = example_comut.figure.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap="Greens"),
                                                 cax=nea_ax, orientation='horizontal')

font_size = 6

nea_colorbar.ax.tick_params(labelsize=font_size)

# remove tick marks and move tick labels slightly inwards. Also remove black border
nea_colorbar.ax.tick_params(size=2)
nea_colorbar.set_ticks([-1.0,0.5,2])
nea_colorbar.set_ticklabels([-1.0,0.5,2])
nea_colorbar.outline.set_visible(False)

# set title of colorbar to line up with other legend elements
nea_colorbar.set_label('log$_2$(NE:A Ratio)', fontsize = font_size)

"""



# color bars must be added manually based on figure coordinates - [left, bottom, width, height]
neem_ax = fig.add_axes([0.85, 0.3, 0.06, 0.01])

# purity ranges from 0 to 1
norm = matplotlib.colors.Normalize(vmin=0, vmax=0.41)

# create the colorbar with colormap used when the continuous data was added (purp_7)
neem_colorbar = example_comut.figure.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap="Blues"),
                                                 cax=neem_ax, orientation='horizontal')

font_size = 6

neem_colorbar.ax.tick_params(labelsize=font_size)

# remove tick marks and move tick labels slightly inwards. Also remove black border
neem_colorbar.ax.tick_params(size=2)
neem_colorbar.set_ticks([0,0.2,0.4])
neem_colorbar.set_ticklabels([0,0.2,0.4])
neem_colorbar.outline.set_visible(False)

# set title of colorbar to line up with other legend elements
neem_colorbar.set_label('chrX NE Promoter\nCGI Methylation (β)', fontsize = font_size)



"""
# color bars must be added manually based on figure coordinates - [left, bottom, width, height]
neam_ax = fig.add_axes([0.93, 0.21, 0.06, 0.01])

# purity ranges from 0 to 1
norm = matplotlib.colors.Normalize(vmin=0.03, vmax=0.11)

# create the colorbar with colormap used when the continuous data was added (purp_7)
neam_colorbar = example_comut.figure.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap="Greys"),
                                                 cax=neam_ax, orientation='horizontal')

font_size = 6

neam_colorbar.ax.tick_params(labelsize=font_size)

# remove tick marks and move tick labels slightly inwards. Also remove black border
neam_colorbar.ax.tick_params(size=2)
neam_colorbar.set_ticks([0.03,0.07,0.11])
neam_colorbar.set_ticklabels([0.03,0.07,0.11])
neam_colorbar.outline.set_visible(False)

# set title of colorbar to line up with other legend elements
neam_colorbar.set_label('Autosome\nMethylation (β)', fontsize = font_size)
"""














example_comut.axes['Same Patient'].set_facecolor("white")
example_comut.axes['Same Patient'].grid(False)

ax = example_comut.axes['XIST Status']
ax.get_xaxis().set_ticks([])

example_comut.figure.savefig('/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PanCan_master_comut_supplement.pdf', dpi = dpi, bbox_inches = 'tight')
example_comut.figure.savefig('/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PanCan_master_comut_supplement.png', dpi = dpi, bbox_inches = 'tight')

