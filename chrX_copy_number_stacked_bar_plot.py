#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 07:37:32 2021

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

plt.rcParams["font.family"] = "Arial"
dpi_set = 72
sns.set(rc={'figure.figsize':(8,5)})
sns.set(font_scale=1)






#MUST add 1 KIRC to XIST- PanCan






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

df4 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/Titan_CN_Xist_neg_pan_can.txt", sep="\t", header=None)

df4.columns = ['Sample', 'Unrounded_chrX_CN', 'Rounded_chrX_CN']

df4 = df4[df4['Sample']!="TCGA-AB-2940"] #XIST+ tumor is -03, CN calculated for -01 -- this is a guess, must confirm

df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

temp_bar = df2['Barcode'].tolist()

end_val = []

for a in temp_bar:
    end_val.append(int(split_advanced(a ,"-", 3)[1]))
    
df2['End_Val'] = end_val

df2 = df2[df2['End_Val']<10]

df4_bar = df4['Sample'].tolist()

new_barcodes = []

for a in df4_bar:
    df_temp = df2[df2['Barcode'].str.contains(a)]
    new_barcodes.append(df_temp['Barcode'].tolist()[0])
    if (len(df_temp.index.tolist())) != 1:
        print("ERROR")

df4['Sample'] = new_barcodes

cn_df = cn_df[cn_df['value']!="Unknown"]

cn_bar = df4['Sample'].tolist() + cn_df['sample'].tolist()
cn_vals = df4['Rounded_chrX_CN'].tolist() + cn_df['value'].tolist()


cn_dict = dict(zip(cn_bar, cn_vals))

print(cn_dict)


df_NSGCT = df2[df2['Secondary_Class']=="TGCT-NS"]

XISTpos_NSGCT = df_NSGCT[df_NSGCT['XIST_TPM']>=3]
XISTneg_NSGCT = df_NSGCT[df_NSGCT['XIST_TPM']<3]

df_SGCT = df2[df2['Secondary_Class']=="TGCT-S"]

XISTpos_SGCT = df_SGCT[df_SGCT['XIST_TPM']>=3]

df_PanCan = df2[df2['Classification']!="TGCT"]

XISTpos_PanCan = df_PanCan[df_PanCan['XIST_TPM']>=3]
XISTneg_PanCan = df_PanCan[df_PanCan['XIST_TPM']<3]


XISTpos_NSGCT_bar = []
XISTpos_NSGCT_cn = []

for a in XISTpos_NSGCT['Barcode'].tolist():
    try:
        XISTpos_NSGCT_cn.append(cn_dict[a])
        XISTpos_NSGCT_bar.append(a)
    except KeyError:
        continue

XISTneg_NSGCT_cn = []
XISTneg_NSGCT_bar = []

for a in XISTneg_NSGCT['Barcode'].tolist():
    try:
        XISTneg_NSGCT_cn.append(cn_dict[a])
        XISTneg_NSGCT_bar.append(a)
    except KeyError:
        continue

XISTpos_SGCT_cn = []
XISTpos_SGCT_bar = []

for a in XISTpos_SGCT['Barcode'].tolist():
    try:
        XISTpos_SGCT_cn.append(cn_dict[a])
        XISTpos_SGCT_bar.append(a)
    except KeyError:
        continue
    
XISTpos_PanCan_cn = []
XISTpos_PanCan_bar = []

for a in XISTpos_PanCan['Barcode'].tolist():
    try:
        XISTpos_PanCan_cn.append(cn_dict[a])
        XISTpos_PanCan_bar.append(a)
    except KeyError:
        continue

XISTneg_PanCan_cn = []
XISTneg_PanCan_bar = []

for a in XISTneg_PanCan['Barcode'].tolist():
    try:
        XISTneg_PanCan_cn.append(cn_dict[a])
        XISTneg_PanCan_bar.append(a)
    except KeyError:
        continue

df_bar = pd.DataFrame([XISTpos_NSGCT_bar, XISTneg_NSGCT_bar, XISTpos_SGCT_bar, XISTpos_PanCan_bar, XISTneg_PanCan_bar]).T
df_bar.columns = ['XISTpos_NSGCT', 'XISTneg_NSGCT', 'XISTpos_SGCT', 'XISTpos_PanCan', 'XISTneg_PanCan']

df_bar.to_csv("/Users/ananthansadagopan/Downloads/temp_output.csv", index=False)


"""
df_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

lineage_dict = dict(zip(df_ref['Barcode'].tolist(), df_ref['Classification'].tolist()))

XISTpos_PanCan_lineage = []
XISTneg_PanCan_lineage = []

for a in XISTpos_PanCan_bar:
    XISTpos_PanCan_lineage.append(lineage_dict[a])

for a in XISTneg_PanCan_bar:
    XISTneg_PanCan_lineage.append(lineage_dict[a])

df_temp = pd.DataFrame([XISTpos_PanCan_lineage, XISTneg_PanCan_lineage]).T
df_temp.columns = ['XIST+', 'XIST-']

df_temp.to_csv("/Users/ananthansadagopan/Downloads/temp_lineage_output.csv", index=False)
"""

q=1
cn_1 = [XISTpos_SGCT_cn.count(q)/len(XISTpos_SGCT_cn), XISTneg_NSGCT_cn.count(q)/len(XISTneg_NSGCT_cn), XISTpos_NSGCT_cn.count(q)/len(XISTpos_NSGCT_cn), XISTneg_PanCan_cn.count(q)/len(XISTneg_PanCan_cn), XISTpos_PanCan_cn.count(q)/len(XISTpos_PanCan_cn)]
q=2
cn_2 = [XISTpos_SGCT_cn.count(q)/len(XISTpos_SGCT_cn), XISTneg_NSGCT_cn.count(q)/len(XISTneg_NSGCT_cn), XISTpos_NSGCT_cn.count(q)/len(XISTpos_NSGCT_cn), XISTneg_PanCan_cn.count(q)/len(XISTneg_PanCan_cn), XISTpos_PanCan_cn.count(q)/len(XISTpos_PanCan_cn)]
q=3
cn_3 = [XISTpos_SGCT_cn.count(q)/len(XISTpos_SGCT_cn), XISTneg_NSGCT_cn.count(q)/len(XISTneg_NSGCT_cn), XISTpos_NSGCT_cn.count(q)/len(XISTpos_NSGCT_cn), XISTneg_PanCan_cn.count(q)/len(XISTneg_PanCan_cn), XISTpos_PanCan_cn.count(q)/len(XISTpos_PanCan_cn)]
q=4
cn_4 = [XISTpos_SGCT_cn.count(q)/len(XISTpos_SGCT_cn), XISTneg_NSGCT_cn.count(q)/len(XISTneg_NSGCT_cn), XISTpos_NSGCT_cn.count(q)/len(XISTpos_NSGCT_cn), XISTneg_PanCan_cn.count(q)/len(XISTneg_PanCan_cn), XISTpos_PanCan_cn.count(q)/len(XISTpos_PanCan_cn)]


"""
print(cn_1)
print(cn_2)
print(cn_3)
print(cn_4)
"""


"""

XISTpos_NSGCT_g2 = XISTpos_NSGCT_cn.count(3) + XISTpos_NSGCT_cn.count(4)
XISTpos_NSGCT_leq2 = XISTpos_NSGCT_cn.count(1) + XISTpos_NSGCT_cn.count(2)

XISTneg_NSGCT_g2 = XISTneg_NSGCT_cn.count(3) + XISTneg_NSGCT_cn.count(4)
XISTneg_NSGCT_leq2 = XISTneg_NSGCT_cn.count(1) + XISTneg_NSGCT_cn.count(2)

fisher_test_list_NSGCT = [[XISTpos_NSGCT_g2, XISTpos_NSGCT_leq2], [XISTneg_NSGCT_g2, XISTneg_NSGCT_leq2]] 

o, p = scipy.stats.fisher_exact(fisher_test_list_NSGCT, alternative='two-sided')

print(o)
print(p)

XISTpos_PanCan_g1 = XISTpos_PanCan_cn.count(2) + XISTpos_PanCan_cn.count(3) + XISTpos_PanCan_cn.count(4)
XISTpos_PanCan_leq1 = XISTpos_PanCan_cn.count(1)
XISTneg_PanCan_g1 = XISTneg_PanCan_cn.count(2) + XISTneg_PanCan_cn.count(3) + XISTneg_PanCan_cn.count(4)
XISTneg_PanCan_leq1 = XISTneg_PanCan_cn.count(1)

fisher_test_list_PanCan = [[XISTpos_PanCan_g1, XISTpos_PanCan_leq1], [XISTneg_PanCan_g1, XISTneg_PanCan_leq1]] 

o, p = scipy.stats.fisher_exact(fisher_test_list_PanCan, alternative='two-sided')

print(o)
print(p)

"""


#XISTpos

print("CHI SQUARED NSGCT")

XISTpos_NSGCT_obs = [cn_1[2]*len(XISTpos_NSGCT_cn), cn_2[2]*len(XISTpos_NSGCT_cn), cn_3[2]*len(XISTpos_NSGCT_cn), cn_4[2]*len(XISTpos_NSGCT_cn)]
XISTpos_NSGCT_exp = [cn_1[1]*len(XISTpos_NSGCT_cn), cn_2[1]*len(XISTpos_NSGCT_cn), cn_3[1]*len(XISTpos_NSGCT_cn), cn_4[1]*len(XISTpos_NSGCT_cn)]

print(XISTpos_NSGCT_obs)
print(XISTpos_NSGCT_exp)

c, p = scipy.stats.chisquare(XISTpos_NSGCT_obs, f_exp=XISTpos_NSGCT_exp, ddof=0, axis=0)

print(c)
print(p)



print("CHI SQUARED PanCan")

XISTpos_PanCan_obs = [cn_1[4]*len(XISTpos_PanCan_cn), cn_2[4]*len(XISTpos_PanCan_cn), cn_3[4]*len(XISTpos_PanCan_cn), cn_4[4]*len(XISTpos_PanCan_cn)]
XISTpos_PanCan_exp = [cn_1[3]*len(XISTpos_PanCan_cn), cn_2[3]*len(XISTpos_PanCan_cn), cn_3[3]*len(XISTpos_PanCan_cn), cn_4[3]*len(XISTpos_PanCan_cn)]

print(XISTpos_PanCan_obs)
print(XISTpos_PanCan_exp)

c, p = scipy.stats.chisquare(XISTpos_PanCan_obs, f_exp=XISTpos_PanCan_exp, ddof=0, axis=0)

print(c)
print(p)








fig, ax = plt.subplots()

ind = np.arange(len(cn_1))
width = 0.4

cn_1 = np.array(cn_1)*100
cn_2 = np.array(cn_2)*100
cn_3 = np.array(cn_3)*100
cn_4 = np.array(cn_4)*100

p1 = ax.bar(ind, cn_1, width, color='b')
p2 = ax.bar(ind, cn_2, width, bottom=cn_1, color='orange')
p3 = ax.bar(ind, cn_3, width, bottom=cn_1+cn_2, color='g')
p4 = ax.bar(ind, cn_4, width, bottom=cn_1+cn_2+cn_3,color='purple')

#labels = ['XIST+\nSGCT\n(n=' + str(len(XISTpos_SGCT_cn)) + ")", 'XIST-\nNSGCT\n(n=' + str(len(XISTneg_NSGCT_cn)) + ")", 'XIST+\nNSGCT\n(n=' + str(len(XISTpos_NSGCT_cn)) + ")", 'XIST-\nPanCan\n(n=' + str(len(XISTneg_PanCan_cn)) + ")", 'XIST+\nPanCan\n(n=' + str(len(XISTpos_PanCan_cn)) + ")"]

labels = ['XIST+ (' + str(len(XISTpos_SGCT_cn)) + ")", 'XIST- (' + str(len(XISTneg_NSGCT_cn)) + ")", 'XIST+ (' + str(len(XISTpos_NSGCT_cn)) + ")", 'XIST- (' + str(len(XISTneg_PanCan_cn)) + ")", 'XIST+ (' + str(len(XISTpos_PanCan_cn)) + ")"]

one_patch = mpatches.Patch(color='b', label='1')
two_patch = mpatches.Patch(color='orange', label='2')
three_patch = mpatches.Patch(color='g', label='3')
four_patch = mpatches.Patch(color='purple', label='4')

leg = ax.legend(handles=[one_patch, two_patch, three_patch, four_patch], facecolor='white', loc="upper left", title="chrX CN", ncol=1, bbox_to_anchor=(1,1))

frame = leg.get_frame()
frame.set_edgecolor("black")
frame.set_linewidth(0.5)

ax.set_xticks(ind)
ax.set_xticklabels(labels, fontsize=14)


ax.set_ylabel("Percent of Samples", size=14)
ax.set_yticklabels([0, 20, 40, 60, 80, 100], size=14)
#ax.set_xlabel("chrX Copy Number", size=14)

plt.ylim([0, 100.2])

ax.grid(False)
ax.set_facecolor("white")
ax.spines['bottom'].set_color('0')
ax.spines['left'].set_color('0')

#ax.spines['bottom'].set_lw(0.9)
#ax.spines['left'].set_lw(0.9)

ax.tick_params(bottom='on', left='on')

fig.tight_layout()

plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/chrX_copy_number_stacked_barplots.png", dpi=dpi_set, bbox_inches = 'tight')
plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/chrX_copy_number_stacked_barplots.pdf", dpi=dpi_set, bbox_inches = 'tight')



