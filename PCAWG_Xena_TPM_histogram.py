#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 19:06:35 2021

@author: ananthansadagopan
"""

import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
import statistics
import seaborn as sns
import collections
from matplotlib.ticker import StrMethodFormatter, NullFormatter, ScalarFormatter, FormatStrFormatter
import matplotlib

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.linewidth'] = 0.3

figsize = (4.11,3.52) # the size of the figure - changes the shape of the squares in the comut
dpi_set = 72 # change the output resolution
sns.set(font_scale=0.6)


df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/XIST_expression_summary.xlsx")

df_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/project_code_sp", sep="\t")

ids_ref = df_ref['icgc_specimen_id'].tolist()
proj_code = df_ref['dcc_project_code'].tolist()

id_dict = dict(zip(ids_ref, proj_code))

df_barcode = df['Barcode'].tolist()

ids_to_remove = []

for a in df_barcode:
    if id_dict[a][-2:] == "US":
        ids_to_remove.append(a)
        
df = df[~(df['Barcode'].isin(ids_to_remove))]

orig_cols = df.columns.tolist()

df_barcode = df['Barcode'].tolist()

temp_vals = df['Specimen_Type'].tolist()

end_barcode = []

for a in temp_vals:
    if "Normal" in a or "normal" in a:
        end_barcode.append(11)
    else:
        end_barcode.append(1)        

df['End_Val'] = end_barcode

df = df[df['End_Val']<10]
df = df[df['Classification'] != "UNKNOWN"]
df = df[df['Gender'] != "UNKNOWN"]

dup_items = list(set([item for item, count in collections.Counter(df['Patient'].tolist()).items() if count > 1]))

new_rows = []

for a in dup_items:
    temp_df = df[df['Patient']==a]
    vals = temp_df['XIST_FPKM_UQ'].tolist()
    
    if len(vals)<2:
        print("ERROR")
    
    temp_average = sum(vals)/len(vals)
    temp_log2 = math.log(temp_average+0.001, 2)
    
    temp_gender = temp_df['Gender'].tolist()[0]
    temp_classification = temp_df['Classification'].tolist()[0]
    
    new_rows.append([str(a)+"-09", temp_log2, temp_average, a, "unknown", "unknown", "unknown", "unknown", temp_gender, temp_classification, "Averaged"]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)

df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Patient'].str.contains(a))]

df = pd.concat([df, df_from_new_rows], ignore_index=True)




df2 = df[df['Gender']=="female"]
df = df[df['Gender']=="male"]

df_class = df['Classification'].tolist()

df2_class = df['Classification'].tolist()

male_exp_unp = df['XIST_FPKM_UQ'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

female_exp_unp = df2['XIST_FPKM_UQ'].tolist()

female_exp = []

for a in female_exp_unp:
    if a>=0:
        female_exp.append(math.log(a+1, 2))
    else:
        female_exp.append(0)

n_bins = 55

fig, axs = plt.subplots(4, 1, sharex=True, tight_layout=True, figsize=figsize)

bins_list = np.arange(0,8,0.1)

#nice_palette = ['#9D2727', '#CFB997', '#E79898', '#9CDBF2', '#C445A1', '#DDDD00', '#36E2EE', '#008BB5', '#DACAFF', '#FFCAFF']





axs[3].hist(male_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")
#axs[3].hist([PanCan_exp, TGCT_exp], density=False, bins=bins_list, color=["#9D2727","#9CDBF2"], linewidth=0.5, edgecolor="none", stacked=True)

axs[1].hist(female_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")

print("COUNTS")

q=0
counter=0

while q<len(male_exp):
    if male_exp[q] >= 0.5:
        counter=counter+1
    q=q+1
        
print(counter)
print(len(male_exp))


q=0
counter=0

while q<len(female_exp):
    if female_exp[q] >= 0.5:
        counter=counter+1
    q=q+1
        
print(counter)
print(len(female_exp))



axs[3].axvline(0.5, c="black", ls="--", lw=0.2)
axs[1].axvline(0.5, c="black", ls="--", lw=0.2)

axs[3].grid(False)
axs[3].set_facecolor("white")
axs[3].spines['bottom'].set_color('0')
axs[3].spines['left'].set_color('0')

axs[3].tick_params(bottom='on', left='on', which="both")


axs[1].grid(False)
axs[1].set_facecolor("white")
axs[1].spines['bottom'].set_color('0')
axs[1].spines['left'].set_color('0')

axs[1].tick_params(bottom='on', left='on', which="both")

























df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/XIST_expression_summary.xlsx")

df_barcode = df['Barcode'].tolist()

ids_to_remove = []

for a in df_barcode:
    if id_dict[a][-2:] == "US":
        ids_to_remove.append(a)
        
df = df[~(df['Barcode'].isin(ids_to_remove))]

orig_cols = df.columns.tolist()

df_barcode = df['Barcode'].tolist()


temp_vals = df['Specimen_Type'].tolist()

end_barcode = []

for a in temp_vals:
    if "Normal" in a or "normal" in a:
        end_barcode.append(11)
    else:
        end_barcode.append(1)        

df['End_Val'] = end_barcode

df = df[df['End_Val']==11]

df = df[df['Classification'] != "UNKNOWN"]
df = df[df['Gender'] != "UNKNOWN"]


"""
dup_items = list(set([item for item, count in collections.Counter(df['Patient'].tolist()).items() if count > 1]))

new_rows = []

for a in dup_items:
    temp_df = df[df['Patient']==a]
    vals = temp_df['XIST_FPKM_UQ'].tolist()
    
    if len(vals)<2:
        print("ERROR")
    
    temp_average = sum(vals)/len(vals)
    temp_log2 = math.log(temp_average+0.001, 2)
    
    temp_gender = temp_df['Gender'].tolist()[0]
    temp_classification = temp_df['Classification'].tolist()[0]
    
    new_rows.append([str(a)+"-09", temp_log2, temp_average, a, temp_gender, temp_classification, "Averaged"]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)

df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Patient'].str.contains(a))]

df = pd.concat([df, df_from_new_rows], ignore_index=True)
"""



df2 = df[df['Gender']=="female"]
df = df[df['Gender']=="male"]

df_class = df['Classification'].tolist()

df2_class = df['Classification'].tolist()

male_exp_unp = df['XIST_FPKM_UQ'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

female_exp_unp = df2['XIST_FPKM_UQ'].tolist()

female_exp = []

for a in female_exp_unp:
    if a>=0:
        female_exp.append(math.log(a+1, 2))
    else:
        female_exp.append(0)

axs[2].hist(male_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")

axs[0].hist(female_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")

"""

sns.distplot(male_exp, hist=False, kde=True, ax=axs[2], rug=False,
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

sns.distplot(female_exp, hist=False, kde=True, ax=axs[0], rug=False,
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

"""


"""
axs[2].text(0.88, 0.98, " TGCT:\n   N=" + str(tgct_yes), transform=axs[2].transAxes,
     fontsize=10, verticalalignment='top', color="blue")

axs[2].text(0.89, 0.8, "Other:\n  N=" + str(other_yes_male), transform=axs[2].transAxes,
     fontsize=10, verticalalignment='top', color="red")

axs[2].text(0.025, 0.98, "TGCT:\n  N=" + str(tgct_no), transform=axs[2].transAxes,
     fontsize=10, verticalalignment='top', color="blue")

axs[2].text(-0.005, 0.8, "    Other:\n   N=" + str(other_no_male), transform=axs[2].transAxes,
     fontsize=10, verticalalignment='top', color="red")
"""

"""
axs[0].text(0.87, 0.98, "N=" + str(other_yes_female), transform=axs[0].transAxes,
     fontsize=10, verticalalignment='top', color="red")

axs[0].text(0.043, 0.98, "N=" + str(other_no_female), transform=axs[0].transAxes,
     fontsize=10, verticalalignment='top', color="red")
"""

# We can set the number of bins with the `bins` kwarg
#axs[0].hist(male_exp, bins=n_bins)
#axs[1].hist(female_exp, bins=n_bins)

plt.xlim([0, 8.1])

print("COUNTS")

q=0
counter=0

while q<len(male_exp):
    if male_exp[q] >= 0.5:
        counter=counter+1
    q=q+1
        
print(counter)
print(len(male_exp))


q=0
counter=0

while q<len(female_exp):
    if female_exp[q] >= 0.5:
        counter=counter+1
    q=q+1
        
print(counter)
print(len(female_exp))


axs[2].axvline(0.5, c="black", ls="--", lw=0.2)
axs[0].axvline(0.5, c="black", ls="--", lw=0.2)

axs[2].grid(False)
axs[2].set_facecolor("white")
axs[2].spines['bottom'].set_color('0')
axs[2].spines['left'].set_color('0')

axs[2].tick_params(bottom='on', left='on', which="both")

axs[0].grid(False)
axs[0].set_facecolor("white")
axs[0].spines['bottom'].set_color('0')
axs[0].spines['left'].set_color('0')

axs[0].tick_params(bottom='on', left='on', which="both")

axs[0].set_yscale('log')
axs[2].set_yscale('log')
axs[1].set_yscale('log')
axs[3].set_yscale('log')


axs[3].set_ylabel("Number of\nSamples")
#axs[1].set_xlabel("XIST Expression (log$_2$(TPM+1))")
axs[1].set_ylabel("Number of\nSamples")
axs[2].set_ylabel("Number of\nSamples")
axs[3].set_xlabel("XIST Expression (log$_2$(FPKM-UQ+1))")
#axs[2].set_xlabel("XIST Expression (log$_2$(TPM+1))")
axs[0].set_ylabel("Number of\nSamples")


axs[0].spines['bottom'].set_lw(0.2)
axs[0].spines['left'].set_lw(0.2)
axs[1].spines['bottom'].set_lw(0.2)
axs[1].spines['left'].set_lw(0.2)
axs[2].spines['bottom'].set_lw(0.2)
axs[2].spines['left'].set_lw(0.2)
axs[3].spines['bottom'].set_lw(0.2)
axs[3].spines['left'].set_lw(0.2)


axs[0].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[0].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[0].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")
axs[0].set_ylim([0, 10.5])

axs[0].set_yticks([2, 3, 4, 5, 6, 7, 8, 9], minor=True)

axs[2].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[2].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[2].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")
axs[2].set_ylim([0, 78])

axs[2].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70], minor=True)

axs[1].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[1].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[1].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")
axs[1].set_ylim([0, 10.5])

axs[1].set_yticks([2, 3, 4, 5, 6, 7, 8, 9], minor=True)

axs[3].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[3].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[3].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")
axs[3].set_ylim([0, 350])

axs[3].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300], minor=True)

axs[0].xaxis.labelpad = 1
axs[0].yaxis.labelpad = 8.5
axs[1].xaxis.labelpad = 1
axs[1].yaxis.labelpad = 8.5
axs[2].xaxis.labelpad = 1
axs[2].yaxis.labelpad = 8.5
axs[3].xaxis.labelpad = 1
axs[3].yaxis.labelpad = 8.5

axs[0].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))

axs[1].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))

axs[2].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))

axs[3].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))



"""
axs[3].set_ylabel("Density (Male Cancers)")
#axs[1].set_xlabel("XIST Expression (log$_2$(TPM+1))")
axs[1].set_ylabel("Density (Female Cancers)")
axs[2].set_ylabel("Density (Normal Males)")
axs[3].set_xlabel("XIST Expression (log$_2$(TPM+1))")
axs[0].set_ylabel("Density (Normal Females)")
"""

#fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PCAWG_OP5_hist_log_normal_and_cancer_XIST_exp_Xena_TCGA_male_and_female.png", dpi=dpi_set, bbox_inches = 'tight')
#fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/PCAWG_OP5_hist_log_normal_and_cancer_XIST_exp_Xena_TCGA_male_and_female.pdf", dpi=dpi_set, bbox_inches = 'tight')




