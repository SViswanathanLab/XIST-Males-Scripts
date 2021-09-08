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
import matplotlib
import matplotlib.patches as mpatches
from matplotlib.ticker import StrMethodFormatter, NullFormatter, ScalarFormatter, FormatStrFormatter
import matplotlib

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.linewidth'] = 0.3
sns.set(font_scale=0.6)

figsize = (3.5,8.3*9/7) # the size of the figure - changes the shape of the squares in the comut
dpi_set = 72 # change the output resolution

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

orig_cols = df.columns.tolist()

df_barcode = df['Barcode'].tolist()

temp_vals = df['XIST_TPM'].tolist()

log2_temp_vals = []

for a in temp_vals:
    log2_temp_vals.append(math.log(a+1, 2))

TPM_dict = dict(zip(df_barcode, log2_temp_vals))

end_barcode = []

for a in df_barcode:
    end_barcode.append(int(split_advanced(a, "-", 3)[1]))

df['End_Val'] = end_barcode

df = df[df['End_Val']==11] #only 11s are present for normal samples
df = df[df['Classification'] != "UNKNOWN"]
df = df[df['Gender'] != "UNKNOWN"]

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df = df[~(df['Barcode'].str.contains(a))]

df2 = df[df['Gender']=="FEMALE"]
df = df[df['Gender']=="MALE"]

df_class = df['Classification'].tolist()
df2_class = df['Classification'].tolist()

male_exp_unp = df['XIST_TPM'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

female_exp_unp = df2['XIST_TPM'].tolist()

female_exp = []

for a in female_exp_unp:
    if a>=0:
        female_exp.append(math.log(a+1, 2))
    else:
        female_exp.append(0)

n_bins = 55

#fix bc these are females

fig, axs = plt.subplots(9, 1, sharex=True, tight_layout=True, figsize=figsize)

sns.distplot(female_exp, hist=False, kde=True, ax=axs[0],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[0].axvline(TPM_dict['TCGA-E9-A1NA-11'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-E9-A1NA-11']-1))

#axs[0].text(0.75, 0.98, TPM, transform=axs[0].transAxes,
#     fontsize=6, verticalalignment='top')

sns.distplot(male_exp, hist=False, kde=True, ax=axs[1],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[1].axvline(TPM_dict['TCGA-22-5482-11'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-22-5482-11']-1))

#axs[1].text(0.75, 0.98, TPM, transform=axs[1].transAxes,
#     fontsize=6, verticalalignment='top')





df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

orig_cols = df.columns.tolist()

df_barcode = df['Barcode'].tolist()

temp_vals = df['XIST_TPM'].tolist()

log2_temp_vals = []

for a in temp_vals:
    log2_temp_vals.append(math.log(a+1, 2))

TPM_dict = dict(zip(df_barcode, log2_temp_vals))

end_barcode = []

for a in df_barcode:
    end_barcode.append(int(split_advanced(a, "-", 3)[1]))

df['End_Val'] = end_barcode

df = df[df['End_Val']<10]
df = df[df['Classification'] != "UNKNOWN"]
df = df[df['Gender'] != "UNKNOWN"]

curr_barcode = df['Barcode'].tolist()

trunc_barcode = []

for a in curr_barcode:
    trunc_barcode.append(split_advanced(a, "-", 3)[0])

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
    
    new_rows.append([str(a)+"-09", temp_log2, temp_average, temp_classification, temp_gender, temp_classification, "UNKNOWN"]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)

df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Barcode'].str.contains(a))]

df = pd.concat([df, df_from_new_rows], ignore_index=True)

invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM', 'TCGA-AB-2872', 'TCGA-B0-4696', 'TCGA-B0-4846', 'TCGA-CJ-4642', 'TCGA-CV-7428', 'TCGA-CZ-4862']

for a in invalid_ids:
    df = df[~(df['Barcode'].str.contains(a))]

df2 = df[df['Gender']=="FEMALE"]
df = df[df['Gender']=="MALE"]

OV_df = df2[df2['Classification']=="OV"]

OV_exp = []
OV_exp_unp = OV_df['XIST_TPM'].tolist()

for a in OV_exp_unp:
    if a>=0:
        OV_exp.append(math.log(a+1, 2))
    else:
        OV_exp.append(0)



NSGCT_df = df[df['Secondary_Class']=="TGCT-NS"]
SGCT_df = df[df['Secondary_Class']=="TGCT-S"]
LUSC_df = df[df['Classification']=="LUSC"]

NSGCT_exp_unp = NSGCT_df['XIST_TPM'].tolist()

NSGCT_exp = []

for a in NSGCT_exp_unp:
    if a>=0:
        NSGCT_exp.append(math.log(a+1, 2))
    else:
        NSGCT_exp.append(0)
   
SGCT_exp_unp = SGCT_df['XIST_TPM'].tolist()

SGCT_exp = []

for a in SGCT_exp_unp:
    if a>=0:
        SGCT_exp.append(math.log(a+1, 2))
    else:
        SGCT_exp.append(0) 

LUSC_exp_unp = LUSC_df['XIST_TPM'].tolist()

LUSC_exp = []

for a in LUSC_exp_unp:
    if a>=0:
        LUSC_exp.append(math.log(a+1, 2))
    else:
        LUSC_exp.append(0)
        
df_class = df['Classification'].tolist()
df2_class = df['Classification'].tolist()

male_exp_unp = df['XIST_TPM'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

female_exp_unp = df2['XIST_TPM'].tolist()

female_exp = []

for a in female_exp_unp:
    if a>=0:
        female_exp.append(math.log(a+1, 2))
    else:
        female_exp.append(0)

sns.distplot(SGCT_exp, hist=False, kde=True, ax=axs[2],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[2].axvline(TPM_dict['TCGA-XY-A8S3-01'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-XY-A8S3-01']-1))

print(TPM)

#axs[2].text(0.75, 0.98, TPM, transform=axs[2].transAxes,
#     fontsize=6, verticalalignment='top')

sns.distplot(NSGCT_exp, hist=False, kde=True, ax=axs[4],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[4].axvline(TPM_dict['TCGA-2G-AAFZ-01'], c="green")

axs[4].set_yticks([0, 0.1, 0.2])
axs[4].set_ylim([0, 0.22])

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-2G-AAFZ-01']-1))

#axs[4].text(0.75, 0.98, TPM, transform=axs[4].transAxes,
#     fontsize=6, verticalalignment='top')

"""
sns.distplot(NSGCT_exp, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[5].axvline(TPM_dict['TCGA-SN-A6IS-01'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-SN-A6IS-01']-1))

#axs[5].text(0.75, 0.98, TPM, transform=axs[5].transAxes,
#     fontsize=6, verticalalignment='top')
"""

sns.distplot(NSGCT_exp, hist=False, kde=True, ax=axs[3],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[3].axvline(TPM_dict['TCGA-2G-AAL5-01'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-2G-AAL5-01']-1))

axs[3].set_yticks([0, 0.1, 0.2])
axs[3].set_ylim([0, 0.22])

#axs[3].text(0.75, 0.98, TPM, transform=axs[3].transAxes,
#     fontsize=6, verticalalignment='top')

sns.distplot(LUSC_exp, hist=False, kde=True, ax=axs[6],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[6].axvline(TPM_dict['TCGA-56-7822-01'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-56-7822-01']-1))





print("hi")

sns.distplot(OV_exp, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[8].axvline(TPM_dict['TCGA-23-2078-01'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-23-2078-01']-1))

print(TPM)

sns.distplot(OV_exp, hist=False, kde=True, ax=axs[7],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[7].axvline(TPM_dict['TCGA-29-1702-01'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-29-1702-01']-1))
print(TPM)


#axs[7].text(0.75, 0.98, TPM, transform=axs[7].transAxes,
#     fontsize=6, verticalalignment='top')


"""
sns.distplot(LUSC_exp, hist=False, kde=True, ax=axs[8],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})


axs[8].axvline(TPM_dict['TCGA-37-4130-01'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-37-4130-01']-1))

#axs[8].text(0.75, 0.98, TPM, transform=axs[8].transAxes,
#     fontsize=6, verticalalignment='top')
"""

sns.distplot(LUSC_exp, hist=False, kde=True, ax=axs[5],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

axs[5].axvline(TPM_dict['TCGA-77-6843-01'], c="green")

TPM = "XIST TPM = " + str("{:.1f}".format(2**TPM_dict['TCGA-77-6843-01']-1))

#axs[6].text(0.75, 0.98, TPM, transform=axs[6].transAxes,
#     fontsize=6, verticalalignment='top')

"""
axs[0].set_ylabel("Density (XIST+ F-LUAD-Normal)")
axs[1].set_ylabel("Density (XIST- M-LUSC-Normal)")
axs[2].set_ylabel("Density (XIST+ M-SGCT)")
axs[3].set_ylabel("Density (XIST+ M-NSGCT)")
axs[4].set_ylabel("Density (XIST+ M-NSGCT)")
axs[5].set_ylabel("Density (XIST- M-NSGCT)")
axs[6].set_ylabel("Density (XIST+ M-LUSC)")
axs[7].set_ylabel("Density (XIST+ M-LUSC)")
axs[8].set_ylabel("Density (XIST- M-LUSC)")
"""

axs[0].set_ylabel("Density")
axs[1].set_ylabel("Density")
axs[2].set_ylabel("Density")
axs[3].set_ylabel("Density")
axs[4].set_ylabel("Density")
axs[5].set_ylabel("Density")
axs[6].set_ylabel("Density")
axs[7].set_ylabel("Density")
axs[8].set_ylabel("Density")

axs[6].set_xlabel("XIST Expression (log$_2$(TPM+1))")

a=0
while a<9:
    axs[a].grid(False)
    axs[a].set_facecolor("white")
    axs[a].spines['bottom'].set_color('0')
    axs[a].spines['left'].set_color('0')
    axs[a].spines['left'].set_lw(0.5)
    axs[a].spines['bottom'].set_lw(0.5)
    axs[a].tick_params(bottom='on', left='on')
    a=a+1

plt.xlim(0, 12.5)


sample_patch = mpatches.Patch(color='green', label='XIST\nExpression\nin Sample')


a=0
while a<9:
    axs[a].yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
    axs[a].xaxis.set_tick_params(width=0.4, size=3, pad=2)
    axs[a].yaxis.set_tick_params(width=0.4, size=3, pad=2)
    a=a+1


#leg = axs[0].legend(handles=[sample_patch], facecolor='white', bbox_to_anchor=(1, 0.89))

#frame = leg.get_frame()
#frame.set_edgecolor("black")
#frame.set_linewidth(0.5)


print(2**statistics.median(SGCT_exp)-1)
print(2**statistics.median(NSGCT_exp)-1)
print(2**statistics.median(LUSC_exp)-1)


fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/select_sample_annotated_XIST_exp_Xena_TCGA_male_and_female_histogram.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/select_sample_annotated_XIST_exp_Xena_TCGA_male_and_female_histogram.pdf", dpi=dpi_set, bbox_inches = 'tight')




