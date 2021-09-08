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

figsize = (3.77,1.62) # the size of the figure - changes the shape of the squares in the comut
dpi_set = 72 # change the output resolution
sns.set(font_scale=0.4)


df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/GTEx/male_and_female_XIST_expression_GTEx_summary.csv")

orig_cols = df.columns.tolist()

df = df[df['Classification']!="UNKNOWN"]
df = df[df['Gender']!="UNKNOWN"]

curr_barcode = df['Barcode'].tolist()

trunc_barcode = []

for a in curr_barcode:
    trunc_barcode.append(split_advanced(a, "-", 3)[0])

dup_items = list(set([item for item, count in collections.Counter(trunc_barcode).items() if count > 1]))

new_rows = []

for a in dup_items:
    temp_df = df[df['Barcode'].str.contains(a)]
    vals = temp_df['XIST_Expression'].tolist()
    
    if len(vals)<2:
        print("ERROR")
    
    temp_average = sum(vals)/len(vals)
    
    temp_gender = temp_df['Gender'].tolist()[0]
    temp_classification = temp_df['Classification'].tolist()[0]
    
    new_rows.append([str(a)+"-AVERAGE", temp_average, temp_classification, temp_gender]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)
df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Barcode'].str.contains(a))]

df = pd.concat([df, df_from_new_rows], ignore_index=True)

df2 = df[df['Gender']=="FEMALE"]
df = df[df['Gender']=="MALE"]

df_class = df['Classification'].tolist()
df2_class = df['Classification'].tolist()

male_exp_unp = df['XIST_Expression'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

female_exp_unp = df2['XIST_Expression'].tolist()

female_exp = []

for a in female_exp_unp:
    if a>=0:
        female_exp.append(math.log(a+1, 2))
    else:
        female_exp.append(0)

n_bins = 40

fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=figsize)

#axs[0].hist(male_exp, density=False, bins=45, color="red", linewidth=0.5, edgecolor="black")

#axs[0].set_yscale('log')

#axs[1].hist(female_exp, density=False, bins=90, color="red", linewidth=0.5, edgecolor="black")

#axs[1].set_yscale('log')

bins_list = np.arange(0,10.01,0.1)

axs[0].hist(male_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")
axs[1].hist(female_exp, density=False, bins=bins_list, color="#9D2727", linewidth=0.5, edgecolor="none")

print(len(male_exp))
print(len(female_exp))

axs[0].set_ylabel("Number of\nSamples")
axs[1].set_ylabel("Number of\nSamples")


"""
sns.distplot(male_exp, hist=False, kde=True, ax=axs[0],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})

sns.distplot(female_exp, hist=False, kde=True, ax=axs[1],
             bins=n_bins, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 0.5})
"""

tgct_yes = 0
tgct_no = 0

other_yes_male = 0
other_no_male = 0

a=0
while a<len(male_exp):
    if df_class[a] == "TGCT":
        if male_exp[a] >= 4:
            tgct_yes = tgct_yes + 1
        else:
            tgct_no = tgct_no + 1
    else:
        if male_exp[a] >= 4:
            other_yes_male = other_yes_male + 1
        else:
            other_no_male = other_no_male + 1
    a=a+1

"""
axs[0].text(0.941, 0.98, "N=" + str(other_yes_male), transform=axs[0].transAxes,
     fontsize=10, verticalalignment='top', color="red")

axs[0].text(0.026, 0.98, "N=" + str(other_no_male), transform=axs[0].transAxes,
     fontsize=10, verticalalignment='top', color="red")
"""

other_yes_female = 0
other_no_female = 0

a=0
while a<len(female_exp):
    if female_exp[a] >= 2:
        other_yes_female = other_yes_female + 1
    else:
        other_no_female = other_no_female + 1
    a=a+1

"""
axs[1].text(0.88, 0.98, "N=" + str(other_yes_female), transform=axs[1].transAxes,
     fontsize=10, verticalalignment='top', color="red")

axs[1].text(0.026, 0.98, "N=" + str(other_no_female), transform=axs[1].transAxes,
     fontsize=10, verticalalignment='top', color="red")
"""

# We can set the number of bins with the `bins` kwarg
#axs[0].hist(male_exp, bins=n_bins)
#axs[1].hist(female_exp, bins=n_bins)

axs[1].set_xlabel("XIST Expression (log$_2$(TPM+1))")

print(str(max(female_exp)))
print(str(max(male_exp)))

plt.xlim(0, 10)

axs[0].grid(False)
axs[0].set_facecolor("white")
axs[0].spines['bottom'].set_color('0')
axs[0].spines['left'].set_color('0')

axs[0].tick_params(bottom='on', left='on', which="both")

axs[1].grid(False)
axs[1].set_facecolor("white")
axs[1].spines['bottom'].set_color('0')
axs[1].spines['left'].set_color('0')

axs[1].tick_params(bottom='on', left='on', which="both")

#axs[0].axvline(4, c="black")
#axs[1].axvline(4, c="black")

axs[0].set_yscale('log')
axs[1].set_yscale('log')


axs[0].spines['bottom'].set_lw(0.2)
axs[0].spines['left'].set_lw(0.2)
axs[1].spines['bottom'].set_lw(0.2)
axs[1].spines['left'].set_lw(0.2)

axs[0].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[0].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[0].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")

axs[1].xaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[1].yaxis.set_tick_params(width=0.2, size=2, pad=2)
axs[1].yaxis.set_tick_params(width=0.2, size=1, pad=2, which="minor")

axs[0].xaxis.labelpad = 1
axs[0].yaxis.labelpad = 1
axs[1].xaxis.labelpad = 1
axs[1].yaxis.labelpad = 1

axs[0].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))

axs[1].yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))

axs[0].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900, 2000, 3000, 4000, 5000, 6000, 7000], minor=True)

axs[0].set_yticks([1, 10, 100, 1000])

axs[1].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300], minor=True)

axs[0].set_ylim([0, 8550])
axs[1].set_ylim([0, 250])


fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XIST_exp_website_GTEx_male_and_female_histogram.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XIST_exp_website_GTEx_male_and_female_histogram.pdf", dpi=dpi_set, bbox_inches = 'tight')




