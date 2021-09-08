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
from matplotlib.ticker import StrMethodFormatter, NullFormatter, ScalarFormatter, FormatStrFormatter
import matplotlib

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.linewidth'] = 0.3
sns.set(font_scale=0.4)


"""

DO NOT USE
"""








figsize = (3.77,1.62) # the size of the figure - changes the shape of the squares in the comut
dpi_set = 72 # change the output resolution

df_annot = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/CCLE_XIST_summary.xlsx")

males_df = df_annot[df_annot['Gender']=="MALE"]

male_exp_unp = males_df['XIST Expression (TPM)'].tolist()

male_exp = []

for a in male_exp_unp:
    if a>=0:
        male_exp.append(math.log(a+1, 2))
    else:
        male_exp.append(0)

females_df = df_annot[df_annot['Gender']=="FEMALE"]

female_exp_unp = females_df['XIST Expression (TPM)'].tolist()

female_exp = []

for a in female_exp_unp:
    if a>=0:
        female_exp.append(math.log(a+1, 2))
    else:
        female_exp.append(0)

n_bins = 25

fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=figsize)

axs[0].hist(male_exp, density=False, bins=90, color="#9D2727", linewidth=0.5, edgecolor="none")
axs[1].hist(female_exp, density=False, bins=90, color="#9D2727", linewidth=0.5, edgecolor="none")

print(len(male_exp))
print(len(female_exp))

axs[0].set_ylabel("Number of\nCell Lines")
axs[1].set_ylabel("Number of\nCell Lines")

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

# We can set the number of bins with the `bins` kwarg
#axs[0].hist(male_exp, bins=n_bins)
#axs[1].hist(female_exp, bins=n_bins)

axs[1].set_xlabel("XIST Expression (log$_2$(TPM+1))")

print(str(max(female_exp)))
print(str(max(male_exp)))

plt.xlim(0, 10.1)

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

axs[0].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200], minor=True)
axs[1].set_yticks([2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90], minor=True)

axs[0].set_ylim([0, 250])
axs[1].set_ylim([0, 110])

axs[0].axvline([2], lw=0.2, c="black", ls="--")
axs[1].axvline([2], lw=0.2, c="black", ls="--")

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XIST_exp_CCLE_male_and_female_histogram.png", dpi=dpi_set, bbox_inches = 'tight')
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/XIST_Males/XIST_exp_CCLE_male_and_female_histogram.pdf", dpi=dpi_set, bbox_inches = 'tight')




