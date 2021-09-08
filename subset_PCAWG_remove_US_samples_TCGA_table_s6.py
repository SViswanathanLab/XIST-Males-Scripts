#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 21:38:02 2021

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

df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/XIST_expression_summary.xlsx")

df_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/project_code_sp", sep="\t")

ids_ref = df_ref['icgc_specimen_id'].tolist()
proj_code = df_ref['dcc_project_code'].tolist()

id_dict = dict(zip(ids_ref, proj_code))

df_barcode = df['Barcode'].tolist()

vals = []

for a in df_barcode:
    if id_dict[a][-2:] == "US":
        vals.append("Yes")
    else:
        vals.append("No")
        
df['In_TCGA'] = vals


df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/table_S6_PCAWG.csv", index=False)