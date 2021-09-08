#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 07:17:18 2021

@author: ananthansadagopan
"""

import pandas as pd

df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/XIST_expression_summary.xlsx")

sample_list = df['Barcode'].tolist()

df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/sp_specimen_type", sep="\t")

df2_samples = df2['icgc_specimen_id'].tolist()

df2_patient = df2['dcc_specimen_type'].tolist()

patient_dict = dict(zip(df2_samples, df2_patient))

new_patient = []

for a in sample_list:
    try:
        new_patient.append(patient_dict[a])
    except:
        new_patient.append("Unknown")
        
df['Specimen_Type'] = new_patient

print(df)

df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/XIST_expression_summary.csv", index=False)