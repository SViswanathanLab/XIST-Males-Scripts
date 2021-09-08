#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 15:57:18 2021

@author: ananthansadagopan
"""

import pandas as pd
import collections

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

df_class = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")
id_ref = df_class['SAMPID'].tolist()
SMTS = df_class['SMTS'].tolist()

class_dict = dict(zip(id_ref, SMTS))

df_gender = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep="\t")

sub_id = df_gender['SUBJID'].tolist()
sex = df_gender['SEX'].tolist()

sex_rev = []

for a in sex:
    if a == 1:
        sex_rev.append("MALE")
    elif a == 2:
        sex_rev.append("FEMALE")

sex_dict = dict(zip(sub_id, sex_rev))


df_new = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/GTEx/XIST_GTEX_RSEM_TPM_expression.txt", sep="\t")

df_new = df_new.T

df_new = df_new.iloc[1:]

df_new.columns = ['XIST_log2TPM']

exp_vals = df_new['XIST_log2TPM'].tolist()

adj_exp = []

for w in exp_vals:
    adj_exp.append(2**w)

df_new['XIST_TPM'] = adj_exp

sample_index = df_new.index.tolist()

gender_vals = []
type_vals = []

for a in sample_index:
    try:
        type_vals.append(class_dict[a])
    except KeyError:
        type_vals.append("UNKNOWN")
    try:
        gender_vals.append(sex_dict[split_advanced(a, "-", 2)[0]])
    except KeyError:
        gender_vals.append("UNKNOWN")

df_new['Classification'] = type_vals
df_new['Gender'] = gender_vals

df_new.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/male_and_female_XIST_expression_GTEx_rev_Xena_TPM.csv", index=True)