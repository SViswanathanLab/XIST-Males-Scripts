#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 30 22:06:05 2021

@author: ananthansadagopan
"""

import pandas as pd
import collections

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv", sep="\t")

cols = df.columns.tolist()

del cols[0]

trunc_cols = []
end_vals = []

for a in cols:
    temp_val = split_advanced(a, "-", 4)
    trunc_cols.append(temp_val[0][:-1])
    end_vals.append(int(temp_val[0][:-1][-2:]))

dup_items = list(set([item for item, count in collections.Counter(trunc_cols).items() if count > 1]))

for a in dup_items:
    temp_cols = []
    for q in cols:
        if a in q:
            temp_cols.append(q)
    
    df[a] = df[temp_cols].mean(axis=1)
    
    for c in temp_cols:
        del df[c]
        
    print(a)
        
cols_new = df.columns.tolist()
del cols_new[0]

trunc_cols_new = ['Composite_Element_REF']

for a in cols_new:
    temp_val = split_advanced(a, "-", 4)
    trunc_cols_new.append(temp_val[0][:-1])
    
df.columns = trunc_cols_new

print(len(trunc_cols_new))
print(len(list(set(trunc_cols_new))))

df_gender = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

df_male = df_gender[df_gender['Gender'] == "FEMALE"]
male_barcodes_unp = df_male['Barcode'].tolist()

male_barcodes = []

for a in male_barcodes_unp:
    if a in trunc_cols_new:
        male_barcodes.append(a)


male_barcodes.insert(0, "Composite_Element_REF")

df = df[male_barcodes]

df_barcode = df.columns.tolist()

del df_barcode[0]

end_barcode = []

print("END BARCODE GENERATION")

for a in df_barcode:
    temp_val = int(split_advanced(a, "-", 3)[1])
    end_barcode.append(temp_val)

subset_1 = ['Composite_Element_REF']

a=0
while a<len(df_barcode):
    if end_barcode[a]<10:
        subset_1.append(df_barcode[a])
    a=a+1

df = df[subset_1]

curr_barcode = df.columns.tolist()

del curr_barcode[0]

trunc_barcode = []

for a in curr_barcode:
    trunc_barcode.append(split_advanced(a, "-", 3)[0])

dup_items = list(set([item for item, count in collections.Counter(trunc_barcode).items() if count > 1]))

uq_dups = list(set(dup_items))

print("DUP AVERAGING")

a=0
while a<len(uq_dups):
    b=0
    rel_barcodes = []
    while b<len(curr_barcode):
        if uq_dups[a] in curr_barcode[b]:
            rel_barcodes.append(curr_barcode[b])
        b=b+1
            
    df[uq_dups[a] + "-09"] = df[rel_barcodes].mean(axis=1)
    
    for c in rel_barcodes:
        del df[c]

    print(a)    
    a=a+1

print(df)

df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_normals_only_averaged_jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.txt", sep="\t", index=False)