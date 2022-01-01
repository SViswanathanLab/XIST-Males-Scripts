#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 20:54:59 2021

@author: ananthansadagopan
"""
import pandas as pd
import collections
import math
import numpy as np
import statistics
import collections

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

df_gender = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

df_male = df_gender[df_gender['Gender'] == "MALE"]
male_barcodes = df_male['Barcode'].tolist()

male_barcodes.insert(0, "sample")

print("READING RNA")

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/tcga_RSEM_gene_tpm", sep="\t")

print("READ RNA")

df = df[male_barcodes]

print(df)

numeric_cols = [col for col in df if df[col].dtype.kind != 'O']

print("CONVERTING VALS")

df[numeric_cols] = df[numeric_cols].apply(lambda x: 2**x)

print("CONVERTED VALS")

print(df)

gene_id = df['sample'].tolist()

orig_index_length = len(gene_id)

print(orig_index_length)

df.drop_duplicates(subset=['sample'], keep=False)

new_index_length = len(df.index)

print(new_index_length)

if orig_index_length != new_index_length:
    print("ERROR")

df_barcode = df.columns.tolist()

del df_barcode[0]

end_barcode = []

print("END BARCODE GENERATION")

for a in df_barcode:
    temp_val = int(split_advanced(a, "-", 3)[1])
    end_barcode.append(temp_val)

subset_1 = ['sample']

a=0
while a<len(df_barcode):
    if end_barcode[a] < 10:
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

df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_tumors_only_averaged_tcga_RSEM_gene_tpm.txt", sep="\t", index=False)