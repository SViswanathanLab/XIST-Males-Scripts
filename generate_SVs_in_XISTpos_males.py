#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 08:01:11 2021

@author: ananthansadagopan
"""

import pandas as pd

df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/XIST_expression_summary.xlsx")

df = df[df['Gender']=="male"]
df = df[df['XIST_FPKM_UQ']>=0.5]

high_barcodes = df['Barcode'].tolist()

df2 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/pcawg_consensus_1.6.161116.somatic_svs.xena.sp", sep="\t")

df2 = df2[df2['sample'].isin(high_barcodes)]

#df2 = df2[df2['chr']=="chrX"]
    
df2.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/PCAWG/SVs_in_XISTpos_males.csv", index=False)


#October_2016_all_patients_2778.snv_mnv_indel.maf.xena.nonUS