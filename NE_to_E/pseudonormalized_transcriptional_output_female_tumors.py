
import pandas as pd
import collections
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from statannot import add_stat_annotation
import statistics
import scipy
import scipy.stats as stats

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

df_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/male_and_female_XIST_expression_TCGA_rev_Xena_TPM.csv")

df_ref = df_ref[df_ref['Gender']=="FEMALE"]
df = df_ref[df_ref['Classification']!="UNKNOWN"]

df_barcode = df['Barcode'].tolist()

end_barcode = []

for a in df_barcode:
    end_barcode.append(int(split_advanced(a, "-", 3)[1]))

df['End_Val'] = end_barcode

df = df[df['End_Val']<10]

del df['End_Val']

orig_cols = df.columns.tolist()

temp_class = []

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
    temp_second_class = temp_df['Secondary_Class'].tolist()[0]
    
    new_rows.append([str(a)+"-09", temp_log2, temp_average, temp_classification, temp_gender, temp_second_class, "UNKNOWN"]) # -9 is an average val, not used elsewhere

df_from_new_rows = pd.DataFrame(new_rows)

df_from_new_rows.columns = orig_cols

for a in dup_items:
    df = df[~(df['Barcode'].str.contains(a))]

df_ref = pd.concat([df, df_from_new_rows], ignore_index=True)

second_class = df_ref['Secondary_Class'].tolist()

a=0
while a<len(df_ref['Classification'].tolist()):
    if df_ref['Classification'].tolist()[a] == "TGCT" and second_class[a] == second_class[a]:
        temp_class.append(second_class[a])
    elif df_ref['Classification'].tolist()[a] == "TGCT" and second_class[a] != second_class[a]:
        temp_class.append("ERROR")        
    else:
        temp_class.append(df_ref['Classification'].tolist()[a])
    a=a+1

df_ref['new_class'] = temp_class
df_ref = df_ref[df_ref['new_class']!="ERROR"]
        
invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

for a in invalid_ids:
    df_ref = df_ref[~(df_ref['Barcode'].str.contains(a))]

df_high = df_ref['Barcode'].tolist()

class_dict = dict(zip(df_ref['Barcode'].tolist(), df_ref['new_class'].tolist()))








print("READING 1st RNA")

df_RNA_temp = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/gene_id_name_and_chr_all_biomart.csv")

print("READ 1st RNA")

gene_id = df_RNA_temp['Gene stable ID'].tolist()
chromosome = df_RNA_temp['Chromosome/scaffold name'].tolist()
gene_symbol = df_RNA_temp['Gene name'].tolist()

chr_dict = dict(zip(gene_id, chromosome))
symbol_dict = dict(zip(gene_id, gene_symbol))

print("READING 2nd RNA")

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/female_tumors_only_averaged_tcga_RSEM_gene_tpm.txt", sep="\t")

print("READ 2nd RNA")
    
gene_ids_temp_unp = df_RNA['sample'].tolist()

gene_ids_temp = []

for a in gene_ids_temp_unp:
    gene_ids_temp.append(a.split(".")[0])

new_chr = []
new_genes = []

for a in gene_ids_temp:
    try:
        new_chr.append(chr_dict[a])
    except KeyError:
        new_chr.append("UNKNOWN")
    try:
        new_genes.append(symbol_dict[a])
    except KeyError:
        new_genes.append("UNKNOWN")
    
df_RNA['chromosome'] = new_chr
df_RNA['gene_symbol'] = new_genes

df_RNA = df_RNA[df_RNA['chromosome']!="UNKNOWN"]
df_RNA = df_RNA[df_RNA['gene_symbol']!="UNKNOWN"]

df_RNA_temp_cols = df_RNA.columns.tolist()

df_high = list(set(df_high) & set(df_RNA_temp_cols))


df_gene_class = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/chrX_gene_classes.xlsx")

inactivated_genes = df_gene_class['Inactivated'].tolist()
escaping_genes = df_gene_class['Escaping'].tolist()

all_chrX_valid_genes = inactivated_genes + escaping_genes

df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]

#df_autosome = df_RNA[~(df_RNA['chromosome'].isin(['X', 'Y']))]

df_autosome = df_RNA[(df_RNA['gene_symbol'].isin(escaping_genes))]

    
#df_e_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/all_redone_e_pseudo_reference_XISTneg_males.csv")
#df_ne_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/all_redone_ne_pseudo_reference_XISTneg_males.csv")

df_e_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/e_pseudo_reference_XISTpos_females.csv")
df_ne_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/ne_pseudo_reference_XISTpos_females.csv")


#NE_list = ['RPL10', 'RPL39', 'RPL36A', 'SSR4', 'FLNA', 'SAT1', 'PGK1', 'SLC25A5', 'BCAP31', 'NONO', 'BGN', 'MORF4L2', 'PLP2', 'NDUFA1', 'UXT', 'MAGED2', 'PRDX4', 'HSD17B10', 'MSN', 'NDUFB11', 'ATP6AP1', 'RBMX', 'MAGED1', 'PDHA1', 'WDR13', 'TCEAL4', 'LAMP2', 'FUNDC2', 'PGRMC1', 'ATP6AP2', 'GDI1', 'TSC22D3', 'FAM50A', 'SMS', 'IRAK1', 'TIMM17B', 'IDH3G', 'HUWE1', 'HNRNPH2', 'IGBP1', 'WDR45', 'CETN2', 'UBE2A', 'LAGE3', 'ARAF', 'HMGB3', 'MCTS1', 'TSR2', 'USP11', 'FAM3A', 'G6PD', 'OTUD5', 'PSMD10', 'ARMCX6', 'BEX4', 'HTATSF1', 'VBP1', 'FTSJ1', 'LDOC1', 'AIFM1', 'SNX12', 'DKC1', 'HPRT1', 'IL13RA1', 'MID1IP1', 'GLA', 'IDS', 'UBL4A', 'MAGT1', 'HDAC6', 'OGT', 'PRAF2', 'PBDC1', 'RBM10', 'SLC35A2', 'ARMCX3', 'THOC2', 'PRPS1', 'APOO', 'TFE3', 'LAS1L', 'FHL1', 'GRIPAP1', 'PRPS2', 'TAZ', 'STAG2', 'ZDHHC9', 'VAMP7', 'DNASE1L1', 'ACOT9', 'SLC10A3', 'HAUS7', 'PIN4', 'NSDHL', 'DYNLT3', 'RBMX2', 'CUL4B', 'CD99L2', 'ELK1', 'MPP1']

def generate_vals(df_inp, vals_ref, sample):
        
        gene_ids = df_inp['gene_symbol'].tolist()
        
        temp_vals = df_inp[sample[0]].tolist()
                        
        temp_sort_list = sorted(vals_ref, reverse=True)
                
        avg_ratio = []
        
        if "CDK16" in gene_ids:
            subset_specific = "Yes"
        else:
            subset_specific = "No"
        
        a=0
        while a<len(gene_ids):
                        
            cutoff1 = 0
            cutoff2 = 0
                
            if subset_specific != "Yes":
                
                """          
                if gene_ids[a] in NE_list[0:45]: #top 10- genes with greatest median expression across all lineage averages in both XIST- males & XIST+ females
                    avg_ratio.append(math.log(temp_vals[a]/vals_ref[a], 2))
                else:
                    avg_ratio.append(float("nan"))       
            
                """
                
                if vals_ref[a] >= cutoff1: #top 5 genes with greatest median expression across all lineage averages
                    avg_ratio.append(math.log(temp_vals[a]/vals_ref[a], 2))
                else:
                    avg_ratio.append(float("nan"))

                    
            else:

                if vals_ref[a] >= cutoff2: #top 5 genes with greatest median expression across all lineage averages
                    avg_ratio.append(math.log(temp_vals[a]/vals_ref[a], 2))
                else:
                    avg_ratio.append(float("nan"))
                    
                """
                
                if gene_ids[a] in ['RPS4X', 'UBA1', 'EIF2S3', 'DDX3X', 'CDK16', 'EIF1AX', 'KDM5C', 'MAOA', 'RAB9A', 'PNPL4A', 'OFD1', 'SMC1A', 'ARSD', 'FUNDC1', 'USP9X', 'GEMIN8', 'MSL3']: #, 'RAB9A']: #, 'PNPL4A', 'OFD1', 'SMC1A', 'ARSD', 'FUNDC1', 'USP9X', 'GEMIN8', 'MSL3']: #top 20 genes with greatest median expression across all lineage averages in both XIST- males & XIST+ females
                    avg_ratio.append(math.log(temp_vals[a]/vals_ref[a], 2))
                else:
                    avg_ratio.append(float("nan"))   
                """ 
                
            
            a=a+1
        
        return avg_ratio

new_NE_rows = []
new_E_rows = []

print(len(df_high))

a=0
while a<len(df_high):
    
    df_e_vals = df_e_ref[class_dict[df_high[a]]].tolist()
    df_ne_vals = df_ne_ref[class_dict[df_high[a]]].tolist()
        
    X_high = generate_vals(df_X, df_ne_vals, [df_high[a]])
    other_high = generate_vals(df_autosome, df_e_vals, [df_high[a]])
    new_NE_rows.append(X_high)
    new_E_rows.append(other_high)
    
    a=a+1
    #print(a)



#Use 10, 150 correct LUSC & HNSC




df_NE_out = pd.DataFrame(new_NE_rows).T
df_NE_out.columns = df_high
df_NE_out.index = df_X['gene_symbol'].tolist()

df_NE_out.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/skewness_female_tumors_normalized_NE_expression_TPM_geq200.csv")

df_E_out = pd.DataFrame(new_E_rows).T
df_E_out.columns = df_high
df_E_out.index = df_autosome['gene_symbol'].tolist()

df_E_out.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/skewness_female_tumors_normalized_E_expression_TPM_geq3.csv")



