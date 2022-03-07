
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

cutoff_value = 0

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(10,10)})
sns.set(font_scale=1)
plt.rcParams.update({'font.size': 14})
plt.rcParams['axes.linewidth'] = 0.3
dpi_set = 72 # change the output resolution

full_class_list = ['PanCan']

www=0
while www<len(full_class_list):

    temp_class = full_class_list[www]
    
    df_XIST = pd.read_csv("../Methylation/Other_Input/male_tumors_averaged_TCGA_Xena_TPM.csv")
    
    invalid_ids = ['TCGA-BP-4974','TCGA-EL-A3T3', 'TCGA-GL-7773', 'TCGA-KO-8403', 'TCGA-M9-A5M8', 'TCGA-98-7454', 'TCGA-G3-A5SM']

    for a in invalid_ids:
        df_XIST = df_XIST[~(df_XIST['Barcode'].str.contains(a))]

        
    """
    
    if temp_class != 'TGCT-NS':
    
        df_XIST = df_XIST[df_XIST['Classification']==temp_class]
        
    else:
        
        df_XIST = df_XIST[df_XIST['Secondary_Class']==temp_class]  
        
    """

    #df_XIST = df_XIST[df_XIST['Secondary_Class']!=temp_class] 
    
    if temp_class == "TGCT-NS":
        df_XIST = df_XIST[df_XIST['Secondary_Class']==temp_class] 
    elif temp_class == "PanCan":
        print("no subsetting")
        #df_XIST = df_XIST[df_XIST['Classification'].isin(['LUSC', 'LUAD', 'SKCM', 'LIHC', 'ESCA', 'STAD', ''])] #Pan-lung, SKCM, LIHC, and Pan-GI (top 4)
    
    df_high = df_XIST[df_XIST['XIST_TPM']>=3]
    df_low = df_XIST[df_XIST['XIST_TPM']<3]
    
    h_samples = df_high['Barcode'].tolist()
    
    l_samples = df_low['Barcode'].tolist()
    
    df_high = h_samples
    
    df_low = l_samples
            
    df2 = pd.read_csv("../Methylation/Other_Input/TCGA_mastercalls.abs_tables_JSedit.fixed.processed.txt", sep="\t")
    
    uq_samples = df2['array'].tolist()
    
    chrX_cn = df2['purity'].tolist()
    
    cn_dict = dict(zip(uq_samples, chrX_cn))
    
    pur_high = []
    
    for a in df_high:
        try:
            pur_high.append(cn_dict[a])
        except KeyError:
            print("KEY ERROR " + str(a))
            continue
            
    pur_low = []
    
    for a in df_low:
        try:
            pur_low.append(cn_dict[a])
        except KeyError:
            print("KEY ERROR " + str(a))
            continue
        
    print(statistics.median(pur_high))
    print(statistics.median(pur_low))
    
    print("READING methylation")
    
    df_RNA = pd.read_csv("../Methylation/Other_Input/male_tumors_only_averaged_jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.txt", sep="\t")
    
    print("READ methylation")
    
    df_annot = pd.read_csv("../Methylation/Other_Input/methylation_annotation_ref.csv")
    
    #df_annot = df_annot[df_annot['Feature_Type'] == "Island"] #Restrict to CpG islands
    
    comp_element_ref = df_annot['Composite Element REF'].tolist()
    chromosome = df_annot['Chromosome'].tolist()
    start = df_annot['Start'].tolist()
    end = df_annot['End'].tolist()
    gene_list = df_annot['Gene_Symbol'].tolist()
    cgi_list = df_annot['CGI_Coordinate'].tolist()
    feature_list = df_annot['Feature_Type'].tolist()
    
    chr_dict = dict(zip(comp_element_ref, chromosome))
    start_dict = dict(zip(comp_element_ref, start))
    end_dict = dict(zip(comp_element_ref, end))
    gene_dict = dict(zip(comp_element_ref, gene_list))
    cgi_dict = dict(zip(comp_element_ref, cgi_list))
    feature_dict = dict(zip(comp_element_ref, feature_list))
    
    df_RNA = df_RNA[df_RNA['Composite_Element_REF'].isin(comp_element_ref)]
    
    orig_comp_elements = df_RNA['Composite_Element_REF'].tolist()
    
    chr_new = []
    start_new = []
    end_new = []
    gene_new = []
    cgi_new = []
    feature_new = []
    
    for a in orig_comp_elements:
        chr_new.append(chr_dict[a])
        start_new.append(start_dict[a])
        end_new.append(end_dict[a])
        gene_new.append(gene_dict[a])
        cgi_new.append(cgi_dict[a])
        feature_new.append(feature_dict[a])
    
    print("inserting columns")
    
    df_RNA.insert(loc=0, column='End', value=end_new)
    df_RNA.insert(loc=0, column='Start', value=start_new)
    df_RNA.insert(loc=0, column='Chromosome', value=chr_new)
    df_RNA.insert(loc=0, column='Gene', value=gene_new)
    df_RNA.insert(loc=0, column='CGI_Coordinate', value=cgi_new)
    df_RNA.insert(loc=0, column='Feature_Type', value=feature_new)
    
    #df_RNA = df_RNA[df_RNA['Feature_Type']=="Island"]
    
    df_RNA_temp_cols = df_RNA.columns.tolist()
    
    df_high = list(set(df_high) & set(df_RNA_temp_cols))
    df_low = list(set(df_low) & set(df_RNA_temp_cols))
    
    cpg_promoter_df = pd.read_csv("../Methylation/Output_Files/chrX_non_escaping_CpG_promoter_islands.csv") #corresponding to only non-escaping TSSs, and within +/-125 bp to the TSS
    
    islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()
    
    #df_X = df_RNA[df_RNA['gene_symbol'].isin(inactivated_genes)]
    
    df_X = df_RNA[df_RNA['Composite_Element_REF'].isin(islands_of_interest)]
    
    df_autosome = df_RNA[~(df_RNA['Chromosome'].isin(['X', 'Y']))]
    
    def DGE(df_inp, temp_sample):
        
        sample_val = []
        
        gene_ids = df_inp['Composite_Element_REF'].tolist()
        
        temp_df = df_inp[temp_sample]
        
        temp_vals = temp_df.values.tolist()
        
        temp_total_TPM = 0
        
        a=0
        while a<len(gene_ids):
            
            if temp_vals[a] >= cutoff_value:
                
                temp_total_TPM = temp_total_TPM + temp_vals[a]
                
                sample_val.append(temp_vals[a])
            
            a=a+1
            
        median_val = statistics.median(sample_val)
    
        return median_val, temp_total_TPM
    
    h_median_list = []
    h_total_TPM_list = []
    
    for a in df_high:
        temp_median, total_TPM = DGE(df_X, a)
        h_median_list.append(temp_median)
        h_total_TPM_list.append(total_TPM)
        
    l_median_list = []
    l_total_TPM_list = []
    
    for a in df_low:
        temp_median, total_TPM = DGE(df_X, a)
        l_median_list.append(temp_median)
        l_total_TPM_list.append(total_TPM)
    
    h_median_list_auto = []
    h_total_TPM_list_auto = []
    
    for a in df_high:
        temp_median, total_TPM = DGE(df_autosome, a)
        h_median_list_auto.append(temp_median)
        h_total_TPM_list_auto.append(total_TPM)
    
    l_median_list_auto = []
    l_total_TPM_list_auto = []
    
    for a in df_low:
        temp_median, total_TPM = DGE(df_autosome, a)
        l_median_list_auto.append(temp_median)
        l_total_TPM_list_auto.append(total_TPM)    
        
    XISTpos_samples_ratio_df = pd.DataFrame([df_high, h_median_list, h_median_list_auto]).T
    XISTpos_samples_ratio_df.to_csv("../Methylation/Output_Files/NEW_male_tumors_sample_level_methylation_XISTpos_samples.csv", index=False)

    low_ratio = []
    
    a=0
    while a<len(l_median_list):
        low_ratio.append((2**l_median_list[a]-1)/(2**l_median_list_auto[a]-1))
        a=a+1       

    XISTpos_samples_ratio_df = pd.DataFrame([df_low, l_median_list, l_median_list_auto]).T
    XISTpos_samples_ratio_df.to_csv("../Methylation/Output_Files/NEW_male_tumors_sample_level_methylation_XISTneg_samples.csv", index=False)

    
    www=www+1
    
