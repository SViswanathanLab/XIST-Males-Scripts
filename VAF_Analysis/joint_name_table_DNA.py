import csv
import pandas as pd
import numpy as np
import shutil
import os
import openpyxl
import matplotlib as mpl
import matplotlib.pylab as plt





excel_file = 'file.xlsx'

df_tot = pd.read_excel(excel_file,engine='openpyxl')

#Wanted columns ['Truncated_Barcodes','Start_Position','Chromosome']
chr_column = df_tot['Chromosome']
chr_start_column = df_tot['Start_Position']
#tumour_sample_column = df_tot['Truncated_Barcodes']

tumour_sample_column = df_tot['Tumor_Sample_Barcode']

unique_id = list(np.unique(tumour_sample_column))

case_id = [elem[:12] for ind, elem in enumerate(unique_id)]



filename = 'tcga_WXS_BAM_HG38_FILES_gdc_sample_sheet.tsv' #This is important

i=0
df = pd.read_csv(filename, sep='\t')
col_list = list(df.columns)  #List column names
sample_type_tum = ['Primary Tumor', 'Metastatic','Primary Blood Derived Cancer - Peripheral Blood']# 
sample_type_normal = ['Blood Derived Normal', 'Solid Tissue Normal'] #Two types of normals we can have

ref_list = df['Case ID']

hg38_dir = '/tcga-artemis/tcga_WXS/tcga_WXS_BAM_HG38_FILES/' #Set relative path


sample_list_for_titan = []
i=0
tum_path = []
norm_path = []
total_index_paths = []
norm_case_id = []
tum_filename = []
norm_filename = []
for i in range(len(case_id)):
    #Get the index that matches the case_id
    #First check if the reference exists - if not add ['0'] so we can filter out in post-processing
    if len(np.where(df['Case ID']==case_id[i])[0])<=1:
        print (i)
        continue

    index_tum = np.where((df['Case ID']==case_id[i]) & ( (df['Sample Type']==sample_type_tum[0]) | (df['Sample Type']==sample_type_tum[1]) | (df['Sample Type']==sample_type_tum[2]) ) )[0] #Get the index for the tumour



    #We now want both the file_id and filename
    file_id_tum = np.array(df['File ID'][index_tum])[0]
    file_name_tum = np.array(df['File Name'][index_tum])[0] #Our BAM file

    # Now get the normal and tumour file paths
    tum_path += [hg38_dir + file_id_tum +'/' + file_name_tum]

    total_index_paths += [hg38_dir + file_id_tum +'/']

    norm_case_id += [case_id[i]]
    tum_filename += [file_name_tum]



'''
#Write out the joing filename + case ID 
output_norm_case_id = 'some_file.txt'

with open(output_norm_case_id, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(norm_case_id,tum_path))
'''    
