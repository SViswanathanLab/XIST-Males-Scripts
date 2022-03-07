import csv
import pandas as pd
import numpy as np
import shutil
import os
import openpyxl
import matplotlib as mpl
import matplotlib.pylab as plt



excel_file = 'input.xlsx'


df_tot = pd.read_excel(excel_file,engine='openpyxl')




#Wanted columns ['Truncated_Barcodes','Start_Position','Chromosome']
chr_column = df_tot['Chromosome']
chr_start_column = df_tot['Start_Position']
tumour_sample_column = df_tot['Tumor_Sample_Barcode']

column_list = len(chr_column)

unique_id = list(np.unique(tumour_sample_column))
case_id = [elem[:12] for ind, elem in enumerate(unique_id)]


filename = 'gdc_sample_sheet.2019-06-24.tsv' #This is important

i=0
df = pd.read_csv(filename, sep='\t')
col_list = list(df.columns)  #List column names

sample_type_tum = ['Primary Tumor', 'Metastatic','Primary Blood Derived Cancer - Peripheral Blood']# 
sample_type_normal = ['Blood Derived Normal', 'Solid Tissue Normal'] #Two types of normals we can have


ref_list = list(df['Case ID'])

hg38_dir = '/tcga/tcga_RNASeq/tcga_RNASeq_BAM_FILES/' #Set relative path


sample_list_for_titan = []
i=0
tum_path = []
norm_path = []
total_index_paths = []
norm_case_id = []
tum_filename = []
norm_filename = []
wanted_extension = '01A'
for i in range(len(case_id)):
    if (case_id[i] in ref_list)==False:
        print (i)
        continue 
    #Get the index that matches the case_id
    index_tum = np.where((df['Case ID']==case_id[i]) & ( (df['Sample Type']==sample_type_tum[0]) | (df['Sample Type']==sample_type_tum[1]) | (df['Sample Type']==sample_type_tum[2]) ) )[0] #ref_list.index(case_id[i]) #Get the index for the tumour
    
    if len(index_tum)==0:
        continue
    

    if len(index_tum)>1:
        index_tum = index_tum[0]
        
    
    index_tum = np.int(index_tum)

    #We now want both the file_id and filename
    file_id_tum = df['File ID'][index_tum]
    file_name_tum = df['File Name'][index_tum] #Our BAM file


    # Now get the normal and tumour file paths
    tum_path += [hg38_dir + file_id_tum +'/' + file_name_tum]
    tum_filename += [file_name_tum]
    norm_case_id += [case_id[i]]

'''
#Write out the joing filename + case ID 
output_norm_case_id = 'some_file.txt'

with open(output_norm_case_id, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(norm_case_id,tum_path))
'''




