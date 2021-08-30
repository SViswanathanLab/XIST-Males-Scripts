import numpy as np
import csv
import json
import pandas as pd


#excel_file = 'LUSC_and_TGCT_purity_0.5_cutoff_non_escaping_chrX_coding_variants_TPM_geq10.xlsx'
#excel_file = 'TGCT_and_PanCan_MAF_DNA_VAF_only_TPM_geq_2.xlsx'
#excel_file = 'male_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2.xlsx'
#excel_file = 'female_non_escaping_chrX_somatic_MAF_10_lineage_TPM_geq2.xlsx'
#csv_file = 'TCGA-HC-7740_complete_MAF_TPM_geq2.csv'
csv_file = 'chrX_muts_REVISED_males_unaveraged_TPM_geq2.csv'

#df = pd.read_excel(excel_file,engine='openpyxl')
df = pd.read_csv(csv_file)

#Wanted columns ['Truncated_Barcodes','Start_Position','Chromosome']
chr_column = np.array(df['Chromosome'])
chr_start_column = np.array(df['Start_Position'])
chr_end_column = np.array(df['End_Position'])
tumour_sample_column = np.array(df['Truncated_Barcodes'])

column_list = len(chr_column)


#Now that we have the required chromosomal positions - we need to figure out how to run ASEreadcounter on these regions
#We first want all of the unique elemts in the tumour_sample_column:
#For each unique ID: we want to get the corresponding start and end position 
#So we can loop over each unique id: retuning their rows and hence we can extract the start and end pos

unique_id = list(np.unique(tumour_sample_column))
print ('n_unique = %s' %len(unique_id))


id_index = []
chrs_start = []
chrs_end = []


#output_norm_case_id = 'maf_sample_id_chr_pos.json'
#sample_id = list(unique_id)

#Now use enumerate to concatenate the 
#sample_id = [new_elem[:-16] for ind, new_elem in enumerate(sample_id)]

for i in range(len(unique_id)):
	id_index = [j for j , elem in enumerate(tumour_sample_column) if elem==unique_id[i]]
	chrs_start += [list(chr_start_column[id_index])]
	chrs_end += [chr_end_column[id_index]]


#Save the dictionary as a json file
#res = {unique_id[i]: chrs_start[i] for i in range(len(unique_id))}

#with open(output_norm_case_id, "w") as outfile: 
#    json.dump(res, outfile)



#Write out all chrom pos to a text file so I can load it into the USCS lib
chr_pos_list = list(chr_start_column)
total_chr_coord = ['chr'+str(chr_column[j])+':'+str(elem1)+'-'+str(elem1) for j, elem1 in enumerate(chr_pos_list)]

'''
hg19_pos = 'hg19_pos_chrX_muts_REVISED_males.txt'
a1 = np.transpose([total_chr_coord])
np.savetxt(hg19_pos,a1,fmt='%s')
'''



