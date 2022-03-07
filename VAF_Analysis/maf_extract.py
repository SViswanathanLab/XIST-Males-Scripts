import numpy as np
import csv
import json
import pandas as pd


csv_file = 'input.csv'


df = pd.read_csv(csv_file)

chr_column = np.array(df['Chromosome'])
chr_start_column = np.array(df['Start_Position'])
chr_end_column = np.array(df['End_Position'])
tumour_sample_column = np.array(df['Truncated_Barcodes'])

column_list = len(chr_column)
unique_id = list(np.unique(tumour_sample_column))
print ('n_unique = %s' %len(unique_id))


id_index = []
chrs_start = []
chrs_end = []




for i in range(len(unique_id)):
	id_index = [j for j , elem in enumerate(tumour_sample_column) if elem==unique_id[i]]
	chrs_start += [list(chr_start_column[id_index])]
	chrs_end += [chr_end_column[id_index]]


#Write out all chrom pos to a text file so I can load it into the USCS lib
chr_pos_list = list(chr_start_column)
total_chr_coord = ['chr'+str(chr_column[j])+':'+str(elem1)+'-'+str(elem1) for j, elem1 in enumerate(chr_pos_list)]


hg19_pos = 'textfile.txt'
a1 = np.transpose([total_chr_coord])
np.savetxt(hg19_pos,a1,fmt='%s')




