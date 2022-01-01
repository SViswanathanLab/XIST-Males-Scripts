import numpy as np
import csv
import json
import pandas as pd

def find_median_string(string):
	median_string = string[round(len(string)/2):]
	return median_string



csv_file = 'file.csv'

df = pd.read_csv(csv_file)

chr_column = df['Chromosome']
chr_start_column = df['Start_Position']
chr_end_column = df['End_Position']
tumour_sample_column = df['Truncated_Barcodes']

#Load in the bed file with hg38 coordinate points: 
hg38_bed_file = list(np.loadtxt('hg38_muts_REVISED_males_unaveraged_TPM_geq2_edit.bed',dtype='str'))

new_hg38_chr_list = []
i=0
for i in range(len(hg38_bed_file)):
	pos = hg38_bed_file[i]
	for chr_len in range(len(pos)):
		if pos[chr_len]=='-':
			print (pos[chr_len+1:])
			new_hg38_chr_list += [pos[chr_len+1:]]

		else:
		    continue	



#We can first convert to a list and remove the chrX:
#Then we can obtain the median position which will always be the '-' and remove that index as well
#new_hg38_chr_list = [find_median_string(elem[6:]) for i, elem in enumerate(hg38_bed_file)]

#Now add this to the pandas data frame as a new column: hg38_start hg38_end
df['hg38_start'] = new_hg38_chr_list
df['hg38_end'] = new_hg38_chr_list

df.to_csv(csv_file)


