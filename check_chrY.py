import pysam
import numpy as np
import os, glob
import matplotlib.pylab as plt
import matplotlib as mpl
from matplotlib.pyplot import cm
from matplotlib import rc

#This will be a file to check if a bam file has a Y chromosome.

#To do this we need a few steps:
#Step 1: Load in all of the Bam files we need
#Step 2:Determine the number of lines within the Bam file
#Step 3: If the number of lines within a bam exceeds a certain number - the bam contains a chrY otherwise chrY is deleted


file_list = np.genfromtxt('TCGA_norm_joint_list.txt',dtype='str') #glob.glob('*.bam') #Gets all files ending with .bam
num_file = len(file_list)
id_list = list(file_list[:,1])


#Loop gets the number of lines in the .bam file for chrY
#We also want to normalise the reads against the total number of reads
autosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']


count_list = np.zeros(num_file)
for i in range(num_file):
	x1 = 0
	samfile = pysam.AlignmentFile(file_list[i,0], "rb")
	chrX_counts = samfile.count() - samfile.count("chrX") - samfile.count("chrY") #samfile.count("chrX") 
	iter = samfile.fetch("chrY")
	for x in iter:
	    if (x.mapping_quality>30.):
                x1 += 1

	
	count_list[i] = x1/chrX_counts	
	print (i)


median_count = np.median(count_list)
count_IQR = np.percentile(count_list,75) - np.percentile(count_list,25)
count_Q1 = np.percentile(count_list,25)
outlier = count_Q1 - (1.5*(count_IQR))
chrY_cond = 2e3 # roughly 2000 for map qual = 20, 150k for the males

good_ind = [] #good meaning containing chrY
bad_ind = []
   
i=0
for i in range(num_file):
	if (count_list[i]>10.*chrY_cond):
		good_ind += [i]
	else:
	    bad_ind += [i]	

'''
#Now plot the histogram
plt.ion()
mpl.rcParams['axes.linewidth'] = 1.5
ms=3
lw = 2.0
axis_labelsize = 15
axis_fontsize = 16
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(5,5))#,figsize=(6,4) inside to change the plot image
(ax1)= axes

n, bins, patches = plt.hist(count_list, 20, facecolor='b', alpha=0.75)
plt.grid(True) # False if we want to remove the grid
plt.xlim(0.0,)
plt.ylabel('N',fontsize=axis_fontsize+2)
plt.xlabel('Normalised chrY reads',fontsize=axis_fontsize+3)
plt.tick_params(which="both",bottom=True, top=True, left=True, right=True,direction='in',width=1.5,labelsize=axis_labelsize)
plt.savefig('fig1.pdf',format='pdf',bbox_inches='tight')
plt.show()
plt.close()
'''

#Print out the filenames with the smallest normalised reads
a1 = sorted(zip(count_list,id_list))
np.savetxt('normalised_chrY_autosome_list.txt', a1, fmt = "%s")


