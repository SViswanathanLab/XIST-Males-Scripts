import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
from matplotlib.pyplot import cm
from matplotlib.gridspec import GridSpec
from matplotlib import rc
import glob
import pandas as pd


#Load in the two file 

#Load in the normalised chrY values
chrY_filename = 'new_normalised_chrY_autosome_list.txt'
norm_chrY_data = pd.read_csv(chrY_filename,delimiter=' ',header=None)
norm_chrY_values = np.array(norm_chrY_data.loc[:,0])
sample_id_chrY = np.array(norm_chrY_data.loc[:,1])

#Load in the het fraction
het_filename = 'new_sample_het_data_corrected.txt'
het_data = pd.read_csv(het_filename,delimiter=' ',header=None)
sample_id_het = np.array(het_data.loc[:,0])
het_values = np.array(het_data.loc[:,1])


ks_samples = ['TCGA-M9-A5M8','TCGA-98-7454','TCGA-G3-A5SM']
fm_samples = ['TCGA-BP-4974','TCGA-EL-A3T3','TCGA-GL-7773','TCGA-KO-8403']

#Now match the values 
plt.ion()
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(5,5))#,figsize=(6,4) inside to change the plot image
(ax1)= axes
mpl.rcParams['axes.linewidth'] = 1.5
ms=1.5
lw = 2.0
axis_labelsize = 12
axis_fontsize = 12
i=0
ks_index = []
female_index = []
for i in range(len(sample_id_chrY)):
	sample = sample_id_chrY[i]
	chrY_reads = norm_chrY_values[i]
	sample_index = np.where(sample_id_het==sample)[0]
	sample_index = np.int(sample_index)
	sample_check = sample_id_het[sample_index]

	if (sample!=sample_check):
		print ('bad')

	#Now get the het values
	het_frac = het_values[sample_index]
	colour = 'black'
	label = 'Xist+ Male'
	if (sample in ks_samples):
		colour = 'red'
		label='Kleinfelter'
		ks_index += [i]
		print (sample)

	if (sample in fm_samples):
	   colour = 'green'	
	   label='Female'
	   female_index += [i]


	#Now we will use a scatter plot
	ax1.loglog(chrY_reads,het_frac,linestyle='',marker='o',ms='8',c=colour)



ax1.set_xlabel('Normalised reads mapped to chrY',fontsize=axis_fontsize)
ax1.set_ylabel('Normalised Heterozygous sites on chrX',fontsize=axis_fontsize)
ax1.minorticks_on()
ax1.tick_params(which="both",bottom=True, top=True, left=True, right=True,direction='in',width=1.5,labelsize=axis_labelsize)
ax1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
#plt.legend(loc='lower left')
plt.tight_layout()
plt.savefig('/home/imi/Documents/harvard_postdoc/plots/year1/week4_august/het_chrY.pdf',format='pdf',bbox_inches='tight')

plt.show()





