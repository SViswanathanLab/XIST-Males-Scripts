import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.pyplot import cm
from matplotlib.gridspec import GridSpec
from matplotlib import rc
import glob
import os
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Arial"

#plt.rcParams['font.serif'] = ['Arial'] + plt.rcParams['font.serif']
ms=1.5
lw = 2.0
axis_labelsize = 12
axis_fontsize = 12



#Load in the normalised reads mapped to the chrY region
male_chrY_filename = 'new_normalised_chrY_autosome_list.txt'
male_norm_chrY_data = np.loadtxt(male_chrY_filename,dtype='str')
male_norm_chrY_values = male_norm_chrY_data[:,0]
male_sample_id_chrY = male_norm_chrY_data[:,1]

female_chrY_filename = 'normalised_female_chrY_autosome_list.txt'
female_norm_chrY_data = np.loadtxt(female_chrY_filename,dtype='str')
female_norm_chrY_values = female_norm_chrY_data[:,0]
female_sample_id_chrY = female_norm_chrY_data[:,1]

#Now plot both distributions using histograms
nbins=40

#ax1.hist(male_norm_chrY_values.astype('float32'),bins=nbins,align='mid',color='r',alpha=0.5,label='Male')
#ax1.hist(female_norm_chrY_values.astype('float32'),bins=nbins,align='mid',color='magenta',alpha=0.5, label='Female')
male_hist, male_bins, _ = plt.hist(male_norm_chrY_values.astype('float32'), bins=nbins)
female_hist, female_bins, _ = plt.hist(female_norm_chrY_values.astype('float32'), bins=nbins)

logbins_male = np.logspace(np.log10(male_bins[0]),np.log10(male_bins[-1]),len(male_bins))
logbins_female = np.logspace(np.log10(female_bins[0]),np.log10(female_bins[-1]),len(female_bins))
plt.close()

fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(5,5))#,figsize=(6,4) inside to change the plot image
(ax1)= axes

ax1.hist(male_norm_chrY_values.astype('float32'),bins=logbins_male,align='mid',color='b',alpha=0.8,label='Xist+ Male, $N=201$')
ax1.hist(female_norm_chrY_values.astype('float32'),bins=logbins_female,align='mid',color='magenta',alpha=0.6, label='Female, $N=201$')
#ax1.axvline(x=2.5e-4, color='black', linestyle='dotted')

ax1.set_xlim(female_norm_chrY_values[0].astype('float32'), male_norm_chrY_values[-1].astype('float32'))
ax1.set_xscale('log')
plt.legend(loc='upper left',fontsize=12)

ax1.set_ylabel('N',fontsize=axis_fontsize+1)
ax1.set_xlabel('ChrY Mapped Reads',fontsize=axis_fontsize+1)
ax1.minorticks_on()
ax1.tick_params(which="both",bottom=False, top=True, left=True, right=True,direction='in',width=1.5,labelsize=axis_labelsize)
ax1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)

#plt.savefig('mapped_reads_comp_hist.pdf',format='pdf',bbox_inches='tight')
plt.show()

