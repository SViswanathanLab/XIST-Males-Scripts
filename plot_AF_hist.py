import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
from matplotlib.pyplot import cm
from matplotlib.gridspec import GridSpec
from matplotlib import rc
import glob
import pandas as pd
from scipy.stats import kde

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Arial"
ms=1.5
lw = 2.0
axis_labelsize = 15
axis_fontsize = 15



joint_file_name = 'TCGA_joint_ASE.txt'
data_joint = np.loadtxt(joint_file_name,dtype='str')
sample_list = data_joint[:,0]
bam_list = data_joint[:,1]
file_list_chrX = data_joint[:,2]
num_file_chrX = len(file_list_chrX)

#Load in the normalised chrY values
chrY_filename = 'normalised_chrY_autosome_list.txt'
norm_chrY_data = np.loadtxt(chrY_filename,dtype='str')
norm_chrY_values = norm_chrY_data[:,0]
sample_id_chrY = norm_chrY_data[:,1]

#Load in the het fraction
het_filename = 'sample_het_data.txt'
het_data = np.genfromtxt(het_filename,dtype='str')
sample_id_het = het_data[:,0]
het_values = het_data[:,1] #np.array([np.float(i) for i, elem in enumerate(list(het_data[:,1]))])


i=0
het_list = []
sample_id_list = []
for i in range(num_file_chrX):
    #Load in the text files
    data_file = np.genfromtxt(file_list_chrX[i])
    sample_id = sample_list[i]
    sample_id_list += [sample_id]
    ref_counts = data_file[1:,5]
    tot_counts = data_file[1:,7]
    #Now we need the corresponding chrY map values
    sample_index = np.where(sample_id_chrY==sample_id)[0] #Find the index corresponding to the chrY values
    norm_chrY_val = norm_chrY_values[sample_index]

    #Find indexes where total counts !=0
    non_zero_ind = np.where(tot_counts > 15.)[0]
    
    ref_counts_non_zero = ref_counts[non_zero_ind]
    tot_counts_non_zero = tot_counts[non_zero_ind]
    chrX_pos = data_file[1:,1][non_zero_ind]/1e6
    chrX_AF = ref_counts_non_zero/tot_counts_non_zero


    het_site = np.where((chrX_AF>0.2) & (chrX_AF<0.8))[0]
    het_val = np.float(len(het_site)/len(chrX_AF))
    het_frac = (len(het_site)/len(chrX_AF)) * 100.0
    het_list += [len(het_site)/len(chrX_AF)]

    #get het values between [0.05,0.95] in order to plot a density distribution
    het_dens_id = np.where((chrX_AF>0.05) & (chrX_AF<0.95))[0]
    het_dens_vals = chrX_AF[het_dens_id]
    het_density = kde.gaussian_kde(het_dens_vals)
    het_x_plot = np.linspace(0.05,0.95,300)
    het_density_plot = het_density(het_x_plot)

    
    fig = plt.figure(figsize=(8,10))
    #plt.ion()
    gs = fig.add_gridspec(3,6,hspace=0.4,wspace=0.7)
    ax1 = fig.add_subplot(gs[1:3,:4])
    ax2 = fig.add_subplot(gs[0,:3])
    ax3 = fig.add_subplot(gs[0,3:6])
    ax4 = fig.add_subplot(gs[1:3,4:6])
    ax1.plot(chrX_pos,chrX_AF,c='b',linestyle='',ms=ms,marker='o')
    ax2.hist(norm_chrY_values.astype('float32'),bins=40,color='blue',align='mid')
    ax2.axvline(norm_chrY_val.astype('float32'),color='red',lw=3)
    ax3.hist(het_values.astype('float32'),bins=40,color='grey')
    ax3.axvline(het_val,color='red',lw=3)
    ax4.hist(het_dens_vals,bins=20,color='green',orientation='horizontal')
    #ax4.plot(het_density_plot,het_x_plot,color='green',lw=3)
    #plt.setp(ax2.get_xticklabels(), visible=False)
    

    #fig, axes = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False, gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(5,5))#,figsize=(6,4) inside to change the plot image
    #(ax1)= axes

    #ax1.plot(chrX_pos,chrX_AF,c='b',linestyle='',ms=ms,marker='o')

    ax1.set_ylim(-0.03,1.03)
    ax1.set_xlim(np.min(chrX_pos), chrX_pos[-1])
    ax1.set_ylabel(r'AF $\left(\frac{ref}{ref + alt}\right)$',fontsize=axis_fontsize+2)
    ax1.set_xlabel('chrX Coordinate (Mb)',fontsize=axis_fontsize+2)
    ax1.minorticks_on()
    ax1.tick_params(which="both",bottom=True, top=True, left=True, right=True,direction='in',width=1.5,labelsize=axis_labelsize)
    ax1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False) #Determines if you label them

    ax2.set_xlim(0.0, norm_chrY_values[-1].astype('float32'))
    ax2.set_ylabel('N',fontsize=axis_fontsize+1)
    ax2.set_xlabel('ChrY Mapped Reads',fontsize=axis_fontsize+1)
    ax2.minorticks_on()
    ax2.tick_params(which="both",bottom=False, top=True, left=True, right=True,direction='in',width=1.5,labelsize=axis_labelsize)
    ax2.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)

    #ax3.set_xlim(0.0, 200.0)
    #ax3.set_ylabel('N',fontsize=axis_fontsize+1)
    ax3.set_xlabel('Het fraction',fontsize=axis_fontsize)
    #ax3.minorticks_on()
    ax3.tick_params(which="both",bottom=False, top=True, left=True, right=True,direction='in',width=1.5,labelsize=axis_labelsize)
    ax3.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)

    #ax4.set_xlabel(r'$\rho_{AF}$',fontsize=axis_fontsize+1)
    ax4.set_xlabel('$N_{het}$',fontsize=axis_fontsize+1)
    ax4.set_xlim(0.0, 200.0)
    ax4.tick_params(which="both",bottom=True, top=True, left=True, right=True,direction='in',width=1.5,labelsize=axis_labelsize)
    ax4.tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=False)


    #plt.tight_layout()
    #plt.show()
    plt.savefig('/czlab/inasim/TCGA_data/X_ist/pdf_plots/revised_'+file_list_chrX[i] + '_chrX_AF_hist.pdf',format='pdf',bbox_inches='tight')
    plt.close()

'''
#Write out the het sited with the corresponding sample_ids
output_het_file = 'sample_het_data.txt'

a1 = sorted(zip(sample_id_list,het_list))
np.savetxt(output_het_file, a1, fmt = "%s")
'''





