# Script to plot CD spectra
# Written by Daniel Chabeda
# April 17, 2023

import numpy as np 
from glob import glob
import os
import sys
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# This script takes a cd.dat file as input and returns high res images to the directory
c = 137.03599
H_to_eV = 1/0.036749308136649
freq_to_wvlen = 2*np.pi*c
Bohr_to_nm = 5.2918e-2

def nm_to_eV(nm):
        eV = Bohr_to_nm*freq_to_wvlen*H_to_eV/nm
        return eV

# Set defaults and command line options
plot_absorp = 0
if 'abs' in sys.argv:
    plot_absorp = 1

xtype = 'eV'
if 'nm' in sys.argv:
    xtype = 'nm'

rs_plot = 0
if 'rs' in sys.argv:
    rs_plot = 1
    
if 'flip' in sys.argv:
    sign = -1
else: sign = 1

if 'rslab' in sys.argv:
    label_rs = 1
else:
    label_rs = ''

font = {'weight' : 'light',
       'size' : 24}
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams["axes.prop_cycle"]
mpl.rcParams['axes.linewidth'] = 2

# Determine which files in the directory need to be plotted
cur_dir = str(os.getcwd())

# Make a list of all 'cd0_perturb.dat' files
cd_files = glob('*cd0_perturb.dat')
#print(cd_files, '\n')
# Remove from the list any file that has a corresponding png file already present
del_ctr = [] # used to keep track of the right indices for the evolving cd_files list
del_idx = 0
for i, cd_file in enumerate(cd_files):
    file_head = cd_file.split('.dat')[0]
    #print(file_head)
    png_dir = file_head[:3] + '_figs'
    #print(f"searching in: {png_dir}")
    
    png_file = cur_dir + '/' + png_dir + '/' + file_head + '.png'
    if os.path.exists(png_file):
        del_ctr.append(i - del_idx)
        del_idx += 1
        #print("File exists")
        

# Remove the already visualized files from the list
for i in del_ctr:
    del cd_files[i]


 
for cd_file in cd_files:
    print('Plotting {}'.format(cd_file))
    # Get the name of the cd file to use as a title
    file_head = cd_file.split('.dat')[0]
    # Get the name of the corresponding 'rs' file.
    # TODO: make this logic prettier
    file_stub = file_head.split('_cd')[0]
    if file_head.split('_cd')[1] == "0":
    	file_nib = '0'  
    	continue # do not plot the unperturbed CD spectra
    elif file_head.split('_cd')[1] == "0_perturb":
    
    	file_nib = '0_perturb'
    	match = re.search('\w(\d+)',file_head.split('_cd')[0])
    	if match:
    		theta = match.group(1)
    rs_file = file_stub + '_rs{}.dat'.format(file_nib)
    
    png_file = file_head + '.png'
    cur_dir = os.getcwd()
    dir_name = file_head[:3] + '_figs'
    
    cd_dat = np.loadtxt(cd_file, dtype=np.float64, skiprows=1)
    nm = cd_dat[:,0]; eV = cd_dat[:,1]; cd = cd_dat[:,2]
    if eV.shape[0] == 0: print('Failed to load file {}'.format(cd_file)); exit();
    
    if xtype == 'nm':
        # Plot within range that BSE is relatively accurate
        #nm_min = nm.max()-90
        #nm_max = nm.max()
        nm_min = 340
        nm_max = 750
        mindex = 0 
        maxdex = np.where(nm < nm_min)[0][0]
        xvals = nm[mindex:maxdex]
        xrange = np.array([nm[mindex],nm[maxdex]])
        xlab = 'Wavelength (nm)'
    elif xtype == 'eV':
        eV_max = eV.min()+0.4
        mindex = 0
        try:
            maxdex = np.where(eV > eV_max)[0][0]
        except IndexError:
            maxdex = np.where(eV.max())[0][0]
            print('The maxdex {}   eV[maxdex] {}  eVmax {}'.format(maxdex, eV[maxdex], eV.max()))
        xvals = eV[mindex:maxdex]
        xrange = np.array([eV[mindex],eV[maxdex]])
        xlab = 'Energy (eV)'
    
    if plot_absorp:
        try:
        	absorp = cd_dat[mindex:maxdex,3] 
        except IndexError:
            print("Warning: No absorption data provided. Absorption will not be plotted.")
            plot_absorp = 0
            pass
            
    # Create plt objects
    fig, ax = plt.subplots()
    fig.set_size_inches(9,6)
    
    leg_lab_list = ['CD']
    
    ax.locator_params(axis='both',nbins=4)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(top=True, bottom=True, left=True, right=True)

    ax.tick_params(axis='both',which='major', direction='in',length=10,width=2.0,color='black',pad=5,labelsize=24,labelcolor='black',
            labelrotation=0)
    ax.tick_params(axis='both',which='minor', right=True, top=True, direction='in',length=5,width=2.0,color='black',pad=15,labelsize=10,labelcolor='black',
            labelrotation=0)

    
    # PLOTTING
    yvals = cd[mindex:maxdex]
    line1, = ax.plot(xvals, sign*yvals, color='k', linewidth=2.5,alpha=1.0, label='CD') # Plot CD spectrum
    ax.text(450,2.5e-4, str(theta) + r'$^\circ$', fontsize=24)
    ax.hlines(y=0.0,xmin=xrange[0],xmax=xrange[1],color='k', linewidth=1.0,alpha=0.7) # Plot horizontal line at 0.0
    leg_art_list = [line1]
    
    if rs_plot:
        rs_dat = np.loadtxt(rs_file, dtype=np.float64)
        rs_lab = rs_dat[:,0]; rs = rs_dat[:,-1]; rs_ene = rs_dat[:,-2]
        
        if xtype == 'nm':
        	rs_ene = rs_ene[rs_ene < Bohr_to_nm*freq_to_wvlen/nm[maxdex]]
        	#print(rs_ene)
        	#rs_ene = rs_ene[:10]
        if xtype == 'eV':
        	rs_ene = rs_ene[rs_ene < eV[maxdex]/H_to_eV]
        rs = rs[:rs_ene.shape[0]]
        #print('max rs, ', rs.max(), '\n')
        yvals_max = max(yvals.min(), yvals.max(), key=abs)
        #print('max y: ', yvals_max, '\n')
        scalef = yvals_max/rs.max()
        if scalef < 0:
            scalef *= -1
        
        for i, ene in enumerate(rs_ene):
            if xtype == 'nm':
            	rs_loc = Bohr_to_nm*freq_to_wvlen/ene
            if xtype == 'eV':
            	rs_loc = H_to_eV * ene
            ax.vlines(rs_loc, scalef*min(0,sign*rs[i]), scalef*max(0,sign*rs[i]), color='k', linewidth=1.0, alpha = 1.0)
            if label_rs:
                ax.text(rs_loc, scalef*max(0,sign*rs[i]), str(int(rs_lab[i])), fontsize=12, verticalalignment='top')
            
    if plot_absorp: # If a plot of absorption is requested, plot on twin x axis
        ax2 = ax.twinx()
        ax2.locator_params(axis='y',nbins=6)
        ax2.tick_params(axis='y',which='major', direction='in',length=10,width=2.0,color='k',pad=5,labelsize=24,labelcolor='k',
               labelrotation=0)
        ax2.tick_params(axis='y',which='minor', direction='in',length=5,width=2.0,color='k',pad=15,labelsize=10,labelcolor='black',
               labelrotation=0)
        if xtype == 'nm':
            line2, = ax2.plot(nm[mindex:maxdex], absorp, color='k', linewidth=1.5, linestyle='dashed', alpha=0.75, label='Absorption')
        if xtype == 'eV':
        	line2, = ax2.plot(eV[mindex:maxdex], absorp, color='k', linewidth=1.5, linestyle='dashed', alpha=0.75, label='Absorption')
        ax2.set_ylabel('Absorption (a.u)', fontsize=24, color='k')
        ax2.set_ylim(-0.002,0.1)
        leg_lab_list.append('Absorption')
        leg_art_list.append(line2)
    
    
    ax.set_title(file_head, pad=20)
    ax.set_xlabel(xlab, fontsize=24)
    ax.set_ylabel('CD (a.u)', fontsize=24, color='k')
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    ax.legend(leg_art_list, leg_lab_list, loc='upper right', frameon=False, fontsize=20)
    ax.set_xlim(430,650)
    ax.set_ylim(-4.2e-4,4.2e-4)

    # Make directory to store images
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    
    os.chdir(dir_name)
    plt.savefig(png_file, dpi=200, format='png', bbox_inches = 'tight', pad_inches=0.2)
    os.chdir(cur_dir)
    