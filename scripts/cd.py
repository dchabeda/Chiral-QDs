# Script to calculate CD spectra from dipole and magnetic dipole output files
# Written by Daniel Chabeda
# March 31, 2023

import numpy as np 
import sys, copy, shutil 
import os as opsys

# This code has two analogous functionalities.
# It can calculate CD spectra for single particle transitions 
# in addition to two-body excitonic states.
# The script will take two files as inputs: 
# rs0.dat and rs.dat.
# The corresponding single particle/excitonic CD will be returned as 
# cd0.dat and cd.dat

# The CD formula is: \Delta \varepsilon = \sum_{ai}' \frac{16 \pi N_A}{3 \hbar c} \omega_{ai} \delta(\omega - \omega_{ai})\text{Im}|\mu_{ai}\cdot m_{ai}|
# The two particle CD is analogous, with a single sum over excitonic states replacing the unique double sum over electrons and holes
if 'fig' in sys.argv:
    filename = sys.argv[sys.argv.index('fig')-1]
    fig_flag = 1
else:
    fig_flag = 0
    
lorentz = lambda w, w0, g : ( g**2/((w - w0)**2 + g**2))  
absorp = lambda w, w0, g, u : w * ( g**2/((w - w0)**2 + g**2))  * u
absorp_g = lambda w, w0, g, u: w * u * np.exp( (-0.5 * (w - w0)**2)/(g**2) ) 

def main():
    '''Main function for computing CD spectra. Usage: python cd.py ['spinor' or 'nospinor']'''
    
    print_flag = 1
    if ('s' in sys.argv) or ('silent' in sys.argv):
        print_flag = 0
        
    if print_flag: print('\nRunning cd.py...')
    
    cwd = opsys.getcwd()
    
    # Set constants used in this program

    c = 137.03604
    N_A = 1
    hbar = 1
    h = 2*np.pi
    eps = 6.2

    H_to_eV = 27.211
    freq_to_wvlen = 2*np.pi*c
    Bohr_to_nm = 5.29177e-2

    sign = 1
    
	
    if 'l' in sys.argv:
        lineshape = 'l'
    else:
        lineshape = 'g'
        if print_flag: print('\nNo lineshape specified. Using default Gaussian lineshape.')

    if 'abs' in sys.argv:
        calc_abs = 1
    else:
        calc_abs = 0
    
    if 'm_abs' in sys.argv:
        m_abs = 1
        calc_abs = 1 
    else:
        m_abs = 0
            
     
    wid = 50 # in meV
    lwid = (wid/1000) / H_to_eV
    prefactor = (16 * np.pi) / (3 * hbar * c)
    

    rs_files = ("rs0.dat", "perturb_rs0.dat")
    
    if calc_abs:
        os_files = ('OS0.dat', 'perturb_OS0.dat')
    
    if m_abs:
        m_files = ('M0.dat', 'perturb_m0.dat')

    for i in range(len(rs_files)):
        if i == 0:
            flag = 'single particle'
            fflag = '0'
            str_fmt = 'Fundamental'
        else:
            flag = 'single particle (perturbatively corrected)'
            fflag = '0_perturb'
            str_fmt = 'Perturbed fundamental'
        header_fmt = '          nm                       eV                        cd'
        if print_flag: print("\nCalculating {} CD spectrum".format(flag))

        # Get the rotational strength (u_x*m_x + u_y*m_y + u_z*m_z) for each state and its energy
        if opsys.path.exists(cwd + '/' + rs_files[i]):
            rs_file = np.loadtxt(rs_files[i], dtype=np.float64)
            len_rs = rs_file.shape[0]
            state = rs_file[:,:2]
            #state = np.concatenate((rs_file[:,0].reshape(len_rs,1),rs_file[:,1].reshape(len_rs,1)), axis=1)
               
            rs = rs_file[:,-2]
            ene = rs_file[:,-3]
            
            if print_flag: print('{} gap: {: .4f} eV'.format(str_fmt, ene.min()*H_to_eV))
            
            # We only want to compute spectra from accurate states near band edge. Apprx. 0.4 Hartree above energy min
            band_edge_idx = ene < ene.min()+(0.4/H_to_eV) 
            new_ene = ene[band_edge_idx]
            rs = rs[band_edge_idx]
            state_cd = state[band_edge_idx]
            
            if print_flag: print('max_rs = {: .6g}'.format(rs.max()))
            
            if calc_abs:
                os_file = np.loadtxt(os_files[i], dtype=np.float64, skiprows=1)
                os = os_file[:,2]**2 # The OS file outputs the sqrt of |u|^2. We want |u|^2 for absorption amplitude
                os_band_edge_idx = ene < ene.min()+(1.5/H_to_eV) # We will plot the absorption at higher energies
                tmp_ene = os_file[:,3][os_band_edge_idx]
                abs_amp = 0
                header_fmt = '          nm                       eV                        cd                      abs'
                
                os = os[os_band_edge_idx]
                state_os = state[os_band_edge_idx]	

            if m_abs:
                m_file = np.loadtxt(m_files[i], dtype=np.float64, skiprows=i)
                m = m_file[:,2]**2 # The M file outputs the sqrt of |m|^2. We want |m|^2 for absorption amplitude
                m_band_edge_idx = ene < ene.min()+(1.5/H_to_eV) # We will plot the absorption at higher energies
                tmp_ene = m_file[:,3][m_band_edge_idx]
                abs_amp = 0
                header_fmt = '          nm                       eV                        cd                      abs                    m_abs'
                
                m = m[m_band_edge_idx]
                state_m = state[m_band_edge_idx]	
            # set up the elements to calculate spectra: Bohr frequencies, x_axis grids...
            Bohr_freq = new_ene # Bohr frequencies in au. E = hbar/w -> w = E/hbar = E (hbar = 1)
            freq = np.linspace(0.06, 0.15, 1500) # Bohr_freq.min()-2*lwid; frequency grid for spectra calculation
            eV_grid = H_to_eV * freq # corresponding energy grid for plotting in eV
            nm_grid = freq_to_wvlen * Bohr_to_nm / freq # lambda grid for plotting in nm
            
            # Generate CD spectra
            
            cd = np.zeros(freq.shape[0])
            for j in range(Bohr_freq.shape[0]):
                if lineshape == 'l':
                    cd += absorp(freq, Bohr_freq[j], lwid, rs[j])
                elif lineshape == 'g':
                    cd += absorp_g(freq, Bohr_freq[j], lwid, rs[j])    
            if calc_abs:
                abs_spec = np.zeros(freq.shape[0])
                for j in range(tmp_ene.shape[0]):
                    abs_spec += absorp(freq, tmp_ene[j], lwid, os[j])
                    
            if m_abs:
                m_abs_spec = np.zeros(freq.shape[0])
                for j in range(tmp_ene.shape[0]):
                    m_abs_spec += absorp(freq, tmp_ene[j], lwid, m[j])
                    
            # Convert amplitude units properly
            #cd = np.asarray(cd, dtype=np.float64)
            cd *= prefactor
            if calc_abs and m_abs:
                abs_spec *= prefactor * (3/2) 
                m_abs_spec *= (prefactor)**2
                cd_out = np.concatenate((nm_grid.reshape(1500,1), eV_grid.reshape(1500,1), cd.reshape(1500,1), abs_spec.reshape(1500,1), m_abs_spec.reshape(1500,1)), axis=1)
            elif calc_abs:
                abs_spec *= prefactor * (3/2)
                cd_out = np.concatenate((nm_grid.reshape(1500,1), eV_grid.reshape(1500,1), cd.reshape(1500,1), abs_spec.reshape(1500,1)), axis=1)
            else:
                cd_out = np.concatenate((nm_grid.reshape(1500,1), eV_grid.reshape(1500,1), cd.reshape(1500,1)), axis=1)
            
            rs_out = np.concatenate((state_cd, new_ene.reshape(new_ene.shape[0],1), rs.reshape(rs.shape[0],1)), axis=1)
            np.savetxt("new_rs{}.dat".format(fflag), rs_out, fmt= '%3d %3d % .8f % .8f')
            if calc_abs:
                os_out = np.concatenate((state_os, tmp_ene.reshape(tmp_ene.shape[0],1), os.reshape(os.shape[0],1)), axis=1)
                np.savetxt("new_os{}.dat".format(fflag), os_out, fmt= '%3d %3d % .8f % .8f')
            if m_abs:
                m_out = np.concatenate((state_m, tmp_ene.reshape(tmp_ene.shape[0],1), m.reshape(m.shape[0],1)), axis=1)
                np.savetxt("new_m{}.dat".format(fflag), m_out, fmt= '%3d %3d % .8f % .8f')
               
            np.savetxt(f"cd{fflag}.dat", cd_out, header=header_fmt)
            
            if print_flag: print('CD generation done!')


def spectrum(w, amp, cntr, wid, lshape = 'l'):
    x = (w - cntr)/(wid/2)
    if lshape == 'g':
        return amp*np.exp( - np.log(2) * x**2 )
    if lshape == 'l':
        return ((amp) / (1 + x**2))
        


main()

if fig_flag:
	target_dir = '/Users/dchabeda/Desktop/Rabani/project/chiral/cdse/final_publication/figs/size_dependence/'
	
	shutil.copy('cd0_perturb.dat', target_dir + filename + '_cd0_perturb.dat')
	shutil.copy('new_rs0_perturb.dat',  target_dir + filename + '_rs0_perturb.dat')
	

