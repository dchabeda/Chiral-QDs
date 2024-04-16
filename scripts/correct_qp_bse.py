import numpy as np

'''This script is for adding energy corrections to charge excitations using diagonal elements
of the Coulomb two-electron integral matrix. The transition energy from quasiparticle state a
to i (a -> i), is \varepsilon_{ia}. This energy is corrected by K^d + 2K_x, the direct and exchange 
Coulomb integrals. The energy of the resulting state is E_{ia} = \varepsilon_{ia} + K^{d}_{ia,ia} + 2K^{x}_{ia,ia}.
The quasiparticle energies are swapped out for the perturbatively corrected energies, E_{ia}, and
the output files accordingly adjusted.'''

# Read in the original quasiparticle output files

os0 = np.loadtxt("OS0.dat", skiprows=1)
rs0 = np.loadtxt("rs0.dat")
m0 = np.loadtxt("M0.dat")

# Resort the entries by electron states so they align with the BSE matrix

bse_align_ind_os = np.lexsort( (os0[:,0], os0[:,1]) )
bse_align_ind_m = np.lexsort( ( m0[:,0], m0[:,1]) )
os0 = os0[bse_align_ind_os]
rs0 = rs0[bse_align_ind_os]
m0 = m0[bse_align_ind_m]

# Read in the HBSmat.dat file containing the corrected energies E_{ia}

hbsmat = np.loadtxt("HBSmat.dat")

# Extract the diagonal elements of this matrix

corrected_energies = np.diagonal(hbsmat)

# Replace the energies of the quasiparticle output with the corrected energies

os0[:,3] = corrected_energies
rs0[:,2] = corrected_energies
m0[:,3] = corrected_energies

# Sort the quasiparticle output by energy, low to high, and print

ene_lowhigh_ind = np.argsort(os0[:,3])

os0 = os0[ene_lowhigh_ind]
rs0 = rs0[ene_lowhigh_ind]
m0 = m0[ene_lowhigh_ind]

np.savetxt("perturb_OS0.dat", os0, fmt=["% 3d", "% 3d", "% .6f", "% .10f", "% .6f", "% .8f", "% .8f", "% .8f"],\
header="i  a   sqrt|u|^2     Ea-Eb    4/3*E*|u|^2    ux          uy          uz")
np.savetxt("perturb_rs0.dat", rs0, fmt=["% 3d", "% 3d", "% .10f", "% .10f"])
np.savetxt("perturb_m0.dat", m0, fmt=["% 3d", "% 3d", "% .6f", "% .10f", "% .6f", "% .8f", "% .8f", "% .8f"],\
header="i  a   sqrt|m|^2     Ea-Eb    16/3*E*|m|^2    mx          my          mz")

#print("\nDone correcting quasiparticle energies!\n")
# Done!