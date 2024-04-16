import numpy as np
import sys
import os

N = int(sys.argv[1])
nstates = 5

dir_stem = os.getcwd()
mag_files = []
dipole_files = []
angle_labels = []

for i in range(N):
    theta = str(i)
    angle_labels.append(int(theta))
    print(theta, '\n')
    os.chdir(dir_stem +'/' + theta + '/bse')
    
    try:
        dipole_files.append(np.loadtxt('perturb_OS0.dat', skiprows=1, dtype=np.float64, max_rows=20))
    except ValueError:
        os.system(r"sed -i '' -e '$ d' perturb_OS0.dat")
        dipole_files.append(np.loadtxt('perturb_OS0.dat', skiprows=1, dtype=np.float64, max_rows=20))
    try:
        mag_files.append(np.loadtxt('perturb_m0.dat', skiprows=1, dtype=np.float64, max_rows = 20))
    except ValueError:
        os.system(r"sed -i '' -e '$ d' perturb_m0.dat")
        mag_files.append(np.loadtxt('perturb_m0.dat', skiprows=1, dtype=np.float64, max_rows=20))

angle_labels = np.array(angle_labels).reshape(-1,1)
mags = np.array(mag_files)
dipoles = np.array(dipole_files)

# Get the state labels
states = ['   homo-'+str(i)+'    ' for i in range(nstates)]

# Get the dipole vectors of the lowest nstates 

m_vec = mags[:,:nstates,-3:]
d_vec = dipoles[:,:nstates,-3:]

# Take the dot product of the states

dotprod = np.sum(np.multiply(d_vec, m_vec), axis=2)

# u . v = |u| |v| cos theta
d_mag = np.sqrt(np.sum(d_vec**2,axis=2))
m_mag = np.sqrt(np.sum(m_vec**2,axis=2))

# Get theta
# theta = arccos( u . v / |u| |v| )

thetas = np.arccos(dotprod/(d_mag*m_mag))
thetas_deg = 180 * thetas / np.pi # Convert to degrees

out = np.concatenate( (angle_labels, thetas_deg), axis=1)

os.chdir(dir_stem)
np.savetxt("orientation_dep_dipole_magnitude.dat", np.concatenate((angle_labels, d_mag),axis=1), fmt=['% 4d', '%.6e', '%.6e', '%.6e', '%.6e', '%.6e'], header = 'ang' + ''.join(states))
np.savetxt("orientation_dep_mag_dipole_magnitude.dat",  np.concatenate((angle_labels, m_mag),axis=1), fmt=['% 4d', '%.6e', '%.6e', '%.6e', '%.6e', '%.6e'], header = 'ang' + ''.join(states))
np.savetxt(f"u.v_angles_{N}.dat", out, fmt=['% 4d', '%.6e', '%.6e', '%.6e', '%.6e', '%.6e'], header = 'ang' + ''.join(states))



