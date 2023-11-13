# Written by Daniel Chabeda, Jul 12, 2023

import sys
import numpy as np
from numpy import linalg as alg


# Create dictionaries that map between atomic symbol and number
atom_label_dict = {'Cd' : 48, 'Se': 34, 'I': 53, 'Cs': 55, 'Pb': 82, 'P1': 84, 'P2': 85, 'P3': 86, 'PC5': 87, 'PC6': 88}
atom_output_label_dict = {48 : 'Cd', 34 : 'Se', 53 : 'I', 55 : 'Cs', 82 : 'Pb', 84 : 'Po', 85 : 'At', 86 : 'Rn', 87 : 'Fr', 88 : 'Ra'}

key_list = list(atom_label_dict.keys())
val_list = list(atom_label_dict.values())

# Set the bond length of each chiral ligand lobe from the central P1 ligand.
# Having all ligands the same distance away enables more accurate cancellation of structural chirality in the model
bond_len = 0.5

def main():
	'''Main driver for chiral passivation code. Usage: python passivate_chiral.py [infile.par] [outfilename] [theta or random] [theta or seed].
	This code can either orient all ligands to angle theta (in degrees) or randomly assign each chiral ligand a random angle on uniform [0,2pi). The 
	infile.par script should contain P1 and P2 passivation, but no chiral ligands. The desired outfile name should be written without a file extension.'''

	# Initialize passivation procedure
	# Take in the unpassivated file as sys.argv[1] and the desired name of the output file as sys.argv[2]
	filename = sys.argv[1]
	foutname = sys.argv[2]

	# Determine if the passivation is for an ordered or random configuration
	if 'theta' in sys.argv and 'rand' not in sys.argv:
		idx = sys.argv.index('theta')
		theta = float(sys.argv[idx + 1])
		orient = 'theta'
	elif 'rand' in sys.argv and 'theta' not in sys.argv:
		idx = sys.argv.index('theta')
		seed = int(sys.argv[idx+1])
		orient = 'rand'
	else:
		print('ERROR: both random and deterministic options have been given.\nusage:python passivate_chiral.py infile.par outfilename [theta angle] or [rand seed]')
		exit(1)

	if not orient:
		print('No orientation type has been specified (rand or theta).\nusage:python passivate_chiral.py infile.par outfilename [theta angle] or [rand seed]')
		exit(1)

	if orient == 'rand':
		wildparam = seed
	elif orient == 'theta':
		wildparam = np.pi*theta/180

	print("\nRUNNING PROGRAM \"passivate_chiral.py\" ON FILE: {}".format(filename))

	# Read in the atoms in filename
	atoms = open_file(filename, 'atomic')
	atoms = np.array(atoms)

	'''Sort the atoms first by z, then y, and then x values. When split in half
	along x = 0, this creates two arrays, ligs_pos and ligs_neg, which have sigma_2 symmetric atoms at 
	ligs_pos[i] and ligs_neg[i]. We can easily give lig_pos[i] and ligs_neg[i] angles theta and minus theta 
	respectively to deal with structural chirality.'''
	ind = np.lexsort( ( abs(atoms[:,3]), abs(atoms[:,2]), atoms[:,1] ) )
	sorted_atoms = atoms[ind]

	# Separate the crystal atoms from the P1 and P2 ligands. They are handled differently in the code
	core_atoms = sorted_atoms[sorted_atoms[:,0]<84]
	ligands = sorted_atoms[sorted_atoms[:,0]>=84]

	# Passivate the nanocrystal with chiral ligands
	passivated_coords = passivate_chiral(core_atoms,ligands, orient, wildparam)

	write_output(foutname, passivated_coords, 'par')
	write_output(foutname, passivated_coords, 'xyz')



def passivate_chiral(core_atoms,ligs, orient, wildparam):
	if orient == 'rand':
		seed = wildparam
		np.random.seed(seed=seed)
	elif orient == 'theta':
		theta = wildparam
		
	# We only passivate the ligand potential on Cd, "P1," 
	# so we select just the cadmium and the P1 ligands.
	coord_p1 = ligs[ligs[:,0]==atom_label_dict['P1']][:,1:]
	coord_cd = core_atoms[core_atoms[:,0]==atom_label_dict['Cd']][:,1:]
	
	p1_pos = coord_p1[coord_p1[:,0] > 0.0]; p1_neg = coord_p1[coord_p1[:,0] < 0.0] # P1 ligands with x > 0 and x < 0, respectively
	cd_pos = coord_cd[coord_cd[:,0] > 0.0]; cd_neg = coord_cd[coord_cd[:,0] < 0.0] # Cd atoms with x < 0 and x < 0, respectively.
	# Sort just the negative ligands so that indexing along y and z lines up with the positive ligands
	# Can check by printing p1_pos and p1_neg that p1_neg = sigma_2*p1_pos where sigma_2 reflects across the yz plane.
	# The structural chirality of atoms on x = 0 cannot be canceled. They will not be passivated.
	lig_ind = np.lexsort((abs(p1_neg[:,2]),abs(p1_neg[:,1]), abs(p1_neg[:,0])))
	core_ind = np.lexsort((abs(cd_neg[:,2]),abs(cd_neg[:,1]), abs(cd_neg[:,0])))
	p1_neg = p1_neg[lig_ind]
	cd_neg = cd_neg[core_ind]
	
	# Container for CLP output
	chiral_ligs = []
	print(theta)
	# Passivation loop
	for i in range(cd_pos.shape[0]):
		cd_now = cd_pos[i]

		for j in range(p1_pos.shape[0]):
			dist = np.sqrt(np.sum((cd_now-p1_pos[j])**2))
			if dist<1.5: # Grab P1 ligands within 1.5 Bohr of the current Cd atom
				z = np.array([0,0,1]) # The "north pole" of the quantum dot.
				n = np.cross(z, cd_now) # Used to check if the Cd atom is right on the pole

				p1_now = p1_pos[j] # Current P1 ligand to be passivated
				
				# Code used to check whether other P1 ligands will be too close for chiral ligand passivation.
				lig_distance = p1_pos - p1_now
				distance_arr = np.sort(np.multiply(lig_distance,lig_distance).sum(1))

				if alg.norm(n)>1e-6 and distance_arr[1] > 4:
					if orient == 'rand':
						theta = 2*np.pi*np.random.rand()
					
					clp1, clp2, clp3 = add_chiral_lig(cd_now, p1_now, theta)
					chiral_ligs.append(clp1); chiral_ligs.append(clp2); chiral_ligs.append(clp3)
				
				# Handle the ligands with x < 0
				n = np.cross(z, cd_neg[i]) # Used to check if the Cd atom is right on the pole
				p1_now = p1_neg[j]
				lig_distance = p1_neg - p1_now
				distance_arr = np.sort(np.multiply(lig_distance,lig_distance).sum(1))

				if alg.norm(n)>1e-6 and distance_arr[1] > 4:
					clp1, clp2, clp3 = add_chiral_lig(cd_neg[i], p1_now, -theta)
					chiral_ligs.append(clp1); chiral_ligs.append(clp2); chiral_ligs.append(clp3)

	out = np.vstack((core_atoms, ligs, np.array(chiral_ligs)))

	return out

def open_file(file, units):
	if units == 'angstrom':
		unit_conv = 1.0
	if units == 'atomic':
		unit_conv = 0.529177
		
	with open(file, 'r') as f:
		conf = f.readlines()
		f.close()
	
	conf = conf[1:]

	labels = []
	coords = []
	for i, line in enumerate(conf):
		if len(line) > 0:
			line = line.split()
			
			labels.append(atom_label_dict[line[0]])
			coords.append([float(line[1])*unit_conv, float(line[2])*unit_conv, float(line[3])*unit_conv])
        

	atoms = []
	for i in range(len(labels)):
		atoms.append([labels[i]]+coords[i])
	
	return atoms

def write_output(filename, coords, opt='xyz'):
	if opt == 'xyz':
		units = 1.0
		filename += '.xyz'
		fmt = '\n\n'
	if opt == 'par':
		units = 1/0.529177 
		filename += '.par'
		fmt = '\n'
	if opt == 'H':
		units = 1.0
		filename += '_H.xyz'
		fmt = '\n\n'
	
	file = open(filename, 'w')
	
	n_atoms = str(len(coords))
	file.write(n_atoms+fmt)
	
	for coord in coords:
		if opt == 'H':
			if coord[0] == 'P1' or coord[0] == 'P2':
				coord[0] = 'H'
		if opt == 'xyz':
			symb = atom_output_label_dict[coord[0]]
		if opt == 'par':
			symb = key_list[val_list.index(coord[0])]
		file.write('{name} {x} {y} {z}\n'.format(name=symb,x=coord[1]*units,y=coord[2]*units,z=coord[3]*units))
		
	file.close()

def add_chiral_lig(cd_coord, p1_now, theta):
	print(theta)
	z = np.array([0,0,1]) # The "north pole" of the quantum dot.

	# First get the vector from p1 to cd
	p1_vector = cd_coord-p1_now; p1_vector /= alg.norm(p1_vector)
	n = np.cross(z, cd_coord)
	n -= np.dot(n,p1_vector)*p1_vector # remove the component of p1_vector along the n direction
	
	# this rotates random norm vector v by amount theta (or at least is supposed to)
	n = np.cos(theta)*n+np.sin(theta)*np.cross(p1_vector,n) + (1-np.cos(theta))*(np.dot(p1_vector,n))*p1_vector
	n /= alg.norm(n)
	nxp1 = np.cross(n,p1_vector)
	b = -p1_vector/3.0+nxp1*2.0*np.sqrt(2.0)/3.0
	r = -(p1_vector+b)*0.5
	c = r+n*np.sqrt(6)/3.0
	d = r-n*np.sqrt(6)/3.0
	atom1 = np.array([atom_label_dict['PC6']] + list(bond_len*b+p1_now))
	atom2 = np.array([atom_label_dict['PC5']] + list(bond_len*c+p1_now))
	atom3 = np.array([atom_label_dict['P3']] + list(bond_len*d+p1_now))
	
	return atom1, atom2, atom3


main()
