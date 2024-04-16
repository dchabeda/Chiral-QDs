# Written by Daniel Chabeda, Nov 9, 2023

import sys
import math
import numpy as np
from numpy import linalg as alg

# Step 1 Read in parameters

atom_label_dict = {'Cd' : 48,
                    'Se': 34,
                    'I': 53,
                    'Cs': 55,
                    'Pb': 82,
                    'P1': 84,
                    'P2': 85,
                    'P3': 86,
                    'PC5': 87,
                    'PC6': 88}

atom_output_label_dict = {48 : 'Cd',
                    34 : 'Se',
                    53 : 'I',
                    55 : 'Cs',
                    82 : 'Pb',
                    84 : 'Po',
                    85 : 'At',
                    86 : 'Rn',
                    87 : 'Fr',
                    88 : 'Ra'}

key_list = list(atom_label_dict.keys())
val_list = list(atom_label_dict.values())
bond_len = 1.0

def main():
    filename = sys.argv[1]
    foutname = sys.argv[2]
    
    if 'theta' in sys.argv and 'rand' not in sys.argv:
        idx = sys.argv.index('theta')
        theta = float(sys.argv[idx + 1])
        orient = 'theta'
    elif 'rand' in sys.argv and 'theta' not in sys.argv:
        idx = sys.argv.index('theta')
        seed = int(sys.argv[idx+1])
        orient = 'rand'
    else:
        print('ERROR: both random and deterministic options have been given.\nusage:python passivate_chiral.py filein fileout [theta angle] or [rand seed]')
        exit(1)
    
    if not orient:
        print('No orientation type has been specified (rand or theta).\nusage:python passivate_chiral.py filein fileout [theta angle] or [rand seed]')
        exit(1)
    
    print("\nRUNNING PROGRAM \"passivate_chiral.py\" ON FILE: {}".format(filename))

    atoms = open_file(filename, 'atomic')
    atoms = np.array(atoms)
    ind = np.lexsort( (atoms[:,3], atoms[:,2], atoms[:,1]) )
    
    sorted_atoms = atoms[ind]
    
    core_atoms = sorted_atoms[sorted_atoms[:,0]<84]
    ligands = sorted_atoms[sorted_atoms[:,0]>=84]
    
    if orient == 'rand':
        wildparam = seed
    elif orient == 'theta':
        wildparam = np.pi*theta/180
    passivated_coords = passivate_chiral(core_atoms,ligands, orient, wildparam)
    
    write_output(foutname, passivated_coords, 'par')
    write_output(foutname, passivated_coords, 'xyz')



def passivate_chiral(core_atoms,ligs, orient, wildparam):
	if orient == 'rand':
		seed = wildparam
	elif orient == 'theta':
		theta = wildparam
    
	half_lig_coords = ligs[ligs[:,1] < 0.0][::-1]
	zero_ligs = ligs[ligs[:,1] == 0.0]
	half_lig_coords_rev = ligs[ligs[:,1] > 0.0]
	
	coord_cd = core_atoms[core_atoms[:,0]==atom_label_dict['Cd']][:,1:]
	coord_p1 = half_lig_coords[half_lig_coords[:,0]==atom_label_dict['P1']][:,1:]
	coord_p1_rev = half_lig_coords_rev[half_lig_coords_rev[:,0]==atom_label_dict['P1']][:,1:]
	
	if orient == 'rand':
		np.random.seed(seed=seed)

	orig_vals = []
	chiral_ligs = []
	theta_list = []
	z = np.array([0,0,1])
	for i in range(coord_cd.shape[0]):
		ligand_now = []
		cd_now = coord_cd[i]

		for j in range(coord_p1.shape[0]):
			
			dist = np.sqrt(np.sum((cd_now-coord_p1[j])**2))
        	
			if dist<1.5:
				ligand_now.append(coord_p1[j])
				
		for j in range(len(ligand_now)):
			p1_now = ligand_now[j]
			# First get the vector from p1 to cd
			diff = cd_now-p1_now
			a = diff/alg.norm(diff)
			# n = numpy.cross(z, a)
			# if alg.norm(n)>-1e-10:
				# print('check 1 dot prod',np.dot(n,a))
			n = np.cross(z, cd_now)
			min_coords = coord_p1 - p1_now
			
			length_array = np.sort(np.multiply(min_coords,min_coords).sum(1))
			# checks to remove other ligands that are too close
			if alg.norm(n)>1e-6 and length_array[1] > 4:
				
				n -= np.dot(n,a)*a
				# check here
				if orient == 'rand':
				    theta= 2 * np.pi* np.random.rand()
				elif orient == 'theta':
				    theta = theta
				
				theta_list.append(theta)
				
				orig_vals.append((list(p1_now), theta))
				
				# this rotates random norm vector v by amount theta (or at least is supposed to)
				n = np.cos(theta)*n+np.sin(theta)*np.cross(a,n) + (1-np.cos(theta))*(np.dot(a,n))*a
				n = n/alg.norm(n)
				nxa = np.cross(n,a)
				b = -a/3.0+nxa*2.0*np.sqrt(2.0)/3.0
				r = -(a+b)*0.5
				c = r+n*np.sqrt(6)/3.0
				d = r-n*np.sqrt(6)/3.0
				atom1 = np.array([atom_label_dict['PC6']] + list(bond_len*b+p1_now))
				atom2 = np.array([atom_label_dict['PC5']] + list(bond_len*c+p1_now))
				atom3 = np.array([atom_label_dict['P3']] + list(bond_len*d+p1_now))
				
				chiral_ligs.append(atom1); chiral_ligs.append(atom2); chiral_ligs.append(atom3)
	
	for i in range(coord_cd.shape[0]):
		ligand_now = []
		cd_now = coord_cd[i]	
		
		for j in range(coord_p1_rev.shape[0]):
			dist = np.sqrt(np.sum((cd_now-coord_p1_rev[j])**2))
			if dist<1.5:
				ligand_now.append(coord_p1_rev[j])
		
		for j in range(len(ligand_now)):
			p1_now = ligand_now[j]
			# First get the vector from p1 to cd
			diff = cd_now-p1_now
			a = diff/alg.norm(diff)
			# n = numpy.cross(z, a)
			# if alg.norm(n)>-1e-10:
				# print('check 1 dot prod',np.dot(n,a))
			n = np.cross(z, cd_now)
			min_coords = coord_p1_rev - p1_now
			
			length_array = np.sort(np.multiply(min_coords,min_coords).sum(1))
			# checks to remove other ligands that are too close
			
			if alg.norm(n)>1e-6 and length_array[1] > 4:
				#print('\nADDING CHIRAL LIGS\n')
				n -= np.dot(n,a)*a
				# check here
				for orig in orig_vals:
					if orig[0][0] == -p1_now[0] and orig[0][1] == p1_now[1] and orig[0][2] == p1_now[2]:
						
						theta = -orig[1]
				
				# this rotates random norm vector v by amount theta (or at least is supposed to)
				n = np.cos(theta)*n+np.sin(theta)*np.cross(a,n) + (1-np.cos(theta))*(np.dot(a,n))*a
				n = n/alg.norm(n)
				nxa = np.cross(n,a)
				b = -a/3.0+nxa*2.0*np.sqrt(2.0)/3.0
				r = -(a+b)*0.5
				c = r+n*np.sqrt(6)/3.0
				d = r-n*np.sqrt(6)/3.0
				atom1 = np.array([atom_label_dict['PC6']] + list(bond_len*b+p1_now))
				atom2 = np.array([atom_label_dict['PC5']] + list(bond_len*c+p1_now))
				atom3 = np.array([atom_label_dict['P3']] + list(bond_len*d+p1_now))
				
				chiral_ligs.append(atom1); chiral_ligs.append(atom2); chiral_ligs.append(atom3)
	
	out = np.vstack((core_atoms, half_lig_coords, zero_ligs, half_lig_coords_rev, np.array(chiral_ligs)))

	return out

def open_file(file, units):
	if units == 'angstrom':
		unit_conv = 1.0
	if units == 'atomic':
		unit_conv = 0.529177
		
	with open(file, 'r') as f:
		conf = f.readlines()
		f.close()
	
	n_atoms = conf[0]
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


main()

