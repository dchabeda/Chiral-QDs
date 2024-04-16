This directory contains the code used to generate chiral ligand configurations.
The script "passivate_chiral.py" takes a NC geometry without chiral ligands
(but MUST contain passivation ligands P1, P2) and decorates the P1 passivation
ligands on Cd with the chiral ligand potentials (CLPs). The CLPs are a tetrahedral
arrangement of Gaussians around the P1 center.

The script runs as
python passivate_chiral.py [infile_name.par] [outfile_name] [theta or rand] [angle or seed]

The input file path must be specified exactly, including the file extension.
The output file name should be specified without an extension, as several files
of different formats will be generated using the output name. The next keyword
is either "theta" or "rand". 

theta - flag to orient all chiral ligands with the same angle
rand - flag to randomly orient each CLP. (structural symetry cancellations are preserved)

The next keyword is a number that sets either the angle or the initial random seed

An example submission line would be

python passivate_chiral.py conf_nochiral.par 2.7nm_chiral theta 60

This submission line would generate a configuration where the ligands are oriented
at 60º and put the output files into 2.7nm_chiral.par, 2.7nm_chiral.xyz, etc

The chiral ligand lobe identities in the .xyz files were changed for visualization 
purposes to Po, Ra, and Rn. This is arbitrary, and you can choose to make the atoms
appear in any convenient form.
