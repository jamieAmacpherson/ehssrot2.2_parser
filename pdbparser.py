#! /usr/bin/python
#
## Parsing script to convert PDB files to input files for EHSSROT2.2
## James Macpherson, Francis Crick Institute London
##
#
#_____________________________________________________________________
## Imports
import numpy
import MDAnalysis as mda
import argparse
import os.path
import sys
import os

#_____________________________________________________________________
## Parse commandline arguments
# check if input file exists
# there are two inputs: trajectory (.dcd) and topology (.pdb)
# if either of those inputs are not supplied, or if the user doesn't invoke the
# help flag the program will display an error message.
def is_valid_file(arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist! Use the --help flag for input options." % arg)
    else:
        return arg

# command line argument parser
parser = argparse.ArgumentParser(description='Parse PDB file for EHSSROT2.2 input')

# the second argument is the topology file (.pdb) supplied after the -s flag
# this is saved an an obect with the variable args.pdbfile
parser.add_argument("-i", dest="pdbfile", required=True,
                    help="Free structure file (format: .pdb)",
                    type=lambda x: is_valid_file(x))

# the arguments are parsed 
args = parser.parse_args()

#_____________________________________________________________________
## Convert PDB file to input for EHSSROT2.2
#
def pdb2ehssrot(pdbfile, traj_per_run, nruns):
	# load the pdb file as a universe object 
	u = mda.Universe(pdbfile)
	# convert the atom positions to a numpy array
	atoms = ref.select_atoms('protein and backbone')
	xyz = np.asarray(all_atoms.positions)
	# convert masses to a numpy array
	masses = all_atoms.masses
	masses = np.reshape(masses, (np.shape(xyz)[0], 1))
	outarray = np.concatenate((xyz, masses), axis=1)
	# add prefixes to EHSSROT file
	natoms = np.shape(xyz)[0]
	prefixes = [traj_per_run, nruns, np.shape(xyz)[0]]
	# save numpy array to file
	np.savetxt('tmp_ehssrotin.txt', outarray, delimiter='\t')
	# add prefixes to parsed file
	np.savetxt('prefixes.txt', prefixes, delimiter='\t')
	os.system('cat prefixes.txt tmp_ehssrotin.txt > %s_ehssrotin.txt' % pdbfile)
	




