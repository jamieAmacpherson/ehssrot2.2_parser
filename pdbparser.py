#! /usr/bin/python
#
## Parsing script to convert PDB files to input files for EHSSROT2.2
## James Macpherson, Francis Crick Institute London
##
#
#_____________________________________________________________________
## Imports
import numpy as np
import MDAnalysis as mda
import argparse
import os.path
import sys
import os
import glob
import time
from joblib import Parallel, delayed
import multiprocessing as mp


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
def pdb2ehssrot(pdbfile):
	# hardcode parameters
	traj_per_run = int(10000)
	nruns = int(30)
	# load the pdb file as a universe object 
	u = mda.Universe(pdbfile)
	# convert the atom positions to a numpy array
	all_atoms = u.select_atoms('protein and backbone')
	xyz = np.asarray(all_atoms.positions)
	# convert masses to a numpy array
	masses = all_atoms.masses
	masses = np.reshape(masses, (np.shape(xyz)[0], 1))
	outarray = np.concatenate((xyz, masses), axis=1)
	# add prefixes to EHSSROT file
	natoms = np.shape(xyz)[0]
	prefixes = [int(traj_per_run), int(nruns), int(np.shape(xyz)[0])]
	# save numpy array to file
	np.savetxt('{0}_tmp_ehssrotin.txt'.format(pdbfile), outarray, delimiter='\t')
	# add prefixes to parsed file
	np.savetxt('{0}_prefixes.txt'.format(pdbfile), prefixes, delimiter='\t', fmt='%i')
	os.system('cat {0}_prefixes.txt {1}_tmp_ehssrotin.txt > {2}_ehssrotin.txt'.format(pdbfile, pdbfile, pdbfile))
#
	# hardcode the src code of ehssrot2.2
	src_dir = '/home/macphej/jm.software/apps/ehssrot2.20/EHSSrot.f'
	os.system('mkdir {0}_calcdir'.format(pdbfile))
	os.system('mv {0}_ehssrotin.txt {1}_calcdir/EHSSin.txt'.format(pdbfile, pdbfile))
	os.chdir('{0}_calcdir'.format(pdbfile))
	os.system('gfortran {0} -o ehssrot'.format(src_dir))
	os.system('./ehssrot')
	os.chdir('../')
	os.system('rm {0}_prefixes.txt'.format(pdbfile))
	os.system('rm {0}_tmp_ehssrotin.txt'.format(pdbfile))


#_____________________________________________________________________
## Launch EHSSROT2.2 job for all pdb files in the directory
#
# what are your inputs, and what operation do you want to 
# perform on each input. For example...

if __name__=='__main__':
    # what are your inputs, and what operation do you want to
    # perform on each input. For example...
    path = '*.pdb'
    file_list = glob.glob(path)
    file_list.sort(key=lambda x: os.path.getmtime(x))
#
    inputs = file_list
    #  removing processes argument makes the code run on all available cores
    pool_size = int(mp.cpu_count() - 2)

    if pool_size == 0:
    	pool_size = int(1)


    pool = mp.Pool(processes=pool_size)
    results = pool.map(pdb2ehssrot, file_list)
    
    pool.close()
    pool.join()


