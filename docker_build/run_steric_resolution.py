'''This script runs the iterative procedure that follows the workflow: steric assessment - alchembed - steric assessment, etc.'''

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-input_coord_file_path", type=str)
parser.add_argument('-index_list', nargs='+', required=True, type=int) #should be list of alternating start and end indices
parser.add_argument('-residue_names_list', nargs='+', required=True, type=str) #should be half as long as index_list
parser.add_argument("-cutoff", type=float, help="cutoff (A)")
parser.add_argument('-list_particles_per_residue', nargs='+', required=True, type=int) #should be half as long as index_list
args = parser.parse_args()

def run_steric_resolution_loop(input_coord_file = args.input_coord_file_path, index_list = args.index_list, residue_names_list = args.residue_names_list, cutoff = args.cutoff, list_particles_per_residue = args.list_particles_per_residue):
    print 'input_coord_file:', input_coord_file
    print 'index_list:', index_list
    print 'residue_names_list:', residue_names_list
    if not len(residue_names_list) == int(len(index_list) / 2.):
        sys.exit('The residue_names_list should be half as long as the index_list as the latter contains start & end indices for each residue.')
    print 'cutoff:', cutoff
    if not len(list_particles_per_residue) == int(len(index_list) / 2.):
        sys.exit('The list_particles_per_residue should be half as long as the index_list as the latter contains start & end indices for each residue.')
    print 'list_particles_per_residue:', list_particles_per_residue
    
if __name__ == '__main__':
    run_steric_resolution_loop()

