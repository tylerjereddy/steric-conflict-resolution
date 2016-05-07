'''This script runs the iterative procedure that follows the workflow: steric assessment - alchembed - steric assessment, etc.'''

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-input_coord_file_path", type=str)
parser.add_argument('-index_list', nargs='+', required=True, type=int) #should be list of alternating start and end indices
args = parser.parse_args()

def run_steric_resolution_loop(input_coord_file = args.input_coord_file_path, index_list = args.index_list):
    print 'input_coord_file:', input_coord_file
    print 'index_list:', index_list
    
if __name__ == '__main__':
    run_steric_resolution_loop()

