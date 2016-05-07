'''This script runs the iterative procedure that follows the workflow: steric assessment - alchembed - steric assessment, etc.'''

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-input_coord_file_path", type=str)
args = parser.parse_args()

def run_steric_resolution_loop(input_coord_file = args.input_coord_file_path):
    print 'input_coord_file:', input_coord_file
    
if __name__ = '__main__':
    run_steric_resolution_loop()

