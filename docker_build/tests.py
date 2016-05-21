import os
import glob
import unittest
import numpy as np
import numpy.testing
import cPickle as pickle
import run_steric_resolution
import subprocess
import mock
import sys

class TestArgparse(unittest.TestCase):
    '''Test argparse code.'''

    def test_parser_simple_dppc_input(self):
        with mock.patch('sys.argv', ['run_steric_resolution.py', '-input_coord_file_path', '/steric_conflict_resolution_work/test_data/dppc_simple_copies/dppc_simple_copies.gro', '-index_list', '1', '24', '-residue_names_list', 'DPPC', '-cutoff', '2.0', '-list_particles_per_residue', '12', '-output_path', '/steric_conflict_resolution_work/dppc_simple_test_argparse', '-topology_filepath', '/steric_conflict_resolution_work/test_data/dppc_simple_copies/sys.top']):
            args = run_steric_resolution.parse_args(sys.argv[1:])
            self.assertEqual(args.input_coord_file_path , '/steric_conflict_resolution_work/test_data/dppc_simple_copies/dppc_simple_copies.gro')
            self.assertEqual(args.index_list , [1, 24])
            self.assertEqual(args.residue_names_list , ['DPPC'])
            self.assertEqual(args.cutoff , 2.0)
            self.assertEqual(args.list_particles_per_residue , [12])
            self.assertEqual(args.output_path , '/steric_conflict_resolution_work/dppc_simple_test_argparse')
            self.assertEqual(args.topology_filepath , '/steric_conflict_resolution_work/test_data/dppc_simple_copies/sys.top')

class TestSimpleCopiesDPPC(unittest.TestCase):
    '''Test that the docker image can completely resolve steric conflicts in the simple case of two superposed DPPC molecules.'''

    def test_resolution_steric_conflicts(self):
        os.mkdir('/steric_conflict_resolution_work/dppc_simple_test') #probably don't need a proper temporary object for now because only running on travis where it will get wiped anyway
        run_steric_resolution.run_steric_resolution_loop(input_coord_file = '/steric_conflict_resolution_work/test_data/dppc_simple_copies/dppc_simple_copies.gro', index_list = [1, 24], residue_names_list = ['DPPC'], cutoff = 2.0, list_particles_per_residue = [12], output_path = '/steric_conflict_resolution_work/dppc_simple_test', topology_filepath = '/steric_conflict_resolution_work/test_data/dppc_simple_copies/sys.top')
        list_steric_conflicts_by_round = []
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/dppc_simple_test/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The dppc_simple_copies.gro file should have all steric conflicts resolved by the last round of steric resolution.")

    def test_resolution_steric_conflicts_argparse(self):
        os.mkdir('/steric_conflict_resolution_work/dppc_simple_test_argparse') #probably don't need a proper temporary object for now because only running on travis where it will get wiped anyway
        '''Test resolution of steric conflicts for simple DPPC test case when run_steric_resolution.py is called via the command line (invoking argparse code).'''
        subprocess.call("python /steric_conflict_resolution/run_steric_resolution.py -input_coord_file_path /steric_conflict_resolution_work/test_data/dppc_simple_copies/dppc_simple_copies.gro -index_list 1 24 -residue_names_list DPPC -cutoff 2.0 -list_particles_per_residue 12 -output_path /steric_conflict_resolution_work/dppc_simple_test_argparse -topology_filepath /steric_conflict_resolution_work/test_data/dppc_simple_copies/sys.top", shell=True)
        list_steric_conflicts_by_round = []
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/dppc_simple_test_argparse/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The dppc_simple_copies.gro file should have all steric conflicts resolved by the last round of steric resolution.")

