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
import MDAnalysis
import collections

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
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/dppc_simple_test/results/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The dppc_simple_copies.gro file should have all steric conflicts resolved by the last round of steric resolution.")

    def test_resolution_steric_conflicts_argparse(self):
        os.mkdir('/steric_conflict_resolution_work/dppc_simple_test_argparse') #probably don't need a proper temporary object for now because only running on travis where it will get wiped anyway
        '''Test resolution of steric conflicts for simple DPPC test case when run_steric_resolution.py is called via the command line (invoking argparse code).'''
        subprocess.call("python /steric_conflict_resolution/run_steric_resolution.py -input_coord_file_path /steric_conflict_resolution_work/test_data/dppc_simple_copies/dppc_simple_copies.gro -index_list 1 24 -residue_names_list DPPC -cutoff 2.0 -list_particles_per_residue 12 -output_path /steric_conflict_resolution_work/dppc_simple_test_argparse -topology_filepath /steric_conflict_resolution_work/test_data/dppc_simple_copies/sys.top", shell=True)
        list_steric_conflicts_by_round = []
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/dppc_simple_test_argparse/results/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The dppc_simple_copies.gro file should have all steric conflicts resolved by the last round of steric resolution.")

class TestMultipleCopiesDPPC(unittest.TestCase):
    '''Test that the docker image can completely resolve steric conflicts in the case of many superposed DPPC molecules.'''

    def test_resolution_steric_conflicts(self):
        os.mkdir('/steric_conflict_resolution_work/dppc_large_test') #probably don't need a proper temporary object for now because only running on travis where it will get wiped anyway
        run_steric_resolution.run_steric_resolution_loop(input_coord_file = '/steric_conflict_resolution_work/test_data/dppc_many_copies/dppc_many_copies.gro', index_list = [1, 240], residue_names_list = ['DPPC'], cutoff = 2.0, list_particles_per_residue = [12], output_path = '/steric_conflict_resolution_work/dppc_large_test', topology_filepath = '/steric_conflict_resolution_work/test_data/dppc_many_copies/sys.top')
        list_steric_conflicts_by_round = []
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/dppc_large_test/results/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The dppc_many_copies.gro file should have all steric conflicts resolved by the last round of steric resolution.")

    def test_resolution_steric_conflicts_argparse(self):
        os.mkdir('/steric_conflict_resolution_work/dppc_large_test_argparse') #probably don't need a proper temporary object for now because only running on travis where it will get wiped anyway
        '''Test resolution of steric conflicts for simple DPPC test case when run_steric_resolution.py is called via the command line (invoking argparse code).'''
        subprocess.call("python /steric_conflict_resolution/run_steric_resolution.py -input_coord_file_path /steric_conflict_resolution_work/test_data/dppc_many_copies/dppc_many_copies.gro -index_list 1 240 -residue_names_list DPPC -cutoff 2.0 -list_particles_per_residue 12 -output_path /steric_conflict_resolution_work/dppc_large_test_argparse -topology_filepath /steric_conflict_resolution_work/test_data/dppc_many_copies/sys.top", shell=True)
        list_steric_conflicts_by_round = []
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/dppc_large_test_argparse/results/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The dppc_many_copies.gro file should have all steric conflicts resolved by the last round of steric resolution.")

class TestLoopArgCheck(unittest.TestCase):
    '''Check that sys.exit() is called appropriately within run_steric_resolution_loop().'''

    def test_residue_names_list(self):
        with self.assertRaises(SystemExit) as context:
            run_steric_resolution.run_steric_resolution_loop(input_coord_file='dummy', index_list = [1, 24,48,60], residue_names_list = ['DPPC'], cutoff = 2.0, list_particles_per_residue = [12, 16], output_path = 'dummy', topology_filepath = 'dummy', alchembed_b_value = 2, alchembed_resolution = 'CG', alchembed_steps = 1000, alchembed_alpha = 0.1, alchembed_dt = 0.01)
        self.assertEqual('The residue_names_list should be half as long as the index_list as the latter contains start & end indices for each residue.' , context.exception.message)

    def test_list_particles_per_residue(self):
        with self.assertRaises(SystemExit) as context:
            run_steric_resolution.run_steric_resolution_loop(input_coord_file='dummy', index_list = [1, 24,48,60], residue_names_list = ['DPPC', 'POPS'], cutoff = 2.0, list_particles_per_residue = [12], output_path = 'dummy', topology_filepath = 'dummy', alchembed_b_value = 2, alchembed_resolution = 'CG', alchembed_steps = 1000, alchembed_alpha = 0.1, alchembed_dt = 0.01)
        self.assertEqual('The list_particles_per_residue should be half as long as the index_list as the latter contains start & end indices for each residue.', context.exception.message)
            
class TestMultipleResiduesCG(unittest.TestCase):
    '''Test that the docker image can completely resolve steric conflicts in the case of several superposed MARTINI lipid residues of different types.'''

    def test_resolution_steric_conflicts(self):
        os.mkdir('/steric_conflict_resolution_work/multi_residue_cg_test') #probably don't need a proper temporary object for now because only running on travis where it will get wiped anyway
        run_steric_resolution.run_steric_resolution_loop(input_coord_file = '/steric_conflict_resolution_work/test_data/multiple_lipids_MARTINI/multiple_MARTINI_residues_steric_conflicts.gro', index_list = [1, 60, 61, 130, 131, 250, 251, 325], residue_names_list = ['DOPC', 'PRPC', 'XNG3', 'XNSM'], cutoff = 2.0, list_particles_per_residue = [12, 14, 24, 15], output_path = '/steric_conflict_resolution_work/multi_residue_cg_test', topology_filepath = '/steric_conflict_resolution_work/test_data/multiple_lipids_MARTINI/sys.top')
        list_steric_conflicts_by_round = []
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/multi_residue_cg_test/results/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The multiple_MARTINI_residues_steric_conflicts.gro file should have all steric conflicts resolved by the last round of steric resolution.")

class TestRedundantCopiesDPPC_tpr_support(unittest.TestCase):
    '''Test that the docker image can completely resolve steric conflicts in the case of two superposed DPPC molecules with the SAME residue numbers (requires .tpr reading support to pass).'''

    def test_resolution_steric_conflicts(self):
        os.mkdir('/steric_conflict_resolution_work/dppc_redundant_test') #probably don't need a proper temporary object for now because only running on travis where it will get wiped anyway
        run_steric_resolution.run_steric_resolution_loop(input_coord_file = '/steric_conflict_resolution_work/test_data/dppc_redundant_copies/dppc_redundant_copies.gro', index_list = [1, 24], residue_names_list = ['DPPC'], cutoff = 2.0, list_particles_per_residue = [12], output_path = '/steric_conflict_resolution_work/dppc_redundant_test', topology_filepath = '/steric_conflict_resolution_work/test_data/dppc_redundant_copies/sys.top', tpr_filepath = '/steric_conflict_resolution_work/test_data/dppc_redundant_copies/topol.tpr')
        list_steric_conflicts_by_round = []
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/dppc_redundant_test/results/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The dppc_redundant_copies.gro file should have all steric conflicts resolved by the last round of steric resolution.")
        

class TestRedundantRestrainedResidueNames(unittest.TestCase):
    '''Testing cases where the RXXX restrained residue names are the same -- which would normally causes issues in the topology unless handled properly by the steric conflict resolution code.'''

    def test_resolution_steric_conflicts(self):
        os.mkdir('/steric_conflict_resolution_work/dopc_popc_naming_redundancy') #probably don't need a proper temporary object for now because only running on travis where it will get wiped anyway
        run_steric_resolution.run_steric_resolution_loop(input_coord_file = '/steric_conflict_resolution_work/test_data/dopc_popc_naming_redundancy/dopc_popc.gro', index_list = [1, 24, 25, 48], residue_names_list = ['DOPC', 'POPC'], cutoff = 2.0, list_particles_per_residue = [12, 12], output_path = '/steric_conflict_resolution_work/dopc_popc_naming_redundancy', topology_filepath = '/steric_conflict_resolution_work/test_data/dopc_popc_naming_redundancy/sys.top')
        list_steric_conflicts_by_round = []
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/dopc_popc_naming_redundancy/results/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The dopc_popc.gro file should have all steric conflicts resolved by the last round of steric resolution.")


class TestRestarts(unittest.TestCase):
    '''Test restarting steric conflict resolution workflows using appropriate pickle files.'''

    def test_resolution_steric_conflicts_from_restart(self):
        os.mkdir('/steric_conflict_resolution_work/dppc_simple_test_restart') #probably don't need a proper temporary object for now because only running on travis where it will get wiped anyway
        run_steric_resolution.run_steric_resolution_loop(input_coord_file = '/steric_conflict_resolution_work/test_data/dppc_simple_copies/dppc_simple_copies.gro', index_list = [1, 24], residue_names_list = ['DPPC'], cutoff = 2.0, list_particles_per_residue = [12], output_path = '/steric_conflict_resolution_work/dppc_simple_test_restart', topology_filepath = '/steric_conflict_resolution_work/test_data/dppc_simple_copies/sys.top', resume_file = '/steric_conflict_resolution_work/test_data/dppc_simple_copies/cumulative_array_per_residue_steric_conflicts_round_1.p')
        list_steric_conflicts_by_round = []
        for steric_conflict_pickle_file in glob.glob('/steric_conflict_resolution_work/dppc_simple_test_restart/results/cumulative_array_per_residue_steric_conflicts_round_*.p'):
            steric_conflict_data_array = pickle.load(open(steric_conflict_pickle_file, 'rb'))
            num_residues_with_steric_conflicts = np.count_nonzero(steric_conflict_data_array)
            list_steric_conflicts_by_round.append(num_residues_with_steric_conflicts)
        self.assertEqual(min(list_steric_conflicts_by_round), 0, "The dppc_simple_copies.gro file should have all steric conflicts resolved by the last round of steric resolution.")

class TestTopologyDataRedundancy(unittest.TestCase):

    def setUp(self):
        self.topology_data_list = []
        self.list_restrained_residue_names = []
        self.output_universe = MDAnalysis.Universe()
        self.u = MDAnalysis.Universe( '/steric_conflict_resolution_work/test_data/dopc_popc_restraint_redundancy/popc_dopc.gro')
        self.residues = self.u.select_atoms('all').residues
        self.ordered_residue_names = [residue.atoms[0].resname for residue in self.residues]
        self.residues_to_restrain = self.residues[numpy.array([0,1,2,3])] #restrain everything in this case
        self.dict_residue_name_mappings = {}
        self.truncated_names_accounted_for = []
        self.list_redundant_residues = []
        for original_residue_name in set(self.ordered_residue_names):
            if original_residue_name[1:] in self.truncated_names_accounted_for:
                restrained_residue_name = 'Z' + original_residue_name[1:]
                self.list_redundant_residues.append(original_residue_name)
            else:
                restrained_residue_name = 'R' + original_residue_name[1:]
            self.dict_residue_name_mappings[restrained_residue_name] = original_residue_name
            self.truncated_names_accounted_for.append(original_residue_name[1:])
        self.dictionary_residues_to_restrain = collections.defaultdict(list)
        for residue_object in self.residues_to_restrain:
            residue_name = residue_object.name
            if residue_name[0] == 'R': #map back to unrestrained version
                residue_name = self.dict_residue_name_mappings[residue_name]
            residue_ag = residue_object.atoms
            self.dictionary_residues_to_restrain[residue_name].append(residue_ag)
        self.dictionary_residues_not_restrained = collections.defaultdict(list) # no unrestrained residues in this case

    def tearDown(self):
        del self.topology_data_list
        del self.list_restrained_residue_names
        del self.output_universe
        del self.u
        del self.ordered_residue_names
        del self.residues_to_restrain
        del self.dict_residue_name_mappings
        del self.truncated_names_accounted_for
        del self.list_redundant_residues
        del self.dictionary_residues_to_restrain
        del self.dictionary_residues_not_restrained

    def test_redundant_names(self):
        '''Test proper handling of topology associations for redundant cases such as POPC and DOPC mapping to ROPC and ZOPC restrained versions. Related to Issue #17'''
        topology_data_list, output_universe = run_steric_resolution.generate_topology_data(topology_data_list = self.topology_data_list, list_restrained_residue_names = self.list_restrained_residue_names, output_universe = self.output_universe, ordered_residue_names = self.ordered_residue_names, dictionary_residues_to_restrain = self.dictionary_residues_to_restrain , dictionary_residues_not_restrained = self.dictionary_residues_not_restrained)
        self.assertEqual(topology_data_list, [('ROPC', 2), ('POPC', 0), ('ZOPC', 2), ('DOPC', 0)])

class TestTopologyDataRedundancyDownstream(unittest.TestCase):

    def setUp(self):
        self.dict_residue_name_mappings_1 = {'RPPC': 'DPPC', 'RHOL':'CHOL', 'ROPS':'POPS', 'RPGS':'DPGS', 'ROPC':'DOPC', 'RPSM':'DPSM', 'ZOPC':'POPC'}
        self.dict_residue_name_mappings_2 = {'RPPC': 'DPPC', 'RHOL':'CHOL', 'ROPS':'POPS', 'RPGS':'DPGS', 'ZOPC':'DOPC', 'RPSM':'DPSM', 'ROPC':'POPC'} # swapping R/Z for DOPC / POPC

    def tearDown(self):
        del self.dict_residue_name_mappings_1
        del self.dict_residue_name_mappings_2

    def test_redundant_names_downstream_dict_1(self):
        '''Test downstream handling of redundant restrained residue names for dictionary 1 -- related to Issue # 17. The dictionary associations for POPC and DOPC should always be respected.'''
        current_resname = 'POPC'
        candidate_restrained_residue_names_used = ['ROPC']
        candidate_restrained_residue_name = run_steric_resolution.generate_candidate_restrained_residue_name(current_residue_name = current_resname, candidate_restrained_residue_names_used = candidate_restrained_residue_names_used, dict_residue_name_mappings = self.dict_residue_name_mappings_1)
        self.assertEqual(self.dict_residue_name_mappings_1[candidate_restrained_residue_name], current_resname)

        candidate_restrained_residue_names_used = ['ZOPC']
        candidate_restrained_residue_name = run_steric_resolution.generate_candidate_restrained_residue_name(current_residue_name = current_resname, candidate_restrained_residue_names_used = candidate_restrained_residue_names_used, dict_residue_name_mappings = self.dict_residue_name_mappings_1)
        self.assertEqual(self.dict_residue_name_mappings_1[candidate_restrained_residue_name], current_resname)

        current_resname = 'DOPC'
        candidate_restrained_residue_names_used = ['ROPC']
        candidate_restrained_residue_name = run_steric_resolution.generate_candidate_restrained_residue_name(current_residue_name = current_resname, candidate_restrained_residue_names_used = candidate_restrained_residue_names_used, dict_residue_name_mappings = self.dict_residue_name_mappings_1)
        self.assertEqual(self.dict_residue_name_mappings_1[candidate_restrained_residue_name], current_resname)

        candidate_restrained_residue_names_used = ['ZOPC']
        candidate_restrained_residue_name = run_steric_resolution.generate_candidate_restrained_residue_name(current_residue_name = current_resname, candidate_restrained_residue_names_used = candidate_restrained_residue_names_used, dict_residue_name_mappings = self.dict_residue_name_mappings_1)
        self.assertEqual(self.dict_residue_name_mappings_1[candidate_restrained_residue_name], current_resname)
       
    def test_redundant_names_downstream_dict_2(self):
        '''Test downstream handling of redundant restrained residue names for dictionary 2 -- related to Issue # 17. The dictionary associations for POPC and DOPC should always be respected.'''
        current_resname = 'POPC'
        candidate_restrained_residue_names_used = ['ROPC']
        candidate_restrained_residue_name = run_steric_resolution.generate_candidate_restrained_residue_name(current_residue_name = current_resname, candidate_restrained_residue_names_used = candidate_restrained_residue_names_used, dict_residue_name_mappings = self.dict_residue_name_mappings_2)
        self.assertEqual(self.dict_residue_name_mappings_2[candidate_restrained_residue_name], current_resname)

        candidate_restrained_residue_names_used = ['ZOPC']
        candidate_restrained_residue_name = run_steric_resolution.generate_candidate_restrained_residue_name(current_residue_name = current_resname, candidate_restrained_residue_names_used = candidate_restrained_residue_names_used, dict_residue_name_mappings = self.dict_residue_name_mappings_2)
        self.assertEqual(self.dict_residue_name_mappings_2[candidate_restrained_residue_name], current_resname)

        current_resname = 'DOPC'
        candidate_restrained_residue_names_used = ['ROPC']
        candidate_restrained_residue_name = run_steric_resolution.generate_candidate_restrained_residue_name(current_residue_name = current_resname, candidate_restrained_residue_names_used = candidate_restrained_residue_names_used, dict_residue_name_mappings = self.dict_residue_name_mappings_2)
        self.assertEqual(self.dict_residue_name_mappings_2[candidate_restrained_residue_name], current_resname)

        candidate_restrained_residue_names_used = ['ZOPC']
        candidate_restrained_residue_name = run_steric_resolution.generate_candidate_restrained_residue_name(current_residue_name = current_resname, candidate_restrained_residue_names_used = candidate_restrained_residue_names_used, dict_residue_name_mappings = self.dict_residue_name_mappings_2)
        self.assertEqual(self.dict_residue_name_mappings_2[candidate_restrained_residue_name], current_resname)
