'''This script runs the iterative procedure that follows the workflow: steric assessment - alchembed - steric assessment, etc.'''

import argparse
import time
import sys
import subprocess
import cPickle as pickle
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-input_coord_file_path", type=str)
parser.add_argument('-index_list', nargs='+', required=True, type=int) #should be list of alternating start and end indices
parser.add_argument('-residue_names_list', nargs='+', required=True, type=str) #should be half as long as index_list
parser.add_argument("-cutoff", type=float, help="cutoff (A)")
parser.add_argument('-list_particles_per_residue', nargs='+', required=True, type=int) #should be half as long as index_list
parser.add_argument("-output_path", type=str) #absolute path for data output (i.e., for writing pickle / plot files)
args = parser.parse_args()

def run_steric_resolution_loop(input_coord_file = args.input_coord_file_path, index_list = args.index_list, residue_names_list = args.residue_names_list, cutoff = args.cutoff, list_particles_per_residue = args.list_particles_per_residue, output_path = args.output_path):
    if not len(residue_names_list) == int(len(index_list) / 2.):
        sys.exit('The residue_names_list should be half as long as the index_list as the latter contains start & end indices for each residue.')
    if not len(list_particles_per_residue) == int(len(index_list) / 2.):
        sys.exit('The list_particles_per_residue should be half as long as the index_list as the latter contains start & end indices for each residue.')

    round_number = 1

    def steric_assessment_all_species(round_number):
        '''Run the steric assessment for all specified residues / species. Should return the pickled aggregate numpy array data for steric violations across all residues.
        round_number should be an integeter representing the current round number in the steric resolution loop'''
        current_residue_species_index = 0
        for residue_name, particles_current_residue in zip(residue_names_list, list_particles_per_residue):
            start_index = index_list[current_residue_species_index]
            end_index = index_list[current_residue_species_index + 1]
            
            subprocess.call("python /steric_analysis/run_assessment.py -start_index {start_index} -end_index {end_index} -coord_filepath {input_coord_file} -particles_per_residue {particles_current_residue} -cutoff {cutoff} -pickle_filename {output_path}/steric_viols.p -plot_filename {output_path}/steric_histogram.png".format(start_index=start_index, end_index=end_index, input_coord_file=input_coord_file, particles_current_residue=particles_current_residue, cutoff=cutoff, output_path=output_path))
            if current_residue_species_index == 0: #initialize cumulative array after first residue type analyzed
                cumulative_array_per_residue_steric_conflicts = pickle.load(open('{output_path}/steric_viols.p'.format(output_path=output_path), 'rb'))
            else: #concatenate new residue steric array data with previous residue data
                new_array_per_residue_steric_conflicts = pickle.load(open('{output_path}/steric_viols.p'.format(output_path=output_path), 'rb'))
                cumulative_array_per_residue_steric_conflicts = np.concatenate((cumulative_array_per_residue_steric_conflicts, new_array_per_residue_steric_conflicts))
            current_residue_species_index += 2
        #the full steric assessment across all residue types in the current round has completed -- so write the cumulative array of per-residue steric conflict counts to a round number-tagged pickle file and return the array proper
        pickle.dump((cumulative_array_per_residue_steric_conflicts, open('cumulative_array_per_residue_steric_conflicts_round_{round_number}.p'.format(round_number = round_number),'wb'))
        return cumulative_array_per_residue_steric_conflicts
        
    percentage_residues_with_steric_conflicts_previous_round = 100
    consecutive_rounds_without_improvement = 0
    while 1:
        print 'starting steric assessment round {round_number}'.format(round_number = round_number)
        start_time = time.time()
        cumulative_array_per_residue_steric_conflicts = steric_assessment_all_species(round_number)
        num_residues_with_steric_conflicts = np.count_nonzero(cumulative_array_per_residue_steric_conflicts)
        percentage_residues_with_steric_conflicts = (float(num_residues_with_steric_conflicts) / float(cumulative_array_per_residue_steric_conflicts.shape[0])) * 100.
        if percentage_residues_with_steric_conflicts >= percentage_residues_with_steric_conflicts_previous_round:
            consecutive_rounds_without_improvement += 1
        end_time = time.time()
        steric_assessment_seconds = end_time - start_time
        steric_assessment_minutes = steric_assessment_seconds / 60.
        steric_assessment_hours = steric_assessment_minutes / 60.
        print 'completed steric assessment round {round_number} in ', steric_assessment_seconds, ' seconds or ', steric_assessment_minutes, ' minutes or ', steric_assessment_hours, ' hours'

        #store the current percentage as the previous now that the comparison has been done
        percentage_residues_with_steric_conflicts_previous_round = percentage_residues_with_steric_conflicts

        #if all steric conflicts have been resolved, exit the while loop
        if percentage_residues_with_steric_conflicts == 0:
            print 'All steric conflicts for the input residues have been resolved -- exiting loop.'
            break

        #likewise, exit the while loop if two consecutive rounds have failed to improve the steric situation
        if consecutive_rounds_without_improvement == 2:
            print 'Exiting loop because two consecutive steric conflict resolution rounds have failed to improve the situation.'
            break

        #if there may still be something to gain with another round of alchembed, run it based on user-specified parameters



    
if __name__ == '__main__':
    run_steric_resolution_loop()

