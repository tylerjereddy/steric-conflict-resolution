'''This script runs the iterative procedure that follows the workflow: steric assessment - alchembed - steric assessment, etc.'''

import glob
import argparse
import time
import sys
import subprocess
import cPickle as pickle
import numpy as np
import generate_mdp
import collections
import MDAnalysis
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import os
assert int(MDAnalysis.__version__.split('.')[1]) > 14, "MDAnalysis version must be > 0.14."

def generate_candidate_restrained_residue_name(current_residue_name, candidate_restrained_residue_names_used, dict_residue_name_mappings):
    for restrained_residue_name, unrestrained_residue_name in dict_residue_name_mappings.iteritems():
        if current_residue_name == unrestrained_residue_name:
            candidate_restrained_residue_name = restrained_residue_name
    return candidate_restrained_residue_name

def generate_topology_data(topology_data_list, list_restrained_residue_names, output_universe, ordered_residue_names, dictionary_residues_to_restrain, dictionary_residues_not_restrained):
    residue_names_accounted_for = []
    for residue_name in ordered_residue_names:
        if residue_name in residue_names_accounted_for:
            continue
        else: 
            residue_names_accounted_for.append(residue_name)
        if residue_name in dictionary_residues_to_restrain.keys():
            new_restrained_residue_name = 'R' + residue_name[1:] #create a unique residue name for the restrained version of the residue (so that separate .itp may be used, etc.)
            if new_restrained_residue_name in list_restrained_residue_names: #in special cases like DOPC and POPC
                new_restrained_residue_name = 'Z' + new_restrained_residue_name[1:]
            list_restrained_residue_names.append(new_restrained_residue_name)
            list_restrained_residue_atomgroups = dictionary_residues_to_restrain[residue_name]
            for ag in list_restrained_residue_atomgroups:
                ag.set_resnames(new_restrained_residue_name)
            if output_universe.select_atoms('all').n_atoms == 0:
                output_universe = MDAnalysis.Merge(*list_restrained_residue_atomgroups)
            else:
                output_universe = MDAnalysis.Merge(output_universe.atoms, *list_restrained_residue_atomgroups)
            output_universe = MDAnalysis.Merge(output_universe.atoms, *dictionary_residues_not_restrained[residue_name])
            topology_data_list.append((new_restrained_residue_name, len(list_restrained_residue_atomgroups)))
            topology_data_list.append((residue_name, len(dictionary_residues_not_restrained[residue_name])))
        else: #just merge in the unrestrained residues if there are no restrained targets for this residue type
            if output_universe.select_atoms('all').n_atoms == 0:
                output_universe = MDAnalysis.Merge(*dictionary_residues_not_restrained[residue_name])
            else:
                output_universe = MDAnalysis.Merge(output_universe.atoms, *dictionary_residues_not_restrained[residue_name])
            topology_data_list.append((residue_name, len(dictionary_residues_not_restrained[residue_name])))
    return topology_data_list, output_universe

def cleanup_files(output_path):
    for filepath in glob.glob(output_path + '/*'):
        if 'alchembed_round' in filepath and not 'gro' in filepath:
            os.remove(filepath)
        elif filepath.split('/')[-1] in ['steric_viols.p', 'steric_histogram.png']:
            os.remove(filepath)

    for filepath in glob.glob('/'.join(output_path.split('/')[:-1]) + '/*'):
        if 'adjusted' in filepath or 'restrained' in filepath:
            os.remove(filepath)

def plot_cumulative_steric_conflicts(cumulative_conflicts_array, outfile, cutoff):
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    matplotlib.pyplot.xticks(rotation=0)
    unique_contact_counts, num_residues_with_those_counts = np.unique(cumulative_conflicts_array, return_counts=True)
    percent_denominator = cumulative_conflicts_array.size / 100.0
    ax.bar(unique_contact_counts, num_residues_with_those_counts / percent_denominator,facecolor = 'green',alpha=0.75,width=0.7, align='center', log = True) #percent histogram
    ax.set_xlim(-1, unique_contact_counts.max() + 1)
    ax.set_xticks(np.arange(-1,unique_contact_counts.max() + 2))
    ax.set_xlabel('# of contacts within {cutoff} $\AA$'.format(cutoff=cutoff))
    ax.set_ylabel('log(%) of residues')
    fig.set_size_inches(6,6)
    fig.subplots_adjust(bottom=0.2, left=0.2)
    fig.savefig(outfile,dpi=300)

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_coord_file_path", type=str)
    parser.add_argument('-index_list', nargs='+', required=True, type=int) #should be list of alternating start and end indices
    parser.add_argument('-residue_names_list', nargs='+', required=True, type=str) #should be half as long as index_list
    parser.add_argument("-cutoff", type=float, help="cutoff (A)")
    parser.add_argument('-list_particles_per_residue', nargs='+', required=True, type=int) #should be half as long as index_list
    parser.add_argument("-output_path", type=str) #absolute path for data output (i.e., for writing pickle / plot files)
    parser.add_argument("-alchembed_b_value", type=int, default = 2)  # b = 2 appears to be optimal based on the Alchembed paper, so that is the default
    parser.add_argument("-alchembed_resolution", type=str, default = 'CG')  # the simulation 'resolution' (AT vs. CG -- currently only have CG .mdp implemented)
    parser.add_argument("-alchembed_steps", type=int, default = 1000)  # code will use this as a starting point, but if it doesn't succeed it will increase the number of steps
    parser.add_argument("-alchembed_alpha", type=float, default = 0.1)  
    parser.add_argument("-alchembed_dt", type=float, default = 0.01)  
    parser.add_argument("-topology_filepath", type=str) # for GROMACS .top file
    parser.add_argument("-tpr_filepath", type=str, default = None) # for GROMACS .tpr file
    parser.add_argument("-resume_file", type=str, default = None) # for path to pickle file used to resume an interrupted run
    return parser.parse_args(args)

def run_steric_resolution_loop(input_coord_file, index_list, residue_names_list, cutoff, list_particles_per_residue, output_path, topology_filepath, alchembed_b_value = 2, alchembed_resolution = 'CG', alchembed_steps = 1000, alchembed_alpha = 0.1, alchembed_dt = 0.01, tpr_filepath = None, resume_file = None):
    if not len(residue_names_list) == int(len(index_list) / 2.):
        sys.exit('The residue_names_list should be half as long as the index_list as the latter contains start & end indices for each residue.')
    if not len(list_particles_per_residue) == int(len(index_list) / 2.):
        sys.exit('The list_particles_per_residue should be half as long as the index_list as the latter contains start & end indices for each residue.')

    round_number = 1

    output_path = output_path + '/results'
    os.mkdir(output_path)

    def steric_assessment_all_species(round_number):
        '''Run the steric assessment for all specified residues / species. Should return the pickled aggregate numpy array data for steric violations across all residues.
        round_number should be an integeter representing the current round number in the steric resolution loop'''
        current_residue_species_index = 0
        for residue_name, particles_current_residue in zip(residue_names_list, list_particles_per_residue):
            start_index = index_list[current_residue_species_index]
            end_index = index_list[current_residue_species_index + 1]
            
            subprocess.call("python /steric_analysis/run_assessment.py -start_index {start_index} -end_index {end_index} -coord_filepath {input_coord_file} -particles_per_residue {particles_current_residue} -cutoff {cutoff} -pickle_filename {output_path}/steric_viols.p -plot_filename {output_path}/steric_histogram.png".format(start_index=start_index, end_index=end_index, input_coord_file=input_coord_file, particles_current_residue=particles_current_residue, cutoff=cutoff, output_path=output_path), shell=True)
            if current_residue_species_index == 0: #initialize cumulative array after first residue type analyzed
                cumulative_array_per_residue_steric_conflicts = pickle.load(open('{output_path}/steric_viols.p'.format(output_path=output_path), 'rb'))
            else: #concatenate new residue steric array data with previous residue data
                new_array_per_residue_steric_conflicts = pickle.load(open('{output_path}/steric_viols.p'.format(output_path=output_path), 'rb'))
                cumulative_array_per_residue_steric_conflicts = np.concatenate((cumulative_array_per_residue_steric_conflicts, new_array_per_residue_steric_conflicts))
            current_residue_species_index += 2
        #the full steric assessment across all residue types in the current round has completed -- so write the cumulative array of per-residue steric conflict counts to a round number-tagged pickle file and return the array proper
        pickle.dump(cumulative_array_per_residue_steric_conflicts, open('{output_path}/cumulative_array_per_residue_steric_conflicts_round_{round_number}.p'.format(round_number = round_number, output_path=output_path),'wb'))
        plot_cumulative_steric_conflicts(cumulative_array_per_residue_steric_conflicts, '{output_path}/cumulative_per_residue_steric_conflicts_round_{round_number}.png'.format(output_path=output_path, round_number=round_number), cutoff = cutoff)
        return cumulative_array_per_residue_steric_conflicts
        
    percentage_residues_with_steric_conflicts_previous_round = 100
    consecutive_rounds_without_improvement = 0
    resume_counter = 0
    while 1:
        print 'starting steric assessment round {round_number}'.format(round_number = round_number)
        start_time = time.time()
        if resume_file is not None and resume_counter == 0:
            cumulative_array_per_residue_steric_conflicts = pickle.load(open(resume_file, 'rb'))
            resume_counter += 1
        else:
            cumulative_array_per_residue_steric_conflicts = steric_assessment_all_species(round_number)
        num_residues_with_steric_conflicts = np.count_nonzero(cumulative_array_per_residue_steric_conflicts)
        # it may be sensible to apply position restraints to those residues with minimal steric conflicts -- as a first step, identify them
        indices_residues_minimal_steric_conflicts = np.where(cumulative_array_per_residue_steric_conflicts == cumulative_array_per_residue_steric_conflicts.min())[0]
        percentage_residues_with_steric_conflicts = (float(num_residues_with_steric_conflicts) / float(cumulative_array_per_residue_steric_conflicts.shape[0])) * 100.
        if percentage_residues_with_steric_conflicts >= percentage_residues_with_steric_conflicts_previous_round and not round_number==1:
            consecutive_rounds_without_improvement += 1
        else:
            consecutive_rounds_without_improvement = 0
        end_time = time.time()
        steric_assessment_seconds = end_time - start_time
        steric_assessment_minutes = steric_assessment_seconds / 60.
        steric_assessment_hours = steric_assessment_minutes / 60.
        print 'completed steric assessment round {round_number} in '.format(round_number=round_number), steric_assessment_seconds, ' seconds or ', steric_assessment_minutes, ' minutes or ', steric_assessment_hours, ' hours'

        #store the current percentage as the previous now that the comparison has been done
        percentage_residues_with_steric_conflicts_previous_round = percentage_residues_with_steric_conflicts

        #if all steric conflicts have been resolved, exit the while loop
        if percentage_residues_with_steric_conflicts == 0:
            print 'All steric conflicts for the input residues have been resolved -- exiting loop.'
            cleanup_files(output_path)
            break

        #likewise, exit the while loop if two consecutive rounds have failed to improve the steric situation
        if consecutive_rounds_without_improvement == 20:
            print 'Exiting loop because 20 consecutive steric conflict resolution rounds have failed to improve the situation.'
            cleanup_files(output_path)
            break

        #if there may still be something to gain with another round of alchembed, run it based on user-specified parameters
        
        # because .mdp files contain some information related to position restraints (i.e., define line), before writing the .mdp file I will want to re-write the coordinate file with separate contiguous groups for the residues involved in indices_residues_minimal_steric_conflicts; I will also need a custom .top file that reflects the associated unique residue names and corresponding forcefield files (basically, the same residue type will need unique names & forcefield files if it falls in the minimal group where some of its molecules are to be restrained and others are not)

        # assuming that we actually have some 'minimum steric conflict' residues, what is the best way to separate them out from the others and re-write the coordinates/ topology / etc?
        # can probably use MDAnalysis to perform the necessary operations on indexed residue objects
        if tpr_filepath is not None:
            u = MDAnalysis.Universe(tpr_filepath, input_coord_file)
        else:
            u = MDAnalysis.Universe(input_coord_file)
        all_selection = u.select_atoms('all')
        residues = all_selection.residues
        residues_to_restrain = residues[indices_residues_minimal_steric_conflicts] # will this work as intended?
        # now, try to rewrite the input coordinate file with restrained and unrestrained residues separated (but might make sense to have the same residue type follow in a restrained / non-restrained topology ordering)

        # if we're going to write a new universe we'll also need information for the residues that are not to be position restrained (in a format that will allow me to access on a per-residue-type basis, because I'll want to rewrite the coordinates with i.e., POPS-restrained, POPS-unrestrained, DOPE-restrained, DOPE-unrestrained topology configuration)
        total_num_residues = all_selection.n_residues
        array_all_residue_numbers = np.arange(total_num_residues)
        mask_residues_not_restrained = np.in1d(array_all_residue_numbers, indices_residues_minimal_steric_conflicts, invert=True)
        unrestrained_residues = residues[mask_residues_not_restrained]

        ordered_residue_names = [residue.atoms[0].resname for residue in residues]
        if round_number == 1:
            original_residue_names = set(ordered_residue_names)
            dict_residue_name_mappings = {}
            truncated_names_accounted_for = []
            list_redundant_residues = []
            for original_residue_name in original_residue_names:
                if original_residue_name[1:] in truncated_names_accounted_for:
                    restrained_residue_name = 'Z' + original_residue_name[1:]
                    list_redundant_residues.append(original_residue_name)
                else:
                    restrained_residue_name = 'R' + original_residue_name[1:]
                dict_residue_name_mappings[restrained_residue_name] = original_residue_name
                truncated_names_accounted_for.append(original_residue_name[1:])
        else: #reset the residue names
            new_list_residue_names = []
            for residue_name in ordered_residue_names:
                if residue_name in list_redundant_residues:
                    restrained_version = 'Z' + residue_name[1:]
                else:
                    restrained_version = 'R' + residue_name[1:]
                new_list_residue_names.append(dict_residue_name_mappings[restrained_version]) #map back unrestrained residue names
            ordered_residue_names = new_list_residue_names[:]

        # need to access the names of the affected residues in order here?
        dictionary_residues_to_restrain = collections.defaultdict(list)
        for residue_object in residues_to_restrain:
            residue_name = residue_object.name
            if residue_name[0] == 'R': #map back to unrestrained version
                residue_name = dict_residue_name_mappings[residue_name]
            residue_ag = residue_object.atoms
            dictionary_residues_to_restrain[residue_name].append(residue_ag)
        # so, dictionary_residues_to_restrain should have a data structure like this: {'POPS': [POPS_1_ag, POPS_2_ag], ... }

        dictionary_residues_not_restrained = collections.defaultdict(list)
        for residue_object in unrestrained_residues:
            residue_name = residue_object.name
            if residue_name[0] == 'R': #map back to unrestrained version
                residue_name = dict_residue_name_mappings[residue_name]
            residue_ag = residue_object.atoms
            residue_ag.set_resnames(residue_name)
            dictionary_residues_not_restrained[residue_name].append(residue_ag)
        # so, dictionary_residues_not_restrained will have a similar data structure to the restrained version

        # iterate through all the residue types in the system, first writing the restrained and then the unrestrained versions of the residues, all while keeping track of the topological details so that I can write the custom .top file later on as well
        topology_data_list = []
        list_restrained_residue_names = []
        output_universe = MDAnalysis.Universe()
        topology_data_list, output_universe = generate_topology_data(topology_data_list, list_restrained_residue_names, output_universe, ordered_residue_names, dictionary_residues_to_restrain, dictionary_residues_not_restrained)

        # write the new input data to a coord file
        forcefield_parent_filepath = '/'.join(topology_filepath.split('/')[:-1]) + '/' # remove the .top file name to obtain the parent path for the input .itp files
        adjusted_coords = forcefield_parent_filepath + 'adjusted_coords.gro'
        output_universe.atoms.write(adjusted_coords)
        #adjust the box vectors (which get set to 0 for whatever reason)
        subprocess.call(['/bin/bash','-i','-c','gmx editconf -f {coord_path} -o {coord_path} -box 300.0 300.0 300.0'.format(coord_path=adjusted_coords)])
        input_coord_file = adjusted_coords #use the new coord file as the algorithm input

        # up next, need to deal with writing new .top and .itp files in preparation for the selective application of position restraints
        # the new topology file will have to include the new .itp files that enable position restraints so I should generate the new .itp files before trying to generate the new .top file
        # however, to obtain the paths of the current .itp files (which I will presumably want to modify programmatically) I will likely want to use the path of the input topology file
        # generate a list of the full filepaths for each .itp file by parsing the topology_filepath
        list_input_itp_filepaths = []
        with open(topology_filepath, 'r') as input_topology:
            for line in input_topology:
                if '#include' in line:
                    itp_filename = line.split(' ')[1].replace('"','').replace("'",'').rstrip()
                    itp_filepath = forcefield_parent_filepath + itp_filename
                    list_input_itp_filepaths.append(itp_filepath)

        # now, try to generate the new .itp files that are needed for application of position restraints (will probably always want to keep copies of the old .itp files as well as we'll likely often have a number of residues that are to remain unrestrained)
        residue_names_to_restrain = dictionary_residues_to_restrain.keys()
        list_new_restrained_itp_files = []
        candidate_restrained_residue_names_used = []
        for input_itp_filepath in list_input_itp_filepaths:
            if 'restrained' in input_itp_filepath:
                continue
            residue_found = 0
            # it will only be necessary to produce special 'restrained' .itp files for the applicable ('minimum conflict') residues
            with open(input_itp_filepath, 'r') as input_itp_file:
                list_input_itp_file_lines = input_itp_file.readlines()
                index_counter = 0
                for line in list_input_itp_file_lines:
                    if 'molname' in line and list_input_itp_file_lines[index_counter + 1].strip().split()[0] in residue_names_to_restrain:
                        name_line_index = index_counter + 1
                        residue_found += 1
                        break
                    else:
                        index_counter += 1
                if residue_found == 0:
                    continue
                current_residue_name = list_input_itp_file_lines[name_line_index].strip().split(' ')[0]
                candidate_restrained_residue_name = generate_candidate_restrained_residue_name(current_residue_name, candidate_restrained_residue_names_used, dict_residue_name_mappings)

                candidate_restrained_residue_names_used.append(candidate_restrained_residue_name)
                if candidate_restrained_residue_name in residue_names_to_restrain or current_residue_name in residue_names_to_restrain:
                    #need to generate a special 'restrained' .itp file (will modify the list_input_itp_file_lines)
                    current_itp_filename = input_itp_filepath.split('/')[-1]
                    new_restrained_itp_filename = ''.join(current_itp_filename.split('.')[:-1]) + '_restrained.itp'
                    list_new_restrained_itp_files.append(new_restrained_itp_filename)
                    new_restrained_itp_filepath = '/'.join(input_itp_filepath.split('/')[:-1]) + '/' + new_restrained_itp_filename
                    with open(new_restrained_itp_filepath, 'w') as output_itp_file:
                        # I'll want to use the new molname in the new (restrained) .itp
                        old_name_line = list_input_itp_file_lines[name_line_index]
                        new_name_line = old_name_line.replace(current_residue_name, candidate_restrained_residue_name)
                        list_input_itp_file_lines[name_line_index] = new_name_line
                        atom_counter = dictionary_residues_to_restrain[current_residue_name][0].n_atoms

                        # also: place an appropriate [ position_restraints ] section within the [ moleculetype ] section
                        # identify the indices of the lines between which the [ position_restraints ] information will be spliced in
                        index_counter = 0
                        for line in list_input_itp_file_lines:
                            if 'atoms' in line and line[0] != ';':
                                atoms_section_start_index = index_counter
                                break
                            else:
                                index_counter += 1

                        index_counter = atoms_section_start_index
                        for line in list_input_itp_file_lines[atoms_section_start_index:]:
                            if line == "\n":
                                posres_start_index = index_counter
                                break
                            else:
                                index_counter += 1

                        list_new_posres_lines = []
                        list_new_posres_lines.append('[ position_restraints ]\n')
                        list_new_posres_lines.append('; ai  funct  fcx    fcy    fcz\n')

                        def generate_posres_line(atom_index):
                            '''Use the atom_index integer argument to generate a strong position restraint line (force of 1000 in all dimensions).
                            Should return the appropriate string for the line.'''
                            if int(atom_index) <= 9:
                                spaces = ' ' * 3
                            else:
                                spaces = ' ' * 2
                            return spaces + '{atom_index}    1    1000   1000   1000\n'.format(atom_index=atom_index)

                        # need to determine the total number of atoms in the molecule in the input itp file to appropriately generate the position restraint data
                        for atom_index in xrange(1,atom_counter + 1):
                            posres_string = generate_posres_line(atom_index)
                            list_new_posres_lines.append(posres_string)

                        # up next, combine the new posres lines with list_input_itp_file_lines (for eventual writing to the new .itp file)
                        list_new_itp_file_lines = list_input_itp_file_lines[:posres_start_index] + list_new_posres_lines + list_input_itp_file_lines[posres_start_index:]

                        # it may also be sensible to replace all instances of the old residue name in the itp file with the new (restrained) residue name (if only to alleviate confusion)
                        list_new_itp_file_lines_renamed = []
                        for line in list_new_itp_file_lines:
                            list_new_itp_file_lines_renamed.append(line.replace(current_residue_name, candidate_restrained_residue_name))

                        output_itp_file.writelines(list_new_itp_file_lines_renamed)

        # use topology_data_list as part of the process to generate the new .top file (and perhaps consider cannibilizing the old / input .top file?)

        molecules_section = 0
        molecule_section_written = False
        list_itps_accounted_for = []
        with open(topology_filepath, 'r') as input_topology:
            topology_filepath = '/'.join(topology_filepath.split('/')[:-1]) + '/adjusted_topology_round_{round_num}.top'.format(round_num = round_number)
            with open(topology_filepath, 'w') as output_topology:
                for line in input_topology:
                    if '#include' in line:
                        if line not in list_itps_accounted_for:
                            output_topology.write(line)
                            list_itps_accounted_for.append(line)
                    elif '[ system ]' in line:
                        for restrained_itp_filename in list_new_restrained_itp_files:
                            include_string = '#include ' + '"' + restrained_itp_filename + '"\n'
                            if include_string not in list_itps_accounted_for:
                                output_topology.write(include_string)
                                list_itps_accounted_for.append(include_string)
                        output_topology.write('[ system ]\n')
                    elif '[ molecules ]' in line:
                        output_topology.write(line)
                        molecules_section += 1
                    elif molecules_section > 0 and not molecule_section_written:
                        for residue_name, num_residues in topology_data_list:
                            if num_residues > 0:
                                output_topology.write(str(residue_name) + ' ' + str(num_residues) + '\n')
                        molecule_section_written = True
                        break
                    else: 
                        output_topology.write(line)

        #we will need to generate an .mdp file for the alchembed simulation
        alchembed_mdp_filename = '{output_path}/alchembed_round_{round_number}.mdp'.format(output_path=output_path, round_number=round_number)
        if consecutive_rounds_without_improvement == 0:
            actual_alchembed_steps = alchembed_steps
            alchembed_delta_lambda = 1. / float(actual_alchembed_steps)
        elif consecutive_rounds_without_improvement == 1:
            actual_alchembed_steps = alchembed_steps * 10 #slow it down if last attempt was problematic
            alchembed_delta_lambda = 1. / float(actual_alchembed_steps)
        else:
            actual_alchembed_steps = alchembed_steps * 10 #slow it down if last attempt was problematic
            #terminate at a progessively lower lambda value if mdrun failures continue to accumulate
            alchembed_delta_lambda = (1.0/float(consecutive_rounds_without_improvement)) / float(actual_alchembed_steps)

        generate_mdp.generate_mdp(resolution=alchembed_resolution, output_filename = alchembed_mdp_filename, b = alchembed_b_value, steps = actual_alchembed_steps, delta_lambda = alchembed_delta_lambda, alpha = alchembed_alpha, dt = alchembed_dt)

        #generate alchembed tpr file
        tpr_filename = alchembed_mdp_filename.replace('mdp','tpr')
        with open(topology_filepath, 'r') as top_file_used:
            subprocess.call(['/bin/bash','-i','-c','gmx grompp -f {mdp_file} -c {input_coords_current_round} -p {top_file} -o {tpr_filename} -maxwarn 99'.format(mdp_file=alchembed_mdp_filename, input_coords_current_round=input_coord_file, top_file=topology_filepath, tpr_filename=tpr_filename)])

        #run the alchembed 'simulation' on a single core
        print 'starting alchembed simulation for round ', round_number
        output_deffnm = alchembed_mdp_filename[:-4]
        mdrun_exit_code = subprocess.call(['/bin/bash','-i','-c','gmx mdrun -v -stepout 100 -s {tpr_filename} -deffnm {output_deffnm} -nt 1'.format(tpr_filename=tpr_filename, round_number=round_number, output_deffnm=output_deffnm)])
        if mdrun_exit_code == 0: #mdrun didn't crash
            input_coord_file = output_deffnm + '.gro' # set coord file for next round
        else: #mdrun crashed for whatever reason, so the input file remains the same in the next round...
            pass


        print 'finished alchembed simulation for round ', round_number
        round_number += 1

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    run_steric_resolution_loop(input_coord_file = args.input_coord_file_path, index_list = args.index_list, residue_names_list = args.residue_names_list, cutoff = args.cutoff, list_particles_per_residue = args.list_particles_per_residue, output_path = args.output_path, alchembed_b_value = args.alchembed_b_value, alchembed_resolution = args.alchembed_resolution, alchembed_steps = args.alchembed_steps, alchembed_alpha = args.alchembed_alpha, alchembed_dt = args.alchembed_dt, topology_filepath = args.topology_filepath, tpr_filepath = args.tpr_filepath, resume_file = args.resume_file)

