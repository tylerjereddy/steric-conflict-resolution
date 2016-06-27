import collections

def generate_mdp(resolution='CG', output_filename = 'alchembed.mdp', b = 2, steps = 1000, delta_lambda = 0.001, alpha = 0.1, dt = 0.01):
    '''Produce custom alchembed / GROMACS .mdp file based on Phil F's github tutorial: https://github.com/philipwfowler/alchembed-tutorial
    b can be either 1 or 2, but should default to 2 based on suggestion from Alchembed paper
    steps defaults to 1000 
    delta_lambda defaults to 0.001'''
    if resolution == 'CG': #Coarse Grained system
        mdp_lines =      collections.OrderedDict([('integrator' , 'md'),
                        ('define'               , '-DFLEXIBLE'),
                        ('tinit'                , '0.0'),
                        ('dt'                   , '0.01'),
                        ('nsteps'               , '1000'),
                        ('nstxout'              , '100000'),
                        ('nstvout'              , '100000'),
                        ('nstfout'              , '10'),
                        ('nstlog'               , '10'),
                        ('nstenergy'            , '10000'),
                        ('nstxtcout'            , '10'),
                        ('xtc_precision'        , '1000'),
                        ('coulombtype'          , 'Reaction_Field'),
                        ('rlist'                , '1.2'),
                        ('rcoulomb_switch'      , '0.0  '),
                        ('rcoulomb'             , '1.2'),
                        ('epsilon_r'            , '15'),
                        ('epsilon_rf'           , '0   '),
                        ('vdw_type'             , 'cutoff '),
                        ('rvdw_switch'          , '0.9'),
                        ('rvdw'                 , '1.2'),
                        ('cutoff-scheme'        , 'verlet'),
                        ('coulomb-modifier'     , 'potential-shift-verlet'),
                        ('vdw-modifier'         , 'potential-shift-verlet'),
                        ('verlet-buffer-drift'  , '0.005'),
                        ('tcoupl'               , 'Berendsen'),
                        ('tc-grps'              , 'SYSTEM'),
                        ('tau_t'                , '1.0 '),
                        ('ref_t'                , '310 '),
                        ('gen_vel'              , 'yes'),
                        ('gen_temp'             , '310'),
                        ('gen-seed'             , '-1'),
                        ('free_energy'          , 'yes'),
                        ('init_lambda'          , '0.00'),
                        ('delta_lambda'         , '1e-3'),
                        ('sc-alpha'             , '0.1000'),
                        ('sc-power'             , '1'),
                        ('sc-r-power'           , '6'),
                        ('couple-moltype'       , 'SYSTEM'),
                        ('couple-lambda0'       , 'none'),
                        ('couple-lambda1'       , 'vdw')])

    #modify mdp values based on user input, if applicable:
    mdp_lines['sc-power'] = str(b) # 2 by default, but can change to 1
    mdp_lines['nsteps'] = str(steps) # 1000 by default, can switch to 10,000 if crashing
    mdp_lines['delta_lambda'] = str(delta_lambda) #usually 1/steps so that lambda = 1 is reached, but sometimes may want to stop before lambda = 1 if system is really tricky
    mdp_lines['sc-alpha'] = str(alpha) 
    mdp_lines['dt'] = str(dt) 

    with open(output_filename, 'w') as output_file:
        for key, value in mdp_lines.iteritems():
            output_file.write(key + ' = ' + value + '\n')

if __name__ == '__main__':
    generate_mdp()
