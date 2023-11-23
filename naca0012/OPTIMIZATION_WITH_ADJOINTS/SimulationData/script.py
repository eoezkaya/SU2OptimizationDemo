#!/home/jrottmay/.conda/envs/su/bin/python3.10
import numpy as np 
import copy
import os, glob, sys, time, shutil
import argparse

import SU2.io as su2io
import SU2.opt as su2opt
import SU2.util as su2util

def find_config_file():
    config_file=""
    for filename in os.listdir():
        if filename.endswith(f'.cfg'):
            config_file = filename

    if config_file=="":
        raise ValueError(f"No configuration file found!")
    
    return config_file

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", dest="filename",help="read config from FILE", metavar="FILE", default="")
    parser.add_argument("-g", "--gradient", dest="gradient", default="DISCRETE_ADJOINT",
                      help="Method for computing the GRADIENT (CONTINUOUS_ADJOINT, DISCRETE_ADJOINT, FINDIFF, NONE)", metavar="GRADIENT")
    parser.add_argument("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")

    args = parser.parse_args()
    # Setup 'Project'
    if not args.filename:
        config_file = find_config_file()
    else:
        config_file = args.filename

    config = su2io.Config(config_file)
    config.NUMBER_PART = args.partitions
    config.GRADIENT_METHOD = args.gradient
    config.NZONES = 1

    state = su2io.State()
    state.find_files(config)

    # Copied from shape...
    if config.get('TIME_DOMAIN', 'NO') == 'YES' and config.get('RESTART_SOL', 'NO') == 'YES' and gradient != 'CONTINUOUS_ADJOINT':
        restart_name = config['RESTART_FILENAME'].split('.')[0]
        restart_filename = restart_name + '_' + str(int(config['RESTART_ITER'])-1).zfill(5) + '.dat'
        if not os.path.isfile(restart_filename): 
            sys.exit("Error: Restart file <" + restart_filename + "> not found.")
        state['FILES']['RESTART_FILE_1'] = restart_filename

        if config.get('TIME_MARCHING', 'NO') == 'DUAL_TIME_STEPPING-2ND_ORDER':
            restart_filename = restart_name + '_' + str(int(config['RESTART_ITER'])-2).zfill(5) + '.dat'
            if not os.path.isfile(restart_filename): 
                sys.exit("Error: Restart file <" + restart_filename + "> not found.")
            state['FILES']['RESTART_FILE_2'] = restart_filename   

    # Create project
    myproject = su2opt.Project(config, state)
    myproject.n_dv = sum(myproject.config['DEFINITION_DV']['SIZE'])
    
    myproject.save()

    # Some useful variables
    n_dv = myproject.n_dv

    lb = np.array([float(myproject.config.OPT_BOUND_LOWER)] * n_dv)
    ub = np.array([float(myproject.config.OPT_BOUND_UPPER)] * n_dv)

    # Random Sample
    sample = np.random.randn(n_dv) * (ub - lb) + ub
    
    with open('dv.dat') as f:
        for i in range(20):
            line = f.readline()
            sample[i] = float(line)
        
    print(sample)    

    try:
        # DIRECT
        myproject.obj_f(sample)
        myproject.con_ceq(sample)
        myproject.con_cieq(sample)
        # ADJOINTS
        myproject.obj_df(sample) 
#        myproject.con_dceq(sample) #
#        myproject.con_dcieq(sample) #

    except:
        print("Simulation did not converge!")
