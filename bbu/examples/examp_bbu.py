#!/usr/bin/env python

from bbu import bbu_main, find_threshold, drscan  #imports bbu package in user python path

#this is the only file a user needs to run bbu

#ALL SETTINGS
#bbu settings
bbu_par = {  \
'lat_filename': "'~/nfs/linux_lib/bsim/bbu/examples/oneturn_lat.bmad'",  
'bunch_freq': 1.3e9,                # Freq in Hz.
'limit_factor': 3,                   # Init_hom_amp * limit_factor = simulation unstable limit  !! Must be >2
'simulation_turns_max': 1000,       # Must be greater than 10
'hybridize': '.true.',                  # Combine non-hom elements to speed up simulation?
'keep_overlays_and_groups': '.true.',  # Keep when hybridizing?
'keep_all_lcavities': '.true.',         # Keep when hybridizing?
'current': 'temp_curr',               # Starting current (amps) set from bbu.py; set current in bbu_params.py
'rel_tol': 1e-3,                    # Final threshold current accuracy
'write_hom_info': '.true.',
'elname': "'T1'",                     # Element to step length for DRSCAN, default is arc
'nstep': 50,                        # Number of steps for DRSCAN
'lat2_filename': "''",                 # Changed for drscan
'nrep': 15,                       # Number of times to repeat threshold calculation
'ran_seed': 5,                      # Set specific seed if desired (0 uses system clock)
'ran_gauss_sigma_cut': 3             # If positive, limit ran_gauss values to within N sigma
}

#py settings
py_par = {  \
'exec_path':'/home/mt728/nfs/linux_lib/production/bin/bbu', 
'ndata_pnts': 65, 
'threshold_start_curr': 1,
'start_dr_arctime': 4.028*10**-9,
'end_dr_arctime': 4.725*10**-9,
'plot_drscan': True   
}

find_threshold.make_init( bbu_par )
bbu_main.single_threshold ( py_par )
#or, for drscan
#bbu_main.drscanner( py_par )

