#!/usr/bin/env python
import os, sys, shutil
from bbu import bbu_main, find_threshold, drscan  #imports bbu package in user python path

#ALL USER SETTINGS:

#bbu settings
bbu_par = {  \
#'lat_filename': "'home/mt728/nfs/linux_lib/bsim/bbu/examples/multicavity_lat.bmad'",  
#'lat_filename': "'~/nfs/linux_lib/bsim/bbu/examples/oneturn_lat.bmad'",  
'lat_filename': "'~/nfs/linux_lib/bsim/bbu/examples/mlc.lat'",  
'bunch_freq': 1.3e9,                # Freq in Hz.
'limit_factor': 3,                   # Init_hom_amp * limit_factor = simulation unstable limit  !! Must be >2
'simulation_turns_max': 1000,       # Must be greater than 10
'hybridize': '.true.',                  # Combine non-hom elements to speed up simulation?
'keep_overlays_and_groups': '.false.',  # Keep when hybridizing?
'keep_all_lcavities': '.true.',         # Keep when hybridizing?
'current': 'temp_curr',               # Starting current (amps) set from bbu.py; set current in bbu_params.py
'rel_tol': 1e-3,                    # Final threshold current accuracy
'elname': "'T1'",                     # Element to step length for DRSCAN, default is arc
'nstep': 50,                        # Number of steps for DRSCAN
'lat2_filename': "''",                 # Changed for drscan
'nrep': 15,                       # Number of times to repeat threshold calculation
'ran_seed': 0,                      # Set specific seed if desired (0 uses system clock)
'ran_gauss_sigma_cut': 3             # If positive, limit ran_gauss values to within N sigma
}

#py settings
py_par = {  \
'exec_path':'/home/mt728/nfs/linux_lib/production/bin/bbu', 
'ndata_pnts': 100, 
'threshold_start_curr': 1,
# For something like the PRSTAB 7, Fig. 3, try startarctime = 4.028*10**-9, endarctime = 4.725*10**-9
'start_dr_arctime': 4.028*10**-9,  
'end_dr_arctime': 4.725*10**-9,  
'plot_drscan': True,   # Creates a python plot
'random_homs': True,    # Will create a lattice file which randomly assigns hom data files to each cavity of the lattice
                        # Set to False for DR scan
'hom_dir': '/home/mt728/nfs/linux_lib/bsim/bbu/test/cut_HOM_lists/',
#'hom_dir': '/home/mt728/nfs/linux_lib/bsim/bbu/test/HOM_lists/', 
'temp_dir': ''
}

# FOR EASE OF USE, KEEP THE COMMENTED LINES BELOW AND INCLUDE THEM WHEN DESIRED
'''
#  For drscan, run the script with no arguments, include either of these two lines:
#bbu_main.drscanner( py_par )    
'''
# This runs the code below from the command line:
#python3 examples/examp_bbu.py #Thresholds #ID '/home/mt728/nfs/linux_lib/bsim/bbu/test/'

def main(argv):

  n_jobs = 1
  if (len(sys.argv) < 2):
    print('Dr scan. No arguments given.')
    n_jobs = 1
    mode = 'dr_scan'
    working_dir = os.getcwd()
  if (len(sys.argv) > 1):  
    n_jobs = int(sys.argv[1])  # Number of times to run 
    f_n  = int(sys.argv[2])  # Output large text file number 
    working_dir = sys.argv[3]  # Location for temp directory
    mode = 'threshold'

  user_lattice = bbu_par['lat_filename'] 
  py_par['temp_dir'] = find_threshold.make_tempdir( 1, working_dir )  
  # cd to temp directory to call bbu
  os.chdir( py_par['temp_dir'])
  bbu_par['lat_filename'] = '\''+os.path.join(py_par['temp_dir'],'temp_lat.lat')+'\''  # make_assign() will include the user-specified lat file and the random_hom latfile
  find_threshold.make_init( bbu_par, py_par['temp_dir'] )
  find_threshold.make_assign( py_par, user_lattice, 'need_names' )  # In threshold mode uts rand_assign_homs.bmad in working dir, but must run bbu to get cavity names in hom_info.txt

  if (mode == 'threshold'):								    # Here, bbu will complain about finding no LR wakes. Its ok. 
    for i in range(n_jobs):
      find_threshold.make_init( bbu_par, py_par['temp_dir'] )
      find_threshold.make_assign( py_par, user_lattice, 'have_names' )  # This will put rand_assign_homs.bmad in the working dir, include this file in the lattice file
      bbu_main.single_threshold ( py_par )  # This loop fills thresholds.txt file

      f2 = open(os.path.join(py_par['temp_dir'],'rand_assign_homs.bmad'), 'r')
      contents2 = f2.readlines()
      f2.close()
      with open('rand_assign_homs_'+str(f_n)+'.bmad', 'a') as myfile2:
        myfile2.write('\nFor threshold # '+ str(i)+ ' the assignments were:\n')
        for line2 in contents2:
          myfile2.write(line2)
        myfile2.close()

    os.chdir(os.path.dirname(working_dir))
    shutil.copyfile(os.path.join(py_par['temp_dir'],'rand_assign_homs_'+str(f_n)+'.bmad'), 'rand_assign_homs_'+str(f_n)+'.bmad')
    f = open(os.path.join(py_par['temp_dir'],'thresholds.txt'), 'r')
    contents = f.readlines()
    f.close()
    with open('bbu_thresholds_'+str(f_n)+'.txt', 'w') as myfile:
      for line in contents:
        myfile.write(line)
      myfile.close()
  
  if(mode == 'dr_scan'):
    bbu_main.drscanner( py_par ) 
    os.chdir(os.path.dirname(working_dir))
    shutil.copyfile(os.path.join(py_par['temp_dir'],'thresh_v_trotb.txt'), 'thresh_v_trotb.txt')


  find_threshold.cleanup_workdir( py_par['temp_dir'] )


    
    
    

# Boilerplate
if __name__ == "__main__":
  print ( sys.argv )
  main(sys.argv[1:]) 
  
