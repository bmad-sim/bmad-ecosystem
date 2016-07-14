#!/usr/bin/env python
import os, sys, shutil, glob, tempfile, time
from bbu import bbu_main, find_threshold, drscan, phase_scan  #imports bbu package in user python path

#ALL USER SETTINGS:

#bbu settings
bbu_par = {  \
#'lat_filename': "'~/nfs/linux_lib/bsim/bbu/examples/oneturn_lat.bmad'",  
#'lat_filename': "'~/nfs/linux_lib/bsim/bbu/lattice/mlc/mlc.lat'",   # Make sure the correct lattice is called 
#'lat_filename': "'~/nfs/linux_lib/bsim/bbu/lattice/1pass_lat.bmad'",   # Make sure the correct lattice is called 
'lat_filename': "'~/nfs/linux_lib/bsim/bbu/lattice/4pass_lat.bmad'", 
#'lat_filename': "'~/nfs/linux_lib/bsim/bbu/lattice/1pass_20160125_hybrid.lat'", 

'bunch_freq': 1.3e9,                # Freq in Hz.
'limit_factor': 3,                  # Init_hom_amp * limit_factor = simulation unstable limit  !! Must be >2
'simulation_turns_max': 250,       # Must be greater than 10
'hybridize': '.true.',                  # Combine non-hom elements to speed up simulation?
#'hybridize': '.true.',                  # Combine non-hom elements to speed up simulation?
'keep_overlays_and_groups': '.false.',  # Keep when hybridizing?
#'keep_all_lcavities': '.true.',         # Keep "cavities without HOM" when hybridizing?
'keep_all_lcavities': '.false.',         # Keep "cavities without HOM" when hybridizing?
'current': 'temp_curr',               # Starting current (amps) set from bbu.py; default is 1A, set in bbu_params.py
'rel_tol': 1e-3,                    # Final threshold current accuracy
'elname': "'T1'",                     # Element to step length for DRSCAN, default is arc
'nstep': 50,                        # Number of steps for DRSCAN
'lat2_filename': "''",                 # Changed for drscan
'nrep': 5,                       # Number of times to repeat threshold calculation
'ran_seed': 100,                      # Set specific seed if desired (0 uses system clock)
'ran_gauss_sigma_cut': 3             # If positive, limit ran_gauss values to within N sigma
}

###############################################################################
# python parameters
py_par = {  \
'exec_path':'/home/wl528/nfs/linux_lib/production/bin/bbu',   # Normal version
#'exec_path':'/home/wl528/nfs/linux_lib/debug/bin/bbu',        # Debug version ( slow )

'threshold_start_curr': 0.1,  # Initial test current for all modes

############## Parameters for DR_SCAN  mode:   #################################

'ndata_pnts_DR': 11,   # integer >=2 required

# For something like the PRSTAB 7, Fig. 3, try startarctime = 4.028E-9, endarctime = 4.725E-9, bunch_freq = 1.3E9
'start_dr_arctime': 4.028*10**-9,  
'end_dr_arctime': 4.725*10**-9,  
#'start_dr_arctime': 4.05*10**-9,  
#'start_dr_arctime': 2.54229*10**-7,  
#'start_dr_arctime': 2.56*10**-7,  
#'end_dr_arctime': 2.56*10**-7,  
#'end_dr_arctime': 4.05*10**-9,  
'plot_drscan': True,   # Create a python plot?

############## Parameters for PHASE_SCAN  mode:   ##################################

'ndata_pnts_PHASE': 1,   # integer >=1 required
'start_phase': 1.00,    # for n_data_pnts >= 2
'end_phase': 2.00,     # for n_data_pnts >= 2
'ONE_phase': 0,       # for n_data_pnts = 1 ONLY
'plot_phase_scan': True,   # Create a python plot ?

######## Parameters for THRESHOLD mode:  ######################################

#'random_homs': True,   # If True, will randomly assign HOMs in 'hom_dir' to the cavities
'random_homs': False,  # Set to False if HOMs are PRE-assigned

# If random_homs is False, hom_dir is not used 
# Make sure hom_dir has the desired HOMs to be RANDOMLY assigned
'hom_dir_number': 250,  # Can be 125,250,500, or 1000 (micrometer). Make sure hom_dir has consistent name!!! 
'hom_dir': '/home/wl528/nfs/linux_lib/bsim/bbu/threshold/HOM_lists_250mm/',
'hom_fixed_file_number': -1,   # Will be the 5th argument from user (if given) to assign all cavities with the same HOMs
'temp_dir': ''   # Will be created, LEAVE IT EMPTY
}

# This runs the code below from the command line:
# python3 .../test_run.py #Thresholds #ID '~/nfs/linux_lib/bsim/bbu/target_directory/'

def main(argv):
  print(time.time())
# Decides which mode the program runs based on the number of arguments
  if (len(sys.argv) == 1):
    print('1 argumnet (including python script) given. DR-SCAN mode.')
    bbu_par['lat_filename']= "'~/nfs/linux_lib/bsim/bbu/examples/oneturn_lat.bmad'"
    mode = 'dr_scan'
    working_dir = os.getcwd() # current directory
    print('WORKING DIR ',os.getcwd())
  
  if (len(sys.argv) == 2):
    print('2 argumnet (including python script) given. PHASE-SCAN mode.')
    mode = 'phase_scan'
    py_par['ONE_phase'] = sys.argv[1]   # If ndata_pnts >=2, ONE_phase is NOT used
    if (py_par['ndata_pnts_PHASE']==1):
      print('Scan for one phase only: ', py_par['ONE_phase'])
    working_dir = os.getcwd() # current directory
    print('WORKING DIR ',os.getcwd())
  
  if (len(sys.argv) >= 4 ):  
    print ('4 or more arguments (including python script) given, threshold (current) mode.')
    n_run = 1
    n_run = int(sys.argv[1])  # Number of times to run
    f_n  = int(sys.argv[2])  # File number to be saved as 
    working_dir = sys.argv[3]  # Location to store output files
    mode = 'threshold'
    if (len(sys.argv) == 5):  
      print ('The 5th argument given =  the HOM_file_number in "hom_dir" used to assign the HOMs for all cavities.')
      print ('Make sure py_par["random_homs"] is TRUE. (Although the assignment is not "random".) ')
      py_par['hom_fixed_file_number'] = int(sys.argv[4]) 
######################################################################

  user_lattice = bbu_par['lat_filename'] 
  print('Lattice name:', bbu_par['lat_filename'])
  # Make a temporary directory (will be removed after program ends properly) 
  # Temporary directory has a randomly-generated name
  py_par['temp_dir'] = make_tempdir( 1, working_dir )  
  os.chdir( py_par['temp_dir'])
  print('Temporary directory created:', py_par['temp_dir']) 
 
  bbu_par['lat_filename'] = '\''+os.path.join(py_par['temp_dir'],'temp_lat.lat')+'\'' 

## creates bbu_template.init which stores all bbu_par
  find_threshold.keep_bbu_param( bbu_par, py_par['temp_dir'] )
  find_threshold.prepare_lat( py_par, user_lattice )  


  if (mode == 'threshold'):						 
    for i in range(n_run):
      find_threshold.keep_bbu_param( bbu_par, py_par['temp_dir'] )
      # This will put rand_assign_homs.bmad in the working dir, include this file in the lattice file
      
      if (py_par['random_homs']):  # If HOMs are not assigned to cavities yet, random HOMs will be assigned for each new job
        find_threshold.prepare_HOM( py_par )  

        # Save (append) the HOM assignments  
        f2 = open(os.path.join(py_par['temp_dir'],'rand_assign_homs.bmad'), 'r')
        contents2 = f2.readlines()
        f2.close()
     
        with open('rand_assign_homs_'+str(f_n)+'.bmad', 'a') as myfile2:
          myfile2.write('\nFor threshold run# '+ str(i)+ ' the (random) assignments were:\n')
          for line2 in contents2:
            myfile2.write(line2)
          myfile2.close()

      bbu_main.single_threshold ( py_par )  # This loop runs BBU and fills thresholds.txt over the runs
    
    # (threshold run(s) end here)
    # Save "bbu_threshold_fn.txt" and "rand_assign_homs_fn.bmad"(if exist) in the working directory
    # os.chdir(os.path.dirname(working_dir))
    os.chdir(os.path.dirname(sys.argv[3]))
    print('Saving the result (threshold current in A) to the working directory...') 
    shutil.copyfile(os.path.join(py_par['temp_dir'],'thresholds.txt'), 'bbu_thresholds_'+str(f_n)+'.txt')
    #shutil.copy(os.path.join(py_par['temp_dir'],'thresholds.txt'), sys.argv[3])
    
    
    if (py_par['random_homs']):
     
      with open('bbu_thresholds_'+str(f_n)+'.txt', 'a') as myfile3:
        myfile3.write('\n')
        myfile3.write('"Random" HOM assignments stored in rand_assign_homs_'+str(f_n)+'.bmad')
        myfile3.close()
      print('Saving rand_assign_homs_', str(f_n),'.bmad to ', str(working_dir) )
      shutil.copyfile(os.path.join(py_par['temp_dir'],'rand_assign_homs_'+str(f_n)+'.bmad'), 'rand_assign_homs_'+str(f_n)+'.bmad')

    else:
      if (not os.path.isfile('assignHOMs.bmad')):
        print('The file with user-assigned-HOMs information was not found!')
        print('The user needs to record the HOMs assigned manually!') 
      else:
        print('Saving user-assigned-HOM information from local assignHOMs.bmad...')
        f3 = open('assignHOMs.bmad', 'r')
        contents3 = f3.readlines()
        f3.close()
     
        with open('bbu_thresholds_'+str(f_n)+'.txt', 'a') as myfile3:
          myfile3.write('\n')
          myfile3.write('HOM assignment from "assignHOMs.bmad": ')
          for line3 in contents3:
            myfile3.write(line3)
          myfile3.close()



    #f = open(os.path.join(py_par['temp_dir'],'thresholds.txt'), 'r')
    #contents = f.readlines()
    #f.close()
    #with open('bbu_thresholds_'+str(f_n)+'.txt', 'a') as myfile:
    #  for line in contents:
    #    myfile.write(line)
    #  myfile.close()
  
  ## for dr scan
  if(mode == 'dr_scan'):
    bbu_main.drscanner( py_par ) 
    #os.chdir(os.path.dirname(working_dir))
    os.chdir(working_dir) # Go back to the working dir from temp dir
    # save the result ( Ith vs tr/tb data)
    print('Copying thresh_v_trotb.txt to ', working_dir) 
    shutil.copyfile(os.path.join(py_par['temp_dir'],'thresh_v_trotb.txt'), 'thresh_v_trotb.txt')

  ## for phase scan
  if(mode == 'phase_scan'):
    bbu_main.phase_scanner( py_par ) 
    print(working_dir)
    # Re-specify the directory to save the files, if necessary
    #working_dir = '~/nfs/linux_lib/bsim/bbu/cbeta_test/Phase_TTT2/' # Go back to the working dir from temp dir
    print("HAHAHA",working_dir)
    os.chdir(working_dir) # Go back to the working dir from temp dir
    # save the result ( Ith vs tr/tb data)
    print('Copying thresh_v_phase.txt to ', working_dir) 
    shutil.copyfile(os.path.join(py_par['temp_dir'],'thresh_v_phase.txt'), 'thresh_v_phase_'+str(py_par['ONE_phase'])+'.txt')
  # clean up the temporary directory for both modes
  # Comment out these two lines if you want to keep the temporary files for debugging 
  print('Deleting temporary directory and its files...') 
  cleanup_workdir( py_par['temp_dir'] )

  print(time.time())



#==========================================================
def make_tempdir ( namecode, dir ):
##################### Makes the temporary directory 
  my_tdir = tempfile.mkdtemp(str(namecode), 'bbu_temp_', dir)
  tdir = os.path.join(dir, my_tdir)
  return tdir


#==========================================================
def cleanup_workdir(tempdir):
# Remove the temporary directory
  if (not os.path.exists(tempdir)):
    print('Error: workdir was already removed!: ', tempdir)
  else:
    shutil.rmtree(tempdir)
    

    

# Boilerplate
if __name__ == "__main__":
  print ( sys.argv )
  print ( sys.argv[0] )
  main(sys.argv[1:]) 
  
