#!/usr/bin/env python
import os, sys, shutil, glob, tempfile
from bbu import bbu_main, find_threshold, drscan  #imports bbu package in user python path

#ALL USER SETTINGS:

#bbu settings
bbu_par = {  \
#'lat_filename': "'~/nfs/linux_lib/bsim/bbu/examples/oneturn_lat.bmad'",  
#'lat_filename': "'~/nfs/linux_lib/bsim/bbu/lattice/mlc/mlc.lat'",   # Make sure the correct lattice is called 
'lat_filename': "'~/nfs/linux_lib/bsim/bbu/lattice/cbeta_lat.bmad'",   # Make sure the correct lattice is called 
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
'nrep': 5,                       # Number of times to repeat threshold calculation
'ran_seed': 0,                      # Set specific seed if desired (0 uses system clock)
'ran_gauss_sigma_cut': 3             # If positive, limit ran_gauss values to within N sigma
}

#py settings
py_par = {  \
'exec_path':'/home/wl528/nfs/linux_lib/production/bin/bbu', 
#'exec_path':'/home/wl528/nfs/linux_lib/debug/bin/bbu',    
#'ndata_pnts': 100, 
'ndata_pnts': 30, 
'threshold_start_curr': 1,
# For something like the PRSTAB 7, Fig. 3, try startarctime = 4.028E-9, endarctime = 4.725E-9, bunch_freq = 1.3E9
'start_dr_arctime': 4.028*10**-9,  
'end_dr_arctime': 4.725*10**-9,  
'plot_drscan': True,   # Creates a python plot
#'random_homs': True,    #  (Threshold mode for now) Will create a lattice file which randomly assigns hom data files to each cavity of the lattice
'random_homs': False,     # Set to False only when HOMs are PRE-assigned

# Make sure hom_dir has the desired HOMs ( to be assigned randomly to cavities )
      # cut_HOM_lists has only 4 dominant lr_wakes for each assignment
      # HOM_lists has 10, 20, 30, or 40 for each assignment 
'hom_dir': '/home/wl528/nfs/linux_lib/bsim/bbu/threshold/cut_HOM_lists/',
#'hom_dir': '/home/wl528/nfs/linux_lib/bsim/bbu/threshold/HOM_lists/',
'temp_dir': ''   # Will be created, leave it empty
}

# FOR EASE OF USE, KEEP THE COMMENTED LINES BELOW AND INCLUDE THEM WHEN DESIRED
# '''
# #  For drscan, run the script with no arguments, include either of these two lines:
# #bbu_main.drscanner( py_par )    
# '''
# This runs the code below from the command line:
#python3 examples/examp_bbu.py #Thresholds #ID '/home/mt728/nfs/linux_lib/bsim/bbu/test/'

def main(argv):

# Decides which mode the program runs based on the number of arguments
  n_jobs = 1
  if (len(sys.argv) < 2):
    print('Dr scan. No arguments given.')
    #n_jobs = 1
    bbu_par['lat_filename']= "'~/nfs/linux_lib/bsim/bbu/examples/oneturn_lat.bmad'"
    #py_par['random_homs']='False'
    mode = 'dr_scan'
    working_dir = os.getcwd() #/bsim/bbu/examples/
    print('WORKING DIR ',os.getcwd())
  if (len(sys.argv) > 1):  
    print ('3 arguments given, threshold (current) mode')
    n_jobs = int(sys.argv[1])  # Number of times to run 
    f_n  = int(sys.argv[2])  # File number to be saved as 
    working_dir = sys.argv[3]  # Location to store output files
    mode = 'threshold'
######################################################################

  user_lattice = bbu_par['lat_filename'] 
  print('Lattice name:', bbu_par['lat_filename'])
  # Make a temporary directory (will be removed after program ends properly) 
  py_par['temp_dir'] = make_tempdir( 1, working_dir )  
  os.chdir( py_par['temp_dir'])
  print('Temporary directory created:', py_par['temp_dir']) 
 
  bbu_par['lat_filename'] = '\''+os.path.join(py_par['temp_dir'],'temp_lat.lat')+'\'' 

## creates bbu_template.init which stores all bbu_par
  find_threshold.make_init( bbu_par, py_par['temp_dir'] )
## Before we assign HOMs to each cavity, we need the cavity names
## To obtain the cavity names, run bbu once ( will fail, as expected, but it's ok )
## bbu program will save the cavity names into hom_info.txt
  find_threshold.make_assign( py_par, user_lattice, 'need_names' )  


  if (mode == 'threshold'):						 
    for i in range(n_jobs):
      find_threshold.make_init( bbu_par, py_par['temp_dir'] )
      # This will put rand_assign_homs.bmad in the working dir, include this file in the lattice file
      find_threshold.make_assign( py_par, user_lattice, 'have_names' )  
      bbu_main.single_threshold ( py_par )  # This loop runs BBU and fills thresholds.txt file

      # Save (append) the HOM assignments over the loops  
      f2 = open(os.path.join(py_par['temp_dir'],'rand_assign_homs.bmad'), 'r')
      contents2 = f2.readlines()
      f2.close()
      with open('rand_assign_homs_'+str(f_n)+'.bmad', 'a') as myfile2:
        myfile2.write('\nFor threshold # '+ str(i)+ ' the assignments were:\n')
        for line2 in contents2:
          myfile2.write(line2)
        myfile2.close()

    # Save "rand_assign_homs_fn.bmad" and "bbu_threshold_fn.txt" in the working directory
    os.chdir(os.path.dirname(working_dir))
    print('Saving rand_assign_homs_', str(f_n),'.bmad to ','working_dir' )
    shutil.copyfile(os.path.join(py_par['temp_dir'],'rand_assign_homs_'+str(f_n)+'.bmad'), 'rand_assign_homs_'+str(f_n)+'.bmad')

    f = open(os.path.join(py_par['temp_dir'],'thresholds.txt'), 'r')
    contents = f.readlines()
    f.close()
    with open('bbu_thresholds_'+str(f_n)+'.txt', 'a') as myfile:
      for line in contents:
        myfile.write(line)
      myfile.close()
  
  ## for dr scan
  if(mode == 'dr_scan'):
    bbu_main.drscanner( py_par ) 
    #os.chdir(os.path.dirname(working_dir))
    os.chdir(working_dir) # Go back to the working dir from temp dir
    # save the result ( Ith vs tr/tb data)
    print('Copying thresh_v_trotb.txt to ', working_dir) 
    shutil.copyfile(os.path.join(py_par['temp_dir'],'thresh_v_trotb.txt'), 'thresh_v_trotb.txt')


  # clean up the temporary directory for both modes
  # Comment out these two lines if you want to keep the temporary files for debugging 
  print('Deleting temporary directory and its files...') 
  cleanup_workdir( py_par['temp_dir'] )




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
  
