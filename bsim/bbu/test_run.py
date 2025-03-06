#!/usr/bin/env python3
import os, sys, shutil, glob, tempfile, time
from bbu import bbu_main, find_threshold, drscan, phase_scan  #imports bbu package in user python path

#ALL USER SETTINGS:

#bbu settings
bbu_par = {  \
# Make sure the correct lattice is called 
'lat_filename':"'$DIST_BASE_DIR/bsim/bbu/examples/oneturn_lat.bmad'",
#'bunch_freq': 1.3e9,               # Freq in Hz.
'bunch_freq': 299792458/0.5/50.0,   # Freq in Hz.
'limit_factor': 3,                  # Init_hom_amp * limit_factor = simulation unstable limit  !! Must be >2
'simulation_turns_max': 30,         # Must be > 10. More turns => more accurate but slower
'hybridize': 'True',                # Combine non-HOM elements to speed up simulation?
'keep_all_lcavities': 'False',      # Keep cavities without HOM when hybridizing (default = false)?
'current': 'temp_curr',             # Test current for the Fortran bbu code. DO NOT MODIFY
#'rel_tol': 1e-2,                   # Final threshold current accuracy. Small => slow

'lat2_filename': "''",              # For DR-scan and phase-scan, LEAVE IT EMPTY
'ran_seed': 100,                    # Set specific seed if desired (0 uses system clock)
'ran_gauss_sigma_cut': 3,           # If positive, limit ran_gauss values to within N sigma
'ramp_on': 'False',
'ramp_n_start': 0,
'ramp_pattern': '1.0, 1.0',
}

###############################################################################
# python parameters
py_par = {  \
'exec_path':'$DIST_BASE_DIR/production/bin/bbu',   # Production version
'temp_dir': '',                # Will be created, LEAVE IT EMPTY
'threshold_start_curr': 0.1,   # Initial test current for all modes
'final_rel_tol': 1e-2,         # Final threshold current accuracy. Small => slow

############## Parameters for DR_SCAN  mode:   #################################

'ndata_pnts_DR': 31,   # integer >=1 required

# For something like the PRSTAB 7, Fig. 3, try startarctime = 4.028E-9, endarctime = 4.725E-9, bunch_freq = 1.3E9
#'start_dr_arctime': 96.5/1.3e9,  
#'end_dr_arctime': 100.5/1.3e9,  
#'start_dr_arctime': 4.028E-9,  
#'end_dr_arctime': 4.725E-9,  
'start_dr_arctime': 99.5*0.5/299792458,  
'end_dr_arctime': 100.5*0.5/299792458,  


'plot_drscan': True,   # Create a python plot?

############## Parameters for PHASE_SCAN  mode:   ##################################

'ndata_pnts_PHASE': 1,   # integer >=1 required
'start_phase': 0.00,     # for n_data_pnts >= 2
'end_phase': 6.28,       # for n_data_pnts >= 2
'ONE_phase': 0,          # for n_data_pnts = 1 ONLY
'plot_phase_scan': False,   # Create a python plot ?


############## Parameters for PHASE_SCAN_XY  mode:   ##################################
'phase_x': 0, 
'phase_y': 0,   
'xy_coupled': 1,        # 1=YES, 0=NO
######## Parameters for THRESHOLD mode:  ######################################

'random_homs': True,   # If True, will (randomly) assign new HOMs in 'hom_dir' to the cavities
#'random_homs': False,  # Set to False if the user wants the PRE-assigned HOMs to be used

# If random_homs is False, hom_dir is not used 
# Make sure hom_dir has the desired HOMs to be RANDOMLY/FIXEDLY assigned
'hom_dir': '$DIST_BASE_DIR/bsim/bbu/hom/',
#'hom_dir_number': 125,  # Can be 125,250,500, or 1000 (micrometer). Make sure hom_dir has consistent name!!! 
'hom_fixed_file_number': -1 # Do not modify 
                            #The 5th argument from user (if given) to assign all cavities with the same HOMs
}

# This runs the code below from the command line:
# python3 .../test_run.py #Thresholds #ID '$DIST_BASE_DIR/bsim/bbu/target_directory/'

def main(argv):
  print(time.time())
  working_dir = os.getcwd() # current directory
  print('WORKING DIR ',os.getcwd())
    
# Decides which mode the program runs based on the number of arguments
  if (len(sys.argv) == 1):
    print('1 argumnet (including python script) given. DR-SCAN mode.')
    bbu_par['lat_filename']= "'$DIST_BASE_DIR/bsim/bbu/examples/oneturn_lat.bmad'"
    mode = 'dr_scan'
    
  if (len(sys.argv) == 2):
    print('2 argumnets (including python script) given. PHASE_SCAN mode.')
    mode = 'phase_scan'
    py_par['ONE_phase'] = sys.argv[1]   # If ndata_pnts >=2, ONE_phase is NOT used
    if (py_par['ndata_pnts_PHASE']==1):
      print('Scan for one phase only: ', py_par['ONE_phase'])
        
  if (len(sys.argv) == 3):
    print('3 argumnets (including python script) given. PHASE_XY_SCAN mode.')
    mode = 'phase_xy_scan'
    py_par['phase_x'] = sys.argv[1]  
    py_par['phase_y'] = sys.argv[2]  
    print('Scan for the XY phase combination: ', py_par['phase_x'], ', ',py_par['phase_y'])

  
  if (len(sys.argv) >= 4 ):  
    print ('4 or more arguments (including python script) given, threshold (current) mode.')
    n_run = int(sys.argv[1])  # Number of times to run
    f_n  = int(sys.argv[2])  # File number to be saved as 
    output_dir = sys.argv[3]  # Location to store output files
 
    if not os.path.exists(output_dir):
       # Create a new output directory because it does not exist
       os.makedirs(output_dir)
      
    mode = 'threshold'
    if (len(sys.argv) == 5):  
      #The 5th argument given =  the HOM_file_number in "hom_dir" used to assign the HOMs for all cavities.
      print ('CAUTION!! All cavities will be assigned with the SAME HOM based on the 5th argument')  
      print ('Make sure py_par["random_homs"] is TRUE. (Although the assignment is not truly "random".) ')
      py_par['hom_fixed_file_number'] = int(sys.argv[4]) 
######################################################################

  user_lattice = bbu_par['lat_filename'] 
  print('Lattice name:', bbu_par['lat_filename'])
  # Create a temp_dir to save all temporary files  
  # The temp_dir has a randomly-generated name
  py_par['temp_dir'] = make_tempdir( 1, working_dir )  
  os.chdir( py_par['temp_dir'])
  print('Temporary directory created:', py_par['temp_dir']) 
 
  bbu_par['lat_filename'] = '\''+os.path.join(py_par['temp_dir'],'temp_lat.lat')+'\'' 

  # creates bbu_template.init which stores all bbu_par
  find_threshold.keep_bbu_param( bbu_par, py_par['temp_dir'] )
  #find_threshold.prepare_lat( py_par, user_lattice, working_dir )  
  find_threshold.prepare_lat( py_par, user_lattice )  


  if (mode == 'threshold'):						 
    for i in range(n_run):
      find_threshold.keep_bbu_param( bbu_par, py_par['temp_dir'] )
     
      # This will put rand_assign_homs.bmad in the working dir, include this file in the lattice file
      if (py_par['random_homs']):  
        # If HOMs are not assigned to cavities yet, random HOMs will be assigned for each new job
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

      else:
        print("Looking for local assignHOMs.bmad...")
        print("If HOMs already assigned with the lattice, leave assignHOMs.bmad blank to avoid over-write \n")
        f_lat2 = open(os.path.join(py_par['temp_dir'],'temp_lat.lat'), 'a')
        f_lat2.write("call, file = \'"+os.path.join(working_dir,'assignHOMs.bmad')+"\'\n")
        f_lat2.close()

      bbu_main.single_threshold ( py_par )  # This loop runs BBU and fills thresholds.txt over the runs
    

    # Save the result ( the final test current attempted ) to the output directory
    file_to_save = os.path.join(os.path.expandvars(py_par['temp_dir']),'thresholds.txt')
    assert os.path.isfile(file_to_save), "The result file is MISSING!"
    file_destination = os.path.join(output_dir, f'bbu_thresholds_N_{n_run}_fn_{f_n}.dat')
    os.chdir( working_dir )
    shutil.move( file_to_save, file_destination )

    # Save the HOM assignments, if available, 
    # The assignments are saved with the result (Ith) in rand_assign_homs_fn.bmad  
    if (py_par['random_homs']):
        file_to_save = os.path.join(py_par['temp_dir'],'rand_assign_homs_'+str(f_n)+'.bmad')
    else:
        file_to_save = os.path.join(working_dir,'assignHOMs.bmad')
    assert os.path.isfile(file_to_save), "The HOM file is MISSING!"
    file_destination = os.path.join(output_dir,  f'HOM_assignment_N_{n_run}_fn_{f_n}.dat')
    os.chdir( working_dir )
    shutil.copy( file_to_save, file_destination )
      
################ End of threshold mode #################################

  ## for DR scan
  if(mode == 'dr_scan'):
    bbu_main.drscanner( py_par ) 
    #os.chdir(os.path.dirname(working_dir))
    os.chdir(working_dir) # Go back to the working dir from temp dir
    # save the result ( Ith vs tr/tb data)
    print('Copying thresh_v_trotb.txt to ', working_dir) 
    shutil.copyfile(os.path.join(py_par['temp_dir'],'thresh_v_trotb.txt'), 'thresh_v_trotb.txt')

  ## for phase scan
  if(mode == 'phase_scan'):
    print('======  PHASE_SCAN MODE ======')
    bbu_main.phase_scanner( py_par ) 
    print(working_dir)
    # Re-specify the directory to save the files, if necessary
    os.chdir(working_dir) # Go back to the working dir from temp dir
    # save the result ( Ith vs phase data)
    print('Copying thresh_v_phase.txt to ', working_dir) 
    shutil.copyfile(os.path.join(py_par['temp_dir'],'thresh_v_phase.txt'), 'thresh_v_phase_'+str(py_par['ONE_phase'])+'.txt')
  
  
  ## for phase_XY scan
  if(mode == 'phase_xy_scan'):
    print('======  PHASE_XY_SCAN MODE ======')
    bbu_main.phase_xy_scanner( py_par ) 
    os.chdir(working_dir) # Go back to the working dir from temp dir
    # save the result ( Ith, phasex, phasey data)
    print('Copying thresh_v_phase_xy.txt to ', working_dir) 
    shutil.copyfile(os.path.join(py_par['temp_dir'],'thresh_v_phase_xy.txt'), 'thresh_v_phase_'+str(py_par['phase_x'])+'_'+str(py_par['phase_y'])+'.txt')
  
  
  # For any mode, clean up the temporary directory
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
  
