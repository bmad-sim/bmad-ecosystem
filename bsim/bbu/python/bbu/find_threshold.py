#!/usr/bin/env python
import subprocess, os, tempfile, shutil
import glob, math, random

#===============================================================
# Report stability of a bbu run 
# In: d[v_gain] and d[lostbool] from dictionary "d" (parsed for_py.txt) 
# Out: 0(unstable) or 1(stable)
def get_stability(gain, lost_bool):
####################### 
  ## This number is critical in determining the stability of the test current 
  ## Ideally this number is 1.0, but there is noise
  ## A stable test current may have voltage noise up to (by observation) 3.0%
  criterion = 1.03
  
  if (lost_bool):  
    is_stable = 0
    print ('Beam is UNstable at this current (Beam LOST, Gain unknown)')
  elif (gain > criterion): 
    ## A stable test current may have voltage noise up to (by observation) 1.0%
    is_stable = 0
    print ('Beam is UNstable at this current (Beam not lost, Gain > '+ str(criterion)+')')
  else:
    is_stable = 1
    print ('Beam is stable at this current (Beam not lost, Gain < '+ str(criterion)+')')

  return is_stable

#===============================================================

# Update dictionary t when finding Ith
# In: dictionary t, stability of the current run, tolerance in finding Ith 
# charge0 is the temporary lower bound, charge1 is the higher bound
def calc_new_charge( t, bool_stable):

  if (not bool_stable):  # Current is unstable at this charge, update upper-bound
    print ('Test current NOT STABLE, reset charge1')
    t['charge1'] = t['bunch_charge']
    t['growth_rate1'] = t['growth_rate']
  else:     # Current is stable at the charge, update lower-bound
    print ('Test current STABLE, reset charge0')
    t['charge0'] = t['bunch_charge']
    t['growth_rate0'] = t['growth_rate']

  if (t['charge1'] > 0):   # Unstable charge has been found
      t['bunch_charge'] = (t['charge0'] + t['charge1']) / 2
  else:  # Still searching for an unstable charge
      t['bunch_charge'] = t['bunch_charge'] * 2

#==========================================================
def keep_bbu_param ( bbu_params, temp_dir ):
# Creates temporary bbu_template.init to store user-defined bbu parameters
# In: the bbu parameters and the temporary directory
#####################################################
  bt = open(os.path.join(temp_dir,'bbu_template.init'),'w')
  bt.write('&bbu_params\n')
  for key in bbu_params:
    bt.write('  bbu_param%'+ key + ' = ' +str( bbu_params[key])+'\n' )
  bt.write('/\n')  # A CR at the end of file is required if the user use ifort compiler
  bt.close()


#==========================================================
def prepare_lat (py_par, lat_file):
### Create "temp_lat.lat" with command lines to call the lattice file 
###########################################################
  print(lat_file)
  f_name = os.path.join(py_par['temp_dir'],'temp_lat.lat')
  f_lat = open(f_name, 'w')
  f_lat.write('call, file = '+lat_file+'\n')
  f_lat.close()

#==========================================================
def prepare_HOM ( py_par ):
# For threshold mode only, assign RANDOM HOMs.
#    - BBU will run once (may fail if no pre-assignment) to output hom_info.txt, which contains cavity names
#    - HOMs are assigned from py_par['hom_dir'] 
#    - Create rand_assign_homs.bmad
#    - Commands to call rand_assign_homs.bmad are added in temp_lat.lat, called by bbu.init
###############################################################################
   
## For threshold mode:
  # Runs bbu ONCE to get the cavity names 
  # The names are stored in hom_info.txt (check bbu_program.f90) 
  print(' Running BBU ONCE to generate hom_info.txt ')
  print(' If no fail, the user has pre-assigned HOMs, and will be overwritten by the random HOMs.')
  print(' To prevent overwritting, set "py_par["random_homs"]" to FALSE before running. ')
  run_bbu( py_par['threshold_start_curr'], py_par, 'threshold' )  

  # Get all cavity names from hom_info.txt
  f = open(os.path.join(py_par['temp_dir'], 'hom_info.txt'),'r')     
  # Make sure hom_dir contains ONLY hom data files ended in "dat"
  hom_files = glob.glob( os.path.join( os.path.expandvars(py_par['hom_dir']), '*dat') ) 
  print('The number of cavity HOM files located:', len(hom_files))
  cav_names = []
  #Skip the first line, which says "cavity_name"
  next(f)
  for line in f:        
    s = line.split()
    seen = False
    if len(s) > 0:
      cav_name = s[0]
      # Exit so no double-count any cavity for multipass 
      for i in range(len(cav_names)):
        if (cav_name == cav_names[i]):
          seen = True
      if (not seen): cav_names.append(cav_name) 


  # Create HOM assignement file
  f_assign = open(os.path.join(py_par['temp_dir'],'rand_assign_homs.bmad'), 'w')
  for r in range(len(cav_names)):            

    # If the user doesn't specifiy the HOM_file_number (the 5th argument )
    if(py_par['hom_fixed_file_number'] == -1):  
      # Choose a file at random for each cavity 
      i = random.randrange(0, len(hom_files))
      cav_string = cav_names[r]+'[lr_wake_file] = '+os.path.join(py_par['hom_dir'],hom_files[i])+'\n' 
      print('Random HOMs being assigned to ' + cav_names[r]) 
    else:
      i = py_par['hom_fixed_file_number'] 
      ################# The default HOM file name is set here 
      ##### Modify this line to change it
      #target_hom_file = os.path.join(py_par['hom_dir'], ('vhoms_'+str(py_par['hom_dir_number'])+'mm_cavity_'+str(i)+'.dat'))
      target_hom_file = os.path.join(py_par['hom_dir'], ('cavity_'+str(i)+'.dat'))
      # Check if the hom_file specificed exists
      if (not os.path.isfile(target_hom_file)):
        print('CANNOT find the HOM file: ' + target_hom_file)
        print('Instead, random homs assigned to '+ cav_names[r])
        i = random.randrange(0, len(hom_files)) 
        cav_string = cav_names[r]+'[lr_wake_file] = '+os.path.join(py_par['hom_dir'],hom_files[i])+'\n' 
      else: 
        print('HOMs in '+target_hom_file+' assigned to '+ cav_names[r]) 
        cav_string = cav_names[r]+'[lr_wake_file] = '+ target_hom_file +'\n' 
    f_assign.write(cav_string)
  f_assign.close()
    
  # APPEND (not overwrite!!) command lines to call the assignment file
  f_lat2 = open(os.path.join(py_par['temp_dir'],'temp_lat.lat'), 'a')
  f_lat2.write("call, file = \'"+os.path.join(py_par['temp_dir'],'rand_assign_homs.bmad')+"\'\n")
    
  f_lat2.close()


#==========================================================
def  run_bbu ( temp_curr, py_par, mode ):
# Prepare/update bbu.init and run the bbu program once for a specific test current 
########################################################
  if ( mode == 'threshold' ):
    
    print ('Updating bbu.init with new temp_curr...')
    # Update bbu.init with bbu_template.init, with new temp current
    template_file = open (os.path.join(py_par['temp_dir'], 'bbu_template.init'), 'r' ) 
    temp_file = open (os.path.join(py_par['temp_dir'],'bbu.init'), 'wt')
    for line in template_file:
      temp_file.write( line.replace( "temp_curr", str(temp_curr)) )  
    template_file.close()
    temp_file.close()
    print ('Running BBU with current ', str(temp_curr), '(A)')

    subprocess.call( py_par['exec_path'], shell = True)  # Run bbu 

  if ( mode == 'drscan' or mode == 'phase_scan' or mode == 'phase_xy_scan'):
    
    print ('Updating bbu.init with new temp_curr and lat2.lat...')
    # Update bbu.init with the new lat2 and temp_curr
    template_file = open(os.path.join(py_par['temp_dir'], 'bbu_template.init'), 'r' ) 
    temp_file = open(os.path.join(py_par['temp_dir'],'bbu.init'), 'wt')   
    ftemp = template_file.read() #import original bbu_par 
    ftemp = ftemp.replace( "temp_curr", str(temp_curr)) 
    ftemp = ftemp.replace( "lat2_filename = ''", "lat2_filename = '" + os.path.join(py_par['temp_dir'],'lat2.lat') + "'")
    temp_file.write(ftemp) # write updated bbu_par in bbu.init 
    template_file.close()
    temp_file.close()

    print ('Subprocess begins!!!  Running BBU with current ', str(temp_curr), '(A)')
    subprocess.call( py_par['exec_path'], shell = True)  # Run bbu
