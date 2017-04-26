#!/usr/bin/env python
import subprocess
import os
import random
import tempfile
import shutil
import glob
import math

#===============================================================
# Report stability of a bbu run 
# In: d[v_gain] and d[lostbool] from dictionary "d" (parsed for_py.txt) 
# Out: 0(unstable) or 1(stable)
def get_stability(gain, lost_bool):
####################### 
  if (lost_bool):  
    success = 0
    print ('Beam is UNstable at this current (Beam LOST, Gain unknown)')
  elif (gain > 1.03): 
    ## A stable test current may have voltage noise up to (by observation) 1.0%
    success = 0
    print ('Beam is UNstable at this current (Beam not lost, Gain > 1.03)')
  else:
    success  = 1
    print ('Beam is stable at this current (Beam not lost, Gain < 1.03)')

  return success

#===============================================================

# Update dictionary t when finding Ith
# In: dictionary t, stability of the current run, tolerance in finding Ith 
# charge0 is the temporary lower bound, charge1 is the higher bound
def calc_new_charge( t, bool_stable, rel_tol ):
####################### Get threshold using the charge and growth rate for two tracking points
  #t['charge_old'] = t['bunch_charge']
  #t['growth_rate_old'] = t['growth_rate']

  ## A stable test current may have voltage noise up to (by observation) 1.0%
  if (t['growth_rate'] > math.log(1.03) or not bool_stable):  # Current is unstable at this charge, update upper-bound
    print ('Test current NOT STABLE, reset charge1')
    t['charge1'] = t['bunch_charge']
    t['growth_rate1'] = t['growth_rate']
  else:     # Current is stable at the charge, update lower-bound
    print ('Test current STABLE, reset charge0')
    t['charge0'] = t['bunch_charge']
    t['growth_rate0'] = t['growth_rate']

  #--- This IF stmt aims to find charge_try, a better approx. for bunch_charge -----------------
  #--- charge_try MUST be updated ( to be a positive value ) in this stmt
  #--- bunch_charge should NOT be updated inside this stmt
#  if (t['charge_old'] > 0): # If we have at least two trackings so far
#    if (t['charge0'] > 0 and t['charge1'] > 0):  # Now we can interpolate for the threshold
#      print('Obtained a non-zero charge0 and finite charge1, guessing test current... ')
#      c0 = t['charge0']  
#      g0 = t['growth_rate0']
#      c1 = t['charge1']
#      g1 = t['growth_rate1']
#
#      if (not t['growth_rate_set0'] or not t['growth_rate_set1'] or g0 == g1):
#        t['charge_threshold'] = (c0 + c1) / 2
#        #print ('One or both growth rates were not set, or the two rates are equal (g0 = g1)')
#      else:
#        t['charge_threshold'] = (c0 * g1 - c1 * g0) / (g1 - g0)
#        #print ('Charge threshold guess is now ', t['charge_threshold'])
#
#      # Current to use in the next tracking must be significantly different from c0 and c1.
#      t['charge_try'] = t['charge_threshold']
#      min_delta = max(c0, c1) * rel_tol
#      dc0 = t['charge_threshold'] - c0
#      dc1 = t['charge_threshold'] - c1
#
#      # If threshold guess is closer to c0...
#      if (abs(dc0) < abs(dc1) and abs(dc0) < min_delta):
#        c = t['charge_threshold'] + min_delta*(dc0/abs(dc0))
#        if (abs(c-c1) < abs(c-c0)): # If moved closer to c1 then just use the average
#          t['charge_try'] = (c1 + c0) / 2
#        else:
#          t['charge_try'] = c
#  
#      # If threshold guess is closer to c1...
#      elif (abs(dc1) < abs(dc0) and abs(dc1) < min_delta):
#        c = t['charge_threshold'] + min_delta*(dc1/abs(dc1))
#        print('Nudging charge')
#        if (abs(c-c0) < abs(c-c1)):   #If moved closer to c0 then just use the average
#          t['charge_try'] = (c1 + c0) / 2
#        else:
#          t['charge_try'] = c
#  
#    # If both charges are stable, no unstable point	
#    elif (t['charge1'] < 0):
#      print('No unstable current found yet... ')
#      c0 = t['charge0']  
#      g0 = t['growth_rate0']
#      c1 = t['charge_old']  
#      g1 = t['growth_rate_old']   
#      if (not t['growth_rate_set0'] or not t['growth_rate_set1'] or g0 == g1): 
#        t['charge_threshold'] = 2*c0
#        print('g0 or g1 not set, or g0=g1... ')
#      else:
#        t['charge_threshold'] = (c0 * g1 - c1 * g0) / (g1 - g0)
#        print(str(c0),str(g0),str(c1),str(g1),str(t['charge_threshold']*1.3e9),str(1.1*c0*1.3e9))
#      # If the threshold is less than charge0 then there must be a lot of noise so
#      # assume we are far from the threshold and increase the current by a factor of 10
#      # In any case, demand at least a 10% change
#      if (t['charge_threshold'] < c0):
#        t['charge_try'] = 10 * c0
#      else:
#        t['charge_try'] = max(t['charge_threshold'], 1.1 * c0)
#  
#    # If both charges are UNstable, no stable point
#    else:      
#      print('No stable current found yet... ')
#      c0 = t['charge_old']
#      g0 = t['growth_rate_old']
#      c1 = t['charge1']
#      g1 = t['growth_rate1']
#    
#      if (not t['growth_rate_set0'] or not t['growth_rate_set1'] or g0 == g1): 
#        t['charge_threshold'] = c0 / 2
#        print('g0 or g1 not set, or g0=g1... ')
#      else:
#        t['charge_threshold'] = (c0 * g1 - c1 * g0) / (g1 - g0)
#        print('Interpolating... ')
#        print(str(c0),str(g0),str(c1),str(g1),str(t['charge_threshold']*1.3e9),str(0.9*c1*1.3e9))
#    
#      # Demand at least a 10% change, but if negative just set to ~0 > 0
#      # t['charge_try'] = min(t['charge_threshold'], 0.9 * c1)
#      # t['charge_threshold'] should never be negtive (?)
#      t['charge_try'] = min(abs(t['charge_threshold']), 0.9 * c1)
#      #if (t['charge_try'] < 0): 
#      #  t['charge_try'] = 1e-30
#      #  print('Warning!! charge being set to very small')
#
#  t['charge_old'] = t['bunch_charge']
#  t['growth_rate_old'] = t['growth_rate']

#  if (t['charge1'] > 0):  # Have already found a stable and an unstable point  ( charge1 starts as -1 ) 
#    print ('charge0 stable, charge1 unstable')
#    t['bunch_charge'] = (t['charge0'] + t['charge1'])/2  # Try the mid-point between charge0 and charge1
#  else:  # Still searching for an unstable point
#    print ('Have NOT found an unstable upper point')
#    t['bunch_charge'] = t['bunch_charge']*2    # Try double charge
  
  # Set new current value to try
#  if (t['charge_try'] >= 0): # Better approx. found using interpolation
    #t['bunch_charge'] = t['charge_try']
#    if (t['charge1'] > 0):   # Have both stable and unstable points
#      t['bunch_charge'] = (t['charge0'] + t['charge1']) / 2
#    else:  # Still searching for an unstable point
#      t['bunch_charge'] = t['bunch_charge'] * 2

#  else:
  if (t['charge1'] > 0):   # Unstable point found
      t['bunch_charge'] = (t['charge0'] + t['charge1']) / 2
  else:  # Still searching for an unstable point
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
  bt.write('/')
  bt.close()


#==========================================================
def prepare_lat (py_par, lat_file):
## Create "temp_lat.lat" with command lines to call the lattice file 
##########################################################
  f_lat = open(os.path.join(py_par['temp_dir'],'temp_lat.lat'), 'w')
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
  run_bbu( py_par['threshold_start_curr'], py_par, 'current' )  

  # Get all cavity names from hom_info.txt
  f = open(os.path.join(py_par['temp_dir'], 'hom_info.txt'),'r')     
  # Make sure hom_dir contains ONLY hom data files ended in "dat"
  hom_files = glob.glob( os.path.join(py_par['hom_dir'], '*dat') ) #hom_dir contains many HOM.dat
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
  if ( mode == 'current' ):
    template_file = open (os.path.join(py_par['temp_dir'], 'bbu_template.init'), 'r' ) 
    temp_file = open (os.path.join(py_par['temp_dir'],'bbu.init'), 'wt')
    # Change the bbu.init's current to the temp_curr
    for line in template_file:
      temp_file.write( line.replace( "temp_curr", str(temp_curr)) )  
    template_file.close()
    temp_file.close()
    print ('Running BBU with current ', str(temp_curr), '(A)')
   
    #print('Running BBU with following parameters : ')
    #temp_bbu_file = open (os.path.join(py_par['temp_dir'],'bbu.init'), 'r')
    #for line in temp_bbu_file:
    #  print(line)
    #temp_bbu_file.close()

    subprocess.call( py_par['exec_path'], shell = True)  # Run bbu 

  if ( mode == 'drscan' or mode == 'phase_scan' or mode == 'phase_xy_scan'):
    print ('Including lat2.lat...')
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
