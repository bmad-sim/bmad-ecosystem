#!/usr/bin/env python
import subprocess
import bbu_params


#===============================================================
def get_stability(gain, lost_bool):
####################### 
  if (lost_bool or (gain > 1)):
    success = 0
    print ('Beam is unstable at current')
  else:
    success  = 1
    print ('Beam is stable at current so will check for threshold')

  return success

  if _name_=='_main_':
    return get_stability()

#===============================================================
def calculate( t, bool_stable, rel_tol ):
####################### Get threshold using the charge and growth rate for two tracking points

  if (t['growth_rate'] > 0 or not bool_stable):  # Current is unstable at this current (~charge)
    print ('This current wasnt stable, bumping up charge1')
    t['charge1'] = t['bunch_charge']
    t['growth_rate1'] = t['growth_rate']
  else:
    t['charge0'] = t['bunch_charge']
    t['growth_rate0'] = t['growth_rate']

  if (t['charge1'] > 0):  # Set with a stable and an unstable point
    t['bunch_charge'] = (t['charge0'] + t['charge1'])/2 # charge1 is unstable, charge0 stable
  else:  # Still searching for an unstable point.
    t['bunch_charge'] = t['bunch_charge']*2   

  if (t['charge_old'] > 0): # If we have at least two trackings so far
  #---------------------------------------------------
    if (t['charge0'] > 0 and t['charge1'] > 0):
      print('Both charge0 and charge1 have been set')
      c0 = t['charge0']  
      g0 = t['growth_rate0']
      c1 = t['charge1']
      g1 = t['growth_rate1']

      if (not t['growth_rate_set0'] or not t['growth_rate_set1'] or g0 == g1):
        t['charge_threshold'] = (c0 + c1) / 2
        print ('One or more growth rates werent set, or g0 = g1')
      else:
        t['charge_threshold'] = (c0 * g1 - c1 * g0) / (g1 - g0)
        print ('Charge threshold guess is now ', str(t['charge_threshold']))

      # Current to use in the next tracking must be significantly different from c0 and c1.
      t['charge_try'] = t['charge_threshold']
      min_delta = max(c0, c1) * rel_tol
      dc0 = t['charge_threshold'] - c0
      dc1 = t['charge_threshold'] - c1
      # If threshold guess is closer to c0...
      if (abs(dc0) < abs(dc1) and abs(dc0) < min_delta):
        c = t['charge_threshold'] + min_delta*(dc0/abs(dc0))
        print('Nudging charge try')
        if (abs(c-c1) < abs(c-c0)): # If moved closer to c1 then just use the average
          t['charge_try'] = (c1 + c0) / 2
        else:
          t['charge_try'] = c
  
      # If threshold guess is closer to c1...
      elif (dc1 < dc0 and dc1 < min_delta):
        c = t['charge_threshold'] + min_delta*(dc1/abs(dc1))
        print('Nudging charge try')
        if (abs(c-c0) < abs(c-c1)):   #If moved closer to c0 then just use the average
          t['charge_try'] = (c1 + c0) / 2
        else:
          t['charge_try'] = c
  
    # If have stable c0, c1 but no unstable point	
    elif (t['charge1'] < 0):
      c0 = t['charge0']  
      g0 = t['growth_rate0']
      c1 = t['charge_old']  
      g1 = t['growth_rate_old']   
      if (not t['growth_rate_set0'] or not t['growth_rate_set1'] or g0 == g1): 
        charge_threshold = 2*c0 #(c0 + c1) / 2
      else:
        charge_threshold = (c0 * g1 - c1 * g0) / (g1 - g0)

      # If the threshold is less than charge0 then there must be a lot of noise so
      # assume we are far from the threshold and increase the current by a factor of 10
      # In any case, demand at least a 10% change
      if (t['charge_threshold'] < c0):
        t['charge_try'] = 10 * c0
      else:
        t['charge_try'] = max(t['charge_threshold'], 1.1 * c0)
  
    # If both charges correspond to unstable points, no stable point
    else:      
      c0 = t['charge_old']
      g0 = t['growth_rate_old']
      c1 = t['charge1']
      g1 = t['growth_rate1']
    
      if (not t['growth_rate_set0'] or not t['growth_rate_set1'] or g0 == g1): 
        t['charge_threshold'] = c0 / 2
      else:
        t['charge_threshold'] = (c0 * g1 - c1 * g0) / (g1 - g0)
	
    
      # Demand at least a 10% change but if negative just set to ~0 > 0
      t['charge_try'] = min(t['charge_threshold'], 0.9 * c1)
      if (t['charge_try'] < 0): t['charge_try'] = 1e-30

  t['charge_old'] = t['bunch_charge']
  t['growth_rate_old'] = t['growth_rate']

  # Set new current value to try
  if (t['charge_try'] >= 0): 
    t['bunch_charge'] = t['charge_try']

  else:
    if (t['charge1'] > 0):   # Have both stable and unstable points
      t['bunch_charge'] = (t['charge0'] + t['charge1']) / 2
    else:  # Still searching for an unstable point
      t['bunch_charge'] = t['bunch_charge'] * 2



#==========================================================
def  run_bbu ( temp_curr ):
##################### Run the bbu program once for a current
  if (bbu_params.variable == "current"):
    print ('Trying current ', str(temp_curr))
    template_file = open (bbu_params.template_file, 'r' ) 
    temp_file = open  ('bbu.init', 'wt')
    # Change the bbu input current to temp_curr
    for line in template_file:
      temp_file.write( line.replace( "temp_curr", str(temp_curr)) )  
    template_file.close()
    temp_file.close()
    subprocess.call(bbu_params.exec_path, shell = True)  # Run bbu 

  if (bbu_params.variable == "drscan"):
    print ('Including lat2.lat')
    template_file = open (bbu_params.template_file, 'r' )
    temp_file = open  ('bbu.init', 'wt')   
    # In current loop, change the parameter in the init file to new current
    # In DR scan mode, include secondary lattice file so arclength may be varied
    ftemp = template_file.read()  
    ftemp = ftemp.replace( "temp_curr", str(temp_curr))
    ftemp = ftemp.replace( "lat2_filename = ''", "lat2_filename = 'lat2.lat'" )
    temp_file.write(ftemp)  
    template_file.close()
    temp_file.close()
    subprocess.call(bbu_params.exec_path, shell = True)  # Run bbu



