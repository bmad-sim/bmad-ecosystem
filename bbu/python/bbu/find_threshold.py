#!/usr/bin/env python
import subprocess

#===============================================================
def get_stability(gain, lost_bool):
####################### 
  if (lost_bool or (gain > 1)):
    success = 0
    print ('Beam is unstable at this current')
  else:
    success  = 1
    print ('Beam is stable at this current')

  return success

#===============================================================
def calc_new_charge( t, bool_stable, rel_tol ):
####################### Get threshold using the charge and growth rate for two tracking points
  if (t['growth_rate'] > 0 or not bool_stable):  # Current is unstable at this current (~charge)
    print ('This current was not stable, bumping up charge1')
    t['charge1'] = t['bunch_charge']
    t['growth_rate1'] = t['growth_rate']
  else:
    print ('This current was stable')
    t['charge0'] = t['bunch_charge']
    t['growth_rate0'] = t['growth_rate']

  if (t['charge1'] > 0):  # Have already found a stable and an unstable point
    print ('charge1 unstable, charge0 stable')
    t['bunch_charge'] = (t['charge0'] + t['charge1'])/2  # charge1 is unstable, charge0 stable
  else:  # Still searching for an unstable point
    print ('Havent found an unstable upper point')
    t['bunch_charge'] = t['bunch_charge']*2   

  #---------------------------------------------------
  if (t['charge_old'] > 0): # If we have at least two trackings so far
    if (t['charge0'] > 0 and t['charge1'] > 0):  # Now we can interpolate for the threshold
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
        print ('Charge threshold guess is now ', t['charge_threshold'])

      # Current to use in the next tracking must be significantly different from c0 and c1.
      t['charge_try'] = t['charge_threshold']
      min_delta = max(c0, c1) * rel_tol
      dc0 = t['charge_threshold'] - c0
      dc1 = t['charge_threshold'] - c1
      # If threshold guess is closer to c0...
      if (abs(dc0) < abs(dc1) and abs(dc0) < min_delta):
        c = t['charge_threshold'] + min_delta*(dc0/abs(dc0))
        if (abs(c-c1) < abs(c-c0)): # If moved closer to c1 then just use the average
          t['charge_try'] = (c1 + c0) / 2
        else:
          t['charge_try'] = c
  
      # If threshold guess is closer to c1...
      elif (abs(dc1) < abs(dc0) and abs(dc1) < min_delta):
        c = t['charge_threshold'] + min_delta*(dc1/abs(dc1))
        print('Nudging charge')
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
        t['charge_threshold'] = 2*c0
      else:
        t['charge_threshold'] = (c0 * g1 - c1 * g0) / (g1 - g0)

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
	
    
      # Demand at least a 10% change, but if negative just set to ~0 > 0
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
def make_init ( bbu_params ):
##################### Makes the template bbu init file from user-set params
  bt = open('bbu_template.init','w')
  bt.write('&bbu_params\n')
  for key in bbu_params:
    bt.write('  bbu_param%'+ key + ' = ' +str( bbu_params[key])+'\n' )
  bt.write('/')
  bt.close()

#==========================================================
def  run_bbu ( temp_curr, exec_path, mode ):
##################### Prepare bbu.init and run the bbu program once for a current
  if ( mode == 'current' ):
    print ('Trying current ', str(temp_curr))
    template_file = open ( 'bbu_template.init', 'r' ) 
    temp_file = open  ('bbu.init', 'wt')
    # Change the bbu.init's current to the temp_curr
    for line in template_file:
      temp_file.write( line.replace( "temp_curr", str(temp_curr)) )  
    template_file.close()
    temp_file.close()
    subprocess.call( exec_path, shell = True)  # Run bbu 

  if ( mode == 'drscan' ):
    print ('Including lat2.lat')
    template_file = open ( 'bbu_template.init', 'r' ) 
    temp_file = open  ('bbu.init', 'wt')   
    # In DR scan mode, include secondary lattice file so arclength may be varied
    ftemp = template_file.read()  
    ftemp = ftemp.replace( "temp_curr", str(temp_curr))
    ftemp = ftemp.replace( "lat2_filename = ''", "lat2_filename = 'lat2.lat'" )
    temp_file.write(ftemp)  
    template_file.close()
    temp_file.close()
    subprocess.call( exec_path, shell = True)  # Run bbu
