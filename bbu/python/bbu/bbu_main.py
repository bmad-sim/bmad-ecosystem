#!/usr/bin/env python

import subprocess
from bbu import find_threshold, drscan
import os
import math

T0  = {'charge0':0,'charge1':-1,'charge_old':-1,'charge_try':-1,'growth_rate0':-1,'growthrate1':-1,'growth_rate_set0':False,'growth_rate_set1':False,'growth_rate':0,'growth_rate_old':0,'bunch_charge':0,'charge_threshold':0}

# Put tracking info from last run of bbu (which was put in for_py.txt) into a dict d
def parse_for_py(filename='for_py.txt'):
  d = {}
  f =  open(filename,'r')
  for line in f:
    line = line.strip()
    key, value = line.replace(' ','').split('=')
    if len(key) > 0:
      if (key == 'lostbool' or key == 'growth_rate_set'):
        if (value == 'True'): val = True
        elif (value == 'False'): val = False
        d[key] = val
      else:
        d[key] = float(value)
  f.close()
  return d

def single_threshold( params ):
  # Define dict t for threshold calculation
  t = T0
  temp_curr = params['threshold_start_curr']
  find_threshold.run_bbu( temp_curr, params['exec_path'], 'current' )   # Run bbu for first tempory current
  d = parse_for_py('for_py.txt')
  t['bunch_charge'] = temp_curr * d['bunch_dt']   # Calculate() will vary bunch charge (~current)
  keep_looking = 1
  while ( keep_looking ):   # Nudge stable current very close to higher, unstable current. Threshold is between the two.
    t['growth_rate_set0'] = d['growth_rate_set']
    t['growth_rate'] = d['growth_rate']
    bool_stable = find_threshold.get_stability(d['v_gain'], d['lostbool'])
    find_threshold.calc_new_charge( t, bool_stable, d['rel_tol'] )
    temp_curr = t['bunch_charge'] / d['bunch_dt']    # Bunch charge updated in calc_new_charge()
    find_threshold.run_bbu( temp_curr, params['exec_path'], 'current' )    # Try tracking new guess for the current
    d.clear()
    d = parse_for_py('for_py.txt')
    t['growth_rate_set1'] = d['growth_rate_set']
    if (abs(t['charge1'] - t['charge0']) < (t['charge1']*d['rel_tol'])): 
      keep_looking = 0
      print('= = = = = = = = =')
      print('	CALCULATED THRESHOLD IS ', temp_curr, ' (A)')
  # Clean up working directory when finished
  if ( os.path.exists("bbu_template.init") ): os.remove("bbu_template.init")    
  if ( os.path.exists("bbu.init") ): os.remove("bbu.init")    

def drscanner( params ):
  # Define dict t for threshold calculation
  t = T0
  if ( os.path.exists("thresh_v_trotb.txt") ): os.remove("thresh_v_trotb.txt")
  my_file = open('thresh_v_trotb.txt','a')
  step_size = (params['end_dr_arctime']-params['start_dr_arctime'])/(params['ndata_pnts']-1)
  for n in range (0, params['ndata_pnts']):
    t.clear()
    t  = {'charge0':0,'charge1':-1,'charge_old':-1,'charge_try':-1,'growth_rate0':-1,'growthrate1':-1,'growth_rate_set0':False,'growth_rate_set1':False,'growth_rate':0,'growth_rate_old':0,'bunch_charge':0,'charge_threshold':0}
    temp_arctime = params['start_dr_arctime'] + n*(step_size)
    # Make lat2 file to vary the arctime (~vary arclength)
    drscan.setup_drscan( temp_arctime ) 
    # Run bbu at this arc_time 
    find_threshold.run_bbu( params['threshold_start_curr'], params['exec_path'], 'drscan' )
    d = parse_for_py('for_py.txt')
    t['bunch_charge'] = params['threshold_start_curr'] * d['bunch_dt']   
    keep_looking = 1
    while ( keep_looking ):   # Nudge stable current very close to the higher, unstable current
      t['growth_rate_set0'] = d['growth_rate_set']
      t['growth_rate'] = d['growth_rate']
      bool_stable = find_threshold.get_stability(d['v_gain'], d['lostbool'])
      find_threshold.calc_new_charge( t, bool_stable, d['rel_tol'] )
      temp_curr = t['bunch_charge'] / d['bunch_dt']    # Bunch charge was just updated in calc_new_charge()
      find_threshold.run_bbu( temp_curr, params['exec_path'], 'drscan' )    # Now try new guess for the current
      d.clear()
      d = parse_for_py('for_py.txt')
      t['growth_rate_set1'] = d['growth_rate_set']
      if ( abs(t['charge1'] - t['charge0']) < abs(t['charge1']*d['rel_tol']) ): 
        keep_looking = 0 
        print('= = = = = = = = =')
        print('	DRSCAN COMPLETE')
        # Append threshold _v_ arctime/tb in txt file for plotting
        trotb = temp_arctime / d['bunch_dt']
        my_file.write(str(trotb)+'	'+str(math.log(temp_curr,10))+'	'+str(temp_curr)+'\n')
    d.clear()
  my_file.close()
  # Remove temporary secondary lattice
  # Clean up working directory when finished
  if ( os.path.exists("lat2.lat") ): os.remove("lat2.lat")    
  if ( os.path.exists("bbu_template.init") ): os.remove("bbu_template.init")    
  if ( os.path.exists("bbu.init") ): os.remove("bbu.init")    
  # Plot after all arctimes are finished
  if (params['plot_drscan']):
    drscan.make_dr_plot()
    		
