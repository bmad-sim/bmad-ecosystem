#!/usr/bin/env python

import subprocess
from bbu import find_threshold, drscan
import os
import math

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

def single_threshold( py_par ):
  if ( os.path.exists(os.path.join(py_par['temp_dir'],"for_py.txt")) ): os.remove(os.path.join(py_par['temp_dir'],"for_py.txt"))    
  # Define dict t for threshold calculation
  t = {'charge0':0,'charge1':-1,'charge_old':-1,'charge_try':-1,'growth_rate0':-1,'growthrate1':-1,'growth_rate_set0':False,'growth_rate_set1':False,'growth_rate':0,'growth_rate_old':0,'bunch_charge':0,'charge_threshold':0}
  temp_curr = py_par['threshold_start_curr']
  find_threshold.run_bbu( temp_curr, py_par, 'current' )   # Run bbu for first tempory current
  d = parse_for_py(os.path.join(py_par['temp_dir'],'for_py.txt'))
  t['bunch_charge'] = temp_curr * d['bunch_dt']   # Calculate() will vary bunch charge (~current)
  keep_looking = 1
  while ( keep_looking ):   # Nudge stable current very close to higher, unstable current. Threshold is between the two.
    t['growth_rate_set0'] = d['growth_rate_set']
    t['growth_rate'] = d['growth_rate']
    bool_stable = find_threshold.get_stability(d['v_gain'], d['lostbool'])
    find_threshold.calc_new_charge( t, bool_stable, d['rel_tol'] )
    temp_curr = t['bunch_charge'] / d['bunch_dt']    # Bunch charge updated in calc_new_charge()
    find_threshold.run_bbu( temp_curr, py_par, 'current' )    # Try tracking new guess for the current
    d.clear()
    d = parse_for_py(os.path.join(py_par['temp_dir'],'for_py.txt'))
    t['growth_rate_set1'] = d['growth_rate_set']
    if (abs(t['charge1'] - t['charge0']) < (t['charge1']*d['rel_tol'])): 
      keep_looking = 0
      histo = open(os.path.join(py_par['temp_dir'],'thresholds.txt'), 'a')
      print('= = = = = = = = =')
      print('	CALCULATED THRESHOLD IS ', temp_curr, ' (A)')
      histo.write(str(temp_curr)+'\n')
      histo.close()
    if ( temp_curr < 10**-15 or temp_curr > 10**5 ): 
      keep_looking = 0
      print('==!!==!!==!!==')
      print('	THRESHOLD DID NOT CONVERGE WITHIN BOUNDS')
      histo = open(os.path.join(py_par['temp_dir'],'thresholds.txt'), 'a')
      histo.write('DID NOT CONVERGE\n')
      histo.close()
  # Clean up working directory when finished
  if ( os.path.exists("bbu_template.init") ): os.remove(os.path.join(py_par['temp_dir'],"bbu_template.init"))    
  if ( os.path.exists(os.path.join(py_par['temp_dir'],"bbu.init")) ): os.remove(os.path.join(py_par['temp_dir'],"bbu.init"))    
  t.clear()



def drscanner( py_par ):
  #keyw = open('volt_v_turn.txt', 'a')
  #keyw.write('NEWCURRENT')
  
  # Define dict t for threshold calculation
  t = {'charge0':0,'charge1':-1,'charge_old':-1,'charge_try':-1,'growth_rate0':-1,'growthrate1':-1,'growth_rate_set0':False,'growth_rate_set1':False,'growth_rate':0,'growth_rate_old':0,'bunch_charge':0,'charge_threshold':0}
  if ( os.path.exists(os.path.join(py_par['temp_dir'],"thresh_v_trotb.txt")) ): os.remove(os.path.join(py_par['temp_dir'],"thresh_v_trotb.txt"))
  my_file = open(os.path.join(py_par['temp_dir'],'thresh_v_trotb.txt'),'a')
  step_size = (py_par['end_dr_arctime']-py_par['start_dr_arctime'])/(py_par['ndata_pnts']-1)
  for n in range (0, py_par['ndata_pnts']):
    t.clear()
    t  = {'charge0':0,'charge1':-1,'charge_old':-1,'charge_try':-1,'growth_rate0':-1,'growthrate1':-1,'growth_rate_set0':False,'growth_rate_set1':False,'growth_rate':0,'growth_rate_old':0,'bunch_charge':0,'charge_threshold':0}
    temp_arctime = py_par['start_dr_arctime'] + n*(step_size)
    # Make lat2 file to vary the arctime (~vary arclength)
    drscan.setup_drscan( temp_arctime, py_par ) 
    # Run bbu at this arc_time 
    find_threshold.run_bbu( py_par['threshold_start_curr'], py_par, 'drscan' )
    d = parse_for_py(os.path.join(py_par['temp_dir'],'for_py.txt'))
    t['bunch_charge'] = py_par['threshold_start_curr'] * d['bunch_dt']   
    keep_looking = 1
    while ( keep_looking ):   # Nudge stable current very close to the higher, unstable current
      t['growth_rate_set0'] = d['growth_rate_set']
      t['growth_rate'] = d['growth_rate']
      bool_stable = find_threshold.get_stability(d['v_gain'], d['lostbool'])
      find_threshold.calc_new_charge( t, bool_stable, d['rel_tol'] )
      temp_curr = t['bunch_charge'] / d['bunch_dt']    # Bunch charge was just updated in calc_new_charge()
      find_threshold.run_bbu( temp_curr, py_par, 'drscan' )    # Now try new guess for the current
      d.clear()
      d = parse_for_py(os.path.join(py_par['temp_dir'],'for_py.txt'))
      t['growth_rate_set1'] = d['growth_rate_set']
      if ( abs(t['charge1'] - t['charge0']) < abs(t['charge1']*d['rel_tol']) ): 
        keep_looking = 0 
        print('= = = = = = = = =')
        print('	STEP IN DRSCAN COMPLETE')
        # Append threshold _v_ arctime/tb in txt file for plotting
        trotb = temp_arctime / d['bunch_dt']
        my_file.write(str(trotb)+'	'+str(math.log(temp_curr,10))+'	'+str(temp_curr)+'\n')
        print('JUST WROTE TO thresh_v_trotb.txt')
    d.clear()
  my_file.close()
  # Remove temporary secondary lattice
  # Clean up working directory when finished
  if ( os.path.exists(os.path.join(py_par['temp_dir'],"lat2.lat")) ): os.remove(os.path.join(py_par['temp_dir'],"lat2.lat"))    
  if ( os.path.exists(os.path.join(py_par['temp_dir'],"bbu_template.init")) ): os.remove(os.path.join(py_par['temp_dir'],"bbu_template.init"))    
  if ( os.path.exists(os.path.join(py_par['temp_dir'],"bbu.init")) ): os.remove(os.path.join(py_par['temp_dir'],"bbu.init"))    
  # Plot after all arctimes are finished
  if (py_par['plot_drscan']):
    drscan.make_dr_plot(py_par)  

