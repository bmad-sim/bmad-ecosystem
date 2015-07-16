#!/usr/bin/env python

import subprocess
import bbu_params
from bbu_params import t0
import bbu_find_threshold
import bbu_drscan
import os
import math

if bbu_params.variable == 'current': 
  t = bbu_params.t   
  temp_curr = bbu_params.threshold_start_curr
  bbu_find_threshold.run_bbu( temp_curr )     # Run bbu for first tempory current
  exec(open( 'for_py.txt' ).read())
  t['bunch_charge'] = temp_curr * bunch_dt   # Calculate() will vary bunch charge (~current)
  keep_looking = 1
  while ( keep_looking ):   # Nudge stable current very close to higher, unstable current. Threshold is between.
    t['growth_rate_set0'] = growth_rate_set
    t['growth_rate'] = growth_rate
		
    bool_stable = bbu_find_threshold.get_stability(v_gain, lostbool)
		
    bbu_find_threshold.calculate( t, bool_stable, rel_tol )
    temp_curr = t['bunch_charge'] / bunch_dt    # bunch charge updated in calculate()
    bbu_find_threshold.run_bbu( temp_curr )    # Try new guess for the current
    exec(open( 'for_py.txt' ).read()) 
    t['growth_rate_set1'] = growth_rate_set
    #t['growth_rate'] = growth_rate
    if (abs(t['charge1'] - t['charge0']) < (t['charge1']*rel_tol)): 
      keep_looking = 0
      print('= = = = = = = = =')
      print('	CALCULATED THRESHOLD IS ', temp_curr, ' (A)')


if bbu_params.variable == 'drscan': 
  if ( os.path.exists("thresh_v_trotb.txt") ): os.remove("thresh_v_trotb.txt")
  my_file = open('thresh_v_trotb.txt','a')
  step_size = (bbu_params.end_dr_arctime-bbu_params.start_dr_arctime)/(bbu_params.ndata_pnts-1)
  t = bbu_params.t0
  for n in range (0,bbu_params.ndata_pnts):
    t.clear()
    t = {'charge0':0,'charge1':-1,'charge_old':-1,'charge_try':-1,'growth_rate0':-1,'growthrate1':-1,'growth_rate_set0':False,'growth_rate_set1':False,'growth_rate':0,'growth_rate_old':0,'bunch_charge':0,'charge_threshold':0}
    temp_arctime = bbu_params.start_dr_arctime + n*(step_size)
    #make lat2 file to vary time to traverse the arc in bbu
    bbu_drscan.setup_dr_scan( temp_arctime ) 
    #run bbu at this arc_time (vary arclength)
    bbu_find_threshold.run_bbu( bbu_params.threshold_start_curr )
    exec(open( 'for_py.txt' ).read())    #'/home/mt728/nfs/linux_lib/bsim/for_py.txt'
    t['bunch_charge'] = bbu_params.threshold_start_curr * bunch_dt   
    keep_looking = 1
    while ( keep_looking ):   # nudge stable current very close to the higher, unstable current
      t['growth_rate_set0'] = growth_rate_set
      t['growth_rate'] = growth_rate
      bool_stable = bbu_find_threshold.get_stability(v_gain, lostbool)
      bbu_find_threshold.calculate( t, bool_stable, rel_tol )
      temp_curr = t['bunch_charge'] / bunch_dt    # bunch charge was just updated in calculate()
      bbu_find_threshold.run_bbu( temp_curr )    # now try new guess for the current
      exec(open( 'for_py.txt' ).read()) 
      t['growth_rate_set1'] = growth_rate_set
      if ( abs(t['charge1'] - t['charge0']) < abs(t['charge1']*rel_tol) ): 
        keep_looking = 0 
        print('= = = = = = = = =')

        #append Threshold _v_ Arctime/tb in txt file
        trotb = temp_arctime / bunch_dt
        my_file.write(str(trotb)+'	'+str(math.log(temp_curr,10))+'	'+str(bunch_dt)+'	'+str(temp_curr)+'\n')
  my_file.close()
    
  #plot after all arctimes are finished
  bbu_drscan.make_dr_plot()


		
