#!/usr/bin/env python

exec_path = '/home/mt728/nfs/linux_lib/production/bin/bbu'
template_file = 'bbu_template.init' # Location of template_init file

ndata_pnts = 65 # Number of steps in current
mode = 'current' # Other mode(s): 'current' for a single threshold; 'drscan' makes plot for PRSTAB 7 (2004) Fig. 3.
threshold_start_curr = 2.0

# Tip for choosing drscan range:  plot x-axis points = arctime / tb, the tb bunch time usually 7.692308e-10
start_dr_arctime = 4.028*10**-9
end_dr_arctime = 4.725*10**-9
plot_drscan = True 

#stable_orbit_anal = True   # Write out cavity orbit data for SIMULATION_TURNS_MAX turns; No thresholds calculated; DRSCAN disabled

# Dictionary for threshold calculation
t0 = {'charge0':0,'charge1':-1,'charge_old':-1,'charge_try':-1,'growth_rate0':-1,'growthrate1':-1,'growth_rate_set0':False,'growth_rate_set1':False,'growth_rate':0,'growth_rate_old':0,'bunch_charge':0,'charge_threshold':0}
