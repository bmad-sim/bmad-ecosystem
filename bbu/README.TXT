BBU PROGRAM __ 2015

GENERAL USAGE NOTES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To run the bbu program, changes need only be made to bbu/scripts/bbu_params.py and bbu/bbu_template.init.  
BBU param settings can not be changed in the bbu.init file, because that is overwritten by the bbu_template.init file.  
BBU is set up to run from the bbu directory, 'python3 scripts/bbu.py'

BBU RUN MODES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'current':  Finds the threshold for a single current, which is set in scripts/bbu_params.py

'drscan':  Scans over a range of arclengths and creates a plot of threshold v. arctime/bunchtime  
           Recreates the plot for PRSTAB 7 (2004) Fig. 3. 

EXAMPLES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bbu/examples contains a simple lattice file 'examples/oneturn_lat.bmad', which calls ::  
-- examples/oneturn_lat_examp.bmad :  A single-turn, single-cavity erl.   
-- examples/oneturn_assignHOMs.bmad : This includes examples/vhoms_cavity_1.dat in the single SRF cavity of the example lattice.
                                      Any number of hom data files can be associated with any cavity of a lattice file.
                                      vhoms_cavity_1.dat has a single wakefield with specs the same as PRSTAB 7 (2004) Fig. 3.

FILES LIST
============== These files are necessary
bbu/bbu_program.f90  ---  Main bmad / tracking code
bbu/bbu_track_mod.f90  ---  BBU tracking and analysis modules
bbu/bbu_template.init  ---  Param file for .f90.  Set lattice file path, hybridize T/F, max turns, n repetitions, etc. 

bbu/scripts/bbu.py  ---  Main python wrapper script
bbu/scripts/bbu_params.py  ---  Param file for .py.  Set run modes, start / end currents, executable file paths, etc.
bbu/scripts/bbu_find_threshold.py  ---  The meat of the threshold calculation; runs bbu_program.f90 
bbu/scripts/bbu_drscan.py  ---  Prepares drscan files and plots. lat2.lat, scripts/thresh_v_trotb.txt will be created.


