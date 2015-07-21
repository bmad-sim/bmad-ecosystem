BBU PROGRAM __ 2015

GENERAL USAGE NOTES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The bbu python package can be run from any user directory after it has been added to their python path (see SET UP).
To run the bbu program, changes need only be made to the file bbu/examples/examp_bbu.py  

BBU RUN MODES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'current':  Finds the threshold for a single current, which is set in scripts/bbu_params.py

'drscan':  Scans over a range of arclengths and creates a plot of current threshold v. arctime/bunchtime  
           Recreates the plot for PRSTAB 7 (2004) Fig. 3. 

SET UP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Get the bbu package:
   Pull bbu/python/bbu/ from the svn and into your python directory (e.g. /home/USERNAME/nfs/linux_lib/python/)
2. Add the bbu package to your python path:
   Add this line in the Python section of your ~.bashrc, with appropriate username,
     export PYTHONPATH=${PYTHONPATH}:/home/USERNAME/nfs/linux_lib/python/bbu
3. Set your parameters in the examples/example_bbu.py:
   --Params for bbu:  Set lattice file path, hybridize T/F, max turns, n repetitions, etc.
   --Params for python: Set run modes, start / end currents, executable file paths, etc.
4. Run examples/examp_bbu.py with python3

EXAMPLES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The examples directory includes everything needed for a test run of bbu  
~~examples/examp_bbu.py : The only file a user needs to edit and run for bbu
--examples/oneturn_lat.bmad : Main lattice file which calls the simple lattice and adds the HOM data files (below)
-- examples/oneturn_lat_examp.bmad :  A single-turn, single-cavity erl   
-- examples/oneturn_assignHOMs.bmad : This includes examples/vhoms_cavity_1.dat in the single SRF cavity of the example lattice.
                                      Any number of hom data files can be associated with any cavity of a lattice file.
                                      vhoms_cavity_1.dat has a single wakefield with specs the same as PRSTAB 7 (2004) Fig. 3.

FILES LIST
============== These files are necessary
In bbu:
bbu_program.f90  ---  Main bmad / tracking code
bbu_track_mod.f90  ---  BBU tracking and analysis modules

In python package:
bbu.py  ---  Main python wrapper script
find_threshold.py  ---  The meat of the threshold calculation; create bbu.init file / runs bbu_program.f90
drscan.py  ---  Prepares drscan files and plots. lat2.lat (temporary) and thresh_v_trotb.txt will be created.


