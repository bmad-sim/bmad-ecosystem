BBU PROGRAM __ 2015

GENERAL USAGE NOTES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The bbu python package can be run from any user directory after it has been added to their python path (see SET UP).
To run the bbu program, changes need only be made to the file bbu/examples/examp_bbu.py  

BBU RUN MODES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'current':  Finds the threshold for a single current, which is specified in the py_par section of the run script.

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
~~ examples/run_examp_bbu.py : The only file a user needs to edit and run for bbu
-- examples/oneturn_lat.bmad : Main lattice file which calls the simple lattice and adds the HOM data files (below)
-- examples/oneturn_lat_examp.bmad :  A single-HOM, single-cavity erl   
-- examples/oneturn_assignHOMs.bmad : This includes examples/vhoms_cavity_1.dat in the single SRF cavity of the example lattice.
                                      Any number of hom data files can be associated with any cavity of a lattice file.
                                      vhoms_cavity_1.dat has a single wakefield with specs the same as PRSTAB 7 (2004) Fig. 3. 
-- This directory also includes the mlc lattice, which contains 6 SRF cavities. The main MLC lattice file is mlc.lat, 
   and the lattice file la.lat.bmad contains the specs. 

Make a directory of cut HOM files
  Many-wake files can slow down simulation. cut_HOMs.py sorts all HOMs in a data file from "worst" to "best", and retains only the "worst" HOMs in output.
  cut_HOM.py accepts a directory of HOM data files, the location for the new cut directory, and the number of HOMs retained per file.
  Command: HOM_cutter( #HOMSTORETAIN, '/home/USER/DIRECTORY OF HOM .DATS', '/home/USER/DESIRED CUT DIRECTORY PATH/')

To randomly assign HOMs to the lattice's cavities, given a directory of hom data files
  Include the line >>> find_threshold.make_assign( py_par,  bbu_par['lat_filename'], 'yes_names' ) <<< after the make_init() function call in the run file.

To perform a DRscan
  The only file a user needs to edit and run for a dr scan. 
  Set the python parameter ``random_homs" to False and set the bbu parameter ``lattice\_filename" to the full path of \texttt{oneturn_lat.bmad} which is in the examples directory.  
  Next, if one runs the script from the command line with no arguments, it will attempt a drscan.

To submit to the grid
  In the examples directory, mlc.sge is the basic submit script for bbu. For bulk submissions, grid_submitter.bash is also included to recursively call mlc.sge. 
  Obviously, if you are not using the mlc lattice, this can be named lattice.sge 

To collect the output files of many jobs
  In the examples directory, collect_thresholds.py takes a directory of threshold output files ('bbu_threshold_#ID.txt') from the bbu package, 
  and combines all of them into a single file of thresholds. It also outputs the number of threshold calculations that did not converge to a threshold
  and the total number of thresholds in the final file.

FILES LIST
============== These files are necessary
In bbu:
bbu_program.f90  ---  Main bmad / tracking code
bbu_track_mod.f90  ---  BBU tracking and analysis modules

In python package:
bbu.py  ---  Main python wrapper script
find_threshold.py  ---  The meat of the threshold calculation; create bbu.init file / runs bbu_program.f90
drscan.py  ---  Prepares drscan files and plots. lat2.lat (temporary) and thresh_v_trotb.txt will be created.
cut_HOM.py  ---  Accepts a directory full of wakefiled data files and a maximum number (N) of desired HOMs per file
                 and creates a new user-specified directory of corresponding files cut to include only the most influential N HOMs

