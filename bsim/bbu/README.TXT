BBU program 2017 (edited by William Lou)
----------------------------------------

The main BBU directory is bsim/bbu/ within the Git repo at <https://github.com/bmad-sim/bmad-ecosystem>
The BBU code and related files are contained in this directory, except for bsim/code/bbu_track_mod.f90. 
The compiled BBU Fortran code is at production/bin/bbu by default.
The latest description of BBU program (in Latex) is under the doc directory.


===========================================================================
Set Up:
-------

1) Edit the environment variable $PYTHONPATH. 
   Add this line to your .bashrc file (with YOUR username):
   "export PYTHONPATH=${PYTHONPATH}:<path-to-Bmad-release>/bsim/bbu/python"
   where <path-to-Bmad-release> is the path to the Bmad release directory.

2) In test_run.py, properly set the BBU AND Python parameters.
   The user is expected to understand what most parameters do. 
   MAKE SURE THE CORRECT LATTICE FILE LOCATION is specified in bbu/test_run.py


==============================================================================
Run BBU:
--------

Note: the output files will be stored in the working directory at which python3 is called.   


1) DR_SCAN mode:  
       -------------------------------------------------------
       Run this command under $DIST_BASE_DIR/bsim/bbu/examples :
       
       python $DIST_BASE_DIR/bsim/bbu/test_run.py 
       -------------------------------------------------------
       Obtain (Ith v.s. tr/tb) for a lattice with one cavity and varying recirculation arclength.
       Can produce the plot for PRSTAB 7 (2004) Fig. 3. (with one dipole HOM assigned, and one-time recirculation)
       This is the simplest mode to test if the BBU program works.
       Output file: thresh_v_trotb.txt

2) THRESHOLD (CURRENT) mode:      
       -------------------------------------------------------------------------
       Command:  python $DIST_BASE_DIR/bsim/bbu/test_run.py  N  fn  output_dir
       -------------------------------------------------------------------------
       For a FIXED lattice with many cavities assigned with (random) HOMs,
       compute the threshold current for N times.
       It is recommended to set N = 1 for parallel computation. 
       ( fn is the file number, and set to the job number ($JOB_ID) during grid submission )
       Output file:  bbu_threshold_fn.txt  and   rand_assign_homs_fn.txt 
       rand_assign_homs_fn.txt will NOT be produced if random_HOM is set false

3) PHASE_SCAN mode:
       -------------------------------------------------------
       Command:  python $DIST_BASE_DIR/bsim/bbu/test_run.py  PHASE
       -------------------------------------------------------
       Obtain (Ith v.s. tr/tb) for a lattice with varying horizontal PHASE ADVANCE.
       Phase is changed via a Taylor element in the lattice, named talylorW.
       The correct Twiss parameters must be extracted from the lattice and set in phase_scan.py.
       Output file: thresh_v_phase_PHASE.txt  (produced at where it's called)


4) PHASE_SCAN_XY mode:
       -------------------------------------------------------
       Command:  python $DIST_BASE_DIR/bsim/bbu/test_run.py  PHASE_X  PHASE_Y
       -------------------------------------------------------
       Obtain (Ith v.s. tr/tb) for a lattice with two free phases (with coupled or decoupled optics)
       Phase advances are changed via a Taylor element in the lattice, named talylorW.
       The correct Twiss parameters must be extracted from the lattice and set in phase_scan.py.
       Output file: thresh_v_phase_PHASE.txt  (produced at where it's called)


==============================================================================
Important Files:
----------------

-> BBU Fortran codes (Core Computation):                 
   1)  .../bsim/bbu/bbu_program.f90     ---  Main bmad / tracking code
   2)  .../bsim/code/bbu_track_mod.f90  ---  BBU tracking and analysis modules
                                     
-> Python main user interface  ( Most upper-level program )
(Location: .../bsim/bbu/*)
   1) test_run.py  --- TOP interface file (formerly as run_examp_bbu.py)
                   --- Define all BBU and python parameters ( can be modified by user )
                   --- Create and Clean up the temporary directory

-> Python package (Call BBU Fortran codes; intermediate calculation, file organization):
(Location: .../bsim/bbu/python/bbu/*)
(Make sure .../bsim/bbu/pyhton/ is included in the $PYTHONPATH in ~/.bashrc)

  1) bbu_main.py        ---  Parse BBU result. Calls other python codes.
  2) find_threshold.py  ---  The core of the threshold calculation; 
                        ---  Prepare temporary files to run bbu / runs bbu_program.f90
  3) drscan.py          ---  Prepares drscan files and plots. lat2.lat (temp) and thresh_v_trotb.txt will be created.
  4) phase_scan.py      ---  Prepares drscan files and plots. lat2.lat (temp) and thresh_v_phase_PH.txt will be created.
  
-> Extra stand-alone programs  
  1) cut_HOM.py  ---  Accepts a directory full of wakefiled data files and a maximum number (N) of desired HOMs per file
                 and creates a new user-specified directory of corresponding files cut to include only the most influential N HOMs

  2) collect_thresholds.py  --- Summarizes the calculated threshold currents stored in the "bbu_thresholds_*" output files
                                Usually used after a complete grid submission
                            --- Output: bbu_combined_thresholds.txt

-> Local HOM assignment file
   1) assignHOMs.bmad   --- If this file exists in the directory which the user calls test_run.py 
                            ( not necessarily the directory where test_run.py locates ),
                            the HOM assignemnt from the original lattice will be over-written.
                            Unless the user has specified py_par['random_homs'] to be true
