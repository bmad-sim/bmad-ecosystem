README for frequency_map
   Revised 2017.02.20 J. Shanks

Theoretical overview: See, for example, the following:
   J. Laskar, "Frequency Analysis for Mult-Dimensional Systems, Global Dynamics and Diffusion", Physica D 67 (1993) 257-281
   J. Laskar, "Frequency Map Analysis and Particle Accelerators", PAC2003 Proceedings

Conceptual overview:
   -For each coordinate in x-y-e space, track a particle for n_turn
   -If the particle remains stable, take the FFT of the first half and second half of the 
     turns separately. Initial and final amplitudes and tunes (x-y-e and Qx-Qy-Qz) are 
     recorded in the output file
   -Results are plotted in two of four projections. In any projection, the color scale is A = log(dQx^2 + dQy^2)
   	-A vs. x-y with fixed e; also projected to Qx-Qy space
	-A vs. x-e with fixed y; also projected to Qx-Qy space

Detailed Procedure:
   -For small scans, use the template .params file and run binary locally.

   -For larger scans, use GRID submission:
   	-If scanning x-y with fixed e (dynamic aperture projection), modify the make_inputs_DA.py script. 
	-If scanning x-e with fixed y (momentum aperture projection), modify make_inputs_MA.py. 
	     -These scripts slice a large 2D job into a large number of 1D scans. After editing, run the 
	       script. This generates all input files necessary for the job.

	-For submission to the GRID, two scripts are necessary: 
	     -q.sh - low-level script which is submitted to the GRID for each individual input file
	     -submit.py - iterates over all input files, submitting one instance of q.sh per input
	  Modify submit.py to point to your local copy of q.sh. Modify q.sh to point to local GRID log 
	  directories and executable. 

	-After generating input files and modifying the submit.py and q.sh scripts, run submit.py in the 
	  directory where all input files are located. All inputs ending in ".in" will be submitted 
	  automatically.

	-Check progress of GRID jobs using "qstat"

	-At completion of all GRID jobs, run combine_outputs.py. (Note that outputs may be combined at any 
	  time while jobs are still running. This is useful if examining the progress of a slow job, or if 
	  looking for gross changes in behavior from a previous run)

   -For x-y scans with fixed e, plot using plot_freq_map_DA.py
   -For x-e scans with fixed y, plot using plot_freq_map_MA.py
        -Both of these plotting scripts assume the original make_inputs_(DA/MA).py script is in the 
	  same directory as the output files, in order to set axis limits, etc.
	-Note that the DA projection is scaled to beam-sigmas; as such, this script requires knowledge
	  of Twiss parameters at the beginning of the lattice and the horizontal and vertical emittance.

