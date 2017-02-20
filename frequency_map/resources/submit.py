#! /usr/bin/python

# Jsh, 2009.10.26
# Python script to submit jobs "*.in" to the grid.

import sys, os, re
import getopt, subprocess

if __name__ == '__main__':

	subList = []

#	opts, params = getopt.getopt(sys.argv[1:],'') 
	
	path = os.getcwd()
	for file in os.listdir(path):
		if re.match('.*\.in$',file):
			subList.append(str(file))
		#endif
	#endfor
	for input in subList:
		#print input
		os.system('qsub -N ' + input + ' -P ilcsim -v inputfile=' + input + ' -v workingdir=' + path + ' /home/shanksj/chess/freq_map_jobs/resources/q.sh')

