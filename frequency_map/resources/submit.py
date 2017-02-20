#! /usr/bin/python

# Jsh, 2009.10.26
# Python script to submit jobs "*.in" to the grid.

import sys, os, re
import getopt, subprocess

if __name__ == '__main__':

	qsh = 'q.sh' # path to q.sh

	subList = []

	path = os.getcwd()
	for file in os.listdir(path):
		if re.match('.*\.in$',file):
			subList.append(str(file))
		#endif
	#endfor
	for input in subList:
		#print input
		os.system('qsub -N ' + input + ' -v inputfile=' + input + ' -v workingdir=' + path + ' ' + qsh)

