#! /home/shanksj/misc/anaconda/bin/ipython

import os, re, sys

outList = []

path = os.getcwd()
for file in os.listdir(path):
	if re.match('.*\.fm',file):
		outList.append(str(file))
	#endif
#endfor

if len(outList) == 0:
	sys.exit("No *.fm files present! Preventing overwrite of combined.txt.")

outList.sort()
combinedFile = open('combined.txt', 'w')
combinedFile.write('! column labels: \n')
combinedFile.write('! x   y   z   Qx_0   Qy_0   Qz_0   Qx_1   Qy_1   Qz_1    delta(Qx)     delta(Qy)    delta(Qz)   \n')
combinedFile.write('!\n!\n!\n!\n!\n!\n!\n')

for ix in outList:
	outfile = open(ix,'r')
	for line in outfile:
		if len(line.strip()) == 0: continue
		combinedFile.write(line)
	outfile.close()


combinedFile.close()
