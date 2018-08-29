#! /usr/bin/python
import os, sys
##########################################################
### Edit the following for your individual job:

x0 = -0.002090
x1 = 0.002090
dx =  0.1000E-04

y0 = 0.0000E+00 
y1 = 0.000012
dy = 0.000001

e0 = -0.00075 
e1 = 0.00075
de = 0.00001

n_turns   =  4096 # total number of turns to track
fft_turns = 2048 # number of turns to FFT at start and end of tracking

target = '/home/cfsd/lovelace/BMAD/bmad_dist_2018_0807/bsim/frequency_map/resources'
lat_file = os.path.join(target,'../example/10GeV.bmad')

###########################################################
########          DO NOT EDIT BELOW HERE           ########
###########################################################

ySlice = int((y1-y0)/dy)+1

for ix in range(ySlice+1):
	yStart = y0+ix*dy
	filename = '%(index)05d' %{'index': round(yStart,5)*100000}
	filename = os.path.join(target,'y_' + filename + '.in')
	#print filename
	outFile = open(filename,'w')
	outFile.write('&parameters\n')
	outFile.write('aperture_limits = .false.\n')
	outFile.write('x0 = '+str(x0)+'\n')
	outFile.write('y0 = '+str(yStart)+'\n') 
	outFile.write('e0 = '+str(e0)+'\n')
	outFile.write('x1 = '+str(x1)+'\n')
	outFile.write('y1 = '+str(yStart)+'\n') 
	outFile.write('e1 = '+str(e1)+'\n')
	outFile.write('dx = '+str(dx)+'\n')
	outFile.write('dy = '+str(dy)+'\n')
	outFile.write('de = '+str(de)+'\n')
	outFile.write('lat_file = "'+lat_file+'"\n')
	outFile.write('out_file_prefix = ' + filename + '\n')
	outFile.write('n_turn = '+str(n_turns)+'\n')
	outFile.write('fft_turns = '+str(fft_turns)+'\n')
	outFile.write('/\n')
	
	outFile.close()

