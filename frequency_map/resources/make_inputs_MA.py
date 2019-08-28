#! /usr/bin/python

x0 = -0.03001
x1 =  0.02999
dx =  0.00016

y0 = 10.e-6
y1 = 10.e-6
dy = 10.e-6

e0 = -0.05001 # check +/-5% aperture
e1 =  0.04999
de =  0.001

fft_turns = 1024
n_turns = 2048

lat_file = "/home/shanksj/chess/lat/chess-u/cu.lat"

#########################################################

sigfig = 1
while True:
	if (int(de*10**(sigfig)) == 10*int(de*10**(sigfig-1)) and int(de*10**(sigfig)) != 0):
		break
	sigfig = sigfig+1
	

eSlice = int((e1-e0)/de+1)

print(eSlice, (e1-e0), de)

for ix in range(eSlice+1):
	eStart = e0+ix*de
	filename = '%(index)04d' %{'index': round(eStart,sigfig)*10000}
	filename = 'e_' + filename + '.in'
	#print filename
	outFile = open(filename,'w')
	outFile.write('&parameters\n')
	outFile.write('aperture_limits = .false.\n')
	outFile.write('x0 = '+str(x0)+'\n')
	outFile.write('y0 = '+str(y0)+'\n') 
	outFile.write('e0 = '+str(eStart)+'\n')
	outFile.write('x1 = '+str(x1)+'\n')
	outFile.write('y1 = '+str(y1)+'\n') 
	outFile.write('e1 = '+str(eStart)+'\n')
	outFile.write('dx = '+str(dx)+'\n')
	outFile.write('dy = '+str(dy)+'\n')
	outFile.write('de = '+str(de)+'\n')
	outFile.write('lat_file = "'+lat_file+'"\n')
	outFile.write('out_file_prefix = ' + filename + '\n')
	outFile.write('n_turn = '+str(n_turns)+'\n')
	outFile.write('fft_turns = '+str(fft_turns)+'\n')
	outFile.write('/\n')
	
	outFile.close()

