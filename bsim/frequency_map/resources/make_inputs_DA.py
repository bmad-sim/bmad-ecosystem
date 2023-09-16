#! /usr/bin/python

x0 = -0.03001
x1 =  0.02999
dx =  0.00016

y0 = 0.00001
y1 = 0.00251
dy = 0.00004

e0 =10.e-6
e1 = 10.e-6
de = 10.e-6

fft_turns = 1024
n_turns = 2048

lat_file = "/home/shanksj/chess/lat/chess-u/cu.lat"

###########################################################

ySlice = int((y1-y0)/dy)+1

for ix in range(ySlice+1):
	yStart = y0+ix*dy
	filename = '%(index)05d' %{'index': round(yStart,5)*100000}
	filename = 'y_' + filename + '.in'
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

