#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os

print "This python script was developed on Python 2.7.10 :: Anaconda 2.3.0 (64-bit)"
print "Which has matplotlib 1.4.3 and numpy 1.9.2"
print "In particular, if you get an Unknown projection '3d' error."

majorFormatter = FormatStrFormatter('%.3f')

data_file = "knob.grid"
lattice_file = "dc12c.bmad"

def poly33(x,y):
	return np.array([x*0+1, y, x, y*y, x*y, x*x, x*y*y, x*x*y, y*y*y, x*x*x]).T

def poly33c(x,y,c):
	return c[0] + c[1]*y + c[2]*x + c[3]*y*y + c[4]*x*y + c[5]*x*x + c[6]*x*y*y + c[7]*x*x*y + c[8]*y*y*y + c[9]*x*x*x

structured_data = np.genfromtxt(data_file, skip_header=0, names=True)
data = structured_data.view(np.float64).reshape(structured_data.shape + (-1,))
# dQx data[:,0]
# dQy data[:,1]
# Var1 data[:,2]
# Var2 data[:,3]
# etc ...

nVars = data.shape[1] - 2
nx = (np.unique(data[:,0].round(decimals=6))).shape[0]
ny = (np.unique(data[:,1].round(decimals=6))).shape[0]

Xm =  np.reshape(data[:,0], (-1,nx))
Ym =  np.reshape(data[:,1], (-1,ny))

coeff = []
for i in range(0,nVars):
	coeff.append(np.linalg.lstsq(poly33(data[:,0],data[:,1]),data[:,i+2])[0])
	fig = plt.figure(figsize=(18,15))
	ax = fig.gca(projection='3d')
	ax.plot_wireframe(Xm, Ym, poly33c(Xm,Ym,coeff[i]), rstride=1, cstride=1)
	ax.scatter(Xm, Ym, data[:,i+2],color='k')
	plt.title('Points are data from matching program.  Wireframe is polynomial fit.')
	ax.set_xlabel('delta Qx')
	ax.set_ylabel('delta Qy')
	ax.set_zlabel('delta '+structured_data.dtype.names[i+2]+' strength')
	plt.show()

format_str = '%0.3e+b*%0.3e+a*%0.3e+b^2*%0.3e+a*b*%0.3e+ &\n a^2*%0.3e+a*b^2*%0.3e+b*a^2*%0.3e+b^3*%0.3e+a^3*%0.3e'

f=open('knob.bmad','w')
f.write('g1 : group = { &\n')
for i in range(0,nVars):
	if(i < nVars-1):
		f.write(structured_data.dtype.names[i+2]+'[k1]:'+format_str%tuple(coeff[i])+', &\n')
	else:
		f.write(structured_data.dtype.names[i+2]+'[k1]:'+format_str%tuple(coeff[i])+'}, var = {a,b}\n')
f.close()

testf=open('test.results','w')
print "Testing knob ..."
npts = 5
a_min = data[0,0]
a_max = data[-1,0]
b_min = data[0,1]
b_max = data[-1,1]
print '   knob a range is '+str(a_min)+' to '+str(a_max)+' ...'
print '   knob b range is '+str(b_min)+' to '+str(b_max)+' ...'
for i in range(npts):
 	a = a_min + (a_max-a_min)/(npts-1)*i
	for j in range(npts):
 		b = b_min + (b_max-b_min)/(npts-1)*j
 		f=open('lat.bmad','w')
	 	f.write('call, file = '+lattice_file+' \n')
	 	f.write('call, file = knob.bmad \n')
	 	f.write('g1[a] = '+str(a)+'\n')
	 	f.write('g1[b] = '+str(b)+'\n')
	 	f.close()
 	 	stdout = os.popen("/afs/psi.ch/user/e/ehrlichman_m/bbin/tunes lat.bmad").read()
 	 	testf.write(str(a)+'  '+str(b)+'   '+stdout.split()[-2:-1][0]+'   '+stdout.split()[-1:][0]+'\n')
testf.close()

testdata = np.genfromtxt('test.results')

fig = plt.figure(figsize=(18,15))
ax = fig.add_subplot(111)
ax.xaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_major_formatter(majorFormatter)
plt.plot(testdata[:,2],testdata[:,3],'.')
for abxy in testdata[:,(0,1,2,3)]:
	ax.annotate('(%.3f,%.3f)' % (abxy[0],abxy[1]), xy=(abxy[2],abxy[3]), textcoords='data')
plt.title('x-y axes are tunes from lattice file.  Point labels are knob settings for a and b.')
plt.xlabel('Qx')
plt.ylabel('Qy')
plt.show()



