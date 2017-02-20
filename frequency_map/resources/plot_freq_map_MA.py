#! /usr/bin/python

#######################################################################
### Plot frequency_map results as a momentum aperture (MA) projection
#######################################################################

import numpy as np
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
import pylab, math
import os

#######################################################################
### Modify these lines according to your job:

filename   = 'combined.txt' # file containing frequency_map output
makeinputs = 'make_inputs_MA.py' # script used to generate GRID inputs

# Set limits for amplitude scale. Recall A = log(sqrt[dQx^2+dQy^2])
minimum = 1.e-12
maximum = 0.1

# Tune plane limits for plotting:
Qxmin = 0.5
Qxmax = 1.0
Qymin = 0.5
Qymax = 1.0







########################################################################
########################################################################
##             Begin Plotting Code - Don't Edit Below Here            ##
########################################################################
########################################################################

xscale = 1.
yscale = 1.

file = open(makeinputs,'r')
for line in file:
        fileline = line.split()
        if len(fileline) == 0: continue
        if   fileline[0] == 'x0':
                x0 = float(fileline[2])
        elif fileline[0] == 'x1':
                x1 = float(fileline[2])
        elif fileline[0] == 'dx':
                dx = float(fileline[2])
        elif fileline[0] == 'e0':
                e0 = float(fileline[2])
        elif fileline[0] == 'e1':
                e1 = float(fileline[2])
        elif fileline[0] == 'de':
                de = float(fileline[2])
        elif fileline[0] == 'for':
                break
file.close()

x = []
y = []
e = []
Qx0 = []
Qy0 = []
Qz0 = []
Qx = []
Qy = []
Qz = []
dQ = []
dQmask = []

steps = int((x1-x0)/(dx))

# populate arrays with float zeroes to start:
for ix in range(steps+2):
	x.append(x0+ix*dx)
	e.append(e0+ix*de)
	Qx0temp = []
	Qy0temp = []
	Qxtemp = []
	Qytemp = []
	dQtemp = []
	dQmasktemp = []
	for jx in range(steps+2):
		Qx0temp.append(minimum)
		Qy0temp.append(minimum)
		Qxtemp.append(minimum)
		Qytemp.append(minimum)
		dQtemp.append(minimum)
		dQmasktemp.append(True)
	Qx0.append(Qx0temp)
	Qy0.append(Qy0temp)
	Qx.append(Qxtemp)
	Qy.append(Qytemp)
	dQ.append(dQtemp)
	dQmask.append(dQmasktemp)

logmin = np.abs(np.log10(minimum))
logmax = np.abs(np.log10(maximum))


tot = -1 # initialize total
file = open(filename, 'r')
for line in file:
	fileLine = line.split()
	if len(fileLine) == 0 or fileLine[0][0] == '!': continue
	
	tot = tot+1

	if any(fileLine[ix] == 'NaN' for ix in range(len(fileLine))):
		continue
	#if ix in range(len(fileline)):
	#	if fileLine[ix] == 'NaN':
	#		fileline[ix] = minimum

	thisx = float(fileLine[0])
	thisy = float(fileLine[1])
	thise = float(fileLine[2])
	thisQx0 = float(fileLine[3])
	thisQy0 = float(fileLine[4])
	thisQz0 = float(fileLine[5])
	thisQx = float(fileLine[6])
	thisQy = float(fileLine[7])
	thisQz = float(fileLine[8])
	thisdQx = float(fileLine[9])
	thisdQy = float(fileLine[10])
	thisdQz = float(fileLine[11])

	xidx = int(round(((thisx-x0) / dx),1))
	eidx = int(round(((thise-e0) / de),1))

	if (xidx > len(Qx0)-1) or (eidx > len(Qx0[0])-1): continue
	
	Qx0[xidx][eidx] = thisQx0
	Qy0[xidx][eidx] = thisQy0
	Qx[xidx][eidx] = thisQx
	Qy[xidx][eidx] = thisQy
	thisdQ = np.sqrt( (np.abs(thisdQx))**2 + (np.abs(thisdQy))**2 )
	if (np.log10(thisdQ) < np.log10(minimum/10000.)):
		dQ[xidx][yidx] = minimum
		continue
	else:
		dQ[xidx][eidx] = np.log10(thisdQ)
		dQmask[xidx][eidx] = False

x = [ix*xscale for ix in x]

xA = np.array(x)
eA = np.array(e)
Qx0A = np.array(Qx0)
Qy0A = np.array(Qy0)
QxA = np.array(Qx)
QyA = np.array(Qy)
dQA = np.array(dQ)

dQmask = np.transpose(dQmask)
dQA = np.transpose(dQA)
Qx0A = np.transpose(Qx0A)
Qy0A = np.transpose(Qy0A)

dQAmasked = ma.masked_array(dQA,dQmask)

#####################
# Plotting:

#plt.clf()
fig = plt.figure(num=1,figsize=(8,6.),dpi=72)
fig.subplots_adjust(left=0.12,right=0.98,top=0.92,bottom=0.10)
ax1 = fig.add_subplot(111)
#ax1.set_aspect('equal')
mappable = ax1.pcolor(xA,eA,dQAmasked, vmin=-logmin, vmax=-logmax)
ax1.set_xlim(x0*xscale,x1*xscale)
ax1.set_ylim(e0,e1)
ax1.set_xlabel('x (m)')
ax1.set_ylabel('$\delta$')
fig.colorbar(mappable, ax=ax1)
ax1.set_title('$\Delta Q$ vs. (x,$\delta$)')
plt.savefig('x_delta_dQ.png')

###########################
# plotting of \Delta Q vs. (Qx,Qy)

fig2 = plt.figure(num=2,figsize=(10.,6.),dpi=100)
fig2.subplots_adjust(left=0.1,right=0.94,top=0.93,bottom=0.1)
ax2 = fig2.add_subplot(111)
ax2.scatter(Qx0A.reshape(-1),Qy0A.reshape(-1), s=8, c=dQAmasked.reshape(-1), marker='s', edgecolors='none')
ax2.set_xlabel('Qx')
ax2.set_ylabel('Qy')
ax2.set_xlim(Qxmin, Qxmax)
ax2.set_ylim(Qymin, Qymax)
ax2.grid()
ax2.set_title('$\Delta Q$ vs. ($Q_x$,$Q_y$)')
plt.savefig('Qx_Qy_dQ_MA.png')
plt.show()

plt.show()
