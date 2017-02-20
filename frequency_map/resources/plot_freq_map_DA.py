#! /usr/bin/python

#######################################################################
### Plot frequency_map results as a dynamic aperture (DA) projection
#######################################################################

import numpy as np
import numpy.ma as ma
import os
import matplotlib
import matplotlib.pyplot as plt
import pylab, math

#######################################################################
### Modify these lines according to your job:

filename   = 'combined.txt' # file containing frequency_map output
makeinputs = 'make_inputs.py' # script used to generate GRID inputs

# Set limits for amplitude scale. Recall A = log(sqrt[dQx^2+dQy^2])
minimum = 1.e-12
maximum = 0.1

# Twiss parameters and emittances at start of lattice: 
betax = 9.253
betay = 2.217
emitx = 30.99e-9
emity = emitx * 0.01

# Tune plane limits for plotting:
Qxmin = 0.50
Qxmax = 1.0
Qymin = 0.50
Qymax = 1.0

# Minimum projected apertures, in beam-sigmas, as projected
# to start of lattice. Only used if plotting these apertures
# onto DA
plotApertures = False
hlim = 34.9  
vlim = 69.65 







########################################################################
########################################################################
##             Begin Plotting Code - Don't Edit Below Here            ##
########################################################################
########################################################################

xscale = 1./np.sqrt(emitx * betax)
yscale = 1./np.sqrt(emity * betay)

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
        elif fileline[0] == 'y0':
                y0 = float(fileline[2])
	elif fileline[0] == 'y1':
                y1 = float(fileline[2])
	elif fileline[0] == 'dy':
                dy = float(fileline[2])
	elif fileline[0] == 'Qz':
                Qs = float(fileline[2])
        elif fileline[0] == 'for':
                break
file.close()

x = []
y = []
z = []
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
	y.append(y0+ix*dx)
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

	thisx = float(fileLine[0])
	thisy = float(fileLine[1])
	thisz = float(fileLine[2])
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
	yidx = int(round(((thisy-y0) / dy),1))

	if (xidx > len(x)-1) or (yidx > len(y)-1 or xidx < 0 or yidx < 0): continue

	Qx0[xidx][yidx] = thisQx0
	Qy0[xidx][yidx] = thisQy0
	Qx[xidx][yidx] = thisQx
	Qy[xidx][yidx] = thisQy
	thisdQ = np.sqrt( (np.abs(thisdQx))**2 + (np.abs(thisdQy))**2 )
	if (np.log10(thisdQ) < np.log10(minimum/10000.)):
		dQ[xidx][yidx] = minimum
		continue
	else:
		dQ[xidx][yidx] = np.log10(thisdQ)
		dQmask[xidx][yidx] = False

### Rescale x and y:
x = [ix*xscale for ix in x]
y = [ix*yscale for ix in y]

xA = np.array(x)
yA = np.array(y)
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

if plotApertures == True:
	xA = ma.masked_outside(xA, -hlim, hlim)
	yA = ma.masked_outside(yA, -vlim, vlim)


#####################
# Plotting:

fig = plt.figure(num=1,figsize=(10.,6.),dpi=100)
fig.subplots_adjust(left=0.1,right=0.98,top=0.93,bottom=0.1)
ax1 = fig.add_subplot(111)
#ax1.set_aspect('equal')
mappable = ax1.pcolor(xA,yA,dQAmasked, vmin=-logmin, vmax=-logmax)
ax1.set_xlim(x0*xscale,x1*xscale)
ax1.set_ylim(y0*yscale,y1*yscale)
ax1.set_xlabel('x/$\sigma_x$')
ax1.set_ylabel('y/$\sigma_y$')

if plotApertures == True:
	ax1.plot((hlim,hlim),(0,vlim), ls='--', c='r',lw=4)
	ax1.plot((-hlim,-hlim),(0,vlim), ls='--', c='r',lw=4)
	ax1.plot((-hlim,hlim),(vlim,vlim), ls='--', c='r',lw=4)
fig.colorbar(mappable, ax=ax1)
ax1.set_title('$\Delta Q$ vs. (x,y)')
plt.savefig('x_y_dQ_DA.png')

###########################
# plotting of \Delta Q vs. (Qx,Qy)

fig2 = plt.figure(num=2,figsize=(10.,6.),dpi=100)
fig2.subplots_adjust(left=0.1,right=0.94,top=0.93,bottom=0.1)
ax2 = fig2.add_subplot(111)
ax2.scatter(Qx0A.reshape(-1),Qy0A.reshape(-1), s=2.5, c=dQAmasked.reshape(-1), marker='s', edgecolors='none')
ax2.set_xlabel('Qx')
ax2.set_ylabel('Qy')
ax2.set_xlim(Qxmin, Qxmax)
ax2.set_ylim(Qymin, Qymax)
ax2.grid()
ax2.set_title('$\Delta Q$ vs. ($Q_x$,$Q_y$)')
plt.savefig('Qx_Qy_dQ_DA.png')
plt.show()
