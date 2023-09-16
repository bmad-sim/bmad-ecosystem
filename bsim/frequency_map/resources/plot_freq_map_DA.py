#! /home/shanksj/misc/anaconda/bin/ipython

import numpy as np
import numpy.ma as ma
import os
import matplotlib
import matplotlib.pyplot as plt
import pylab, math
import matplotlib.colors as colors

matplotlib.rcParams["savefig.directory"] = os.getcwd()
matplotlib.rcParams['image.cmap'] = 'jet'
#matplotlib.rcParams['image.cmap'] = 'terrain'
#matplotlib.rcParams['image.cmap'] = 'CMRmap'

filename = 'combined.txt'
makeinputs = 'make_inputs.py'

masked = False



dQxcutoff = 2.5e-5
dQycutoff = 2.5e-3
minimum = 1.e-12
maximum = 0.1
plotTitle = ""

## evaluate at und_mid
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


xscale = 1.#/np.sqrt(22.74e-9  * 11.187)
yscale = 1.#/np.sqrt(22.74e-11 *  2.520)

Qxmin = 0.50
Qxmax = 1.0
Qymin = 0.50
Qymax = 1.0

#Qxmin = 0.535
#Qxmax = 0.550
#Qymin = 0.625
#Qymax = 0.6325

#Qxmin = 0.5325
#Qxmax = 0.5575
#Qymin = 0.626
#Qymax = 0.6335

hlim = 26.5 # in beam sigmas
vlim = 28.4 # in beam sigmas


##########################
# For analytic tune plane:
PlotAnalytic = False

pMax = 2
qMax = 2
rMax = 4
nMax = 10
pqrMax = 6
pqMax = 4
##########################

#################################################
## Begin Plotting Code - Don't Edit Below Here ##
#################################################

windXmin = Qxmin
windXmax = Qxmax
windYmin = Qymin
windYmax = Qymax

#def analyticTunePlane(ax):
#	resLines = []
#	for n in range(nMax+1):
#		for r in range(-(rMax),(rMax+1)):
#			for p in range(-(pMax),(pMax+1)):
#				q = 0  #q=0 lines are vertical and handled differently
#				if ((abs(p)+abs(r)) <= pqrMax): #q is zero here
#					if (abs(p) <= pqMax): #q is zero here
#						if (n==0):
#							if (p<0):
#								continue;
#						if (p != 0):
#							x1 = (n - r*Qs)/p
#							ax.plot([x1,x1],[0.0,1.0], linestyle='solid', color='0.0', lw=0.5)
#
#							#Save line in y=mx+b format, with m = inf, b = nan, and x = x-intercept
#							m = float("inf")  #line is vertical
#							b = float("nan")  #no y intercept
#							x = x1
#							resLines.append([m,b,x,"(%1i,%1i,%1i,%1i)"%(p,q,r,n)]) # [slope,y-int,x-int] x-int is used for vertical lines
#				qrange = list(range(-(qMax),(qMax+1)))
#				qrange.remove(0)
#				for q in qrange:
#					if ((abs(p)+abs(q)+abs(r)) <= pqrMax):
#						if ((abs(p)+abs(q)) <= pqMax):
#							if (n==0):
#								if (p<0):
#									continue;
#								if (p==0) & (q<0):
#									continue;
#							x1 = 0.0
#							x2 = 1.0
#							y1 = (n-p*x1-r*Qs)/q
#							y2 = (n-p*x2-r*Qs)/q
#							ax.plot([x1,x2],[y1,y2], linestyle='solid', color='0.0', lw=0.5)
#
#							#Save line in y=mx+b format. third element, x, is used for vertical lines only
#							m = -float(p)/q
#							b = (n-r*Qs)/q
#							resLines.append([m,b,0.0,"(%1i,%1i,%1i,%1i)"%(p,q,r,n)])
#	
#	#Add window borders as lines
#	resLines.append([float('inf'),float('nan'),1.00001*windXmin,'lb'])  #left border
#	resLines.append([float('inf'),float('nan'),0.99999*windXmax,'rb'])  #right border
#	resLines.append([0.0,1.00001*windYmin,0.0,'bb']) #bottom border
#	resLines.append([0.0,0.99999*windYmax,0.0,'tp']) #top border
#
#	#Find point on each line that is furthest from any other lines.
#	#This point is where the label will be placed
#	for i in range(len(resLines)-4):  #Subtract 4 because the last 4 lines are borders, and we do not wish to label borders.
#		[m1,b1,x1,labelText] = resLines[i]
#
#		#print "Processing line [m,b,x,label] = ", m1, b1, x1, labelText
#
#		locs = []  #Will be populated with intersections
#
#		otherLines = list(range(len(resLines)))
#		otherLines.remove(i)
#
#		#Calculate where line i intersects each other line
#		if (m1 != float('inf')):
#			for j in otherLines:
#				[m2,b2,x2,notUsed] = resLines[j]
#				if (m2 != float('inf')):
#					#neither m1 nor m2 is vertical
#					if b1 != b2:
#						if m1 != m2:
#							x0 = (b2-b1)/(m1-m2)
#							y0 = m1*x0+b1
#							#Check if intersection is inside window
#							if (x0<windXmax) & (x0>windXmin):
#								if (y0<windYmax) & (y0>windYmin):
#									#Use y-intercept of line i as origin
#									loc = np.sqrt(x0*x0+(y0-b1)*(y0-b1))
#									locs.append(loc)
#				else:
#					#m1 is not vertical, m2 is vertical
#						x0 = x2
#						y0 = m1*x0+b1
#						if (x0<windXmax) & (x0>windXmin):
#							if (y0<windYmax) & (y0>windYmin):
#								#Use y-intercept of line i as origin
#								loc = np.sqrt(x0*x0+(y0-b1)*(y0-b1))
#								locs.append(loc)
#		else: #m1 is vertical
#			for j in otherLines:
#				[m2,b2,x2,notUsed] = resLines[j]
#				if (m2 != float('inf')):
#					#m1 is vertical, m2 is not vertical
#					x0 = x1
#					y0 = m2*x0+b2
#					if (x0<windXmax) & (x0>windXmin):
#						if (y0<windYmax) & (y0>windYmin):
#							#Use x-intercept of line i as origin
#							loc = np.sqrt((x0-x1)*(x0-x1)+y0*y0)
#							locs.append(loc)
#
#		if (len(locs) > 0):
#			#Determine best location for label as location furthest any other line intersection
#			locs.sort() #Contains locations of intersections, relative to y-intercept or x-intercept
#			dists = []  #Will contain distances between neighboring points
#			for i in range(len(locs)-1):
#				dists.append( locs[i+1] - locs[i] )
#
#			idx = dists.index(max(dists))
#			labelLoc = locs[idx] + dists[idx]/2.0  #Distance along line i from y-intercept or x-intercept to place label
#
#			if (m1 != float('inf')):
#				# line i is not vertical
#				theta = math.atan(m1)
#				xlabel = labelLoc * math.cos(theta)
#				ylabel = labelLoc * math.sin(theta) + b1
#			else:
#				# line i is vertical
#				theta = 3.14159/2.0
#				xlabel = x1
#				ylabel = labelLoc
#
#			thetaDeg = math.degrees(theta)
#			#Transform angle from plot to screen coordinate system
#			coords = pylab.array((xlabel,ylabel))
#			trans_angle = pylab.gca().transData.transform_angles(pylab.array((thetaDeg,)),coords.reshape((1,2)))[0]
#
#			#Plot the label
#			ax.text(xlabel,ylabel,labelText,size=10.0,ha='center',va='bottom',rotation_mode='anchor',rotation=trans_angle)


###########################################
#######      Primary routine      #########
###########################################

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

xsteps = int((x1-x0)/(dx))
ysteps = int((y1-y0)/(dy))

# populate arrays with float zeroes to start:

Qx0    = np.zeros([xsteps,ysteps])
Qy0    = np.zeros([xsteps,ysteps])
Qx     = np.zeros([xsteps,ysteps])
Qy     = np.zeros([xsteps,ysteps])
dQ     = np.zeros([xsteps,ysteps])
dQmask = np.zeros([xsteps,ysteps])

Qx0.fill(minimum)
Qx0.fill(minimum)
Qx.fill(minimum)
Qx.fill(minimum)
dQ.fill(minimum)
dQmask.fill(True)

x = np.array([x0+ix*dx for ix in range(xsteps)])
y = np.array([y0+iy*dy for iy in range(ysteps)])

logmin = np.abs(np.log10(minimum))
logmax = np.abs(np.log10(maximum))


tot = -1 # initialize total
file = open(filename, 'r')
for line in file:
	fileLine = line.split()
	if len(fileLine) == 0 or fileLine[0][0] == '!': continue
	
	tot = tot+1

	if any(fileLine[ix] == '*' for ix in range(len(fileLine))):
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
y = [iy*yscale for iy in y]

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

if masked: 
	xA = ma.masked_outside(xA, -hlim, hlim)
	yA = ma.masked_outside(yA, -vlim, vlim)


#####################
# Plotting:

# manually reset limits for plotting:
x0 = -30./xscale
x1 = 30./xscale
y0 = 0./yscale
y1 = 30./yscale


#plt.clf()
#fig = plt.figure(num=1,figsize=(20.,6),dpi=72)
fig = plt.figure(num=1,figsize=(6.67,4.),dpi=150)
fig.subplots_adjust(left=0.1,right=0.98,top=0.93,bottom=0.1)
ax1 = fig.add_subplot(111)
#ax1.set_aspect('equal')
mappable = ax1.pcolor(xA,yA,dQAmasked, vmin=-logmin, vmax=-logmax)
ax1.set_xlim(x0*xscale,x1*xscale)
ax1.set_ylim(y0*yscale,y1*yscale)
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.plot((hlim,hlim),(0,vlim), ls='--', c='r',lw=4)
ax1.plot((-hlim,-hlim),(0,vlim), ls='--', c='r',lw=4)
ax1.plot((-hlim,hlim),(vlim,vlim), ls='--', c='r',lw=4)
fig.colorbar(mappable, ax=ax1)
ax1.set_title('$\Delta Q$ vs. (x,y)')
fig.tight_layout()
plt.savefig('x_y_dQ.png')
plt.savefig('x_y_dQ.eps')


###########################
# plotting of \Delta Q vs. (Qx,Qy)

fig2 = plt.figure(num=2,figsize=(6.67,4.),dpi=150)
fig2.subplots_adjust(left=0.1,right=0.94,top=0.93,bottom=0.1)
ax2 = fig2.add_subplot(111)
ax2.scatter(Qx0A.reshape(-1),Qy0A.reshape(-1), s=2.5, c=dQAmasked.reshape(-1), norm=colors.Normalize(vmin=-logmin, vmax=-logmax, clip=True), marker='s', edgecolors='none')
#if PlotAnalytic == True:
#	analyticTunePlane(ax2)
ax2.set_xlabel('Qx')
ax2.set_ylabel('Qy')
ax2.axis('equal')
ax2.set_xlim(Qxmin, Qxmax)
ax2.set_ylim(Qymin, Qymax)
ax2.grid()
ax2.set_title('$\Delta Q$ vs. ($Q_x$,$Q_y$)')
fig2.tight_layout()
plt.savefig('Qx_Qy_dQ.png')
plt.savefig('Qx_Qy_dQ.eps')
plt.show()
