#!/usr/bin/env python

import subprocess
import pylab
import os
import matplotlib.pyplot as plt
import numpy as np
import math

#######################
def setup_phase_scan ( phase, py_par ):
#######################
  # Make lat2 txt file
  # print("Filling lat2.txt with the new phase (Taylor) ")
  # Calculate the arclength for this given recirculation time

  bW = 13.58259507     # 1-pass 2016_01_25, after LA.END.MAR\1  
  aW = -0.78076767    
  #bW = 34.54001603    # 4-pass 2016_05_17, after LA.END.MAR\1
  #aW = -1.68026439
  #bW = 30.36578721     # 4-pass 2016_09_01, after LA.END.MAR\1
  #aW = -1.21024368


  gW = (1+aW*aW)/bW

  pW = float(phase)

  m11W = math.cos(pW) + aW*math.sin(pW)
  m12W = bW*math.sin(pW)
  m21W = -gW*math.sin(pW)
  m22W = math.cos(pW) - aW*math.sin(pW)

  #taylorW: Taylor, {1: m11W | 1}, {1:  m12W | 2}, {2: m21W | 1}, {2:  m22W | 2}

  my_file = open(os.path.join(py_par['temp_dir'],'lat2.lat'),'w')
    
  print("Test phase included in lat2: ", pW)
  my_file.write('taylorW[tt11]='+ str(m11W)+'\n')
  my_file.write('taylorW[tt12]='+ str(m12W)+'\n')
  my_file.write('taylorW[tt21]='+ str(m21W)+'\n')
  my_file.write('taylorW[tt22]='+ str(m22W)+'\n')

  my_file.close()

#######################
def setup_phase_xy_scan ( py_par ):
#######################
  # Make lat2 txt file
  # print("Filling lat2.txt with the new phase (Taylor) ")
  # Calculate the arclength for this given recirculation time
    
  phasex = py_par['phase_x']
  phasey = py_par['phase_y']

  #1-pass 2016_01_25, after LA.END.MAR\1
  #bxW = 13.58259507  
  #axW = -0.78076767    
  #byW = 13.48215194
  #ayW = -0.75675121
   
  # 4-pass 2016_05_17, after LA.END.MAR\1
  #bxW = 34.54001603   
  #axW = -1.68026439
  #byW = 34.11944669
  #ayW = -1.65547205
    
  # 4-pass 2016_09_01, after LA.END.MAR\1
  #bxW = 30.36578721   
  #axW = -1.21024368
  #byW = 29.99476043
  #ayW = -1.19124485
  
  # 4-pass 2016_11_17, after LA.MAR.END\1
  # Also for 1-pass 2016_12_12
  bxW = 31.40988174
  axW = -1.35641640
  byW = 31.02811471
  ayW = -1.33553551
  
  pxW = float(phasex)
  pyW = float(phasey)
  
  if (py_par['xy_coupled'] == 0):

    gxW = (1+axW*axW)/bxW
    gyW = (1+ayW*ayW)/byW
    
    m11W = math.cos(pxW) + axW*math.sin(pxW)
    m12W = bxW*math.sin(pxW)
    m21W = -gxW*math.sin(pxW)
    m22W = math.cos(pxW) - axW*math.sin(pxW)
    m33W = math.cos(pyW) + ayW*math.sin(pyW)
    m34W = byW*math.sin(pyW)
    m43W = -gyW*math.sin(pyW)
    m44W = math.cos(pyW) - ayW*math.sin(pyW)

    my_file = open(os.path.join(py_par['temp_dir'],'lat2.lat'),'w')
    
    print("DECOUPLED X-Y phase included in lat2: ", pxW, pyW)
    my_file.write('taylorW[tt11]='+ str(m11W)+'\n')
    my_file.write('taylorW[tt12]='+ str(m12W)+'\n')
    my_file.write('taylorW[tt21]='+ str(m21W)+'\n')
    my_file.write('taylorW[tt22]='+ str(m22W)+'\n')
    my_file.write('taylorW[tt33]='+ str(m33W)+'\n')
    my_file.write('taylorW[tt34]='+ str(m34W)+'\n')
    my_file.write('taylorW[tt43]='+ str(m43W)+'\n')
    my_file.write('taylorW[tt44]='+ str(m44W)+'\n')
    
    my_file.close()

  elif (py_par['xy_coupled'] == 1):
   
    m13W = math.sqrt(bxW/byW)*(math.cos(pxW) + ayW*math.sin(pxW))
    m14W = math.sqrt(bxW*byW)*math.sin(pxW)
    m23W = ((ayW-axW)*math.cos(pxW)-(1+ayW*axW)*math.sin(pxW))/math.sqrt(bxW*byW)
    m24W = math.sqrt(byW/bxW)*(math.cos(pxW) - axW*math.sin(pxW))

    m31W = math.sqrt(byW/bxW)*(math.cos(pyW) + axW*math.sin(pyW))
    m32W = math.sqrt(byW*bxW)*math.sin(pyW)
    m41W = ((axW-ayW)*math.cos(pyW)-(1+axW*ayW)*math.sin(pyW))/math.sqrt(byW*bxW)
    m42W = math.sqrt(bxW/byW)*(math.cos(pyW) - ayW*math.sin(pyW))

    # Lattice file syntax
    # taylorW: Taylor, {1: m11W | 1}, {1:  m12W | 2}, {2: m21W | 1}, {2:  m22W | 2}

    my_file = open(os.path.join(py_par['temp_dir'],'lat2.lat'),'w')
    
    print(" COUPLED X-Y phase included in lat2: ", pxW, ", ", pyW)
    
    my_file.write('taylorW[tt13]='+ str(m13W)+'\n')
    my_file.write('taylorW[tt14]='+ str(m14W)+'\n')
    my_file.write('taylorW[tt23]='+ str(m23W)+'\n')
    my_file.write('taylorW[tt24]='+ str(m24W)+'\n')
    
    my_file.write('taylorW[tt31]='+ str(m31W)+'\n')
    my_file.write('taylorW[tt32]='+ str(m32W)+'\n')
    my_file.write('taylorW[tt41]='+ str(m41W)+'\n')
    my_file.write('taylorW[tt42]='+ str(m42W)+'\n')

    my_file.write('taylorW[tt11]= 0.0 '+'\n')
    my_file.write('taylorW[tt22]= 0.0 '+'\n')
    my_file.write('taylorW[tt33]= 0.0 '+'\n')
    my_file.write('taylorW[tt44]= 0.0 '+'\n')
    
    my_file.close()

  else:
    print("py_par['xy_coupled'] must be 0 or 1 !!")
#######################
def make_phase_plot ( py_par ):
####################### 
  x = []
  y = []
  f = []
  lines = []
  p = []
  my_plotfile = "thresh_v_phase.txt"
  f = open(os.path.join(py_par['temp_dir'],my_plotfile), 'r')
  lines = f.readlines()
  f.close()

  for line in lines:
    p = line.split()
    x.append(float(p[0]))
    y.append(float(p[1]))

    xv = np.array(x)
    yv = np.array(y)

  plt.scatter(xv, yv, marker = 'o', color = 'b')
  #plt.title("Phase Scan for Q=10^-4, R/Q=100Ohm, f=2E9Hz")
  plt.rcParams.update({'font.size': 20})
  plt.xlabel("Phase")
  plt.ylabel("Ith(A)")
  #plt.text(.15, .9, 'PRSTAB 7 (2004) Fig. 3.')
  #fig = plt.figure()
  #ax = fig.add_subplot(111)  
  #ax.annotate('PRSTAB 7 (2004) Fig. 3.',xy=(0,0))

  plt.yscale('log')
  plt.show()

