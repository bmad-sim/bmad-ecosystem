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

    #bW = 13.58259507   # 1-pass, after LA.END.MAR  
    #aW = -0.78076767    
    bW = 34.54001603    # 4-pass, after LA.END.MAR
    aW = -1.68026439
    
    gW = (1+aW*aW)/bW

    pW = float(phase)

    m11W = math.cos(pW) + aW*math.sin(pW)
    m12W = bW*math.sin(pW)
    m21W = -gW*math.sin(pW)
    m22W = math.cos(pW) - aW*math.sin(pW)

    #taylorW: Taylor, {1: m11W | 1}, {1:  m12W | 2}, {2: m21W | 1}, {2:  m22W | 2}

    my_file = open(os.path.join(py_par['temp_dir'],'lat2.lat'),'w')
    
    #my_file.write('taylorW: Taylor, {1:'+str(m11W)+'|1}, {1:'+ str(m12W)+'|2}, {2:'+str(m21W)+'|1}, {2:'+ str(m22W)+'|2}')
    #my_file.write('arc: Taylor, {1:'+str(m11W)+'|1}, {1:'+ str(m12W)+'|2}, {2:'+str(m21W)+'|1}, {2:'+ str(m22W)+'|2}')
    
    #my_file.write('TaylorW[L]= 10')
    #my_file.write('tt: taylor, {1: 10 | 1}')
    #my_file.write('tt: taylor, tt13 = 0.5')
    print("Test phase included in lat2: ", pW)
    #print(m11W,m12W,m21W,m22W)
    my_file.write('taylorW[tt11]='+ str(m11W)+'\n')
    my_file.write('taylorW[tt12]='+ str(m12W)+'\n')
    my_file.write('taylorW[tt21]='+ str(m21W)+'\n')
    my_file.write('taylorW[tt22]='+ str(m22W)+'\n')
    #my_file.write('arc[tt12]='+ str(phase)+'\n')

    my_file.close()

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
    y.append(float(p[2]))

    xv = np.array(x)
    yv = np.array(y)

  plt.scatter(xv, yv, marker = 'o', color = 'b')
  plt.title("Phase Scan for Q=10^-4, R/Q=100Ohm, f=2E9Hz")
  plt.rcParams.update({'font.size': 20})
  plt.xlabel("Phase")
  plt.ylabel("Ith")
  #plt.text(.15, .9, 'PRSTAB 7 (2004) Fig. 3.')
  fig = plt.figure()
  ax = fig.add_subplot(111)  
  ax.annotate('PRSTAB 7 (2004) Fig. 3.',xy=(0,0))

  plt.show()

