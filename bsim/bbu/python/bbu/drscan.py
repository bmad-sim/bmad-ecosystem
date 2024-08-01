#!/usr/bin/env python3

import subprocess
import pylab
import os
import matplotlib.pyplot as plt
import numpy as np

#######################
def setup_drscan ( arc_time, py_par ):
#######################
    # Make lat2 txt file
    print("Filling lat2.txt with the new arclength")
    # Calculate the arclength for this given recirculation time
    arc_l = calc_arcl ( arc_time )
    my_file = open(os.path.join(py_par['temp_dir'],'lat2.lat'),'w')
    my_file.write('arc[L] = '+str(arc_l))
    my_file.close()

#######################
def calc_arcl ( arctime ):
#######################
  c = 299792458
  # Find the new arc length for given arc time
  arc_l = arctime * c  
  return arc_l

#######################
def make_dr_plot ( py_par ):
####################### 
  x = []
  y = []
  f = []
  lines = []
  p = []
  my_plotfile = "thresh_v_trotb.txt"
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
#  plt.title("DR Scan for Q=10^-4, R/Q=100Ohm, f=2E9Hz, m12=10")
  plt.rcParams.update({'font.size': 20})
  plt.xlabel("Arc Time / Bunch Time")
  plt.ylabel("Ith (A)")
  #plt.text(.15, .9, 'PRSTAB 7 (2004) Fig. 3.')
  #fig = plt.figure()
  #ax = fig.add_subplot(2,1,1)  
  #ax.annotate('PRSTAB 7 (2004) Fig. 3.',xy=(0,0))
  plt.yscale('log')
  plt.show()

