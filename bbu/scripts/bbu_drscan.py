#!/usr/bin/env python
import bbu_find_threshold
import bbu_params
import subprocess
import pylab
import os
import matplotlib.pyplot as plt
import numpy as np

#######################
def setup_dr_scan ( arc_time ):
#######################
    # Make lat2 txt file
    print("Filling LAT2 txt file now")
    # Calculate the arclength for this given recirculation time
    arc_l = calc_arcl ( arc_time )
    my_file = open('lat2.lat','w')
    my_file.write('arc[L] = '+str(arc_l))
    my_file.close()


#######################
def calc_arcl ( arctime ):
#######################
  c = 2.99792458*10**8
  beta = 0.999987
  # Find the new arc length for given arc time
  arc_l = arctime*beta*c    # Beta for speed of bunch (varies by lattice, can just use approximation too)
  #arc_l = arctime * c  # Approximation
  return arc_l


#######################
def make_dr_plot ( ):
####################### 
  x = []
  y = []
  f = []
  lines = []
  p = []
  my_plotfile = "scripts/thresh_v_trotb.txt"
  f = open(my_plotfile, 'r')
  lines = f.readlines()
  f.close()

  for line in lines:
    p = line.split()
    x.append(float(p[0]))
    y.append(float(p[1]))

    xv = np.array(x)
    yv = np.array(y)


  plt.scatter(xv, yv, marker = 'o', color = 'b')
  plt.title("BBI paper reproduction")
  plt.xlabel("tr/tb")
  plt.ylabel("Log( HOM Voltage )")
  #plt.text(.15, .9, 'PRSTAB 7 (2004) Fig. 3.')
  plt.show()

