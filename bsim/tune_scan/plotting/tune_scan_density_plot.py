#!/usr/bin/env python

# To execute this file from within the Python interpreter:
#   execfile('this_file_name')

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

T = True    # For converting from data file header T/F notation
F = False

# Check version

if sys.version_info[0] == 2 or sys.version_info[1] < 6: sys.exit("MUST USE PYTHON 3.6 OR GREATER!")

# Help

def print_help():
  print ('''
Usage:
  tune_scan_density_plot.py {-cmap <color_map>} {-column <col_to_plot>} {-min <min_val>} {-max <max_val>} {-z <z_index>} {<data_file_name>}

Defaults:
  <color_map>    = "gnuplot2_r" # Color map. Google "matplotlib color maps" for more info.
  <col_to_plot>  = 9            # Column index (0 = first column, 9 = Amp_b_rms).
  <max_val>      = -1           # If positive: Set all data values = min(max_val, data_val).
  <min_val>      = 0            # Set all data values = max(min_val, data_val).
  <z_index>      = 0            # Qz slice index (0 = first slice at Qx = Q_z0).
  <data_file_name> = tune_scan.dat
''')
  exit() 

# Defaults

dat_file_name = 'tune_scan.dat'
dat_col = 8
z_index = 0
min_val = 0
max_val = -1
cmap = 'gnuplot2_r'

# Command line arguments

i = 1
while i < len(sys.argv):
  n = len(sys.argv[i])
  if sys.argv[i] == '-':
    print_help()

  elif sys.argv[i] == '-column'[:n]:
    dat_col = int(sys.argv[i+1])
    i += 1

  elif sys.argv[i] == '-max'[:n]:
    max_val = float(sys.argv[i+1])
    i += 1

  elif sys.argv[i] == '-min'[:n]:
    min_val = float(sys.argv[i+1])
    i += 1

  elif sys.argv[i] == '-z'[:n]:
    z_index = int(sys.argv[i+1])
    i += 1

  elif sys.argv[i][0] == '-':
    print_help()
    
  else:
    dat_file_name = sys.argv[i]

  i += 1

# Read data file parameters in header lines and data file data

dat_file = open (dat_file_name)

for n_header in range(1, 1000):
  line = dat_file.readline()
  if line[0:2] == '#-': break
  print (line.strip())
  exec (line[1:].strip())

col_label = line.split()[1:]
dat_file.close()

# [jx, jy, jz, n_turn, dat]
pix_dat = np.loadtxt(dat_file_name, usecols=(0,1,2,6,dat_col), skiprows = n_header)

# Create density matrix

na_min = 0
nb_min = 0

print (f'na min/max: {na_min}, {na_max}')
print (f'nb min/max: {nb_min}, {nb_max}')

pix_mat = np.zeros((na_max+1, nb_max+1))
if max_val < 0: max_val = pix_dat[:,4].max()

for pix in pix_dat:
  if round(pix[2]) != z_index: continue
  if pix[3] == 0:
    pix_mat[int(round(pix[0])), int(round(pix[1]))] = max_val
  else:
    pix_mat[int(round(pix[0])), int(round(pix[1]))] = max(min(pix[4], max_val), min_val)

# And plot

x_min = Q_a0
x_max = Q_a1
y_min = Q_b0
y_max = Q_b1

if 'ix_plot' not in locals(): ix_plot = 0
ix_plot = ix_plot + 1

fig = plt.figure(ix_plot)
ax = fig.add_subplot(111)

if x_max-x_min > y_max - y_min:
  it = max(1, int((y_max-y_min) * 8 / (x_max-x_min)))
  ax.yaxis.set_major_locator(ticker.MaxNLocator(it))
else:
  it = max(1, int((x_max-x_min) * 8 / (y_max-y_min)))
  ax.xaxis.set_major_locator(ticker.MaxNLocator(it))

dens = ax.imshow(np.transpose(pix_mat), origin = 'lower', extent = (x_min, x_max, y_min, y_max))
dens.set_cmap(cmap)
ax.set_xlabel('Qa')
ax.set_ylabel('Qb')
if rf_on:
  ax.set_title(f'{col_label[dat_col]}  Qz: {Q_z0+z_index*dQ_z}')
else:
  ax.set_title(f'{col_label[dat_col]}  pz: {pz0+z_index*dpz}')


fig.colorbar(dens)

plt.show()
