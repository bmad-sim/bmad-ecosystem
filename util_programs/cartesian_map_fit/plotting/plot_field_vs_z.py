import numpy as np
import matplotlib.pyplot as plt
import time
# Input parameters

dat_file_name = 'bmad_field.table'
fit_file_name = 'fit.table'
x_val = 0.1    # Must correspond to a grid x-position value
y_val = 0.2    # Must correspond to a grid y-position value
plot_type = 'Bz'

# Read header

dat_file = open(dat_file_name, 'r')

vals = dat_file.readline().replace(',', ' ').split('!')[0].split()
pos_scale = float(vals[0])

vals = dat_file.readline().replace(',', ' ').split('!')[0].split()
field_scale = float(vals[0])

vals = dat_file.readline().replace(',', ' ').split('!')[0].split()
nx_min = float(vals[0])
nx_max = float(vals[1])

vals = dat_file.readline().replace(',', ' ').split('!')[0].split()
ny_min = float(vals[0])
ny_max = float(vals[1])

vals = dat_file.readline().replace(',', ' ').split('!')[0].split()
nz_min = float(vals[0])
nz_max = float(vals[1])

vals = dat_file.readline().replace(',', ' ').split('!')[0].split()
del_x = float(vals[0])
del_y = float(vals[1])
del_z = float(vals[2])

vals = dat_file.readline().replace(',', ' ').split('!')[0].split()
x_off = float(vals[0])
y_off = float(vals[1])
z_off = float(vals[2])

dat_file.close()

# Read data and fit tables

if plot_type == 'Bx':
  b_row = 3
elif plot_type == 'By':
  b_row = 4
else:
  b_row = 5

x_dat = []
y_dat = []

dat_table = np.loadtxt(dat_file_name, skiprows = 7)
for row in dat_table:
  if abs(row[0] - x_val) > 1e-5*del_x: continue
  if abs(row[1] - y_val) > 1e-5*del_y: continue
  x_dat.append(row[2])
  y_dat.append(row[b_row])


x_fit = []
y_fit = []

fit_table = np.loadtxt(fit_file_name, skiprows = 7)
for row in fit_table:
  if abs(row[0] - x_val) > 1e-5*del_x: continue
  if abs(row[1] - y_val) > 1e-5*del_y: continue
  x_fit.append(row[2])
  y_fit.append(row[b_row])

# Plot

plt.plot(x_dat, y_dat, 'ro', markersize = 1, label = dat_file_name)
plt.plot(x_fit, y_fit, 'go', markersize = 1, label = fit_file_name)
plt.plot([], [], ' ', label = 'x = ' + str(x_val))  # For creating the legend
plt.plot([], [], ' ', label = 'y = ' + str(y_val))

plt.xlabel('Z')
plt.ylabel(plot_type)
plt.legend()

plt.show()
