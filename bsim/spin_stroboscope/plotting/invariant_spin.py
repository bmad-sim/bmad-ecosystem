#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

who = 'z'
if len(sys.argv) == 3: who = sys.argv[2]


dat_file = open(sys.argv[1], 'r')
line = dat_file.readline()  # Skip header

xm = []
ym = []

for line in dat_file.readlines():
  lst = line.split()

  if who == 'x':
    x = float(lst[ 9])
    y = float(lst[10])
    z = float(lst[11])
  elif who == 'y':
    x = float(lst[12])
    y = float(lst[13])
    z = float(lst[14])
  elif who == 'z':
    x = float(lst[15])
    y = float(lst[16])
    z = float(lst[17])
  else:
    print ('BAD WHO')
    sys.exit()

  xm.append(np.arctan2(x, z))
  ym.append(np.arcsin(y))

#

fig = plt.figure()

plt.subplot(111, projection="mollweide", facecolor ='#F6FFFF') # Very light cyan
plt.grid(True, linestyle = 'dashed')

plt.plot(xm, ym, 'bo', markersize = 3)
plt.title("y")
plt.text(3.2, 0.0, "x", size = 'large')


plt.show()
