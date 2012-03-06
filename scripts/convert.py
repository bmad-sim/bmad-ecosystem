#!/usr/bin/python

import sys
import shutil
import os

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open('temp.temp', 'w')
  found = False

  for line in in_file.readlines(): 

    if 'single-mode' in line:
      line = line.replace('single-mode', 'single_mode')
      found = True

    if 'end-file' in line:
      line = line.replace('end-file', 'end_file')
      found = True

    if 'xy-scale' in line:
      line = line.replace('xy-scale', 'xy_scale')
      found = True

    if 'x-scale' in line:
      line = line.replace('x-scale', 'x_scale')
      found = True

    if '{x-axis}' in line:
      line = line.replace('{x-axis}', '{x_axis}')
      found = True

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
