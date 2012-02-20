#!/usr/bin/python

import sys
import shutil
import os

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open('temp.temp', 'w')
  found = False

  for line in in_file.readlines(): 

    if 'lat_compute_reference_energy' in line and 'if' not in line: 
      line = line.replace('lat_compute_reference_energy', 'lat_compute_ref_energy_and_time')
      found = True

    if 'compute_ele_reference_energy' in line and 'if' not in line: 
      line = line.replace('compute_ele_reference_energy', 'ele_compute_ref_energy_and_time')
      found = True

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
