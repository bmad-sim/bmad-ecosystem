#!/usr/bin/python

import sys
import shutil
import os

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open('temp.temp', 'w')
  found = False

  for line in in_file.readlines(): 

    if 'genplt mpm_utils cesr_utils sim_utils mpmnet' in line:
      line = line.replace('genplt mpm_utils cesr_utils sim_utils mpmnet', 'genplt mpm_utils cesr_utils sim_utils mpmnet cbi_net c_utils')
      found = True

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
