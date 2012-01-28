#!/usr/bin/python

import sys
import shutil
import os

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open('temp.temp', 'w')
  found = False

  for line in in_file.readlines(): 

    if 'auto_rf_phase_and_amp' in line: 
      line = line.replace('auto_rf_phase_and_amp', 'rf_auto_phase_and_amp_correct')
      found = True

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
