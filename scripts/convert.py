#!/usr/bin/python

import sys
import shutil
import os

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open('temp.temp', 'w')
  found = False

  for line in in_file.readlines(): 

    if 'GENERAL1' in line:
      line = line.replace('GENERAL1', 'CUSTOM1')
      found = True

    if 'GENERAL2' in line:
      line = line.replace('GENERAL2', 'CUSTOM2')
      found = True

    if 'GENERAL3' in line:
      line = line.replace('GENERAL3', 'CUSTOM3')
      found = True

    if 'general1' in line:
      line = line.replace('general1', 'custom1')
      found = True

    if 'general2' in line:
      line = line.replace('general2', 'custom2')
      found = True

    if 'general3' in line:
      line = line.replace('general3', 'custom3')
      found = True

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
