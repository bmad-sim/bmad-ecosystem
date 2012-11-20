#!/usr/bin/python

import sys
import shutil
import os

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open('temp.temp', 'w')
  found = False
  skip_blank = False

  for line in in_file.readlines(): 

#    if '#include "CESR_platform.inc' in line:
#      line = line.replace('scaler', 'scalar')
#      found = True

    if '#include "CESR_platform.inc' in line:
      found = True
      skip_blank = True
      continue

    if '#include "CESR_platform.h' in line:
      found = True
      skip_blank = True
      continue

    if skip_blank and line == '\n': 
      skip_blank = False
      continue

    skip_blank = False
    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
