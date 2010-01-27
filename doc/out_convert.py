#!/usr/bin/python

import sys, os

for arg in sys.argv[1:]:
  print arg
  file_in = open(arg)
  file_out = open(arg + 'x', mode = 'w')

  for line in file_in.readlines():

    while True:
      ix1 = line.find('hyperref[')
      if ix1 == -1: break
      ix2 = line.find(']', ix1)
      line = line[0:ix1] + 'Hyperref{' + line[ix1+9:ix2] + '}' + line[ix2+1:]
      ## print line

    file_out.write(line)
    
  file_in.close()
  file_out.close()
  os.rename(arg + 'x', arg)
