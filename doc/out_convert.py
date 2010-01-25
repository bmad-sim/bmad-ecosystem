#!/usr/bin/python

import sys, os

file_in = open('subroutines.tex')
file_out = open('out.tex', mode = 'w')

for arg in sys.argv[1:]:
  print arg
  file_in = open(arg)
  file_out = open(arg + 'x', mode = 'w')

  for line in file_in.readlines():
    ixh = line.find(r']{routine!')
    if ixh != -1:
      line = line[0:ixh+2] + line[ixh+10:]
      print line
    file_out.write(line)
    
  file_in.close()
  file_out.close()
  os.rename(arg + 'x', arg)
