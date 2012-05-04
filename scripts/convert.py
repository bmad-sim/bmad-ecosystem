#!/usr/bin/python

import sys
import shutil
import os
import re

g = re.compile(r'(\s*)(call)( element_at_s.*?,.*?, )(\w*)(.*)')

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open(r'temp.temp', 'w')
  found = False

  for line in in_file.readlines(): 

    gg = g.match(line)
    if gg:
      line = gg.group(1) + gg.group(4) + ' =' + gg.group(3) + '.true.' + gg.group(5) + '\n'
      found = True

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
