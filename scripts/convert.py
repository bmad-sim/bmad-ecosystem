#!/usr/bin/python

import sys
import shutil
import os
import re

g = re.compile(r'(\s*)(call)( element_at_s.*?,.*?,) *(\w*)(.*)', re.IGNORECASE)

print 'Number of files: ' + str(len(sys.argv[1:]))

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open(r'temp.temp', 'w')
  found = False

  for line in in_file.readlines(): 

    if 'n_attrib_maxx' in line:
      found = True
      line = line.replace('n_attrib_maxx', 'num_ele_attrib$')

    if 'n_attrib_special_maxx' in line:
      found = True
      line = line.replace('n_attrib_special_maxx', 'num_ele_attrib_extended$')

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
