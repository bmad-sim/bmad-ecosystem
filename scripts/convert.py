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

    if 'ix_track_end' in line:
      found = True
      line = line.replace('ix_track_end', 'track_state')

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
