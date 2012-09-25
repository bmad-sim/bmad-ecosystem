#!/usr/bin/python

import sys
import shutil
import os

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open('temp.temp', 'w')
  found = False

  for line in in_file.readlines(): 

    if 'bmad_status%type_out' in line:
      line = line.replace('bmad_status%type_out', 'global_com%type_out')
      found = True

    if 'bmad_status%exit_on_error' in line:
      line = line.replace('bmad_status%exit_on_error', 'global_com%exit_on_error')
      found = True

    if 'call err_exit' in line:
      line = line.replace('call err_exit', 'if (global_com%exit_on_error) call err_exit')
      found = True

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
