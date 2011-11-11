#!/usr/bin/python

#+
# Script to run regression tests for the Bmad libraries.
#-

import re
import os
import subprocess

# First make a list of tests to run by looking in Makefile for
# the M_FILE_LIST list.


re_this = re.compile("^ *M_FILE_LIST *:= *(.*)")

makefile = open ('Makefile', 'r')
for line in makefile.readlines():
  match_obj = re_this.match(line)
  if match_obj:
    programs = match_obj.group(1).replace('M.', '').split()
    break

makefile.close()

# Print the list of programs

print 'Programs:'
for prog in programs:
  print '    ' + prog 

# Run the programs

for prog in programs:
  os.chdir(prog)
  print '/n'
  print 'Running: ' + prog
  process = subprocess.Popen(['../../bin/' + prog], shell = True)

  # Compare the output of the program "output.now" to the expected output "output.correct"

  now_file = open('output.now', 'r')  
  correct_file = open('output.correct', 'r')

  for now_line in now_file.readline():
    if now_line.strip[0] == '!': continue
    if not now_line.strip: continue blank line


  if len(now_line) == 0: 


  os.chdir('..')
  

