#!/usr/bin/python

import re
import glob
import os

# Script to make a VMS mms file with dependencies.

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Function to convert form Unix directory syntax to VMS directory syntax

def vms_file_syntax (unix_file):
  ix_s = unix_file.rfind('/')
  if ix_s == -1: return unix_file
  vms_file = '[.' + unix_file[0:ix_s] + ']' + unix_file[ix_s+1:]
  return vms_file.replace('/', '.')

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Function to break a long string into lines of a certain length.
# "indent" is the number of spaces to indent after the first line.

def break_into_lines (string, line_len = 80, end_char = '-', indent = 5):

  out = ''
  string = string.strip()

  while len(string) > line_len:
    ix = string.rfind(' ', indent, line_len)
    if ix == -1: ix = string.find(' ', line_len)
    if ix == -1: break
    out += string[:ix].rstrip() + ' ' + end_char + '\n'
    string = ' ' * indent + string[ix+1:].strip()

  out += string
  
  return out

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Get the current directory

root_dir = os.getcwd().split('/')[-1]

# Read the Makefile and look for the LIB_SRC_DIRS subdirectories

p = re.compile("^ *LIB_SRC_DIRS *:= *(.*)")

for line in open ("Makefile", 'r'):
  m = p.match(line)
  if (m):
    libs =  m.group(1).split()
    print libs

# Make a list of all routine files in the subdirectories

full_file_list = []   # List of files with directory. EG: './modules.bmad_struct.f90'
file_list = []        # List of files w/o directory. EG: 'bmad_struct.f90'
prefix_list = []      # List of files w/o suffix. EG: 'bmad_struct'
mod_list = []         # List with 'xxx.mod' for .f90 files and 'xxx.cpp' for .cpp files

for subdir in libs:
  full_file_list.extend(glob.glob(subdir + '/*.f90')) 
  full_file_list.extend(glob.glob(subdir + '/*.cpp'))

for file in full_file_list:
  file = file.split('/')[-1] # EG: .split -> ['.', 'modules', 'bmad.f90']
  file_list.append(file)
  prefix_list.append(file.split('.')[0])
  if '.f90' in file:
    mod_list.append(file.replace('.f90', '.mod'))
  else:
    mod_list.append(file)

# Make a list of the source files

print break_into_lines('file_source = ' + ', '.join(file_list))
print '\n!-----------------------------------------------------\n'


# For each routine file: Find the corresponding dependency file

for i in range(len(prefix_list)):
  prefix = prefix_list[i]
  depend_file = '.depend/' + prefix + '.d'
  
  # Now open the dependency file and read in the dependencies

  depend_upon_list = []

  for line in open (depend_file, 'r'):

    if line.isspace(): continue
    if '# DO NOT DELETE' in line: continue
    if 'CESR_platform.inc' in line: continue
    if '/lib/' in line: continue
    if '/usr/' in line: continue

    if line.find('$(LIBRARY)(') != 0:
      print 'Problem dependency: ' + line
      stop

    # Everything after the ':' is the list of dependencies we are after 
    d_line = line.split(':', 1)[1]
    depend_upon_list.extend(d_line.split())

  # Throw out dependencies on files that are not part of this library.

  mini_depend_list = []
  for j in range(len(depend_upon_list)):
    d1 = depend_upon_list[j]
    d1 = d1.replace('.mod:', '.mod')
    d1 = d1.split('/')[-1]
    if d1 in mod_list: mini_depend_list.append(d1)

  # Print the dependencies

  d_line = root_dir + '(' + file_list[i] + ') : ' 
  d_line += vms_file_syntax(full_file_list[i].replace('.', root_dir, 1))
  for d in mini_depend_list:
    if '.mod' in d:
      d_line += ", " + vms_file_syntax("modules/" + d.replace('.mod', '.f90$mod'))
    else:
      d_line += ", " + vms_file_syntax(d)

  print break_into_lines(d_line)
