#!/usr/bin/env python

# Script to read in Fortran structures and create:
#   Fortran bind(c) version.
#   C version
#   Translators between Fortran original structure and bind(c) version.
#   Translators between C version and C++ version

import sys
import shutil
import os
import copy

#---------------------------------------------------------
# Class defs

class f_var_class:

  def __init__(self):
    self.name = ''
    self.type = ''
    self.pointer_type = '-'    # 'pointer', 'allocatable'
    self.array = ''            # [':', ':'] or ['6']
    self.init_value = ''
    self.comment = ''
    
  def __repr__(self):
    return '[%s, %s, %s, %s, %s]' % (self.type, self.pointer_type, self.name, self.array, self.init_value)

class f_struct_def_class:
  def __init__(self):
    self.name = ''
    self.var = []

  def __repr__(self):
    return '[name: %s, #var: %i]' % (self.name, len(self.var))

#---------------------------------------------------------
# Get the list of structs

f_struct_list_file = open('fortran_structs.list')

f_struct_def = {}
f_struct_list = []
f_module_files = []

for line in f_struct_list_file:
  line = line.strip()
  if len(line) == 0: continue
  split_line = line.split()

  if split_line[0] == 'FILE:':
    f_module_files.append(split_line[1])
  else:
    f_struct_list.append(split_line)
    f_struct_def[line.split()[0]] = f_struct_def_class()

f_struct_list_file.close()

### print f_struct_list

#--------------------------------------------------------
# Parse structure definitions

for file_name in f_module_files:
  f_module_file = open('../../' + file_name)

  for line in f_module_file:
    split_line = line.split()
    if len(split_line) < 2: continue
    if split_line[0] != 'type': continue
    if not split_line[1] in f_struct_def: continue

    struct_name = split_line[1]
    f_struct_def[struct_name].name = struct_name

    # Now collect the struct variables

    for line in f_module_file:
      var = f_var_class()

      part = line.partition('!')
      var.comment = part[2].strip()

      part = part[0].partition('=')
      var.init_value = part[2].strip()
      if var.init_value[:1] == '>': var.init_value = ''  # Ignore pointer "=> null()" init.

      split_line = part[0].replace(',', ' ').replace('::', ' ').split()
      if len(split_line) == 0: continue
      if split_line[0:2] == ['end', 'type']: break

      var.type = split_line.pop(0)
      if var.type == 'type': var.type = 'type:' + split_line.pop(0)[1:-1]

      if split_line[0] == 'allocatable' or split_line[0] == 'pointer':
        var.pointer_type = split_line.pop(0)

      # Handle lines like: "real(rp) kx, ky, kz"

      for name in split_line:
        var = copy.copy(var)
        var.name = name
        part = var.name.partition('(')
        var.name = part[0]
        var.array = part[2][:-1]     # Remove last ')'
        f_struct_def[struct_name].var.append(var)
        
  f_module_file.close()

#--------------------------------------------------------
# Write results to file

f_out = open('f_structs.parsed', 'w')
for struct in f_struct_def:
  f_out.write('******************************************\n')
  f_out.write (struct + '    ' + str(len(f_struct_def[struct].var)) + '\n')
  for var in f_struct_def[struct].var:
    f_out.write ('    ' + str(var) + '\n')

#--------------------------------------------------------
# Create Fortran bind(c) struct version


#--------------------------------------------------------
# Create C struct version

#--------------------------------------------------------
# Create C++ class version

#--------------------------------------------------------
# Create Fortran to C++ code

#--------------------------------------------------------
# Create C++ to Fortran code

#--------------------------------------------------------
# Create Fortran struct equality check code

#--------------------------------------------------------
# Create C++ class equality check code

#--------------------------------------------------------
# Create interface check code

