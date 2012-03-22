#!/usr/bin/env python

# Script to read in Fortran structures and create:
#   Corresponding C++ class
#   Translator between Fortran structure and C++ class
#   Routines to check for equality between instances of a given fortran structure.
#   Routines to check for equality between instances of a given C++ class
#   Program to check the Fortran / C++ translator

import sys
import shutil
import os
import copy

##################################################################################
##################################################################################
# f_var_class

class f_var_class:

  def __init__(self):
    self.name = ''
    self.type = ''
    self.pointer_type = '-'    # 'pointer', 'allocatable'
    self.array = ''            # [':', ':'] or ['6']
    self.init_value = '-'
    self.comment = ''
    self.n_array = 1

  def __repr__(self):
    return '[%s, %s, %s, (%s), %s]' % (self.type, self.pointer_type, self.name, self.array, self.init_value)

# f_struct_def_class

class f_struct_def_class:
  def __init__(self):
    self.name = ''
    self.var = []

  def __repr__(self):
    return '[name: %s, #var: %i]' % (self.name, len(self.var))

# f_arg_class

class f_arg_class:
  def __init__(self, name):
    self.c_name = name
    self.num = 0
    self.input = []

  def __repr__(self):
    return '[%s, %s]' % (self.c_arg, self.num)


##################################################################################
##################################################################################
# Translations

f_type_to_bind_c = {
            'real(rp)'    : 'real(c_double)',
            'complex(rp)' : 'complex(c_double_complex)',
            'integer'     : 'integer(c_int)',
            'logical'     : 'logical(c_bool)'}

##################################################################################
##################################################################################
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

##################################################################################
##################################################################################
# Parse structure definitions

# Current restrictions: Syntax to avoid:
#   Multiple inits: "real a = 5, b = 7"
#   No space in type parens: "type( abc_struct )" or "real( rp )"
#   Space in array def: "arr( : )"

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
      if split_line[1][0:1] == '(': var.type = var.type + split_line.pop(0)

      if var.type == 'type(':
        var.type = 'type:' + var.type[5:-1]

      if var.type == 'type': var.type = 'type:' + split_line.pop(0)[1:-1]

      if split_line[0] == 'allocatable' or split_line[0] == 'pointer':
        var.pointer_type = split_line.pop(0)

      # Handle lines like: "real(rp) kx, ky, kz"

      for ix, name in enumerate(split_line):
        var = copy.copy(var)
        var.name = name

        part = var.name.partition('(')
        var.name = part[0]
        var.array = part[2][:-1]     # Remove last ')'

        if ix+1 < len(split_line) and split_line[ix+1][0:1] == '(': 
          var.array = split_line.pop(ix+1)[:-1]

        f_struct_def[struct_name].var.append(var)
        
  f_module_file.close()

##################################################################################
##################################################################################
# As a check, write results to file. 

f_out = open('f_structs.parsed', 'w')

for struct in f_struct_def:
  f_out.write('******************************************\n')
  f_out.write (struct + '    ' + str(len(f_struct_def[struct].var)) + '\n')
  for var in f_struct_def[struct].var:
    f_out.write ('    ' + str(var) + '\n')

f_out.close()

##################################################################################
##################################################################################
# Create Fortran side of interface...

# First the header

f_out = open('../cpp_interface/bmad_and_cpp_mod.f90', 'w')

f_out.write ('''
!+
! Fortran side of the Bmad / C++ interface.
!
! File Generated by: create_interface.py
! Do not edit this file directly! 
!-

module bmad_and_cpp_mod

use bmad_struct
use bmad_interface
use fortran_and_cpp_mod

contains

''')


##############
# Loop over all the structures...

for struct in f_struct_def:

  ##############
  # Sort the variables into the argument list

  s_name = struct[0:-7]  # Strip off ending "_struct"

  # Form argument list

  args = [f_arg_class('c_real'), f_arg_class('c_int')]  # real and int arrays

  for var in f_struct_def[struct]:

    if var.pointer_type == 'allocatable' or var.pointer_type == 'pointer':

    elif var.type == 'character':

    elif var.array /= '':
      args.append(f_arg_class('var.name'))
      args[-1].num = 1
      for d1 in var.array.split(','):
        d1_apart = d1.partition(':')
        if d1_apart[2] == '':
          args[-1].num = args[-1].num * int(d1_apart[0])
        else:
          args[-1].num = args[-1].num * (int(d1_apart[2]) - int(d1_apart[0]))

    elif var.name[0:5] == 'type:':

    else:
      if var.type == 'real(rp)':
        args[0].num += 1
        args[0].input.append(var.name)
      elif var.type == 'complex(rp)':
        args[0].num += 1
        args[0].input.append(var.name)
      elif var.type == 'integer':
        args[1].num += 1
        args[1].input.append(var.name)
      elif var.type == 'logical':
        args[1].num += 1
        args[1].input.append(var.name)

  if args[1].num == 0: args.pop[1]
  if args[0].num == 0: args.pop[0]

  ##############
  # Write out structure

  f_out.write ('''
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine zzz_to_c (f_zzz, c_zzz) bind(c)
!
! Routine to convert a Bmad zzz_struct to a C++ C_zzz
!
! Input:
!   f_zzz -- type (c_ptr), value :: Input Bmad zzz_struct structure.
!
! Output:
!   c_zzz -- type(c_ptr), value :: Output C++ C_zzz struct.
!-

subroutine zzz_to_c (f_zzz, c_zzz)

implicit none

interface
'''.replace('zzz', s_name))

  f_out.write ('  subroutine zzz_to_c2 (c_zzz'.replace('zzz', s_name) + (', ').join([a.c_name for a in args]) + ')\n')
  f_out.write ('  import fortran_and_cpp_mod\n')
  f_out.write ('  type (c_ptr), value :: c_zzz\n'.replace('zzz', s_name) + '\n')
  for arg in args:
    f_out.write ('  type (c_ptr), value :: ' + arg + '\n')
  f_out.write ('''
  end subroutine
end interface

type (zzz_struct), pointer :: f_zzz
type (c_ptr), value :: c_zzz
'''.replace('zzz', s_name))

  for arg in args:
    if arg.num == 1:
      f_out.write ('type (c_ptr), value :: ' + arg.c_name + '(' + arg.num + ')\n')
    else:
      f_out.write ('type (c_ptr), value :: ' + arg.c_name + '\n')

  f_out.write ('\n')




  f_out.write ('''

!

c_real_arr = [f_zzz%, ...]
c_in_arr = [f_zzz%, ...]

call zzz_to_c2 (c_zzz, c_real_arr, c_int_arr)

end subroutine zzz_to_c
'''.replace('zzz', s_name))

########################
# End stuff

f_out.write('end module\n')

##################################################################################
##################################################################################
# Create C++ class

##################################################################################
##################################################################################
# Create C++ side of interface

##################################################################################
##################################################################################
# Create Fortran struct equality check code

##################################################################################
##################################################################################
# Create C++ class equality check code

##################################################################################
##################################################################################
# Create interface check code

