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
# struct_def_class
# Class for a structure

class struct_def_class:
  def __init__(self):
    self.short_name = '' # Struct name without trailing '_struct'. Note: C name is 'C_<name>'
    self.var = []        # List of structrure components. Array of var_class.
    self.dim_var = []    # Array dimensions

  def __repr__(self):
    return '[name: %s, #var: %i]' % (self.short_name, len(self.var))

# var_class

class var_class:

  def __init__(self):
    self.name = ''             # Name of variable
    self.type = ''             # Fortran type. EG: 'real(rp)', 'type'
    self.sub_type = ''         # For: components that are structures: the structure name. 
                               #      character strings: The number of characters
    self.pointer_type = '-'    # '-', 'PTR', 'ALLOC'
    self.array = []            # EG: [':', ':'] or ['6']
    self.init_value = '-'      # Initialization value
    self.comment = ''          # Comment with Fortran structure def.
    self.f_side = 0
    self.pointer_trans = 0
    self.c_side = 0

  def __repr__(self):
    return '[%s, %s, %s, [%s], "%s"]' % (self.type, self.pointer_type, self.name, self.array, self.init_value)

##################################################################################
##################################################################################
# Translations

NOT = '-'
PTR = 'PTR'
ALLOC = 'ALLOC'
T = True
F = False

class f_side_trans_class:

  def __init__(self, to_c2_arg, bindc_type, bindc_name, to_f2_trans):
    self.to_c2_arg = to_c2_arg
    self.bindc_type = bindc_type
    self.bindc_name = bindc_name
    self.bindc_const = ''
    self.to_f2_trans = to_f2_trans

  def __repr__(self):
    return '%s,  %s,  %s :: %s' % (self.to_c2_arg, self.bindc_type, self.bindc_name, self.to_f2_trans)


#                 Dim  Has_P                     to_c2_arg             bindc_type                    bindc_name   to_f2_trans
f_side_trans = {
  ('real(rp)',     0,   F) : f_side_trans_class('name',               'real(c_double)',             'name',      'f_zzz%name = name'),
  ('real(rp)',     1,   F) : f_side_trans_class('name',               'real(c_double)',             'name(*)',   'f_zzz%name = name'),
  ('real(rp)',     2,   F) : f_side_trans_class('mat2arr(name)',      'real(c_double)',             'name(*)',   'f_zzz%name = arr2mat(name, name_dim)'),
  ('real(rp)',     3,   F) : f_side_trans_class('tensor2arr(name)',   'real(c_double)',             'name(*)',   'f_zzz%name = arr2tensor(name, name_dim)'),
  ('complex(rp)',  0,   F) : f_side_trans_class('name',               'complex(c_double_complex)',  'name',      'f_zzz%name = name'),
  ('complex(rp)',  1,   F) : f_side_trans_class('name',               'complex(c_double_complex)',  'name(*)',   'f_zzz%name = name'),
  ('integer',      0,   F) : f_side_trans_class('name',               'int(c_int)',                 'name',      'f_zzz%name = name'),
  ('integer',      1,   F) : f_side_trans_class('name',               'int(c_int)',                 'name(*)',   'f_zzz%name = name'),
  ('integer',      2,   F) : f_side_trans_class('imat2arr(name)',     'int(c_int)',                 'name(*)',   'f_zzz%name = arr2imat(name, name_dim)'),
  ('logical',      0,   F) : f_side_trans_class('c_logic(name)',      'int(c_bool)',                'name',      'f_zzz%name = f_logic(name)'),
  ('logical',      1,   F) : f_side_trans_class('c_logic(name)',      'int(c_bool)',                'name(*)',   'f_zzz%name = f_logic(name)'),
  ('character',    0,   F) : f_side_trans_class('c_str(name)',        'character(c_char)',          'name(*)',   'f_zzz%name = f_str(name)'),
  ('type',         0,   F) : f_side_trans_class('c_loc(name)',        'type(c_ptr), value',         'name',      'call zzz_to_f(name, f_zzz%name)', ),
  ('type',         1,   F) : f_side_trans_class('c_loc(name)',        'type(c_ptr), value',         'name(*)',   'do i = 1, name_dim; call zzz_to_f(name(i), f_zzz%name(i)); enddo')}

for key, f in f_side_trans.items(): f.bindc_const = f.bindc_type.partition('(')[2].partition(')')[0]

#############################################################

class pointer_trans_class:

  def __init__(self, type, name, arg, to_c):
    self.type = type
    self.name = name
    self.arg = arg
    self.to_c = to_c

  def __repr__(self):
    return '%s,  %s,  %s,  %s' % (self.type, self.name, self.arg, self.to_c)

# Dim                      type         name            arg           to_c
pointer_trans = {
  0 : pointer_trans_class('integer',   'name_dim',     'name_dim',   'name_dim = 0; if (associated(f_zzz%name)) name_dim = 1'),
  1 : pointer_trans_class('integer',   'name_dim',     'name_dim',   'name_dim = 0; if (associated(f_zzz%name)) name_dim = size(f_zzz%name)'),
  2 : pointer_trans_class('integer',   'name_dim(2)',  'name_dim',   'name_dim = 0; if (associated(f_zzz%name)) name_dim = size(f_zzz%name)'),
  3 : pointer_trans_class('integer',   'name_dim(3)',  'name_dim',   'name_dim = 0; if (associated(f_zzz%name)) name_dim = size(f_zzz%name)')}

#############################################################

class c_side_trans_class:

  def __init__(self, c_class, to_f2_arg, to_f2_call, to_c2_arg):
    self.c_class = c_class
    self.to_f2_arg = to_f2_arg
    self.to_f2_call = to_f2_call
    self.to_c2_arg = to_c2_arg

  def __repr__(self):
    return '%s,  %s,  %s,  %s' % (self.c_class, self.to_f2_arg, self.to_f2_call, self.to_c2_arg)



#                  Dim  Has_P                      c_class          to_f2_arg            to_f2_call             to_c2_arg
c_side_trans = {
  ('real(rp)',     0,    F) : c_side_trans_class('double',        'Re&',               'c_zzz.name',          'Re& name'),
  ('real(rp)',     1,    F) : c_side_trans_class('Real_Array',    'ReArr',             '&c_zzz.name[0]',      'ReArr name'),
  ('real(rp)',     2,    F) : c_side_trans_class('Real_Matrix',   'ReArr',             'c_zzz.name',          'ReArr name'),
  ('complex(rp)',  0,    F) : c_side_trans_class('Complx',        'CComplx',           'c_zzz.name',          'CComplx name'),
  ('complex(rp)',  1,    F) : c_side_trans_class('Complx_Array',  'const ComplxArr',   'c_zzz.name',          'const ComplxArr name'),
  ('integer',      0,    F) : c_side_trans_class('int',           'Int&',              'c_zzz.name',          'Int& name'),
  ('integer',      1,    F) : c_side_trans_class('Int_Array',     'IntArr',            'c_zzz.name',          'intarr name'),
  ('integer',      2,    F) : c_side_trans_class('Int_Matrix',    'IntArr',            'c_zzz.name',          'IntArr name'),
  ('logical',      0,    F) : c_side_trans_class('bool',          'Int',               'c_zzz.name',          'Int& name'),
  ('logical',      1,    F) : c_side_trans_class('Bool_Array',    'IntArr',            'c_zzz.name',          'IntArr name'),
  ('character',    0,    F) : c_side_trans_class('string',        'Char',              'c_zzz.name',          'Char name'),
  ('type',         0,    F) : c_side_trans_class('C_zzz',         'const C_zzz&',      '&c_zzz.name',         'const C_zzz&'),
  ('type',         1,    F) : c_side_trans_class('C_zzz_array',   'const C_zzz&',      'c_zzz.name[0]',       'const C_zzz&')}

##################################################################################
##################################################################################
# Get the list of structs

f_struct_list_file = open('scripts/fortran_structs.list')

struct_def = {}
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
    struct_def[line.split()[0]] = struct_def_class()

f_struct_list_file.close()

##################################################################################
##################################################################################
# Parse structure definitions

# Current restrictions: Syntax to avoid:
#   Multiple inits: "real a = 5, b = 7"
#   No space in type parens: "type( abc_struct )" or "real( rp )"
#   Space in array def: "arr( : )"

for file_name in f_module_files:
  f_module_file = open('../' + file_name)

  for line in f_module_file:
    split_line = line.split()
    if len(split_line) < 2: continue
    if split_line[0] != 'type': continue
    if not split_line[1] in struct_def: continue

    struct_name = split_line[1]
    struct = struct_def[struct_name]
    struct.short_name = struct_name[:-7]   # Remove '_struct' suffix
    
    # Now collect the struct variables

    for line in f_module_file:

      var = var_class()

      part = line.partition('!')

      var.comment = part[2].strip()

      part = part[0].replace('  ', ' ').replace(' (', '(').partition('=')
      var.init_value = part[2].strip()
      if var.init_value[:1] == '>': var.init_value = ''  # Ignore pointer "=> null()" init.

      split_line = part[0].replace(',', ' ').replace('::', ' ').split()
      if len(split_line) == 0: continue
      if split_line[0:2] == ['end', 'type']: break

      var.type = split_line.pop(0)
      if split_line[0][0:1] == '(': var.type = var.type + split_line.pop(0)

      if var.type == 'character(':
        var.sub_type = var.type[10:-1]
        var.type = 'character'

      if var.type == 'type(':
        var.sub_type = var.type[5:-1]
        var.type = 'type' 

      if split_line[0] == 'allocatable':
        var.pointer_type = ALLOC
        split_line.pop(0)

      if split_line[0] == 'pointer':
        var.pointer_type = PTR
        split_line.pop(0)

      # Handle lines like: "real(rp) kx, ky, kz"

      for ix, name in enumerate(split_line):
        var = copy.copy(var)
        var.name = name

        part = var.name.partition('(')
        var.name = part[0]
        var.array = part[2][:-1]     # Remove last ')'

        if ix+1 < len(split_line) and split_line[ix+1][0:1] == '(': 
          var.array = split_line.pop(ix+1)[:-1]

        # F side translation

        n_dim = len(var.array)
        has_p = F     ### has_p = (var.pointer_type != '-')

        if (var.type, n_dim, has_p) not in f_side_trans:
          print 'NO TRANSLATION FOR: ' + struct.short_name + '%' + var.name
          continue

        var.f_side = f_side_trans[var.type, n_dim, has_p]
        var.c_side = c_side_trans[var.type, n_dim, has_p]

        # Dimensionality

        var.pointer_trans = pointer_trans[n_dim]

        struct.var.append(var)        

        if len(var.array) != 0:
          if var.array[0] == ':':
            for n in range(1, len(var.array)):
              dim_var = var_class()
              dim_var.name = 'n' + str(n) + '_' + var.name
              dim_var.type = 'integer'
              dim_var.f_side = f_side_class('integer', 0, F)
              struct.dim_var.append(dim_var)

  f_module_file.close()

##################################################################################
##################################################################################
# As a check, write results to file. 

f_out = open('f_structs.parsed', 'w')

for struct in struct_def:
  f_out.write('******************************************\n')
  f_out.write (struct + '    ' + str(len(struct_def[struct].var)) + '\n')
  for var in struct_def[struct].var:
    f_out.write ('    ' + str(var) + '\n')

f_out.close()

##################################################################################
##################################################################################
# Create Fortran side of interface...

# First the header

f_face = open('code/bmad_cpp_convert_mod.f90', 'w')

f_face.write ('''
!+
! Fortran side of the Bmad / C++ interface.
!
! File Generated by: create_interface.py
! Do not edit this file directly! 
!-

module bmad_cpp_convert_mod

use bmad_struct
use bmad_interface
use fortran_and_cpp_mod
use, intrinsic :: iso_c_binding

contains
''')


##############
# Loop over all the structures...

for struct_key, struct in struct_def.items():

  s_name = struct.short_name

  f_face.write ('''
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

  import_set = set(['c_ptr'])
  to_c2_arg_def = {'type (c_ptr), value' : ['c_' + s_name]}

  for var in struct.var: 
    f_side = var.f_side
    import_set.add(f_side.bindc_const)
    if not f_side.bindc_type in to_c2_arg_def: to_c2_arg_def[f_side.bindc_type] = []
    to_c2_arg_def[f_side.bindc_type].append(f_side.bindc_name.replace('name', var.name))

  for dim_var in struct.dim_var: 
    f_side = dim_var.f_side
    import_set.add(f_side.bindc_const)
    if not f_side.bindc_type in to_c2_arg_def: to_c2_arg_def[f_side.bindc_type] = []
    to_c2_arg_def[f_side.bindc_type].append(f_side.bindc_name.replace('name', var.name))

  f_face.write ('  subroutine zzz_to_c2 (c_zzz'.replace('zzz', s_name))
  for var in struct.var: f_face.write (', ' + var.name)
  for var in struct.dim_var: f_face.write (', ' + var.name)
  f_face.write (')\n')
  f_face.write ('    import ' + ', '.join(import_set) + '\n')

  for arg_type, args in to_c2_arg_def.items():
    f_face.write ('    ' + arg_type + ' :: ' + ', '.join(args) + '\n')

  f_face.write ('''  end subroutine
end interface

type (zzz_struct), pointer :: f_zzz
type (c_ptr), value :: c_zzz

call zzz_to_c2 (c_zzz'''.replace('zzz', s_name))

  for var in struct.var:
    f_face.write (', ' + var.f_side.to_c2_arg.replace('name', var.name))

  f_face.write(''')

end subroutine zzz_to_c
'''.replace('zzz', struct.short_name))

########################
# End stuff

f_face.write('end module\n')
f_face.close()

##################################################################################
##################################################################################
# Create Fortran struct equality check code

f_equ = open('code/bmad_equality.f90', 'w')

f_equ.write ('''
module bmad_equality

use bmad_struct

interface operator (==)
''')

for i in range(0, len(struct_def), 5):
  f_equ.write ('  module procedure ' + ', '.join('eq_' + 
                              f.short_name for f in struct_def.values()[i:i+5]) + '\n')

f_equ.write ('''end interface

contains
''')

for key, struct in struct_def.items():
  f_equ.write ('''
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

elemental function eq_zzz (f1, f2) result (is_eq)

implicit none

type (zzz_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = .true.
'''.replace('zzz', struct.short_name))

  for var in struct.var:
    f_equ.write ('is_eq = is_eq .and. (f1%xxx == f2%xxx)'.replace('xxx', var.name) + '\n')

  f_equ.write ('''
end function eq_zzz
'''.replace('zzz', struct.short_name))
  
f_equ.write ('end module\n')
f_equ.close()

##################################################################################
##################################################################################
# Create code check main program

f_test = open('test/bmad_cpp_test.f90', 'w')
f_test.write('''
program bmad_cpp_test

use bmad_cpp_test_mod

logical ok, all_ok

!

all_ok = .true.
''')

for struct in struct_def.values():
  f_test.write ('call test1_f_' + struct.short_name + '(ok); if (.not. ok) all_ok = .false.\n')

f_test.write('''
if (all_ok) then
  print *, 'Bottom Line: Everything OK!'
else
  print *, 'BOTTOM LINE: PROBLEMS FOUND!'
endif

end program
''')

f_test.close()

##################################################################################
##################################################################################
# Create Fortran side check code

f_ftest = open('test/bmad_cpp_test_mod.f90', 'w')
f_ftest.write('''
module bmad_cpp_test_mod

use bmad_cpp_convert_mod

contains
''')

for struct_key, struct in struct_def.items():
  f_ftest.write('''
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_zzz (ok)

implicit none

type (zzz_struct), target :: f_zzz, f2_zzz
integer c_ok
logical ok

interface
  subroutine test_c_zzz (c_zzz, c_ok) bind(c)
    import c_ptr, c_int
    type (c_ptr), value :: c_zzz
    integer(c_int) c_ok
  end subroutine
end interface

!

ok = .true.
f_zzz = zzz_test_pattern (1)

call test_c_zzz(c_loc(f_zzz), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

f2_zzz = zzz_test_pattern (1)
if (f_zzz == f2_zzz) then
  print *, zzz: C_side_convert C->F: Good'
else
  print *, zzz: C_side_convert C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_zzz

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_zzz (c_zzz, c_ok) bind(c)

implicit  none

type (c_ptr), value ::  c_zzz
type (zzz_struct), target :: f_zzz, f2_zzz
integer(c_int) c_ok

!

call zzz_to_f (c_zzz, c_loc(f_zzz))

f2_zzz = zzz_test_pattern (2)
if (f_zzz == f2_zzz) then
  print *, zzz: F_side_convert C->F: Good'
else
  print *, zzz: F_side_convert C->F: FAILED!'
  ok = .false.
endif

call zzz_to_c (c_loc(f2_zzz), c_zzz)

end subroutine test2_f_zzz

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

Function zzz_test_pattern (ix_patt) result (f_zzz)

implicit none

type (zzz_struct) f_zzz
integer ix_patt, offset

!

offset = 100 * ix_patt
'''.replace('zzz', struct.short_name))

for i, var in enumerate(struct.var, 1):
  short_type = var.type.partition('(')[0]  # EG: 'real(rp)' -> 'real'
  if var.pointer_type == 'PTR':
  if (len(var.array)
  if var.type == 'type':
  elif var.type == 'character':
  else

  f_ftest.write ('f_' + struct.short_name + '%' + var.name + ' = ' + str(i) + ' + offset \n')

f_ftest.write('''
end function zzz_test_pattern
'''.replace('zzz', struct.short_name))

f_ftest.write ('''
end module
''')

f_ftest.close()

##################################################################################
##################################################################################
# Create C++ class

f_class = open('include/cpp_bmad_classes.h', 'w')
f_class.write('''

#ifndef CPP_BMAD_CLASSES

#include <string>
#include <string.h>
#include <valarray>
#include <complex>
#include "bmad_parameters.h"

using namespace std;

typedef const double    Re;
typedef const int       Int;
typedef const char*     Char;
typedef const bool      Bool;
typedef const double*   ReArr;
typedef const int*      IntArr;

typedef complex<double>                 Complx;
typedef const complex<double>           CComplx;

typedef valarray<double>                Real_Array;
typedef valarray<Complx>                Complx_Array;
typedef valarray<bool>                  Bool_Array;
typedef valarray<int>                   Int_Array;

typedef valarray<Real_Array>            Real_Matrix;
typedef valarray<Bool_Array>            Bool_Matrix;

typedef valarray<Real_Matrix>           Real_Tensor;

const Real_Array V2_array(double(0), 2);
const Real_Array V3_array(double(0), 3);
const Real_Array V6_array(double(0), 6);
const Real_Matrix M2_mat(V2_array, 2);
const Real_Matrix M3_mat(V3_array, 3);
const Real_Matrix M6_mat(V6_array, 6);
const double       V0(0);
const Complx       C0(0);

''')

for s_name, struct in struct_def.items():
  f_class.write('''
//--------------------------------------------------------------------
// zzz

class zzz_struct {};

class C_zzz {
public:
'''.replace('zzz', s_name))

  for var in struct.var:
    f_class.write('  ' + var.c_side.c_class + ' ' + var.name  + ';\n')

  f_class.write ('''
  };   // End Class

  extern "C" void zzz_to_c (zzz_struct*, C_zzz&);
  extern "C" void zzz_to_f (C_zzz&, zzz_struct*);

  typedef valarray<C_zzz>    C_zzz_array;

  '''.replace('zzz', struct.short_name))

f_class.write('''
#define CPP_AND_BMAD
#endif
''')

f_class.close()

##################################################################################
##################################################################################
# Create C++ side of interface

f_face = open('code/cpp_bmad_convert.cpp', 'w')
f_face.write('''
#include "cpp_bmad_classes.h"
#include <iostream>
''')

for s_name, struct in struct_def.items():
  f_face.write('''
//--------------------------------------------------------------------
//--------------------------------------------------------------------

extern "C" void zzz_to_c (zzz_struct*, C_zzz&);
extern "C" void zzz_to_f2 (zzz_struct*, IntArr);

extern "C" void zzz_to_f (C_zzz& c_zzz, zzz_struct* f_zzz) {
'''.replace('zzz', struct.short_name))


  f_face.write('''
  zzz_to_f2 (f_zzz, c_int);
}

extern "C" void zzz_to_c2 (C_zzz& c_zzz, '''.replace('zzz', struct.short_name))



f_face.write ('}')


##################################################################################
##################################################################################
# Create C++ class equality check code

##################################################################################
##################################################################################
# Create C++ side code check

