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
    self.name = ''     # Fortran name. C name is 'C_<name>'
    self.var = []      # List of structrure components. Array of var_class.
    self.args = []     # Argument list

  def __repr__(self):
    return '[name: %s, #var: %i]' % (self.name, len(self.var))

# var_class

class var_class:

  def __init__(self):
    self.name = ''
    self.type = ''
    self.c_type = ''
    self.pointer_type = '-'    # 'pointer', 'allocatable'
    self.array = ''            # [':', ':'] or ['6']
    self.init_value = '-'
    self.comment = ''
    self.n_array = 1

  def __repr__(self):
    return '[%s, %s, %s, (%s), %s]' % (self.type, self.pointer_type, self.name, self.array, self.init_value)

# arg_class

class arg_class:
  def __init__(self, name):
    self.name = name
    self.c_type = ''
    self.f_type = ''
    self.num = 1

  def __repr__(self):
    return '[%s, %s]' % (self.name, self.num)


##################################################################################
##################################################################################
# Translations

T = True
F = False

class f_side_trans_class:

  def __init__(self, to_c_arg, bindc_type, bindc_name, to_f2_trans):
    self.to_c_arg = to_c_arg
    self.bindc_type = bindc_type
    self.bindc_name = bindc_name
    self.to_f2_trans = to_f2_trans

  def __repr__(self):
    return '%s,  %s,  %s :: %s' % (self.to_c_arg, self.bindc_type, self.bindc_name, self.to_f2_trans)


#                  Dim  Fix                       to_c_arg              bindc_type                    bindc_name  to_f2_trans
f_side_trans = {
  ('real(rp)',      0,   T) : f_side_trans_class('name',               'real(c_double)',             'name',      'f_zzz%name = name'),
  ('real(rp)',      1,   T) : f_side_trans_class('name',               'real(c_double)',             'name(*)',   'f_zzz%name = name'),
  ('real(rp)',      2,   T) : f_side_trans_class('mat2arr(name)',      'real(c_double)',             'name(*)',   'f_zzz%name = arr2mat(name, name_dim)'),
  ('real(rp)',      3,   T) : f_side_trans_class('tensor2arr(name)',   'real(c_double)',             'name(*)',   'f_zzz%name = arr2tensor(name, name_dim)'),
  ('complex(rp)',   0,   T) : f_side_trans_class('name',               'complex(c_double_complex)',  'name',      'f_zzz%name = name'),
  ('complex(rp)',   1,   T) : f_side_trans_class('name',               'complex(c_double_complex)',  'name(*)',   'f_zzz%name = name'),
  ('integer',       0,   T) : f_side_trans_class('name',               'int(c_int)',                 'name',      'f_zzz%name = name'),
  ('integer',       1,   T) : f_side_trans_class('name',               'int(c_int)',                 'name(*)',   'f_zzz%name = name'),
  ('integer',       2,   T) : f_side_trans_class('imat2arr(name)',     'int(c_int)',                 'name(*)',   'f_zzz%name = arr2imat(name, name_dim)'),
  ('logical',       0,   T) : f_side_trans_class('c_logic(name)',      'int(c_bool)',                'name',      'f_zzz%name = f_logic(name)'),
  ('logical',       1,   T) : f_side_trans_class('c_logic(name)',      'int(c_bool)',                'name(*)',   'f_zzz%name = f_logic(name)'),
  ('character(n)',  0,   T) : f_side_trans_class('c_str(name)',        'character(c_char)',          'name(*)',   'f_zzz%name = f_str(name)'),
  ('type:zzz',      0,   T) : f_side_trans_class('c_loc(name)',        'type(c_ptr), value',         'name',      'call zzz_to_f(name, f_zzz%name)', ),
  ('type:zzz',      1,   T) : f_side_trans_class('c_loc(name)',        'type(c_ptr), value',         'name(*)',   'do i = 1, name_dim; call zzz_to_f(name(i), f_zzz%name(i)); enddo')}

#############################################################

class pointer_trans_class:

  def __init__(self, type, def, arg, to_c):
    self.type = type
    self.def = def
    self.arg = arg
    self.to_c = to_c

  def __repr__(self):
    return '%s,  %s,  %s,  %s' % (self.type, self.def, self.arg, self.to_c)

# Dim                      type         def             arg           to_c
Pointer_trans = {
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



#                  Dim  Fix                        c_class          to_f2_arg            to_f2_call             to_c2_arg
c_side_trans = {
  ('real(rp)',      0,    T) : c_side_trans_class('double',        'Re&',               'c_zzz.name',          'Re& name'),
  ('real(rp)',      1,    T) : c_side_trans_class('Real_Array',    'ReArr',             '&c_zzz.name[0]',      'ReArr name'),
  ('real(rp)',      2,    T) : c_side_trans_class('Real_Matrix',   'ReArr',             'c_zzz.name',          'ReArr name'),
  ('complex(rp)',   0,    T) : c_side_trans_class('Complx',        'CComplx',           'c_zzz.name',          'CComplx name'),
  ('complex(rp)',   1,    T) : c_side_trans_class('Complx_Array',  'const ComplxArr',   'c_zzz.name',          'const ComplxArr name'),
  ('integer',       0,    T) : c_side_trans_class('int',           'Int&',              'c_zzz.name',          'Int& name'),
  ('integer',       1,    T) : c_side_trans_class('Int_Array',     'IntArr',            'c_zzz.name',          'intarr name'),
  ('integer',       2,    T) : c_side_trans_class('Int_Matrix',    'IntArr',            'c_zzz.name',          'IntArr name'),
  ('logical',       0,    T) : c_side_trans_class('bool',          'Int',               'c_zzz.name',          'Int& name'),
  ('logical',       1,    T) : c_side_trans_class('Bool_Array',    'IntArr',            'c_zzz.name',          'IntArr name'),
  ('character(n)',  0,    T) : c_side_trans_class('string',        'Char',              'c_zzz.name',          'Char name'),
  ('type:zzz',      0,    T) : c_side_trans_class('C_zzz',         'const C_zzz&',      '&c_zzz.name',         'const C_zzz&'),
  ('type:zzz',      1,    T) : c_side_trans_class('C_zzz_array',   'const C_zzz&',      'c_zzz.name[0]',       'const C_zzz&')}

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
    struct_def[struct_name].name = struct_name

    # Now collect the struct variables

    for line in f_module_file:
      var = var_class()

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

        struct_def[struct_name].var.append(var)
        
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
# Sort the variables into the argument list

for struct in struct_def:

  s_name = struct[0:-7]  # Strip off ending "_struct"

  # Form argument list

  for var in struct_def[struct]:

    struct.args.append(arg_class('var.name'))
    arg = args[-1]

    if var.pointer_type == 'allocatable' or var.pointer_type == 'pointer':

    elif var.type == 'character':

    elif var.array /= '':
      for d1 in var.array.split(','):
        d1_apart = d1.partition(':')
        if d1_apart[2] == '':
          arg.num = arg.num * int(d1_apart[0])
        else:
          arg.num = arg.num * (int(d1_apart[2]) - int(d1_apart[0]))

    elif var.name[0:5] == 'type:':

    elif var.type == 'real(rp)':
      arg.f_type = 'real(c_double)'
      arg.c_type = 'double'

    elif var.type == 'complex(rp)':
      arg.f_type = 'complex(c_double_complex)'
      arg.c_type = 'Complx'

    elif var.type == 'integer':
      arg.f_type = 'integer(c_int)'
      arg.c_type = 'int'

    elif var.type == 'logical':
      arg.f_type = 'integer(c_bool)'
      arg.c_type = 'int'

    else:
      print 'CONVERSION TO ARGUMENT FAILED:' + var.name

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
use intrinsic :: iso_c_binding

contains

''')


##############
# Loop over all the structures...

for struct in struct_def:

  args = struct.args

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

  f_face.write ('  subroutine zzz_to_c2 (c_zzz'.replace('zzz', s_name) + (', ').join([a.f_name for a in args]) + ')\n')
  f_face.write ('  import c_ptr\n')
  f_face.write ('  type (c_ptr), value :: c_zzz\n'.replace('zzz', s_name) + '\n')
  for arg in args:
    f_face.write ('  type (c_ptr), value :: ' + arg + '\n')
  f_face.write ('''
  end subroutine
end interface

type (zzz_struct), pointer :: f_zzz
type (c_ptr), value :: c_zzz
'''.replace('zzz', s_name))

  for arg in args:
    if arg.num == 1:
      f_face.write ('type (c_ptr), value :: ' + arg.c_name + '(' + arg.num + ')\n')
    else:
      f_face.write ('type (c_ptr), value :: ' + arg.c_name + '\n')

  f_face.write ('\n')


  f_face.write ('''

!

c_real_arr = [f_zzz%, ...]
c_in_arr = [f_zzz%, ...]

call zzz_to_c2 (c_zzz, c_real_arr, c_int_arr)

end subroutine zzz_to_c
'''.replace('zzz', s_name))


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
'''

for i in range(0, len(struct_def), 5):
  f_equ.write ('  module procedure ' + ', '.join('eq_' + f.name for f in struct_def[i:i+5] + '\n')

f_equ.write ('''
end interface

contains
'''

for struct in struct_def:
  f_equ.write ('''
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

elemental function eq_zzz (f1, f2) result (is_eq)

implicit none

type (zzz_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = .true.
'''.replace('zzz', struct.name)

  for var in struct.var:
    f_equ.write ('is_eq = is_eq .and. (f1%xxx == f2%xxx)'.replace('xxx', var.name)

  f_equ.write ('''
end function eq_zzz
'''.replace('zzz', struct.name)
  
f_equ.write ('end module')
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
'''

for struct in struct_def:
  f_test.write ('call test1_f_zzz(ok); if (.not. ok) all_ok = .false.\n')

f_test.write('''
if (all_ok) then
  print *, 'Bottom Line: Everything OK!'
else
  print *, 'BOTTOM LINE: PROBLEMS FOUND!'

end program
'''

f_test.close()

##################################################################################
##################################################################################
# Create Fortran side check code

f_ftest = open('test/bmad_cpp_test_mod.f90', 'w')
f_ftest.write('''
module bmad_cpp_test_mod

use bmad_cpp_convert_mod

contains

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
call set_zzz_test_pattern (f_zzz, 1)

call test_c_zzz(c_loc(f_zzz), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_zzz_test_pattern (f2_zzz, 1)
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

call set_zzz_test_pattern (f2_zzz, 2)
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

subroutine set_zzz_test_pattern (f_zzz, ix_patt)

implicit none

type (zzz_struct) f_zzz
integer ix_patt

!

'''.replace('zzz', struct.name))

for i, var in enumerate(struct.var):
  f_ftest.write ('f_' + struct.name + '%' + var.name + ' = ' + str(i) + ' * ix_patt\n' 

f_ftest.write('''
end subroutine set_zzz_test_pattern
'''.replace('zzz', struct.name))

f_ftest.write ('''
end module
'''

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

for struct in struct_def:
  f_class.write('''
//--------------------------------------------------------------------
// zzz

class zzz_struct {};

class C_zzz {
public:
'''.replace('zzz', struct.name))

  for var in struct:
    f_class.write('  ' + f_type_to_cpp(var.type) + ' ' + var.name  + ';\n')

  f_class.write ('''
  };   // End Class

  extern "C" void zzz_to_c (zzz_struct*, C_zzz&);
  extern "C" void zzz_to_f (C_zzz&, zzz_struct*);

  typedef valarray<C_zzz>    C_zzz_array;

  '''.replace('zzz', struct.name))

f_class.write('''
#define CPP_AND_BMAD
#endif
'''

f_class.close()

##################################################################################
##################################################################################
# Create C++ side of interface

f_face = open('code/cpp_bmad_convert.cpp', 'w')
f_face.write('''
#include "cpp_bmad_classes.h"
#include <iostream>
''')

for struct in struct_def:
  f_face.write('''
//--------------------------------------------------------------------
//--------------------------------------------------------------------

extern "C" void zzz_to_c (zzz_struct*, C_zzz&);
extern "C" void zzz_to_f2 (zzz_struct*, IntArr);

extern "C" void zzz_to_f (C_zzz& c_zzz, zzz_struct* f_zzz) {
'''.replace('zzz', struct.name))

  for arg in struct.args:


  f_face.write('''
  zzz_to_f2 (f_zzz, c_int);
}

extern "C" void zzz_to_c2 (C_zzz& c_zzz, '''.replace('zzz', struct.name))



f_face.write ('}')


##################################################################################
##################################################################################
# Create C++ class equality check code

##################################################################################
##################################################################################
# Create C++ side code check

