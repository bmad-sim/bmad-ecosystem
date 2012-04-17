#!/usr/bin/env python

# Note: Run this script in the cpp_bmad_interface directory.

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
    self.type = ''             # Fortran type. EG: 'real(rp)', 'type', 'character'
    self.full_type = ''        # Same as type except: Structures: 'type(<struct_name>)'
                               #                      Strings: 'character(<len>)'
    self.pointer_type = '-'    # '-', 'PTR', 'ALLOC'
    self.array = []            # EG: [':', ':'] or ['0:6', '3']
    self.full_array = ''       # EG: '(:,:)', '(0:6, 3)'
    self.lbound = []
    self.ubound = []
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
REAL  = 'real(rp)'
CMPLX = 'complex(rp)'
INT   = 'integer'
LOGIC = 'logical'
CHAR  = 'character'
TYPE  = 'type'

# Fortran side translation

class f_side_trans_class:

  def __init__(self, to_c2_arg, bindc_type, bindc_name):
    self.to_c2_arg = to_c2_arg
    self.bindc_type = bindc_type
    self.bindc_name = bindc_name
    self.bindc_const = ''
    self.equal_test = '(f1%name == f2%name)'
    self.to_f2_trans = 'FP%name = name'
    self.test_pat = 'FF%name = XXX + offset'

  def __repr__(self):
    return '%s,  %s,  %s :: %s' % (self.to_c2_arg, self.bindc_type, self.bindc_name, self.to_f2_trans)


#       Dim  Has_P                    to_c2_arg             bindc_type                    bindc_name   equal_test
f_side_trans = {
  (REAL,  0, F) : f_side_trans_class('name',               'real(c_double)',             'name'),
  (REAL,  1, F) : f_side_trans_class('name',               'real(c_double)',             'name(*)'),
  (REAL,  2, F) : f_side_trans_class('mat2arr(name)',      'real(c_double)',             'name(*)'),
  (REAL,  3, F) : f_side_trans_class('tensor2arr(name)',   'real(c_double)',             'name(*)'),
  (CMPLX, 0, F) : f_side_trans_class('name',               'complex(c_double_complex)',  'name'),
  (CMPLX, 1, F) : f_side_trans_class('name',               'complex(c_double_complex)',  'name(*)'),
  (INT,   0, F) : f_side_trans_class('name',               'integer(c_int)',             'name'),
  (INT,   1, F) : f_side_trans_class('name',               'integer(c_int)',             'name(*)'),
  (INT,   2, F) : f_side_trans_class('imat2arr(name)',     'integer(c_int)',             'name(*)'),
  (LOGIC, 0, F) : f_side_trans_class('c_logic(name)',      'logical(c_bool)',            'name'),
  (LOGIC, 1, F) : f_side_trans_class('c_logic(name)',      'logical(c_bool)',            'name(*)'),
  (CHAR,  0, F) : f_side_trans_class('c_str(name)',        'character(c_char)',          'name(*)'),
  (TYPE,  0, F) : f_side_trans_class('c_loc(name)',        'type(c_ptr), value',         'name'),
  (TYPE,  1, F) : f_side_trans_class('c_loc(name)',        'type(c_ptr), value',         'name(*)')} 

f_side_trans[REAL,  2, F].to_f2_trans = 'FP%name = arr2tensor(name, name_dim)'
f_side_trans[INT,   2, F].to_f2_trans = 'FP%name = arr2imat(name, name_dim)'

f_side_trans[CMPLX, 1, F].test_pat = 'FF%name = [(cmplx(100 + jd1 + XXX + offset, 200 + jd1 + XXX + offset), jd1 = 1, size(FF%name))]'
f_side_trans[REAL,  1, F].test_pat = 'FF%name = [(100 + jd1 + XXX + offset, jd1 = 1, size(FF%name))]'
f_side_trans[REAL,  2, F].test_pat = 'FF%name = reshape([100 + jd1 + XXX + offset, jd1 = 1, size(FF%name)], shape(FF%name))'

for key, f in f_side_trans.items(): 
  f.bindc_const = f.bindc_type.partition('(')[2].partition(')')[0]
  if key[1] != 0: f.equal_test = 'all(f1%name == f2%name)'

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
    self.to_c2_set = 'C.name = name;'
    self.constructor = 'name(0)'
    self.equal_test = '(x.name == y.name)'
    self.test_pat = 'C.name = XXX + offset;'

  def __repr__(self):
    return '%s,  %s,  %s,  %s' % (self.c_class, self.to_f2_arg, self.to_f2_call, self.to_c2_arg)



#        Dim  Has_P                    c_class          to_f2_arg          to_f2_call             to_c2_arg
c_side_trans = {
  (REAL,   0,  F) : c_side_trans_class('double',         'Real&',           'C.name',          'Real& name'),
  (REAL,   1,  F) : c_side_trans_class('Real_Array',     'RealArr',         '&C.name[0]',      'RealArr name'),
  (REAL,   2,  F) : c_side_trans_class('Real_Matrix',    'RealArr',         'C.name',          'RealArr name'),
  (CMPLX,  0,  F) : c_side_trans_class('Dcomplex',       'Dcomplexd&',      'C.name',          'Dcomplex name'),
  (CMPLX,  1,  F) : c_side_trans_class('Dcomplex_Array', 'DcomplexArr',     '&C.name[0]',      'DcomplexArr name'),
  (INT,    0,  F) : c_side_trans_class('int',            'Int&',            'C.name',          'Int& name'),
  (INT,    1,  F) : c_side_trans_class('Int_Array',      'IntArr',          '&C.name[0]',      'IntArr name'),
  (INT,    2,  F) : c_side_trans_class('Int_Matrix',     'IntArr',          'C.name',          'IntArr name'),
  (LOGIC,  0,  F) : c_side_trans_class('bool',           'Bool&',           'C.name',          'Bool& name'),
  (LOGIC,  1,  F) : c_side_trans_class('Bool_Array',     'BoolArr',         'C.name',          'BoolArr name'),
  (CHAR,   0,  F) : c_side_trans_class('string',         'Char',            'C.name',          'Char name'),
  (TYPE,   0,  F) : c_side_trans_class('C_zzz',          'const C_zzz&',    '&C.name',         'const C_zzz& name'),
  (TYPE,   1,  F) : c_side_trans_class('C_zzz_array',    'const C_zzz&',    'C.name[0]',       'const C_zzz& name')}

c_side_trans[REAL,  1, F].constructor = 'name(0.0, DIM1)'
c_side_trans[REAL,  1, F].to_c2_set   = 'C.name = Real_Array(name, DIM1);'
c_side_trans[REAL,  1, F].test_pat    = 'for (int i = 0; i < C.name.size(); i++) C.name[i] = 101 + i + XXX + offset;'

c_side_trans[CMPLX, 1, F].constructor = 'name(0.0, DIM1)'
c_side_trans[CMPLX, 1, F].to_c2_set   = 'C.name = Dcomplex_Array(name, DIM1);'
c_side_trans[CMPLX, 1, F].test_pat    = 'for (int i = 0; i < C.name.size(); i++) C.name[i] = Dcomplex(101 + i + XXX + offset, 201 + i + XXX + offset);'

c_side_trans[INT,   1, F].constructor = 'name(DIM1)'
c_side_trans[INT,   1, F].to_c2_set   = 'C.name = Int_Array(name, DIM1);'
c_side_trans[INT,   1, F].test_pat    = 'for (int i = 0; i < C.name.size(); i++) C.name[i] = 101 + i + XXX + offset;'

c_side_trans[LOGIC, 1, F].constructor = 'name(DIM1)'
c_side_trans[LOGIC, 1, F].to_c2_set   = 'C.name = Int_Array(name, DIM1);'
c_side_trans[LOGIC, 1, F].test_pat    = 'for (int i = 0; i < C.name.size(); i++) C.name[i] = 101 + i + XXX + offset;'

for key, c in c_side_trans.items(): 
  if key[1] != 0: c.equal_test = 'is_all_equal(x.name, y.name)'

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

      var.full_type = split_line.pop(0)
      if split_line[0][0:1] == '(': var.full_type = var.full_type + split_line.pop(0)

      var.type = var.full_type

      if var.full_type[0:10] == 'character(':
        var.type = 'character'

      if var.full_type[0:5] == 'type(':
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

        if len(var.array) > 0 and var.array[0] != ':':
          for dim in var.array:
            if ':' in dim:
              var.lbound.append(dim.partition(':')[0])
              var.ubound.append(dim.partition(':')[2])
            else:
              var.lbound.append('1')
              var.ubound.append(dim)

        # F side translation

        n_dim = len(var.array)
        has_p = F     ### has_p = (var.pointer_type != '-')

        if (var.type, n_dim, has_p) not in f_side_trans:
          print 'NO TRANSLATION FOR: ' + struct.short_name + '%' + var.name + ' [', var.type + ', ' + str(n_dim) + ', ' + str(has_p) + ']'
          continue

        var.f_side = f_side_trans[var.type, n_dim, has_p]
        var.c_side = c_side_trans[var.type, n_dim, has_p]

        # Dimensionality

        var.pointer_trans = pointer_trans[n_dim]

        struct.var.append(var)        

        if len(var.array) != 0:
          var.full_array = '(' + ', '.join(var.array) + ')'
          if var.array[0] == ':':
            for n in range(1, len(var.array)):
              dim_var = var_class()
              dim_var.name = 'n' + str(n) + '_' + var.name
              dim_var.type = 'integer'
              dim_var.f_side = f_side_class('integer', 0, F)
              struct.dim_var.append(dim_var)

  f_module_file.close()

# Make some name substitutions

for struct in struct_def.values():
  for var in struct.var:
    if len(var.array) == 1: 
      var.c_side.constructor = var.c_side.constructor.replace('DIM1', str(1 + int(var.ubound[0]) - int(var.lbound[0])))
      var.c_side.to_c2_set = var.c_side.to_c2_set.replace('DIM1', str(1 + int(var.ubound[0]) - int(var.lbound[0])))

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

if not os.path.exists('code'): os.makedirs('code')
f_face = open('code/bmad_cpp_convert_mod.f90', 'w')

f_face.write ('''
!+
! Fortran side of the Bmad / C++ structure interface.
!
! File Generated by: create_interface.py
! Do not edit this file directly! 
!-

module bmad_cpp_convert_mod

use bmad_struct
use bmad_interface
use fortran_cpp_utils
use, intrinsic :: iso_c_binding
''')

##############
# zzz_to_f interface

for struct in struct_def.values():
  f_face.write ('''
!--------------------------------------------------------------------------

interface 
  subroutine zzz_to_f (CC, FF) bind(c)
    import c_ptr
    type (c_ptr), value :: CC, FF
  end subroutine
end interface
'''.replace('zzz', struct.short_name))


f_face.write ('contains\n')


##############
# zzz_to_c definitions

for struct_key, struct in struct_def.items():

  s_name = struct.short_name

  f_face.write ('''
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine zzz_to_c (FF, CC) bind(c)
!
! Routine to convert a Bmad zzz_struct to a C++ C_zzz structure
!
! Input:
!   FF -- type (c_ptr), value :: Input Bmad zzz_struct structure.
!
! Output:
!   CC -- type(c_ptr), value :: Output C++ C_zzz struct.
!-

subroutine zzz_to_c (FF, CC) bind(c)

implicit none

interface
'''.replace('zzz', s_name))

  import_set = set(['c_ptr'])
  to_c2_arg_def = {}

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

  f_face.write ('  subroutine zzz_to_c2 (CC'.replace('zzz', s_name))
  for var in struct.var: f_face.write (', ' + var.name)
  for var in struct.dim_var: f_face.write (', ' + var.name)
  f_face.write (') bind(c)\n')
  f_face.write ('    import ' + ', '.join(import_set) + '\n')
  f_face.write ('    type (c_ptr), value :: CC\n')
  for arg_type, args in to_c2_arg_def.items():
    f_face.write ('    ' + arg_type + ' :: ' + ', '.join(args) + '\n')

  f_face.write ('''  end subroutine
end interface

type (c_ptr), value :: FF
type (c_ptr), value :: CC
type (zzz_struct), pointer :: FP

!

call c_f_pointer (FF, FP)

call zzz_to_c2 (CC'''.replace('zzz', s_name))

  for var in struct.var:
    f_face.write (', FP%' + var.f_side.to_c2_arg.replace('name', var.name))

  f_face.write(''')

end subroutine zzz_to_c

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine zzz_to_f2 (FF, ...etc...) bind(c)
!
! Routine used in converting a C++ C_zzz structure to a Bmad zzz_struct structure.
! This routine is called by zzz_to_c and is not meant to be called directly.
!
! Input:
!   ...etc... -- Components of the structure. See the zzz_to_f2 code for more details.
!
! Output:
!   FF -- type (c_ptr), value :: Bmad zzz_struct structure.
!-

subroutine zzz_to_f2 (FF'''.replace('zzz', struct.short_name))

  for var in struct.var:
    f_face.write(', ' + var.name.replace('name', var.name))

  f_face.write(''') bind(c)\n

implicit none

type (c_ptr), value :: FF
type (zzz_struct), pointer :: FP
'''.replace('zzz', struct.short_name))

  f2_arg_list = {}
  for var in struct.var:
    if not var.full_type in f2_arg_list: f2_arg_list[var.full_type] = []
    f2_arg_list[var.full_type].append(var.name + var.full_array)

  for arg_type, arg_list in f2_arg_list.items():
    f_face.write(arg_type + ' ' + ', '.join(arg_list) + '\n')

  f_face.write('''
call c_f_pointer (FF, FP)
''')

  for var in struct.var:
    f_face.write (var.f_side.to_f2_trans.replace('name', var.name) + '\n')

  f_face.write ('''
end subroutine zzz_to_f2
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
    f_equ.write ('is_eq = is_eq .and. ' + var.f_side.equal_test.replace('name', var.name) + '\n')

  f_equ.write ('''
end function eq_zzz
'''.replace('zzz', struct.short_name))
  
f_equ.write ('end module\n')
f_equ.close()

##################################################################################
##################################################################################
# Create code check main program

if not os.path.exists('interface_test'): os.makedirs('interface_test')
f_test = open('interface_test/bmad_cpp_test.f90', 'w')

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

f_test = open('interface_test/bmad_cpp_test_mod.f90', 'w')
f_test.write('''
module bmad_cpp_test_mod

use bmad_cpp_convert_mod
use bmad_equality

contains
''')

for struct_key, struct in struct_def.items():
  f_test.write('''
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_zzz (ok)

implicit none

type (zzz_struct), target :: f_zzz, f2_zzz
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_zzz (c_zzz, c_ok) bind(c)
    import c_ptr, c_bool
    type (c_ptr), value :: c_zzz
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call zzz_test_pattern (f2_zzz, 1)

call test_c_zzz(c_loc(f2_zzz), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call zzz_test_pattern (f_zzz, 4)
if (f_zzz == f2_zzz) then
  print *, 'zzz: C side convert C->F: Good'
else
  print *, 'zzz: C side convert C->F: FAILED!'
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

c_ok = c_logic(.true.)
call zzz_to_f (c_zzz, c_loc(f_zzz))

call zzz_test_pattern (f2_zzz, 2)
if (f_zzz == f2_zzz) then
  print *, 'zzz: F side convert C->F: Good'
else
  print *, 'zzz: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call zzz_test_pattern (f2_zzz, 3)
call zzz_to_c (c_loc(f2_zzz), c_zzz)

end subroutine test2_f_zzz

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine zzz_test_pattern (FF, ix_patt)

implicit none

type (zzz_struct) FF
integer ix_patt, offset, jd1, jd2, jd3

!

offset = 100 * ix_patt

'''.replace('zzz', struct.short_name))

for i, var in enumerate(struct.var, 1):
  f_test.write (var.f_side.test_pat.replace('XXX', str(i)).replace('name', var.name) + '\n')

f_test.write('''
end subroutine zzz_test_pattern
'''.replace('zzz', struct.short_name))

f_test.write ('''
end module
''')

f_test.close()

##################################################################################
##################################################################################
# Create C++ class

if not os.path.exists('include'): os.makedirs('include')
f_class = open('include/cpp_bmad_classes.h', 'w')
f_class.write('''
//+
// C++ classes definitions for Bmad / C++ structure interface.
//
// File Generated by: create_interface.py
// Do not edit this file directly! 
//-

#ifndef CPP_BMAD_CLASSES

#include <string>
#include <string.h>
#include <valarray>
#include <complex>
#include "bmad_enums.h"

using namespace std;

typedef complex<double>          dcomplex;

typedef const bool               Bool;
typedef const dcomplex           Dcomplex;
typedef const char*              Char;
typedef const double             Real;
typedef const int                Int;

typedef const dcomplex*          DcomplexArr;
typedef const double*            RealArr;
typedef const int*               IntArr;

typedef valarray<bool>           Bool_Array;
typedef valarray<dcomplex>       Dcomplex_Array;
typedef valarray<double>         Real_Array;
typedef valarray<int>            Int_Array;

typedef valarray<Real_Array>     Real_Matrix;
typedef valarray<Bool_Array>     Bool_Matrix;

typedef valarray<Real_Matrix>    Real_Tensor;

''')

for key, struct in struct_def.items():
  f_class.write('''
//--------------------------------------------------------------------
// C_zzz

class zzz_struct {};  // Opaque class for pointers to corresponding fortran structs.

class C_zzz {
public:
'''.replace('zzz', struct.short_name))

  for var in struct.var:
    f_class.write('  ' + var.c_side.c_class.replace('zzz', struct.short_name) + ' ' + var.name  + ';\n')

  f_class.write ('''
  C_zzz() :
'''.replace('zzz', struct.short_name))

  construct_list = []
  for var in struct.var:
    construct_list.append(var.c_side.constructor.replace('name', var.name))

  f_class.write ('    ' + ', '.join(construct_list) + '\n')

  f_class.write('''    {}

};   // End Class

extern "C" void zzz_to_c (zzz_struct*, C_zzz&);
extern "C" void zzz_to_f (C_zzz&, zzz_struct*);

bool operator== (const C_zzz&, const C_zzz&);

typedef valarray<C_zzz>    C_zzz_array;

'''.replace('zzz', struct.short_name))

f_class.write('''
//--------------------------------------------------------------------

#define CPP_BMAD_CLASSES
#endif
''')

f_class.close()

##################################################################################
##################################################################################
# Create C++ side of interface

f_cpp = open('code/cpp_bmad_convert.cpp', 'w')
f_cpp.write('''
//+
// C++ side of the Bmad / C++ structure interface.
//
// File Generated by: create_interface.py
// Do not edit this file directly! 
//-

#ifndef CPP_BMAD_CONVERT

#include <iostream>
#include "cpp_bmad_classes.h"

//---------------------------------------------------------------------------

template <class T> void operator<< (valarray<T>& arr, const T* ptr) {
  int n = arr.size();
  for (int i = 0; i < n; i++) arr[i] = ptr[i];
}

template <class T> void operator<< (valarray< valarray<T> >& mat, const T* ptr) {
  int n1 = mat.size();
  if (n1 == 0) return;
  int n2 = mat[0].size();
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      mat[i][j] = ptr[i*n2+j];
    }
  }
}

template <class T> void operator<< (valarray< valarray< valarray<T> > >& tensor, const T* ptr) {
  int n1 = tensor.size();
  if (n1 == 0) return;
  int n2 = tensor[0].size();
  int n3 = tensor[0][0].size();
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      for (int k = 0; k < n3; k++) {
        tensor[i][j][k] = ptr[i*n2*n3 + j*n3 + k];
      }
    }
  }
}

template <class T> void operator<< (valarray<T>& arr1, const valarray<T>& arr2) {
  int n1 = arr1.size(), n2 = arr2.size();
  if (n1 != n2) arr1.resize(n2);
  arr1 = arr2;
}

template <class T> void operator<< (valarray< valarray<T> >& mat1, 
                              const valarray< valarray<T> >& mat2) {
  int n1_1 = mat1.size(), n2_1 = mat2.size();
  int n1_2 = 0, n2_2 = 0;
  if (n1_1 > 0) n1_2 = mat1[0].size();
  if (n2_1 > 0) n2_2 = mat2[0].size();
  if (n1_1 != n2_1) mat1.resize(n2_1);
  if (n1_2 != n2_2) {for (int i = 0; i < n1_1; i++) mat1[i].resize(n2_2);}
  mat1 = mat2;
}

//---------------------------------------------------------------------------
// Instantiate needed template instances.

template void operator<< (Real_Array&,  const double*);
template void operator<< (Real_Matrix&, const double*);
template void operator<< (Real_Tensor&, const double*);
template void operator<< (Int_Array&,   const int*);

template void operator<< (Real_Array&,  const Real_Array&);
template void operator<< (Real_Matrix&, const Real_Matrix&);
template void operator<< (Int_Array&,   const Int_Array&);

//---------------------------------------------------------------------------

void matrix_to_vector (const Real_Matrix& mat, double* vec) {
  int n1 = mat.size();
  if (n1 == 0) return;
  int n2 = mat[0].size();
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      vec[i*n2+j] = mat[i][j];
    }
  }
}

//---------------------------------------------------------------------------

void tensor_to_vector (const Real_Tensor& tensor, double* vec) {
  int n1 = tensor.size();
  if (n1 == 0) return;
  int n2 = tensor[0].size();
  int n3 = tensor[0][0].size();
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      for (int k = 0; k < n3; k++) {
        vec[i*n2*n3 + j*n3 + k] = tensor[i][j][k];
      }
    }
  }
}

''')

for key, struct in struct_def.items():

  # zzz_to_f2
  f_cpp.write('''
//--------------------------------------------------------------------
//--------------------------------------------------------------------
// zzz

extern "C" void zzz_to_c (zzz_struct*, C_zzz&);

extern "C" void zzz_to_f2 (zzz_struct*'''.replace('zzz', struct.short_name))

  for var in struct.var:
    f_cpp.write (', ' + var.c_side.to_f2_arg.replace('zzz', struct.short_name))

  f_cpp.write(');\n')

  # zzz_to_f

  f_cpp.write('''
extern "C" void zzz_to_f (C_zzz& C, zzz_struct* F) {
  zzz_to_f2 (F'''.replace('zzz', struct.short_name))

  for var in struct.var:
    f_cpp.write (', ' + var.c_side.to_f2_call.replace('name', var.name))

  f_cpp.write(');\n')
  f_cpp.write('}\n')

  # zzz_to_c2

  f_cpp.write('\n')
  f_cpp.write('extern "C" void zzz_to_c2 (C_zzz& C'.replace('zzz', struct.short_name))

  for var in struct.var:
    f_cpp.write (', ' + var.c_side.to_c2_arg.replace('name', var.name).replace('zzz', struct.short_name))

  f_cpp.write(') {\n')

  for var in struct.var:
    f_cpp.write ('  ' + var.c_side.to_c2_set.replace('name', var.name) + '\n')

  f_cpp.write('}\n')

f_cpp.write('''
#define CPP_BMAD_CONVERT
#endif
''')

f_cpp.close()

##################################################################################
##################################################################################
# Create C++ class equality check code

f_eq = open('code/cpp_equality.cpp', 'w')

f_eq.write('''
//+
// C++ equality functions for Bmad / C++ structure interface.
//
// File Generated by: create_interface.py
// Do not edit this file directly! 
//-

#include <iostream>
#include <stdlib.h>
#include "cpp_bmad_classes.h"

using namespace std;

//---------------------------------------------------

template <class T> bool is_all_equal (const valarray<T>& v1, const valarray<T>& v2) {
  bool is_eq = true;
  T t1, t2;
  if (v1.size() != v2.size()) return false;
  for (int i = 0; i < v1.size(); i++) {
    t1 = v1[i];
    t2 = v2[i];
    is_eq = is_eq && (v1[i] == v2[i]);
  }
  return is_eq;
}

//---------------------------------------------------

bool is_all_equal (const Real_Matrix& mat1, const Real_Matrix& mat2) {
  bool is_eq = true;
  if (mat1.size() != mat2.size()) return false;
  for (int i = 0; i < mat1.size(); i++) {
    if (mat1[i].size() != mat2[i].size()) return false;
    for (int j = 0; j < mat1[i].size(); j++) {
      is_eq = is_eq && (mat1[i][j] == mat2[i][j]);
    }
  }
  return is_eq;
};

bool is_all_equal (const Real_Tensor& tensor1, const Real_Tensor& tensor2) {
  bool is_eq = true;
  if (tensor1.size() != tensor2.size()) return false;
  for (int i = 0; i < tensor1.size(); i++) {
    if (tensor1[i].size() != tensor2[i].size()) return false;
    for (int j = 0; j < tensor1[i].size(); j++) {
      if (tensor1[i][j].size() != tensor2[i][j].size()) return false;
      for (int k = 0; k < tensor1[i][j].size(); k++) {
        is_eq = is_eq && (tensor1[i][j][k] == tensor2[i][j][k]);
      }
    }
  }
  return is_eq;
};

//---------------------------------------------------

template bool is_all_equal (const Dcomplex_Array&, const Dcomplex_Array&);
template bool is_all_equal (const Real_Array&,     const Real_Array&);
template bool is_all_equal (const Int_Array&,      const Int_Array&);
''')

for key, struct in struct_def.items():
  f_eq.write ('\n//--------------------------------------------------------------\n\n')
  f_eq.write ('bool operator== (const C_zzz& x, const C_zzz& y) {'.replace('zzz', struct.short_name) + '\n')
  f_eq.write ('  bool is_eq = true;\n')

  for var in struct.var:
    f_eq.write ('  is_eq = is_eq && ' + var.c_side.equal_test.replace('name', var.name)  + ';\n')

  f_eq.write ('  return is_eq;\n')
  f_eq.write ('};\n\n')

  f_eq.write ('template bool is_all_equal (const C_zzz_array&, const C_zzz_array&);\n'.replace('zzz', struct.short_name))

f_eq.close()

##################################################################################
##################################################################################
# Create C++ side code check

f_test = open('interface_test/cpp_bmad_test.cpp', 'w')
f_test.write('''
//+
// C++ classes definitions for Bmad / C++ structure interface.
//
// File Generated by: create_interface.py
// Do not edit this file directly! 
//-

#include <stdio.h>
#include <iostream>
#include "cpp_bmad_classes.h"

using namespace std;
''')

for key, struct in struct_def.items():
  f_test.write ('''
//--------------------------------------------------------------
//--------------------------------------------------------------

extern "C" void test2_f_zzz (C_zzz&, bool&);

void C_zzz_test_pattern (C_zzz& C, int ix_patt) {

int offset = 100 * ix_patt;

'''.replace('zzz', struct.short_name))

for i, var in enumerate(struct.var, 1):
  f_test.write (var.c_side.test_pat.replace('XXX', str(i)).replace('name', var.name) + '\n')

f_test.write('''
}

//--------------------------------------------------------------

extern "C" void test_c_zzz (zzz_struct* F, bool& c_ok) {

  C_zzz C, C2;

  c_ok = true;

  zzz_to_c (F, C);
  C_zzz_test_pattern (C2, 1);

  if (C == C2) {
    cout << " zzz: C side convert F->C: Good" << endl;
  } else {
    cout << " zzz: C SIDE CONVERT C->F: FAILED!" << endl;
    c_ok = false;
  }

  C_zzz_test_pattern (C2, 2);
  bool c_ok2;
  test2_f_zzz (C2, c_ok2);
  if (!c_ok2) c_ok = false;

  C_zzz_test_pattern (C, 3);
  if (C == C2) {
    cout << " zzz: F side convert F->C: Good" << endl;
  } else {
    cout << " zzz: F SIDE CONVERT C->F: FAILED!" << endl;
    c_ok = false;
  }

  C_zzz_test_pattern (C2, 4);
  zzz_to_f (C2, F);

}
'''.replace('zzz', struct.short_name))

f_test.close()
