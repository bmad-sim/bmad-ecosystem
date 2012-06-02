#!/usr/bin/env python

# Note: Run this script in the cpp_bmad_interface directory.

# Script to read in Fortran structures and create:
#   Corresponding C++ class
#   Translator between Fortran structure and C++ class
#   Routines to check for equality between instances of a given fortran structure.
#   Routines to check for equality between instances of a given C++ class
#   Program to check the Fortran / C++ translator

# Note: The corresponding C++ class component for a pointer or allocatable Fortran 
# scaler struct component is an array whose length is zero if the Fortran component
# is nullified and whose length is 1 otherwise.

import sys
import shutil
import os
import copy
import re

##################################################################################
##################################################################################
# For printing of intermediate steps, etc.

debug = False # Change to True to enable printout
 
def print_debug (line):
  if (debug): print line

def indent (string, numspace):
  x = ' ' * numspace
  if string[-1] == '\n':
    return x + string[:-1].replace('\n', '\n' + x) + '\n'
  else:
    return x + string.replace('\n', '\n' + x)

##################################################################################
##################################################################################
# struct_def_class
# Class for a structure
# Note: The .arg array will contain, in addtion to the structure components, the
# array bounds for allocatable and pointer structure components.

class struct_def_class:
  def __init__(self, f_name = ''):
    self.f_name = f_name # Struct name on Fortran side
    self.short_name = '' # Struct name without trailing '_struct'. Note: C name is 'C_<name>'
    self.arg = []        # Array of arg_class. List of structrure components + array bound dimensions. 

  def __repr__(self):
    return '[name: %s, #arg: %i]' % (self.short_name, len(self.arg))

# arg_class.

class arg_class:

  def __init__(self):
    self.is_component = True   # Is a structure component? If not, then will be array bound.
    self.name = ''             # Name of argument
    self.type = ''             # Fortran type without '(...)'. EG: 'real', 'type', 'character', etc.
    self.kind = ''             # Fortran kind. EG: '', 'rp', 'coord_struct', etc.
    self.pointer_type = 'NOT'  # '-', 'PTR', 'ALLOC'
    self.array = []            # EG: [':', ':'] or ['0:6', '3']
    self.full_array = ''       # EG: '(:,:)', '(0:6, 3)'
    self.lbound = []
    self.ubound = []
    self.init_value = ''       # Initialization value
    self.comment = ''          # Comment with Fortran structure def.
    self.f_side = 0
    self.c_side = 0

  def __repr__(self):
    return '["%s(%s)", "%s", "%s", %s, "%s"]' % (self.type, self.kind, self.pointer_type, self.name, self.array, self.init_value)

  def full_repr(self):
    return '["%s(%s)", "%s", "%s", %s, "%s" %s %s "%s"]' % (self.type, 
              self.kind, self.pointer_type, self.name, self.array, self.full_array, 
              self.lbound, self.ubound, self.init_value)

##################################################################################
##################################################################################
# Translations

NOT = 'NOT'
PTR = 'PTR'
ALLOC = 'ALLOC'

T = True
F = False

REAL  = 'real'
CMPLX = 'complex'
INT   = 'integer'
LOGIC = 'logical'
CHAR  = 'character'
TYPE  = 'type'
SIZE  = 'size'

# Fortran side translation

class f_side_trans_class:

  def __init__(self):
    self.to_c2_call = ''
    self.bindc_type = ''
    self.bindc_name = ''
    self.equality_test = 'is_eq = is_eq .and. all(f1%NAME == f2%NAME)\n'
    self.to_c2_f2_sub_arg = 'z_NAME'
    self.to_f2_trans = 'F%NAME = z_NAME'
    self.to_f2_extra_var_name = ''
    self.to_f2_extra_var_type = ''
    self.test_pat = 'rhs = XXX + offset; F%NAME = NNN\n'
    self.to_c_var = ''
    self.to_c_trans = ''
    self.size_var = []              # For communicating the size of allocatable and pointer variables

  def __repr__(self):
    return '%s,  %s,  %s :: %s' % (self.to_c2_call, self.bindc_type, self.bindc_name, self.to_f2_trans)

#------------------------

x2 = ' ' * 2
x4 = ' ' * 4
x6 = ' ' * 6

to_f2_trans_pointer = \
'''
if (associated(F%NAME)) then
  if (n1_NAME == 0 .or. any(shape(F%NAME) /= [DIMS])) deallocate(F%NAME)
endif
if (n1_NAME /= 0) then
  call c_f_pointer (z_NAME, f_NAME, [TOTDIM])
  if (.not. associated(F%NAME)) allocate(F%NAME(DIMS))
  SET
else
  if (associated(F%NAME)) deallocate(F%NAME)
endif'''

equality_test_pointer = \
'''
is_eq = is_eq .and. (associated(f1%NAME) .eqv. associated(f2%NAME))
if (.not. is_eq) return
if (associated(f1%NAME)) is_eq = all(shape(f1%NAME) == shape(f2%NAME))
if (.not. is_eq) return
if (associated(f1%NAME)) is_eq = all(f1%NAME == f2%NAME)'''

jd1_loop = 'do jd1 = 1, size(F%NAME,1); lb1 = lbound(F%NAME,1) - 1\n'
jd2_loop = 'do jd2 = 1, size(F%NAME,2); lb2 = lbound(F%NAME,2) - 1\n'
jd3_loop = 'do jd3 = 1, size(F%NAME,3); lb3 = lbound(F%NAME,3) - 1\n'

rhs1 = '  rhs = 100 + jd1 + XXX + offset\n'
rhs2 = '  rhs = 100 + jd1 + 10*jd2 + XXX + offset\n'
rhs3 = '  rhs = 100 + jd1 + 10*jd2 + 100*jd3 + XXX + offset\n'

set1 = '  F%NAME(jd1+lb1) = NNN\n' 
set2 = '  F%NAME(jd1+lb1,jd2+lb2) = NNN\n' 
set3 = '  F%NAME(jd1+lb1,jd2+lb2,jd3+lb3) = NNN\n' 

test_pat1 = jd1_loop + rhs1 + set1 + 'enddo\n'
test_pat2 = jd1_loop + jd2_loop + rhs2 + set2 + 'enddo; enddo\n'
test_pat3 = jd1_loop + jd2_loop + jd3_loop + rhs3 + set3 + 'enddo; enddo; enddo\n'

n_str = {0:'', 1:'1', 2:'2', 3:'3'}

f_side_trans = {}

for type in [REAL, CMPLX, INT, LOGIC, TYPE, SIZE]:
  for dim in range(4):
    f_side_trans[type, dim, NOT] = copy.deepcopy(f_side_trans_class())
    f = f_side_trans[type, dim, NOT]

    if type == REAL:
      f.bindc_type = 'real(c_double)'
      test_value = 'rhs'

    if type == CMPLX:
      f.bindc_type = 'complex(c_double_complex)'
      test_value = 'cmplx(rhs, 100+rhs)'

    if type == INT:
      f.bindc_type = 'integer(c_int)'
      test_value = 'rhs'

    if type == LOGIC:
      f.bindc_type = 'logical(c_bool)'
      test_value = '(modulo(rhs, 2) == 0)'

    if type == TYPE:
      f.bindc_type = 'type(c_ptr)'
      test_value = 'NNN'

    if type == SIZE:
      f.to_c2_call       = 'NAME'
      f.bindc_type       = 'integer(c_int), value'
      f.bindc_name       = 'NAME'
      f.to_f2_trans      = ''
      f.to_c2_f2_sub_arg = 'NAME'
      f.to_c_var         = 'integer(c_int) :: NAME'
      f.test_value       = ''
      continue

    #-----------------------------------------------

    if dim == 0: 
      f.bindc_name    = 'z_NAME'
      f.equality_test = 'is_eq = is_eq .and. (f1%NAME == f2%NAME)\n'
      f.to_c2_call    = 'F%NAME'
      f.test_pat      = f.test_pat
      if type == LOGIC:
        f.to_c2_call  = 'c_logic(F%NAME)'
        f.to_f2_trans = 'F%NAME = f_logic(z_NAME)'
        f.equality_test = f.equality_test.replace('==', '.eqv.')
      if type == TYPE:
        f.bindc_type = 'type(c_ptr), value'
        f.to_c2_call = 'c_loc(F%NAME)'
        f.to_f2_trans = 'call KIND_to_f(z_NAME, c_loc(F%NAME))'
        f.test_pat    = 'call set_KIND_test_pattern (F%NAME, ix_patt)\n'

    #-----------------------------------------------

    if dim == 1:
      f.to_c2_call  = 'fvec2vec(F%NAME, DIM1)'
      f.bindc_name  = 'z_NAME(*)'
      f.to_f2_trans = 'F%NAME = z_NAME(1:DIM1)'
      f.test_pat    = test_pat1
      if type == LOGIC:
        f.to_f2_trans = 'call vec2fvec (z_NAME, F%NAME)'
        f.equality_test = f.equality_test.replace('==', '.eqv.')
      if type == TYPE:
        f.to_c2_call = 'z_NAME'
        f.to_f2_trans = jd1_loop + '  call KIND_to_f(z_NAME(jd1), c_loc(F%NAME(jd1+lb1)))\nenddo'
        f.test_pat    = jd1_loop + rhs1 + '  call set_KIND_test_pattern (F%NAME(jd1+lb1), ix_patt+jd1)\n' + 'enddo\n'
        f.to_c_var    = 'type(c_ptr) :: z_NAME(DIM1)'
        f.to_c_trans  = jd1_loop + '  z_NAME(jd1) = c_loc(F%NAME(jd1+lb1))\nenddo\n\n'

    #-----------------------------------------------

    if dim == 2:
      f.to_f2_trans = 'call vec2mat(z_NAME, F%NAME)'
      f.to_c2_call  = 'mat2vec(F%NAME, SIZE2)'
      f.bindc_name  = 'z_NAME(*)'
      f.test_pat    = test_pat2
      if type == LOGIC:
        f.equality_test = f.equality_test.replace('==', '.eqv.')
      if type == TYPE:
        f.to_c2_call = 'z_NAME'
        f.to_f2_trans = jd1_loop + jd2_loop + \
                    '  call KIND_to_f(z_NAME(DIM2*(jd1-1) + jd2), c_loc(F%NAME(jd1+lb1,jd2+lb2)))\n' + 'enddo; enddo\n'
        f.test_pat    = jd1_loop + jd2_loop + rhs2 + \
                    '  call set_KIND_test_pattern (F%NAME(jd1+lb1,jd2+lb2), ix_patt+jd1+10*jd2)\n' + 'enddo; enddo\n'
        f.to_c_var    = 'type(c_ptr) :: z_NAME(DIM1*DIM2)'
        f.to_c_trans  = jd1_loop + jd2_loop + '  z_NAME(DIM2*(jd1-1) + jd2) = c_loc(F%NAME(jd1+lb1,jd2+lb2))\n' + \
                                              'enddo; enddo\n\n'

    #-----------------------------------------------

    if dim == 3:
      f.to_f2_trans = 'call vec2tensor(z_NAME, F%NAME)'
      f.to_c2_call  = 'tensor2vec(F%NAME, SIZE3)'
      f.bindc_name  = 'z_NAME(*)'
      f.test_pat    = test_pat3
      if type == LOGIC:
        f.equality_test = f.equality_test.replace('==', '.eqv.')
      if type == TYPE:
        f.to_c2_call = 'z_NAME'
        f.to_f2_trans = jd1_loop + jd2_loop + jd3_loop + \
              '  call KIND_to_f(z_NAME(DIM3*DIM2*(jd1-1) + DIM3*(jd2-1) + jd3), c_loc(F%NAME(jd1+lb1,jd2+lb2,jd3+lb3)))\n' + \
              'enddo; enddo; enddo\n'
        f.test_pat    = jd1_loop + jd2_loop + jd3_loop + rhs3 + \
              '  call set_KIND_test_pattern (F%NAME(jd1+lb1,jd2+lb2,jd3+lb3), ix_patt+jd1+10*jd2+100*jd3)\n' + \
              'enddo; enddo; enddo\n'
        f.to_c_var    = 'type(c_ptr) :: z_NAME(DIM1*DIM2*DIM3)'
        f.to_c_trans  = jd1_loop + jd2_loop + jd3_loop + \
              '  z_NAME(DIM3*DIM2*(jd1-1) + DIM3*(jd2-1) + jd3) = c_loc(F%NAME(jd1+lb1,jd2+lb2,jd3+lb3))\n' + \
              'enddo; enddo; enddo\n\n'

    #-------------------------

    f.test_pat    = f.test_pat.replace('NNN', test_value)

    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------
    # Pointers

    f_side_trans[type, dim, PTR] = f_side_trans_class()
    fp = f_side_trans[type, dim, PTR]
    fp.to_f2_extra_var_name = 'f_NAME(:)'
    fp.to_f2_extra_var_type = f.bindc_type + ', pointer'
    fp.bindc_type = 'type(c_ptr), value'
    fp.bindc_name = 'z_NAME'

    #---------------------
    # dim = 0

    if dim == 0:
      fp.to_c2_call = 'c_loc(F%NAME)'
      fp.to_f2_extra_var_name = 'f_NAME'
      fp.to_f2_trans = '''
if (n_NAME == 0) then
  if (associated(F%NAME)) deallocate(F%NAME)
else
  call c_f_pointer (z_NAME, f_NAME)
  if (.not. associated(F%NAME)) allocate(F%NAME)
  F%NAME = f_NAME
endif'''

      fp.to_c_trans = '''\
n_NAME = 0
if (associated(F%NAME)) n_NAME = 1
'''

      fp.equality_test = '''
is_eq = is_eq .and. (associated(f1%NAME) .eqv. associated(f2%NAME))
if (.not. is_eq) return
if (associated(f1%NAME)) is_eq = (f1%NAME == f2%NAME)'''

      fp.test_pat    = '''
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME)
  rhs = XXX + offset
  F%NAME = NNN
endif
'''

      if type == LOGIC:
        fp.equality_test = fp.equality_test.replace('== f', '.eqv. f')
        fp.to_f2_trans = fp.to_f2_trans.replace('= f_NAME', '= f_logic(f_NAME)')
        fp.to_c2_call = 'c_loc(fscaler2scaler(F%NAME, n_NAME))'

      if type == TYPE:
        fp.to_f2_extra_var_type = 'type(KIND_struct), pointer'
        fp.test_pat = fp.test_pat.replace('F%NAME = NNN', 'call set_KIND_test_pattern (F%NAME, ix_patt)')

    #---------------------
    # Pointer, dim = 1

    if dim == 1:

      fp.to_c2_call = 'c_loc(fvec2vec(F%NAME, n1_NAME))'
      fp.to_f2_trans = to_f2_trans_pointer.replace('DIMS', \
                    'n1_NAME').replace('TOTDIM', 'n1_NAME').replace('SET', 'F%NAME = f_NAME(1:n1_NAME)')
      fp.equality_test = equality_test_pointer

      fp.to_c_trans  = \
'''n1_NAME = 0
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME, 1)
endif
'''
      fp.test_pat    = '''
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME(-1:1))
''' + x2 + jd1_loop + x2 + rhs1 + x2 + set1 + '  enddo\n' + 'endif\n'

      if type == LOGIC:
        fp.equality_test = fp.equality_test.replace('== f', '.eqv. f')
        fp.to_f2_trans = fp.to_f2_trans.replace('F%NAME = f_NAME(1:n1_NAME)', 'call vec2fvec (f_NAME, F%NAME)')

    #---------------------
    # Pointer, dim = 2

    if dim == 2:
      fp.to_c2_call = 'c_loc(mat2vec(F%NAME, n1_NAME*n2_NAME))'
      fp.to_f2_trans = to_f2_trans_pointer.replace('DIMS', 'n1_NAME, n2_NAME')\
                                          .replace('TOTDIM', 'n1_NAME*n2_NAME')\
                                          .replace('SET', 'call vec2mat(f_NAME, F%NAME)')
      fp.equality_test = equality_test_pointer

      fp.to_c_trans  = \
'''n1_NAME = 0; n2_NAME = 0
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME, 1)
  n2_NAME = size(F%NAME, 2)
endif
'''
      fp.test_pat    = '''
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME(-1:1, 2))
''' + x2 + jd1_loop + x2 + jd2_loop + x2 + rhs2 + \
      x2 + set2 + '  enddo; enddo\n' + 'endif\n'

      if type == LOGIC:
        fp.equality_test = fp.equality_test.replace('== f', '.eqv. f')

    #---------------------
    # Pointer, dim = 3

    if dim == 3:
      fp.to_c2_call = 'c_loc(tensor2vec(F%NAME, n1_NAME*n2_NAME*n3_NAME))'
      fp.to_f2_trans = to_f2_trans_pointer.replace('DIMS', 'n1_NAME, n2_NAME, n3_NAME')\
                                          .replace('TOTDIM', 'n1_NAME*n2_NAME*n3_NAME')\
                                          .replace('SET', 'call vec2tensor(f_NAME, F%NAME)')
      fp.equality_test = equality_test_pointer

      fp.to_c_trans  = \
'''n1_NAME = 0; n2_NAME = 0; n3_NAME = 0
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME, 1)
  n2_NAME = size(F%NAME, 2)
  n3_NAME = size(F%NAME, 3)
endif
'''
      fp.test_pat    = '''
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME(-1:1, 2, 1))
''' + x2 + jd1_loop + x2 + jd2_loop + x2 + jd3_loop + x2 + rhs3 + \
      x2 + set3 + '  enddo; enddo; enddo\n' + 'endif\n'

      if type == LOGIC:
        fp.equality_test = fp.equality_test.replace('== f', '.eqv. f')

    #----------

    fp.test_pat = fp.test_pat.replace('NNN', test_value)

#---------------------------
# CHAR

f_side_trans[CHAR, 0, NOT] = f_side_trans_class()
f_side_trans[CHAR, 0, NOT].bindc_type = 'character(c_char)'
f_side_trans[CHAR, 0, NOT].bindc_name = 'z_NAME(*)'
f_side_trans[CHAR, 0, NOT].to_c2_call = 'trim(F%NAME) // c_null_char'
f_side_trans[CHAR, 0, NOT].equality_test = 'is_eq = is_eq .and. (f1%NAME == f2%NAME)\n'
f_side_trans[CHAR, 0, NOT].test_pat    = \
        'do jd1 = 1, len(F%NAME)\n  F%NAME(jd1:jd1) = char(ichar("a") + modulo(100+XXX+offset+jd1, 26))\nenddo\n'
f_side_trans[CHAR, 0, NOT].to_f2_trans = 'call to_f_str(z_NAME, F%NAME)'

# Allocatable components are very similar to pointer components
# with the simple replacement of 'allocated' for 'associated'.

for trans in f_side_trans.keys():
  if trans[2] == PTR:
    trans_alloc = (trans[0], trans[1], ALLOC)
    f_side_trans[trans_alloc] = copy.deepcopy(f_side_trans[trans])
    t = f_side_trans[trans_alloc]
    t.to_f2_trans   = t.to_f2_trans.replace('associated', 'allocated')
    t.to_c_trans    = t.to_c_trans.replace('associated', 'allocated')
    t.test_pat      = t.test_pat.replace('associated', 'allocated')
    t.equality_test = t.equality_test.replace('associated', 'allocated')

#############################################################

class c_side_trans_class:

  def __init__(self):
    self.c_class = ''
    self.to_f_setup = ''
    self.to_f_cleanup = ''
    self.to_f2_arg = ''
    self.to_f2_call = ''
    self.to_c2_arg = ''
    self.to_c2_set = '  C.NAME = z_NAME;'
    self.constructor = 'NAME(0)'
    self.destructor = ''
    self.equality_test = '  is_eq = is_eq && (x.NAME == y.NAME);\n'
    self.test_pat = '  rhs = XXX + offset; C.NAME = NNN;\n'

  def __repr__(self):
    return '%s,  %s,  %s,  %s' % (self.c_class, self.to_f2_arg, self.to_f2_call, self.to_c2_arg)
    self.size_var = []       # For communicating the size of allocatable and pointer variables

#------------------

to_c2_set_pointer = \
'''
  if (n1_NAME == 0) 
    delete C.NAME;
  else {
    C.NAME = new KIND;
    SET
  }
'''

test_pat_pointer = \
'''  if (ix_patt < 3) 
    C.NAME == NULL;
  else {
'''

equality_test_pointer = \
'''
  is_eq = is_eq && ((x.NAME == NULL) == (y.NAME == NULL));
  if (!is_eq) return false;
  if (x.NAME != NULL) is_eq = TEST;
'''

for1 = '  for (int i = 0; i < C.NAME.size(); i++)'
for2 = ' for (int j = 0; j < C.NAME[0].size(); j++) ' 
for3 = ' for (int k = 0; k < C.NAME[0][0].size(); k++)'

test_pat1 = for1 + '\n    {int rhs = 101 + i + XXX + offset; C.NAME[i] = NNN;}'
test_pat2 = for1 + for2 + '\n    {int rhs = 101 + i + 10*(j+1) + XXX + offset; C.NAME[i][j] = NNN;}'
test_pat3 = for1 + for2 + for3 + '\n    {int rhs = 101 + i + 10*(j+1) + 100*(k+1) + XXX + offset; C.NAME[i][j][k] = NNN;}'

c_side_trans = {}

for type in [REAL, CMPLX, INT, LOGIC, TYPE, SIZE]:
  for dim in range(4):
    c_side_trans[type, dim, NOT] = c_side_trans_class()
    c = c_side_trans[type, dim, NOT]

    if type == REAL:
      c_type = 'Real'
      c_arg  = 'c_Real'
      test_value = 'rhs'
      construct_value = '0.0'

    if type == CMPLX:
      c_type = 'Complex'
      c_arg  = 'c_Complex'
      test_value = 'Complex(rhs, 100+rhs)'
      construct_value = '0.0'

    if type == INT:
      c_type = 'Int'
      c_arg  = 'c_Int'
      test_value = 'rhs'
      construct_value = '0'

    if type == LOGIC:
      c_type = 'Bool'
      c_arg  = 'c_Bool'
      test_value = '(rhs % 2 == 0)'
      construct_value = 'false'

    if type == TYPE:
      c_type = 'C_KIND'
      c_arg  = 'const C_KIND'
      test_value = ''
      construct_value = ''

    if type == SIZE:
      c.to_f2_arg  = 'Int'
      c.to_f2_call = 'NAME'
      c.to_c2_arg  = 'Int NAME'
      continue

    #-------------------------------------------------------

    if dim == 0: 
      c.c_class      = c_type
      c.to_f2_arg    = c_arg + '&'
      c.to_f2_call   = 'C.NAME'
      c.to_c2_arg    = c_arg + '& z_NAME'
      c.test_pat     = c.test_pat

      if type == TYPE:
        c.constructor = 'NAME()'
        c.to_c2_set   = '  KIND_to_c(z_NAME, C.NAME);' 
        c.test_pat    = '  set_C_KIND_test_pattern(C.NAME, ix_patt);\n'

    #-------------------------------------------------------

    if dim == 1:
      c.c_class      = c_type + '_Array'
      c.to_f2_arg    = c_arg + 'Arr'
      c.to_f2_call   = '&C.NAME[0]'
      c.to_c2_arg    = c_arg + 'Arr z_NAME'
      c.constructor  = 'NAME(' + construct_value + ', DIM1)'
      c.to_c2_set    = '  C.NAME = ' + c_type + '_Array(z_NAME, DIM1);'
      c.test_pat    = test_pat1
      c.equality_test = '  is_eq = is_eq && is_all_equal(x.NAME, y.NAME);\n'

      if type == TYPE:
        c.constructor = 'NAME(C_KIND_Array(C_KIND(), DIM1))'
        c.to_c2_set   = for1 + ' KIND_to_c(z_NAME[i], C.NAME[i]);' 
        c.test_pat    = test_pat1.replace('C.NAME[i] = NNN', 'set_C_KIND_test_pattern(C.NAME[i], ix_patt+i+1);')
        c.to_f_setup  = '  const C_KIND* z_NAME[DIM1];\n' + for1 + ' z_NAME[i] = &C.NAME[i];\n'

    #-------------------------------------------------------

    if dim == 2:
      c.c_class      = c_type + '_Matrix'
      c.to_f2_arg    = c_arg + 'Arr'
      c.to_f2_call   = 'z_NAME'
      c.to_c2_arg    = c_arg + 'Arr z_NAME'
      c.constructor  = 'NAME(' + c_type + '_Array(' + construct_value + ', DIM2), DIM1)'
      c.to_c2_set    = '  C.NAME << z_NAME;'
      c.test_pat     = test_pat2
      c.to_f_setup  = '  ' + c_type + ' z_NAME[DIM1*DIM2]; matrix_to_vec(C.NAME, z_NAME);\n'
      c.equality_test = '  is_eq = is_eq && is_all_equal(x.NAME, y.NAME);\n'

      if type == TYPE:
        c.constructor = 'NAME(C_KIND_Array(C_KIND(), DIM2), DIM1)'
        c.to_c2_set   = for1 + for2 + '\n    {int m = DIM2*i + j; KIND_to_c(z_NAME[m], C.NAME[i][j]);}' 
        c.test_pat    = test_pat2.replace('C.NAME[i][j] = NNN', 'set_C_KIND_test_pattern(C.NAME[i][j], ix_patt+i+1+10*(j+1))')
        c.to_f_setup  = '  const C_KIND* z_NAME[DIM1*DIM2];\n' + \
                  for1 + for2 +  '\n    {int m = DIM2*i + j; z_NAME[m] = &C.NAME[i][j];}\n'

    #-------------------------------------------------------

    if dim == 3:
      c.c_class      = c_type + '_Tensor'
      c.to_f2_arg    = c_arg + 'Arr'
      c.to_f2_call   = 'z_NAME'
      c.to_c2_arg    = c_arg + 'Arr z_NAME'
      c.constructor  = 'NAME(' + c_type + '_Matrix(' + c_type + '_Array(' + construct_value + ', DIM3), DIM2), DIM1)'
      c.to_c2_set    = '  C.NAME << z_NAME;'
      c.test_pat     = test_pat3
      c.to_f_setup   = '  ' + c_type + ' z_NAME[DIM1*DIM2*DIM3]; tensor_to_vec(C.NAME, z_NAME);\n'
      c.equality_test = '  is_eq = is_eq && is_all_equal(x.NAME, y.NAME);\n'

      if type == TYPE:
        c.constructor = 'NAME(C_KIND_Matrix(C_KIND_Array(C_KIND(), DIM3), DIM2), DIM1)'
        c.to_c2_set   = for1 + for2 + for3 + '\n    {int m = DIM3*DIM2*i + DIM3*j + k; KIND_to_c(z_NAME[m], C.NAME[i][j][k]);}' 
        c.test_pat    = test_pat3.replace('C.NAME[i][j][k] = NNN', \
                                          'set_C_KIND_test_pattern(C.NAME[i][j][k], ix_patt+i+1+10*(j+1)+100*(k+1))')
        c.to_f_setup  = '  const C_KIND* z_NAME[DIM1*DIM2*DIM3];\n' + \
                  for1 + for2 + for3 + '\n    {int m = DIM3*DIM2*i + DIM3*j + k; z_NAME[m] = &C.NAME[i][j][k];}\n'

    #-------------------------------------------------------

    c.test_pat    = c.test_pat.replace('NNN', test_value)

    if type == TYPE:
      c.to_c2_arg = 'const KIND_struct* z_NAME'
      if dim > 0: 
        c.to_f2_arg  = c.to_f2_arg.replace('Arr', '**')
        c.to_c2_arg  = 'const KIND_struct** z_NAME'
        c.to_f2_call = 'z_NAME'

    #-------------------------------------------------------
    # Pointers

    c_side_trans[type, dim, PTR] = copy.deepcopy(c)
    cp = c_side_trans[type, dim, PTR] 
    cp.c_class       = c.c_class + '*'
    cp.constructor   = 'NAME(NULL)'
    cp.destructor    = 'delete NAME;'
    cp.equality_test = equality_test_pointer.replace('TEST', 'is_all_equal(*x.NAME, *y.NAME)')

    if type != TYPE:
      cp.to_f2_arg     = c_arg + 'Arr'
      cp.to_c2_arg     = cp.to_f2_arg + ' z_NAME'
    else:
      cp.to_c2_arg     = 'KIND_struct* z_NAME'

    #------------------------------
    # Pointer, dim = 0

    if dim == 0:
      cp.to_c2_set     = to_c2_set_pointer.replace('KIND', c_type).replace('n1_NAME', 'n_NAME')
      cp.equality_test = equality_test_pointer.replace('TEST', '(*x.NAME == *y.NAME)')
      cp.test_pat      = test_pat_pointer + '    C.NAME = new ' + c_type + ';\n' + \
                                    indent(c.test_pat.replace('C.NAME', '(*C.NAME)'), 2) + '  }\n'
      cp.to_f_setup    = '  int n_NAME = 0; if (C.NAME != NULL) n_NAME = 1;\n'
 
      if type == TYPE:
        cp.test_pat   = cp.test_pat
        cp.to_f2_call = '*C.NAME' 
        cp.to_c2_set = cp.to_c2_set.replace('SET', 'KIND_to_c(z_NAME, *C.NAME);')
      else:
        cp.to_c2_set = cp.to_c2_set.replace('SET', '*C.NAME = *z_NAME;')

    #------------------------------
    # Pointer, dim = 1

    if dim == 1:
      cp.to_f2_call  = 'z_NAME' 
      cp.to_c2_set   = to_c2_set_pointer.replace('KIND', c_type + '_Array').replace('SET', \
                                  '*C.NAME = ' + c_type + '_Array(z_NAME, n1_NAME);')
      cp.test_pat    = test_pat_pointer + '    C.NAME = new ' + c_type + '_Array(3);\n' + \
                                    indent(c.test_pat.replace('C.NAME', '(*C.NAME)'), 2) + '  }\n'
      cp.to_f_setup  = '''
  int n1_NAME = 0;
  c_TYPEArr z_NAME = NULL;
  if (C.NAME != NULL) {
    n1_NAME = (*C.NAME).size();
    z_NAME = &(*C.NAME)[0];
  }
'''.replace('TYPE', c_type)

      if type == TYPE:
        cp.test_pat   = cp.test_pat.replace('C.NAME', 'C_KIND_test_pattern(C.NAME').replace(' = rhs', ', rhs)')
        cp.to_f2_call = '*C.NAME' 
  
    #------------------------------
    # Pointer, dim = 2

    if dim == 2:
      s = '(*C.NAME).resize(n1_NAME);\n' + x2 + for1.replace('C.NAME', '(*C.NAME)') + '\n' + \
          '      (*C.NAME)[i].resize(n2_NAME);\n' +  '    *C.NAME << z_NAME;\n'

      cp.to_c2_set    = to_c2_set_pointer.replace('KIND', c_type + '_Matrix')\
                                         .replace('SET', s)
      cp.test_pat     = test_pat_pointer + '    C.NAME = new ' + c_type + '_Matrix(3);\n' + \
                        indent((for1 + '\n    C.NAME[i].resize(2);\n' + c.test_pat)\
                                .replace('C.NAME', '(*C.NAME)'), 2) + '  }\n'
      cp.to_f_cleanup = '  delete z_NAME;\n' 
      cp.to_f_setup   = '''
  int n1_NAME = 0, n2_NAME = 0;
  TYPE* z_NAME = NULL;
  if (C.NAME != NULL) {
    n1_NAME = (*C.NAME).size();
    n2_NAME = (*C.NAME)[0].size();
    z_NAME = new TYPE [(*C.NAME).size()*(*C.NAME)[0].size()];
    matrix_to_vec (*C.NAME, z_NAME);
  }
'''.replace('TYPE', c_type)

      if type == TYPE:
        cp.test_pat   = cp.test_pat.replace('C.NAME', 'C_KIND_test_pattern(C.NAME').replace(' = rhs', ', rhs)')
        cp.to_f2_call = '*C.NAME' 
  
    #------------------------------
    # Pointer, dim = 3

    if dim == 3:
      s = '(*C.NAME).resize(n1_NAME);\n' + x2 + for1.replace('C.NAME', '(*C.NAME)') + '{\n' + \
          '      (*C.NAME)[i].resize(n2_NAME);\n' +  x4 + for2.replace('C.NAME', '(*C.NAME)') + '\n' +\
          x6 + '(*C.NAME)[i][j].resize(n3_NAME);\n' + x4 + '}\n' + x4 + '*C.NAME << z_NAME;\n'

      cp.to_c2_set = '''\
  if (n1_NAME == 0) 
    delete C.NAME;
  else {
    C.NAME = new TYPE_Tensor;
    (*C.NAME).resize(n1_NAME);
    for (int i = 0; i < (*C.NAME).size(); i++) {
      (*C.NAME)[i].resize(n2_NAME);
      for (int j = 0; j < (*C.NAME)[0].size(); j++)
        (*C.NAME)[i][j].resize(n3_NAME);
    }
    *C.NAME << z_NAME;
  }
'''.replace('TYPE', c_type)

      cp.test_pat   = '''\
  if (ix_patt < 3) 
    C.NAME == NULL;
  else {
    C.NAME = new TYPE_Tensor(3);
    for (int i = 0; i < (*C.NAME).size(); i++) {
      (*C.NAME)[i].resize(2);
      for (int j = 0; j < (*C.NAME)[0].size(); j++) {
        (*C.NAME)[i][j].resize(1);
        for (int k = 0; k < (*C.NAME)[0][0].size(); k++)
          {int rhs = 101 + i + 10*(j+1) + 100*(k+1) + XXX + offset; (*C.NAME)[i][j][k] = NNN;}  
      }
    }
  }
'''.replace('TYPE', c_type)

      cp.to_f_cleanup = '  delete z_NAME;\n' 
      cp.to_f_setup   = '''
  int n1_NAME = 0, n2_NAME = 0, n3_NAME = 0;
  TYPE* z_NAME = NULL;
  if (C.NAME != NULL) {
    n1_NAME = (*C.NAME).size();
    n2_NAME = (*C.NAME)[0].size();
    n3_NAME = (*C.NAME)[0][0].size();
    z_NAME = new TYPE [(*C.NAME).size()*(*C.NAME)[0].size()*(*C.NAME)[0][0].size()];
    tensor_to_vec (*C.NAME, z_NAME);
  }
'''.replace('TYPE', c_type)

    #------------------------------

    cp.test_pat = cp.test_pat.replace('NNN', test_value)

#----------------------------------------------------------------------
# CHAR

c_side_trans[CHAR,  0, NOT] = c_side_trans_class()
c_side_trans[CHAR,  0, NOT].c_class    = 'string'
c_side_trans[CHAR,  0, NOT].to_f2_arg  = 'c_Char'
c_side_trans[CHAR,  0, NOT].to_f2_call = 'C.NAME.c_str()'
c_side_trans[CHAR,  0, NOT].to_c2_arg  = 'c_Char z_NAME'
c_side_trans[CHAR,  0, NOT].test_pat    = 'C.NAME.resize(STR_LEN);\n' + test_pat1.replace('NNN', "'a' + rhs % 26")
c_side_trans[CHAR,  0, NOT].constructor = 'NAME()'

# Allocatable components on the C side are the same as pointer components.

for trans in c_side_trans.keys():
  if trans[2] == PTR:
    trans_alloc = (trans[0], trans[1], ALLOC)
    c_side_trans[trans_alloc] = copy.deepcopy(c_side_trans[trans])

##################################################################################
##################################################################################
# Get the list of structs

struct_list_file = 'fortran_structs'
if len(sys.argv) > 1: struct_list_file = sys.argv[1]

module = __import__(struct_list_file)

struct_definitions = []
for name in module.structs:
  struct_definitions.append(struct_def_class(name))

##################################################################################
##################################################################################
# Parse structure definitions

# Examples: 
#  1) "type(abc), pointer :: a(:,:),b(7) = 23 ! Comment"
#  2) "integer ZZZ"
# Notice that only in example 2 is space significant.

# Current restrictions. That is, syntax to avoid:
#   1) Line continuations: '&'
#   2) Dimensions: "integer, dimension(7) :: abc"
#   3) Kind: "integer(kind = 8) abc"
#   4) Variable inits using "," or "(" characters: "real ZZZ(2) = [1, 2]"

re_end_type = re.compile('^\s*end\s*type')  # Match to: 'end type'
re_match1 = re.compile('([,(]|::|\s+)')     # Match to: ',', '::', '(', ' '
re_match2 = re.compile('([=[,(]|::)')       # Match to: ',', '::', '(', '[', '='

for file_name in module.files:
  f_module_file = open(file_name)

  for line in f_module_file:
    split_line = line.lower().split()
    if len(split_line) < 2: continue
    if split_line[0] != 'type': continue

    found = False
    for struct in struct_definitions:
      if struct.f_name != split_line[1]: continue
      found = True
      break

    if not found: continue

    struct.short_name = struct.f_name[:-7]   # Remove '_struct' suffix

    # Now collect the struct components

    for line in f_module_file:
      if re_end_type.match(line): break
      print_debug('\nStart: ' + line.strip())
      base_arg = arg_class()

      part = line.partition('!')

      base_arg.comment = part[2].strip()
      line = part[0].strip()
      if len(line) == 0: continue   # Blank line.
      print_debug('P1: ' + line.strip())

      # Get base_arg.type

      split_line = re_match1.split(line, 1)
      print_debug('P2: ' + str(split_line))
      base_arg.type = split_line.pop(0)

      if split_line[0][0] == ' ': 
        split_line = re_match2.split(split_line[1], 1)
        if split_line[0] == '': split_line.pop(0)
 
      print_debug('P3: ' + str(split_line))

      # Now split_line[0] is a delimiter or component name
      # Add type information if there is more...

      if split_line[0] == '(':
        split_line = split_line[1].partition(')')
        base_arg.kind = split_line[0].strip()
        split_line = re_match2.split(split_line[2].lstrip(), 1)
        if split_line[0] == '': split_line.pop(0)   # EG: "real(rp) :: ..."

      print_debug('P4: ' + str(split_line))

      if split_line[0] == ',':
        split_line = split_line[1].partition('::')

        if split_line[0].strip() == 'allocatable':
          base_arg.pointer_type = ALLOC
        elif split_line[0].strip() == 'pointer':
          base_arg.pointer_type = PTR

        split_line = [split_line[2].lstrip()]

      if split_line[0] == '::': split_line.pop(0)

      # Join split_line into one string so that we are starting from a definite state.

      if len(split_line) > 1: split_line = [''.join(split_line)]
      print_debug('P5: ' + str(split_line))

      # Now len(split_line) = 1 and the first word in split_line[0] is the structure component name.
      # There may be multiple components defined so loop over all instances.

      while True:

        print_debug('L1: ' + str(split_line))

        if len(split_line) > 1:
          print 'Confused parsing of struct component: ' + line.strip()

        split_line = re_match2.split(split_line[0], 1)

        print_debug('L2: ' + str(split_line))

        arg = copy.deepcopy(base_arg)
        arg.name = split_line.pop(0).strip()

        if len(split_line) == 0: 
          struct.arg.append(arg)        
          break

        # Get array bounds

        if split_line[0] == '(':
          split_line = split_line[1].lstrip().partition(')')
          arg.full_array = '(' + split_line[0].strip().replace(' ', '') + ')'
          arg.array = arg.full_array[1:-1].split(',')
          print_debug('L2p1: ' + str(split_line))
          split_line = re_match2.split(split_line[2].lstrip(), 1)
          print_debug('L2p2: ' + str(split_line))
          if split_line[0] == '': split_line.pop(0)  # Needed for EG: "integer aaa(5)"

          if arg.array[0] != ':':   # If has explicit bounds...
            for dim in arg.array:
              if ':' in dim:
                arg.lbound.append(dim.partition(':')[0])
                arg.ubound.append(dim.partition(':')[2])
              else:
                arg.lbound.append('1')
                arg.ubound.append(dim)

        print_debug('L3: ' + str(split_line))

        if len(split_line) == 0: 
          struct.arg.append(arg)        
          break

        # Get initial value

        if split_line[0] == '=':
          split_line = re_match2.split(split_line[1].lstrip(), 1)
          print_debug('L3p1: ' + str(split_line))

          # Combine back if have "(...)" or "[...]" construct as part of init string.

          if len(split_line) > 1 and split_line[1] == '(':
            s2 = split_line[2].partition(')')
            split_line = [split_line[0] + '(' + s2[0] + ')', s2[2]]

          if len(split_line) > 1 and split_line[1] == '[':
            s2 = split_line[2].partition(']')
            split_line = [split_line[0] + '[' + s2[0] + ']', s2[2]]

          arg.init_value = split_line[0]
          if len(split_line) == 1:
            split_line[0] = ''
          else:
            split_line.pop(0)

        print_debug('L4: ' + str(split_line))

        struct.arg.append(arg)        
        if len(split_line) == 0 or split_line[0] == '': break

        if split_line[0] != ',':
          print 'Confused parsing of struct2: ' + line.strip()

        split_line.pop(0)

  # End of parsing

  f_module_file.close()

##################################################################################
##################################################################################
# Add Fortran and C++ side translation info

for struct in struct_definitions:
  for arg in struct.arg:

    # F side translation

    n_dim = len(arg.array)
    p_type = arg.pointer_type

    if (arg.type, n_dim, p_type) not in f_side_trans:
      print 'NO TRANSLATION FOR: ' + struct.short_name + '%' + arg.name + ' [', arg.type + ', ' + str(n_dim) + ', ' + str(p_type) + ']'
      continue

    arg.f_side = copy.deepcopy(f_side_trans[arg.type, n_dim, p_type])
    arg.c_side = copy.deepcopy(c_side_trans[arg.type, n_dim, p_type])

##################################################################################
##################################################################################
# Add array bound info for pointer structure components

for struct in struct_definitions:

  ia = 0
  while ia < len(struct.arg):
    arg = struct.arg[ia]
    ia += 1
    if arg.pointer_type == NOT: continue

    struct.arg.insert(ia, arg_class())
    arg1 = struct.arg[ia]
    ia += 1

    arg1.is_component = False
    arg1.type = 'integer'

    if len(arg.array) == 0: 
      arg1.f_side = copy.deepcopy(f_side_trans[SIZE, 1, NOT])
      arg1.c_side = copy.deepcopy(c_side_trans[SIZE, 1, NOT])
      arg1.name = 'n_' + arg.name
      continue

    if len(arg.array) >= 1:
      arg1.f_side = copy.deepcopy(f_side_trans[SIZE, 1, NOT])
      arg1.c_side = copy.deepcopy(c_side_trans[SIZE, 1, NOT])
      arg1.name = 'n1_' + arg.name

    if len(arg.array) >= 2:
      struct.arg.insert(ia, copy.deepcopy(struct.arg[ia-1]))
      arg2 = struct.arg[ia]
      ia += 1
      arg2.f_side = copy.deepcopy(f_side_trans[SIZE, 2, NOT])
      arg2.c_side = copy.deepcopy(c_side_trans[SIZE, 2, NOT])
      arg2.name = 'n2_' + arg.name

    if len(arg.array) >= 3:
      struct.arg.insert(ia, copy.deepcopy(struct.arg[ia-1]))
      arg3 = struct.arg[ia]
      ia += 1
      arg3.f_side = copy.deepcopy(f_side_trans[SIZE, 3, NOT])
      arg3.c_side = copy.deepcopy(c_side_trans[SIZE, 3, NOT])
      arg3.name = 'n3_' + arg.name

##################################################################################
##################################################################################
# Make name substitutions

for struct in struct_definitions:
  for arg in struct.arg:

    # F side translation

    n_dim = len(arg.array)
    p_type = arg.pointer_type

    arg.c_side.test_pat = arg.c_side.test_pat.replace('STR_LEN', arg.kind)

    if arg.type == 'type':
      kind = arg.kind[:-7]
      arg.f_side.to_f2_trans = arg.f_side.to_f2_trans.replace('KIND', kind)
      arg.f_side.to_f2_extra_var_type = arg.f_side.to_f2_extra_var_type.replace('KIND', kind)
      arg.f_side.test_pat    = arg.f_side.test_pat.replace('KIND', kind)
      arg.c_side.test_pat    = arg.c_side.test_pat.replace('KIND', kind)
      arg.c_side.c_class     = arg.c_side.c_class.replace('KIND', kind)
      arg.c_side.to_c2_set   = arg.c_side.to_c2_set.replace('KIND', kind)
      arg.c_side.to_f_setup  = arg.c_side.to_f_setup.replace('KIND', kind)
      arg.c_side.to_f2_arg   = arg.c_side.to_f2_arg.replace('KIND', kind)
      arg.c_side.to_c2_arg   = arg.c_side.to_c2_arg.replace('KIND', kind)
      arg.c_side.constructor = arg.c_side.constructor.replace('KIND', kind)

    if len(arg.array) >= 1 and p_type == NOT:
      d1 = 1 + int(arg.ubound[0]) - int(arg.lbound[0])
      dim1 = str(d1)
      arg.f_side.to_f2_trans = arg.f_side.to_f2_trans.replace('DIM1', dim1)
      arg.f_side.test_pat    = arg.f_side.test_pat.replace('DIM1', dim1)
      arg.f_side.to_c_var    = arg.f_side.to_c_var.replace('DIM1', dim1)
      arg.f_side.to_c_trans  = arg.f_side.to_c_trans.replace('DIM1', dim1)
      arg.f_side.to_c2_call  = arg.f_side.to_c2_call.replace('DIM1', dim1)
      arg.c_side.constructor = arg.c_side.constructor.replace('DIM1', dim1)
      arg.c_side.to_c2_set   = arg.c_side.to_c2_set.replace('DIM1', dim1)
      arg.c_side.to_f_setup  = arg.c_side.to_f_setup.replace('DIM1', dim1)

    if len(arg.array) >= 2 and p_type == NOT:
      d2 = 1 + int(arg.ubound[1]) - int(arg.lbound[1])
      dim2 = str(d2)
      arg.f_side.to_f2_trans = arg.f_side.to_f2_trans.replace('DIM2', dim2)
      arg.f_side.test_pat    = arg.f_side.test_pat.replace('DIM2', dim2)
      arg.f_side.to_c_var    = arg.f_side.to_c_var.replace('DIM2', dim2)
      arg.f_side.to_c_trans  = arg.f_side.to_c_trans.replace('DIM2', dim2)
      arg.f_side.to_c2_call  = arg.f_side.to_c2_call.replace('SIZE2', str(d1*d2))
      arg.c_side.constructor = arg.c_side.constructor.replace('DIM2', dim2)
      arg.c_side.to_c2_set   = arg.c_side.to_c2_set.replace('DIM2', dim2)
      arg.c_side.to_f_setup  = arg.c_side.to_f_setup.replace('DIM2', dim2)

    if len(arg.array) >= 3 and p_type == NOT:
      d3 = 1 + int(arg.ubound[2]) - int(arg.lbound[2])
      dim3 = str(d3)
      arg.f_side.to_f2_trans = arg.f_side.to_f2_trans.replace('DIM3', dim3)
      arg.f_side.test_pat    = arg.f_side.test_pat.replace('DIM3', dim3)
      arg.f_side.to_c_var    = arg.f_side.to_c_var.replace('DIM3', dim3)
      arg.f_side.to_c_trans  = arg.f_side.to_c_trans.replace('DIM3', dim3)
      arg.f_side.to_c2_call  = arg.f_side.to_c2_call.replace('SIZE3', str(d1*d2*d3))
      arg.c_side.constructor = arg.c_side.constructor.replace('DIM3', dim3)
      arg.c_side.to_c2_set   = arg.c_side.to_c2_set.replace('DIM3', dim3)
      arg.c_side.to_f_setup  = arg.c_side.to_f_setup.replace('DIM3', dim3)

    arg.f_side.bindc_name           = arg.f_side.bindc_name.replace('NAME', arg.name)
    arg.f_side.to_c2_f2_sub_arg     = arg.f_side.to_c2_f2_sub_arg.replace('NAME', arg.name)
    arg.f_side.to_c_var             = arg.f_side.to_c_var.replace('NAME', arg.name)
    arg.f_side.to_c_trans           = arg.f_side.to_c_trans.replace('NAME', arg.name)
    arg.f_side.to_c2_call           = arg.f_side.to_c2_call.replace('NAME', arg.name)
    arg.f_side.to_c2_f2_sub_arg     = arg.f_side.to_c2_f2_sub_arg.replace('NAME', arg.name)
    arg.f_side.bindc_name           = arg.f_side.bindc_name.replace('NAME', arg.name)
    arg.f_side.to_f2_extra_var_name = arg.f_side.to_f2_extra_var_name.replace('NAME', arg.name)
    arg.f_side.to_f2_trans          = arg.f_side.to_f2_trans.replace('NAME', arg.name)
    arg.f_side.equality_test        = arg.f_side.equality_test.replace('NAME', arg.name)
    arg.f_side.test_pat             = arg.f_side.test_pat.replace('NAME', arg.name)
    arg.c_side.constructor          = arg.c_side.constructor.replace('NAME', arg.name)
    arg.c_side.destructor           = arg.c_side.destructor.replace('NAME', arg.name)
    arg.c_side.to_f_setup           = arg.c_side.to_f_setup.replace('NAME', arg.name)
    arg.c_side.to_f_cleanup         = arg.c_side.to_f_cleanup.replace('NAME', arg.name)
    arg.c_side.to_f2_call           = arg.c_side.to_f2_call.replace('NAME', arg.name)
    arg.c_side.to_c2_arg            = arg.c_side.to_c2_arg.replace('NAME', arg.name)
    arg.c_side.to_c2_set            = arg.c_side.to_c2_set.replace('NAME', arg.name)
    arg.c_side.equality_test        = arg.c_side.equality_test.replace('NAME', arg.name)
    arg.c_side.test_pat             = arg.c_side.test_pat.replace('NAME', arg.name)

##################################################################################
##################################################################################
# As a check, write results to file. 

if debug: 
  f_out = open('f_structs.parsed', 'w')
  for struct in struct_definitions:
    f_out.write('******************************************\n')
    f_out.write (struct.f_name + '    ' + str(len(struct.arg)) + '\n')
    for arg in struct.arg:
      f_out.write ('    ' + arg.full_repr() + '\n')

  f_out.close()

#------

n_found = 0
for struct in struct_definitions:
  if struct.short_name != '': n_found = n_found + 1

print 'Number of structs in input list: ' + str(len(struct_definitions))
print 'Number of structs found:         ' + str(n_found)

if len(struct_definitions) != n_found:
  sys.exit('COULD NOT FIND ALL THE STRUCTS! STOPPING HERE!')  

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

''')

f_face.write ('\n'.join(module.use))

f_face.write ('''
use fortran_cpp_utils
use, intrinsic :: iso_c_binding
''')

##############
# ZZZ_to_f interface

for struct in struct_definitions:
  f_face.write ('''
!--------------------------------------------------------------------------

interface 
  subroutine ZZZ_to_f (C, Fp) bind(c)
    import c_ptr
    type(c_ptr), value :: C, Fp
  end subroutine
end interface
'''.replace('ZZZ', struct.short_name))


f_face.write ('contains\n')


##############
# ZZZ_to_c definitions

for struct in struct_definitions:

  s_name = struct.short_name

  f_face.write ('''
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine ZZZ_to_c (Fp, C) bind(c)
!
! Routine to convert a Bmad ZZZ_struct to a C++ C_ZZZ structure
!
! Input:
!   Fp -- type(c_ptr), value :: Input Bmad ZZZ_struct structure.
!
! Output:
!   C -- type(c_ptr), value :: Output C++ C_ZZZ struct.
!-

subroutine ZZZ_to_c (Fp, C) bind(c)

implicit none

interface
'''.replace('ZZZ', s_name))

  import_set = set(['c_ptr'])
  to_c2_call_def = {}

  for arg in struct.arg: 
    bindc_const =  arg.f_side.bindc_type.partition('(')[2].partition(')')[0]
    import_set.add(bindc_const)
    if not arg.f_side.bindc_type in to_c2_call_def: to_c2_call_def[arg.f_side.bindc_type] = []
    to_c2_call_def[arg.f_side.bindc_type].append(arg.f_side.bindc_name)

  f_face.write ('  subroutine ZZZ_to_c2 (C'.replace('ZZZ', s_name))
  for arg in struct.arg: f_face.write (', ' + arg.f_side.to_c2_f2_sub_arg)
  f_face.write (') bind(c)\n')
  f_face.write ('    import ' + ', '.join(import_set) + '\n')
  f_face.write ('    type(c_ptr), value :: C\n')
  for arg_type, args in to_c2_call_def.items():
    f_face.write ('    ' + arg_type + ' :: ' + ', '.join(args) + '\n')

  f_face.write ('''  end subroutine
end interface

type(c_ptr), value :: Fp
type(c_ptr), value :: C
type(ZZZ_struct), pointer :: F
integer jd1, jd2, jd3, lb1, lb2, lb3
'''.replace('ZZZ', s_name))

  for arg in struct.arg:
    if arg.f_side.to_c_var != '': f_face.write (arg.f_side.to_c_var + '\n')

  f_face.write (
'''
!

call c_f_pointer (Fp, F)
''')

  for arg in struct.arg:
    f_face.write (arg.f_side.to_c_trans)

  f_face.write ('call ZZZ_to_c2 (C'.replace('ZZZ', s_name))

  for arg in struct.arg:
    f_face.write (', ' + arg.f_side.to_c2_call)

  f_face.write(''')

end subroutine ZZZ_to_c

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine ZZZ_to_f2 (Fp, ...etc...) bind(c)
!
! Routine used in converting a C++ C_ZZZ structure to a Bmad ZZZ_struct structure.
! This routine is called by ZZZ_to_c and is not meant to be called directly.
!
! Input:
!   ...etc... -- Components of the structure. See the ZZZ_to_f2 code for more details.
!
! Output:
!   Fp -- type(c_ptr), value :: Bmad ZZZ_struct structure.
!-

subroutine ZZZ_to_f2 (Fp'''.replace('ZZZ', struct.short_name))

  for arg in struct.arg:
    f_face.write(', ' + arg.f_side.to_c2_f2_sub_arg)

  f_face.write(''') bind(c)\n

implicit none

type(c_ptr), value :: Fp
type(ZZZ_struct), pointer :: F
integer jd1, jd2, jd3, lb1, lb2, lb3
'''.replace('ZZZ', struct.short_name))

  f2_arg_list = {}
  for arg in struct.arg:
    if not arg.f_side.bindc_type in f2_arg_list: f2_arg_list[arg.f_side.bindc_type] = []
    f2_arg_list[arg.f_side.bindc_type].append(arg.f_side.bindc_name)
    if arg.f_side.to_f2_extra_var_name == '': continue
    if not arg.f_side.to_f2_extra_var_type in f2_arg_list: f2_arg_list[arg.f_side.to_f2_extra_var_type] = []
    f2_arg_list[arg.f_side.to_f2_extra_var_type].append(arg.f_side.to_f2_extra_var_name)

  for arg_type, arg_list in f2_arg_list.items():
    f_face.write(arg_type + ' :: ' + ', '.join(arg_list) + '\n')

  f_face.write('''
call c_f_pointer (Fp, F)
''')

  for arg in struct.arg:
    f_face.write (arg.f_side.to_f2_trans + '\n')

  f_face.write ('''
end subroutine ZZZ_to_f2
'''.replace('ZZZ', struct.short_name))

########################
# End stuff

f_face.write('end module\n')
f_face.close()

##################################################################################
##################################################################################
# Create Fortran struct equality check code

f_equ = open('code/bmad_equality.f90', 'w')

f_equ.write ('module bmad_equality\n\n')
f_equ.write ('\n'.join(module.use))

f_equ.write ('''

interface operator (==)
''')

for i in range(0, len(struct_definitions), 5):
  f_equ.write ('  module procedure ' + ', '.join('eq_' + f.short_name for f in struct_definitions[i:i+5]) + '\n')

f_equ.write ('''end interface

contains
''')

for struct in struct_definitions:
  f_equ.write ('''
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

elemental function eq_ZZZ (f1, f2) result (is_eq)

implicit none

type(ZZZ_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = .true.
'''.replace('ZZZ', struct.short_name))

  for arg in struct.arg:
    if not arg.is_component: continue
    f_equ.write (arg.f_side.equality_test)

  f_equ.write ('\n' + 'end function eq_ZZZ\n'.replace('ZZZ', struct.short_name))
  
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

for struct in struct_definitions:
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

for struct in struct_definitions:
  f_test.write('''
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_ZZZ (ok)

implicit none

type(ZZZ_struct), target :: f_ZZZ, f2_ZZZ
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_ZZZ (c_ZZZ, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_ZZZ
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_ZZZ_test_pattern (f2_ZZZ, 1)

call test_c_ZZZ(c_loc(f2_ZZZ), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_ZZZ_test_pattern (f_ZZZ, 4)
if (f_ZZZ == f2_ZZZ) then
  print *, 'ZZZ: C side convert C->F: Good'
else
  print *, 'ZZZ: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_ZZZ

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_ZZZ (c_ZZZ, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_ZZZ
type(ZZZ_struct), target :: f_ZZZ, f2_ZZZ
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call ZZZ_to_f (c_ZZZ, c_loc(f_ZZZ))

call set_ZZZ_test_pattern (f2_ZZZ, 2)
if (f_ZZZ == f2_ZZZ) then
  print *, 'ZZZ: F side convert C->F: Good'
else
  print *, 'ZZZ: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_ZZZ_test_pattern (f2_ZZZ, 3)
call ZZZ_to_c (c_loc(f2_ZZZ), c_ZZZ)

end subroutine test2_f_ZZZ

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_ZZZ_test_pattern (F, ix_patt)

implicit none

type(ZZZ_struct) F
integer ix_patt, offset, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

'''.replace('ZZZ', struct.short_name))

  for i, arg in enumerate(struct.arg, 1):
    if not arg.is_component: continue
    f_test.write (arg.f_side.test_pat.replace('XXX', str(i)))

  f_test.write('''
end subroutine set_ZZZ_test_pattern
'''.replace('ZZZ', struct.short_name))

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

typedef bool               Bool;
typedef complex<double>    Complex;
typedef double             Real;
typedef int                Int;
typedef char*              Char;

typedef const bool               c_Bool;
typedef const Complex            c_Complex;
typedef const double             c_Real;
typedef const int                c_Int;
typedef const char*              c_Char;

typedef const bool*              c_BoolArr;
typedef const Complex*           c_ComplexArr;
typedef const double*            c_RealArr;
typedef const int*               c_IntArr;

typedef valarray<bool>           Bool_Array;
typedef valarray<Complex>        Complex_Array;
typedef valarray<double>         Real_Array;
typedef valarray<int>            Int_Array;

typedef valarray<Bool_Array>     Bool_Matrix;
typedef valarray<Complex_Array>  Complex_Matrix;
typedef valarray<Real_Array>     Real_Matrix;
typedef valarray<Int_Array>      Int_Matrix;

typedef valarray<Bool_Matrix>      Bool_Tensor;
typedef valarray<Complex_Matrix>   Complex_Tensor;
typedef valarray<Real_Matrix>      Real_Tensor;
typedef valarray<Int_Matrix>       Int_Tensor;

''')

for struct in struct_definitions:
  f_class.write('''
//--------------------------------------------------------------------
// C_ZZZ

class ZZZ_struct {};  // Opaque class for pointers to corresponding fortran structs.

class C_ZZZ {
public:
'''.replace('ZZZ', struct.short_name))

  for arg in struct.arg:
    if not arg.is_component: continue
    f_class.write('  ' + arg.c_side.c_class.replace('ZZZ', struct.short_name) + ' ' + arg.name  + ';\n')

  f_class.write ('''
  C_ZZZ() :
'''.replace('ZZZ', struct.short_name))

  # Constructor

  construct_list = []
  for arg in struct.arg:
    if not arg.is_component: continue
    construct_list.append(arg.c_side.constructor)

  f_class.write ('    ' + ',\n    '.join(construct_list) + '\n')
  f_class.write('    {}\n\n')

  # Destructor

  f_class.write ('  ~C_ZZZ() {\n'.replace('ZZZ', struct.short_name))
  for arg in struct.arg:
    if arg.c_side.destructor == '': continue
    f_class.write ('    ' + arg.c_side.destructor + '\n')

  f_class.write ('  }\n')

  # End class

  f_class.write('''
};   // End Class

extern "C" void ZZZ_to_c (const ZZZ_struct*, C_ZZZ&);
extern "C" void ZZZ_to_f (const C_ZZZ&, ZZZ_struct*);

bool operator== (const C_ZZZ&, const C_ZZZ&);

typedef valarray<C_ZZZ>          C_ZZZ_Array;
typedef valarray<C_ZZZ_Array>    C_ZZZ_Matrix;
typedef valarray<C_ZZZ_Matrix>   C_ZZZ_Tensor;

'''.replace('ZZZ', struct.short_name))

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

template <class T> void matrix_to_vec (const valarray< valarray<T> >& mat, T* vec) {
  int n1 = mat.size();
  if (n1 == 0) return;
  int n2 = mat[0].size();
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      vec[i*n2+j] = mat[i][j];
    }
  }
}

template <class T> void tensor_to_vec (const valarray< valarray< valarray<T> > >& tensor, T* vec) {
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

//---------------------------------------------------------------------------
// Instantiate instances for conversion from array to C++ structure.

template void operator<< (Bool_Array&,  c_Bool*);
template void operator<< (Bool_Matrix&, c_Bool*);

template void operator<< (Real_Array&,  c_Real*);
template void operator<< (Real_Matrix&, c_Real*);
template void operator<< (Real_Tensor&, c_Real*);

template void operator<< (Complex_Array&,  c_Complex*);
template void operator<< (Complex_Matrix&, c_Complex*);
template void operator<< (Complex_Tensor&, c_Complex*);

template void operator<< (Int_Array&,  c_Int*);
template void operator<< (Int_Matrix&, c_Int*);
template void operator<< (Int_Tensor&, c_Int*);

//---------------------------------------------------------------------------
// Instantiate instances for transfer

template void operator<< (Real_Array&,  const Real_Array&);
template void operator<< (Real_Matrix&, const Real_Matrix&);
template void operator<< (Real_Tensor&, const Real_Tensor&);

template void operator<< (Complex_Array&,  const Complex_Array&);
template void operator<< (Complex_Matrix&, const Complex_Matrix&);
template void operator<< (Complex_Tensor&, const Complex_Tensor&);

template void operator<< (Int_Array&,  const Int_Array&);
template void operator<< (Int_Matrix&, const Int_Matrix&);
template void operator<< (Int_Tensor&, const Int_Tensor&);

//---------------------------------------------------------------------------

template void matrix_to_vec (const Bool_Matrix&,     Bool*);
template void matrix_to_vec (const Complex_Matrix&,  Complex*);
template void matrix_to_vec (const Real_Matrix&,     Real*);
template void matrix_to_vec (const Int_Matrix&,      Int*);

template void tensor_to_vec (const Complex_Tensor&,  Complex*);
template void tensor_to_vec (const Real_Tensor&,     Real*);
template void tensor_to_vec (const Int_Tensor&,      Int*);

//---------------------------------------------------------------------------

void void_matrix_to_vec (const valarray< valarray< void** > >& mat, void** vec) {
  int n1 = mat.size();
  if (n1 == 0) return;
  int n2 = mat[0].size();
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      vec[i*n2+j] = mat[i][j];
    }
  }
}

void void_tensor_to_vec (const valarray< valarray< valarray< void** > > >& tensor, void** vec) {
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

for struct in struct_definitions:

  # ZZZ_to_f2
  f_cpp.write('''
//--------------------------------------------------------------------
//--------------------------------------------------------------------
// C_ZZZ

extern "C" void ZZZ_to_c (const ZZZ_struct*, C_ZZZ&);

extern "C" void ZZZ_to_f2 (ZZZ_struct*'''.replace('ZZZ', struct.short_name))

  for arg in struct.arg:
    f_cpp.write (', ' + arg.c_side.to_f2_arg.replace('ZZZ', struct.short_name))

  f_cpp.write(');\n\n')

  # ZZZ_to_f

  f_cpp.write('extern "C" void ZZZ_to_f (const C_ZZZ& C, ZZZ_struct* F) {\n'.replace('ZZZ', struct.short_name))

  for arg in struct.arg:
    f_cpp.write (arg.c_side.to_f_setup)

  f_cpp.write('\n  ZZZ_to_f2 (F'.replace('ZZZ', struct.short_name))

  for arg in struct.arg:
    f_cpp.write (', ' + arg.c_side.to_f2_call)

  f_cpp.write(');\n\n')

  for arg in struct.arg:
    f_cpp.write (arg.c_side.to_f_cleanup)

  f_cpp.write('}\n')

  # ZZZ_to_c2

  f_cpp.write('\n')
  f_cpp.write('extern "C" void ZZZ_to_c2 (C_ZZZ& C'.replace('ZZZ', struct.short_name))

  for arg in struct.arg:
    f_cpp.write (', ' + arg.c_side.to_c2_arg.replace('ZZZ', struct.short_name))

  f_cpp.write(') {\n')

  for arg in struct.arg:
    if not arg.is_component: continue
    f_cpp.write (arg.c_side.to_c2_set + '\n')

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

template <class T> bool is_all_equal (const valarray<T>& vec1, const valarray<T>& vec2) {
  bool is_eq = true;
  if (vec1.size() != vec2.size()) return false;
  for (int i = 0; i < vec1.size(); i++) {
    is_eq = is_eq && (vec1[i] == vec2[i]);
  }
  return is_eq;
}

template <class T> bool is_all_equal (const valarray< valarray<T> >& mat1, const valarray< valarray<T> >& mat2) {
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

template <class T> bool is_all_equal (const valarray< valarray< valarray<T> > >& tensor1, const valarray< valarray< valarray<T> > >& tensor2) {
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

template bool is_all_equal (const Bool_Array&,     const Bool_Array&);
template bool is_all_equal (const Complex_Array&,  const Complex_Array&);
template bool is_all_equal (const Real_Array&,     const Real_Array&);
template bool is_all_equal (const Int_Array&,      const Int_Array&);

template bool is_all_equal (const Bool_Matrix&,     const Bool_Matrix&);
template bool is_all_equal (const Complex_Matrix&,  const Complex_Matrix&);
template bool is_all_equal (const Real_Matrix&,     const Real_Matrix&);
template bool is_all_equal (const Int_Matrix&,      const Int_Matrix&);

template bool is_all_equal (const Complex_Tensor&,  const Complex_Tensor&);
template bool is_all_equal (const Real_Tensor&,     const Real_Tensor&);
template bool is_all_equal (const Int_Tensor&,      const Int_Tensor&);

''')

for struct in struct_definitions:
  f_eq.write ('\n//--------------------------------------------------------------\n\n')
  f_eq.write ('bool operator== (const C_ZZZ& x, const C_ZZZ& y) {'.replace('ZZZ', struct.short_name) + '\n')
  f_eq.write ('  bool is_eq = true;\n')

  for arg in struct.arg:
    if not arg.is_component: continue
    f_eq.write (arg.c_side.equality_test)

  f_eq.write ('  return is_eq;\n')
  f_eq.write ('};\n\n')

  f_eq.write ('template bool is_all_equal (const C_ZZZ_Array&, const C_ZZZ_Array&);\n'.replace('ZZZ', struct.short_name))
  f_eq.write ('template bool is_all_equal (const C_ZZZ_Matrix&, const C_ZZZ_Matrix&);\n'.replace('ZZZ', struct.short_name))

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

for struct in struct_definitions:
  f_test.write ('''
//--------------------------------------------------------------
//--------------------------------------------------------------

extern "C" void test2_f_ZZZ (C_ZZZ&, bool&);

void set_C_ZZZ_test_pattern (C_ZZZ& C, int ix_patt) {

  int rhs, offset = 100 * ix_patt;

'''.replace('ZZZ', struct.short_name))

  for i, arg in enumerate(struct.arg, 1):
    if not arg.is_component: continue
    f_test.write (arg.c_side.test_pat.replace('XXX', str(i)) + '\n')

  f_test.write('''
}

//--------------------------------------------------------------

extern "C" void test_c_ZZZ (ZZZ_struct* F, bool& c_ok) {

  C_ZZZ C, C2;

  c_ok = true;

  ZZZ_to_c (F, C);
  set_C_ZZZ_test_pattern (C2, 1);

  if (C == C2) {
    cout << " ZZZ: C side convert F->C: Good" << endl;
  } else {
    cout << " ZZZ: C SIDE CONVERT F->C: FAILED!" << endl;
    c_ok = false;
  }

  set_C_ZZZ_test_pattern (C2, 2);
  bool c_ok2;
  test2_f_ZZZ (C2, c_ok2);
  if (!c_ok2) c_ok = false;

  set_C_ZZZ_test_pattern (C, 3);
  if (C == C2) {
    cout << " ZZZ: F side convert F->C: Good" << endl;
  } else {
    cout << " ZZZ: F SIDE CONVERT F->C: FAILED!" << endl;
    c_ok = false;
  }

  set_C_ZZZ_test_pattern (C2, 4);
  ZZZ_to_f (C2, F);

}
'''.replace('ZZZ', struct.short_name))

f_test.close()
