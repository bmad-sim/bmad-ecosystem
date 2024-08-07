#!/usr/bin/env python3

# Note: Run this script in the cpp_bmad_interface directory.

# Script to read in Fortran structures and create:
#   Corresponding C++ class
#   Translator between Fortran structure and C++ class
#   Routines to check for equality between instances of a given fortran structure.
#   Routines to check for equality between instances of a given C++ class
#   Program to check the Fortran / C++ translator

# Note: The corresponding C++ class component for a pointer or allocatable Fortran 
# scalar struct component is an array whose length is zero if the Fortran component
# is nullified and whose length is 1 otherwise.

import sys
import shutil
import os
import copy
import re
import textwrap

##################################################################################
##################################################################################
# Init

## master_input_file = 'test_interface_input'   # Used for testing
master_input_file = 'interface_input_params'

n_char_max = 95
debug = False # Change to True to enable printout

##################################################################################
##################################################################################
# Deciding if something is a number

def is_number(s):
  try:
    float(s.replace('d', 'e').replace('D', 'e'))
    return True
  except ValueError:
    return False

##################################################################################
##################################################################################
# For printing of intermediate steps, etc.

def wrap_line (line, indent, cont_char):
  lines = textwrap.wrap(line, n_char_max, initial_indent = indent, subsequent_indent = indent + '    ')
  for i in range(len(lines)-1):
    lines[i] += cont_char + '\n'
  lines[-1] += '\n'
  return ''.join(lines)
 
def print_debug (line):
  if (debug): print(line)

def indent (string, numspace):
  x = ' ' * numspace
  if string[-1] == '\n':
    return x + string[:-1].replace('\n', '\n' + x) + '\n'
  else:
    return x + string.replace('\n', '\n' + x)

##################################################################################
##################################################################################
# Constants

NOT = 'NOT'
PTR = 'PTR'
ALLOC = 'ALLOC'

T = True
F = False

REAL   = 'real'
CMPLX  = 'complex'
INT    = 'integer'
INT8   = 'integer8'
LOGIC  = 'logical'
CHAR   = 'character'
STRUCT = 'type'
SIZE   = 'size'

##################################################################################
##################################################################################
# struct_def_class
# Class for a structure
# Note: The .arg array will contain, in addtion to the structure components, the
# array bounds for allocatable and pointer structure components.

class struct_def_class:
  def __init__(self, f_name = ''):
    self.f_name = f_name # Struct name on Fortran side
    self.short_name = '' # Struct name without trailing '_struct'. Note: C++ name is 'CPP_<short_name>'
    self.cpp_class  = '' # C++ name.
    self.arg = []        # Array of arg_class. List of structrure components + array bound dimensions. 
    self.c_constructor_arg_list = ''
    self.c_constructor_body = '{}'  # Body of the C++ constructor
    self.c_extra_methods = ''       # Additional custom methods

  def __repr__(self):
    return '[name: %s, #arg: %i]' % (self.short_name, len(self.arg))

# arg_class.

class arg_class:

  def __init__(self):
    self.is_component = True   # Is a structure component? If not, then will be array bound.
    self.f_name = ''           # Fortran side name of argument. Will be lower case
    self.c_name = ''           # C++ side name of argument. May be mangled to avoid reserved word conflicts.
    self.type = ''             # Fortran type without '(...)'. EG: 'real', 'type', 'character', etc.
    self.kind = ''             # Fortran kind. EG: '', 'rp', 'coord_struct', etc.
    self.pointer_type = NOT    # NOT, PTR, or ALLOC
    self.array = []            # EG: [':', ':'] or ['0:6', '3']
    self.full_array = ''       # EG: '(:,:)', '(0:6, 3)'
    self.lbound = []
    self.ubound = []
    self.init_value = ''       # Initialization value
    self.comment = ''          # Comment with Fortran structure def.
    self.f_side = 0
    self.c_side = 0

  def __repr__(self):
    return '["%s(%s)", "%s", "%s", %s, "%s"]' % (self.type, self.kind, self.pointer_type, self.f_name, self.array, self.init_value)

  def full_repr(self):
    return '["%s(%s)", "%s", "%s", %s, "%s" %s %s "%s"]' % (self.type, 
              self.kind, self.pointer_type, self.f_name, self.array, self.full_array, 
              self.lbound, self.ubound, self.init_value)

##################################################################################
##################################################################################
# Fortran side translation

class f_side_trans_class:

  def __init__(self):
    self.to_c2_call = ''
    self.to_c2_type = ''
    self.to_c2_name = ''
    self.to_f2_type = ''
    self.to_f2_name = ''
    self.equality_test = 'is_eq = is_eq .and. all(f1%NAME == f2%NAME)\n'
    self.test_pat = 'rhs = XXX + offset; F%NAME = NNN\n'
    self.to_c2_f2_sub_arg = 'z_NAME'
    self.to_f2_trans = 'F%NAME = z_NAME'
    self.to_f2_var = []
    self.to_c_var = []
    self.to_c_trans = ''
    self.size_var = []              # For communicating the size of allocatable and pointer variables

  def __repr__(self):
    return '%s,  %s,  %s :: %s' % (self.to_c2_call, self.to_c2_type, self.to_c2_name, 
                                   self.to_f2_trans, self.to_f2_type, self.to_f2_name)

#--------------------------------------

x2 = ' ' * 2
x4 = ' ' * 4
x6 = ' ' * 6
x8 = ' ' * 8

to_f2_trans_pointer = '''\
if (associated(F%NAME)) then
  if (n1_NAME == 0 .or. any(shape(F%NAME) /= [DIMS])) deallocate(F%NAME)
  if (any(lbound(F%NAME) /= LBOUND)) deallocate(F%NAME)
endif
if (n1_NAME /= 0) then
  call c_f_pointer (z_NAME, f_NAME, [TOTDIM])
  if (.not. associated(F%NAME)) allocate(F%NAME(DIMS))
  SET
else
  if (associated(F%NAME)) deallocate(F%NAME)
endif
'''

equality_test_pointer = '''\
is_eq = is_eq .and. (associated(f1%NAME) .eqv. associated(f2%NAME))
if (.not. is_eq) return
if (associated(f1%NAME)) is_eq = all(shape(f1%NAME) == shape(f2%NAME))
if (.not. is_eq) return
if (associated(f1%NAME)) is_eq = all(f1%NAME == f2%NAME)
'''

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

for type in [REAL, CMPLX, INT, INT8, LOGIC, STRUCT, SIZE]:
  for dim in range(4):
    f_side_trans[type, dim, NOT] = copy.deepcopy(f_side_trans_class())
    f = f_side_trans[type, dim, NOT]

    if type == REAL:
      f.to_c2_type = 'real(c_double)'
      test_value = 'rhs'

    if type == CMPLX:
      f.to_c2_type = 'complex(c_double_complex)'
      test_value = 'cmplx(rhs, 100+rhs)'

    if type == INT:
      f.to_c2_type = 'integer(c_int)'
      test_value = 'rhs'

    if type == INT8:
      f.to_c2_type = 'integer(c_long)'
      test_value = 'rhs'

    if type == LOGIC:
      f.to_c2_type = 'logical(c_bool)'
      test_value = '(modulo(rhs, 2) == 0)'

    if type == STRUCT:
      f.to_c2_type = 'type(c_ptr)'
      test_value = 'NNN'

    if type == SIZE:
      f.to_c2_call       = 'NAME'
      f.to_c2_type       = 'integer(c_int), value'
      f.to_c2_name       = 'NAME'
      f.to_f2_type       = 'integer(c_int), value'
      f.to_f2_name       = 'NAME'
      f.to_f2_trans      = ''
      f.to_c2_f2_sub_arg = 'NAME'
      f.to_c_var         = ['integer(c_int) :: NAME']
      f.test_value       = ''
      continue

    #-----------------------------------------------

    if dim == 0: 
      f.to_c2_name    = 'z_NAME'
      f.equality_test = 'is_eq = is_eq .and. (f1%NAME == f2%NAME)\n'
      f.to_c2_call    = 'F%NAME'
      f.test_pat      = f.test_pat
      if type == LOGIC:
        f.to_c2_call  = 'c_logic(F%NAME)'
        f.to_f2_trans = 'F%NAME = f_logic(z_NAME)'
        f.equality_test = f.equality_test.replace('==', '.eqv.')
      if type == STRUCT:
        f.to_c2_type = 'type(c_ptr), value'
        f.to_c2_call = 'c_loc(F%NAME)'
        f.to_f2_trans = 'call KIND_to_f(z_NAME, c_loc(F%NAME))'
        f.test_pat    = 'call set_KIND_test_pattern (F%NAME, ix_patt)\n'

    #-----------------------------------------------

    if dim == 1:
      f.to_c2_call  = 'fvec2vec(F%NAME, DIM1)'
      f.to_c2_name  = 'z_NAME(*)'
      f.to_f2_trans = 'F%NAME = z_NAME(1:DIM1)'
      f.test_pat    = test_pat1
      if type == LOGIC:
        f.to_f2_trans = 'call vec2fvec (z_NAME, F%NAME)'
        f.equality_test = f.equality_test.replace('==', '.eqv.')
      if type == STRUCT:
        f.to_c2_call = 'z_NAME'
        f.to_f2_trans = jd1_loop + '  call KIND_to_f(z_NAME(jd1), c_loc(F%NAME(jd1+lb1)))\nenddo'
        f.test_pat    = jd1_loop + rhs1 + '  call set_KIND_test_pattern (F%NAME(jd1+lb1), ix_patt+jd1)\n' + 'enddo\n'
        f.to_c_var    = ['type(c_ptr) :: z_NAME(DIM1)']
        f.to_c_trans  = jd1_loop + '  z_NAME(jd1) = c_loc(F%NAME(jd1+lb1))\nenddo\n'

    #-----------------------------------------------

    if dim == 2:
      f.to_f2_trans = 'call vec2mat(z_NAME, F%NAME)'
      f.to_c2_call  = 'mat2vec(F%NAME, DIM2)'
      f.to_c2_name  = 'z_NAME(*)'
      f.test_pat    = test_pat2
      if type == LOGIC:
        f.equality_test = f.equality_test.replace('==', '.eqv.')
      if type == STRUCT:
        f.to_c2_call = 'z_NAME'
        f.to_f2_trans = jd1_loop + jd2_loop + \
                    '  call KIND_to_f(z_NAME(DIM2*(jd1-1) + jd2), c_loc(F%NAME(jd1+lb1,jd2+lb2)))\n' + 'enddo; enddo\n'
        f.test_pat    = jd1_loop + jd2_loop + rhs2 + \
                    '  call set_KIND_test_pattern (F%NAME(jd1+lb1,jd2+lb2), ix_patt+jd1+10*jd2)\n' + 'enddo; enddo\n'
        f.to_c_var    = ['type(c_ptr) :: z_NAME(DIM1*DIM2)']
        f.to_c_trans  = jd1_loop + jd2_loop + '  z_NAME(DIM2*(jd1-1) + jd2) = c_loc(F%NAME(jd1+lb1,jd2+lb2))\n' + \
                                              'enddo; enddo\n'

    #-----------------------------------------------

    if dim == 3:
      f.to_f2_trans = 'call vec2tensor(z_NAME, F%NAME)'
      f.to_c2_call  = 'tensor2vec(F%NAME, DIM3)'
      f.to_c2_name  = 'z_NAME(*)'
      f.test_pat    = test_pat3
      if type == LOGIC:
        f.equality_test = f.equality_test.replace('==', '.eqv.')
      if type == STRUCT:
        f.to_c2_call = 'z_NAME'
        f.to_f2_trans = jd1_loop + jd2_loop + jd3_loop + \
              '  call KIND_to_f(z_NAME(DIM3*DIM2*(jd1-1) + DIM3*(jd2-1) + jd3), c_loc(F%NAME(jd1+lb1,jd2+lb2,jd3+lb3)))\n' + \
              'enddo; enddo; enddo\n'
        f.test_pat    = jd1_loop + jd2_loop + jd3_loop + rhs3 + \
              '  call set_KIND_test_pattern (F%NAME(jd1+lb1,jd2+lb2,jd3+lb3), ix_patt+jd1+10*jd2+100*jd3)\n' + \
              'enddo; enddo; enddo\n'
        f.to_c_var    = ['type(c_ptr) :: z_NAME(DIM1*DIM2*DIM3)']
        f.to_c_trans  = jd1_loop + jd2_loop + jd3_loop + \
              '  z_NAME(DIM3*DIM2*(jd1-1) + DIM3*(jd2-1) + jd3) = c_loc(F%NAME(jd1+lb1,jd2+lb2,jd3+lb3))\n' + \
              'enddo; enddo; enddo\n'

    #-------------------------

    f.test_pat    = f.test_pat.replace('NNN', test_value)
    if f.to_f2_type == '': f.to_f2_type = f.to_c2_type
    if f.to_f2_name == '': f.to_f2_name = f.to_c2_name

    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------
    # Pointers

    f_side_trans[type, dim, PTR] = f_side_trans_class()
    fp = f_side_trans[type, dim, PTR]
    fp.to_c2_type = f.to_c2_type
    fp.to_c2_name = 'z_NAME(*)'
    fp.to_f2_type = 'type(c_ptr), value'
    fp.to_f2_name = 'z_NAME'
    fp.to_f2_var = [f.to_f2_type + ', pointer :: f_NAME(:)'] 

    #---------------------
    # Pointer, dim = 0

    if dim == 0:
      fp.to_f2_var = [f.to_f2_type + ', pointer :: f_NAME'] 
      fp.to_c2_name = 'z_NAME'
      fp.to_c2_call = 'F%NAME'
      fp.to_f2_trans = '''\
if (n_NAME == 0) then                                                                                  
  if (associated(F%NAME)) deallocate(F%NAME)                                                           
else                                                                                                   
  call c_f_pointer (z_NAME, f_NAME)                                                                    
  if (.not. associated(F%NAME)) allocate(F%NAME)                                                       
  F%NAME = f_NAME
endif                                                                                                  
'''

      fp.to_c_trans = '''\
n_NAME = 0
if (associated(F%NAME)) n_NAME = 1
'''

      fp.equality_test = '''
is_eq = is_eq .and. (associated(f1%NAME) .eqv. associated(f2%NAME))
if (.not. is_eq) return
if (associated(f1%NAME)) is_eq = (f1%NAME == f2%NAME)
'''

      test_pat    = '''\
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME)
  rhs = XXX + offset
  SET
endif
'''

      fp.test_pat = test_pat.replace('SET', 'F%NAME = ' + test_value)

      if type == LOGIC:
        fp.equality_test = fp.equality_test.replace('== f', '.eqv. f')
        fp.to_f2_trans = fp.to_f2_trans.replace('= f_NAME', '= f_logic(f_NAME)')
        fp.to_c2_call = 'fscalar2scalar(F%NAME, n_NAME)'
        fp.to_c2_type = 'logical(c_bool)'

      if type == STRUCT:
        fp.to_c2_call = 'c_loc(F%NAME)'
        fp.to_c2_type = 'type(c_ptr), value'
        fp.to_f2_var = ['type(KIND_struct), pointer :: f_NAME'] 
        fp.test_pat = test_pat.replace('SET', 'call set_KIND_test_pattern (F%NAME, ix_patt)')
        fp.to_f2_trans = '''\
if (n_NAME == 0) then
  if (associated(F%NAME)) deallocate(F%NAME)
else
  if (.not. associated(F%NAME)) allocate(F%NAME)
  call KIND_to_f (z_NAME, c_loc(F%NAME))
endif
'''

    #---------------------
    # Pointer, dim = 1

    if dim == 1:

      fp.to_c2_call = 'fvec2vec(F%NAME, n1_NAME)'
      fp.to_f2_trans = to_f2_trans_pointer.replace('DIMS', 'n1_NAME')\
                                          .replace('TOTDIM', 'n1_NAME')\
                                          .replace('SET', 'F%NAME = f_NAME(1:n1_NAME)')
      fp.equality_test = equality_test_pointer

      fp.to_c_trans  = \
'''\
n1_NAME = 0
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME, 1)
endif
'''
      tp1 = '''
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME(-1:1))
'''

      fp.test_pat    = tp1 + x2 + jd1_loop + x2 + rhs1 + x2 + set1.replace('NNN', test_value) + '  enddo\n' + 'endif\n'

      if type == LOGIC:
        fp.equality_test = fp.equality_test.replace('== f', '.eqv. f')
        fp.to_f2_trans = fp.to_f2_trans.replace('F%NAME = f_NAME(1:n1_NAME)', 'call vec2fvec (f_NAME, F%NAME)')

      if type == STRUCT:
        fp.to_c2_call  = 'z_NAME'
        fp.to_c2_type  = 'type(c_ptr)'
        fp.to_c2_name  = 'z_NAME(*)'
        fp.to_f2_type = fp.to_c2_type
        fp.to_f2_name = fp.to_c2_name
        fp.to_c_var    = ['type(c_ptr), allocatable :: z_NAME(:)']
        fp.to_f2_var = []
        fp.test_pat    = tp1 + x2 + jd1_loop + x4 + \
                        'call set_KIND_test_pattern (F%NAME(jd1+lb1), ix_patt+jd1)\n' + '  enddo\n' + 'endif\n'
        ## fp.equality_test = fp.equality_test.replace
        fp.to_c_trans  = ''' \
n1_NAME = 0
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME); lb1 = lbound(F%NAME, 1) - 1
  allocate (z_NAME(n1_NAME))
  do jd1 = 1, n1_NAME
    z_NAME(jd1) = c_loc(F%NAME(jd1+lb1))
  enddo
endif
'''
        fp.to_f2_trans = '''\
if (n1_NAME == 0) then
  if (associated(F%NAME)) deallocate(F%NAME)
else
  if (associated(F%NAME)) then
    if (n1_NAME == 0 .or. any(shape(F%NAME) /= [n1_NAME])) deallocate(F%NAME)
    if (any(lbound(F%NAME) /= LBOUND)) deallocate(F%NAME)
  endif
  if (.not. associated(F%NAME)) allocate(F%NAME(LBOUND:n1_NAME+LBOUND-1))
  do jd1 = 1, n1_NAME
    call KIND_to_f (z_NAME(jd1), c_loc(F%NAME(jd1+LBOUND-1)))
  enddo
endif
'''

    #---------------------
    # Pointer, dim = 2

    if dim == 2:
      fp.to_c2_call = 'mat2vec(F%NAME, n1_NAME*n2_NAME)'
      fp.to_f2_trans = to_f2_trans_pointer.replace('DIMS', 'n1_NAME, n2_NAME')\
                                          .replace('TOTDIM', 'n1_NAME*n2_NAME')\
                                          .replace('SET', 'call vec2mat(f_NAME, F%NAME)')
      fp.equality_test = equality_test_pointer

      fp.to_c_trans  = \
'''\
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME, 1)
  n2_NAME = size(F%NAME, 2)
else
  n1_NAME = 0; n2_NAME = 0
endif
'''
      tp2 = '''
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME(-1:1, 2))
''' 
      fp.test_pat    = tp2 + x2 + jd1_loop + x2 + jd2_loop + x2 + rhs2 + \
                       x2 + set2.replace('NNN', test_value) + '  enddo; enddo\n' + 'endif\n'

      if type == LOGIC:
        fp.equality_test = fp.equality_test.replace('== f', '.eqv. f')

      if type == STRUCT:
        fp.to_c2_call  = 'z_NAME'
        fp.to_c2_type  = 'type(c_ptr)'
        fp.to_c2_name  = 'z_NAME(*)'
        fp.to_f2_type = fp.to_c2_type
        fp.to_f2_name = fp.to_c2_name
        fp.to_c_var    = ['type(c_ptr), allocatable :: z_NAME(:)']
        fp.to_f2_var   = []
        fp.test_pat    = tp2 + x2 + jd1_loop + x2 + jd2_loop + x4 + \
                        'call set_KIND_test_pattern (F%NAME(jd1+lb1,jd2+lb2), ix_patt+jd1+2*jd2)\n' + \
                         '  enddo\n' + '  enddo\n' + 'endif\n'
        ## fp.equality_test = fp.equality_test.replace
        fp.to_c_trans  = '''\
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME, 1); lb1 = lbound(F%NAME, 1) - 1
  n2_NAME = size(F%NAME, 2); lb2 = lbound(F%NAME, 2) - 1
  allocate (z_NAME(n1_NAME * n2_NAME))
  do jd1 = 1, n1_NAME; do jd2 = 1, n2_NAME
    z_NAME(n2_NAME*(jd1-1) + jd2) = c_loc(F%NAME(jd1+lb1, jd2+lb2))
  enddo;  enddo
else
  n1_NAME = 0; n2_NAME = 0
endif
'''
        fp.to_f2_trans = '''\
if (n1_NAME == 0) then
  if (associated(F%NAME)) deallocate(F%NAME)
else
  if (associated(F%NAME)) then
    if (n1_NAME == 0 .or. any(shape(F%NAME) /= [n1_NAME, n2_NAME])) deallocate(F%NAME)
    if (any(lbound(F%NAME) /= LBOUND)) deallocate(F%NAME)
  endif
  if (.not. associated(F%NAME)) allocate(F%NAME(LBOUND:n1_NAME+LBOUND-1, LBOUND:n2_NAME+LBOUND-1))
  do jd1 = 1, n1_NAME
  do jd2 = 1, n2_NAME
    call KIND_to_f (z_NAME(n2_NAME*(jd1-1) + jd2), c_loc(F%NAME(jd1+LBOUND-1,jd2+LBOUND-1)))
  enddo
  enddo
endif
'''

    #---------------------
    # Pointer, dim = 3

    if dim == 3:
      fp.to_c2_call = 'tensor2vec(F%NAME, n1_NAME*n2_NAME*n3_NAME)'
      fp.to_f2_trans = to_f2_trans_pointer.replace('DIMS', 'n1_NAME, n2_NAME, n3_NAME')\
                                          .replace('TOTDIM', 'n1_NAME*n2_NAME*n3_NAME')\
                                          .replace('SET', 'call vec2tensor(f_NAME, F%NAME)')
      fp.equality_test = equality_test_pointer

      fp.to_c_trans  = \
'''\
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME, 1)
  n2_NAME = size(F%NAME, 2)
  n3_NAME = size(F%NAME, 3)
else
  n1_NAME = 0; n2_NAME = 0; n3_NAME = 0
endif
'''

      tp3 = '''\
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME(-1:1, 2, 1))
'''
      fp.test_pat    = tp3 + x2 + jd1_loop + x2 + jd2_loop + x2 + jd3_loop + x2 + rhs3 + \
                       x2 + set3.replace('NNN', test_value) + '  enddo; enddo; enddo\n' + 'endif\n'

      if type == LOGIC:
        fp.equality_test = fp.equality_test.replace('== f', '.eqv. f')

      if type == STRUCT:
        fp.to_c2_call  = 'z_NAME'
        fp.to_c2_type  = 'type(c_ptr)'
        fp.to_c2_name  = 'z_NAME(*)'
        fp.to_f2_type = fp.to_c2_type
        fp.to_f2_name = fp.to_c2_name
        fp.to_c_var    = ['type(c_ptr), allocatable :: z_NAME(:)']
        fp.to_f2_var   = []
        fp.test_pat    = tp3 + x2 + jd1_loop + x2 + jd2_loop + x2 + jd3_loop + x4 + \
                        'call set_KIND_test_pattern (F%NAME(jd1+lb1,jd2+lb2,jd3+lb3), ix_patt+jd1+2*jd2+3*jd3)\n' + \
                         '  enddo\n' + '  enddo\n' + '  enddo\n' + 'endif\n'
        ## fp.equality_test = fp.equality_test.replace
        fp.to_c_trans  = '''\
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME, 1); lb1 = lbound(F%NAME, 1) - 1
  n2_NAME = size(F%NAME, 2); lb2 = lbound(F%NAME, 2) - 1
  n3_NAME = size(F%NAME, 3); lb3 = lbound(F%NAME, 3) - 1
  allocate (z_NAME(n1_NAME * n2_NAME * n3_NAME))
  do jd1 = 1, n1_NAME; do jd2 = 1, n2_NAME; do jd3 = 1, n3_NAME
    z_NAME(n3_NAME*n2_NAME*(jd1-1) + n3_NAME*(jd2-1) + jd3) = c_loc(F%NAME(jd1+lb1, jd2+lb2, jd3+lb3))
  enddo;  enddo; enddo
else
  n1_NAME = 0; n2_NAME = 0; n3_NAME = 0
endif
'''
        fp.to_f2_trans = '''\
if (n1_NAME == 0) then
  if (associated(F%NAME)) deallocate(F%NAME)
else
  if (associated(F%NAME)) then
    if (n1_NAME == 0 .or. any(shape(F%NAME) /= [n1_NAME, n2_NAME, n3_NAME])) deallocate(F%NAME)
    if (any(lbound(F%NAME) /= LBOUND)) deallocate(F%NAME)
  endif
  if (.not. associated(F%NAME)) allocate(F%NAME(LBOUND:n1_NAME+LBOUND-1, LBOUND:n2_NAME+LBOUND-1, LBOUND:n3_NAME+LBOUND-1))
  do jd1 = 1, n1_NAME;  do jd2 = 1, n2_NAME;  do jd3 = 1, n3_NAME
    call KIND_to_f (z_NAME(n3_NAME*n2_NAME*(jd1-1) + n3_NAME*(jd2-1) + jd3), c_loc(F%NAME(jd1+LBOUND-1,jd2+LBOUND-1,jd3+LBOUND-1)))
  enddo;  enddo;  enddo
endif
'''

#---------------------------
# CHAR 0 NOT

f_side_trans[CHAR, 0, NOT] = f_side_trans_class()
fc = f_side_trans[CHAR, 0, NOT] 
fc.to_c2_call    = 'trim(F%NAME) // c_null_char'
fc.to_c2_type    = 'character(c_char)'
fc.to_c2_name    = 'z_NAME(*)'
fc.to_f2_type    = fc.to_c2_type
fc.to_f2_name    = fc.to_c2_name
fc.equality_test = 'is_eq = is_eq .and. (f1%NAME == f2%NAME)\n'
fc.test_pat      = \
        'do jd1 = 1, len(F%NAME)\n  F%NAME(jd1:jd1) = char(ichar("a") + modulo(100+XXX+offset+jd1, 26))\nenddo\n'
fc.to_f2_trans   = 'call to_f_str(z_NAME, F%NAME)'

# CHAR 0 PTR

f_side_trans[CHAR, 0, PTR] = copy.deepcopy(f_side_trans[INT, 0, PTR])
fc = f_side_trans[CHAR, 0, PTR] 
fc.to_c2_type       = 'character(c_char)'
fc.to_c2_name       = 'z_NAME(*)'
fc.to_f2_type    = fc.to_c2_type
fc.to_f2_name    = fc.to_c2_name
fc.to_f2_trans      = '''\
if (n_NAME == 0) then
  if (associated(F%NAME)) deallocate(F%NAME)
else
  if (.not. associated(F%NAME)) allocate(F%NAME)
  call to_f_str(z_NAME, F%NAME)
endif
'''

fc.to_c_var        = ['character(STR_LEN+1), target :: f_NAME']
fc.to_c_trans      = '''\
n_NAME = 0
if (associated(F%NAME)) then
  n_NAME = 1
  f_NAME = trim(F%NAME) // c_null_char 
endif
'''

fc.to_c2_call = 'f_NAME'

fc.test_pat        = '''\
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME)
  do jd1 = 1, len(F%NAME)
    F%NAME(jd1:jd1) = char(ichar("a") + modulo(100+XXX+offset+jd1, 26))
  enddo
endif
'''

# CHAR 1 NOT

f_side_trans[CHAR, 1, NOT] = copy.deepcopy(f_side_trans[STRUCT, 1, NOT])
fc = f_side_trans[CHAR, 1, NOT]
fc.to_c2_name    = 'z_NAME(*)'
fc.to_f2_type    = fc.to_c2_type
fc.to_f2_name    = fc.to_c2_name
fc.to_f2_var     = ['character(c_char), pointer :: f_NAME']
fc.test_pat      = '''\
do jd1 = lbound(F%NAME, 1), ubound(F%NAME, 1)
  do jd = 1, len(F%NAME(jd1))
    F%NAME(jd1)(jd:jd) = char(ichar("a") + modulo(100+XXX+offset+10*jd+jd1, 26))
  enddo
enddo
'''
fc.to_f2_trans = jd1_loop + '''\
  call c_f_pointer (z_NAME(jd1), f_NAME)
  call to_f_str(f_NAME, F%NAME(jd1+lb1))
enddo
''' 
fc.to_c_trans  = jd1_loop + '''\
  a_NAME(jd1) = trim(F%NAME(jd1+lb1)) // c_null_char
  z_NAME(jd1) = c_loc(a_NAME(jd1))
enddo
'''
fc.to_c_var += ['character(STR_LEN+1), target :: a_NAME(DIM1)']

# CHAR 1 PTR

f_side_trans[CHAR, 1, PTR] = copy.deepcopy(f_side_trans[STRUCT, 1, PTR])
fc = f_side_trans[CHAR, 1, PTR] 
fc.to_c2_type       = 'type(c_ptr)'
fc.to_c2_name       = 'z_NAME(*)'
fc.to_f2_type    = fc.to_c2_type
fc.to_f2_name    = fc.to_c2_name
fc.to_f2_var        = ['character(c_char), pointer :: f_NAME']
fc.to_f2_trans      = fc.to_f2_trans.replace('ZZZ', 'string')
fc.test_pat         = '''\
if (ix_patt < 3) then
  if (associated(F%NAME)) deallocate (F%NAME)
else
  if (.not. associated(F%NAME)) allocate (F%NAME(3))
  do jd1 = 1, 3
  do jd = 1, len(F%NAME)
    F%NAME(jd1)(jd:jd) = char(ichar("a") + modulo(100+XXX+offset+10*jd+jd1, 26))
  enddo; enddo
endif
'''

fc.to_f2_trans = '''\
if (n1_NAME == 0) then
  if (associated(F%NAME)) deallocate(F%NAME)
else
  if (associated(F%NAME)) then
    if (n1_NAME == 0 .or. any(shape(F%NAME) /= [n1_NAME])) deallocate(F%NAME)
    if (any(lbound(F%NAME) /= LBOUND)) deallocate(F%NAME)
  endif
  if (.not. associated(F%NAME)) allocate(F%NAME(LBOUND:n1_NAME+LBOUND-1))
  do jd1 = 1, n1_NAME
    call c_f_pointer (z_NAME(jd1), f_NAME)
    call to_f_str(f_NAME, F%NAME(jd1+LBOUND-1))
  enddo
endif
'''

fc.to_c_trans  = '''\
n1_NAME = 0
if (associated(F%NAME)) then
  n1_NAME = size(F%NAME); lb1 = lbound(F%NAME, 1) - 1
  allocate (a_NAME(n1_NAME))
  allocate (z_NAME(n1_NAME))
  do jd1 = 1, n1_NAME
    a_NAME(jd1) = trim(F%NAME(jd1+lb1)) // c_null_char
    z_NAME(jd1) = c_loc(a_NAME(jd1))
  enddo
endif
'''

fc.to_c_var += ['character(STR_LEN+1), allocatable, target :: a_NAME(:)']

#--------------------------------------------------------------------------------------
# Allocatable components are very similar to pointer components
# with the simple replacement of 'allocated' for 'associated'.

for trans in list(f_side_trans.keys()):
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
    self.c_class             = ''        # EG: 'CPP_ele_Array'
    self.c_class_suffix      = ''        # EG: '*'
    self.to_f_setup          = ''
    self.to_f_cleanup        = ''
    self.to_f2_arg           = ''
    self.to_f2_call          = ''
    self.to_c2_arg           = ''
    self.to_c2_set           = '  C.NAME = z_NAME;'
    self.constructor         = 'NAME(VALUE)'
    self.construct_value     = '0'
    self.destructor          = ''
    self.equality_test       = '  is_eq = is_eq && (x.NAME == y.NAME);\n'
    self.test_pat            = '  rhs = XXX + offset; C.NAME = NNN;\n'

  def __repr__(self):
    return '%s,  %s,  %s,  %s' % (self.c_class, self.to_f2_arg, self.to_f2_call, self.to_c2_arg)
    self.size_var = []       # For communicating the size of allocatable and pointer variables

#------------------

test_pat_pointer0 = '''\
  if (ix_patt < 3) 
    C.NAME = NULL;
  else {
'''

test_pat_pointer1 = '''\
  if (ix_patt < 3) 
    C.NAME.resize(0);
  else {
    C.NAME.resize(3);
'''

equality_test_pointer = '''\
  is_eq = is_eq && ((x.NAME == NULL) == (y.NAME == NULL));
  if (!is_eq) return false;
  if (x.NAME != NULL) is_eq = TEST;
'''

for1 = '  for (unsigned int i = 0; i < C.NAME.size(); i++)'
for2 = '  for (unsigned int j = 0; j < C.NAME[0].size(); j++) ' 
for3 = '  for (unsigned int k = 0; k < C.NAME[0][0].size(); k++)'

test_pat1 = for1 + '\n    {int rhs = 101 + i + XXX + offset; C.NAME[i] = NNN;}'
test_pat2 = for1 + for2 + '\n    {int rhs = 101 + i + 10*(j+1) + XXX + offset; C.NAME[i][j] = NNN;}'
test_pat3 = for1 + for2 + for3 + '\n' + x4 +\
                  '{int rhs = 101 + i + 10*(j+1) + 100*(k+1) + XXX + offset; C.NAME[i][j][k] = NNN;}'

c_side_trans = {}

for type in [REAL, CMPLX, INT, INT8, LOGIC, STRUCT, SIZE]:
  for dim in range(4):
    c_side_trans[type, dim, NOT] = c_side_trans_class()
    c = c_side_trans[type, dim, NOT]

    if type == REAL:
      c_type = 'Real'
      c_arg  = 'c_Real'
      test_value = 'rhs'
      c.construct_value = '0.0'

    if type == CMPLX:
      c_type = 'Complex'
      c_arg  = 'c_Complex'
      test_value = 'Complex(rhs, 100+rhs)'
      c.construct_value = '0.0'

    if type == INT:
      c_type = 'Int'
      c_arg  = 'c_Int'
      test_value = 'rhs'
      c.construct_value = '0'

    if type == INT8:
      c_type = 'Int8'
      c_arg  = 'c_Int8'
      test_value = 'rhs'
      c.construct_value = '0'

    if type == LOGIC:
      c_type = 'Bool'
      c_arg  = 'c_Bool'
      test_value = '(rhs % 2 == 0)'
      c.construct_value = 'false'

    if type == STRUCT:
      c_type = 'CPP_KIND'
      c_arg  = 'const CPP_KIND'
      test_value = ''
      c.construct_value = ''

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

      if type == STRUCT:
        c.constructor = 'NAME()'
        c.to_c2_set   = '  KIND_to_c(z_NAME, C.NAME);' 
        c.test_pat    = '  set_CPP_KIND_test_pattern(C.NAME, ix_patt);\n'

    #-------------------------------------------------------

    if dim == 1:
      c.c_class      = c_type + '_ARRAY'
      c.to_f2_arg    = c_arg + 'Arr'
      c.to_f2_call   = '&C.NAME[0]'
      c.to_c2_arg    = c_arg + 'Arr z_NAME'
      c.constructor  = 'NAME(VALUE, DIM1)'
      c.to_c2_set    = '  C.NAME << z_NAME;'
      c.test_pat    = test_pat1
      c.equality_test = '  is_eq = is_eq && is_all_equal(x.NAME, y.NAME);\n'

      if type == STRUCT:
        c.constructor = 'NAME(CPP_KIND_ARRAY(CPP_KIND(), DIM1))'
        c.to_c2_set   = for1 + ' KIND_to_c(z_NAME[i], C.NAME[i]);' 
        c.test_pat    = test_pat1.replace('C.NAME[i] = NNN', 'set_CPP_KIND_test_pattern(C.NAME[i], ix_patt+i+1)')
        c.to_f_setup  = '''\
  const CPP_KIND* z_NAME[DIM1];
  for (int i = 0; i < DIM1; i++) {z_NAME[i] = &C.NAME[i];}
'''

    #-------------------------------------------------------

    if dim == 2:
      c.c_class      = c_type + '_MATRIX'
      c.to_f2_arg    = c_arg + 'Arr'
      c.to_f2_call   = 'z_NAME'
      c.to_c2_arg    = c_arg + 'Arr z_NAME'
      c.constructor  = 'NAME(' + c_type + '_ARRAY(VALUE, DIM2), DIM1)'
      c.to_c2_set    = '  C.NAME << z_NAME;'
      c.test_pat     = test_pat2
      c.to_f_setup  = '  ' + c_type + ' z_NAME[DIM1*DIM2]; matrix_to_vec(C.NAME, z_NAME);\n'
      c.equality_test = '  is_eq = is_eq && is_all_equal(x.NAME, y.NAME);\n'

      if type == STRUCT:
        c.constructor = 'NAME(CPP_KIND_ARRAY(CPP_KIND(), DIM2), DIM1)'
        c.to_c2_set   = for1 + for2 + '\n    {int m = DIM2*i + j; KIND_to_c(z_NAME[m], C.NAME[i][j]);}' 
        c.test_pat    = test_pat2.replace('C.NAME[i][j] = NNN', 'set_CPP_KIND_test_pattern(C.NAME[i][j], ix_patt+i+1+10*(j+1))')
        c.to_f_setup  = '  const CPP_KIND* z_NAME[DIM1*DIM2];\n' + \
                  for1 + for2 +  '\n    {int m = DIM2*i + j; z_NAME[m] = &C.NAME[i][j];}\n'

    #-------------------------------------------------------

    if dim == 3:
      c.c_class      = c_type + '_TENSOR'
      c.to_f2_arg    = c_arg + 'Arr'
      c.to_f2_call   = 'z_NAME'
      c.to_c2_arg    = c_arg + 'Arr z_NAME'
      c.constructor  = 'NAME(' + c_type + '_MATRIX(' + c_type + '_ARRAY(VALUE, DIM3), DIM2), DIM1)'
      c.to_c2_set    = '  C.NAME << z_NAME;'
      c.test_pat     = test_pat3
      c.to_f_setup   = '  ' + c_type + ' z_NAME[DIM1*DIM2*DIM3]; tensor_to_vec(C.NAME, z_NAME);\n'
      c.equality_test = '  is_eq = is_eq && is_all_equal(x.NAME, y.NAME);\n'

      if type == STRUCT:
        c.constructor = 'NAME(CPP_KIND_MATRIX(CPP_KIND_ARRAY(CPP_KIND(), DIM3), DIM2), DIM1)'
        c.to_c2_set   = for1 + for2 + for3 + '\n    {int m = DIM3*DIM2*i + DIM3*j + k; KIND_to_c(z_NAME[m], C.NAME[i][j][k]);}' 
        c.test_pat    = test_pat3.replace('C.NAME[i][j][k] = NNN', \
                                          'set_CPP_KIND_test_pattern(C.NAME[i][j][k], ix_patt+i+1+10*(j+1)+100*(k+1))')
        c.to_f_setup  = '  const CPP_KIND* z_NAME[DIM1*DIM2*DIM3];\n' + \
                  for1 + for2 + for3 + '\n    {int m = DIM3*DIM2*i + DIM3*j + k; z_NAME[m] = &C.NAME[i][j][k];}\n'

    #-------------------------------------------------------

    c.test_pat    = c.test_pat.replace('NNN', test_value)

    if type == STRUCT:
      c.to_c2_arg = 'const Opaque_KIND_class* z_NAME'
      if dim > 0: 
        c.to_f2_arg  = c.to_f2_arg.replace('Arr', '**')
        c.to_c2_arg  = 'const Opaque_KIND_class** z_NAME'
        c.to_f2_call = 'z_NAME'

    #-------------------------------------------------------
    # Pointers

    c_side_trans[type, dim, PTR] = copy.deepcopy(c)
    cp = c_side_trans[type, dim, PTR] 
    cp.destructor    = ''

    if type == STRUCT:
      cp.to_c2_arg     = 'Opaque_KIND_class** z_NAME'
    else:
      cp.to_f2_arg     = c_arg + 'Arr'
      cp.to_c2_arg     = cp.to_f2_arg + ' z_NAME'

    #------------------------------
    # Pointer, dim = 0

    if dim == 0:
      cp.c_class_suffix = '*'
      cp.constructor    = 'NAME(NULL)'
      cp.destructor     = 'delete NAME;'
      cp.test_pat       = test_pat_pointer0 + '    C.NAME = new ' + c_type + ';\n' + \
                                    indent(c.test_pat.replace('C.NAME', '(*C.NAME)'), 2) + '  }\n'
      cp.to_f_setup     = '  unsigned int n_NAME = 0; if (C.NAME != NULL) n_NAME = 1;\n'
      cp.equality_test  = '''\
  is_eq = is_eq && ((x.NAME == NULL) == (y.NAME == NULL));
  if (!is_eq) return false;
  if (x.NAME != NULL) is_eq = (*x.NAME == *y.NAME);
'''


      cp.to_c2_set     = '''\
  if (n_NAME == 0)
    delete C.NAME;
  else {
    C.NAME = new KIND;
    SET
  }
'''.replace('KIND', c_type)

 
      if type == STRUCT:
        cp.test_pat      = cp.test_pat
        cp.to_f2_call    = '*C.NAME' 
        cp.to_c2_arg     = 'Opaque_KIND_class* z_NAME'
        cp.to_c2_set     = cp.to_c2_set.replace('SET', 'KIND_to_c(z_NAME, *C.NAME);')
      else:
        cp.to_c2_set     = cp.to_c2_set.replace('SET', '*C.NAME = *z_NAME;')

    #------------------------------
    # Pointer, dim = 1

    if dim == 1:
      cp.constructor = cp.constructor.replace('DIM1', '0')
      cp.to_f2_call  = 'z_NAME' 
      cp.to_c2_set   = '''
  C.NAME.resize(n1_NAME);
  C.NAME << z_NAME;
'''
      cp.test_pat    = test_pat_pointer1 + indent(c.test_pat, 2) + '  }\n'
      cp.to_f_setup  = '''\
  int n1_NAME = C.NAME.size();
  c_TYPEArr z_NAME = NULL;
  if (n1_NAME > 0) {
    z_NAME = &C.NAME[0];
  }
'''.replace('TYPE', c_type)

      if type == STRUCT:
        cp.constructor = 'NAME(CPP_KIND_ARRAY(CPP_KIND(), 0))'
        cp.test_pat    = test_pat_pointer1 + x2 + for1 + \
            '  {set_CPP_KIND_test_pattern(C.NAME[i], ix_patt+i+1);}\n' + '  }\n'
        cp.to_f_setup  = c.to_f_setup  
        cp.to_f_setup  = '''\
  int n1_NAME = C.NAME.size();
  const CPP_KIND** z_NAME = NULL;
  if (n1_NAME != 0) {
    z_NAME = new const CPP_KIND*[n1_NAME];
    for (int i = 0; i < n1_NAME; i++) z_NAME[i] = &C.NAME[i];
  }
'''
        cp.to_c2_set   = '''\
  C.NAME.resize(n1_NAME);
  for (int i = 0; i < n1_NAME; i++) KIND_to_c(z_NAME[i], C.NAME[i]);
'''
        cp.to_f_cleanup = ' delete[] z_NAME;\n'

    #------------------------------
    # Pointer, dim = 2

    if dim == 2:
      cp.constructor = cp.constructor.replace('DIM1', '0').replace('DIM2', '0')
      cp.to_c2_set   = '''\
  C.NAME.resize(n1_NAME);
  for (int i = 0; i < n1_NAME; i++) C.NAME[i].resize(n2_NAME);
  C.NAME << z_NAME;
'''
      cp.test_pat     = test_pat_pointer1 + \
                        indent((for1 + '\n    C.NAME[i].resize(2);\n' + c.test_pat), 2) + '  }\n'
      cp.to_f_cleanup = '  delete z_NAME;\n' 
      cp.to_f_setup   = '''\
  int n1_NAME = C.NAME.size(), n2_NAME = 0;
  TYPE* z_NAME = NULL;
  if (n1_NAME > 0) {
    n2_NAME = C.NAME[0].size();
    z_NAME = new TYPE [n1_NAME*n2_NAME];
    matrix_to_vec (C.NAME, z_NAME);
  }
'''.replace('TYPE', c_type)
      cp.to_f_cleanup = '  delete[] z_NAME;\n'


      if type == STRUCT:
        cp.constructor = 'NAME(CPP_KIND_ARRAY(CPP_KIND(), 0), 0)'
        cp.test_pat   = test_pat_pointer1 + '''\
    for (unsigned int i = 0; i < C.NAME.size(); i++) {
      C.NAME[i].resize(2);\n
      for (unsigned int j = 0; j < C.NAME[0].size(); j++) {
        set_CPP_KIND_test_pattern(C.NAME[i][j], ix_patt+i+2*j+3);
      }
    }
  }
'''
        cp.to_c2_set   = '''\
  C.NAME.resize(n1_NAME);
  for (int i = 0; i < n1_NAME; i++) {
    C.NAME[i].resize(n2_NAME);
    for (int j = 0; j < n2_NAME; j++) KIND_to_c(z_NAME[n2_NAME*i+j], C.NAME[i][j]);
  }
'''
        cp.to_f_setup   = '''
  int n1_NAME = C.NAME.size(), n2_NAME = 0;
  const TYPE** z_NAME = NULL;
  if (n1_NAME > 0) {
    n2_NAME = C.NAME[0].size();
    z_NAME = new const TYPE* [n1_NAME*n2_NAME];
    for (int i = 0; i < n1_NAME; i++) {
      for (int j = 0; j < n2_NAME; j++) z_NAME[i*n2_NAME + j] = &C.NAME[i][j];}
  }
'''.replace('TYPE', c_type)
  
    #------------------------------
    # Pointer, dim = 3

    if dim == 3:
      cp.constructor = cp.constructor.replace('DIM1', '0').replace('DIM2', '0').replace('DIM3', '0')
      s = 'C.NAME.resize(n1_NAME);\n' + x2 + for1 + '{\n' + \
          '      C.NAME[i].resize(n2_NAME);\n' +  x4 + for2 + '\n' +\
          x6 + 'C.NAME[i][j].resize(n3_NAME);\n' + x4 + '}\n' + x4 + 'C.NAME << z_NAME;\n'

      cp.to_c2_set = '''\
  C.NAME.resize(n1_NAME);
  for (unsigned int i = 0; i < C.NAME.size(); i++) {
    C.NAME[i].resize(n2_NAME);
    for (unsigned int j = 0; j < C.NAME[0].size(); j++)
      C.NAME[i][j].resize(n3_NAME);
  }
  C.NAME << z_NAME;
'''

      cp.test_pat   = '''\
  if (ix_patt < 3) 
    C.NAME.resize(0);
  else {
    C.NAME.resize(3);
    for (unsigned int i = 0; i < C.NAME.size(); i++) {
      C.NAME[i].resize(2);
      for (unsigned int j = 0; j < C.NAME[0].size(); j++) {
        C.NAME[i][j].resize(1);
        for (unsigned int k = 0; k < C.NAME[0][0].size(); k++) {
          int rhs = 101 + i + 10*(j+1) + 100*(k+1) + XXX + offset; C.NAME[i][j][k] = NNN;
        }
      }
    }
  }
'''

      cp.to_f_cleanup = '  delete z_NAME;\n' 
      cp.to_f_setup   = '''
  int n1_NAME = C.NAME.size(), n2_NAME = 0, n3_NAME = 0;
  TYPE* z_NAME = NULL;
  if (n1_NAME > 0) {
    n2_NAME = C.NAME[0].size();
    n3_NAME = C.NAME[0][0].size();
    z_NAME = new TYPE [C.NAME.size()*C.NAME[0].size()*C.NAME[0][0].size()];
    tensor_to_vec (C.NAME, z_NAME);
  }
'''.replace('TYPE', c_type)
      cp.to_f_cleanup = '  delete[] z_NAME;\n'

      if type == STRUCT:
        cp.constructor = 'NAME(CPP_KIND_MATRIX(CPP_KIND_ARRAY(CPP_KIND(), 0), 0), 0)'
        cp.to_c2_set   = '''
  C.NAME.resize(n1_NAME);
  for (int i = 0; i < n1_NAME; i++) {
    C.NAME[i].resize(n2_NAME);
    for (int j = 0; j < n2_NAME; j++) {
      C.NAME[i][j].resize(n3_NAME);
      for (int k = 0; k < n3_NAME; k++) {
        KIND_to_c(z_NAME[n3_NAME*n2_NAME*i+n3_NAME*j+k], C.NAME[i][j][k]);
    } } }
'''
        cp.test_pat   = '''\
  if (ix_patt < 3) 
    C.NAME.resize(0);
  else {
    C.NAME.resize(3);
    for (unsigned int i = 0; i < C.NAME.size(); i++) {
      C.NAME[i].resize(2);
      for (unsigned int j = 0; j < C.NAME[0].size(); j++) {
        C.NAME[i][j].resize(1);
        for (unsigned int k = 0; k < C.NAME[0][0].size(); k++) {
          set_CPP_KIND_test_pattern(C.NAME[i][j][k], ix_patt+i+2*j+3*k+6);
    } } }
  }
'''
        cp.to_f_setup   = '''
  int n1_NAME = C.NAME.size(), n2_NAME = 0, n3_NAME = 0;
  const TYPE** z_NAME = NULL;
  if (n1_NAME > 0) {
    n2_NAME = C.NAME[0].size();
    n3_NAME = C.NAME[0][0].size();
    z_NAME = new const TYPE* [n1_NAME*n2_NAME*n3_NAME];
    for (int i = 0; i < n1_NAME; i++) {
      for (int j = 0; j < n2_NAME; j++) {
        for (int k = 0; k < n3_NAME; k++) {
          z_NAME[i*n2_NAME*n3_NAME + j*n3_NAME + k] = &C.NAME[i][j][k];
    } } }
  }
'''.replace('TYPE', c_type)


    #------------------------------

    cp.test_pat = cp.test_pat.replace('NNN', test_value)

#----------------------------------------------------------------------
# CHAR, 0, NOT

c_side_trans[CHAR,  0, NOT] = c_side_trans_class()
c_side_trans[CHAR,  0, NOT].c_class     = 'string'
c_side_trans[CHAR,  0, NOT].to_f2_arg   = 'c_Char'
c_side_trans[CHAR,  0, NOT].to_f2_call  = 'C.NAME.c_str()'
c_side_trans[CHAR,  0, NOT].to_c2_arg   = 'c_Char z_NAME'
c_side_trans[CHAR,  0, NOT].test_pat    = '  C.NAME.resize(STR_LEN);\n' + test_pat1.replace('NNN', "'a' + rhs % 26")
c_side_trans[CHAR,  0, NOT].constructor = 'NAME()'

# CHAR, 0, PTR

c_side_trans[CHAR,  0, PTR] = copy.deepcopy(c_side_trans[STRUCT, 0, PTR])
cc = c_side_trans[CHAR, 0, PTR] 
cc.c_class       = 'string'
cc.constructor   = 'NAME(NULL)'
cc.destructor    = 'delete NAME;'
cc.to_f2_call    = 'z_NAME'
cc.to_f2_arg     = 'c_Char'
cc.to_f_setup    = '''\
  unsigned int n_NAME = 0;
  const char* z_NAME = NULL;  
  if (C.NAME != NULL) {
    z_NAME = C.NAME->c_str();
    n_NAME = 1;
  }
'''
cc.to_c2_arg     = 'c_Char z_NAME'
cc.to_c2_set     = '''\
  if (n_NAME == 0) 
    delete C.NAME;
  else {
    C.NAME = new string;
    *(C.NAME) = z_NAME;
  }
'''
cc.test_pat    = '''\
  if (ix_patt < 3) 
    C.NAME == NULL;
  else {
    C.NAME = new string(STR_LEN, ' ');
    for (unsigned int i = 0; i < C.NAME->size(); i++) {
      (*C.NAME)[i] = 'a' + (101 + i + XXX + offset) % 26; }
  }
'''

# CHAR, 1, NOT

c_side_trans[CHAR,  1, NOT] = c_side_trans_class()
cc = c_side_trans[CHAR, 1, NOT] 
cc.c_class       = 'String_ARRAY'
cc.to_f2_arg     = 'c_Char*'
cc.to_f_setup    = '''\
  c_Char z_NAME[DIM1];
  for (int i = 0; i < DIM1; i++) {z_NAME[i] = C.NAME[i].c_str();}
'''
cc.to_c2_arg     = 'c_Char* z_NAME'
cc.test_pat      = for1 + ''' {
    C.NAME[i].resize(STR_LEN);
    for (unsigned int j = 0; j < C.NAME[i].size(); j++) 
      {C.NAME[i][j] = 'a' + (101 + i + 10*(j+1) + XXX + offset) % 26;}
  }
'''
cc.constructor   = 'NAME(String_ARRAY(string(), DIM1))'
cc.equality_test = '  is_eq = is_eq && is_all_equal(x.NAME, y.NAME);\n'
cc.to_f2_call    = c_side_trans[STRUCT, 1, NOT].to_f2_call
cc.to_c2_set     = for1 + ' C.NAME[i] = z_NAME[i];'

# CHAR, 1, PTR

c_side_trans[CHAR,  1, PTR] = copy.deepcopy(c_side_trans[STRUCT, 1, PTR])
cc = c_side_trans[CHAR, 1, PTR] 
cc.c_class       = 'String_ARRAY'
cc.constructor   = 'NAME(String_ARRAY(string(), 0))'
cc.destructor    = ''
cc.equality_test = '  is_eq = is_eq && is_all_equal(x.NAME, y.NAME);\n'
cc.to_f2_arg     = 'c_Char*'
cc.to_c2_arg     = 'c_Char* z_NAME'
cc.to_f_setup  = '''\
  int n1_NAME = C.NAME.size();
  c_Char* z_NAME = NULL;
  if (n1_NAME != 0) {
    z_NAME = new c_Char[n1_NAME];
    for (int i = 0; i < n1_NAME; i++) z_NAME[i] = C.NAME[i].c_str();
  }
'''
cc.to_c2_set     = '''\
  C.NAME.resize(n1_NAME);
  for (int i = 0; i < n1_NAME; i++) C.NAME[i] = z_NAME[i];
'''
cc.test_pat    = test_pat_pointer1 + x2 + for1 + '{\n' + \
                x6 + 'C.NAME[i].resize(STR_LEN);\n' + \
                x4 + for2 + "{\n" + \
                x8 + "C.NAME[i][j] = 'a' + (101 + i + 10*(j+1) + XXX + offset) % 26;\n" + \
                x4 + '} }\n' + x2 + '}\n'


# Allocatable components on the C side are the same as pointer components.

for trans in list(c_side_trans.keys()):
  if trans[2] == PTR:
    trans_alloc = (trans[0], trans[1], ALLOC)
    c_side_trans[trans_alloc] = copy.deepcopy(c_side_trans[trans])

##################################################################################
##################################################################################
# Get the list of structs
# See test_interface_input.py (or whatever file is used).

if len(sys.argv) > 1: master_input_file = sys.argv[1]
print(('Input file: ' + master_input_file))

params = __import__(master_input_file)

struct_definitions = []
for name in params.struct_list:
  struct_definitions.append(struct_def_class(name))

##################################################################################
##################################################################################
# Parse structure definitions

# Examples: 
#  1) "type(abc), pointer :: a(:,:),b(7) = 23 ! Comment"
#  2) "integer abc"
# Notice that only in example 2 is space significant.

# Current restrictions. That is, syntax to avoid:
#   1) Line continuations: '&'
#   2) Dimensions: "integer, dimension(7) :: abc"
#   3) Kind: "integer(kind = 8) abc"
#   4) Variable inits using "," or "(" characters: "real abc(2) = [1, 2]"

re_end_type = re.compile('^\s*end\s*type')  # Match to: 'end type'
re_match1 = re.compile('([,(]|::|\s+)')     # Match to: ',', '::', '(', ' '
re_match2 = re.compile('([=[,(]|::)')       # Match to: ',', '::', '(', '[', '='
re_contains = re.compile('^\s*contains')    # Match to: 'contains' (indicating type procedures are defined.)

for file_name in params.struct_def_files:
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
    struct.cpp_class = 'CPP_' + struct.short_name

    # Now collect the struct components

    found_contains_statement = False

    for line in f_module_file:
      if re_end_type.match(line): break
      if re_contains.match(line): found_contains_statement = True
      if found_contains_statement: continue

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
      if base_arg.type == 'integer' and split_line[0][0] == '(': base_arg.type = 'integer8'

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
          print('Confused parsing of struct component: ' + line.strip() + ' in: ' + struct.f_name)

        split_line = re_match2.split(split_line[0], 1)

        print_debug('L2: ' + str(split_line))

        arg = copy.deepcopy(base_arg)
        arg.f_name = split_line.pop(0).strip().lower()

        # Sometimes must avoid reserved words on the C++ side. 
        # This is handled on a case-by-case basis by params.c_side_name_translation

        full_name = struct.f_name + '%' + arg.f_name
        if full_name in params.c_side_name_translation:
          arg.c_name = params.c_side_name_translation[full_name]
        else:
          arg.c_name = arg.f_name

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

          # If have EG: "b(2) = [3, 4], c => null()" need to 
          # combine back "(...)" or "[...]" construct which is part of init string.

          if len(split_line) > 1 and (split_line[1] == '(' or split_line[1] == '['):
            split0 = split_line[0] + split_line[1]
            n_parens = 1
            for ix, char in enumerate(split_line[2]):
              split0 = split0 + char
              if char == '(' or char == '[': n_parens = n_parens + 1
              if char == ')' or char == ']': n_parens = n_parens - 1
              if n_parens == 0: break
            split1 = split_line[2][ix+1:]
            if split1 == '':
              split_line = [split0]
            elif split1[0] == ',':
              split_line = [split0, ',', split1[1:]]
            else:
              print ('?????')
              sys.exit()

          print_debug('L3p2: ' + str(split_line))

          arg.init_value = split_line[0]
          if len(split_line) == 1:
            split_line[0] = ''
          else:
            split_line.pop(0)

        print_debug('L4: ' + str(split_line))

        struct.arg.append(arg)        
        if len(split_line) == 0 or split_line[0] == '': break

        if split_line[0] != ',':
          print('Expected "," while parsing: ' + line.strip()  + ' in: ' + struct.f_name)

        split_line.pop(0)

  # End of parsing

  f_module_file.close()

##################################################################################
##################################################################################
# Add Fortran and C++ side translation info.
# Also throw out any sub-structures that are not to be translated.

for struct in struct_definitions:

  ia = 0
  while ia < len(struct.arg):

    if struct.arg[ia].kind in params.component_no_translate_list or \
                 struct.f_name + '%' + struct.arg[ia].f_name in params.component_no_translate_list: 
      struct.arg.pop(ia)
      continue

    ia += 1

  #--------

  for arg in struct.arg:

      
    # F side translation

    n_dim = len(arg.array)
    p_type = arg.pointer_type

    if (arg.type, n_dim, p_type) not in f_side_trans:
      print('NO TRANSLATION FOR: ' + struct.short_name + '%' + arg.f_name + ' [', arg.type + ', ' + str(n_dim) + ', ' + str(p_type) + ']')
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

    if len(arg.array) == 0: 
      if 'n_' in arg.c_side.to_f_setup: 
        struct.arg.insert(ia, arg_class())
        arg1 = struct.arg[ia]
        ia += 1
        arg1.is_component = False
        arg1.type = 'integer'
        arg1.f_side = copy.deepcopy(f_side_trans[SIZE, 1, NOT])
        arg1.c_side = copy.deepcopy(c_side_trans[SIZE, 1, NOT])
        arg1.f_name = 'n_' + arg.f_name
        arg1.c_name = 'n_' + arg.c_name
      continue

    if len(arg.array) >= 1:
      if 'n1_' not in arg.c_side.to_f_setup: continue
      struct.arg.insert(ia, arg_class())
      arg1 = struct.arg[ia]
      ia += 1
      arg1.is_component = False
      arg1.type = 'integer'
      arg1.f_side = copy.deepcopy(f_side_trans[SIZE, 1, NOT])
      arg1.c_side = copy.deepcopy(c_side_trans[SIZE, 1, NOT])
      arg1.f_name = 'n1_' + arg.f_name
      arg1.c_name = 'n1_' + arg.c_name

    if len(arg.array) >= 2:
      struct.arg.insert(ia, copy.deepcopy(struct.arg[ia-1]))
      arg2 = struct.arg[ia]
      ia += 1
      arg2.f_side = copy.deepcopy(f_side_trans[SIZE, 2, NOT])
      arg2.c_side = copy.deepcopy(c_side_trans[SIZE, 2, NOT])
      arg2.f_name = 'n2_' + arg.f_name
      arg2.c_name = 'n2_' + arg.c_name

    if len(arg.array) >= 3:
      struct.arg.insert(ia, copy.deepcopy(struct.arg[ia-1]))
      arg3 = struct.arg[ia]
      ia += 1
      arg3.f_side = copy.deepcopy(f_side_trans[SIZE, 3, NOT])
      arg3.c_side = copy.deepcopy(c_side_trans[SIZE, 3, NOT])
      arg3.f_name = 'n3_' + arg.f_name
      arg3.c_name = 'n3_' + arg.c_name

##################################################################################
##################################################################################
# Make name substitutions

for struct in struct_definitions:
  print_debug ('\nStruct: ' + str(struct))
  for arg in struct.arg:
    print_debug ('Arg: ' + str(arg))
    n_dim = len(arg.array)
    p_type = arg.pointer_type

    arg.c_side.test_pat            = arg.c_side.test_pat.replace('STR_LEN', arg.kind)
    arg.f_side.to_c_var            = [var.replace('STR_LEN', arg.kind) for var in arg.f_side.to_c_var]

    id_name = struct.short_name + '%' + arg.f_name
    lbound = params.f_side_lbound(id_name)
    arg.f_side.to_f2_trans = arg.f_side.to_f2_trans.replace('LBOUND', lbound)

    if arg.type == 'type':
      kind = arg.kind[:-7]
      arg.f_side.to_f2_trans = arg.f_side.to_f2_trans.replace('KIND', kind)
      arg.f_side.to_f2_var   = [var.replace('KIND', kind) for var in arg.f_side.to_f2_var]
      arg.f_side.test_pat    = arg.f_side.test_pat.replace('KIND', kind)
      arg.c_side.test_pat    = arg.c_side.test_pat.replace('KIND', kind)
      arg.c_side.c_class     = arg.c_side.c_class.replace('KIND', kind)
      arg.c_side.to_c2_set   = arg.c_side.to_c2_set.replace('KIND', kind)
      arg.c_side.to_f_setup  = arg.c_side.to_f_setup.replace('KIND', kind)
      arg.c_side.to_f2_arg   = arg.c_side.to_f2_arg.replace('KIND', kind)
      arg.c_side.to_c2_arg   = arg.c_side.to_c2_arg.replace('KIND', kind)
      arg.c_side.constructor = arg.c_side.constructor.replace('KIND', kind)

    if len(arg.array) >= 1 and p_type == NOT:
      if arg.ubound[0][-1] == '$':
        f_dim1 = arg.ubound[0]
        c_dim1 = 'Bmad::' + arg.ubound[0][0:-1].upper()
        if arg.lbound[0] != '1':
          print('lbound not "1" with parameter upper bound!')
          sys.exit('STOPPING HERE')
      else:
        f_dim1 = str(1 + int(arg.ubound[0]) - int(arg.lbound[0]))
        c_dim1 = f_dim1

      arg.c_side.to_f_setup          = arg.c_side.to_f_setup.replace('DIM1', c_dim1)
      arg.f_side.to_f2_trans         = arg.f_side.to_f2_trans.replace('DIM1', f_dim1)
      arg.f_side.test_pat            = arg.f_side.test_pat.replace('DIM1', f_dim1)
      arg.f_side.to_c_var            = [var.replace('DIM1', f_dim1) for var in arg.f_side.to_c_var]
      arg.f_side.to_c_trans          = arg.f_side.to_c_trans.replace('DIM1', f_dim1)
      arg.f_side.to_c2_call          = arg.f_side.to_c2_call.replace('DIM1', f_dim1)
      arg.c_side.to_c2_set           = arg.c_side.to_c2_set.replace('DIM1', c_dim1)
      arg.c_side.constructor         = arg.c_side.constructor.replace('DIM1', c_dim1)

    if len(arg.array) >= 2 and p_type == NOT:
      d2 = 1 + int(arg.ubound[1]) - int(arg.lbound[1])
      dim2 = str(d2)
      arg.c_side.to_f_setup          = arg.c_side.to_f_setup.replace('DIM2', dim2)
      arg.f_side.to_f2_trans         = arg.f_side.to_f2_trans.replace('DIM2', dim2)
      arg.f_side.test_pat            = arg.f_side.test_pat.replace('DIM2', dim2)
      arg.f_side.to_c_var            = [var.replace('DIM2', dim2) for var in arg.f_side.to_c_var]
      arg.f_side.to_c_trans          = arg.f_side.to_c_trans.replace('DIM2', dim2)
      arg.f_side.to_c2_call          = arg.f_side.to_c2_call.replace('DIM2', f_dim1+'*'+dim2)
      arg.c_side.to_c2_set           = arg.c_side.to_c2_set.replace('DIM2', dim2)
      arg.c_side.constructor         = arg.c_side.constructor.replace('DIM2', dim2)

    if len(arg.array) >= 3 and p_type == NOT:
      d3 = 1 + int(arg.ubound[2]) - int(arg.lbound[2])
      dim3 = str(d3)
      arg.c_side.to_f_setup          = arg.c_side.to_f_setup.replace('DIM3', dim3)
      arg.f_side.to_f2_trans         = arg.f_side.to_f2_trans.replace('DIM3', dim3)
      arg.f_side.test_pat            = arg.f_side.test_pat.replace('DIM3', dim3)
      arg.f_side.to_c_var            = [var.replace('DIM3', dim3) for var in arg.f_side.to_c_var]
      arg.f_side.to_c_trans          = arg.f_side.to_c_trans.replace('DIM3', dim3)
      arg.f_side.to_c2_call          = arg.f_side.to_c2_call.replace('DIM3', f_dim1+'*'+dim2+'*'+dim3)
      arg.c_side.to_c2_set           = arg.c_side.to_c2_set.replace('DIM3', dim3)
      arg.c_side.constructor         = arg.c_side.constructor.replace('DIM3', dim3)

    arg.f_side.to_c_var             = [var.replace('NAME', arg.f_name) for var in arg.f_side.to_c_var]
    arg.f_side.to_c_trans           = arg.f_side.to_c_trans.replace('NAME', arg.f_name)
    arg.f_side.to_c2_call           = arg.f_side.to_c2_call.replace('NAME', arg.f_name)
    arg.f_side.to_c2_f2_sub_arg     = arg.f_side.to_c2_f2_sub_arg.replace('NAME', arg.f_name)
    arg.f_side.to_c2_name           = arg.f_side.to_c2_name.replace('NAME', arg.f_name)
    arg.c_side.to_c2_arg            = arg.c_side.to_c2_arg.replace('NAME', arg.c_name)
    arg.c_side.to_c2_set            = arg.c_side.to_c2_set.replace('NAME', arg.c_name)

    arg.c_side.to_f_setup           = arg.c_side.to_f_setup.replace('NAME', arg.c_name)
    arg.c_side.to_f_cleanup         = arg.c_side.to_f_cleanup.replace('NAME', arg.c_name)
    arg.c_side.to_f2_call           = arg.c_side.to_f2_call.replace('NAME', arg.c_name)

    arg.f_side.to_f2_var            = [var.replace('NAME', arg.f_name) for var in arg.f_side.to_f2_var]
    arg.f_side.to_f2_trans          = arg.f_side.to_f2_trans.replace('NAME', arg.f_name)
    arg.f_side.to_f2_name           = arg.f_side.to_f2_name.replace('NAME', arg.f_name)

    arg.f_side.equality_test        = arg.f_side.equality_test.replace('NAME', arg.f_name)
    arg.f_side.test_pat             = arg.f_side.test_pat.replace('NAME', arg.f_name)

    arg.c_side.equality_test        = arg.c_side.equality_test.replace('NAME', arg.c_name)
    arg.c_side.test_pat             = arg.c_side.test_pat.replace('NAME', arg.c_name)

    arg.c_side.constructor          = arg.c_side.constructor.replace('NAME', arg.c_name)
    arg.c_side.destructor           = arg.c_side.destructor.replace('NAME', arg.c_name)

    # On Fortran side "complex abc(2) = 0" is allowed but on C++ side want "0.0" for init value.
    # Therefore, ignore "0" as an init value.

    if arg.init_value == '':
      pass
    elif arg.init_value == '0':
      pass
    elif arg.init_value[0] == '>':   # Pointer: '=> null()'
      pass 
    elif '_rp' in arg.init_value:
      arg.init_value = arg.init_value.replace('_rp', '')
    elif arg.init_value == '.true.':
      arg.c_side.construct_value = 'true'
    elif arg.init_value == '.false.':
      arg.c_side.construct_value = 'false'
    elif '$' in arg.init_value:
      arg.c_side.construct_value = 'Bmad::' + arg.init_value[:-1].upper()
    elif ('d' in arg.init_value or 'D' in arg.init_value) and is_number(arg.init_value):
      arg.c_side.construct_value = arg.init_value.replace('d', 'e').replace('D', 'e')
    else:
      arg.c_side.construct_value = arg.init_value

    # If there is an array of values, just use first one.
    if len(arg.c_side.construct_value) > 0 and arg.c_side.construct_value[0] == '[':       
      arg.c_side.construct_value = arg.c_side.construct_value[1:].split(',')[0]      

    arg.c_side.constructor = arg.c_side.constructor.replace('VALUE', arg.c_side.construct_value)

##################################################################################
##################################################################################
# Customize the interface code

params.customize(struct_definitions)

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
  if struct.short_name == '': 
    print('NOT FOUND: ' + struct.f_name)
  else:
    n_found = n_found + 1

print('Number of structs in input list: ' + str(len(struct_definitions)))
print('Number of structs found:         ' + str(n_found))

if len(struct_definitions) != n_found:
  sys.exit('COULD NOT FIND ALL THE STRUCTS! STOPPING HERE!')  

struct_names = set([struct.f_name for struct in struct_definitions])

err = False

for struct in struct_definitions:
  for arg in struct.arg:
    if arg.type != STRUCT: continue
    if arg.kind in params.structs_defined_externally: continue
    if arg.kind not in struct_names:
      print(('NO DEFINITION OF STRUCTURE: ' + arg.kind + ' WHICH IS A COMPONENT OF: ' + struct.short_name))
      err = True

if err: sys.exit()

##################################################################################
##################################################################################
# Create Fortran side of interface...

# First the header

f_face = open(params.code_dir + '/bmad_cpp_convert_mod.f90', 'w')

f_face.write ('''
!+
! Fortran side of the Bmad / C++ structure interface.
!
! This file is generated by the Bmad/C++ interface code generation.
! The code generation files can be found in cpp_bmad_interface.
!
! DO NOT EDIT THIS FILE DIRECTLY! 
!-

module bmad_cpp_convert_mod

''')

f_face.write ('\n'.join(params.conversion_use_statements))

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


f_face.write ('\ncontains\n')


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
! Routine to convert a Bmad ZZZ_struct to a C++ CPP_ZZZ structure
!
! Input:
!   Fp -- type(c_ptr), value :: Input Bmad ZZZ_struct structure.
!
! Output:
!   C -- type(c_ptr), value :: Output C++ CPP_ZZZ struct.
!-

subroutine ZZZ_to_c (Fp, C) bind(c)

implicit none

interface
'''.replace('ZZZ', s_name))

  to_c2_call_def = {}

  for arg in struct.arg: 
    if not arg.f_side.to_c2_type in to_c2_call_def: to_c2_call_def[arg.f_side.to_c2_type] = []
    to_c2_call_def[arg.f_side.to_c2_type].append(arg.f_side.to_c2_name)

  line = 'subroutine ZZZ_to_c2 (C'.replace('ZZZ', s_name)
  for arg in struct.arg: 
      line += ', ' + arg.f_side.to_c2_f2_sub_arg
  line += ') bind(c)\n'

  f_face.write ('  !! f_side.to_c2_f2_sub_arg\n')
  f_face.write (wrap_line(line, '  ', ' &'))
  f_face.write ('    import c_bool, c_double, c_ptr, c_char, c_int, c_long, c_double_complex\n')
  f_face.write ('    !! f_side.to_c2_type :: f_side.to_c2_name\n')
  f_face.write ('    type(c_ptr), value :: C\n')
  for arg_type, args in list(to_c2_call_def.items()):
    for i in range(1+(len(args)-1)//7): 
      f_face.write ('    ' + arg_type + ' :: ' + ', '.join(args[i*7:i*7+7]) + '\n')

  f_face.write ('''  end subroutine
end interface

type(c_ptr), value :: Fp
type(c_ptr), value :: C
type(ZZZ_struct), pointer :: F
integer jd, jd1, jd2, jd3, lb1, lb2, lb3
'''.replace('ZZZ', s_name))

  f_face.write ('!! f_side.to_c_var\n')
  for arg in struct.arg:
    for var in arg.f_side.to_c_var: f_face.write (var + '\n')

  f_face.write ( \
'''
!

call c_f_pointer (Fp, F)

''')

  for arg in struct.arg:
    if arg.f_side.to_c_trans == '': continue
    f_face.write ('!! f_side.to_c_trans[' + arg.type + \
                  ', ' + str(len(arg.array)) + ', ' +  arg.pointer_type + ']\n')
    f_face.write (arg.f_side.to_c_trans)

  f_face.write ('\n' + '!! f_side.to_c2_call\n')

  line = 'call ZZZ_to_c2 (C'.replace('ZZZ', s_name)
  for arg in struct.arg:
      line += ', ' + arg.f_side.to_c2_call
  line += ')'
  f_face.write(wrap_line(line, '', ' &'))

  f_face.write('''
end subroutine ZZZ_to_c

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine ZZZ_to_f2 (Fp, ...etc...) bind(c)
!
! Routine used in converting a C++ CPP_ZZZ structure to a Bmad ZZZ_struct structure.
! This routine is called by ZZZ_to_c and is not meant to be called directly.
!
! Input:
!   ...etc... -- Components of the structure. See the ZZZ_to_f2 code for more details.
!
! Output:
!   Fp -- type(c_ptr), value :: Bmad ZZZ_struct structure.
!-

'''.replace('ZZZ', struct.short_name))

  f_face.write ('!! f_side.to_c2_f2_sub_arg\n')
  line = 'subroutine ZZZ_to_f2 (Fp'.replace('ZZZ', struct.short_name)
  for arg in struct.arg:
    line += ', ' + arg.f_side.to_c2_f2_sub_arg
  line += ') bind(c)'
  f_face.write (wrap_line (line, '', ' &'))

  f_face.write('''

implicit none

type(c_ptr), value :: Fp
type(ZZZ_struct), pointer :: F
integer jd, jd1, jd2, jd3, lb1, lb2, lb3
'''.replace('ZZZ', struct.short_name))

  f2_arg_list = {}
  for arg in struct.arg:
    if not arg.f_side.to_f2_type in f2_arg_list: f2_arg_list[arg.f_side.to_f2_type] = []
    f2_arg_list[arg.f_side.to_f2_type].append(arg.f_side.to_f2_name)
    for var in arg.f_side.to_f2_var:
      var_type = var.split('::')[0].strip()
      var_name = var.split('::')[1].strip()
      if not var_type in f2_arg_list: f2_arg_list[var_type] = []
      f2_arg_list[var_type].append(var_name)

  f_face.write ('!! f_side.to_f2_var && f_side.to_f2_type :: f_side.to_f2_name\n')
  for arg_type, arg_list in list(f2_arg_list.items()):
    for i in range(1+(len(arg_list)-1)//7): 
      f_face.write(arg_type + ' :: ' + ', '.join(arg_list[i*7:i*7+7]) + '\n')

  f_face.write('''
call c_f_pointer (Fp, F)

''')

  for arg in struct.arg:
    if arg.f_side.to_f2_trans == '': continue
    f_face.write ('!! f_side.to_f2_trans[' + arg.type + \
                  ', ' + str(len(arg.array)) + ', ' +  arg.pointer_type + ']\n')
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

f_equ = open(params.equality_mod_dir + '/' + params.equality_mod_file + '.f90', 'w')

f_equ.write ('''\
!+
! Module XXX
!
! This module defines a set of functions which overload the equality operator ("==").
! These functions test for equality between instances of a given structure. 
!
! This file is generated as a by product of the Bmad/C++ interface code generation
! The code generation files can be found in cpp_bmad_interface.
!
! DO NOT EDIT THIS FILE DIRECTLY! 
!- 

module XXX
'''.replace('XXX', params.equality_mod_file))

f_equ.write ('\n'.join(params.equality_use_statements))

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
    if struct.f_name + '%' + arg.f_name in params.interface_ignore_list: continue 
    f_equ.write ('!! f_side.equality_test[' + arg.type + ', ' + str(len(arg.array)) + ', ' +  arg.pointer_type + ']\n')
    f_equ.write (arg.f_side.equality_test)

  f_equ.write ('\n' + 'end function eq_ZZZ\n'.replace('ZZZ', struct.short_name))
  
f_equ.write ('end module\n')
f_equ.close()

##################################################################################
##################################################################################
# Create code check main program

if not os.path.exists(params.test_dir): 
    sys.exit ('DIRECTORY DOES NOT EXIST: ' + params.test_dir)
f_test = open(params.test_dir + '/main.f90', 'w')

f_test.write('''
program cpp_bmad_interface_test

use bmad_cpp_test_mod

logical ok, all_ok

!

all_ok = .true.
''')

for struct in struct_definitions:
  f_test.write ('call test1_f_' + struct.short_name + '(ok); if (.not. ok) all_ok = .false.\n')

f_test.write('''
print *
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

f_test = open(params.test_dir + '/bmad_cpp_test_mod.f90', 'w')
f_test.write('''
module bmad_cpp_test_mod

use bmad_cpp_convert_mod
use XXX
'''.replace('XXX', params.equality_mod_file))

f_test.write ('\n'.join(params.test_use_statements) + '\n\n')

f_test.write('contains\n\n')

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
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

'''.replace('ZZZ', struct.short_name))

  for i, arg in enumerate(struct.arg, 1):
    if not arg.is_component: continue
    if struct.f_name + '%' + arg.f_name in params.interface_ignore_list: continue 
    f_test.write ('!! f_side.test_pat[' + arg.type + \
                  ', ' + str(len(arg.array)) + ', ' +  arg.pointer_type + ']\n')

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
// This file is generated as part of the Bmad/C++ interface code generation.
// The code generation files can be found in cpp_bmad_interface.
//
// DO NOT EDIT THIS FILE DIRECTLY! 
//-

#ifndef CPP_BMAD_CLASSES

#include <string>
#include <valarray>
#include <complex>
''')

for line in params.include_header_files:
  f_class.write(line + '\n')

for struct in struct_definitions:
  f_class.write('''
class CPP_ZZZ;
typedef valarray<CPP_ZZZ>          CPP_ZZZ_ARRAY;
typedef valarray<CPP_ZZZ_ARRAY>    CPP_ZZZ_MATRIX;
typedef valarray<CPP_ZZZ_MATRIX>   CPP_ZZZ_TENSOR;
'''.replace('ZZZ', struct.short_name))

#

for struct in struct_definitions:
  f_class.write('''
//--------------------------------------------------------------------
// CPP_ZZZ

class Opaque_ZZZ_class {};  // Opaque class for pointers to corresponding fortran structs.

class CPP_ZZZ {
public:
'''.replace('ZZZ', struct.short_name))

  for arg in struct.arg:
    if not arg.is_component: continue
    f_class.write('  ' + arg.c_side.c_class.replace('ZZZ', struct.short_name) + arg.c_side.c_class_suffix + \
                  ' ' + arg.c_name  + ';\n')

  # Extra methods

  f_class.write(struct.c_extra_methods)

  f_class.write ('''
  CPP_ZZZ(AAA) :
'''.replace('ZZZ', struct.short_name).replace('AAA', struct.c_constructor_arg_list))

  # Constructor

  construct_list = []
  for arg in struct.arg:
    if not arg.is_component: continue
    construct_list.append(arg.c_side.constructor)

  f_class.write ('    ' + ',\n    '.join(construct_list) + '\n')
  f_class.write('    ' + struct.c_constructor_body + '\n\n')

  # Destructor

  f_class.write ('  ~CPP_ZZZ() {\n'.replace('ZZZ', struct.short_name))
  for arg in struct.arg:
    if arg.c_side.destructor == '': continue
    if struct.f_name + '%' + arg.f_name in params.interface_ignore_list: continue 
    f_class.write ('    ' + arg.c_side.destructor + '\n')

  f_class.write ('  }\n')

  # End class

  f_class.write('''
};   // End Class

extern "C" void ZZZ_to_c (const Opaque_ZZZ_class*, CPP_ZZZ&);
extern "C" void ZZZ_to_f (const CPP_ZZZ&, Opaque_ZZZ_class*);

bool operator== (const CPP_ZZZ&, const CPP_ZZZ&);

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

f_cpp = open(params.code_dir + '/cpp_bmad_convert.cpp', 'w')
f_cpp.write('''
//+
// C++ side of the Bmad / C++ structure interface.
//
// This file is generated as part of the Bmad/C++ interface code generation.
// The code generation files can be found in cpp_bmad_interface.
//
// DO NOT EDIT THIS FILE DIRECTLY! 
//-

#include <iostream>
#include "converter_templates.h"
#include "cpp_bmad_classes.h"

''')

for struct in struct_definitions:

  # ZZZ_to_f2
  f_cpp.write('''
//--------------------------------------------------------------------
//--------------------------------------------------------------------
// CPP_ZZZ

extern "C" void ZZZ_to_c (const Opaque_ZZZ_class*, CPP_ZZZ&);

'''.replace('ZZZ', struct.short_name))

  f_cpp.write ('// c_side.to_f2_arg\n')

  line = 'extern "C" void ZZZ_to_f2 (Opaque_ZZZ_class*'.replace('ZZZ', struct.short_name)
  for arg in struct.arg:
    line += ', ' + arg.c_side.to_f2_arg.replace('ZZZ', struct.short_name)
  line += ');'

  f_cpp.write (wrap_line(line, '', ''))


  # ZZZ_to_f

  f_cpp.write('\n')
  f_cpp.write('extern "C" void ZZZ_to_f (const CPP_ZZZ& C, Opaque_ZZZ_class* F) {\n'.replace('ZZZ', struct.short_name))

  for arg in struct.arg:
    if arg.c_side.to_f_setup == '': continue
    f_cpp.write ('  // c_side.to_f_setup[' + arg.type + \
                  ', ' + str(len(arg.array)) + ', ' +  arg.pointer_type + ']\n')
    f_cpp.write (arg.c_side.to_f_setup)

  f_cpp.write('\n')
  f_cpp.write('  // c_side.to_f2_call\n')

  line = 'ZZZ_to_f2 (F'.replace('ZZZ', struct.short_name)
  for arg in struct.arg:
    line += ', ' + arg.c_side.to_f2_call
  line += ');'
  f_cpp.write(wrap_line(line, '  ', ''))

  f_cpp.write('\n')

  for arg in struct.arg:
    if arg.c_side.to_f_cleanup == '': continue
    f_cpp.write ('  // c_side.to_f_cleanup[' + arg.type + \
                  ', ' + str(len(arg.array)) + ', ' +  arg.pointer_type + ']\n')
    f_cpp.write (arg.c_side.to_f_cleanup)

  f_cpp.write('}\n')

  # ZZZ_to_c2

  f_cpp.write('\n')
  f_cpp.write('// c_side.to_c2_arg\n')

  line = 'extern "C" void ZZZ_to_c2 (CPP_ZZZ& C'.replace('ZZZ', struct.short_name)
  for arg in struct.arg:
    line += ', ' + arg.c_side.to_c2_arg.replace('ZZZ', struct.short_name)
  line += ') {'
  f_cpp.write(wrap_line(line, '', ''))

  f_cpp.write('\n')
  for arg in struct.arg:
    if not arg.is_component: continue
    f_cpp.write ('  // c_side.to_c2_set[' + arg.type + \
                  ', ' + str(len(arg.array)) + ', ' +  arg.pointer_type + ']\n')
    f_cpp.write (arg.c_side.to_c2_set + '\n')

  f_cpp.write('}\n')

f_cpp.close()

##################################################################################
##################################################################################
# Create C++ class equality check code

f_eq = open(params.code_dir + '/cpp_equality.cpp', 'w')

f_eq.write('''
//+
// C++ equality functions for Bmad / C++ structure interface.
//
// This file is generated as part of the Bmad/C++ interface code generation.
// The code generation files can be found in cpp_bmad_interface.
//
// DO NOT EDIT THIS FILE DIRECTLY! 
//-

#include <iostream>
#include <stdlib.h>
#include "cpp_bmad_classes.h"

using namespace std;

//---------------------------------------------------

template <class T> bool is_all_equal (const valarray<T>& vec1, const valarray<T>& vec2) {
  bool is_eq = true;
  if (vec1.size() != vec2.size()) return false;
  for (unsigned int i = 0; i < vec1.size(); i++) {
    is_eq = is_eq && (vec1[i] == vec2[i]);
  }
  return is_eq;
}

template <class T> bool is_all_equal (const valarray< valarray<T> >& mat1, const valarray< valarray<T> >& mat2) {
  bool is_eq = true;
  if (mat1.size() != mat2.size()) return false;
  for (unsigned int i = 0; i < mat1.size(); i++) {
    if (mat1[i].size() != mat2[i].size()) return false;
    for (unsigned int j = 0; j < mat1[i].size(); j++) {
      is_eq = is_eq && (mat1[i][j] == mat2[i][j]);
    }
  }
  return is_eq;
};

template <class T> bool is_all_equal (const valarray< valarray< valarray<T> > >& tensor1, const valarray< valarray< valarray<T> > >& tensor2) {
  bool is_eq = true;
  if (tensor1.size() != tensor2.size()) return false;
  for (unsigned int i = 0; i < tensor1.size(); i++) {
    if (tensor1[i].size() != tensor2[i].size()) return false;
    for (unsigned int j = 0; j < tensor1[i].size(); j++) {
      if (tensor1[i][j].size() != tensor2[i][j].size()) return false;
      for (unsigned int k = 0; k < tensor1[i][j].size(); k++) {
        is_eq = is_eq && (tensor1[i][j][k] == tensor2[i][j][k]);
      }
    }
  }
  return is_eq;
};

//---------------------------------------------------

template bool is_all_equal (const Bool_ARRAY&,     const Bool_ARRAY&);
template bool is_all_equal (const Complex_ARRAY&,  const Complex_ARRAY&);
template bool is_all_equal (const Real_ARRAY&,     const Real_ARRAY&);
template bool is_all_equal (const Int_ARRAY&,      const Int_ARRAY&);
template bool is_all_equal (const String_ARRAY&,   const String_ARRAY&);

template bool is_all_equal (const Bool_MATRIX&,     const Bool_MATRIX&);
template bool is_all_equal (const Complex_MATRIX&,  const Complex_MATRIX&);
template bool is_all_equal (const Real_MATRIX&,     const Real_MATRIX&);
template bool is_all_equal (const Int_MATRIX&,      const Int_MATRIX&);

template bool is_all_equal (const Complex_TENSOR&,  const Complex_TENSOR&);
template bool is_all_equal (const Real_TENSOR&,     const Real_TENSOR&);
template bool is_all_equal (const Int_TENSOR&,      const Int_TENSOR&);

''')

for struct in struct_definitions:
  f_eq.write ('\n//--------------------------------------------------------------\n\n')
  f_eq.write ('bool operator== (const CPP_ZZZ& x, const CPP_ZZZ& y) {'.replace('ZZZ', struct.short_name) + '\n')
  f_eq.write ('  bool is_eq = true;\n')

  for arg in struct.arg:
    if not arg.is_component: continue
    if struct.f_name + '%' + arg.f_name in params.interface_ignore_list: continue 
    f_eq.write (arg.c_side.equality_test)

  f_eq.write ('  return is_eq;\n')
  f_eq.write ('};\n\n')

  f_eq.write ('template bool is_all_equal (const CPP_ZZZ_ARRAY&, const CPP_ZZZ_ARRAY&);\n'.replace('ZZZ', struct.short_name))
  f_eq.write ('template bool is_all_equal (const CPP_ZZZ_MATRIX&, const CPP_ZZZ_MATRIX&);\n'.replace('ZZZ', struct.short_name))

f_eq.close()

##################################################################################
##################################################################################
# Create C++ side code check

f_test = open(params.test_dir + '/cpp_bmad_test.cpp', 'w')
f_test.write('''
//+
// C++ classes definitions for Bmad / C++ structure interface.
//
// This file is generated as part of the Bmad/C++ interface code generation.
// The code generation files can be found in cpp_bmad_interface.
//
// DO NOT EDIT THIS FILE DIRECTLY! 
//-

#include <stdio.h>
#include <iostream>
#include "cpp_bmad_classes.h"

using namespace std;
''')

for struct in params.structs_defined_externally:
  head = struct.replace('_struct', '')
  f_test.write('void set_CPP_' + head + '_test_pattern (CPP_' + head + '& C, int ix_patt);\n')

for struct in struct_definitions:
  f_test.write ('''
//--------------------------------------------------------------
//--------------------------------------------------------------

extern "C" void test2_f_ZZZ (CPP_ZZZ&, bool&);

void set_CPP_ZZZ_test_pattern (CPP_ZZZ& C, int ix_patt) {

  int rhs, offset = 100 * ix_patt;

'''.replace('ZZZ', struct.short_name))

  for i, arg in enumerate(struct.arg, 1):
    if not arg.is_component: continue
    if struct.f_name + '%' + arg.f_name in params.interface_ignore_list: continue 
    f_test.write ('  // c_side.test_pat[' + arg.type + \
                  ', ' + str(len(arg.array)) + ', ' +  arg.pointer_type + ']\n')

    f_test.write (arg.c_side.test_pat.replace('XXX', str(i)) + '\n')

  f_test.write('''
}

//--------------------------------------------------------------

extern "C" void test_c_ZZZ (Opaque_ZZZ_class* F, bool& c_ok) {

  CPP_ZZZ C, C2;

  c_ok = true;

  ZZZ_to_c (F, C);
  set_CPP_ZZZ_test_pattern (C2, 1);

  if (C == C2) {
    cout << " ZZZ: C side convert F->C: Good" << endl;
  } else {
    cout << " ZZZ: C SIDE CONVERT F->C: FAILED!" << endl;
    c_ok = false;
  }

  set_CPP_ZZZ_test_pattern (C2, 2);
  bool c_ok2;
  test2_f_ZZZ (C2, c_ok2);
  if (!c_ok2) c_ok = false;

  set_CPP_ZZZ_test_pattern (C, 3);
  if (C == C2) {
    cout << " ZZZ: F side convert F->C: Good" << endl;
  } else {
    cout << " ZZZ: F SIDE CONVERT F->C: FAILED!" << endl;
    c_ok = false;
  }

  set_CPP_ZZZ_test_pattern (C2, 4);
  ZZZ_to_f (C2, F);

}
'''.replace('ZZZ', struct.short_name))

f_test.close()
