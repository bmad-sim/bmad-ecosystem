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
import re

##################################################################################
##################################################################################

def print_debug (line):
  ## print line
  pass

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
    self.type = ''             # Fortran type without '(...)'. EG: 'real', 'type', 'character', etc.
    self.full_type = ''        # Fortran type. EG: 'real(rp)', 'type(coord_struct)', etc.
                               #                      Strings: 'character(<len>)'
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
    return '["%s", "%s", "%s", %s, "%s"]' % (self.full_type, self.pointer_type, self.name, self.array, self.init_value)

  def full_repr(self):
    return '["%s" "%s", "%s", "%s", %s, "%s" %s %s "%s"]' % (self.type, 
              self.full_type, self.pointer_type, self.name, self.array, self.full_array, 
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

# Fortran side translation

class f_side_trans_class:

  def __init__(self, to_c2_arg, bindc_type, bindc_name):
    self.to_c2_arg = to_c2_arg
    self.bindc_type = bindc_type
    self.bindc_name = bindc_name
    self.bindc_const = ''
    self.equal_test = '(f1%NAME == f2%NAME)'
    self.to_f2_trans = 'FP%NAME = NAME'
    self.test_pat = 'FF%NAME = XXX + offset'

  def __repr__(self):
    return '%s,  %s,  %s :: %s' % (self.to_c2_arg, self.bindc_type, self.bindc_name, self.to_f2_trans)


#       Dim  P_type                    to_c2_arg                bindc_type                    bindc_name   equal_test
f_side_trans = {
  (REAL,  0, NOT) : f_side_trans_class('FP%NAME',               'real(c_double)',             'NAME'),
  (REAL,  1, NOT) : f_side_trans_class('FP%NAME',               'real(c_double)',             'NAME(*)'),
  (REAL,  2, NOT) : f_side_trans_class('mat2arr(FP%NAME)',      'real(c_double)',             'NAME(*)'),
  (REAL,  3, NOT) : f_side_trans_class('tensor2arr(FP%NAME)',   'real(c_double)',             'NAME(*)'),
  (CMPLX, 0, NOT) : f_side_trans_class('FP%NAME',               'complex(c_double_complex)',  'NAME'),
  (CMPLX, 1, NOT) : f_side_trans_class('FP%NAME',               'complex(c_double_complex)',  'NAME(*)'),
  (INT,   0, NOT) : f_side_trans_class('FP%NAME',               'integer(c_int)',             'NAME'),
  (INT,   1, NOT) : f_side_trans_class('FP%NAME',               'integer(c_int)',             'NAME(*)'),
  (INT,   2, NOT) : f_side_trans_class('imat2arr(FP%NAME)',     'integer(c_int)',             'NAME(*)'),
  (LOGIC, 0, NOT) : f_side_trans_class('c_logic(FP%NAME)',      'logical(c_bool)',            'NAME'),
  (LOGIC, 1, NOT) : f_side_trans_class('c_logic(FP%NAME)',      'logical(c_bool)',            'NAME(*)'),
  (TYPE,  0, NOT) : f_side_trans_class('c_loc(FP%NAME)',        'type(c_ptr), value',         'NAME'),
  (TYPE,  1, NOT) : f_side_trans_class('c_loc(FP%NAME)',        'type(c_ptr), value',         'NAME(*)'), 
  (CHAR,  1, NOT) : f_side_trans_class('c_str(FP%NAME)',        'character(c_char)',          'NAME(*)')
  }

f_side_trans[REAL,  1, NOT].test_pat    = 'FF%NAME = [(100 + jd1 + XXX + offset, jd1 = 1, size(FF%NAME))]'

f_side_trans[REAL,  2, NOT].to_f2_trans = 'FP%NAME = arr2tensor(NAME, size(FP%NAME, 1), size(FP%NAME, 2))'
f_side_trans[REAL,  2, NOT].test_pat    = 'FF%NAME = reshape([100 + jd1 + XXX + offset, jd1 = 1, size(FF%NAME)], shape(FF%NAME))'

f_side_trans[INT,   1, NOT].test_pat    = 'FF%NAME = [(100 + jd1 + XXX + offset, jd1 = 1, size(FF%NAME))]'

f_side_trans[INT,   2, NOT].to_f2_trans = 'FP%NAME = arr2imat(NAME, size(FP%NAME, 1), size(FP%NAME, 2))'
f_side_trans[INT,   2, NOT].test_pat    = '''do jd1 = lbound(FF%NAME, 1), ubound(FF%NAME, 1); do jd2 = lbound(FF%NAME, 2), ubound(FF%NAME, 2)
  FF%NAME = 100 + jd1 + 10*jd2 + XXX + offset
enddo; enddo'''

f_side_trans[CMPLX, 1, NOT].test_pat    = 'FF%NAME = [(cmplx(100 + jd1 + XXX + offset, 200 + jd1 + XXX + offset), jd1 = 1, size(FF%NAME))]'

for key, f in f_side_trans.items(): 
  f.bindc_const = f.bindc_type.partition('(')[2].partition(')')[0]
  if key[1] != 0: f.equal_test = 'all(f1%NAME == f2%NAME)'

#############################################################

class c_side_trans_class:

  def __init__(self, c_class, to_f2_arg, to_f2_call, to_c2_arg):
    self.c_class = c_class
    self.to_f2_arg = to_f2_arg
    self.to_f2_call = to_f2_call
    self.to_c2_arg = to_c2_arg
    self.to_f_setup = ''
    self.to_c2_set = 'C.NAME = NAME;'
    self.constructor = 'NAME(0)'
    self.equal_test = '(x.NAME == y.NAME)'
    self.test_pat = 'C.NAME = XXX + offset;'

  def __repr__(self):
    return '%s,  %s,  %s,  %s' % (self.c_class, self.to_f2_arg, self.to_f2_call, self.to_c2_arg)



#        Dim  P_type                    c_class            to_f2_arg          to_f2_call         to_c2_arg
c_side_trans = {
  (REAL,   0, NOT) : c_side_trans_class('double',          'Real&',           'C.NAME',          'Real& NAME'),
  (REAL,   1, NOT) : c_side_trans_class('Real_Array',      'RealArr',         '&C.NAME[0]',      'RealArr NAME'),
  (REAL,   2, NOT) : c_side_trans_class('Real_Matrix',     'RealArr',         'NAME',            'RealArr NAME'),
  (CMPLX,  0, NOT) : c_side_trans_class('Dcomplex',        'Dcomplexd&',      'C.NAME',          'Dcomplex& NAME'),
  (CMPLX,  1, NOT) : c_side_trans_class('Dcomplex_Array',  'DcomplexArr',     '&C.NAME[0]',      'DcomplexArr NAME'),
  (CMPLX,  2, NOT) : c_side_trans_class('Dcomplex_Matrix', 'DcomplexArr',     'NAME',            'DcomplexArr NAME'),
  (INT,    0, NOT) : c_side_trans_class('int',             'Int&',            'C.NAME',          'Int& NAME'),
  (INT,    1, NOT) : c_side_trans_class('Int_Array',       'IntArr',          '&C.NAME[0]',      'IntArr NAME'),
  (INT,    2, NOT) : c_side_trans_class('Int_Matrix',      'IntArr',          'NAME',            'IntArr NAME'),
  (LOGIC,  0, NOT) : c_side_trans_class('bool',            'Bool&',           'C.NAME',          'Bool& NAME'),
  (LOGIC,  1, NOT) : c_side_trans_class('Bool_Array',      'BoolArr',         'C.NAME',          'BoolArr NAME'),
  (TYPE,   0, NOT) : c_side_trans_class('C_zzz',           'const C_zzz&',    '&C.NAME',         'const C_zzz& NAME'),
  (TYPE,   1, NOT) : c_side_trans_class('C_zzz_array',     'const C_zzz&',    'C.NAME[0]',       'const C_zzz& NAME'),
  (CHAR,   1, NOT) : c_side_trans_class('string',          'Char',            'C.NAME',          'Char NAME')
  }

c_side_trans[REAL,  1, NOT].constructor = 'NAME(0.0, DIM1)'
c_side_trans[REAL,  1, NOT].to_c2_set   = 'C.NAME = Real_Array(NAME, DIM1);'
c_side_trans[REAL,  1, NOT].test_pat    = 'for (int i = 0; i < C.NAME.size(); i++) C.NAME[i] = 101 + i + XXX + offset;'

c_side_trans[CMPLX, 1, NOT].constructor = 'NAME(0.0, DIM1)'
c_side_trans[CMPLX, 1, NOT].to_c2_set   = 'C.NAME = Dcomplex_Array(NAME, DIM1);'
c_side_trans[CMPLX, 1, NOT].test_pat    = 'for (int i = 0; i < C.NAME.size(); i++) C.NAME[i] = Dcomplex(101 + i + XXX + offset, 201 + i + XXX + offset);'

c_side_trans[INT,   1, NOT].constructor = 'NAME(DIM1)'
c_side_trans[INT,   1, NOT].to_c2_set   = 'C.NAME = Int_Array(NAME, DIM1);'
c_side_trans[INT,   1, NOT].test_pat    = 'for (int i = 0; i < C.NAME.size(); i++) C.NAME[i] = 101 + i + XXX + offset;'

c_side_trans[INT,   2, NOT].constructor = 'NAME(Int_Array(0, DIM2), DIM1)'
c_side_trans[INT,   2, NOT].to_c2_set   = 'C.NAME = Int_Array(NAME, DIM1);'
c_side_trans[INT,   2, NOT].test_pat    = 'for (int i = 0; i < C.NAME.size(); i++); for (int j = 0; i < C.NAME[0].size(); J++ C.NAME[i][j] = 101 + i + 10*(j+1) + XXX + offset;'
c_side_trans[INT,   2, NOT].to_f_setup  = 'int NAME[DIM1]; matrix_to_vec(C.NAME, NAME);\n'

c_side_trans[LOGIC, 1, NOT].constructor = 'NAME(DIM1)'
c_side_trans[LOGIC, 1, NOT].to_c2_set   = 'C.NAME = Int_Array(NAME, DIM1);'
c_side_trans[LOGIC, 1, NOT].test_pat    = 'for (int i = 0; i < C.NAME.size(); i++) C.NAME[i] = 101 + i + XXX + offset;'

for key, c in c_side_trans.items(): 
  if key[1] == 1: c.equal_test = 'is_all_equal_vec(x.NAME, y.NAME)'
  if key[1] == 2: c.equal_test = 'is_all_equal_mat(x.NAME, y.NAME)'
  if key[1] == 3: c.equal_test = 'is_all_equal_tensor(x.NAME, y.NAME)'

##################################################################################
##################################################################################
# Get the list of structs

struct_list_file = 'scripts/fortran_structs.list'
if len(sys.argv) > 1: struct_list_file = sys.argv[1]

f_struct_list_file = open(struct_list_file)

struct_def = {}
f_struct_list = []
f_module_files = []
f_use_list = []

for line in f_struct_list_file:
  line = line.strip()
  if len(line) == 0: continue
  split_line = line.split()

  if split_line[0] == 'FILE:':
    f_module_files.append(split_line[1])
  elif split_line[0] == 'USE:':
    f_use_list.append('use ' + split_line[1])
  else:
    f_struct_list.append(split_line)
    struct_def[line.split()[0]] = struct_def_class()

f_struct_list_file.close()

##################################################################################
##################################################################################
# Parse structure definitions

# Examples: 
#  1) "type (abc), pointer :: a(:,:),b(7) = 23 ! Comment"
#  2) "integer zzz"
# Notice that only in example 2 is space significant.

# Current restrictions. That is, syntax to avoid:
#   1) Line continuations: '&'
#   2) Dimensions: "integer, dimension(7) :: abc"
#   3) Kind: "integer(kind = 8) abc"
#   4) Variable inits using "," or "(" characters: "real zzz(2) = [1, 2]"

re_end_type = re.compile('^\s*end\s*type')  # Match to: 'end type'
re_match1 = re.compile('([,(]|::|\s+)')     # Match to: ',', '::', '(', ' '
re_match2 = re.compile('([,(]|::|=)')       # Match to: ',', '::', '(', '='

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
      if re_end_type.match(line): break
      print_debug('\nStart: ' + line.strip())
      base_var = var_class()

      part = line.partition('!')

      base_var.comment = part[2].strip()
      line = part[0].strip()
      if len(line) == 0: continue   # Blank line.
      print_debug('P1: ' + line.strip())

      # Get base_var.type

      split_line = re_match1.split(line, 1)
      print_debug('P2: ' + str(split_line))
      base_var.type = split_line.pop(0)
      base_var.full_type = base_var.type

      if split_line[0][0] == ' ': 
        split_line = re_match2.split(split_line[1], 1)
        if split_line[0] == '': split_line.pop(0)
 
      print_debug('P3: ' + str(split_line))

      # Now split_line[0] is a delimiter or variable name
      # Add type information if there is more...

      if split_line[0] == '(':
        split_line = split_line[1].partition(')')
        base_var.full_type = base_var.type + '(' + split_line[0].strip() + ')'
        split_line = re_match2.split(split_line[2].lstrip(), 1)
        if split_line[0] == '': split_line.pop(0)   # EG: "real(rp) :: ..."

      print_debug('P4: ' + str(split_line))

      if split_line[0] == ',':
        split_line = split_line[1].partition('::')

        if split_line[0].strip() == 'allocatable':
          base_var.pointer_type = ALLOC
        elif split_line[0].strip() == 'pointer':
          base_var.pointer_type = PTR

        split_line = [split_line[2].lstrip()]

      if split_line[0] == '::': split_line.pop(0)

      print_debug('P5: ' + str(split_line))

      # Now the first word in split_line[0] is the variable name
      # There may be multiple variables defined so loop over all instances.

      while True:

        print_debug('L1: ' + str(split_line))

        if len(split_line) > 1:
          print 'Confused parsing of struct component: ' + line.strip()

        split_line = re_match2.split(split_line[0], 1)

        print_debug('L2: ' + str(split_line))

        var = copy.deepcopy(base_var)
        var.name = split_line.pop(0).strip()

        if len(split_line) == 0: 
          struct.var.append(var)        
          break

        # Get array bounds

        if split_line[0] == '(':
          split_line = split_line[1].lstrip().partition(')')
          var.full_array = '(' + split_line[0].strip().replace(' ', '') + ')'
          var.array = var.full_array[1:-1].split(',')
          print_debug('L2p1: ' + str(split_line))
          split_line = re_match2.split(split_line[2].lstrip(), 1)
          print_debug('L2p2: ' + str(split_line))
          if split_line[0] == '': split_line.pop(0)  # Needed for EG: "integer aaa(5)"

          if var.array[0] != ':':   # If has explicit bounds...
            for dim in var.array:
              if ':' in dim:
                var.lbound.append(dim.partition(':')[0])
                var.ubound.append(dim.partition(':')[2])
              else:
                var.lbound.append('1')
                var.ubound.append(dim)

        print_debug('L3: ' + str(split_line))

        if len(split_line) == 0: 
          struct.var.append(var)        
          break

        # Get initial value

        if split_line[0] == '=':
          split_line = re_match2.split(split_line[1].lstrip(), 1)
          var.init_value = split_line[0]
          if len(split_line) == 1:
            split_line[0] = ''
          else:
            split_line.pop(0)

        print_debug('L4: ' + str(split_line))

        struct.var.append(var)        
        if len(split_line) == 0 or split_line[0] == '': break

        if split_line[0] != ',':
          print 'Confused parsing of struct2: ' + line.strop()

        split_line.pop(0)

  # End of parsing

  f_module_file.close()

##################################################################################
##################################################################################
# As a check, write results to file. 

f_out = open('f_structs.parsed', 'w')

n_found = 0
for key, struct in struct_def.items():
  if struct.short_name != '': n_found = n_found + 1
  f_out.write('******************************************\n')
  f_out.write (key + '    ' + str(len(struct.var)) + '\n')
  for var in struct.var:
    f_out.write ('    ' + var.full_repr() + '\n')

f_out.close()


print 'Number of structs in input list: ' + str(len(struct_def))
print 'Number of structs found:         ' + str(n_found)

##################################################################################
##################################################################################
# Add translation info and make some name substitutions

for struct in struct_def.values():
  for var in struct.var:

    # F side translation

    n_dim = len(var.array)
    p_type = var.pointer_type

    if (var.type, n_dim, p_type) not in f_side_trans:
      print 'NO TRANSLATION FOR: ' + struct.short_name + '%' + var.name + ' [', var.type + ', ' + str(n_dim) + ', ' + str(p_type) + ']'
      continue

    var.f_side = f_side_trans[var.type, n_dim, p_type]
    var.c_side = c_side_trans[var.type, n_dim, p_type]

    # If allocatable or pointer type then add the dimensional varaibles to the struct.dim_var list.

    if len(var.array) != 0 and var.array[0] == ':':
      for n in range(1, len(var.array)):
        dim_var = var_class()
        dim_var.name = 'n' + str(n) + '_' + var.name
        dim_var.type = 'integer'
        dim_var.f_side = f_side_trans['integer', 0, NOT]
        struct.dim_var.append(dim_var)

    if len(var.array) > 0 and p_type == NOT:
      dim1 = str(1 + int(var.ubound[0]) - int(var.lbound[0]))
      var.c_side.constructor = var.c_side.constructor.replace('DIM1', dim1)
      var.c_side.to_c2_set   = var.c_side.to_c2_set.replace('DIM1', dim1)
      var.c_side.to_f_setup  = var.c_side.to_f_setup.replace('DIM1', dim1)
      var.f_side.to_f2_trans = var.f_side.to_f2_trans.replace('DIM1', dim1)
      var.f_side.test_pat    = var.f_side.test_pat.replace('DIM1', dim1)

    if len(var.array) > 1 and p_type == NOT:
      dim2 = str(1 + int(var.ubound[1]) - int(var.lbound[1]))
      var.c_side.constructor = var.c_side.constructor.replace('DIM2', dim2)
      var.c_side.to_c2_set   = var.c_side.to_c2_set.replace('DIM2', dim2)
      var.c_side.to_f_setup  = var.c_side.to_f_setup.replace('DIM2', dim2)
      var.f_side.to_f2_trans = var.f_side.to_f2_trans.replace('DIM2', dim2)
      var.f_side.test_pat    = var.f_side.test_pat.replace('DIM2', dim2)


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

f_face.write ('\n'.join(f_use_list))

f_face.write ('''
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
    to_c2_arg_def[f_side.bindc_type].append(f_side.bindc_name.replace('NAME', var.name))

  for dim_var in struct.dim_var: 
    f_side = dim_var.f_side
    import_set.add(f_side.bindc_const)
    if not f_side.bindc_type in to_c2_arg_def: to_c2_arg_def[f_side.bindc_type] = []
    to_c2_arg_def[f_side.bindc_type].append(f_side.bindc_name.replace('NAME', var.name))

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
    f_face.write (', ' + var.f_side.to_c2_arg.replace('NAME', var.name))

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
    f_face.write(', ' + var.name.replace('NAME', var.name))

  f_face.write(''') bind(c)\n

implicit none

type (c_ptr), value :: FF
type (zzz_struct), pointer :: FP
'''.replace('zzz', struct.short_name))

  f2_arg_list = {}
  for var in struct.var:
    if not var.f_side.bindc_type in f2_arg_list: f2_arg_list[var.f_side.bindc_type] = []
    f2_arg_list[var.f_side.bindc_type].append(var.f_side.bindc_name.replace('NAME', var.name))

  for arg_type, arg_list in f2_arg_list.items():
    f_face.write(arg_type + ' ' + ', '.join(arg_list) + '\n')

  f_face.write('''
call c_f_pointer (FF, FP)
''')

  for var in struct.var:
    f_face.write (var.f_side.to_f2_trans.replace('NAME', var.name) + '\n')

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

f_equ.write ('module bmad_equality\n\n')
f_equ.write ('\n'.join(f_use_list))

f_equ.write ('''

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
    f_equ.write ('is_eq = is_eq .and. ' + var.f_side.equal_test.replace('NAME', var.name) + '\n')

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
  print *, 'zzz: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_zzz

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_zzz (c_zzz, c_ok) bind(c)

implicit  none

type (c_ptr), value ::  c_zzz
type (zzz_struct), target :: f_zzz, f2_zzz
logical(c_bool) c_ok

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
  f_test.write (var.f_side.test_pat.replace('XXX', str(i)).replace('NAME', var.name) + '\n')

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

typedef const bool*              BoolArr;
typedef const dcomplex*          DcomplexArr;
typedef const double*            RealArr;
typedef const int*               IntArr;

typedef valarray<bool>           Bool_Array;
typedef valarray<dcomplex>       Dcomplex_Array;
typedef valarray<double>         Real_Array;
typedef valarray<int>            Int_Array;

typedef valarray<Bool_Array>     Bool_Matrix;
typedef valarray<Dcomplex_Array> Dcomplex_Matrix;
typedef valarray<Real_Array>     Real_Matrix;
typedef valarray<Int_Array>      Int_Matrix;

typedef valarray<Real_Matrix>      Real_Tensor;
typedef valarray<Dcomplex_Matrix>  Dcomplex_Tensor;

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
    construct_list.append(var.c_side.constructor.replace('NAME', var.name))

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

//---------------------------------------------------------------------------
// Instantiate instances for conversion from array to C++ structure.

template void operator<< (Bool_Array&,  const bool*);
template void operator<< (Bool_Matrix&, const bool*);

template void operator<< (Real_Array&,  const double*);
template void operator<< (Real_Matrix&, const double*);
template void operator<< (Real_Tensor&, const double*);

template void operator<< (Dcomplex_Array&,  const dcomplex*);
template void operator<< (Dcomplex_Matrix&, const dcomplex*);
template void operator<< (Dcomplex_Tensor&, const dcomplex*);

template void operator<< (Int_Array&,   const int*);
template void operator<< (Int_Matrix&,  const int*);

//---------------------------------------------------------------------------
// Instantiate instances for .

template void operator<< (Real_Array&,  const Real_Array&);
template void operator<< (Real_Matrix&, const Real_Matrix&);
template void operator<< (Real_Tensor&, const Real_Tensor&);

template void operator<< (Dcomplex_Array&,  const Dcomplex_Array&);
template void operator<< (Dcomplex_Matrix&, const Dcomplex_Matrix&);
template void operator<< (Dcomplex_Tensor&, const Dcomplex_Tensor&);

template void operator<< (Int_Array&,   const Int_Array&);
template void operator<< (Int_Matrix&,  const Int_Matrix&);

template void matrix_to_vec (const Bool_Matrix&,     bool*);
template void matrix_to_vec (const Dcomplex_Matrix&, dcomplex*);
template void matrix_to_vec (const Real_Matrix&,     double*);
template void matrix_to_vec (const Int_Matrix&,      int*);

//---------------------------------------------------------------------------

void tensor_to_vec (const Real_Tensor& tensor, double* vec) {
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

  f_cpp.write(');\n\n')

  # zzz_to_f

  f_cpp.write('extern "C" void zzz_to_f (C_zzz& C, zzz_struct* F) {\n'.replace('zzz', struct.short_name))

  for var in struct.var:
    f_cpp.write (var.c_side.to_f_setup.replace('NAME', var.name))

  f_cpp.write('  zzz_to_f2 (F'.replace('zzz', struct.short_name))

  for var in struct.var:
    f_cpp.write (', ' + var.c_side.to_f2_call.replace('NAME', var.name))

  f_cpp.write(');\n')
  f_cpp.write('}\n')

  # zzz_to_c2

  f_cpp.write('\n')
  f_cpp.write('extern "C" void zzz_to_c2 (C_zzz& C'.replace('zzz', struct.short_name))

  for var in struct.var:
    f_cpp.write (', ' + var.c_side.to_c2_arg.replace('NAME', var.name).replace('zzz', struct.short_name))

  f_cpp.write(') {\n')

  for var in struct.var:
    f_cpp.write ('  ' + var.c_side.to_c2_set.replace('NAME', var.name) + '\n')

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

template <class T> bool is_all_equal_vec (const valarray<T>& vec1, const valarray<T>& vec2) {
  bool is_eq = true;
  if (vec1.size() != vec2.size()) return false;
  for (int i = 0; i < vec1.size(); i++) {
    is_eq = is_eq && (vec1[i] == vec2[i]);
  }
  return is_eq;
}

template <class T> bool is_all_equal_mat (const valarray< valarray<T> >& mat1, const valarray< valarray<T> >& mat2) {
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

template <class T> bool is_all_equal_tensor (const valarray< valarray< valarray<T> > >& tensor1, const valarray< valarray< valarray<T> > >& tensor2) {
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

template bool is_all_equal_vec (const Bool_Array&,     const Bool_Array&);
template bool is_all_equal_vec (const Dcomplex_Array&, const Dcomplex_Array&);
template bool is_all_equal_vec (const Real_Array&,     const Real_Array&);
template bool is_all_equal_vec (const Int_Array&,      const Int_Array&);

template bool is_all_equal_mat (const Bool_Matrix&,     const Bool_Matrix&);
template bool is_all_equal_mat (const Dcomplex_Matrix&, const Dcomplex_Matrix&);
template bool is_all_equal_mat (const Real_Matrix&,     const Real_Matrix&);
template bool is_all_equal_mat (const Int_Matrix&,      const Int_Matrix&);

template bool is_all_equal_tensor (const Dcomplex_Tensor&, const Dcomplex_Tensor&);
template bool is_all_equal_tensor (const Real_Tensor&,     const Real_Tensor&);

''')

for key, struct in struct_def.items():
  f_eq.write ('\n//--------------------------------------------------------------\n\n')
  f_eq.write ('bool operator== (const C_zzz& x, const C_zzz& y) {'.replace('zzz', struct.short_name) + '\n')
  f_eq.write ('  bool is_eq = true;\n')

  for var in struct.var:
    f_eq.write ('  is_eq = is_eq && ' + var.c_side.equal_test.replace('NAME', var.name)  + ';\n')

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
  f_test.write (var.c_side.test_pat.replace('XXX', str(i)).replace('NAME', var.name) + '\n')

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
