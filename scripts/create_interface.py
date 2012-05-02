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

debug = False

def print_debug (line):
  if (debug): print line

##################################################################################
##################################################################################
# struct_def_class
# Class for a structure

class struct_def_class:
  def __init__(self, f_name = ''):
    self.f_name = f_name # Struct name on Fortran side
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

# Fortran side translation

class f_side_trans_class:

  def __init__(self, to_c2_arg, bindc_type, bindc_name):
    self.to_c2_arg = to_c2_arg
    self.bindc_type = bindc_type
    self.bindc_name = bindc_name
    self.bindc_const = ''
    self.equal_test = '(f1%NAME == f2%NAME)'
    self.to_f2_trans = 'F%NAME = z_NAME'
    self.test_pat = 'F%NAME = XXX + offset'
    self.to_c_var = ''
    self.to_c_trans = ''

  def __repr__(self):
    return '%s,  %s,  %s :: %s' % (self.to_c2_arg, self.bindc_type, self.bindc_name, self.to_f2_trans)


#       Dim  P_type                    to_c2_arg                bindc_type                    bindc_name   equal_test
f_side_trans = {
  (REAL,  0, NOT) : f_side_trans_class('F%NAME',               'real(c_double)',             'z_NAME'),
  (REAL,  1, NOT) : f_side_trans_class('F%NAME',               'real(c_double)',             'z_NAME(*)'),
  (REAL,  2, NOT) : f_side_trans_class('mat2vec(F%NAME)',      'real(c_double)',             'z_NAME(*)'),
  (REAL,  3, NOT) : f_side_trans_class('tensor2vec(F%NAME)',   'real(c_double)',             'z_NAME(*)'),
  (CMPLX, 0, NOT) : f_side_trans_class('F%NAME',               'complex(c_double_complex)',  'z_NAME'),
  (CMPLX, 1, NOT) : f_side_trans_class('F%NAME',               'complex(c_double_complex)',  'z_NAME(*)'),
  (CMPLX, 2, NOT) : f_side_trans_class('mat2vec(F%NAME)',      'complex(c_double_complex)',  'z_NAME(*)'),
  (CMPLX, 3, NOT) : f_side_trans_class('tensor2vec(F%NAME)',   'complex(c_double_complex)',  'z_NAME(*)'),
  (INT,   0, NOT) : f_side_trans_class('F%NAME',               'integer(c_int)',             'z_NAME'),
  (INT,   1, NOT) : f_side_trans_class('F%NAME',               'integer(c_int)',             'z_NAME(*)'),
  (INT,   2, NOT) : f_side_trans_class('mat2vec(F%NAME)',      'integer(c_int)',             'z_NAME(*)'),
  (LOGIC, 0, NOT) : f_side_trans_class('c_logic(F%NAME)',      'logical(c_bool)',            'z_NAME'),
  (LOGIC, 1, NOT) : f_side_trans_class('c_logic(F%NAME)',      'logical(c_bool)',            'z_NAME(*)'),
  (LOGIC, 2, NOT) : f_side_trans_class('mat2vec(F%NAME)',      'logical(c_bool)',            'z_NAME(*)'),
  (TYPE,  0, NOT) : f_side_trans_class('c_loc(F%NAME)',        'type(c_ptr), value ::',      'z_NAME'),
  (TYPE,  1, NOT) : f_side_trans_class('z_NAME',               'type(c_ptr) ::',             'z_NAME(*)'), 
  (TYPE,  2, NOT) : f_side_trans_class('z_NAME',               'type(c_ptr) ::',             'z_NAME(*)'), 
  (TYPE,  3, NOT) : f_side_trans_class('z_NAME',               'type(c_ptr) ::',             'z_NAME(*)'), 
  (CHAR,  0, NOT) : f_side_trans_class('trim(F%NAME) // c_null_char', 'character(c_char)',   'z_NAME(*)')
  }

jd1_loop = 'do jd1 = 1, size(F%NAME,1); lb1 = lbound(F%NAME,1) - 1\n'
jd2_loop = 'do jd2 = 1, size(F%NAME,2); lb2 = lbound(F%NAME,2) - 1\n'
jd3_loop = 'do jd3 = 1, size(F%NAME,3); lb3 = lbound(F%NAME,3) - 1\n'

rhs1 = '  rhs = 100 + jd1 + XXX + offset\n'
rhs2 = '  rhs = 100 + jd1 + 10*jd2 + XXX + offset\n'
rhs3 = '  rhs = 100 + jd1 + 10*jd2 + 100*jd3 + XXX + offset\n'

test_pat1 = jd1_loop + rhs1 + '  F%NAME(jd1+lb1) = NNN\n' + 'enddo\n'
test_pat2 = jd1_loop + jd2_loop + rhs2 + '  F%NAME(jd1+lb1,jd2+lb2) = NNN\n' + 'enddo; enddo\n'
test_pat3 = jd1_loop + jd2_loop + jd3_loop + rhs3 + '  F%NAME(jd1+lb1,jd2+lb2,jd3+lb3) = NNN\n' + 'enddo; enddo; enddo\n'

f_side_trans[CHAR,  0, NOT].test_pat    = \
        'do jd1 = 1, len(F%NAME)\n  F%NAME(jd1:jd1) = char(ichar("a") + modulo(100+XXX+offset+jd1, 26))\nenddo'
f_side_trans[CHAR,  0, NOT].to_f2_trans = 'call to_f_str(z_NAME, F%NAME)'

f_side_trans[LOGIC, 0, NOT].test_pat    = 'F%NAME = (modulo(XXX + offset, 2) == 0)'

f_side_trans[TYPE,  0, NOT].to_f2_trans = 'call KIND_to_f(z_NAME, c_loc(F%NAME))'
f_side_trans[TYPE,  0, NOT].test_pat    = 'call KIND_test_pattern (F%NAME, ix_patt)'

f_side_trans[TYPE,  1, NOT].to_f2_trans = jd1_loop + \
                              '  call KIND_to_f(z_NAME(jd1), c_loc(F%NAME(jd1+lb1)))\nenddo'
f_side_trans[TYPE,  1, NOT].test_pat    = \
                        jd1_loop + rhs1 + '  call KIND_test_pattern (F%NAME(jd1+lb1), ix_patt+jd1)\n' + 'enddo\n'
f_side_trans[TYPE,  1, NOT].to_c_var    = 'type(c_ptr) :: z_NAME(DIM1)'
f_side_trans[TYPE,  1, NOT].to_c_trans  = jd1_loop + '  z_NAME(jd1) = c_loc(F%NAME(jd1+lb1))\nenddo\n\n'

f_side_trans[REAL,  1, NOT].to_f2_trans = 'F%NAME = z_NAME(1:size(F%NAME))'
f_side_trans[REAL,  1, NOT].test_pat    = test_pat1.replace('NNN', 'rhs')

f_side_trans[INT,   1, NOT].to_f2_trans = 'F%NAME = z_NAME(1:size(F%NAME))'
f_side_trans[INT,   1, NOT].test_pat    = test_pat1.replace('NNN', 'rhs')

f_side_trans[LOGIC, 1, NOT].to_f2_trans = 'F%NAME = f_logic(z_NAME(1:size(F%NAME)))'
f_side_trans[LOGIC, 1, NOT].test_pat    = test_pat1.replace('NNN', 'modulo(rhs, 2) == 0')

f_side_trans[CMPLX, 1, NOT].to_f2_trans = 'F%NAME = z_NAME(1:size(F%NAME))'
f_side_trans[CMPLX, 1, NOT].test_pat    = test_pat1.replace('NNN', 'cmplx(rhs, 100 + rhs)')

f_side_trans[TYPE,  2, NOT].to_f2_trans = jd1_loop + jd2_loop + \
                        '  call KIND_to_f(z_NAME(DIM2*(jd1-1) + jd2), c_loc(F%NAME(jd1+lb1,jd2+lb2)))\n' + 'enddo; enddo\n'
f_side_trans[TYPE,  2, NOT].test_pat    = jd1_loop + jd2_loop + rhs2 + \
                    '  call KIND_test_pattern (F%NAME(jd1+lb1,jd2+lb2), ix_patt+jd1+10*jd2)\n' + 'enddo; enddo\n'
f_side_trans[TYPE,  2, NOT].to_c_var    = 'type(c_ptr) :: z_NAME(DIM1*DIM2)'
f_side_trans[TYPE,  2, NOT].to_c_trans  = jd1_loop + jd2_loop + \
                    '  z_NAME(DIM2*(jd1-1) + jd2) = c_loc(F%NAME(jd1+lb1,jd2+lb2))\nenddo; enddo\n\n'

f_side_trans[REAL,  2, NOT].to_f2_trans = 'F%NAME = vec2mat(z_NAME, size(F%NAME, 1), size(F%NAME, 2))'
f_side_trans[REAL,  2, NOT].test_pat    = test_pat2.replace('NNN', 'rhs')

f_side_trans[INT,   2, NOT].to_f2_trans = 'F%NAME = vec2mat(z_NAME, size(F%NAME, 1), size(F%NAME, 2))'
f_side_trans[INT,   2, NOT].test_pat    = test_pat2.replace('NNN', 'rhs')

f_side_trans[LOGIC, 2, NOT].to_f2_trans = 'F%NAME = vec2mat(z_NAME, size(F%NAME, 1), size(F%NAME, 2))'
f_side_trans[LOGIC, 2, NOT].test_pat    = test_pat2.replace('NNN', '(modulo(rhs, 2) == 0)')

f_side_trans[CMPLX, 2, NOT].to_f2_trans = 'F%NAME = vec2mat(z_NAME, size(F%NAME, 1), size(F%NAME, 2))'
f_side_trans[CMPLX, 2, NOT].test_pat    = test_pat2.replace('NNN', 'cmplx(rhs, 100+rhs)')

f_side_trans[REAL,  3, NOT].to_f2_trans = 'F%NAME = vec2tensor(z_NAME, size(F%NAME, 1), size(F%NAME, 2), size(F%NAME, 3))'
f_side_trans[REAL,  3, NOT].test_pat    = test_pat3.replace('NNN', 'rhs')

f_side_trans[CMPLX, 3, NOT].to_f2_trans = 'F%NAME = vec2tensor(z_NAME, size(F%NAME, 1), size(F%NAME, 2), size(F%NAME, 3))'
f_side_trans[CMPLX, 3, NOT].test_pat    = test_pat3.replace('NNN', 'cmplx(rhs, 100+rhs)')

f_side_trans[TYPE,  3, NOT].to_f2_trans = jd1_loop + jd2_loop + jd3_loop + \
                   '  call KIND_to_f(z_NAME(DIM3*DIM2*(jd1-1) + DIM3*(jd2-1) + jd3), c_loc(F%NAME(jd1+lb1,jd2+lb2,jd3+lb3)))\n' + \
                   'enddo; enddo; enddo\n'
f_side_trans[TYPE,  3, NOT].test_pat    = jd1_loop + jd2_loop + jd3_loop + rhs3 + \
                    '  call KIND_test_pattern (F%NAME(jd1+lb1,jd2+lb2,jd3+lb3), ix_patt+jd1+10*jd2+100*jd3)\n' + 'enddo; enddo; enddo\n'
f_side_trans[TYPE,  3, NOT].to_c_var    = 'type(c_ptr) :: z_NAME(DIM1*DIM2*DIM3)'
f_side_trans[TYPE,  3, NOT].to_c_trans  = jd1_loop + jd2_loop + jd3_loop + \
                    '  z_NAME(DIM3*DIM2*(jd1-1) + DIM3*(jd2-1) + jd3) = c_loc(F%NAME(jd1+lb1,jd2+lb2,jd3+lb3))\n' + 'enddo; enddo; enddo\n\n'

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
    self.to_c2_set = 'C.NAME = z_NAME;'
    self.constructor = 'NAME(0)'
    self.equal_test = '(x.NAME == y.NAME)'
    self.test_pat = 'C.NAME = XXX + offset;'

  def __repr__(self):
    return '%s,  %s,  %s,  %s' % (self.c_class, self.to_f2_arg, self.to_f2_call, self.to_c2_arg)


#        Dim  P_type                    c_class            to_f2_arg          to_f2_call         to_c2_arg
c_side_trans = {
  (REAL,   0, NOT) : c_side_trans_class('double',          'Real&',           'C.NAME',          'Real& z_NAME'),
  (REAL,   1, NOT) : c_side_trans_class('Real_Array',      'RealArr',         '&C.NAME[0]',      'RealArr z_NAME'),
  (REAL,   2, NOT) : c_side_trans_class('Real_Matrix',     'RealArr',         'z_NAME',          'RealArr z_NAME'),
  (REAL,   3, NOT) : c_side_trans_class('Real_Tensor',     'RealArr',         'z_NAME',          'RealArr z_NAME'),
  (CMPLX,  0, NOT) : c_side_trans_class('dcomplex',        'Dcomplex&',       'C.NAME',          'Dcomplex& z_NAME'),
  (CMPLX,  1, NOT) : c_side_trans_class('Dcomplex_Array',  'DcomplexArr',     '&C.NAME[0]',      'DcomplexArr z_NAME'),
  (CMPLX,  2, NOT) : c_side_trans_class('Dcomplex_Matrix', 'DcomplexArr',     'z_NAME',          'DcomplexArr z_NAME'),
  (CMPLX,  3, NOT) : c_side_trans_class('Dcomplex_Tensor', 'DcomplexArr',     'z_NAME',          'DcomplexArr z_NAME'),
  (INT,    0, NOT) : c_side_trans_class('int',             'Int&',            'C.NAME',          'Int& z_NAME'),
  (INT,    1, NOT) : c_side_trans_class('Int_Array',       'IntArr',          '&C.NAME[0]',      'IntArr z_NAME'),
  (INT,    2, NOT) : c_side_trans_class('Int_Matrix',      'IntArr',          'z_NAME',          'IntArr z_NAME'),
  (LOGIC,  0, NOT) : c_side_trans_class('bool',            'Bool&',           'C.NAME',          'Bool& z_NAME'),
  (LOGIC,  1, NOT) : c_side_trans_class('Bool_Array',      'BoolArr',         '&C.NAME[0]',      'BoolArr z_NAME'),
  (LOGIC,  2, NOT) : c_side_trans_class('Bool_Matrix',     'BoolArr',         'z_NAME',          'BoolArr z_NAME'),
  (TYPE,   0, NOT) : c_side_trans_class('C_KIND',          'const C_KIND&',   'C.NAME',          'const KIND_struct* z_NAME'),
  (TYPE,   1, NOT) : c_side_trans_class('C_KIND_Array',    'const C_KIND**',  'z_NAME',          'const KIND_struct** z_NAME'),
  (TYPE,   2, NOT) : c_side_trans_class('C_KIND_Matrix',   'const C_KIND**',  'z_NAME',          'const KIND_struct** z_NAME'),
  (TYPE,   3, NOT) : c_side_trans_class('C_KIND_Tensor',   'const C_KIND**',  'z_NAME',          'const KIND_struct** z_NAME'),
  (CHAR,   0, NOT) : c_side_trans_class('string',          'Char',            'C.NAME.c_str()',  'Char z_NAME')
  }

for1 = '  for (int i = 0; i < C.NAME.size(); i++)'
for2 = ' for (int j = 0; j < C.NAME[0].size(); j++) ' 
for3 = ' for (int k = 0; k < C.NAME[0][0].size(); k++)'

test_pat1 = for1 + '\n    {int rhs = 101 + i + XXX + offset; C.NAME[i] = NNN;}'
test_pat2 = for1 + for2 + '\n    {int rhs = 101 + i + 10*(j+1) + XXX + offset; C.NAME[i][j] = NNN;}'
test_pat3 = for1 + for2 + for3 + '\n    {int rhs = 101 + i + 10*(j+1) + 100*(k+1) + XXX + offset; C.NAME[i][j][k] = NNN;}'

c_side_trans[CHAR,  0, NOT].test_pat    = 'C.NAME.resize(STR_LEN);\n' + test_pat1.replace('NNN', "'a' + rhs % 26")
c_side_trans[CHAR,  0, NOT].constructor = "NAME()"

c_side_trans[LOGIC, 0, NOT].test_pat    = 'C.NAME = ((XXX + offset % 2) == 0);'

c_side_trans[TYPE,  0, NOT].constructor = "NAME()"
c_side_trans[TYPE,  0, NOT].to_c2_set   = 'KIND_to_c(z_NAME, C.NAME);' 
c_side_trans[TYPE,  0, NOT].test_pat    = 'C_KIND_test_pattern(C.NAME, ix_patt);'

c_side_trans[TYPE,  1, NOT].constructor = "NAME(C_KIND_Array(C_KIND(), DIM1))"
c_side_trans[TYPE,  1, NOT].to_c2_set   = for1 + ' KIND_to_c(z_NAME[i], C.NAME[i]);' 
c_side_trans[TYPE,  1, NOT].test_pat    = test_pat1.replace('C.NAME[i] = NNN', 'C_KIND_test_pattern(C.NAME[i], ix_patt+i+1);')
c_side_trans[TYPE,  1, NOT].to_f_setup  = '  const C_KIND* z_NAME[DIM1];\n' + for1 + ' z_NAME[i] = &C.NAME[i];\n'

c_side_trans[REAL,  1, NOT].constructor = 'NAME(0.0, DIM1)'
c_side_trans[REAL,  1, NOT].to_c2_set   = 'C.NAME = Real_Array(z_NAME, DIM1);'
c_side_trans[REAL,  1, NOT].test_pat    = test_pat1.replace('NNN', 'rhs')

c_side_trans[CMPLX, 1, NOT].constructor = 'NAME(0.0, DIM1)'
c_side_trans[CMPLX, 1, NOT].to_c2_set   = 'C.NAME = Dcomplex_Array(z_NAME, DIM1);'
c_side_trans[CMPLX, 1, NOT].test_pat    = test_pat1.replace('NNN', 'Dcomplex(rhs, 100+rhs)')

c_side_trans[INT,   1, NOT].constructor = 'NAME(DIM1)'
c_side_trans[INT,   1, NOT].to_c2_set   = 'C.NAME = Int_Array(z_NAME, DIM1);'
c_side_trans[INT,   1, NOT].test_pat    = test_pat1.replace('NNN', 'rhs')

c_side_trans[LOGIC, 1, NOT].constructor = 'NAME(DIM1)'
c_side_trans[LOGIC, 1, NOT].to_c2_set   = 'C.NAME = Bool_Array(z_NAME, DIM1);'
c_side_trans[LOGIC, 1, NOT].test_pat    = test_pat1.replace('NNN', '((rhs % 2) == 0)')

c_side_trans[REAL,  2, NOT].constructor = 'NAME(Real_Array(0.0, DIM2), DIM1)'
c_side_trans[REAL,  2, NOT].to_c2_set   = 'C.NAME << z_NAME;'
c_side_trans[REAL,  2, NOT].test_pat    = test_pat2.replace('NNN', 'rhs')
c_side_trans[REAL,  2, NOT].to_f_setup  = '  double z_NAME[DIM1*DIM2]; matrix_to_vec(C.NAME, z_NAME);\n'

c_side_trans[INT,   2, NOT].constructor = 'NAME(Int_Array(0, DIM2), DIM1)'
c_side_trans[INT,   2, NOT].to_c2_set   = 'C.NAME << z_NAME;'
c_side_trans[INT,   2, NOT].test_pat    = test_pat2.replace('NNN', 'rhs')
c_side_trans[INT,   2, NOT].to_f_setup  = '  int z_NAME[DIM1*DIM2]; matrix_to_vec(C.NAME, z_NAME);\n'

c_side_trans[LOGIC, 2, NOT].constructor = 'NAME(Bool_Array(false, DIM2), DIM1)'
c_side_trans[LOGIC, 2, NOT].to_c2_set   = 'C.NAME << z_NAME;'
c_side_trans[LOGIC, 2, NOT].test_pat    = test_pat2.replace('NNN', '((rhs % 2) == 0)')
c_side_trans[LOGIC, 2, NOT].to_f_setup  = '  bool z_NAME[DIM1*DIM2]; matrix_to_vec(C.NAME, z_NAME);\n'

c_side_trans[CMPLX, 2, NOT].constructor = 'NAME(Dcomplex_Array(0.0, DIM2), DIM1)'
c_side_trans[CMPLX, 2, NOT].to_c2_set   = 'C.NAME << z_NAME;'
c_side_trans[CMPLX, 2, NOT].test_pat    = test_pat2.replace('NNN', 'Dcomplex(rhs, 100+rhs)')
c_side_trans[CMPLX, 2, NOT].to_f_setup  = '  dcomplex z_NAME[DIM1*DIM2]; matrix_to_vec(C.NAME, z_NAME);\n'

c_side_trans[TYPE,  2, NOT].constructor = "NAME(C_KIND_Array(C_KIND(), DIM2), DIM1)"
c_side_trans[TYPE,  2, NOT].to_c2_set   = for1 + for2 + '\n    {int m = DIM2*i + j; KIND_to_c(z_NAME[m], C.NAME[i][j]);}' 
c_side_trans[TYPE,  2, NOT].test_pat    = test_pat2.replace('C.NAME[i][j] = NNN', 'C_KIND_test_pattern(C.NAME[i][j], ix_patt+i+1+10*(j+1))')
c_side_trans[TYPE,  2, NOT].to_f_setup  = \
    '  const C_KIND* z_NAME[DIM1*DIM2];\n' + for1 + for2 +  '\n    {int m = DIM2*i + j; z_NAME[m] = &C.NAME[i][j];}\n'

c_side_trans[REAL,  3, NOT].constructor = 'NAME(Real_Matrix(Real_Array(0.0, DIM3), DIM2), DIM1)'
c_side_trans[REAL,  3, NOT].to_c2_set   = 'C.NAME << z_NAME;'
c_side_trans[REAL,  3, NOT].test_pat    = test_pat3.replace('NNN', 'rhs')
c_side_trans[REAL,  3, NOT].to_f_setup  = '  double z_NAME[DIM1*DIM2*DIM3]; tensor_to_vec(C.NAME, z_NAME);\n'

c_side_trans[CMPLX, 3, NOT].constructor = 'NAME(Dcomplex_Matrix(Dcomplex_Array(0.0, DIM3), DIM2), DIM1)'
c_side_trans[CMPLX, 3, NOT].to_c2_set   = 'C.NAME << z_NAME;'
c_side_trans[CMPLX, 3, NOT].test_pat    = test_pat3.replace('NNN', 'Dcomplex(rhs, 100+rhs)')
c_side_trans[CMPLX, 3, NOT].to_f_setup  = '  dcomplex z_NAME[DIM1*DIM2*DIM3]; tensor_to_vec(C.NAME, z_NAME);\n'

c_side_trans[TYPE,  3, NOT].constructor = "NAME(C_KIND_Matrix(C_KIND_Array(C_KIND(), DIM3), DIM2), DIM1)"
c_side_trans[TYPE,  3, NOT].to_c2_set   = for1 + for2 + for3 + '\n    {int m = DIM3*DIM2*i + DIM3*j + k; KIND_to_c(z_NAME[m], C.NAME[i][j][k]);}' 
c_side_trans[TYPE,  3, NOT].test_pat    = test_pat3.replace('C.NAME[i][j][k] = NNN', 'C_KIND_test_pattern(C.NAME[i][j][k], ix_patt+i+1+10*(j+1)+100*(k+1))')
c_side_trans[TYPE,  3, NOT].to_f_setup  = \
    '  const C_KIND* z_NAME[DIM1*DIM2*DIM3];\n' + for1 + for2 + for3 + '\n    {int m = DIM3*DIM2*i + DIM3*j + k; z_NAME[m] = &C.NAME[i][j][k];}\n'

for key, c in c_side_trans.items(): 
  if key[1] == 1: c.equal_test = 'is_all_equal (x.NAME, y.NAME)'
  if key[1] == 2: c.equal_test = 'is_all_equal (x.NAME, y.NAME)'
  if key[1] == 3: c.equal_test = 'is_all_equal (x.NAME, y.NAME)'

##################################################################################
##################################################################################
# Get the list of structs

##  struct_list_file = 'scripts/fortran_structs.list'
struct_list_file = 'test.list'
if len(sys.argv) > 1: struct_list_file = sys.argv[1]

f_struct_list_file = open(struct_list_file)

struct_def = []
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
    struct_def.append(struct_def_class(line.split()[0]))

f_struct_list_file.close()

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
re_match2 = re.compile('([,(]|::|=)')       # Match to: ',', '::', '(', '='

for file_name in f_module_files:
  f_module_file = open('../' + file_name)

  for line in f_module_file:
    split_line = line.split()
    if len(split_line) < 2: continue
    if split_line[0] != 'type': continue

    found = False
    for struct in struct_def:
      if struct.f_name != split_line[1]: continue
      found = True
      break

    if not found: continue

    struct.short_name = struct.f_name[:-7]   # Remove '_struct' suffix

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

      if split_line[0][0] == ' ': 
        split_line = re_match2.split(split_line[1], 1)
        if split_line[0] == '': split_line.pop(0)
 
      print_debug('P3: ' + str(split_line))

      # Now split_line[0] is a delimiter or variable name
      # Add type information if there is more...

      if split_line[0] == '(':
        split_line = split_line[1].partition(')')
        base_var.kind = split_line[0].strip()
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

      # Join split_line into one string so that we are starting from a definite state.

      if len(split_line) > 1: split_line = [''.join(split_line)]
      print_debug('P5: ' + str(split_line))

      # Now len(split_line) = 1 and the first word in split_line[0] is the variable name
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
for struct in struct_def:
  if struct.short_name != '': n_found = n_found + 1
  f_out.write('******************************************\n')
  f_out.write (struct.f_name + '    ' + str(len(struct.var)) + '\n')
  for var in struct.var:
    f_out.write ('    ' + var.full_repr() + '\n')

f_out.close()


print 'Number of structs in input list: ' + str(len(struct_def))
print 'Number of structs found:         ' + str(n_found)

if len(struct_def) != n_found:
  sys.exit('COULD NOT FIND ALL THE STRUCTS! STOPPING HERE!')  

##################################################################################
##################################################################################
# Add translation info and make some name substitutions

for struct in struct_def:
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

    var.c_side.test_pat = var.c_side.test_pat.replace('STR_LEN', var.kind)

    if var.type == 'type':
      kind = var.kind[:-7]
      var.f_side.to_f2_trans = var.f_side.to_f2_trans.replace('KIND', kind)
      var.f_side.test_pat    = var.f_side.test_pat.replace('KIND', kind)
      var.c_side.test_pat    = var.c_side.test_pat.replace('KIND', kind)
      var.c_side.c_class     = var.c_side.c_class.replace('KIND', kind)
      var.c_side.to_c2_set   = var.c_side.to_c2_set.replace('KIND', kind)
      var.c_side.to_f_setup  = var.c_side.to_f_setup.replace('KIND', kind)
      var.c_side.to_f2_arg   = var.c_side.to_f2_arg.replace('KIND', kind)
      var.c_side.to_c2_arg   = var.c_side.to_c2_arg.replace('KIND', kind)
      var.c_side.constructor = var.c_side.constructor.replace('KIND', kind)

    if len(var.array) >= 1 and p_type == NOT:
      dim1 = str(1 + int(var.ubound[0]) - int(var.lbound[0]))
      var.c_side.constructor = var.c_side.constructor.replace('DIM1', dim1)
      var.c_side.to_c2_set   = var.c_side.to_c2_set.replace('DIM1', dim1)
      var.c_side.to_f_setup  = var.c_side.to_f_setup.replace('DIM1', dim1)
      var.f_side.to_f2_trans = var.f_side.to_f2_trans.replace('DIM1', dim1)
      var.f_side.test_pat    = var.f_side.test_pat.replace('DIM1', dim1)
      var.f_side.to_c_var    = var.f_side.to_c_var.replace('DIM1', dim1)
      var.f_side.to_c_trans  = var.f_side.to_c_trans.replace('DIM1', dim1)

    if len(var.array) >= 2 and p_type == NOT:
      dim2 = str(1 + int(var.ubound[1]) - int(var.lbound[1]))
      var.c_side.constructor = var.c_side.constructor.replace('DIM2', dim2)
      var.c_side.to_c2_set   = var.c_side.to_c2_set.replace('DIM2', dim2)
      var.c_side.to_f_setup  = var.c_side.to_f_setup.replace('DIM2', dim2)
      var.f_side.to_f2_trans = var.f_side.to_f2_trans.replace('DIM2', dim2)
      var.f_side.test_pat    = var.f_side.test_pat.replace('DIM2', dim2)
      var.f_side.to_c_var    = var.f_side.to_c_var.replace('DIM2', dim2)
      var.f_side.to_c_trans  = var.f_side.to_c_trans.replace('DIM2', dim2)

    if len(var.array) >= 3 and p_type == NOT:
      dim3 = str(1 + int(var.ubound[2]) - int(var.lbound[2]))
      var.c_side.constructor = var.c_side.constructor.replace('DIM3', dim3)
      var.c_side.to_c2_set   = var.c_side.to_c2_set.replace('DIM3', dim3)
      var.c_side.to_f_setup  = var.c_side.to_f_setup.replace('DIM3', dim3)
      var.f_side.to_f2_trans = var.f_side.to_f2_trans.replace('DIM3', dim3)
      var.f_side.test_pat    = var.f_side.test_pat.replace('DIM3', dim3)
      var.f_side.to_c_var    = var.f_side.to_c_var.replace('DIM3', dim3)
      var.f_side.to_c_trans  = var.f_side.to_c_trans.replace('DIM3', dim3)

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
# ZZZ_to_f interface

for struct in struct_def:
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

for struct in struct_def:

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

  f_face.write ('  subroutine ZZZ_to_c2 (C'.replace('ZZZ', s_name))
  for var in struct.var: f_face.write (', z_' + var.name)
  for var in struct.dim_var: f_face.write (', ' + var.name)
  f_face.write (') bind(c)\n')
  f_face.write ('    import ' + ', '.join(import_set) + '\n')
  f_face.write ('    type(c_ptr), value :: C\n')
  for arg_type, args in to_c2_arg_def.items():
    f_face.write ('    ' + arg_type + ' ' + ', '.join(args) + '\n')

  f_face.write ('''  end subroutine
end interface

type(c_ptr), value :: Fp
type(c_ptr), value :: C
type(ZZZ_struct), pointer :: F
integer jd1, jd2, jd3, lb1, lb2, lb3
'''.replace('ZZZ', s_name))

  for var in struct.var:
    if var.f_side.to_c_var != '': f_face.write (var.f_side.to_c_var.replace('NAME', var.name) + '\n')

  f_face.write (
'''
!

call c_f_pointer (Fp, F)
''')

  for var in struct.var:
    f_face.write (var.f_side.to_c_trans.replace('NAME', var.name))

  f_face.write ('call ZZZ_to_c2 (C'.replace('ZZZ', s_name))

  for var in struct.var:
    f_face.write (', ' + var.f_side.to_c2_arg.replace('NAME', var.name))

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

  for var in struct.var:
    f_face.write(', z_' + var.name.replace('NAME', var.name))

  f_face.write(''') bind(c)\n

implicit none

type(c_ptr), value :: Fp
type(ZZZ_struct), pointer :: F
integer jd1, jd2, jd3, lb1, lb2, lb3
'''.replace('ZZZ', struct.short_name))

  f2_arg_list = {}
  for var in struct.var:
    if not var.f_side.bindc_type in f2_arg_list: f2_arg_list[var.f_side.bindc_type] = []
    f2_arg_list[var.f_side.bindc_type].append(var.f_side.bindc_name.replace('NAME', var.name))

  for arg_type, arg_list in f2_arg_list.items():
    f_face.write(arg_type + ' ' + ', '.join(arg_list) + '\n')

  f_face.write('''
call c_f_pointer (Fp, F)
''')

  for var in struct.var:
    f_face.write (var.f_side.to_f2_trans.replace('NAME', var.name) + '\n')

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
f_equ.write ('\n'.join(f_use_list))

f_equ.write ('''

interface operator (==)
''')

for i in range(0, len(struct_def), 5):
  f_equ.write ('  module procedure ' + ', '.join('eq_' + 
                              f.short_name for f in struct_def[i:i+5]) + '\n')

f_equ.write ('''end interface

contains
''')

for struct in struct_def:
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

  for var in struct.var:
    f_equ.write ('is_eq = is_eq .and. ' + var.f_side.equal_test.replace('NAME', var.name) + '\n')

  f_equ.write ('''
end function eq_ZZZ
'''.replace('ZZZ', struct.short_name))
  
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

for struct in struct_def:
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

for struct in struct_def:
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
call ZZZ_test_pattern (f2_ZZZ, 1)

call test_c_ZZZ(c_loc(f2_ZZZ), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call ZZZ_test_pattern (f_ZZZ, 4)
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

call ZZZ_test_pattern (f2_ZZZ, 2)
if (f_ZZZ == f2_ZZZ) then
  print *, 'ZZZ: F side convert C->F: Good'
else
  print *, 'ZZZ: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call ZZZ_test_pattern (f2_ZZZ, 3)
call ZZZ_to_c (c_loc(f2_ZZZ), c_ZZZ)

end subroutine test2_f_ZZZ

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine ZZZ_test_pattern (F, ix_patt)

implicit none

type(ZZZ_struct) F
integer ix_patt, offset, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

'''.replace('ZZZ', struct.short_name))

  for i, var in enumerate(struct.var, 1):
    f_test.write (var.f_side.test_pat.replace('XXX', str(i)).replace('NAME', var.name) + '\n')

  f_test.write('''
end subroutine ZZZ_test_pattern
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

for struct in struct_def:
  f_class.write('''
//--------------------------------------------------------------------
// C_ZZZ

class ZZZ_struct {};  // Opaque class for pointers to corresponding fortran structs.

class C_ZZZ {
public:
'''.replace('ZZZ', struct.short_name))

  for var in struct.var:
    f_class.write('  ' + var.c_side.c_class.replace('ZZZ', struct.short_name) + ' ' + var.name  + ';\n')

  f_class.write ('''
  C_ZZZ() :
'''.replace('ZZZ', struct.short_name))

  construct_list = []
  for var in struct.var:
    construct_list.append(var.c_side.constructor.replace('NAME', var.name))

  f_class.write ('    ' + ',\n'.join(construct_list) + '\n')

  f_class.write('''    {}

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
// Instantiate instances for transfer

template void operator<< (Real_Array&,  const Real_Array&);
template void operator<< (Real_Matrix&, const Real_Matrix&);
template void operator<< (Real_Tensor&, const Real_Tensor&);

template void operator<< (Dcomplex_Array&,  const Dcomplex_Array&);
template void operator<< (Dcomplex_Matrix&, const Dcomplex_Matrix&);
template void operator<< (Dcomplex_Tensor&, const Dcomplex_Tensor&);

template void operator<< (Int_Array&,   const Int_Array&);
template void operator<< (Int_Matrix&,  const Int_Matrix&);

//---------------------------------------------------------------------------

template void matrix_to_vec (const Bool_Matrix&,     bool*);
template void matrix_to_vec (const Dcomplex_Matrix&, dcomplex*);
template void matrix_to_vec (const Real_Matrix&,     double*);
template void matrix_to_vec (const Int_Matrix&,      int*);

template void tensor_to_vec (const Dcomplex_Tensor&, dcomplex*);
template void tensor_to_vec (const Real_Tensor&,     double*);

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

for struct in struct_def:

  # ZZZ_to_f2
  f_cpp.write('''
//--------------------------------------------------------------------
//--------------------------------------------------------------------
// C_ZZZ

extern "C" void ZZZ_to_c (const ZZZ_struct*, C_ZZZ&);

extern "C" void ZZZ_to_f2 (ZZZ_struct*'''.replace('ZZZ', struct.short_name))

  for var in struct.var:
    f_cpp.write (', ' + var.c_side.to_f2_arg.replace('ZZZ', struct.short_name))

  f_cpp.write(');\n\n')

  # ZZZ_to_f

  f_cpp.write('extern "C" void ZZZ_to_f (const C_ZZZ& C, ZZZ_struct* F) {\n'.replace('ZZZ', struct.short_name))

  for var in struct.var:
    f_cpp.write (var.c_side.to_f_setup.replace('NAME', var.name))

  f_cpp.write('  ZZZ_to_f2 (F'.replace('ZZZ', struct.short_name))

  for var in struct.var:
    f_cpp.write (', ' + var.c_side.to_f2_call.replace('NAME', var.name))

  f_cpp.write(');\n')
  f_cpp.write('}\n')

  # ZZZ_to_c2

  f_cpp.write('\n')
  f_cpp.write('extern "C" void ZZZ_to_c2 (C_ZZZ& C'.replace('ZZZ', struct.short_name))

  for var in struct.var:
    f_cpp.write (', ' + var.c_side.to_c2_arg.replace('NAME', var.name).replace('ZZZ', struct.short_name))

  f_cpp.write(') {\n')

  for var in struct.var:
    f_cpp.write (var.c_side.to_c2_set.replace('NAME', var.name) + '\n')

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
template bool is_all_equal (const Dcomplex_Array&, const Dcomplex_Array&);
template bool is_all_equal (const Real_Array&,     const Real_Array&);
template bool is_all_equal (const Int_Array&,      const Int_Array&);

template bool is_all_equal (const Bool_Matrix&,     const Bool_Matrix&);
template bool is_all_equal (const Dcomplex_Matrix&, const Dcomplex_Matrix&);
template bool is_all_equal (const Real_Matrix&,     const Real_Matrix&);
template bool is_all_equal (const Int_Matrix&,      const Int_Matrix&);

template bool is_all_equal (const Dcomplex_Tensor&, const Dcomplex_Tensor&);
template bool is_all_equal (const Real_Tensor&,     const Real_Tensor&);

''')

for struct in struct_def:
  f_eq.write ('\n//--------------------------------------------------------------\n\n')
  f_eq.write ('bool operator== (const C_ZZZ& x, const C_ZZZ& y) {'.replace('ZZZ', struct.short_name) + '\n')
  f_eq.write ('  bool is_eq = true;\n')

  for var in struct.var:
    f_eq.write ('  is_eq = is_eq && ' + var.c_side.equal_test.replace('NAME', var.name)  + ';\n')

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

for struct in struct_def:
  f_test.write ('''
//--------------------------------------------------------------
//--------------------------------------------------------------

extern "C" void test2_f_ZZZ (C_ZZZ&, bool&);

void C_ZZZ_test_pattern (C_ZZZ& C, int ix_patt) {

int offset = 100 * ix_patt;

'''.replace('ZZZ', struct.short_name))

  for i, var in enumerate(struct.var, 1):
    f_test.write (var.c_side.test_pat.replace('XXX', str(i)).replace('NAME', var.name) + '\n')

  f_test.write('''
}

//--------------------------------------------------------------

extern "C" void test_c_ZZZ (ZZZ_struct* F, bool& c_ok) {

  C_ZZZ C, C2;

  c_ok = true;

  ZZZ_to_c (F, C);
  C_ZZZ_test_pattern (C2, 1);

  if (C == C2) {
    cout << " ZZZ: C side convert F->C: Good" << endl;
  } else {
    cout << " ZZZ: C SIDE CONVERT C->F: FAILED!" << endl;
    c_ok = false;
  }

  C_ZZZ_test_pattern (C2, 2);
  bool c_ok2;
  test2_f_ZZZ (C2, c_ok2);
  if (!c_ok2) c_ok = false;

  C_ZZZ_test_pattern (C, 3);
  if (C == C2) {
    cout << " ZZZ: F side convert F->C: Good" << endl;
  } else {
    cout << " ZZZ: F SIDE CONVERT C->F: FAILED!" << endl;
    c_ok = false;
  }

  C_ZZZ_test_pattern (C2, 4);
  ZZZ_to_f (C2, F);

}
'''.replace('ZZZ', struct.short_name))

f_test.close()
