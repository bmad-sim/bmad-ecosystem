#!/usr/bin/env python

# Note: Run this script in the cpp_bmad_interface directory.

import re
import os

def searchit (file):

  re_p = re.compile('INTEGER, *PARAMETER *:: *')
  re_a = re.compile('\[')

  params_here = False

  f_in = open(file)
  for line in f_in:
    line = line.partition('!')[0]   # Strip off comment
    line = line.upper()
    if not re_p.match(line) and not params_here: continue
    if '[' in line: continue                              # Skip parameter arrays

    line = re_p.sub('const int ', line)
    line = line.replace('$', '')

    if '&' in line:
      line = line.replace('&', '')
      params_here = True
      line = '  ' + line.rstrip() + '\n'
    else:
      params_here = False
      line = '  ' + line.rstrip() + ';\n'

    f_out.write(line)

#---------------------------------------

if not os.path.exists('include'): os.makedirs('include')
f_out = open('include/bmad_enums.h', 'w')

f_out.write('''
#ifndef BMAD_ENUMS

namespace Bmad {
''')

searchit('../bmad/modules/bmad_struct.f90')
searchit('../bmad/modules/basic_bmad_mod.f90')
searchit('../sim_utils/io/output_mod.f90')

f_out.write('''
}

#define BMAD_ENUMS
#endif
''')

print 'Created: include/bmad_enums.h'

