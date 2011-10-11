#+
# Script to check that all library symbols are unique.
# Script checks the libraries:
#   vms libraries
#   $ACC_RELEASE_DIR/lib/mpm_utils
#   $ACC_RELEASE_DIR/lib/cesr_utils
#   $ACC_RELEASE_DIR/lib/sim_utils
#
# How to run:
#   1) Make a listing of the stulib routines in [cesr.axlib]stulib.olb:
#         lib/list stulib
#
#   2) Store this list along with listings from any other VMS library in the linux file: 
#         util_programs/scripts/vms_lib.symbols
# 
#   3) Run python program in util_programs/scripts (must be in this directory)
#         python beam_lib_routine_check
#
# Note: to get a merged list of symbols, use the -dump argument on the command line:
#         python beam_lib_routine_check -dump
#-

import subprocess
import sys

#

symbols = dict()
duplicate = False
dump = False

if len(sys.argv) > 1:
  if sys.argv[1] == '-dump':
    dump = True
  else:
    print 'ERROR: UNKNOWN ARGUMENT: ' + sys.argv[1]

# stulib

stu_file = open('vms_lib.symbols', 'r')
for line in stu_file.readlines():
  line = line.rstrip() 
  if line == '': continue
  symbols[line.lower()] = 'stulib'

# mpm_utils

process = subprocess.Popen(['ar t $ACC_RELEASE_DIR/lib/libmpm_utils.a'], shell = True, stdout = subprocess.PIPE)
mpm_list = process.communicate()[0].replace('_DBL.o', '').split()

for sym in mpm_list:
  if symbols.has_key(sym):
    print 'Duplicate symbol: ' + sym 
    print '   In: mpm_utils, ' + symbols[sym]
    duplicate = True
  else:
    symbols[sym] = 'mpm_utils'
  
# cesr_utils

process = subprocess.Popen(['ar t $ACC_RELEASE_DIR/lib/libcesr_utils.a'], shell = True, stdout = subprocess.PIPE)
cesr_list = process.communicate()[0].replace('_DBL.o', '').split()

for sym in cesr_list:
  if symbols.has_key(sym):
    print 'Duplicate symbol: ' + sym 
    print '   In: cesr_utils, ' + symbols[sym]
    duplicate = True
  else:
    symbols[sym] = 'cesr_utils'
  
# sim_utils

process = subprocess.Popen(['ar t $ACC_RELEASE_DIR/lib/libsim_utils.a'], shell = True, stdout = subprocess.PIPE)
sim_list = process.communicate()[0].replace('_DBL.o', '').split()

for sym in sim_list:
  if symbols.has_key(sym):
    print 'Duplicate symbol: ' + sym 
    print '   In: sim_utils, ' + symbols[sym]
    duplicate = True
  else:
    symbols[sym] = 'sim_utils'
  
#

if not duplicate: 
  print 'Everything OK!'

if dump:
  for key in sorted(symbols.keys()):
    print '%-30s  %s' % (key, symbols[key])
