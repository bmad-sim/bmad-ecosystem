#+
# Script to check that all library symbols are unique.
# Script checks the libraries:
# On the VMS side:
#   stulib
#   vmpmlib
#   grlib
# On the linux side:
#   $ACC_RELEASE_DIR/lib/mpm_utils
#   $ACC_RELEASE_DIR/lib/cesr_utils
#   $ACC_RELEASE_DIR/lib/sim_utils
#
# How to run:
#   1) Run the script
#         get_vms_lib_symbols
#   This will run a script on the VMS side and pull the two files that are generated over to the Linux side: 
#         vms_lib.symbols     
#         vms_lib_deb.symbols
# The second file has the symbols from the *_deb versions of the VMS libraries.
#
#   2) Run the python program:
#         python beam_lib_routine_check
#
# Note: to get a merged list of symbols, use the -dump argument on the command line:
#         python beam_lib_routine_check -dump
#-

import subprocess
import sys

#

vms_sym = dict()
vms_deb_sym = dict()

duplicate = False
dump = False

if len(sys.argv) > 1:
  if sys.argv[1] == '-dump':
    dump = True
  else:
    print 'ERROR: UNKNOWN ARGUMENT: ' + sys.argv[1]

# vms libraries

vms_file = open('vms_lib.symbols', 'r')
for line in vms_file:
  line = line.rstrip() 

  # '!!!' marks the start of a library listing
  if line[0:3] == '!!!': 
    lib_name = line[4:]
    # skip the next 7 lines
    for i in range(7):
      vms_file.next()
    continue

  if line == '': continue      # Ignore blank lines
  vms_sym[line.lower()] = lib_name

vms_file.close()

# vms deb libraries

vms_file = open('vms_lib_deb.symbols', 'r')
for line in vms_file:
  line = line.rstrip() 

  # '!!!' marks the start of a library listing
  if line[0:3] == '!!!': 
    lib_name = line[4:]
    # skip the next 7 lines
    for i in range(7):
      vms_file.next()
    continue

  if line == '': continue      # Ignore blank lines
  vms_deb_sym[line.lower()] = lib_name

vms_file.close()

# mpm_utils

symbols = vms_sym.copy()

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
  print 'No Duplicates!'

for sym in vms_sym:
  if not vms_deb_sym.has_key(sym):
    print 'Symbol in production library but not in debug library: ' + sym
    print '   Library: ' + vms_sym[sym]

for sym in vms_deb_sym:
  if not vms_sym.has_key(sym):
    print 'Symbol in debug library but not in production library: ' + sym
    print '   Library: ' + vms_deb_sym[sym]

if dump:
  for key in sorted(symbols.keys()):
    print '%-30s  %s' % (key, symbols[key])
