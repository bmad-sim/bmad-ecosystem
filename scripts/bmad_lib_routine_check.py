#!/usr/bin/env python


#+
# Script to check that all library symbols are unique.
# Script checks the libraries:
# On the VMS side:
#   stulib
#   vmpmlib
#   grlib
# On the linux side:
#   $ACC_RELEASE_DIR/production/lib/mpm_utils
#   $ACC_RELEASE_DIR/production/lib/cesr_utils
#   $ACC_RELEASE_DIR/production/lib/sim_utils
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

#----------------------------------------------------------------------
# get_vms_symbols function

def get_vms_symbols (file_name):

  sym_dict = dict()
  vms_file = open(file_name, 'r')
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
    sym_dict[line.lower()] = lib_name

  vms_file.close()
  return sym_dict

#----------------------------------------------------------------------
# get_linux_symbols function

def get_linux_symbols (lib_file, lib_name):
  global duplicate, symbols
  process = subprocess.Popen(['ar t ' + lib_file], shell = True, stdout = subprocess.PIPE)
  mpm_list = process.communicate()[0].replace('_DBL.o', '').split()

  for sym in mpm_list:
    if symbols.has_key(sym):
      print 'Duplicate symbol: ' + sym 
      print '   In: ' + lib_name + ', ' + symbols[sym]
      duplicate = True
    else:
      symbols[sym] = lib_name

#----------------------------------------------------------------------
# main program

duplicate = False
dump = False

if len(sys.argv) > 1:
  if sys.argv[1] == '-dump':
    dump = True
  else:
    print 'ERROR: UNKNOWN ARGUMENT: ' + sys.argv[1]

# vms libraries

vms_sym = get_vms_symbols ('vms_lib.symbols')
vms_deb_sym = get_vms_symbols ('vms_lib_deb.symbols')

# mpm_utils

symbols = vms_sym.copy()
get_linux_symbols('$ACC_RELEASE_DIR/production/lib/libmpm_utils.a', 'mpm_utils') 
get_linux_symbols('$ACC_RELEASE_DIR/production/lib/libcesr_utils.a', 'cesr_utils') 
get_linux_symbols('$ACC_RELEASE_DIR/production/lib/libsim_utils.a', 'sim_utils') 

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
