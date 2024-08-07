#!/usr/bin/env python3

#+
# Script to run regression tests for the Bmad libraries.
#-

import re
import os
import sys
import time
import math
from shutil import which

num_tests = 0
num_failures = 0
num_flow_failures = 0
num_programs = 0
pass_all_tests = True

warning_color = '\033[91m\033[1m'   # Red + Bold
normal_color = '\033[0m'

#----------------------------------------------------------

def print_all(string, terminate = False, color = False, failing = False):
  global pass_all_tests
  if failing: pass_all_tests = False

  results.write(string + '\n')
  if color:
    print(warning_color + string + normal_color)
  else:
    print(string)

  if terminate: 
    string2 = '     Flow Failure. Stopping here for this regression.'
    results.write(string2 + '\n')
    print(string2)
    global num_flow_failures
    num_flow_failures += 1

#----------------------------------------------------------
def print_help():
  print('''
Usage:
   run_test.py {-bin <bin_dir>} {-debug} {-test <test_dir>} {-list <test_list_file>}
Note: Do not use -debug with -bin
Defaults:
   <bin_dir>  = "../production/bin" ! Relative to current directory.
              = "../debug/bin"      ! If -debug switch is present
   <test_dir> = ""                  ! For running a single test. Overrides test.list list.
   <test_list_file> = "test.list"   ! For running multiple tests.''')
  exit()

#----------------------------------------------------------
# List of tests is in "test.list".

results = open('regression.results', 'w')

bin_dir = '../production/bin/'
test_dir_list = []
test_list_file = 'TESTS.LIST'
time0 = time.time()

i = 1
while i < len(sys.argv):
  if sys.argv[i] == '-bin':
    bin_dir = sys.argv[i+1]
    i += 1
  elif sys.argv[i] == '-test':
    test_dir_list = [sys.argv[i+1]]
    i += 1
  elif sys.argv[i] == '-list':
    test_list_file = [sys.argv[i+1]]
    i += 1
  elif sys.argv[i] == '-debug':
    bin_dir = '../debug/bin'
  else:
    print_help()

  i += 1

if bin_dir[0] != '/' and bin_dir[0] != '$': bin_dir = '../' + bin_dir
if bin_dir[-1] != '/': bin_dir = bin_dir + '/'
if len(test_dir_list) == 1 and test_dir_list[0] == 'all': test_dir_list = []

if len(test_dir_list) == 0:
  dir_file = open (test_list_file, 'r')
  test_dir_list = dir_file.readlines()

#-------------------------------------------------------------

for test_dir in test_dir_list:
  time0_test = time.time()
  test_dir = test_dir.strip()
  ix = test_dir.find('!')
  if ix != -1: test_dir = test_dir[:ix]
  if len(test_dir) == 0: continue

  # Is this a note:

  if test_dir[:5] == 'NOTE:':
    print_all ('Note in TESTS.LIST file: ' + test_dir, False, True, False)
    continue

  #-----------------------------------------------------------
  # Run the programs

  dir_split = test_dir.split()
  num_programs += 1
  num_local_tests = 0
  num_local_failures = 0

  if len(dir_split) > 2:
    print_all ('\nExtra stuff on line in "TESTS.LIST": ' + dir_split, True, True, True)
    continue

  max_fail = 0
  if len(dir_split) == 2: max_fail = int(dir_split[1])

  subdir = dir_split[0]
  if subdir[-1] == "/": subdir = subdir[:-1]

  if not os.path.exists(subdir):
    print_all ('\nNon-existant subdirectory given in "TESTS.LIST": ' + subdir, True, True, True)
    continue

  os.chdir(subdir)

  print_all ('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print_all ('Starting testing in subdirectory: ' + subdir)

  # Remove output.now

  if os.path.exists('output.now'):
    os.system('rm -f output.now')
#  if pathlib.Path('output.now'.exists('output.now'):
#     os.system('rm output.now')
#  else:
#    print ('did not find output.now')

  # Run process and make sure output.now has been created

  program = subdir
  if program[-1] == '/': program = program[:-1]

  # run.py
  if os.path.exists('run.py'):
    if which('python3'):
      print_all ('     Found run.py. Running this script with python3.')        
      os.system('python3 run.py ' + bin_dir)
    else:
      print_all ('     Found run.py. Running this script with python(2.#).')
      os.system('python run.py ' + bin_dir)

  else:
    program = bin_dir + program
    print_all ('     Running program: ' + program)
  
    if not os.path.isfile(program):
      print_all ('     !!! Program does not exist!', True, True, True)
      os.chdir('..')
      continue

    os.system(program)

  # Look for output

  if not os.path.isfile('output.now'):
    print_all ('    ' + subdir + ': !!! Program failed to create "output.now" file', True, True, True)
    os.chdir('..')
    continue

  if not os.path.isfile('output.correct'):
    print_all ('    ' + subdir + ': !!! No "output.correct" file', True, True, True)
    os.chdir('..')
    continue

  # Compare the output of the program "output.now" to the expected output "output.correct"

  now_file = open('output.now', 'r')  
  correct_file = open('output.correct', 'r')
  test_count = 0

  while True:

    while True:
      now_line = now_file.readline()
      if len(now_line) == 0: break                # End of file
      if len(now_line.strip()) == 0: continue     # Skip blank line
      if now_line.strip()[0] == '!': continue     # Skip comment line
      break

    while True:
      correct_line = correct_file.readline()
      if len(correct_line) == 0: break                # End of file
      if len(correct_line.strip()) == 0: continue     # Skip blank line
      if correct_line.strip()[0] == '!': continue     # Skip comment line
      break

    if len(now_line) == 0 or len(correct_line) == 0: 
      print_all ('')
      if len(now_line) != 0:
        print_all ('    ' + subdir + ': Confusion! End of "output.correct" reached before End of "output.now"', True, True, True)
      if len(correct_line) != 0:
        print_all ('    ' + subdir + ': Confusion! End of "output.now" reached before End of "output.correct"', True, True, True)
      break

    now_line = now_line.strip()
    correct_line = correct_line.strip()

    now_split = now_line.split('"', 2)
    correct_split = correct_line.split('"', 2)
    
    if now_split[0] != '' or len(now_split) != 3:
      print_all ('    ' + subdir + ': Cannot parse line from "output.now": ' + now_line, True, True, True)
      break

    if correct_split[0] != '' or len(correct_split) != 3:
      print_all ('    ' + subdir + ': Cannot parse line from "output.correct": ' + correct_line, True, True, True)
      break

    if now_split[1] != correct_split[1]:
      print_all ('    ' + subdir + ': Identification string for a line in "output.now":    ' + now_split[1], False, True, True)
      print_all ('    ' + subdir + ': Does not match corresponding ID in "output.correct": ' + correct_split[1], True, True, True)

    now_end = now_split[2].strip().split()

    #----------------------------------------------
    # String test

    num_local_tests += 1

    if now_end[0] == 'STR':
      now2_split = now_split[2].split('"')
      correct2_split = correct_split[2].split('"')[1:]

      if len(now2_split) < 2:
        print_all ('    ' + subdir + ': Bad line line "output.now": ' + now_line, True, True, True)
        break

      now2_split.pop(0)    # Get rid of STR item.

      if len(now2_split) != len(correct2_split):
        print_all ('    ' + subdir + ': Number of components in "output.now" line: ' + now_line, False, True, True)
        print_all ('    ' + subdir + ': Does not match number in "output.correct:  ' + correct_line, True, True, True)
        break

      for ix, (now1, correct1) in enumerate(list(zip(now2_split, correct2_split))):
        if now1 != correct1:
          print_all ('')
          if len(now2_split) == 2:     # Will always have blank item in list.
            print_all ('    ' + subdir + ': Regression test failed:', color = True)
          else:
            print_all ('    ' + subdir + ': Regression test failed for datum number: ' + str(ix+1), color = True)

          print_all ('          Line from "output.now": ' + now_line, color = True)
          print_all ('          Line from "output.correct": ' + correct_line, color = True)
          num_local_failures += 1
          break

    #----------------------------------------------
    # Real test

    elif now_end[0] == 'ABS' or now_end[0] == 'REL' or now_end[0] == 'VEC_REL':
      now2_split = now_split[2].strip().split()
      correct2_split = correct_split[2].strip().split()[2:]   # [2:] -> Throw away EG: "ABS 2E-7"
      
      if len(now2_split) < 3:
        print_all ('    ' + subdir + ': Bad line in "output.now": ' + now_line, True, True, True)
        break

      tol_type = now2_split.pop(0)           # Pop REL or ABS item.
      tol_val  = float(now2_split.pop(0))    # Pop tollerance

      if len(now2_split) != len(correct2_split):
        print_all ('    ' + subdir + ': Number of components in "output.now" line: ' + now_line, False, True, True)
        print_all ('    ' + subdir + ': Does not match number in "output.correct:  ' + correct_line, True, True, True)
        break

      bad_at = -1
      bad_diff_val  = 0

      if tol_type == 'VEC_REL':
        vec_amp = 0
        for ix, (now1, correct1) in enumerate(list(zip(now2_split, correct2_split))):
          now_val = float(now1)
          correct_val = float(correct1)
          vec_amp += ((abs(now_val) + abs(correct_val)) / 2) ** 2
        vec_amp = math.sqrt(vec_amp)

      for ix, (now1, correct1) in enumerate(list(zip(now2_split, correct2_split))):
        try:
          now_val = float(now1)
        except:
          now_val = 1e100
        correct_val = float(correct1)
        diff_val = abs(now_val - correct_val)
        abs_val = (abs(now_val) + abs(correct_val)) / 2
        factor = 1
        if tol_type == 'REL': factor = abs_val
        if tol_type == 'VEC_REL': factor = vec_amp

        if diff_val > factor * tol_val and diff_val > bad_diff_val: 
          bad_at = ix
          bad_diff_val = diff_val
          bad_abs_val = abs_val

      if bad_at > -1:
        print_all ('')
        if now_end[0] == 'STR':
          print_all ('    ' + subdir + ': Regression test failed for: "' + now_split[1] + '"', color = True)
        else:
          print_all ('    ' + subdir + ': Regression test failed for: "' + now_split[1] + '"   ' + now_end[0] + '   ' + now_end[1], color = True)

        if len(now2_split) != 1: 
          print_all ('     Regression test failed for datum number: ' + str(bad_at+1), color = True)

        print_all ('        Data from "output.now":     ' + str(now2_split), color = True)
        print_all ('        Data from "output.correct": ' + str(correct2_split), color = True)
        print_all ('        Diff: ' + str(bad_diff_val) + '  Diff/Val: ' + str(abs(bad_diff_val) / bad_abs_val), color = True)
        num_local_failures += 1

    #----------------------------------------------
    # Error test

    else:
      print_all ('     Bad data ID string in "output.now" file: ' + now_line, False, True, True)
      print_all ('     Should be one of: STR, REL, or ABS.', True, True, True)
      break

  #------------------

  print_all ('    ' + subdir + ': Number of tests:        ' + str(num_local_tests))
  print_all ('     Number of failed tests: ' + str(num_local_failures), False, color = (num_local_failures != 0))
  print_all ('     Duration of test (sec): ' + str(time.time() - time0_test))
  print_all ('     Maximum allowed failed tests: ' + str(max_fail))
  if num_local_failures > max_fail: 
    print_all ('     Grade for tests in subdirectory ' + subdir + ': FAILED!', False, True, True)
  else:
    print_all ('     Grade for tests in subdirectory ' + subdir + ': Passed.')


  num_tests += num_local_tests
  num_failures += num_local_failures

  os.chdir('..')
  
#------------------------------------------------------------

print_all ('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print_all ('Total number of tests:           ' + str(num_tests))
print_all ('Total number of failed tests:    ' + str(num_failures), color = (num_failures != 0))
print_all ('Number of Program flow failures: ' + str(num_flow_failures), color = (num_flow_failures != 0))
print_all ('Duration of all tests (sec): %5.2f' % (time.time() - time0))

print('Results file: regression.results')

if pass_all_tests:
  print_all ('\nBottom line for all tests: The code PASSES regression testing.')
  exit(0)
else:
  print_all ('\nBottom line for all tests: The code FAILS regression testing.', color = True)
  exit(1)

results.close()

