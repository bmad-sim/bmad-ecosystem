#!/usr/bin/python

#+
# Script to run regression tests for the Bmad libraries.
#-

import re
import os
import sys

num_tests = 0
num_failures = 0
num_flow_failures = 0
num_programs = 0

#----------------------------------------------------------

def print_all(str, terminate = False):
  results.write(str + '\n')
  print str
  if terminate: 
    str2 = '     Stopping here for this regression.'
    results.write(str2 + '\n')
    print str2
    global num_flow_failures
    num_flow_failures += 1

#----------------------------------------------------------
# List of tests is in "test.list".

if len(sys.argv) > 2: 
  print 'Usage:'
  print '   run_test.py {<program_dir>}'
  exit()

results = open('regression.results', 'w')


program_dir = '../../bin/'
if len(sys.argv) == 2: program_dir = sys.argv[1]
if program_dir[-1] != '/': program_dir = program_dir + '/'
dir_file = open ('tests.list', 'r')
dir_list = dir_file.readlines()

for line in dir_list:
  line = line.strip()
  if len(line) == 0: continue

  #-----------------------------------------------------------
  # Run the programs

  dir_split = line.split()
  num_programs += 1
  num_local_tests = 0
  num_local_failures = 0

  if len(dir_split) > 2:
    print_all ('     Extra stuff on line in "tests.list": ' + dir_split, True)
    continue

  if not os.path.exists(dir_split[0]):
    print_all ('     Non-existant subdirectory given in "tests.list": ' + dir_split[0], True)
    continue

  os.chdir(dir_split[0])

  # Remove output.now

  os.system('rm output.now')

  # Run process and make sure output.now has been created

  program = dir_split[0]
  if len(dir_split) == 2: program = dir_split[1]
  program = program_dir + program

  print_all ('\nStarting testing in subdirectory: ' + dir_split[0])
  print_all ('     Running program: ' + program)
  
  if not os.path.isfile(program):
    print_all ('    !!! Program does not exist!', True)
    continue

  os.system(program)

  if not os.path.isfile('output.now'):
    print_all ('    !!! Program failed to create "output.now" file', True)
    continue

  if not os.path.isfile('output.correct'):
    print_all ('    !!! No "output.correct" file', True)
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
      if len(now_line) != 0:
        print_all ('     Confusion! End of "output.correct" reached before End of "output.now"', True)
      if len(correct_line) != 0:
        print_all ('     Confusion! End of "output.now" reached before End of "output.correct"'), True
      break

    now_line = now_line.strip()
    correct_line = correct_line.strip()

    now_split = now_line.split('"', 2)
    correct_split = correct_line.split('"', 2)
    
    if now_split[0] != '' or len(now_split) != 3:
      print_all ('     Cannot parse line from "output.now": ' + now_line, True)
      break

    if correct_split[0] != '' or len(correct_split) != 3:
      print_all ('     Cannot parse line from "output.correct": ' + correct_line, True)
      break

    if now_split[1] != correct_split[1]:
      print_all ('     Identification string for a line in "output.now":    ' + now_split[1])
      print_all ('     Does not match corresponding ID in "output.correct": ' + corect_split[1], True)

    now_end = now_split[2].strip().split()

    #----------------------------------------------
    # String test

    num_local_tests += 1

    if now_end[0] == 'STR':
      now2_split = now_split[2].split('"')
      correct2_split = correct_split[2].split('"')[1:]

      if len(now2_split) < 2:
        print_all ('     Bad line line "output.now": ' + now_line, True)
        break

      now2_split.pop(0)    # Get rid of STR item.

      if len(now2_split) != len(correct2_split):
        print_all ('     Number of components in "output.now" line: ' + now_line)
        print_all ('     Does not match number in "output.correct:  ' + correct_line, True)
        break

      for ix, (now1, correct1) in enumerate(zip(now2_split, correct2_split)):
        if now1 != correct1:
          if len(now2_split) == 2:     # Will always have blank item in list.
            print_all ('     Regression test failed:')
          else:
            print_all ('     Regression test failed for datum number: ' + str(ix))
          print_all ('          Line from "output.now": ' + now_line)
          print_all ('          Line from "output.correct": ' + correct_line)
          num_local_failures += 1
          break



    #----------------------------------------------
    # Real test

    elif now_end[0] == 'ABS' or now_end[0] == 'REL':
      now2_split = now_split[2].strip().split()
      correct2_split = correct_split[2].strip().split()
      
      if len(now2_split) < 3:
        print_all ('     Bad line line "output.now": ' + now_line, True)
        break

      tol_type = now2_split.pop(0)           # Pop REL or ABS item.
      tol_val  = float(now2_split.pop(0))    # Pop tollerance

      if len(now2_split) != len(correct2_split):
        print_all ('     Number of components in "output.now" line: ' + now_line)
        print_all ('     Does not match number in "output.correct:  ' + correct_line, True)
        break

      for ix, (now1, correct1) in enumerate(zip(now2_split, correct2_split)):
        now_val = float(now1)
        correct_val = float(correct1)
        diff_val = now_val - correct_val
        ave_abs_val = (abs(now_val) + abs(correct_val)) / 2

        ok = True
        if tol_type == 'REL':
          if abs(diff_val) > ave_abs_val * tol_val: ok = False
        else:
          if abs(diff_val) > tol_val: ok = False

        if ok == False:
          if now_end[0] == 'STR':
            print_all ('     Regression test failed for: "' + now_split[1] + '"')
          else:
            print_all ('     Regression test failed for: "' + now_split[1] + '"   ' + now_end[0] + '   ' + now_end[1])
          if len(now2_split) != 1: 
            print_all ('     Regression test failed for datum number: ' + str(ix))
          print_all ('        Data from "output.now":     ' + str(now2_split))
          print_all ('        Data from "output.correct": ' + str(correct2_split))
          print_all ('        Diff: ' + str(diff_val) + '  Diff/Val: ' + str(abs(diff_val) / ave_abs_val))
          num_local_failures += 1
          break

    #----------------------------------------------
    # Error test

    else:
      print_all ('     Bad data ID string in "correct.now" file: ' + now_line)
      print_all ('     Should be one of: STR, REL, or ABS.', True)
      break

  #------------------

  print_all ('     Number of tests:        ' + str(num_local_tests))
  print_all ('     Number of failed tests: ' + str(num_local_failures))

  num_tests += num_local_tests
  num_failures += num_local_failures

  os.chdir('..')
  
#------------------------------------------------------------

print_all ('\n')
print_all ('Total number of tests:           ' + str(num_tests))
print_all ('Total number of failed tests:    ' + str(num_failures))
print_all ('Number of Program flow failures: ' + str(num_flow_failures))

results.close()
print '\nResults file: regression.results\n'
