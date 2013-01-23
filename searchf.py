#!/usr/bin/env python

import os
import sys
import re

# The idea is to look for a local copy of the library to search.
# We have found a local copy when we find one specific file that we know 
# is in the library.

release_dir   = os.environ['ACC_RELEASE_DIR'] + '/'
dist_dir  = os.environ['DIST_BASE_DIR'] + '/'

#-----------------------------------------

def choose_path (base_dir, base_file, dist_sub_dir):
  if os.path.isfile(base_dir + base_file):                  return base_dir
  if os.path.isfile('../' + base_dir + base_file):          return '../' + base_dir
  if os.path.isfile('../../' + base_dir + base_file):       return '../../' + base_dir
  if os.path.isfile(release_dir + base_dir + base_file):    return release_dir + base_dir
  return dist_dir + dist_sub_dir + base_dir

#-----------------------------------------

def print_help_message ():
  print '''
  Usage:
    getf {options} <search_string>
  Options:
     -a          # Search Numerical recipes, forest, and varies program directories as well.
     -d <dir>    # Search files in <dir> and sub-directories for matches.

  Standard Libraries searched:
     bmad
     sim_utils
     cesr_utils
     mpm_utils
     tao
     bmadz
'''
  sys.exit()

#-----------------------------------------

re_routine = re.compile('^(subroutine|recursive subroutine|elemental subroutine|' + \
      'function|recursive function|real\(rp\) *function|' +  \
      'integer *function|logical *function|interface) ')

re_routine_name = re.compile('\w+')

def routine_here (line2, routine_name):
  match = re_routine.match(line2)
  if match:
    routine_name[0] = re_routine_name.match(line2[match.end(1):].lstrip()).group(0)
    return True
  else:
    return False

#-----------------------------------------

re_interface_end      = re.compile('^end +interface')
re_type               = re.compile(r'^type *\(')
re_module_begin       = re.compile('^module')
re_module_end         = re.compile('^contains')
re_parameter          = re.compile(' parameter.*::')
re_parameter1         = re.compile(r'^\s*([\$\w]+)[\(-:\) ]*=')  # match to: "charge_of(-3:3) = "
re_type_interface_end = re.compile('^end +(type|interface)')
re_end                = re.compile('^end')
re_routine_name_here  = re.compile('^(subroutine|function|interface)')

def search_f90 (file_name, match_str, found_one):

  re_match_str  = re.compile('^' + match_str.lower() + '$')
  re_type_interface_match = re.compile('^(type|interface) +' + match_str.lower() + '\s') 

  found_one_in_this_file = False
  have_printed_file_name = False
  in_module_header = False
  routine_name = ['']

  comments = []

  f90_file = open(file_name)
  while True:
    line = f90_file.readline()
    if line == '': return
    line2 = line.lstrip().lower()
    if line2 == '': continue

    # Skip blank interface blocks

    if line2 == 'interface':
      while True:
        line = f90_file.readline()
        if line == '': return
        line2 = line.lstrip().lower()
        if re_interface_end.search(line2): break

    # Skip "type (" constructs and separator comments.

    if re_type.search(line2): continue
    if line2[0] == '#': continue
    if line2[0:10] == '!---------': continue   # ignore separator comment

    # In the header section of a module

    if re_module_begin.search(line2): in_module_header = True
    if re_module_end.search(line2): in_module_header = False
    
    # Search for parameters

    if in_module_header:
      if re_parameter.search(line2):
        chunks = re_parameter.split(line2)[1].split(',')
        for chunk in chunks:
          chunk_match = re_parameter1.search(chunk)
          if chunk_match:
            param = chunk_match.group(1)
            if re_match_str.search(param) or \
               (param[-1] == '$' and re_match_str.search(param[:-1])):
              found_one[0] = True
              found_one_in_this_file = True
              if not have_printed_file_name:
                print '\nFile:', file_name
                have_printed_file_name = True
              print '    ' + line.rstrip()

      # Add to comment block if a comment

      if line2[0] == '!':
        if full_doc: comments.append(line)
        continue

      # Match to type or interface statement
      # These we type the whole definition

      if re_type_interface_match.search(line2):
        found_one[0] = True
        found_one_in_this_file = True
        print '\nFile:', file_name
        if full_doc:
          for com in comments: print com.rstrip()
          print ''
          print line.rstrip()
          while True:
            line = f90_file.readline()
            if line == '': return
            line2 = line.lstrip().lower()
            print line.rstrip()
            if re_type_interface_end.search(line2): break
        else:
          print line.rstrip()
        comments = []
        continue

      # match to subroutine, function, etc.

      if routine_here(line2, routine_name):
        if re_match_str.match(routine_name[0]):
          found_one[0] = True
          found_one_in_this_file = True
          print '\nFile:', file_name
          if full_doc:
            for com in comments: print com.rstrip()          
          else:
            print '    ', line.rstrip()

        # Skip rest of routine including contained routines

        count = 1
        while True:
          line = f90_file.readline()
          if line == '': return
          line2 = line.lstrip().lower()

          if re_end.search(line2):
            if re_routine_name_here.search(line2[4:].lstrip()):
              count -= 1
          elif routine_here(line2, routine_name):
              count += 1

          if count == 0: break

      #

      comments = []


#-----------------------------------------

def search_c (file_name, match_str, found_one):
  return

#-----------------------------------------

def searchit (this_dir, match_str, found_one):

  # Loop over all directories

  for root, dirs, files in os.walk(this_dir):

    # Remove from searching hidden directories plus "production" and "debug" derectories
    i = 0
    while i < len(dirs):
      if dirs[i] == 'production' or dirs[i] == 'debug' or dirs[i][0] == '.': 
        del dirs[i]
      else:
        i += 1

    # Now loop over all files

    for this_file in files:
      if re.search (this_file, '#'): continue
      file_name = os.path.join(root, this_file)
      if this_file[-4:] == '.f90' or this_file[-4:] == '.inc': search_f90(file_name, match_str, found_one)
      if this_file[-4:] == '.cpp' or this_file[-2:] == '.h' or this_file[-2:] == '.c': search_c(file_name, match_str, found_one)

  # End

  return

#-----------------------------------------

bmad_dir          = choose_path ('bmad', '/modules/bmad_struct.f90', '')
cesr_utils_dir    = choose_path ('cesr_utils', '/modules/cesr_utils.f90', '')
sim_utils_dir     = choose_path ('sim_utils', '/interfaces/sim_utils.f90', '')
mpm_utils_dir     = choose_path ('mpm_utils', '/code/butout.f90', '')
recipes_dir       = choose_path ('recipes_f-90_LEPP', '/lib_src/nr.f90', '')
forest_dir        = choose_path ('forest', '/code/i_tpsa.f90', '/packages')
tao_dir           = choose_path ('tao', '/code/tao_struct.f90', '')
bmadz_dir         = choose_path ('bmadz', '/modules/bmadz_struct.f90', '')
nonlin_bpm_dir    = choose_path ('nonlin_bpm', '/code/nonlin_bpm_init.f90', '')
recipes_dir       = choose_path ('recipes_f-90_LEPP', '/lib_src/nr.f90', '')
bsim_dir          = choose_path ('bsim', '/code/bsim_interface.f90', '')
bsim_cesr_dir     = choose_path ('bsim_cesr', '/modules/bsim_cesr_interface.f90', '')
cesr_programs_dir = choose_path ('cesr_programs', '/bmad_to_ing_knob/bmad_to_ing_knob.f90', '')
cesrv_dir         = choose_path ('cesrv', '/code/cesrv_struct.f90', '')
util_programs_dir = choose_path ('util_programs', '/bmad_to_mad_and_xsif/bmad_to_mad_and_xsif.f90', '')
examples_dir      = choose_path ('examples', '/simple_bmad_program/simple_bmad_program.f90', '')

#-----------------------------------------------------------
# Look for arguments

extra_dir = ''
search_all = False

i = 1
while i+1 < len(sys.argv):
  i += 1
  arg = sys.argv[i]

  if arg[0] != '-': break

  if arg == '-d':
    extra_dir = sys.argv[i+1]
    i += 1
    continue

  if arg == '-a':
    search_all = True
    continue

  if arg == '-h':
    print_help_message ()

#----------------------------------------------------------
# Search for a match.

if i >= len(sys.argv): print_help_message()

match_str_in = sys.argv[i]
match_str = match_str_in.replace('*', '\w*') 
found_one = [False]
full_doc = True

if extra_dir != '':
  print 'Searching also: extra_dir\n'
  searchit (extra_dir, match_str, found_one)

bmad_dir = '/home/dcs16/linux_lib/getf_dir'
searchit (bmad_dir, match_str, found_one)
#searchit (sim_utils_dir, match_str, found_one)
#searchit (tao_dir, match_str, found_one)
#searchit (cesr_utils_dir, match_str, found_one)
#searchit (mpm_utils_dir, match_str, found_one)
#searchit (bmadz_dir, match_str, found_one)

if search_all:
  searchit (recipes_dir, match_str, found_one)
  searchit (forest_dir, match_str, found_one)
  searchit (bsim_dir, match_str, found_one)
  searchit (bsim_cesr_dir, match_str, found_one)
  searchit (cesr_programs_dir, match_str, found_one)
  searchit (cesrv_dir, match_str, found_one)
  searchit (util_programs_dir, match_str, found_one)
  searchit (examples_dir, match_str, found_one)

if not found_one[0]:
  print 'Cannot match String! match_str\n'
  print 'Use "-h" command line option to list options.\n'
else:
  print '\n'

