#!/usr/bin/env python

import os
import sys
import re
from multiprocessing import Pool, Process

# The idea is to look for a local copy of the library to search.
# We have found a local copy when we find one specific file that we know 
# is in the library.

release_dir = ''
dist_dir = ''

if 'ACC_RELEASE_DIR' in os.environ: release_dir   = os.environ['ACC_RELEASE_DIR'] + '/'
if 'DIST_BASE_DIR'   in os.environ: dist_dir      = os.environ['DIST_BASE_DIR'] + '/'

class search_com_class:
  def __init__(self):
    self.found_one      = False
    self.doc            = 'FULL'
    self.match_str      = ''
    self.case_sensitive = False
    self.namelist_file  = ''
    self.file_name_rel_root = ''   # File name relative to the root search directory

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# choose_path function

def choose_path (base_dir, base_file, dist_sub_dir):
  if os.path.isfile(base_dir + base_file):                  return base_dir
  if os.path.isfile('../' + base_dir + base_file):          return '../' + base_dir
  if os.path.isfile('../../' + base_dir + base_file):       return '../../' + base_dir
  if os.path.isfile(release_dir + base_dir + base_file):    return release_dir + base_dir
  if os.path.isfile(dist_dir + dist_sub_dir + base_dir):    return dist_dir + dist_sub_dir + base_dir
  # If release_dir is defined then we should have found the directory.
  if release_dir != '': print 'CANNOT FIND DIRECTORY FOR SEARCHING:', base_dir
  return ''

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# print_help_message function

def print_help_message ():
  print '''
  Usage for getf and listf:
    getf  {options} <search_string>
    listf {options} <search_string>
  Options:
     -a          # Search Numerical recipes, forest, and varies program directories as well.
     -d <dir>    # Search files in <dir> and sub-directories for matches.
     -c          # Case sensitive search.

  Standard Libraries searched:
     bmad
     sim_utils
     cesr_utils
     mpm_utils
     tao
     bmadz
'''
  sys.exit()

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# routine_here function

re_routine = re.compile('(subroutine|recursive subroutine|elemental subroutine|' + \
      'function|recursive function|real\(rp\) *function|' +  \
      'integer *function|logical *function|interface) ')

re_routine_name = re.compile('\w+')

def routine_here (line2, routine_name):
  match = re_routine.match(line2)
  if match:
    name_match = re_routine_name.match(line2[match.end(1):].lstrip())
    if name_match: routine_name[0] = name_match.group(0)
    return True
  else:
    return False

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# search_f90 function

re_blank_interface_begin = re.compile('interface\s*$')
re_interface_end         = re.compile('end +interface')
re_type                  = re.compile(r'type *\(')
re_module_begin          = re.compile('module')
re_module_header_end     = re.compile('contains')
re_parameter             = re.compile(' parameter.*::')
re_parameter1            = re.compile(r'\s*([\$\w]+)[\(-:\) ]*=')  # match to: "charge_of(-3:3) = "
re_type_interface_end    = re.compile('end +(type|interface)')
re_end                   = re.compile('end')
re_routine_name_here     = re.compile('subroutine|function|interface')

def search_f90 (file_name, search_com):

  re_match_str  = re.compile(search_com.match_str.lower() + '$')
  re_type_interface_match = re.compile('^(type|interface) +' + search_com.match_str.lower() + '\s') 

  found_one_in_this_file = False
  have_printed_file_name = False
  in_module_header = False
  routine_name = ['']
  blank_line_found = False

  comments = []

  f90_file = open(file_name)
  while True:
    line = f90_file.readline()
    if line == '': return
    line2 = line.lstrip().lower()
    if line2.rstrip() == '': 
      blank_line_found = True
      continue

    # Skip blank interface blocks

    if re_blank_interface_begin.match(line2):
      while True:
        line = f90_file.readline()
        if line == '': return
        line2 = line.lstrip().lower()
        if re_interface_end.match(line2): break

    # Skip "type (" constructs and separator comments.

    if re_type.match(line2): continue
    if line2[0] == '#': continue
    if line2[0:10] == '!---------': continue   # ignore separator comment

    # In the header section of a module

    if re_module_begin.match(line2): in_module_header = True
    if re_module_header_end.match(line2): in_module_header = False
    
    # Search for parameters

    if in_module_header:
      if re_parameter.search(line2):
        chunks = re_parameter.split(line2)[1].split(',')
        for chunk in chunks:
          chunk_match = re_parameter1.match(chunk)
          if chunk_match:
            if search_com.doc == 'LIST':
              if not have_printed_file_name: search_com.namelist_file.write('\nFile: '  + search_com.file_name_rel_root + '\n')
              have_printed_file_name = True
              search_com.namelist_file.write(chunk_match.group(1) + '\n')
            else:
              param = chunk_match.group(1)
              if re_match_str.match(param) or \
                 (param[-1] == '$' and re_match_str.match(param[:-1])):
                search_com.found_one = True
                found_one_in_this_file = True
                if not have_printed_file_name:
                  print '\nFile:', file_name
                  have_printed_file_name = True
                print '    ' + line.rstrip()
                break

    # Add to comment block if a comment

    if line2[0] == '!':
      if blank_line_found:
        comments = []
        blank_line_found = False
      if search_com.doc == 'FULL': comments.append(line)
      continue

    # Match to type or interface statement
    # These we type the whole definition

    match = re_type_interface_match.match(line2)
    if match:
      search_com.found_one = True
      found_one_in_this_file = True
      if search_com.doc == 'LIST':
        if not have_printed_file_name: search_com.namelist_file.write('\nFile: '  + search_com.file_name_rel_root + '\n')
        have_printed_file_name = True
        search_com.namelist_file.write(match.group(2) + '\n')
      elif search_com.doc == 'FULL':
        print '\nFile:', file_name
        for com in comments: print com.rstrip()
        if len(comments) > 0: print ''
        print line.rstrip()
        while True:
          line = f90_file.readline()
          if line == '': return
          line2 = line.lstrip().lower()
          print line.rstrip()
          if re_type_interface_end.match(line2): break
      else:
        print '\nFile:', file_name
        print '    ', line.rstrip()
      comments = []
      continue

    # match to subroutine, function, etc.

    if routine_here(line2, routine_name):
      if re_match_str.match(routine_name[0]):
        search_com.found_one = True
        found_one_in_this_file = True
        if search_com.doc == 'LIST':
          if not have_printed_file_name: search_com.namelist_file.write('\nFile: '  + search_com.file_name_rel_root + '\n')
          have_printed_file_name = True
          search_com.namelist_file.write(routine_name[0] + '\n')
        elif search_com.doc == 'FULL':
          print '\nFile:', file_name
          for com in comments: print com.rstrip()          
        else:
          print '\nFile:', file_name
          print '    ', line.rstrip()

      # Skip rest of routine including contained routines

      count = 1
      while True:
        line = f90_file.readline()
        if line == '': return
        line2 = line.lstrip().lower()

        if re_end.match(line2):
          if re_routine_name_here.match(line2[4:].lstrip()):
            count -= 1
        elif routine_here(line2, routine_name):
            count += 1

        if count == 0: break

    #

    comments = []


#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# search_c function

re_quote            = re.compile('"|\'')

def search_c (file_name, search_com):

  found_one_in_this_file = False
  in_extended_comment = False
  blank_line_here = False
  n_curly = 0
  comments = []
  lines_after_comments = []
  function_line = ''
  have_printed_file_name = False

  c_file = open(file_name)
  while True:
    line = c_file.readline()
    if line == '': return
    line2 = line.lstrip()
    if line2.rstrip() == '':
      blank_line_here = True
      continue

    # Ignore preprocessor lines

    if line[0] == '#': continue

    # Throw out quoted substrings

    line2.replace(r'\"', '').replace(r"\'", '')
    while True:
      match = re_quote.search(line2)
      if not match: break
      char = match.group(0)      
      ix = line2.find(char, match.end(0))
      if ix == -1: break
      line2 = line2[0:match.start(0)] + line2[ix+1:]

    # Look For multiline comment "/* ... */" construct and remove if present.

    if n_curly == 0:
      if line2[0:2] == '//' or line2[0:2] == '/*' or in_extended_comment: 
        if blank_line_here:
          comments = []
          blank_line_here = False
        comments.append(line)
        lines_after_comments = []
      else:
        lines_after_comments.append(line)


    while True:
      ix_save = 0
      if in_extended_comment:
        ix = line2.find('*/')
        if ix == -1: break
        in_extended_comment = False
        line2 = line2[0:ix_save] + line2[ix+2:]
      else:
        ix = line2.find('/*')
        if ix == -1: break
        in_extended_comment = True
        ix_save = ix

    if line2[0:2] == '//': continue
    if line2.strip() == '': continue

    ix = line2.find('//')
    if ix > -1: line2 = line2[0:ix]

    # Count curly brackets

    for char in line2:

      if n_curly == 0: 
        function_line = function_line + char
        if char == ';': 
          function_line = ''
          comments = []
          lines_after_comments = []

      if char == '{':
        n_curly += 1
        if n_curly == 1:
          if search_com.case_sensitive:
            is_match = re.search(search_com.match_str + '\s*(\(.*\))\s*{', function_line)
          else:
            is_match = re.search(search_com.match_str + '\s*(\(.*\))\s*{', function_line, re.I)
          if is_match:
            search_com.found_one = True
            if search_com.doc == 'LIST':
              if not have_printed_file_name: search_com.namelist_file.write('\nFile: '  + search_com.file_name_rel_root + '\n')
              have_printed_file_name = True
              search_com.namelist_file.write(is_match.group(1) + '\n')
            elif search_com.doc == 'FULL':
              print '\nFile:', file_name
              for com in comments: print com.rstrip()
              for com in lines_after_comments: print com.rstrip()
            else:
              if not found_one_in_this_file: 
                print '\nFile:', file_name
                found_one_in_this_file = True
              for com in lines_after_comments: print '    ', com.rstrip()

      elif char == '}':
        n_curly -= 1
        if n_curly == 0: 
          function_line = ''
          comments = []
          lines_after_comments = []
  return

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# search_file function

def search_file (search_root_dir, file_dir, file_name, search_com):
  if re.search ('#', file_name): return
  if file_name[0] == '.': return
  full_file_name = os.path.join(file_dir, file_name)
  search_com.file_name_rel_root = full_file_name.replace(search_root_dir, '', 1)
  if file_name[-4:] == '.f90' or file_name[-4:] == '.inc': search_f90(full_file_name, search_com)
  if file_name[-4:] == '.cpp' or file_name[-2:] == '.h' or file_name[-2:] == '.c': search_c(full_file_name, search_com)

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# search_tree function

def search_tree (search_root_dir, search_com):

  if search_root_dir == '': return    # Directory not found by choose_path
  if search_root_dir[-1] != '/': search_root_dir = search_root_dir + '/'
  namelist_file = search_root_dir + 'searchf.namelist'

  # Open file for namelist output if needed

  if search_com.doc == 'LIST':
    if os.access(namelist_file, os.W_OK):
      print 'Opening:', namelist_file
      search_com.namelist_file = open(namelist_file, 'w')
    else:
      print 'CANNOT WRITE TO:', namelist_file
      return

  # If there is an existing searchf.namelist file then use this to see if there are matches.

  if search_com.doc != 'LIST' and os.path.isfile(namelist_file):

    f_namelist = open(namelist_file)
    have_searched_file = False

    for line in f_namelist:
      if line == '': return
      if line.strip() == '': continue

      if line[0:5] == 'File:':
        file = line[6:].strip().rsplit('/', 1)
        have_searched_file = False
        if len(file) == 1:     # No directory spec
          this_root_dir = search_root_dir
          this_file = file[0]
        else:
          this_root_dir = search_root_dir + file[0]
          this_file = file[1]
        continue
      elif have_searched_file:
        continue

      if re.search(search_com.match_str, line):
        search_file (search_root_dir, this_root_dir, this_file, search_com)
        have_searched_file = True

    return

  # No searchf.namelist: Loop over all directories

  for this_root_dir, sub_dirs, files in os.walk(search_root_dir):

    # Remove from searching hidden directories plus "production" and "debug" derectories
    i = 0
    while i < len(sub_dirs):
      if sub_dirs[i] == 'production' or sub_dirs[i] == 'debug' or sub_dirs[i][0] == '.': 
        del sub_dirs[i]
      else:
        i += 1

    # Now loop over all files

    for this_file in files: search_file (search_root_dir, this_root_dir, this_file, search_com)

  # End

  return

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Main routine

def search_all (doc_type):

  search_com = search_com_class()
  search_com.found_one = False
  search_com.doc = doc_type

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

  i = 0
  while i < len(sys.argv):
    i += 1
    if i >= len(sys.argv): break
    arg = sys.argv[i]

    if arg[0] != '-': break

    if arg == '-d':
      extra_dir = sys.argv[i+1]
      i += 1
      continue

    if arg == '-a':
      search_all = True
      continue

    if arg == '-c':
      search_com.case_sensitive = True
      continue

    if arg == '-h':
      print_help_message ()

    print '!!! UNKNOWN ARGUMENT:', arg
    print_help_message ()

  #----------------------------------------------------------
  # Setup dir_list list, etc

  dir_list = [bmad_dir, sim_utils_dir, tao_dir, cesr_utils_dir, mpm_utils_dir, bmadz_dir]
  if search_all:
    dir_list.extend([recipes_dir, forest_dir, bsim_dir, bsim_cesr_dir, cesr_programs_dir, 
                     cesrv_dir, util_programs_dir, examples_dir])
  if extra_dir != '':
    print 'Searching directory:', extra_dir
    dir_list = [extra_dir]

  if search_com.doc == 'LIST':
    search_com.match_str = '(\w+)'
    if i > 0 and i < len(sys.argv): dir_list = [sys.argv[i]]
  else:
    if i == 0 or i >= len(sys.argv): 
      print 'NO SEARCH STRING FOUND!'
      print_help_message()  # Nothing to match to
    match_str_in = sys.argv[i]
    search_com.match_str = match_str_in.replace('*', '\w*') 

  # Search for a match.

  for dir in dir_list:
    search_tree (dir, search_com)

  # And finish

  if not search_com.found_one:
    print 'Cannot match String:',  match_str_in
    print 'Use "-h" command line option to list options.'
  else:
    print ''

