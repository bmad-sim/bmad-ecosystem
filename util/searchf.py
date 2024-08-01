#!/usr/bin/env python3

#+
# searchf.py is the file that holds the code for the scripts:
#   getf
#   listf
#   create_searchf_namelist
#
# See the Bmad manual for a description of listf and getf.
# See create_searchf_namelist for documentation on the searchf.namelist files
#-

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
    self.doc_type       = 'FULL'   # (for getf), 'SHORT' (for listf), 'LIST' (for create_searchf_namelist), or 'RAW'
    self.match_str      = ''
    self.case_sensitive = False
    self.namelist_file  = ''
    self.file_name_rel_root = ''   # File name relative to the root search directory
    self.search_only_for = ''

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# choose_path function

def choose_path (dir_list, root_dir, base_dir, base_file, release_sub_dir):
  this_dir = ''

  if root_dir == '':
    if os.path.isfile(base_dir + base_file):
      this_dir = base_dir

    elif os.path.isfile('../' + base_dir + base_file):
      this_dir = '../' + base_dir

    elif os.path.isfile('../../' + base_dir + base_file):
      this_dir = '../../' + base_dir

    elif os.path.isfile(release_sub_dir + base_dir + base_file):
      this_dir = release_sub_dir + base_dir

    elif os.path.isfile(release_dir + release_sub_dir + base_dir + base_file):
      this_dir = release_dir + release_sub_dir + base_dir

    elif os.path.isfile(dist_dir + base_dir + base_file):
      this_dir = dist_dir + base_dir

    # If release_dir is defined then we should have found the directory.
    elif release_dir != '': 
      print ('CANNOT FIND DIRECTORY FOR SEARCHING: ', base_dir)


  else:
    if root_dir[-1] != '/': root_dir = root_dir + '/'

    if os.path.isfile(root_dir + base_dir + base_file):
      this_dir = root_dir + base_dir

    elif os.path.isfile(root_dir + release_sub_dir + base_dir + base_file):
      this_dir = root_dir + release_sub_dir + base_dir

    else:
      print ('CANNOT FIND DIRECTORY FOR SEARCHING: ', base_dir)

  if this_dir != '': dir_list.append(this_dir)
  return

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# print_help_message function

def print_help_message ():
  print ('''
  Usage for getf and listf:
    getf  {options} <search_string>
    listf {options} <search_string>

  Options:
     -c          # Case sensitive search when searching C/C++ files.
     -d <s_dir>  # Use <s_dir> as the search directory. Will not search standard directories. 
     -h          # Print this help message.
     -r <r_dir>  # Use <r_dir> as the root directory to search for the search directories.
     -s <what>   # Search only for: <what> = "struct", "routine", "parameter", or "module".

  Explanation: getf/listf will search the "Search directories" and any sub-directories
  for files of the type:
      *.f90    *.inc    *.cpp    *.h    *.c
  Within each of these files getf/listf will search for any routine, struct, parameter or 
  module that matches <search_string>. Wild cards "*" and "." may be used. See the Bmad
  manual for more details.

  Note: getf/listf look for search directories locally and then, if not found, look for the
  search directories in a release or distribution. The exception is that if the "-r <r_dir>" option
  is used, getf/listf will only look at the subdirectories of <r_dir> for the search directories.

  Standard Search directories:
      bmad                recipes_f-90_LEPP      
      bsim                sim_utils         
      code_examples       tao
      forest              util_programs

''')
  sys.exit()

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# routine_here function

re_routine = re.compile('(program|subroutine|recursive subroutine|elemental subroutine|' + \
      'function|recursive function|elemental function|real\(rp\) *function|' +  \
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
re_type_var              = re.compile(r'type *\(')
re_type_def              = re.compile(r'type +\w')
re_type_def_end          = re.compile('end +type')
re_module_begin          = re.compile('module')
re_module_header_end     = re.compile('contains')
re_parameter             = re.compile(' parameter\s*::')
re_parameter1            = re.compile(r'\s*([\$\w]+)\s*(\(.+\)|)?\s*=')  # match to: "charge_of(-3:n_charge$) = "
re_type_interface_end    = re.compile('end +(type|interface)')
re_end                   = re.compile('end')
re_routine_name_here     = re.compile('program|subroutine|function|interface')

def search_f90 (file_name, search_com):

  re_match_str  = re.compile(search_com.match_str.lower() + '$')
  re_type_interface_match = re.compile('^(type|interface) +' + search_com.match_str.lower() + '\s') 

  found_one_in_this_file = False
  have_printed_file_name = False
  in_module_header = False
  in_type_def = False
  routine_name = ['']
  blank_line_found = False

  comments = []

  try:
    if sys.version_info[0] < 3:
      f90_file = open(file_name, 'r')
    else:
      f90_file = open(file_name, 'r', encoding = 'ISO-8859-1')
  except:
    print ('Note: Cannot open: ' + file_name)
    return

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

    if re_type_var.match(line2): continue
    if line2[0] == '#': continue
    if line2[0:10] == '!---------': continue   # ignore separator comment
    if line2[:11] == 'recursive &': 
      comments.append(line)
      continue

    line_list = [line2]
    while line2.rstrip()[-1] == '&':
      aline = f90_file.readline()
      line_list.append(aline)
      line2 = line2.rstrip()[:-1] + aline

    # In the header section of a module

    match = re_module_begin.match(line2)
    if match:
      in_module_header = True
      name_match = re.match('\w+', line2[match.end(0):].lstrip())
      if name_match and 'module'.startswith(search_com.search_only_for):
        module_name = name_match.group(0)
        if re_match_str.match(module_name):
          search_com.found_one = True
          found_one_in_this_file = True
          if search_com.doc_type == 'LIST':
            if not have_printed_file_name: search_com.namelist_file.write('\nFile: '  + search_com.file_name_rel_root + '\n')
            have_printed_file_name = True
            search_com.namelist_file.write(module_name + '\n')
          elif search_com.doc_type == 'FULL':
            print ('\nFile: ', file_name)
            for com in comments: print (com.rstrip())
          elif search_com.doc_type == 'SHORT':
            print ('\nFile: ' + file_name)
            print ('    ' + line.rstrip())

    if not in_type_def and re_module_header_end.match(line2): in_module_header = False
    
    # Search for parameters

    if in_module_header:
      if re_parameter.search(line2) and 'parameter'.startswith(search_com.search_only_for):

        chunks = re_parameter.split(line2)[1].split(',')
        for chunk in chunks:
          chunk_match = re_parameter1.match(chunk)
          if chunk_match:
            if search_com.doc_type == 'LIST':
              if not have_printed_file_name: search_com.namelist_file.write('\nFile: '  + search_com.file_name_rel_root + '\n')
              have_printed_file_name = True
              search_com.namelist_file.write(chunk_match.group(1) + '\n')
            elif search_com.doc_type != 'RAW':
              param = chunk_match.group(1)
              if re_match_str.match(param) or (param[-1] == '$' and re_match_str.match(param[:-1])):
                search_com.found_one = True
                found_one_in_this_file = True
                if not have_printed_file_name:
                  print ('\nFile: ' + file_name)
                  have_printed_file_name = True
                for bline in line_list:
                  print ('    ' + bline.rstrip())             
                break
              # endif
            # endif
          # endif
        # endfor
      # endif
    # endif

    # Add to comment block if a comment

    if line2[0] == '!':
      if blank_line_found:
        comments = []
        blank_line_found = False
      if search_com.doc_type == 'FULL': comments.append(line)
      continue

    # Match to type or interface statement
    # These we type the whole definition

    if in_type_def and re_type_def_end.match(line2): in_type_def = False
    if not in_type_def and re_type_def.match(line2): in_type_def = True

    match = re_type_interface_match.match(line2)
    if match and 'struct'.startswith(search_com.search_only_for):
      search_com.found_one = True
      found_one_in_this_file = True

      if search_com.doc_type == 'LIST':
        if not have_printed_file_name: search_com.namelist_file.write('\nFile: '  + search_com.file_name_rel_root + '\n')
        have_printed_file_name = True
        search_com.namelist_file.write(match.group(2) + '\n')
      elif search_com.doc_type == 'FULL':
        print ('\nFile: ' + file_name)
        for com in comments: print (com.rstrip())
        if len(comments) > 0: print ('')
        print (line.rstrip())
        while True:
          line = f90_file.readline()
          if line == '': return
          line2 = line.lstrip().lower()
          print (line.rstrip())
          if re_type_interface_end.match(line2): break
          in_type_def = False
      elif search_com.doc_type == 'SHORT':
        print ('\nFile: ' + file_name)
        print ('    ' + line.rstrip())
      elif search_com.doc_type == 'RAW':
        while True:
          line = f90_file.readline()
          if line == '': return
          line2 = line.lstrip().lower()
          if re_type_interface_end.match(line2): break
          print (line.rstrip())
          in_type_def = False

      comments = []
      continue

    # match to subroutine, function, etc. --------

    if routine_here(line2, routine_name):
      if re_match_str.match(routine_name[0]) and 'routine'.startswith(search_com.search_only_for):
        search_com.found_one = True
        found_one_in_this_file = True
        if search_com.doc_type == 'LIST':
          if not have_printed_file_name: search_com.namelist_file.write('\nFile: '  + search_com.file_name_rel_root + '\n')
          have_printed_file_name = True
          search_com.namelist_file.write(routine_name[0] + '\n')
        elif search_com.doc_type == 'FULL':
          print ('\nFile: ' + file_name)
          for com in comments: print (com.rstrip())
          for aline in line_list: print (aline.rstrip())
        else:
          print ('\nFile: ' + file_name)
          print ('    ' + line.rstrip())

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
    try:
      line = c_file.readline()
    except:
      continue    # If line contains a non-ascii character.
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
            is_match = re.search(' ' + search_com.match_str + '_?\s*(\(.*\))\s*{', function_line)
          else:
            is_match = re.search(' ' + search_com.match_str + '_?\s*(\(.*\))\s*{', function_line, re.I)
          if is_match and 'routine'.startswith(search_com.search_only_for):
            search_com.found_one = True
            if search_com.doc_type == 'LIST':
              if not have_printed_file_name: search_com.namelist_file.write('\nFile: '  + search_com.file_name_rel_root + '\n')
              have_printed_file_name = True
              search_com.namelist_file.write(is_match.group(1) + '\n')
            elif search_com.doc_type == 'FULL':
              print ('\nFile: ' + file_name)
              for com in comments: print (com.rstrip())
              for com in lines_after_comments: print (com.rstrip())
            else:
              if not found_one_in_this_file: 
                print ('\nFile: ' + file_name)
                found_one_in_this_file = True
              for com in lines_after_comments: print ('    ' + com.rstrip())

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

def search_file (search_base_dir, file_dir, file_name, search_com):
  if re.search ('#', file_name): return
  if file_name[0] == '.': return
  full_file_name = os.path.join(file_dir, file_name)
  search_com.file_name_rel_root = full_file_name.replace(search_base_dir, '', 1)
  if file_name[-4:] == '.f90' or file_name[-4:] == '.inc': search_f90(full_file_name, search_com)
  if file_name[-4:] == '.cpp' or file_name[-2:] == '.h' or file_name[-2:] == '.c': search_c(full_file_name, search_com)

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# search_tree function

def search_tree (search_base_dir, search_com):

  if search_base_dir == '': return    # Directory not found by choose_path
  if search_base_dir[-1] != '/': search_base_dir = search_base_dir + '/'
  namelist_file = search_base_dir + 'searchf.namelist'

  # Open file for namelist output if needed

  if search_com.doc_type == 'LIST':
    if os.access(search_base_dir, os.W_OK):
      print ('Creating: ' + namelist_file)
      search_com.namelist_file = open(namelist_file, 'w')
    else:
      print ('CANNOT WRITE TO: ' + namelist_file)
      return

  # If there is an existing searchf.namelist file then use this to see if there are matches.

  if search_com.doc_type != 'LIST' and os.path.isfile(namelist_file):

    f_namelist = open(namelist_file)
    have_searched_file = False

    for line in f_namelist:
      if line == '': return
      if line.strip() == '': continue

      if line[0:5] == 'File:':
        file = line[6:].strip().rsplit('/', 1)
        have_searched_file = False
        if len(file) == 1:     # No directory spec
          this_search_base_dir = search_base_dir
          this_file = file[0]
        else:
          this_search_base_dir = search_base_dir + file[0]
          this_file = file[1]
        continue
      elif have_searched_file:
        continue

      if re.search(search_com.match_str, line):
        search_file (search_base_dir, this_search_base_dir, this_file, search_com)
        have_searched_file = True

    return

  # No searchf.namelist: Loop over all directories

  for this_search_base_dir, sub_dirs, files in os.walk(search_base_dir):

    # Remove from searching hidden directories plus "production" and "debug" derectories
    i = 0
    while i < len(sub_dirs):
      if sub_dirs[i] == 'production' or sub_dirs[i] == 'debug' or sub_dirs[i][0] == '.': 
        del sub_dirs[i]
      else:
        i += 1

    # Now loop over all files

    for this_file in files: search_file (search_base_dir, this_search_base_dir, this_file, search_com)

  # End

  return

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Main routine

def search_all (doc_type):

  search_com = search_com_class()
  search_com.found_one = False
  search_com.doc_type = doc_type

  #-----------------------------------------------------------
  # Look for arguments

  dir_list = []
  root_dir = ''
  search_all = False
  search_com.search_only_for = ''

  i = 0
  while i < len(sys.argv):
    i += 1
    if i >= len(sys.argv): break
    arg = sys.argv[i]

    if i == 1 and len(arg) == 0:
      print_help_message ()

    if arg[0] != '-': break

    if arg == '-a':
      search_all = True    # Not used.
      continue

    if arg == '-c':
      search_com.case_sensitive = True
      continue

    if arg == '-d':
      dir_list = [sys.argv[i+1]]
      i += 1
      continue

    if arg == '-h':
      print_help_message ()

    if arg == '-r':
      root_dir = sys.argv[i+1]
      i += 1
      continue

    if arg == '-s':
      s = sys.argv[i+1] 
      search_com.search_only_for = s
      if not 'struct'.startswith(s) and not 'routine'.startswith(s) and \
         not 'parameter'.startswith(s) and not 'module'.startswith(s):
        print ('-s ARGUMENT NOT CORRECT.')
        print_help_message ()
      i += 1
      continue

    print ('!!! UNKNOWN ARGUMENT: ' + arg)
    print_help_message ()

  #----------------------------------------------------------
  # Setup dir_list list, etc

  if len(dir_list) == 0:    # If no -d command line arg
    choose_path (dir_list, root_dir, 'util_programs', '/mad_to_bmad/madx_to_bmad.py', '')
    choose_path (dir_list, root_dir, 'forest', '/code/i_tpsa.f90', '')
    choose_path (dir_list, root_dir, 'bsim', '/code/bsim_interface.f90', '')
    choose_path (dir_list, root_dir, 'code_examples', '/simple_bmad_program/simple_bmad_program.f90', '')
    choose_path (dir_list, root_dir, 'sim_utils', '/interfaces/sim_utils.f90', '')
    choose_path (dir_list, root_dir, 'tao', '/code/tao_struct.f90', '')
    choose_path (dir_list, root_dir, 'bmad', '/modules/bmad_struct.f90', '')

  if search_com.doc_type == 'LIST':
    search_com.match_str = '(\w+)'
    if i > 0 and i < len(sys.argv): dir_list = [sys.argv[i]]
  else:
    if i == 0 or i >= len(sys.argv): 
      print ('NO SEARCH STRING FOUND!')
      print_help_message()  # Nothing to match to
    match_str_in = sys.argv[i]
    search_com.match_str = match_str_in.replace('*', '\w*') 

  # Search for a match.

  for dir in dir_list:
    search_tree (dir, search_com)

  # And finish

  if search_com.doc_type != 'LIST':
    if not search_com.found_one:
      print ('Cannot match String: ' + match_str_in)
      print ('Use "-h" command line option to list options.')
    else:
      print ('')

