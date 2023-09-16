#!/usr/bin/python

#+
# sxf_to_bmad.py
# 
# script to convert SXF lattice files to Bmad format.
# SXF lattice files can be constructed from MADX using the SXFWRITE command.
#-

import sys, re, math
import time


start_time = time.time()

class ele_param_struct:
  def __init__(self, name = ''):
    self.name = name
    self.type = ''
    self.value = 0

class ele_struct:
  def __init__(self, name = ''):
    self.name = name
    self.type = ''
    self.param = dict()

class sequence_struct:
  def __init__(self, name = ''):
    self.name = name
    self.length = 0
    self.line = []

#------------------------------------------------------------------
#------------------------------------------------------------------

def param_to_string(param, ele):

  if param.name in ['l', 'e1', 'e2', 'fint', 'hgap', 'ks']:
    return ', ' + param.name + ' = ' + param.value

  if param.name in ['at', 'lrad']: return ''

  if param.name in ['tag']: return ', type = ' + param.value
  if param.name == 'arc' and ele.type == 'rbend': return ', l_arc = ' + param.value
  if param.name == 'arc' and ele.type == 'sbend': return ', l = ' + param.value
  if param.name == 'volt': return ', l = 1e6 * ' + param.value

  length = 0
  if ele.type in ['rbend', 'sbend'] and 'arc' in ele.param: 
    length = float(ele.param['arc'].value)
  elif 'l' in ele.param: 
    length = float(ele.param['l'].value)

  if param.name == 'kl':
    kl_arr = param.value
    line = ''

    if ele.type in ['rbend', 'sbend']:
      line = ', angle = ' + kl_arr[0]
      kl_arr[0] = 0

    if ele.type in ['quadrupole', 'rbend', 'sbend'] and len(kl_arr) > 1 and kl_arr[1] != 0:
      line = line + ', k1 = ' + str(float(kl_arr[1])/length)
      kl_arr[1] = 0
    
    if ele.type in ['sextupole', 'rbend', 'sbend'] and len(kl_arr) > 2 and kl_arr[2] != 0:
      line = line + ', k2 = ' + str(float(kl_arr[2])/length)
      kl_arr[2] = 0
    
    if ele.type in ['octupole'] and len(kl_arr) > 3 and kl_arr[3] != 0:
      line = line + ', k3 = ' + str(float(kl_arr[3])/length)
      kl_arr[3] = 0
    
    for ix, kl in enumerate(kl_arr):
      if float(kl) == 0: continue
      if ele.type in ['multipole']:
        line = line + ', k' + str(ix) + ' = ' + kl
      else:
        line = line + ', b' + str(ix) + ' = ' + str(float(kl)/math.factorial(ix))

    return line


  error_exit('UNKNOWN ELEMENT PARAMETER: ' + param.name)

#------------------------------------------------------------------
#------------------------------------------------------------------

def print_help():
  print (''' \
Syntax: 
  sxf_to_bmad.py <xsf_lattice_file>
''')
  sys.exit()

#-------------------------------------------------------------------
#------------------------------------------------------------------

def token_is_name_check(token, err_str, line):
  if re.match('^[\w.]+$', token): return
  error_exit (err_str, line)

#-------------------------------------------------------------------
#------------------------------------------------------------------

def token_matches_check(token, match_to, err_str, line):
  if token == match_to: return
  error_exit (err_str, line)

#-------------------------------------------------------------------
#------------------------------------------------------------------

def error_exit (err_str, line = ''):
  print (err_str)
  if line != '': print ('ERROR AT OR NEAR LINE: ' + line)
  sys.exit()

#-------------------------------------------------------------------
#------------------------------------------------------------------

def pop_token (token_list):
  while True:
    if len(token_list) == 0: return None
    token = token_list.pop(0)
    if token == '': continue   # Happens when two delims are next to one another: "(["
    if token == ' ': continue
    return token

#-------------------------------------------------------------------
#------------------------------------------------------------------

def WrapWrite(f_out, line):
  MAXLEN = 120
  tab = ''

  while True:
    if len(line) <= MAXLEN:
      f_out.write(tab + line + '\n')
      return

    ix = line[:MAXLEN].rfind(',')

    if ix == -1: 
      ix = line[:MAXLEN].rfind(' ')
      f_out.write(tab + line[:ix+1] + ' &\n')
    else:
      f_out.write(tab + line[:ix+1] + '\n')  # Don't need '&' after a comma

    tab = '         '
    line = line[ix+1:]

#-------------------------------------------------------------------
#------------------------------------------------------------------
# Main program.

if len(sys.argv) == 2:
  sxf_lat_file = sys.argv[1]
elif len(sys.argv) > 2:
  print_help()

f_in = open(sxf_lat_file, 'r')

#------------------------------------------------------------------
# Read in sxf file.

token_list = []
sequence_status = 'start'
ele_status = 'start'
param_status = 'start'
param_stack = []

reading_file = True
seq = sequence_struct

while True:

  while len(token_list) < 20 and reading_file:
    line = f_in.readline()
    if len(line) == 0: reading_file = False
    line = line.partition('//')[0].rstrip()   # Remove comments
    token_list.extend(re.split(r'(,|=|\(|\)|\[|\]|;\{|\}| )\s*', line))

  token = pop_token(token_list)
  if token == None: break

  if sequence_status == 'end':
    error_exit ('EXTRA STUFF AT END', line)

  if sequence_status == 'start':
    seq = sequence_struct(token)
    token_matches_check (pop_token(token_list), 'sequence', 'NO SEQUENCE WORD FOUND AT BEGINNING OF FILE!', line)
    token_matches_check (pop_token(token_list), '{', 'NO BEGINNING "{" FOUND FOR SEQUENCE', line)
    sequence_status = 'in_body'

  elif ele_status == 'start' and token == 'endsequence':
    token_matches_check (pop_token(token_list), 'at', 'NO "at" FOUND FOLLOWING "endsequence"', line)
    token_matches_check (pop_token(token_list), '=', 'NO "=" FOUND FOLLOWING "endsequence at"', line)
    seq.length = pop_token(token_list)
    token_matches_check (pop_token(token_list), '}', 'NO "}" FOUND AT END OF SEQUENCE', line)
    sequence_status = 'end'

  elif ele_status == 'start':
    token_is_name_check (token, 'EXPECTING AN ELEMENT NAME BUT GOT: ' + token, line)
    ele = ele_struct(token)
    ele.type = pop_token(token_list)
    seq.line.append(ele)
    ele_status == 'found_type'
    token_matches_check (pop_token(token_list), '{', 'NO BEGINNING "{" FOUND FOR ELEMENT: ' + ele.name + ' GOT: ' + token, line)
    ele_status = 'in_body'

  elif param_status == 'start' and token == '}':
    if len(param_stack) == 0:
      token_matches_check (pop_token(token_list), ';', 'NO CLOSING ";" FOUND AT END OF ELEMENT: ' + ele.name, line)
      ele_status = 'start'

    else:
      parameter = param_stack.pop()
      param0_dict = parameter.value

  elif param_status == 'start':
    parameter = ele_param_struct(token)
    if len(param_stack) == 0:
      param0_dict = ele.param
    else:
      param0_dict = param_stack[-1].value

    token_matches_check (pop_token(token_list), '=', 'NO "=" FOR PARAMETER: ' + parameter.name + ' IN ELEMENT: ' + ele.name + ' GOT: ' + token, line)
    param_status = 'in_value'

  elif param_status == 'in_value':
    param_status = 'start'

    if token == '{':
      param0_dict[parameter.name] = parameter
      parameter.value = dict()
      parameter.type = 'container'
      param0_dict = parameter.value
      param_stack.append(parameter)

    elif token == '[':
      value_array = []
      while True:
        token = pop_token(token_list)
        if token != ']':
          value_array.append(token)
          continue

        if parameter.name in param0_dict:
          v0_array = param0_dict[parameter.name].value
          for i in range(len(value_array)):
            if i == len(v0_array):
              v0_array.append(value_array[i])
            elif value_array[i] != '0':
              v0_array[i] = value_array[i]

        else:
          parameter.value = value_array
          parameter.type = 'array'
          param0_dict[parameter.name] = parameter

        break

    else:
      parameter.value = token
      parameter.type = 'scalar'
      param0_dict[parameter.name] = parameter

#-----------------------------------------------------------
print_seq = False

if print_seq:
  print ('')
  print ('Sequence name: ' + seq.name)
  print ('Number elements: ' + str(len(seq.line)))
  print ('Sequence length: ' + str(seq.length))
  print ('')

  for ele in seq.line:
    print ('Element: ' + ele.name + '  ' + ele.type)
    for name, param in ele.param.items():
      if param.type == 'scalar' or param.type == 'array':
        print ('  ' + param.name + ' = ' + str(param.value))
      else:
        print ('  ' + param.name + ':')
        for name2, param2 in param.value.items():
          print ('    ' + param2.name + ' = ' + str(param2.value))


#-----------------------------------------------------------
# Write Bmad file

if sxf_lat_file.find('sxf') != -1:
  bmad_lat_file = sxf_lat_file.replace('sxf', 'bmad')
elif sxf_lat_file.find('Sxf') != -1:
  bmad_lat_file = sxf_lat_file.replace('Sxf', 'bmad')
elif sxf_lat_file.find('SXF') != -1:
  bmad_lat_file = sxf_lat_file.replace('SXF', 'bmad')
else:
  bmad_lat_file = sxf_lat_file + '.bmad'

f_out = open(bmad_lat_file, 'w')
f_out.write ('! Translated from SXF file: ' + sxf_lat_file + '\n\n')

f_out.write ('sequence_drift: drift, l = ' + seq.length + '\n')
f_out.write ('machine: line = (sequence_drift)\n')
f_out.write ('use, machine\n')

for ele in seq.line:
  line = ele.name + ': ' + ele.type
  if 'at' in ele.param:
    line = line + ', superimpose, offset = ' + ele.param['at'].value
  else:
    line = line + ', superimpose, ref = ' + old_ele.name + 'ref_origin = end, ele_origin = beginning'

  for name, param in ele.param.items():
    if param.type == 'container':
      for name2, param2 in param.value.items():
        line = line + param_to_string(param2, ele)
    else:
      line = line + param_to_string(param, ele)

  f_out.write('\n')
  WrapWrite(f_out, line)
  old_ele = ele
