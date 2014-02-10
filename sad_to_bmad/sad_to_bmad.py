#!/usr/bin/python

import sys, getopt, re, textwrap

class ele_struct:
  def __init__(self):
    self.name = ''
    self.type = ''
    self.param = dict()
    self.printed = False

class lat_line_struct:
  def __init__(self):
    self.name = ''
    self.element = []

# ------------------------------------------------------------------

def WrapWrite(line):
  MAXLINELENGTH = 120
  lines = textwrap.wrap(line,MAXLINELENGTH)
  tab = False
  for line in lines[:-1]:
    if tab:
      f_out.write('         '+line+' &\n')
    else:
      f_out.write(line+' &\n')
      tab = True
  if tab:
    f_out.write('         '+lines[-1]+'\n')
  else:
    f_out.write(lines[-1]+'\n')

# ------------------------------------------------------------------

ele_type_to_bmad = {
  'drift': 'drift',
  'bend': 'rbend',
  'quad': 'quadrupole',
  'sext': 'sextupole',
  'oct': 'octupole',
  'mult': 'sad_mult',
  'sol': 'solenoid',
  'cavi': 'rfcavity',
  'moni': 'marker',
  'mark': 'marker',
  'beambeam': 'marker',
  'apert': 'marker',
}

ele_param_translate = {
    'l': 'l',
    'dx': 'x_offset',
    'dy': 'y_offset',
    'dz': 'z_offset',
    'radius': 'aperture',
    'angle': 'angle',
    'e1': 'e1',
    'e2': 'e2',
    'bz': 'bs_field',
    'f1': 'f1',
    'f2': 'f2',
    'eps': 'eps_step_scale',
    'freq': 'rf_frequency',
    'volt': 'voltage',
    'bound': 'bound',
    'geo': 'geo',
    'rotate': 'tilt', 
    'chi1': 'x_pitch',
    'chi2': 'y_pitch',
    'chi3': 'tilt',
    'k0': 'b0', 'k1': 'b1', 'k2': 'b2', 'k3': 'b3', 'k4': 'b4', 'k5': 'b5', 'k6': 'b6', 'k7': 'b7', 
    'k8': 'b8', 'k9': 'b9', 'k10': 'b10', 'k11': 'b11', 'k12': 'b12', 'k13': 'b13', 'k14': 'b14', 
    'k15': 'b15', 'k16': 'b16', 'k17': 'b17', 'k18': 'b18', 'k19': 'b19', 'k20': 'b20', 'k21'
    'sk0': 'a0', 'sk1': 'a1', 'sk2': 'a2', 'sk3': 'a3', 'sk4': 'a4', 'sk5': 'a5', 'sk6': 'a6', 'sk7': 'a7', 
    'sk8': 'a8', 'sk9': 'a9', 'sk10': 'a10', 'sk11': 'a11', 'sk12': 'a12', 'sk13': 'a13', 'sk14': 'a14', 
    'sk15': 'a15', 'sk16': 'a16', 'sk17': 'a17', 'sk18': 'a18', 'sk19': 'a19', 'sk20': 'a20', 'sk21': 'a21 ',
}

ignore_sad_param = ['ldev', 'fringe', 'disfrin', 'disrad', 'r1', 'r2', 'r3', 'r4', 'betax', 'betay']

# ------------------------------------------------------------------

def sad_ele_to_bmad (sad_ele, bmad_ele, inside_sol, bz):

  bmad_ele.name = sad_ele.name

  if not sad_ele.type in ele_type_to_bmad:
    print ('TYPE OF ELEMENT NOT RECOGNIZED: ' + sad_ele.type + '\n' + '    FOR ELEMENT: ' + sad_ele.name)
    return

  bmad_ele.type = ele_type_to_bmad[sad_ele.type]

  if inside_sol: 
    bmad_ele.type = 'sad_mult'
    bmad_ele.param['bs_field'] = bz

  for sad_name in sad_ele.param:
    if sad_name in ignore_sad_param: continue
    if sad_name not in ele_param_translate:
      print ('SAD PARAMETER NOT RECOGNIZED: ' + sad_name + '\n' + '    DEFINED IN ELEMENT: ' + sad_ele.name)
      continue

    bmad_name = ele_param_translate[sad_name]
    if sad_ele.type == 'bend' and sad_name == 'rotate': bmad_name = 'ref_tilt'
    if bmad_ele.type != 'sad_mult' and sad_name == 'k1': bmad_name = 'k1'
    if bmad_ele.type != 'sad_mult' and sad_name == 'k2': bmad_name = 'k2'

    value = sad_ele.param[sad_name]
    if 'deg' in value: value.replace('deg', '* degrees')
    if 'kev' in value: value.replace('kev', '* 1e3')
    if 'mev' in value: value.replace('mev', '* 1e6')
    if 'gev' in value: value.replace('gev', '* 1e9')

    bmad_ele.param[bmad_name] = value


  fringe = '0'   # default 
  if 'fringe' in sad_ele.param: fringe = sad_ele.param['fringe']

  if fringe == '1':
    bmad_ele.param['fringe_at'] = 'entrance_end'
  elif fringe == '2':
    bmad_ele.param['fringe_at'] = 'exit_end'

  disfrin = '0'
  if 'disfrin' in sad_ele.param: disfrin = sad_ele.param['disfrin']

  if bmad_ele.type == 'sad_mult':
    if fringe == '0':
      if disfrin != '0': bmad_ele.param['fringe_kind'] = 'none'
    elif disfrin == '0': 
      bmad_ele.param['fringe_kind'] = 'full'
    else:
      bmad_ele.param['fringe_kind'] = 'linear'

# ------------------------------------------------------------------

def parse_line(definitions, lat_line_list):
  this_line = lat_line_struct()
  line_name, e_sign, line_list = definitions.partition("=")
  line_name = line_name.strip()
  line_list = line_list.strip()
  line_list = line_list[1:-1]    # Remove parenteses from "(...)"
  
  this_line.name = line_name
  for elename in re.split('\W+', line_list):
    this_line.element.append(elename)

  lat_line_list[line_name] = this_line

# ------------------------------------------------------------------

regexDef = r"(\S*?)\s*=\s*\((.*?)\)"
regexParam = r"(\S*)\s*=\s*(\S*\s*(deg)?)"

def parse_ele (head, definitions, sad_ele_list):
  for name_and_params in re.finditer(regexDef, definitions):
    ele = ele_struct()
    ele.name = name_and_params.group(1).strip()
    ele.type = head
    params = name_and_params.group(2).strip()

    for param_and_val in re.finditer(regexParam, params):
      ele.param[param_and_val.group(1).strip()] = param_and_val.group(2).strip()

    sad_ele_list[ele.name] = ele

# ------------------------------------------------------------------

def parse_param (head, definitions, sad_param_list):
  definitions = definitions.strip()
  if definitions[0] == '=': definitions = definitions[1:].strip()  # Strip off '='
  sad_param_list[head] = definitions

# ------------------------------------------------------------------

sad_param_names = ['momentum', 'use', 'betax', 'betay']

def parse_directive(directive, sad_ele_list, lat_line_list, sad_param_list):

  directive = directive.strip()  # Remove leading and trailing blanks.
  head, blank, definitions = directive.partition(" ")
  if head == 'ffs': head, blank, definitions = definitions.partition(" ") # Strip off FFS if present

  if head in sad_param_names:
    parse_param (head, definitions, sad_param_list)

  elif head == 'line':
    parse_line(definitions, lat_line_list)

  elif head in sad_names:
    parse_ele(head, definitions, sad_ele_list)

# ------------------------------------------------------------------

def print_help():
  print ('Syntax: sad_to_bmad.py {-open} {-closed} <input-sad-file>')
  sys.exit()

# ------------------------------------------------------------------
# ------------------------------------------------------------------
# ------------------------------------------------------------------

mark_open = False
mark_closed = False

i = 1
while i < len(sys.argv):
  if sys.argv[i] == '-open': 
    mark_open = True
  elif sys.argv[i] == '-closed': 
    mark_closed = True
  elif sys.argv[i][0] == '-':
    print_help()
  else:
    inputfile = sys.argv[i]

  i += 1

if inputfile.find('sad') != -1:
  outputfile = inputfile.replace('sad', 'bmad')
elif inputfile.find('Sad') != -1:
  outputfile = inputfile.replace('Sad', 'Bmad')
elif inputfile.find('SAD') != -1:
  outputfile = inputfile.replace('SAD', 'BMAD')
else:
  outputfile = inputfile + '.bmad'

print ('Input file is ', inputfile)
print ('Output file is ', outputfile)

f_in = open(inputfile, 'r')
f_out = open(outputfile, 'w')

sad_names = ("drift", "bend", "quad", "sext", "oct", "mult", "sol", "cavi", "moni", "line", "beambeam", "apert", "mark")


# ------------------------------------------------------------------
# Read in SAD file line-by-line.  Assemble lines into directives, which are delimited by a ; (colon).
# Call parse_directive whenever an entire directive has been obtained.

sad_ele_list = {}
sad_param_list = {}
lat_line_list = {}

directive = ''
in_comment = False

for line in f_in:
  line = line.strip()              # Remove leading and trailing blanks.
  line = line.lower()              # All letters to lower case.
  line = line.partition('!')[0]    # Remove comments

  if in_comment:
    ix2 = line.find('*)')
    if ix2 == -1: continue    # Next line
    line = line[ix2+2:]
    in_comment = False

  ix = line.find('(*')
  if ix != -1:
    ix2 = line.find('*)')
    if ix2 != -1:
      line = line[:ix] + line[ix2+2:]
    else:
      line = line[:ix]
      in_comment = True

  directive = directive + line + " "
  
  while True:
    ix = directive.find(';')
    if ix == -1: break
    parse_directive(directive[:ix], sad_ele_list, lat_line_list, sad_param_list)
    directive = directive[ix+1:]

# ------------------------------------------------------------------
# Get root lattice line

line0_name = sad_param_list['use']
print ('Using line: ' + line0_name)

if line0_name not in lat_line_list:
  print ('USED LINE NOT FOUND. STOPPING HERE.')
  sys.exit()

this_line = lat_line_list[line0_name]

# For betax and betay translations

ele0_name = this_line.element[0]
if ele0_name in sad_ele_list:
  ele0 = sad_ele_list[ele0_name]
  for key in ele0.param:
    if key in sad_param_names:
      sad_param_list[key] = ele0.param[key]

# ------------------------------------------------------------------
# Translate and write parameters

global_param_translate = {
  'momentum': 'beginning[p0c]',
  'betax':    'beginning[beta_a]',
  'betay':    'beginning[beta_b]'
}

if mark_open: f_out.write ('parameter[geometry] = open\n')
if mark_closed: f_out.write ('parameter[geometry] = closed\n')

for name in sad_param_list:
  if name not in global_param_translate: continue
  f_out.write(global_param_translate[name] + ' = ' + sad_param_list[name] + '\n')

# ------------------------------------------------------------------
# Translate and write element defs

inside_sol = False
bz = '0'

f_out.write ('\n')

for ele_name in this_line.element:

  if not ele_name in sad_ele_list:
    print ('No definition found for element name: ' + ele_name)
    continue

  s_ele = sad_ele_list[ele_name]

  if s_ele.type == 'sol':
    bz = s_ele.param.get('bz')  
    if s_ele.param.get('bound') == '1': inside_sol = not inside_sol
    continue

  if s_ele.printed == True: continue

  b_ele = ele_struct()
  sad_ele_to_bmad (s_ele, b_ele, inside_sol, bz)

  bmadline = b_ele.name + ': ' + b_ele.type
  for param in iter(b_ele.param):
    bmadline += ', ' + param + ' = ' + b_ele.param[param]

  WrapWrite(bmadline)
  b_ele.printed = True

# ------------------------------------------------------------------
# Write lat_line_list line

f_out.write ('\n')

lat_ele = this_line.element
bmadline = this_line.name + ': line = ('

for name in lat_ele: 
  if not ele_name in sad_ele_list: continue  # Already have printed error message
  if sad_ele_list[name].type == 'sol': continue
  bmadline += name + ', '


bmadline = bmadline[:-2] + ')'
WrapWrite(bmadline)


f_out.write ('\n')
f_out.write ('use, ' + line0_name + '\n')

f_in.close()
f_out.close()

