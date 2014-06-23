#!/usr/bin/python

import sys, getopt, re, math
from collections import *

class ele_struct:
  def __init__(self):
    self.name = ''
    self.type = ''
    self.param = dict()
    self.printed = False
    self.instances = 0


class line_item_struct:
  def __init__(self, name = '', reversed = False):
    self.name = name
    self.reversed = reversed

class lat_line_struct:
  def __init__(self):
    self.name = ''
    self.list = []

# ------------------------------------------------------------------

def print_help():
  print (''' \
Syntax: 
  sad_to_bmad.py {-open} {-closed} {-ignore_marker_offsets} {-include <bmad_header_file>} <input-sad-file>

Options:
  -open                   # Put "parameter[geometry] = open" line in the bmad lattice file
  -closed                 # Put "parameter[geometry] = closed" line in the bmad lattice file
  -ignore_marker_offsets  # SAD mark elements which have an offset are translated to a marker element that 
                          #   is superimposed on the lattice. The -ignore_marker_offsets switch means 
                          #   that the offset is ignored and no superposition is done. 
  -include <bmad_header_file>  # Include lines from this file in the Bmad lattice file.

Output:
  If the input SAD file has "sad" in the name, the output bmad lattice file name will substitute "bmad".
  If "sad" is not in the input file name, the output file name will be the input file name with ".bmad" appended.

Please see the NOTES file for known limitations of this translation script.
''')

  sys.exit()

# ------------------------------------------------------------------

def WrapWrite(line):
  MAXLEN = 120
  tab = ''

  while True:
    if len(line) <= MAXLEN:
      f_out.write(tab + line + '\n')
      return

    ix = line[:MAXLEN].rfind(',')
    if ix == -1: ix = line[:MAXLEN].rfind(' ')

    f_out.write(tab + line[:ix+1] + ' &\n')

    tab = '         '
    line = line[ix+1:]

# ------------------------------------------------------------------

ele_type_to_bmad = {
  'drift': 'drift',
  'bend': 'sbend',
  'quad': 'quadrupole',
  'sext': 'sextupole',
  'oct': 'octupole',
  'mult': 'sad_mult',
  'sol': 'marker',
  'cavi': 'rfcavity',
  'moni': 'monitor',
  'mark': 'marker',
  'beambeam': 'marker',
  'apert': 'marker',
}

# Translation rule: Specific translations (of the form 'element:parameter') take
# precedence over generaic translations (of the form 'parameter').

ele_param_translate = {
    'bend:k0': ['k0', ' / @l@'],
    'bend:k1': ['k1', ' / @l@'],
    'quad:k1': ['k1', ' / @l@'],
    'sext:k2': ['k2', ' / @l@'],
    'bend:rotate': ['ref_tilt', ' * -1'],
    'bend:k0': ['g_err', ' / @l@'],
    'l': 'l',
    'radius': 'aperture',
    'offset': 'offset',
    'angle': 'angle',
    'ae1': 'ae1',
    'ae2': 'ae2',
    'e1': ['e1', ' * @angle@'],
    'e2': ['e2', ' * @angle@'],
    'bz': 'bs_field',
    'f1': 'f1',
    'f2': 'f2',
    'eps': 'eps_step_scale',
    'freq': 'rf_frequency',
    'phi': ['phi0', ' / twopi'],
    'volt': 'voltage',
    'bound': 'bound',
    'geo': 'geo',
    'mult:dx': 'x_offset_mult',
    'mult:dy': 'y_offset_mult',
    'dx': 'x_offset',
    'dy': 'y_offset',
    'dz': 'z_offset',
    'sol:chi1': ['x_pitch', ' * -1'],
    'sol:chi2': ['y_pitch', ' * -1'],
    'sol:dz': ['t_offset', ' / -c_light'],
    'mult:chi1': 'x_pitch_mult',
    'mult:chi2': 'y_pitch_mult',
    'chi1': 'x_pitch',
    'chi2': 'y_pitch',
    'chi3': ['tilt', ' * -1'],
    'rotate': ['tilt', ' * -1'], 
    'k0': 'b0', 'k1': 'b1', 'k2': ['b2', ' / factorial(2)'], 'k3': ['b3', ' / factorial(3)'], 
    'k4': ['b4', ' / factorial(4)'], 'k5': ['b5', ' / factorial(5)'], 'k6': ['b6', ' / factorial(6)'], 
    'k7': ['b7', ' / factorial(7)'], 'k8': ['b8', ' / factorial(8)'], 'k9': ['b9', ' / factorial(9)'], 
    'k10': ['b10', ' / factorial(10)'], 'k11': ['b11', ' / factorial(11)'], 'k12': ['b12', ' / factorial(12)'], 
    'k13': ['b13', ' / factorial(13)'], 'k14': ['b14', ' / factorial(14)'], 'k15': ['b15', ' / factorial(15)'], 
    'k16': ['b16', ' / factorial(16)'], 'k17': ['b17', ' / factorial(17)'], 'k18': ['b18', ' / factorial(18)'], 
    'k19': ['b19', ' / factorial(19)'], 'k20': ['b20', ' / factorial(20)'], 'k21': ['b21', ' / factorial(21)'],
    'sk0': 'a0', 'sk1': 'a1', 'sk2': ['a2', ' / factorial(2)'], 'sk3': ['a3', ' / factorial(3)'], 
    'sk4': ['a4', ' / factorial(4)'], 'sk5': ['a5', ' / factorial(5)'], 'sk6': ['a6', ' / factorial(6)'], 
    'sk7': ['a7', ' / factorial(7)'], 'sk8': ['a8', ' / factorial(8)'], 'sk9': ['a9', ' / factorial(9)'], 
    'sk10': ['a10', ' / factorial(10)'], 'sk11': ['a11', ' / factorial(11)'], 'sk12': ['a12', ' / factorial(12)'], 
    'sk13': ['a13', ' / factorial(13)'], 'sk14': ['a14', ' / factorial(14)'], 'sk15': ['a15', ' / factorial(15)'], 
    'sk16': ['a16', ' / factorial(16)'], 'sk17': ['a17', ' / factorial(17)'], 'sk18': ['a18', ' / factorial(18)'], 
    'sk19': ['a19', ' / factorial(19)'], 'sk20': ['a20', ' / factorial(20)'], 'sk21': ['a21 ', ' / factorial(21)'],
}

ignore_sad_param = ['ldev', 'fringe', 'disfrin', 'disrad', 'r1', 'r2', 'r3', 'r4', 'betax', 'betay',
                  'sol:f1', 'sol:bz', 'geo', 'bound', 'index', 'ex', 'ey', 'ax', 'ay', 'bx', 'by', 
                  'epx', 'epy', 'dpx', 'dpy', 'emitx', 'emity', 'dp', 'psix', 'psiy', 'psiz',
                  'sigx', 'sigy', 'sigz', 'slice', 'sturn', 'xangle', 'np']

# ------------------------------------------------------------------

def sad_ele_to_bmad (sad_ele, bmad_ele, inside_sol, bz):

  bmad_ele.name = sad_ele.name

  if not sad_ele.type in ele_type_to_bmad:
    print ('TYPE OF ELEMENT NOT RECOGNIZED: ' + sad_ele.type + '\n' + '    FOR ELEMENT: ' + sad_ele.name)
    return

  bmad_ele.type = ele_type_to_bmad[sad_ele.type]

  # SAD sol with misalignments becomes a Bmad patch

  if sad_ele.type == 'sol':
    if 'dx' in sad_ele.param or 'dy' in sad_ele.param or 'dz' in sad_ele.param or 'chi1' in sad_ele.param or \
                                'chi2' in sad_ele.param or 'chi3' in sad_ele.param or 'rotate' in sad_ele.param:
      if sad_ele.param.get('geo') == '1':
        print ('MISALIGNMENTS IN SOL ELEMENT WITH GEO = 1 NOT YET IMPLEMENTED!')
      else:
        bmad_ele.type = 'patch'

  # Handle case when inside solenoid

  if inside_sol and bmad_ele.type != 'marker' and bmad_ele.type != 'monitor' and bmad_ele.type != 'patch': 
    bmad_ele.type = 'sad_mult'
    bmad_ele.param['bs_field'] = bz

  for sad_param_name in sad_ele.param:
    full_param_name = sad_ele.type + ':' + sad_param_name 
    if sad_param_name in ignore_sad_param: continue
    if full_param_name in ignore_sad_param: continue

    # Use more specific translation first

    if full_param_name in ele_param_translate:
      result = ele_param_translate[full_param_name]
    elif sad_param_name in ele_param_translate:
      result = ele_param_translate[sad_param_name]
    else:
      print ('SAD PARAMETER NOT RECOGNIZED: ' + sad_param_name + '\n' + '    DEFINED IN ELEMENT: ' + sad_ele.name)
      continue

    if type(result) is list:
      bmad_name = result[0]
      value_suffix = result[1]
      if '@' in value_suffix:
        val_parts = value_suffix.split('@')
        if val_parts[1] in sad_ele.param:
          value_suffix = val_parts[0] + sad_ele.param[val_parts[1]] + val_parts[2]
        else:
          print ('SAD ELEMENT: ' + sad_ele.name + '\n' + 
                 '  DOES NOT HAVE A ' + val_parts[1] + ' NEEDED FOR CONVERSION: ' + sad_param_name)
          value_suffix = val_parts[0] + '???' + val_parts[2]
    else:
      bmad_name = result
      value_suffix = ''

    value = sad_ele.param[sad_param_name]
    if 'deg' in value: value = value.replace('deg', '* degrees')
    if 'kev' in value: value = value.replace('kev', '* 1e3')
    if 'mev' in value: value = value.replace('mev', '* 1e6')
    if 'gev' in value: value = value.replace('gev', '* 1e9')

    if sad_param_name == 'radius':   
      if value[0] == '-': value = value[1:]   # Remove negative sign from radius

    bmad_ele.param[bmad_name] = value + value_suffix

  # Correct patch signs

  if bmad_ele.type == 'patch':
    # If entering solenoid 
    if inside_sol:
      if 'x_offset' in bmad_ele.param: bmad_ele.param['x_offset'] = bmad_ele.param['x_offset'] + ' * -1'
      if 'y_offset' in bmad_ele.param: bmad_ele.param['y_offset'] = bmad_ele.param['y_offset'] + ' * -1'
      if 'z_offset' in bmad_ele.param: bmad_ele.param['z_offset'] = bmad_ele.param['z_offset'] + ' * -1'
    # If exiting solenoid 
    else:
      if 'z_offset' in bmad_ele.param: bmad_ele.param['z_offset'] = bmad_ele.param['z_offset'] + ' * -1'

  # Fringe 

  fringe = '0'   # default 
  if 'fringe' in sad_ele.param: fringe = sad_ele.param['fringe']

  disfrin = '0'    # default
  if 'disfrin' in sad_ele.param: disfrin = sad_ele.param['disfrin']

  # Mult and quad fringe

  if sad_ele.type == 'mult' or sad_ele.type == 'quad':
    if fringe == '1':
      bmad_ele.param['fringe_at'] = 'entrance_end'
    elif fringe == '2':
      bmad_ele.param['fringe_at'] = 'exit_end'

    if disfrin == '0':
      if fringe == '0':
        bmad_ele.param['fringe_type'] = 'sad_nonlin_only'
      else:
        bmad_ele.param['fringe_type'] = 'sad_full'

    else:   # disfrin != '0' 
      # fringe == '0' --> Default: bmad_ele.param['fringe_type'] = 'none'
      if fringe != '0':
        bmad_ele.param['fringe_type'] = 'sad_linear'

  # Bend fringe

  elif sad_ele.type == 'bend':
    # fringe == '0' & disfrin != '0' ==> default fringe_type = basic_bend
    if fringe == '0' and disfrin == '0':
      bmad_ele.param['fringe_type'] = 'sad_nonlin_only'
    elif fringe != '0' and disfrin == '0':
      bmad_ele.param['fringe_type'] = 'sad_full'
    elif fringe != '0' and disfrin != '0':
      bmad_ele.param['fringe_type'] = 'sad_linear'

  # Cavi fringe

  elif sad_ele.type == 'cavi':
    if fringe == '1':
      bmad_ele.param['fringe_at'] = 'entrance_end'
    elif fringe == '2':
      bmad_ele.param['fringe_at'] = 'exit_end'

    if disfrin == '0':
      bmad_ele.param['fringe_type'] = 'sad_full'

  # All other fringes  

  elif sad_ele.type == 'deca' or sad_ele.type == 'dodeca' or \
       sad_ele.type == 'oct' or sad_ele.type == 'sext':
    if disfrin == '0':
      bmad_ele.param['fringe_type'] = 'sad_full'

  # bend edge

  if bmad_ele.type == 'sbend':
    if 'ae1' in bmad_ele.param:
      if 'e1' in bmad_ele.param:
        bmad_ele.param['e1'] = bmad_ele.param['e1'] + ' + ' + bmad_ele.param['ae1']
      else:
        bmad_ele.param['e1'] = bmad_ele.param['ae1']
      del bmad_ele.param['ae1']

    if 'ae2' in bmad_ele.param:
      if 'e2' in bmad_ele.param:
        bmad_ele.param['e2'] = bmad_ele.param['e2'] + ' + ' + bmad_ele.param['ae2']
      else:
        bmad_ele.param['e2'] = bmad_ele.param['ae2']
      del bmad_ele.param['ae2']

# ------------------------------------------------------------------

def parse_line(rest_of_line, lat_line_list):
  sad_line = lat_line_struct()
  line_name, e_sign, line_list = rest_of_line.partition("=")
  line_name = line_name.strip()
  line_list = line_list.strip()
  line_list = line_list[1:-1].strip()    # Remove parenteses from "(...)"
  
  sad_line.name = line_name
  for elename in re.split('\s+', line_list):
    if elename[0] == '-':
      line_item = line_item_struct(elename[1:], True)
    else:
      line_item = line_item_struct(elename)
    sad_line.list.append(line_item)

  lat_line_list[line_name] = sad_line

# ------------------------------------------------------------------

regexDef = r"(\S*?)\s*=\s*\((.*?)\)"
regexParam = r"(\S*)\s*=\s*(\S*\s*(deg)?)"

def parse_ele (head, rest_of_line, sad_ele_list):
  for name_and_params in re.finditer(regexDef, rest_of_line):
    ele = ele_struct()
    ele.name = name_and_params.group(1).strip()
    ele.type = head
    params = name_and_params.group(2).strip()

    for param_and_val in re.finditer(regexParam, params):
      value = param_and_val.group(2).strip()
      if value[-1] == ',': value = value[:-1]
      ele.param[param_and_val.group(1).strip()] = value

    # Ignore MARK element offsets?

    if ignore_marker_offsets and ele.type == 'mark' and 'offset' in ele.param: del ele.param['offset']

    sad_ele_list[ele.name] = ele


# ------------------------------------------------------------------

def parse_param (head, rest_of_line, sad_param_list):
  rest_of_line = rest_of_line.strip()
  if len(rest_of_line) > 0 and rest_of_line[0] == '=': rest_of_line = rest_of_line[1:].strip()  # Strip off '='
  sad_param_list[head] = rest_of_line

# ------------------------------------------------------------------
# Rule: After a "calc" command, the only thing left to parse is a possible "initialorbit" setting.

sad_param_names = ['momentum', 'use', 'betax', 'betay', 'nocod']

def parse_directive(directive, sad_ele_list, lat_line_list, sad_param_list):

  global calc_command_found

  directive = directive.strip()  # Remove leading and trailing blanks.
  head, blank, rest_of_line = directive.partition(" ")
  if head == 'ffs': head, blank, rest_of_line = rest_of_line.partition(" ") # Strip off FFS if present

  if head in sad_param_names:
    if calc_command_found: return
    parse_param (head, rest_of_line, sad_param_list)

  elif head == 'line':
    if calc_command_found: return
    parse_line(rest_of_line, lat_line_list)

  elif head in sad_ele_type_names:
    if calc_command_found: return
    parse_ele(head, rest_of_line, sad_ele_list)

  elif head == 'calc':
    calc_command_found = True

  elif 'initialorbit' in rest_of_line:
    line = rest_of_line.partition('initialorbit')[2]
    line = line.partition('{')[2].partition('}')[0]
    orbit = line.replace(',', ' ').split()
    sad_param_list['x_orb']  = orbit[0]
    sad_param_list['px_orb'] = orbit[1]
    sad_param_list['y_orb']  = orbit[2]
    sad_param_list['py_orb'] = orbit[3]
    sad_param_list['z_orb']  = orbit[4]
    sad_param_list['pz_orb'] = orbit[5]

  elif '=' in directive:
    head, blank, rest_of_line = directive.partition('=')
    head = head.strip()
    if head == 'dxi':  sad_param_list['x_orb']  = rest_of_line
    if head == 'dpxi': sad_param_list['px_orb'] = rest_of_line
    if head == 'dyi':  sad_param_list['y_orb']  = rest_of_line
    if head == 'dpyi': sad_param_list['py_orb'] = rest_of_line
    if head == 'dzi':  sad_param_list['z_orb']  = rest_of_line
    if head == 'ddpi': sad_param_list['pz_orb'] = rest_of_line


# ------------------------------------------------------------------
# ------------------------------------------------------------------
# ------------------------------------------------------------------

header_file = ''
mark_open = False
mark_closed = False
inputfile = None
sad_sol_to_marker = True           # False -> No element in bmad lattice.
ignore_marker_offsets = False

i = 1
while i < len(sys.argv):
  n = len(sys.argv[i])
  if sys.argv[i] == '-open': 
    mark_open = True
  elif sys.argv[i] == '-closed': 
    mark_closed = True
  elif n > 1 and '-ignore_marker_offsets'[:n] == sys.argv[i]:
    ignore_marker_offsets = True
  elif sys.argv[i] == '-include':
    i += 1
    header_file = sys.argv[i]
  elif sys.argv[i][0] == '-':
    print_help()
  else:
    inputfile = sys.argv[i]

  i += 1

#

if inputfile == None:
  print_help()

if inputfile.find('sad') != -1:
  outputfile = inputfile.replace('sad', 'bmad')
elif inputfile.find('Sad') != -1:
  outputfile = inputfile.replace('Sad', 'Bmad')
elif inputfile.find('SAD') != -1:
  outputfile = inputfile.replace('SAD', 'BMAD')
else:
  outputfile = inputfile + '.bmad'

print ('Input file is:  ' + inputfile)
print ('Output file is: ' + outputfile)

f_in = open(inputfile, 'r')
f_out = open(outputfile, 'w')

sad_ele_type_names = ("drift", "bend", "quad", "sext", "oct", "mult", "sol", "cavi", "moni", "line", "beambeam", "apert", "mark")


# ------------------------------------------------------------------
# Read in SAD file line-by-line.  Assemble lines into directives, which are delimited by a ; (colon).
# Call parse_directive whenever an entire directive has been obtained.

sad_ele_list = OrderedDict()
sad_param_list = OrderedDict()
lat_line_list = OrderedDict()
calc_command_found = False

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

sad_line = lat_line_list[line0_name]

# For betax and betay translations

ele0_name = sad_line.list[0].name
if ele0_name in sad_ele_list:
  ele0 = sad_ele_list[ele0_name]
  for key in ele0.param:
    if key in sad_param_names:
      sad_param_list[key] = ele0.param[key]

# ------------------------------------------------------------------
# Header

if header_file != '':
  h_file = open(header_file, 'r')
  for line in h_file:
    f_out.write (line)

# ------------------------------------------------------------------
# Translate and write parameters

global_param_translate = {
  'momentum': 'beginning[p0c]',
  'betax':    'beginning[beta_a]',
  'betay':    'beginning[beta_b]',
  'nocod':    'parameter[geometry] = open',
  'x_orb':    'beam_start[x]',
  'px_orb':   'beam_start[px]',
  'y_orb':    'beam_start[y]',
  'py_orb':   'beam_start[py]',
  'z_orb':    'beam_start[z]',
  'pz_orb':   'beam_start[pz]',
}

if mark_open: f_out.write ('parameter[geometry] = open\n')
if mark_closed: f_out.write ('parameter[geometry] = closed\n')

for name in sad_param_list:
  if name not in global_param_translate: continue
  if sad_param_list[name] == '':
    f_out.write(global_param_translate[name] + '\n')
  else:
    f_out.write(global_param_translate[name] + ' = ' + sad_param_list[name] + '\n')

# ------------------------------------------------------------------
# Translate and write element defs

inside_sol = False
bz = '0'
bmad_line = []

f_out.write ('\n')

for ix_s_ele, sad_line_list in enumerate(sad_line.list):

  ele_name = sad_line_list.name

  if not ele_name in sad_ele_list:
    print ('No definition found for element name: ' + ele_name)
    continue

  s_ele = sad_ele_list[ele_name]

  # sol element

  if s_ele.type == 'sol':
    bz = s_ele.param.get('bz')  
    if s_ele.param.get('bound') == '1': inside_sol = not inside_sol

    # misalignment with geo = 1 or patch to orbit if not.

    if 'dx' in s_ele.param or 'dy' in s_ele.param or 'dz' in s_ele.param or 'chi1' in s_ele.param or \
                                'chi2' in s_ele.param or 'chi3' in s_ele.param or 'rotate' in s_ele.param:
      if s_ele.param.get('geo') == '1': 
        print ('MISALIGNMENTS IN SOL ELEMENT WITH GEO = 1 NOT YET IMPLEMENTED!')
        if not sad_sol_to_marker: continue
    
    else:
      if not sad_sol_to_marker: continue

  # A MARK element with an offset gets translated to a marker superimpsed with respect to a null_ele

  if s_ele.type == 'mark' and 'offset' in s_ele.param:
    null_ele_name = 'null_' + s_ele.name + '#' + str(ix_s_ele)  # Guaranteed unique
    bmad_line.append (null_ele_name)                            # Put null_ele in the line
    WrapWrite(null_ele_name + ': null_ele')                     # Define the null_ele

    # Now define the marker element
    sad_offset = float(s_ele.param['offset'])
    int_off = int(math.floor(sad_offset))
    frac_off = sad_offset - int_off

    offset = 0
    direc = 1
    if int_off < 0: direc = -1
    for ix in range(0, int_off, direc):
      ss_ele = sad_ele_list[sad_line.list[ix_s_ele+ix].name]
      if 'l' in ss_ele.param: offset += direc * float(ss_ele.param['l'])
    ss_ele = sad_ele_list[sad_line.list[ix_s_ele+int_off].name]
    if 'l' in ss_ele.param: offset += frac_off * float(ss_ele.param['l'])

    if s_ele.instances == 0:
      suffix = ''
    else:
      suffix = '.' + str(s_ele.instances)
    bmad_ele_def = s_ele.name + suffix + ': marker, superimpose, ref = ' + null_ele_name + ', offset = ' + str(offset)
    WrapWrite(bmad_ele_def)
    s_ele.printed = True
    s_ele.instances += 1
    continue

  # Regular element not getting superimposed

  bmad_line.append(s_ele.name)

  if s_ele.printed == True: continue

  b_ele = ele_struct()
  sad_ele_to_bmad (s_ele, b_ele, inside_sol, bz)

  bmad_ele_def = b_ele.name + ': ' + b_ele.type
  for param in iter(b_ele.param):
    bmad_ele_def += ', ' + param + ' = ' + b_ele.param[param]

  WrapWrite(bmad_ele_def)
  s_ele.printed = True

# ------------------------------------------------------------------
# Write lat_line_list line

f_out.write ('\n')

bmad_line_str = sad_line.name + ': line = ('

for name in bmad_line: 
  bmad_line_str += name + ', '


bmad_line_str = bmad_line_str[:-2] + ')'
WrapWrite(bmad_line_str)


f_out.write ('\n')
f_out.write ('use, ' + line0_name + '\n')

f_in.close()
f_out.close()

