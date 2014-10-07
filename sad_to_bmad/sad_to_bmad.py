#!/usr/bin/python

import sys, getopt, re, math
from collections import *
import time

start_time = time.time()
ix_null = 0

class ele_struct:
  def __init__(self):
    self.name = ''
    self.type = ''
    self.param = dict()
    self.printed = False
    self.instances = 0

class line_item_struct:
  def __init__(self, name = '', sign = '', multiplyer = '1'):
    self.name = name
    self.sign = ''
    self.multiplyer = multiplyer

class lat_line_struct:
  def __init__(self):
    self.name = ''
    self.printed = False
    self.list = []

class sad_info_struct:
  def __init__(self):
    self.lat_line_list = OrderedDict()
    self.ele_list = OrderedDict()
    self.param_list = OrderedDict()
    self.var_list = OrderedDict()

#------------------------------------------------------------------
#------------------------------------------------------------------

def add_units (line):
  if 'deg' in line: line = re.sub(' deg$', ' * degrees', line)
  if 'kev' in line: line = re.sub(' kev$', ' * 1e3', line)
  if 'mev' in line: line = re.sub(' mev$', ' * 1e6', line)
  if 'gev' in line: line = re.sub(' gev$', ' * 1e9', line)
  return line

#------------------------------------------------------------------
#------------------------------------------------------------------

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

#-------------------------------------------------------------------
#------------------------------------------------------------------

def WrapWrite(line):
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

#------------------------------------------------------------------
# Adds parenteses around expressions with '+' or '-' operators.
# Otherwise just returns the expression.
# Eg: '7+3' ->  '(7+3)'
#     '7*3  ->  '7*3'

def add_parens (str):
  for ix in range(1, len(str)):
    if str[ix] != '+' and str[ix] != '-': continue
    if str[ix-1] == 'e' or str[ix-1] == 'E': continue  # '+' in '3.0e+7' is not an operator
    return '(' + str + ')'

  # No +/- op found so just return the expression.
  return str    

#------------------------------------------------------------------
#------------------------------------------------------------------

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
  'apert': 'rcollimator',
}

# Translation rule: Specific translations (of the form 'element:parameter') take
# precedence over generaic translations (of the form 'parameter').

# Something like ['k0', ' / @l@'] is translated to "<k0> / <l>" where <k0> and <l> are the values
# of the k0 and l parameters.

ele_param_translate = {
    'apert:dx': 'x_offset',
    'apert:dy': 'y_offset',
    'bend:k0': ['g_err', ' / @l@'],
    'bend:k1': ['k1', ' / @l@'],
    'bend:fb1': ['hgap', '/6, fint = 0.5'],
    'bend:fb2': ['hgapx', '/6, fintx = 0.5'],
    'quad:k1': ['k1', ' / @l@'],
    'sext:k2': ['k2', ' / @l@'],
    'oct:k3': ['k3', ' / @l@'],
    'bend:rotate': ['ref_tilt', ' * -1'],
    'l': 'l',
    'radius': 'aperture',
    'offset': 'offset',
    'angle': 'angle',
    'ae1': 'ae1',
    'ae2': 'ae2',
    'e1': ['e1', ' * @angle@'],
    'e2': ['e2', ' * @angle@'],
    'bz': 'bs_field',
    'sol:f1': 'f1',
    'f1': 'fq1',       # For mult and quad elements
    'f2': 'fq2',
    'eps': 'eps_step_scale',
    'freq': 'rf_frequency',
    'phi': ['phi0', ' / twopi'],
    'dphi': ['phi0_err', ' / twopi'],
    'volt': 'voltage',
    'bound': 'bound',
    'harm': 'harmon',
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
                  'sigx', 'sigy', 'sigz', 'slice', 'sturn', 'xangle', 'np', 'ddp', 
                  'pex', 'pepx', 'pey', 'pepy', 'trx', 'try', 'leng', 'ax', 'ay', 'dx1', 'dx2', 'dy1', 'dy2']

sad_reversed_params = {
      'ae1': 'ae2',
      'ae2': 'ae1',
      'e1': 'e2',
      'e2': 'e1',
      'fb1': 'fb2',
      'fb2': 'fb1'
}
  
#------------------------------------------------------------------
#------------------------------------------------------------------

def output_line (sad_line, sad_info, inside_sol, bz):

  global ix_null

  f_out.write ('\n')

  bmad_line = []
  sad_line.printed = True

  for ix_s_ele, sad_ele in enumerate(sad_line.list):

    ele_name = sad_ele.name

    # If the line element is itself a line then print this line info.

    if ele_name in sad_info.lat_line_list:
      if not sad_info.lat_line_list[ele_name].printed: 
        output_line(sad_info.lat_line_list[ele_name], sad_info, inside_sol, bz)
      bmad_line.append(sad_ele)
      continue

    if not ele_name in sad_info.ele_list:
      print ('No definition found for element name: ' + ele_name)
      continue

    s_ele = sad_info.ele_list[ele_name]

    # sol element

    if s_ele.type == 'sol':
      bz = s_ele.param.get('bz', '0')
      try:
        b_z = float(bz)
        if b_z == 0: bz = '0'    # EG convert '0.0' to '0'
      except ValueError:
        pass
      if s_ele.param.get('bound') == '1': inside_sol = not inside_sol

    # A MARK element with an offset gets translated to a marker superimpsed with respect to a null_ele

    if s_ele.type == 'mark' and 'offset' in s_ele.param:
      if ignore_marker_offsets:
        del s_ele.param['offset']
      else:
        ix_null += 1
        null_ele_name = 'null_' + s_ele.name + '#' + str(ix_null)   # Guaranteed unique
        bmad_line.append (line_item_struct(null_ele_name))          # Put null_ele in the line
        WrapWrite(null_ele_name + ': null_ele')                     # Define the null_ele
  
        # Now define the marker element
        sad_offset = float(s_ele.param['offset'])
        int_off = int(math.floor(sad_offset))
        frac_off = sad_offset - int_off
  
        offset = 0
        direc = 1
        if int_off < 0: direc = -1
        for ix in range(0, int_off, direc):
          this_name = sad_line.list[ix_s_ele+ix].name
          if this_name in sad_info.ele_list:
            ss_ele = sad_info.ele_list[this_name]
            if 'l' in ss_ele.param: offset += direc * eval(ss_ele.param['l'])
          else:   # Must be a line
            for sub_ele in sad_info.lat_line_list[this_name].list:
              ss_ele = sad_info.ele_list[sub_ele.name]
              if 'l' in ss_ele.param: offset += direc * eval(ss_ele.param['l'])
  
        ss_ele = sad_info.ele_list[sad_line.list[ix_s_ele+int_off].name]
        if 'l' in ss_ele.param: offset += frac_off * eval(ss_ele.param['l'])
  
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

    bmad_line.append(sad_ele)

    if s_ele.printed == True: continue

    b_ele = ele_struct()
    sad_ele_to_bmad (s_ele, b_ele, inside_sol, bz, sad_ele.sign == '-')

    bmad_ele_def = b_ele.name + ': ' + b_ele.type
    for param in iter(b_ele.param):
      bmad_ele_def += ', ' + param + ' = ' + b_ele.param[param]

    WrapWrite(bmad_ele_def)
    s_ele.printed = True

  #---------------------------------------------
  # Write line

  f_out.write ('\n')

  bmad_line_str = sad_line.name + ': line = ('

  for ele in bmad_line:
    if ele.multiplyer == '1':
      bmad_line_str += ele.sign + ele.name + ', '
    else:
      bmad_line_str += ele.sign + ele.multiplyer + '*' + ele.name + ', '

  bmad_line_str = bmad_line_str[:-2] + ')'
  WrapWrite(bmad_line_str)

#------------------------------------------------------------------
#------------------------------------------------------------------

def sad_ele_to_bmad (sad_ele, bmad_ele, inside_sol, bz, reversed):

  bmad_ele.name = sad_ele.name

  if not sad_ele.type in ele_type_to_bmad:
    print ('TYPE OF ELEMENT NOT RECOGNIZED: ' + sad_ele.type + '\n' + '    FOR ELEMENT: ' + sad_ele.name)
    return

  bmad_ele.type = ele_type_to_bmad[sad_ele.type]

  # SAD sol with misalignments becomes a Bmad patch

  if sad_ele.type == 'sol':
    if len(set(['dx', 'dy', 'dz', 'chi1', 'chi2', 'chi3', 'rotate']).intersection(sad_ele.param)) > 0:
      if sad_ele.param.get('geo') == '1':
        print ('MISALIGNMENTS IN SOL ELEMENT '+ sad_ele.name + ' WITH GEO = 1 NOT YET IMPLEMENTED! WILL BE IGNORED!')
        print ('  IF MISALIGNMENTS ARE SMALL, OR BZ = 0 THROUGH THE SOLENOID, THIS IS NOT A PROBLEM:')
        for param in sad_ele.param:
          if param in ['dx', 'dy', 'dz', 'chi1', 'chi2', 'chi3', 'rotate']:
            print ('  ' + param + ' = ' + sad_ele.param[param])
      else:
        bmad_ele.type = 'patch'

  # Handle case when inside solenoid

  if inside_sol and bmad_ele.type != 'marker' and bmad_ele.type != 'monitor' and \
                    bmad_ele.type != 'patch' and bz != '0': 
    bmad_ele.type = 'sad_mult'
    bmad_ele.param['bs_field'] = bz

  # convert apert

  if sad_ele.type == 'apert':
    ax = float(sad_ele.param.get('ax', '0'))
    ay = float(sad_ele.param.get('ay', '0'))
    dx1 = float(sad_ele.param.get('dx1', '0'))
    dx2 = float(sad_ele.param.get('dx2', '0'))
    dy1 = float(sad_ele.param.get('dy1', '0'))
    dy2 = float(sad_ele.param.get('dy2', '0'))

    if ax == 0 and ay == 0:
      bmad_ele.param['x1_limit'] = str(-min(dx1, dx2))
      bmad_ele.param['x2_limit'] = str(max(dx1, dx2))
      bmad_ele.param['y1_limit'] = str(-min(dy1, dy2))
      bmad_ele.param['y2_limit'] = str(max(dy1, dy2))
    elif dx1 == dx2 and dy1 == dy2:
      bmad_ele.type = 'ecollimator'
      bmad_ele.param['x_limit'] = str(ax)
      bmad_ele.param['y_limit'] = str(ay)
    else:
      print ('Combined collimators not yet implemented. Please contact David Sagan')

  # For a bend, f1 mut be added to fb1 and fb2

  if sad_ele.type == 'bend' and 'f1' in sad_ele.param:
    if 'fb1' in sad_ele.param:
      sad_ele.param['fb1'] = sad_ele.param['fb1'] + sad_ele.param['f1']
    else:
      sad_ele.param['fb1'] = sad_ele.param['f1']

    if 'fb2' in sad_ele.param:
      sad_ele.param['fb2'] = sad_ele.param['fb2'] + sad_ele.param['f1']
    else:
      sad_ele.param['fb2'] = sad_ele.param['f1']

    del sad_ele.param['f1']

    # If fb1 == fb2 then don't need fb2 (Bmad will take them as equal if fb2 is not present).

    if 'fb1' in sad_ele.param:
      if 'fb2' in sad_ele.param:
        if sad_ele.param['fb1'] == sad_ele.param['fb2']: del sad_ele.param['fb2']
      else:
        sad_ele.param['fb2'] = '0'

  # Loop over all parameters

  for sad_param_name in sad_ele.param:

    full_param_name = sad_ele.type + ':' + sad_param_name 
    if sad_param_name in ignore_sad_param: continue
    if full_param_name in ignore_sad_param: continue

    value = sad_ele.param[sad_param_name]

    # For reversed elements, ae1 -> ae2, etc.

    if reversed and sad_param_name in sad_reversed_params: 
      sad_param_name = sad_reversed_params[sad_param_name]
      full_param_name = sad_ele.type + ':' + sad_param_name 

    # Use more specific translation first

    if full_param_name in ele_param_translate:
      result = ele_param_translate[full_param_name]
    elif sad_param_name in ele_param_translate:
      result = ele_param_translate[sad_param_name]
    else:
      print ('SAD PARAMETER NOT RECOGNIZED: ' + sad_param_name + '\n' + '    DEFINED IN ELEMENT: ' + sad_ele.name)
      continue

    # 

    if sad_param_name == 'k1' and sad_ele.type == 'quad' and 'l' not in sad_ele.param:
      bmad_ele.type = 'multipole'
      bmad_name = 'k1l'
      value_suffix = ''

    elif sad_param_name == 'k2' and sad_ele.type == 'sext' and 'l' not in sad_ele.param:
      bmad_ele.type = 'multipole'
      bmad_name = 'k2l'
      value_suffix = ''

    elif sad_param_name == 'k3' and sad_ele.type == 'oct' and 'l' not in sad_ele.param:
      bmad_ele.type = 'multipole'
      bmad_name = 'k3l'
      value_suffix = ''

    elif sad_param_name == 'rotate' and sad_ele.type == 'quad' and 'l' not in sad_ele.param:
      bmad_ele.type = 'multipole'
      bmad_name = 't1'
      value_suffix = ''

    elif sad_param_name == 'rotate' and sad_ele.type == 'sext' and 'l' not in sad_ele.param:
      bmad_ele.type = 'multipole'
      bmad_name = 't2'
      value_suffix = ''

    elif sad_param_name == 'rotate' and sad_ele.type == 'oct' and 'l' not in sad_ele.param:
      bmad_ele.type = 'multipole'
      bmad_name = 't3'
      value_suffix = ''

    elif type(result) is list:   # EG: result = ['k1', ' / @l@']
      bmad_name = result[0]
      value_suffix = result[1]
      value = add_parens(value)
      if '@' in value_suffix:
        val_parts = value_suffix.split('@')
        if val_parts[1] in sad_ele.param:
          value_suffix = val_parts[0] + add_parens(sad_ele.param[val_parts[1]]) + val_parts[2]
        elif '/' not in val_parts[0]:      # Assume zero
          value_suffix = val_parts[0] + '0' + val_parts[2]
        else:
          print ('SAD ELEMENT: ' + sad_ele.name + '\n' + 
                 '  DOES NOT HAVE A ' + val_parts[1] + ' NEEDED FOR CONVERSION: ' + sad_param_name)
          value_suffix = val_parts[0] + '???' + val_parts[2]
    else:
      bmad_name = result
      value_suffix = ''

    #

    if reversed:
      if sad_param_name == 'offset': value = '1 - ' + value
      if sad_param_name in ['chi1', 'chi2', 'dz']: 
        if value[0] == '-': 
          value = value[1:]
        else:
          value = '-' + value 

    if sad_param_name == 'radius':   
      if value[0] == '-': value = value[1:]   # Remove negative sign from radius

    bmad_ele.param[bmad_name] = value + value_suffix

  # End of parameter loop
  #-----------------------

  # f1 -> fq1 translation for mult and quad elements involves some work

  if 'fq1' in bmad_ele.param:
    if bmad_ele.param['fq1'][0:1] == '-':
      bmad_ele.param['fq1'] = '(' + bmad_ele.param['fq1'][1:] + ')^2 / 24'
    else:
      bmad_ele.param['fq1'] = '-(' + bmad_ele.param['fq1'] + ')^2 / 24'

  # If a SAD mult element has acceleration then it become an rfcavity or lcavity.

  if sad_ele.type == 'mult' and 'volt' in sad_ele.param: bmad_ele.type = 'rfcavity'

  # If the SAD cavi has a nonzero phi then the reference particle's energy is changing
  # and so the corresponding Bmad element must be an lcavity.
  # Also remember that an lcavity has a differnt phase convention.

  if bmad_ele.type == 'rfcavity':
    if 'phi0' in bmad_ele.param and 'harmon' not in bmad_ele.param: 
      bmad_ele.type = 'lcavity'
      bmad_ele.param['phi0'] = '0.25 - ' + add_parens(bmad_ele.param['phi0'])
      if 'phi0_err' in bmad_ele.param: bmad_ele.param['phi0_err'] = '-' + add_parens(bmad_ele.param['phi0_err'])

    elif 'phi0_err' in bmad_ele.param:
      if 'phi0' in bmad_ele.param:
        bmad_ele.param['phi0'] = bmad_ele.param['phi0'] + ' + ' + bmad_ele.param['phi0_err']
      else:
        bmad_ele.param['phi0'] = bmad_ele.param['phi0_err']
      del bmad_ele.param['phi0_err']

  # Correct patch signs

  if bmad_ele.type == 'patch':
    # If entering solenoid 
    if inside_sol:
      if 'x_offset' in bmad_ele.param: bmad_ele.param['x_offset'] = add_parens(bmad_ele.param['x_offset']) + ' * -1'
      if 'y_offset' in bmad_ele.param: bmad_ele.param['y_offset'] = add_parens(bmad_ele.param['y_offset']) + ' * -1'
      if 'z_offset' in bmad_ele.param: bmad_ele.param['z_offset'] = add_parens(bmad_ele.param['z_offset']) + ' * -1'
    # If exiting solenoid 
    else:
      if 'z_offset' in bmad_ele.param: bmad_ele.param['z_offset'] = add_parens(bmad_ele.param['z_offset']) + ' * -1'

  # Fringe 

  fringe = '0'   # default 
  if 'fringe' in sad_ele.param: fringe = sad_ele.param['fringe']

  if reversed:
    if fringe == 1:
      fringe = 2
    elif fringe == 2:
      fringe = 1

  disfrin = '0'    # default
  if 'disfrin' in sad_ele.param: disfrin = sad_ele.param['disfrin']

  # Bmad multipoles [created since sad element had zero length], do not have a fringe

  if bmad_ele.type == 'multipole':
    pass    # No fringe

  # Mult and quad fringe

  elif sad_ele.type == 'mult' or sad_ele.type == 'quad':
    if fringe == '1':
      bmad_ele.param['fringe_at'] = 'entrance_end'
    elif fringe == '2':
      bmad_ele.param['fringe_at'] = 'exit_end'

    if disfrin == '0':
      if fringe == '0':
        bmad_ele.param['fringe_type'] = 'hard_edge_only'
      else:
        bmad_ele.param['fringe_type'] = 'full'

    else:   # disfrin != '0' 
      # fringe == '0' --> Default: bmad_ele.param['fringe_type'] = 'none'
      if fringe != '0':
        bmad_ele.param['fringe_type'] = 'soft_edge_only'

  # Bend fringe

  elif sad_ele.type == 'bend':
    # fringe == '0' & disfrin != '0' ==> default fringe_type = basic_bend
    if fringe == '0' and disfrin == '0':
      bmad_ele.param['fringe_type'] = 'hard_edge_only'
    elif fringe != '0' and disfrin == '0':
      bmad_ele.param['fringe_type'] = 'sad_full'
    elif fringe != '0' and disfrin != '0':
      bmad_ele.param['fringe_type'] = 'sad_soft_edge_only'

  # Cavi fringe

  elif sad_ele.type == 'cavi':
    if fringe == '1':
      bmad_ele.param['fringe_at'] = 'entrance_end'
    elif fringe == '2':
      bmad_ele.param['fringe_at'] = 'exit_end'

    if disfrin == '0':
      bmad_ele.param['fringe_type'] = 'full'

  # All other fringes  

  elif sad_ele.type == 'deca' or sad_ele.type == 'dodeca' or \
       sad_ele.type == 'oct' or sad_ele.type == 'sext':
    if disfrin == '0':
      bmad_ele.param['fringe_type'] = 'full'

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

#------------------------------------------------------------------
#------------------------------------------------------------------

def parse_line(rest_of_line, sad_info):

  parse_status = 'init'

  for token in re.split('\s+|(=|\)|\(|\+|\*|-|,)', rest_of_line):
    if token == '': continue
    if token == None: continue
    if token == ',': continue

    if parse_status == 'init':
      if token in '+-*=()':
        print ('ERROR PARSING LINE: ' + rest_of_line)
        sys.exit()
      sad_line = lat_line_struct()
      sad_line.name = token
      parse_status = 'got line name'
      continue

    if parse_status == 'got line name':
      if token != '=':
        print ('ERROR PARSING LINE: ' + rest_of_line)
        sys.exit()
      parse_status = 'got ='
      continue

    if parse_status == 'got =':
      if token != '(':
        print ('ERROR PARSING LINE: ' + rest_of_line)
        sys.exit()
      parse_status = 'got ('
      sign = ''
      multiplyer = '1'  
      continue

    #

    if token == '-':
      sign = '-'

    elif token == '*':
      continue

    elif token == ')':
      sad_info.lat_line_list[sad_line.name] = sad_line
      sign = ''
      multiplyer = '1'  
      parse_status = 'init'

    elif token == '(':
      print ('ERROR PARSING LINE: ' + rest_of_line)
      sys.exit()

    elif token == '+':
      continue

    elif token.isdigit():
      multiplyer = token

    else:
      line_item = line_item_struct(token, sign, multiplyer)
      if line_item.name == 'end': line_item.name = 'end_ele'  # 'end' is a reserved name
      sad_line.list.append(line_item)
      sign = ''
      multiplyer = '1'

  #

  if parse_status != 'init':
    print ('ERROR PARSING LINE: ' + rest_of_line)
    sys.exit()

#------------------------------------------------------------------
#------------------------------------------------------------------

regexDef = r"(\S*?)\s*=\s*\((.*?)\)"
regexParam = r"(\S*)\s*=\s*(\S*\s*(deg)?)"

def parse_ele (head, rest_of_line, sad_info):

  # Since the length of rest_of_line may be huge (~ 1M characters for KEKB mult elements.), break
  # this string into chunks to make the script *much* faster.

  len_str =  9000
  len_min =  3000

  rest_of_line = rest_of_line.strip()

  if len(rest_of_line) > len_str:
    line = rest_of_line[:len_str]
    line_saved = rest_of_line[len_str:]
  else:
    line = rest_of_line
    line_saved = ''


  while True:

    if len(line) < len_min and len(line_saved) != 0:
      if len(line_saved) > len_str:
        line = line + line_saved[:len_str]
        line_saved = line_saved[len_str:]
      else:
        line = (line + line_saved)
        line_saved = ''

    if len(line) == 0: break

    for ix in range(len(line)):
      if line[ix] == ' ': continue

      if ix > len(line) - 3:
        print ('MALFORMED ELEMENT DEFINITION: ' + rest_of_line)
        sys.exit()

      if line[ix] == '=':
        ele = ele_struct()
        ele.type = head
        ele.name = line[:ix].strip()
        ## print ('Name: "' + ele.name + '"')
        line = line[ix+1:].lstrip()
        break

    if line[0] != '(':
      print ('MALFORMED ELEMENT DEFINITION: ' + rest_of_line)
      sys.exit()
    line = line[1:].lstrip()

    # parameter loop

    # First look for first parameter name

    for ix in range(len(line)):
      if line[ix] == '=' or line[ix] == ')':
        param_name = line[:ix].strip()
        delim = line[ix]
        line = line[ix+1:]
        break
    
    if delim != ')':

      # get parameter value and name of next parameter
      # Find next '=' or ')'

      ix0 = 0
      n_parens = 0
      ## print ('Line: "' + line + '"')

      for ix in range(len(line)):
        if line[ix] == '(':
          n_parens += 1
          continue

        if line[ix] == ')':
          if n_parens == 0:
            ## print ('Found ): "' + line[ix0:ix] + '"')
            ele.param[param_name] = add_units(line[ix0:ix].strip())
            line = line[ix+1:]
            break

          n_parens -= 1
          continue

        if line[ix] == '=':
          ## print ('Found =: "' + line[ix0:ix] + '" + "' + line[ix:] + '"')
          sub_str = line[ix0:ix].strip()
          j = max(sub_str.rfind(' '), sub_str.rfind(','))
          ## print ('sub_str: "' + sub_str + '"  ' + str(j))
          ele.param[param_name] = add_units(sub_str[:j].strip())
          param_name = sub_str[j+1:].strip()
          delim = line[ix]
          ix0 = ix + 1

    # Put element in list

    if ele.name == 'end': ele.name = 'end_ele'  # 'end' is a reserved name
    sad_info.ele_list[ele.name] = ele

#------------------------------------------------------------------
#------------------------------------------------------------------

def parse_param (head, line, sad_info):
  line = line.strip()
  if len(line) > 0 and line[0] == '=': line = line[1:].strip()  # Strip off '='
  sad_info.param_list[head] = add_units(line)

#------------------------------------------------------------------
#------------------------------------------------------------------
# Rule: After a "calc" command, the only thing left to parse is a possible "initialorbit" setting.

sad_param_names = ['momentum', 'use', 'betax', 'betay', 'nocod']

global_param_translate = {
  'momentum': 'beginning[p0c]',
  'betax':    'beginning[beta_a]',
  'betay':    'beginning[beta_b]',
  'nocod':    'parameter[geometry] = open',
  'cod':      'parameter[geometry] = closed',
  'x_orb':    'beam_start[x]',
  'px_orb':   'beam_start[px]',
  'y_orb':    'beam_start[y]',
  'py_orb':   'beam_start[py]',
  'z_orb':    'beam_start[z]',
  'pz_orb':   'beam_start[pz]',
  'dxi':      'beam_start[x]',
  'dpxi':     'beam_start[px]',
  'dyi':      'beam_start[y]',
  'dpyi':     'beam_start[py]',
  'dzi':      'beam_start[z]',
  'ddpi':     'beam_start[pz]',
  'npara':    ''
}

def parse_directive(directive, sad_info):

  global calc_command_found

  directive = directive.strip()  # Remove leading and trailing blanks.
  head, blank, rest_of_line = directive.partition(" ")
  if head == 'ffs': head, blank, rest_of_line = rest_of_line.partition(" ") # Strip off FFS if present

  if ',' in head or '=' in head:
    if ',' in head: 
      p1, delim, p2 = head.partition(',')
    else:
      p1, delim, p2 = head.partition('=')

    head = p1
    rest_of_line = delim + p2 + rest_of_line 

  if head in global_param_translate or head == 'use':
    if calc_command_found: return
    parse_param (head, rest_of_line, sad_info)

  elif head == 'line':
    if calc_command_found: return
    parse_line(rest_of_line, sad_info)

  elif head in sad_ele_type_names:
    if calc_command_found: return
    parse_ele(head, rest_of_line, sad_info)

  elif head == 'calc' or head == 'cal':
    calc_command_found = True

  elif 'initialorbit' in rest_of_line:
    line = rest_of_line.partition('initialorbit')[2]
    line = line.partition('{')[2].partition('}')[0]
    orbit = line.replace(',', ' ').split()
    sad_info.param_list['x_orb']  = orbit[0]
    sad_info.param_list['px_orb'] = orbit[1]
    sad_info.param_list['y_orb']  = orbit[2]
    sad_info.param_list['py_orb'] = orbit[3]
    sad_info.param_list['z_orb']  = orbit[4]
    sad_info.param_list['pz_orb'] = orbit[5]

  elif not calc_command_found and len(rest_of_line) > 1 and rest_of_line[0] == '=':  # Parameter
    # Might be a variable def (EG "xxx = 7"). 
    # But if there are any special characters then ignore
    for c in '"[$,@{\'>=': 
      if c in head or c in rest_of_line[1:]: return
    if head in sad_info.var_list:
      print ('WARNING!! TRANSLATOR CANNOT HANDLE VARIABLE REDEFINITION OF: ' + head + '\n' + 
             '     WILL IGNORE FIRST DEFINITION!\n' +
             '     THIS CAN CAUSE THE BMAD LATTICE TO NOT BE EQUIVALENT TO THE SAD LATTICE!!\n' + 
             '     YOU HAVE BEEN WARNED!!')
    sad_info.var_list[head] = add_units(rest_of_line[1:])

#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
# Main program.

header_file = ''
mark_open = False
mark_closed = False
inputfile = None
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


#------------------------------------------------------------------
# Read in SAD file line-by-line.  Assemble lines into directives, which are delimited by a ; (colon).
# Call parse_directive whenever an entire directive has been obtained.

sad_info = sad_info_struct()
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
    parse_directive(directive[:ix], sad_info)
    directive = directive[ix+1:]

#------------------------------------------------------------------
# Get root lattice line

if 'use' not in sad_info.param_list:
  print ('NO USE STATEMENT FOUND!')
  sys.exit()

line0_name = sad_info.param_list['use']

if line0_name not in sad_info.lat_line_list:
  print ('USED LINE NOT FOUND. STOPPING HERE.')
  sys.exit()

sad_line = sad_info.lat_line_list[line0_name]

# For betax and betay translations

ele0_name = sad_line.list[0].name
if ele0_name in sad_info.ele_list:
  ele0 = sad_info.ele_list[ele0_name]
  for key in ele0.param:
    if key in sad_param_names:
      sad_info.param_list[key] = ele0.param[key]

#------------------------------------------------------------------
# Header

if header_file != '':
  h_file = open(header_file, 'r')
  for line in h_file:
    f_out.write (line)

#------------------------------------------------------------------
# Translate and write parameters

if mark_open: f_out.write ('parameter[geometry] = open\n')
if mark_closed: f_out.write ('parameter[geometry] = closed\n')

for name in sad_info.param_list:
  if name not in global_param_translate: continue
  if global_param_translate[name] != '':
    if '=' in global_param_translate[name]:   # EG: 'parameter[geometry] = open'
      f_out.write(global_param_translate[name] + '\n')
    else:
      f_out.write(global_param_translate[name] + ' = ' + sad_info.param_list[name] + '\n')

#------------------------------------------------------------------
# Write variable definitions

f_out.write ('\n')

for var in sad_info.var_list:
  f_out.write (var + ' = ' + sad_info.var_list[var] + '\n')

#------------------------------------------------------------------
# Translate and write element defs

inside_sol = False
bz = '0'

output_line (sad_line, sad_info, inside_sol, bz)

#-------------------------------------------------------------------

f_out.write ('\n')
f_out.write ('use, ' + line0_name + '\n')

print ('Execution time: ' + str(time.time() - start_time))

f_in.close()
f_out.close()

