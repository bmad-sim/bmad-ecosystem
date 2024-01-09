#!/usr/bin/python

import sys, getopt, re, math, copy
from collections import *
import time
import subprocess

start_time = time.time()

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
    self.sign = sign
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
    self.ix_null = 0           # Index used for generating unique null_ele names

#------------------------------------------------------------------
#------------------------------------------------------------------

def add_units (line):
  if len(line) < 4: return line
  if line[-4] not in ' 0123456789': return line
  if line[-3:] == 'deg': return line[:-3] + ' * degrees'
  if line[-3:] == 'kev': return line[:-3] + ' * 1e3'
  if line[-3:] == 'mev': return line[:-3] + ' * 1e6'
  if line[-3:] == 'gev': return line[:-3] + ' * 1e9'
  return line

#------------------------------------------------------------------
#------------------------------------------------------------------

def print_help():
  print (''' \
Syntax: 
  sad_to_bmad.py <parameter_file> <sad_lattice_file>
Defaults: 
  <parameter_file>   = "sad_to_bmad.params"
  <sad_lattice_file> = As specified in the parameter file.
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
#     '7*3  ->  '7*3'     If operand = '/'
#     '7*3  ->  '(7*3)'   If operand != '/'

def add_parens (str, operand):

  if operand == '-':
    if str.lstrip()[0] == '-':
      return '(' + str + ')'
    else:
      return str

  #
  slash_here = (operand == '/')

  for ix in range(1, len(str)):
    if slash_here and (str[ix] == '*' or str[ix] == '/'): return '(' + str + ')'
    if str[ix] != '+' and str[ix] != '-': continue
    if str[ix-1] == 'e' or str[ix-1] == 'E': continue  # '+' in '3.0e+7' is not an operator
    return '(' + str + ')'

  # No +/- op found so just return the expression.
  return str    

#------------------------------------------------------------------
#------------------------------------------------------------------

ele_type_to_bmad = {
  'drift':     'drift',
  'bend':      'sbend',
  'quad':      'quadrupole',
  'sext':      'sextupole',
  'oct':       'octupole',
  'mult':      'sad_mult',
  'sol':       'marker',
  'cavi':      'rfcavity',
  'moni':      'monitor',
  'map':       'marker',
  'mark':      'marker',
  'beambeam':  'beambeam',
  'apert':     'rcollimator',
  'coord':     'patch',
}

# Translation rule: Specific translations (of the form 'element:parameter') take
# precedence over generaic translations (of the form 'parameter').

# Something like ['k0', ' / @l@'] is translated to "<k0> / <l>" where <k0> and <l> are the values
# of the k0 and l parameters.

ele_param_translate = {
    'apert:dx': 'x_offset',
    'apert:dy': 'y_offset',
    'bend:k0': ['dg', ' / @l@'],
    'bend:k1': ['k1', ' / @l@'],
    'bend:fb1': ['hgap', '/6, fint = 0.5'],
    'bend:fb2': ['hgapx', '/6, fintx = 0.5'],
    'quad:k1': ['k1', ' / @l@'],
    'sext:k2': ['k2', ' / @l@'],
    'oct:k3': ['k3', ' / @l@'],
    'sad_mult:k1': ['b1', ' / @l@'],            # In case a SAD quad -> Bmad sad_mult
    'sad_mult:k2': ['b2', ' / (2 * @l@)'],      # In case a SAD sext -> Bmad sad_mult
    'sad_mult:k3': ['b3', ' / (6 * @l@)'],      # In case a SAD oct  -> Bmad sad_mult
    'bend:rotate': ['ref_tilt', ' * -1'],
    'bend:drotate': ['roll', ' * -1'],
    'l': 'l',
    'radius': 'aperture',
    'offset': 'offset',
    'angle': 'angle',
    'ae1': 'ae1',
    'ae2': 'ae2',
    'e1': ['e1', ' * @angle@'],
    'e2': ['e2', ' * @angle@'],
    'sol:bz': 'sad_bz',
    'sol:f1': 'sad_f1',
    'geo': 'sad_geo',
    'bound': 'sad_bound',
    'bz': 'bs_field',
    'f1': 'fq1',       # For mult and quad elements
    'f2': 'fq2',
    'eps': 'eps_step_scale',
    'freq': 'rf_frequency',
    'phi': ['phi0', ' / twopi'],
    'dphi': ['phi0_err', ' / twopi'],
    'volt': 'voltage',
    'harm': 'harmon',
    'mult:dx': 'x_offset_mult',
    'mult:dy': 'y_offset_mult',
    'sol:chi1': ['x_pitch', ' * -1'],
    'sol:chi2': ['y_pitch', ' * -1'],
    'sol:dz': ['t_offset', ' / c_light'],
    'dx': 'x_offset',
    'dy': 'y_offset',
    'dz': 'z_offset',
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

# Stuff to ignore or stuff that must be handled specially.

ignore_sad_param = ['ldev', 'fringe', 'disfrin', 'disrad', 'r1', 'r2', 'r3', 'r4', 'betax', 'betay',
                  'index', 'ex', 'ey', 'ax', 'ay', 'bx', 'by', 'detr', 'sigmaz', 'sige', 'emitz', 
                  'epx', 'epy', 'dpx', 'dpy', 'emitx', 'emity', 'dp', 'psix', 'psiy', 'psiz',
                  'sigx', 'sigy', 'sigz', 'slice', 'sturn', 'xangle', 'np', 'ddp', 
                  'pex', 'pepx', 'pey', 'pepy', 'trx', 'try', 'leng', 'ax', 'ay', 'dx1', 'dx2', 'dy1', 'dy2',
                  'zpx', 'zx', 'az', 'mark:bz']

#

sad_reversed_param = {
      'ae1': 'ae2',
      'e1':  'e2',
      'fb1': 'fb2',
}

sad_reverse_sign_flip_param = ['offset', 'bz', 'dz']

#------------------------------------------------------------------
#------------------------------------------------------------------

def output_lattice_line (sad_line, sad_info, sol_status, bz, rf_list):

  f_out.write ('\n')

  bmad_line = []
  sad_line.printed = True

  # If last element is end_marker then ignore this since Bmad will naturally put in an end marker.

  if sad_line.list[-1].name == 'end_marker': del sad_line.list[-1]

  # Loop over all SAD elements

  for ix_s_ele, sad_line_ele in enumerate(sad_line.list):

    ele_name = sad_line_ele.name

    # If the line element is itself a line then print this line info.

    if ele_name in sad_info.lat_line_list:
      if not sad_info.lat_line_list[ele_name].printed: 
        output_lattice_line(sad_info.lat_line_list[ele_name], sad_info, sol_status, bz, rf_list)
      bmad_line.append(sad_line_ele)
      continue

    if not ele_name in sad_info.ele_list:
      print ('No definition found for element name: ' + ele_name)
      continue

    sad_ele_def = sad_info.ele_list[ele_name]

    # Reversed and not longitudinally symmetric?
    # If so create a new reversed element

    if sad_line_ele.sign == '-':
      symmetric = True
      for pname in sad_reversed_param:
        if sad_ele_def.param.get(pname, '0') != sad_ele_def.param.get(sad_reversed_param[pname], '0'): symmetric = False
      for pname in sad_reverse_sign_flip_param:
        if pname in sad_ele_def.param: symmetric = False

      if not symmetric:
        ele_name = ele_name + '_inverse'
        sad_line_ele.name = ele_name

        if ele_name in sad_info.ele_list:
          sad_ele_def = sad_info.ele_list[ele_name]

        else:
          sad_ele_def = copy.deepcopy(sad_ele_def)
          sad_ele_def.name = ele_name
          sad_ele_def.printed = False

          for pname in sad_reversed_param:
            rname = sad_reversed_param[pname]
            if pname in sad_ele_def.param and rname in sad_ele_def.param:
              sad_ele_def.param[pname], sad_ele_def.param[rname] = sad_ele_def.param[rname], sad_ele_def.param[pname]
            elif pname in sad_ele_def.param:
              sad_ele_def.param[rname] = sad_ele_def.param[pname]
              del sad_ele_def.param[pname]
            elif rname in sad_ele_def.param:
              sad_ele_def.param[pname] = sad_ele_def.param[rname]
              del sad_ele_def.param[rname]

          for pname in sad_reverse_sign_flip_param:
            if pname in sad_ele_def.param: 
              if sad_ele_def.param[pname][0] == '-':
                sad_ele_def.param[pname] = sad_ele_def.param[pname][1:]
              else:
                sad_ele_def.param[pname] = '-' + sad_ele_def.param[pname]

          sad_info.ele_list[ele_name] = sad_ele_def

    # sol element

    if sad_ele_def.type == 'sol':
      if sad_ele_def.param.get('bound') == '1': 
        if sol_status == 0:            # If was outside solenoid..
          if sad_line_ele.sign == '':
            sol_status = 1             # Now inside
          else:
            sol_status = -1            # Now inside a reversed solenoid
        else:                          # If was inside solenoid...
          sol_status = 0               # Now outside solenoid

    if sad_ele_def.type == 'sol':
      if sol_status == 1:
        bz = sad_ele_def.param.get('bz', '0')
      elif sol_status == -1:
        for ise in range(ix_s_ele+1, len(sad_line.list)):
          s2_ele = sad_line.list[ise]
          s2_type = sad_info.ele_list[s2_ele.name]
          if s2_type.type == 'sol':
            bz = s2_type.param.get('bz', '0')
            break

      if sol_status == 1 or sol_status == -1:
        try:
          b_z = float(bz)
          if b_z == 0: bz = '0'    # EG convert '0.0' to '0'
        except ValueError:
          pass

      if sol_status == -1:
        if bz[0] == '-':
          bz = bz[1:]
        else:
          bz = '-' + bz

    # A MARK element with an offset gets translated to a marker superimpsed with respect to a null_ele

    if sad_ele_def.type == 'mark' and 'offset' in sad_ele_def.param:
      if ignore_marker_offsets:
        del sad_ele_def.param['offset']
      else:
        sad_info.ix_null += 1
        null_ele_name = 'null_' + sad_ele_def.name + '#' + str(sad_info.ix_null)   # Guaranteed unique
        bmad_line.append (line_item_struct(null_ele_name))          # Put null_ele in the line
        WrapWrite(null_ele_name + ': null_ele')                     # Define the null_ele
  
        # Now define the marker element
        sad_offset = float(sad_ele_def.param['offset'])
        int_off = int(math.floor(sad_offset))
        frac_off = sad_offset - int_off
        direc = 1
        if int_off < 0: 
          direc = -1
          frac_off = frac_off - 1

        offset = 0
        for ix in range(0, int_off, direc):
          this_name = sad_line.list[ix_s_ele+ix].name
          if this_name in sad_info.ele_list:
            sad_ele_def2 = sad_info.ele_list[this_name]
            if 'l' in sad_ele_def2.param: offset += direc * eval(sad_ele_def2.param['l'])
          else:   # Must be a line
            for sub_ele in sad_info.lat_line_list[this_name].list:
              sad_ele_def2 = sad_info.ele_list[sub_ele.name]
              if 'l' in sad_ele_def2.param: offset += direc * eval(sad_ele_def2.param['l'])
  
        sad_ele_def2 = sad_info.ele_list[sad_line.list[ix_s_ele+int_off].name]
        if 'l' in sad_ele_def2.param: offset += frac_off * eval(sad_ele_def2.param['l'])
  
        if sad_ele_def.instances == 0:
          suffix = ''
        else:
          suffix = '.' + str(sad_ele_def.instances)
        bmad_ele_def = sad_ele_def.name + suffix + ': marker, superimpose, ref = ' + null_ele_name + ', offset = ' + str(offset)
        WrapWrite(bmad_ele_def)
        sad_ele_def.printed = True
        sad_ele_def.instances += 1
        continue

    # Regular element not getting superimposed

    bmad_line.append(sad_line_ele)
    if sad_ele_def.type == 'cavi': rf_list.append(sad_ele_def.name)

    if sad_ele_def.printed == True: continue

    b_ele = ele_struct()
    sad_ele_to_bmad (sad_ele_def, b_ele, sol_status, bz, sad_line_ele.sign == '-')

    bmad_ele_def = b_ele.name + ': ' + b_ele.type
    for param in iter(b_ele.param):
      try:
        val = float(b_ele.param[param])
        if val == 0: continue
      except:
        pass
      bmad_ele_def += ', ' + param + ' = ' + b_ele.param[param]

    WrapWrite(bmad_ele_def)
    sad_ele_def.printed = True

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

def sad_ele_to_bmad (sad_ele, bmad_ele, sol_status, bz, reversed):

  bmad_ele.name = sad_ele.name

  #

  if 'l' in sad_ele.param:
    try:
      zero_length = (float(sad_ele.param['l']) == 0)
      if zero_length: sad_ele.param['l'] = '0'
    except ValueError:
      zero_length = False
  else:
    zero_length = True

  if 'bz' in sad_ele.param:
    try:
      zero_field = (float(sad_ele.param['bz']) == 0)
      if zero_field: sad_ele.param['bz'] = '0'
    except ValueError:
      zero_field = False
  else:
    zero_field = True

  #

  if not sad_ele.type in ele_type_to_bmad:
    print ('TYPE OF ELEMENT NOT RECOGNIZED: ' + sad_ele.type + '\n' + '    FOR ELEMENT: ' + sad_ele.name)
    return

  bmad_ele.type = ele_type_to_bmad[sad_ele.type]

  # SAD sol with misalignments becomes a Bmad patch

  if sad_ele.type == 'sol' and sad_ele.param.get('bound') == '1':
    if len(set(['dx', 'dy', 'dz', 'chi1', 'chi2', 'chi3', 'rotate']).intersection(sad_ele.param)) > 0:
      bmad_ele.type = 'patch'

  # Handle case when inside solenoid

  if sol_status != 0 and bmad_ele.type != 'marker' and bmad_ele.type != 'monitor' and \
                    bmad_ele.type != 'patch' and bz != '0' and not zero_length:
    bmad_ele.param['bs_field'] = bz
    if bmad_ele.type == 'drift':
      bmad_ele.type = 'solenoid'
    else:
      bmad_ele.type = 'sad_mult'

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

  if sad_ele.type == 'bend':

    if 'f1' in sad_ele.param and 'fb1' in sad_ele.param:
      bmad_ele.param['hgap'] = '(' + sad_ele.param['f1'] + ' + ' + sad_ele.param['fb1'] + ')/6, fint = 0.5'
    elif 'f1' in sad_ele.param:
      bmad_ele.param['hgap'] = '(' + sad_ele.param['f1'] + ')/6, fint = 0.5'
    elif 'fb1' in sad_ele.param:
      bmad_ele.param['hgap'] = '(' + sad_ele.param['fb1'] + ')/6, fint = 0.5'

    # If fb1 == fb2 then don't need fb2 
    if sad_ele.param.get('fb1', '0') != sad_ele.param.get('fb2', '0'):
      if 'f1' in sad_ele.param and 'fb2' in sad_ele.param:
        bmad_ele.param['hgapx'] = '(' + sad_ele.param['f1'] + ' + ' + sad_ele.param['fb2'] + ')/6, fintx = 0.5'
      elif 'f1' in sad_ele.param:
        bmad_ele.param['hgapx'] = '(' + sad_ele.param['f1'] + ')/6, fintx = 0.5'
      elif 'fb2' in sad_ele.param:
        bmad_ele.param['hgapx'] = '(' + sad_ele.param['fb2'] + ')/6, fintx = 0.5'

  # Loop over all parameters

  for sad_param_name in sad_ele.param:

    if sad_ele.type == 'bend' and sad_param_name == 'f1': continue
    if sad_ele.type == 'bend' and sad_param_name == 'fb1': continue
    if sad_ele.type == 'bend' and sad_param_name == 'fb2': continue

    value = sad_ele.param[sad_param_name]

    try:
      zero_value = (float(value) == 0)
      if zero_value: value = '0'
    except ValueError:
      zero_value = False

#   # geo and bound parameters
#
#   if sad_param_name == 'geo' or sad_param_name == 'bound':
#     bmad_ele.param['sad_' + sad_param_name] = value
#
#     continue

    #

    full_param_name = sad_ele.type + ':' + sad_param_name 
    if sad_param_name in ignore_sad_param: continue
    if full_param_name in ignore_sad_param: continue

    if bmad_ele.type == 'sad_mult' and sad_ele.type != 'mult':    # EG: In a solenoid field, SAD quad -> Bmad sad_mult
      full_param_name = bmad_ele.type + ':' + sad_param_name 

    # Use more specific translation first

    if full_param_name in ele_param_translate:
      result = ele_param_translate[full_param_name]
    elif sad_param_name in ele_param_translate:
      result = ele_param_translate[sad_param_name]
    else:
      print ('SAD PARAMETER NOT RECOGNIZED: ' + sad_param_name + '\n' + '    DEFINED IN ELEMENT: ' + sad_ele.name)
      continue

    # 

    value_suffix = ''
    bmad_name = sad_param_name

    if sad_param_name == 'k0' and sad_ele.type == 'bend' and zero_length:
      bmad_ele.type = 'multipole'
      bmad_name = 'k0l'

    elif sad_param_name == 'k1' and sad_ele.type == 'bend' and zero_length:
      bmad_ele.type = 'multipole'
      bmad_name = 'k1l'

    elif sad_param_name == 'k1' and sad_ele.type == 'quad' and zero_length:
      bmad_ele.type = 'multipole'
      bmad_name = 'k1l'

    elif sad_param_name == 'k2' and sad_ele.type == 'sext' and zero_length:
      bmad_ele.type = 'multipole'
      bmad_name = 'k2l'

    elif sad_param_name == 'k3' and sad_ele.type == 'oct' and zero_length:
      bmad_ele.type = 'multipole'
      bmad_name = 'k3l'

    elif sad_param_name == 'rotate' and sad_ele.type == 'quad' and zero_length:
      bmad_ele.type = 'multipole'
      bmad_name = 't1'
      value_suffix = ''

    elif sad_param_name == 'rotate' and sad_ele.type == 'sext' and zero_length:
      bmad_ele.type = 'multipole'
      bmad_name = 't2'
      value_suffix = ''

    elif sad_param_name == 'rotate' and sad_ele.type == 'oct' and zero_length:
      bmad_ele.type = 'multipole'
      bmad_name = 't3'
      value_suffix = ''

    elif type(result) is list:   # EG: result = ['k1', ' / @l@']
      bmad_name = result[0]
      value_suffix = result[1]
      value = add_parens(value, '')
      if '@' in value_suffix:
        val_parts = value_suffix.split('@')
        if val_parts[1] in sad_ele.param:
          operand = '/' if '/' in val_parts[0] else ''
          value_suffix = val_parts[0] + add_parens(sad_ele.param[val_parts[1]], operand) + val_parts[2]
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

    if reversed and sad_ele.type != 'sol':
      if sad_param_name == 'offset': value = '1 - ' + value
      if sad_param_name in ['chi1', 'chi2', 'dz', 'bz']: 
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

  if sad_ele.type == 'mult' and ('freq' in sad_ele.param or 'phi' in sad_ele.param): bmad_ele.type = 'rfcavity'

  # Use ptc_standard type cavities

  if bmad_ele.type == 'rfcavity' or bmad_ele.type == 'lcavity':
    bmad_ele.param['cavity_type'] = 'ptc_standard'

  # If the SAD cavi has a nonzero phi then the reference particle's energy is changing and
  # in this case the corresponding Bmad element must be an lcavity.
  # Also remember that an lcavity has a differnt phase convention.

  if bmad_ele.type == 'rfcavity':
    if 'phi0' in bmad_ele.param and 'harmon' not in bmad_ele.param: # If has harmon then must be in a ring
      bmad_ele.type = 'lcavity'
      bmad_ele.param['phi0'] = '0.25 + ' + add_parens(bmad_ele.param['phi0'], '')
      if 'phi0_err' in bmad_ele.param: bmad_ele.param['phi0_err'] = '+' + add_parens(bmad_ele.param['phi0_err'], '')

    elif 'phi0_err' in bmad_ele.param:
      if 'phi0' in bmad_ele.param:
        bmad_ele.param['phi0'] = bmad_ele.param['phi0'] + ' + ' + bmad_ele.param['phi0_err']
      else:
        bmad_ele.param['phi0'] = bmad_ele.param['phi0_err']
      del bmad_ele.param['phi0_err']

  # Correct patch signs. 
  # And SAD applies pitches and offsets in reverse order to Bmad which is why the trig is needed.
  # Note: When a SOL at the solenoid boundary is translated into a patch, the dz will be translated to a t_offset and
  #  the z_offset, before rotation by a finite x_pitch or y_pitch, will be 0.

  if bmad_ele.type == 'patch':

    # If exiting solenoid or not in a solenoid.
    if sol_status == 0:
      if 'z_offset' in bmad_ele.param: bmad_ele.param['z_offset'] = str(eval(bmad_ele.param['z_offset'] + ' * -1'))

    # If entering solenoid 
    else:
      if 'x_offset' in bmad_ele.param: bmad_ele.param['x_offset'] = str(eval(bmad_ele.param['x_offset'] + ' * -1'))
      if 'y_offset' in bmad_ele.param: bmad_ele.param['y_offset'] = str(eval(bmad_ele.param['y_offset'] + ' * -1'))
      if 'z_offset' in bmad_ele.param: bmad_ele.param['z_offset'] = str(eval(bmad_ele.param['z_offset'] + ' * -1'))


    zo = eval(bmad_ele.param.get('z_offset', '0'))
    xo = eval(bmad_ele.param.get('x_offset', '0'))
    xp = eval(bmad_ele.param.get('x_pitch', '0'))
    yo = eval(bmad_ele.param.get('y_offset', '0'))
    yp = eval(bmad_ele.param.get('y_pitch', '0'))

    if xp != 0 and xo != 0:
      bmad_ele.param['z_offset'] = str(zo * math.cos(xp) - xo * math.sin(xp))
      bmad_ele.param['x_offset'] = str(zo * math.sin(xp) + xo * math.cos(xp))
      zo = eval(bmad_ele.param.get('z_offset', '0'))
    if yp != 0 and yo != 0:
      bmad_ele.param['z_offset'] = str(zo * math.cos(xp) - yo * math.sin(yp))
      bmad_ele.param['y_offset'] = str(zo * math.sin(xp) + yo * math.cos(yp))

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
      bmad_ele.param['fringe_type'] = 'soft_edge_only'

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
      if line_item.name == 'end': line_item.name = 'end_marker'  # 'end' is a reserved name
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

      if ix > len(line) - 2:
        print ('MALFORMED ELEMENT DEFINITION: ' + rest_of_line)
        sys.exit()

      if line[ix] == '=' or line[ix] == '(':    # Looking for "ename = (..." or "emane (..."
        ele = ele_struct()
        ele.type = head
        ele.name = line[:ix].strip()
        ## print ('Name: "' + ele.name + '"')
        if line[ix] == '(':
          line = line[ix:].lstrip()
        else:
          line = line[ix+1:].lstrip()
        break

    if line[0] != '(':
      print ('MALFORMED ELEMENT DEFINITION. EXPECTING "(": ' + rest_of_line)
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

# For translating parameters of marker at beginning of lattice.

sad_ele0_param_names = ['momentum', 'use', 'bx', 'by', 'ax', 'ay', 
        'ex', 'epx', 'ey', 'epy', 'emitx', 'emity', 'nocod']  

#
global_param_translate = {
  'momentum': 'beginning[p0c]',
  'bx':       'beginning[beta_a]',
  'by':       'beginning[beta_b]',
  'ax':       'beginning[alpha_a]',
  'ay':       'beginning[alpha_b]',
  'ex':       'beginning[eta_x]',
  'epx':      'beginning[etap_x]',
  'ey':       'beginning[eta_y]',
  'epy':      'beginning[etap_y]',
  'emitx':    'particle_start[emittance_a]',
  'emity':    'particle_start[emittance_b]',
  'nocod':    'parameter[geometry] = open',
  'cod':      'parameter[geometry] = closed',
  'x_orb':    'particle_start[x]',
  'px_orb':   'particle_start[px]',
  'y_orb':    'particle_start[y]',
  'py_orb':   'particle_start[py]',
  'z_orb':    'particle_start[z]',
  'pz_orb':   'particle_start[pz]',
  'dxi':      'particle_start[x]',
  'dpxi':     'particle_start[px]',
  'dyi':      'particle_start[y]',
  'dpyi':     'particle_start[py]',
  'dzi':      'particle_start[z]',
  'ddpi':     'particle_start[pz]',
  'npara':    ''
}

#------------------------------------------------------------------
#------------------------------------------------------------------

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

if sys.version_info.major != 3:
  print ('This script requires Python 3!\n')
  sys.exit(1)

# Read the parameter file specifying the SAD lattice file, etc.

param_file = "sad_to_bmad.params"

if len(sys.argv) == 1:
  exec (open(param_file).read())

if len(sys.argv) == 2:
  param_file = sys.argv[1]
  exec (open(param_file).read())

elif len(sys.argv) == 3:
  param_file = sys.argv[1]
  exec (open(param_file).read())
  sad_lattice_file = sys.argv[2]

elif len(sys.argv) > 3:
  print_help()

# Construct the bmad lattice file name

if bmad_lattice_file == '':
  if sad_lattice_file.find('sad') != -1:
    bmad_lattice_file = sad_lattice_file.replace('sad', 'bmad')
  elif sad_lattice_file.find('Sad') != -1:
    bmad_lattice_file = sad_lattice_file.replace('Sad', 'bmad')
  elif sad_lattice_file.find('SAD') != -1:
    bmad_lattice_file = sad_lattice_file.replace('SAD', 'bmad')
  else:
    bmad_lattice_file = sad_lattice_file + '.bmad'

print ('Input lattice file is:  ' + sad_lattice_file)
print ('Output lattice file is: ' + bmad_lattice_file)

# Open files for reading and writing

f_in = open(sad_lattice_file, 'r')
f_out = open(bmad_lattice_file, 'w')
f_out.write ('! Translated from SAD file: ' + sad_lattice_file + "\n\n")

sad_ele_type_names = ("drift", "bend", "quad", "sext", "oct", "mult", "sol", "cavi", "map", "moni", "line", "beambeam", "apert", "mark", "coord")

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
for i in range(100):
  if ele0_name not in sad_info.lat_line_list: break
  ele0_name = sad_info.lat_line_list[ele0_name].list[0].name

ele0 = sad_info.ele_list[ele0_name]
for key in ele0.param:
  if key in sad_ele0_param_names:
    sad_info.param_list[key] = ele0.param[key]

#------------------------------------------------------------------
# Header

f_out.write (header_lines + '\n')

#------------------------------------------------------------------
# Translate and write parameters

if lattice_geometry != '': f_out.write ('parameter[geometry] = ' + lattice_geometry + '\n')

for name in sad_info.param_list:
  if name not in global_param_translate: continue
  if global_param_translate[name] != '':
    if global_param_translate[name][:19] == 'parameter[geometry]':
      if lattice_geometry != '': f_out.write(global_param_translate[name] + '\n')
    elif '=' in global_param_translate[name]: 
      f_out.write(global_param_translate[name] + '\n')
    else:
      f_out.write(global_param_translate[name] + ' = ' + sad_info.param_list[name] + '\n')

# The SuperKEK-B sler lattice may need PTC_exact_model = True

f_out.write ('parameter[ptc_exact_model] = true\n')

# If there is a SOL element with an F1 attribute. See the DOC file for more info.

f_out.write('''
! Save SAD SOL F1 and other info in a custom attribute in case lattice is back translated to to SAD
parameter[custom_attribute1] = "marker::sad_f1"
parameter[custom_attribute1] = "patch::sad_f1"
parameter[custom_attribute2] = "marker::sad_geo"  
parameter[custom_attribute2] = "patch::sad_geo"   
parameter[custom_attribute3] = "marker::sad_bound"
parameter[custom_attribute3] = "patch::sad_bound" 
parameter[custom_attribute4] = "marker::sad_bz"
parameter[custom_attribute4] = "patch::sad_bz"
parameter[custom_attribute5] = "patch::sad_fshift"
''')

# If the first element is a marker with Twiss parameters...


#------------------------------------------------------------------
# Write variable definitions

if patch_for_fshift != 'MAYBE' and patch_for_fshift != 'TRUE' and patch_for_fshift != 'FALSE':
  print ('Possible settings for patch_for_fshift are: "MAYBE", "TRUE", or "FALSE".')
  print ('I suspect you are using an old version of of the sad_to_bmad.params file.')
  sys.exit()

if patch_for_fshift == 'MAYBE':
  if 'fshift' in sad_info.var_list:
    if float(sad_info.var_list['fshift']) == 0: 
      patch_for_fshift = 'FALSE'
    else:
      patch_for_fshift = 'TRUE'
  else:
    patch_for_fshift = 'FALSE'

if patch_for_fshift == 'TRUE' and 'fshift' not in sad_info.var_list: sad_info.var_list['fshift'] = '0'

f_out.write ('\n')

for var in sad_info.var_list:
  f_out.write (var + ' = ' + sad_info.var_list[var] + '\n')

#------------------------------------------------------------------
# Translate and write element defs

sol_status = 0
bz = '0'

rf_list = []
output_lattice_line (sad_line, sad_info, sol_status, bz, rf_list)

#-------------------------------------------------------------------

f_out.write ('\n')
f_out.write ('use, ' + line0_name + '\n')

#------------------------------------------------------------------
# Footer

f_out.write ('\n' + footer_lines)

print ('Execution time: ' + str(time.time() - start_time))

#-------------------------------------------------------------------
# Insert patches for finite fshift 

if patch_for_fshift == 'TRUE':
  f_out.write ('\n' + 'expand_lattice\n')
  f_out.write ('t_scale = 1\n')

  fshift = sad_info.var_list.get('fshift', '1e-30') # Default is just some small non-zero number 

  rf_dict = {}
  for rf_name in rf_list:
    if rf_name in rf_dict:
      rf_dict[rf_name] = rf_dict[rf_name] + 1
    else:
      rf_dict[rf_name] = 1

    ns = str(rf_dict[rf_name])
    full_rf_name = rf_name + '##' + ns
    patch_name = rf_name + '_patch' + ns
    f_out.write ('t_' + patch_name + ' = 0  ! Will be replaced by sad_to_bmad_postprocess\n')
    f_out.write (patch_name + ': patch, superimpose, ref_origin = beginning, ref = ' + full_rf_name +
                 ',\n    sad_fshift = ' + fshift + ', t_offset = t_scale * t_' + patch_name + '\n')

  patch_name = 'last_rf_time_patch'
  f_out.write ('t_' + patch_name + ' = 0  ! Will be replaced by sad_to_bmad_postprocess\n')
  f_out.write (patch_name + ': patch, superimpose, ref_origin = end, ref = ' + full_rf_name +
                 ',\n    sad_fshift = ' + fshift + ', t_offset = t_scale * t_' + patch_name + '\n')

#-------------------------------------------------------------------

f_in.close()
f_out.close()

if patch_for_fshift == 'TRUE':
  command = sad_to_bmad_postprocess_exe + ' ' + bmad_lattice_file + ' ' + calc_fshift_for
  print (f'\nRunning sad_to_bmad_postprocess to complete the translation. Command is:\n   {command}')
  subprocess.call (command, shell = True)
