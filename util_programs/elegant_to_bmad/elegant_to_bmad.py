#!/usr/bin/env python3

#+
# Script to convert from Elegant lattice format to Bmad lattice format.
# See the README file for more details.
#-

import sys, re, argparse, time
import math as m

from collections import OrderedDict

if sys.version_info[0] < 3 or sys.version_info[1] < 6:
  raise Exception("Must be using Python 3.6+")

#------------------------------------------------------------------

class ele_struct:
  def __init__(self, name):
    self.name = name
    self.elegant_type    = ''    # Is quadrupole, etc.
    self.bmad_type    = ''       # Is quadrupole, etc.
    self.at = '0'                # Used if element is in a sequence
    self.param = OrderedDict()

class common_struct:
  def __init__(self):
    self.debug = False
    self.one_file = True
    self.add_constants = False       # Add elegant defined constants to lattice file?
    self.ele_dict = {}               # Dict of elements
    self.f_in = []         # Elegant input files
    self.f_out = []        # Bmad output files
    self.beam_line_name = ''
    self.command = ''

#------------------------------------------------------------------
#------------------------------------------------------------------

ele_type_translate = {
  'bumper':      'ac_kicker',
  'cepl':        'ac_kicker',
  'mbumper':     'ac_kicker',
  'mrfdf':       'ac_kicker',
  'rmdf':        'ac_kicker',
  'twpl':        'ac_kicker',
  'beambeam':    'beambeam',
  'rfdf':        'crab_cavity',
  'rftm110':     'crab_cavity',
  'bggexp':      'drift',
  'bmapxy':      'drift',
  'bmxyz':       'drift',
  'boffaxe':     'drift',
  'corgpipe':    'drift',
  'csrdrift':    'drift',
  'drif':        'drift',
  'edrift':      'drift',
  'fmult':       'drift',
  'gkickmap':    'drift',
  'hkpoly':      'drift',
  'lscdrift':    'drift',
  'matter':      'drift',
  'peppot':      'drift',
  'polyseries':  'drift',
  'qufringe':    'drift',
  'stray':       'drift',
  'apcontour':   'ecollimator',
  'ecol':        'ecollimator',
  'maxamp':      'ecollimator',
  'speedbump':   'ecollimator',
  'taperapc':    'ecollimator',
  'taperape':    'ecollimator',
  'floor':       'floor_shift',
  'ehkick':      'hkicker',
  'hkick':       'hkicker',
  'ekicker':     'kicker',
  'kicker':      'kicker',
  'lthinlens':   'lens',
  'slice':       'maker',
  'alph':        'marker',
  'branch':      'marker',
  'charge':      'marker',
  'clean':       'marker',
  'cpickup':     'marker',
  'dscatter':    'marker',
  'emittance':   'marker',
  'frfmode':     'marker',
  'ftrfmode':    'marker',
  'histogram':   'marker',
  'ibscatter':   'marker',
  'ioneffects':  'marker',
  'lrwake':      'marker',
  'magnify':     'marker',
  'mark':        'marker',
  'mhistogram':  'marker',
  'pfilter':     'marker',
  'rampp':       'marker',
  'recirc':      'marker',
  'reflect':     'marker',
  'remcor':      'marker',
  'rfmode':      'marker',
  'rimult':      'marker',
  'scatter':     'marker',
  'scmult':      'marker',
  'script':      'marker',
  'shrfdf':      'marker',
  'sreffects':   'marker',
  'tfbdriver':   'marker',
  'tfbpickup':   'marker',
  'trcount':     'marker',
  'trfmode':     'marker',
  'trwake':      'marker',
  'tscatter':    'marker',
  'twiss':       'marker',
  'wake':        'marker',
  'watch':       'marker',
  'zlongit':     'marker',
  'ztransverse': 'marker',
  'ckicker':     'marker',
  'ilmatrix':    'match',
  'lmirror':     'mirror',
  'hmon':        'monitor',
  'moni':        'monitor',
  'vmon':        'monitor',
  'mult':        'multipole',
  'octu':        'octuple',
  'koct':        'octupole',
  'center':      'patch',
  'energy':      'patch',
  'malign':      'gkicker',
  'rotate':      'marker', #  [exclude_floor != 0, exclude_optics != 0]
                           # taylor      [exclude_floor != 0, exclude_optics == 0]
                           # floor_shift [exclude_floor == 0, exclude_optics != 0]
                           # patch       [exclude_floor == 0, exclude_optics == 0]
  'kquad':       'quadrupole',
  'kquse':       'quadrupole',
  'quad':        'quadrupole',
  'rben':        'rbend',
  'tubend':      'rbend',
  'rcol':        'rcollimator',
  'scraper':     'rcollimator',
  'taperapr':    'rcollimator',
  'rftmez0':     'rf_cavity',
  'modrf':       'rfcavity',
  'ramprf':      'rfcavity',
  'rfca':        'rfcavity',
  'rfcw':        'rfcavity',
  'tmcf':        'rfcavity',
  'twla':        'rfcavity',
  'twmta':       'rfcavity',
  'sample':      'sample',
  'brat':        'sbend',
  'ccbend':      'sbend',
  'csbend':      'sbend',
  'csrcsbend':   'sbend',
  'nibend':      'sbend',
  'nisept':      'sbend',
  'sben':        'sbend',
  'ftable':      'sbend',
  'sext':        'sextuple',
  'ksext':       'sextupole',
  'mapsolenoid': 'solenoid',
  'sole':        'solenoid',
  'ematrix':     'taylor',
  'matr':        'taylor',
  'kpoly':       'taylor',
  'lsrmdltr':    'undulator',
  'ukickmap':    'undulator',
  'evkick':      'vhicker',
  'vkick':       'vkicker',
  'cwiggler':    'wiggler',
  'gfwiggler':   'wiggler',
  'wiggler':     'wiggler',
}

problematical_translation_list = ['corgpipe', 'hkpoly', 'lscdrift', 'matter', 'peppot', 'qufringe', 'maxamp', 
          'slice', 'alph', 'branch', 'charge', 'clean', 'cpickup', 'emittance', 'frfmode', 'ftrfmode',
          'ibscatter', 'ioneffects', 'lrwake', 'mhistogram', 'pfilter', 'recirc', 'reflect',
          'remcor', 'rfmode', 'rimult', 'scatter', 'scmult', 'script', 'shrfdf', 'sreffects',
          'tfbdriver', 'tfbpickup', 'trcount', 'trfmode', 'trwake', 'tscatter', 'twiss',
          'wake', 'zlongit', 'ztransverse', 'center', 'rftm110', 'bumper', 'mbumper',
          'bggexp', 'bmapxy', 'bmxyz', 'boffaxe', 'fmult', 'gkickmap', 'polyseries', 'apcontour',
          'dscatter', 'rampp', 'brat', 'ftable', 'mapsolenoid', 'matr', 'gfwiggler', 'cepl',
          'mrfdf', 'rmdf', 'twpl', 'stray', 'speedbump', 'taperapc', 'taperape', 'lthinlens',
          'histogram', 'ckicker', 'tubend', 'scraper', 'taperapr', 'rftmez0', 'modrf', 'ramprf',
          'tmcf', 'twla', 'twmta', 'ccbend', 'csrcsbend', 'nibend', 'nisept', 'kpoly',
          'lsrmdltr', 'ukickmap', 'cwiggler']

bmad_param_name = {
  'dx':           'x_offset',
  'dy':           'y_offset',
  'dz':           'z_offset',
  'l':            'l',
  'angle':        'angle',
  'b':            'b_field',
  'e1':           'e1',
  'e2':           'e2',
  'ks':           'ks',
  'k1':           'k1',
  'k2':           'k2',
  'k3':           'k3',
  'tilt':         'tilt',  # bend ref_tilt handled in bmad_param routine
  'pitch':        'y_pitch',
  'yaw':          'x_pitch',
  'h1':           'h1',
  'h2':           'h2',
  'hgap':         'hgap',
  'hgapx':        'hgapx',
  'fint':         'fint',
  'fintx':        'fintx',
  'fint1':        'fint',
  'fint2':        'fintx',
  'xkick':        'hkick',
  'hkick':        'hkick',
  'ykick':        'vkick',
  'vkick':        'vkick',
  'kick':         'kick',
  'xcenter':      'x_offset',
  'ycenter':      'y_offset',
  'xsize':        'sig_x',
  'ysize':        'sig_y',
  'b_max':        'b_max',
  'periods':      'n_period',
  'vertical':     'tilt = pi/2',
  'helical':      'field_calc = helical_model',
  'x_max':        'x_limit',
  'y_max':        'y_limit',
  'frequency':    'rf_frequency',
  'freq':         'rf_frequency',
  'phase':        'phi0',
  'volt':         'voltage',
}

#------------------------------------------------------------------
#------------------------------------------------------------------
# Routine to parse a namelist

def namelist_dict(dlist):
  nl_name = dlist.pop(0)   # Remove "&namelist-name" and "&end"
  dlist.pop(-1)
  ndir = {}

  while True:
    if len(dlist) == 0: return ndir
    if len(dlist) < 3 or dlist[1] != '=':
      print (f'ERROR READING NAMELIST: {nl_name}')
      return ndir

    name = dlist.pop(0).lower()
    dlist.pop(0)   # Pop equal sign

    try: 
      ixe = dlist.index('=')
      ndir[name] = ''.join(dlist[:ixe-1])
      if ndir[name][-1] == ',': ndir[name] = ndir[name][:-1]
      dlist = dlist[ixe-1:]
    except:
      ndir[name] = ''.join(dlist)
      return ndir

#------------------------------------------------------------------
#------------------------------------------------------------------
# Routines for postfix to infix converter

def infixer(tokens, ix0, toplevel=True):
  ix0 -= 1
  #print (f'In: {ix0} {tokens}')
  if ix0 == -1:
    return tokens, -1

  elif tokens[ix0].upper() in ['ABS', 'TAN', 'DTAN', 'SIN', 'DSIN', 
                  'COS', 'DCOS', 'REC', 'RTOD', 'DTOR', 'SINH', 'COSH', 'TANH', 
                  'ASIN', 'ACOS', 'ATAN', 'ACOSH', 'ASINH', 'ATANH', 'LOG', 'HYPOT',
                  'MAX2', 'MIN2', 'SQR', 'SQRT']:
    fn = tokens.pop(ix0)
    tokens, ixn = infixer(tokens, ix0, False)

    if fn.upper() in ['ABS', 'TAN', 'SIN', 'COS', 'SINH', 'COSH', 'TANH', 'ASIN', 'ACOS', 'ATAN', 
              'ACOSH', 'ASINH', 'ATANH', 'LOG', 'SQRT']:
      tokens[ixn] = f'{fn}({tokens[ixn]})'

    elif fn.upper() in ['MAX2', 'MIN2']:
      arg2 = tokens.pop(ixn)
      tokens, ixn = infixer(tokens, ixn, False)
      tokens[ixn] = f'{fn}({tokens[ixn]}, arg2)'

    elif fn.upper() in ['DTAN', 'DSIN', 'DCOS']:
      tokens[ixn] = f'{fn[1:]}({tokens[ixn]}*degrees)'

    elif fn.upper() == 'REC':
      tokens[ixn] = f'1 / {tokens[ixn]}'

    elif fn.upper() == 'RTOD':
      tokens[ixn] = f'{tokens[ixn]}*raddeg'

    elif fn.upper() == 'DTOR':
      tokens[ixn] = f'{tokens[ixn]}*degrees'

    elif fn.upper() == 'HYPOT':
      arg2 = tokens.pop(ixn)
      tokens, ixn = infixer(tokens, ixn, False)
      tokens[ixn] = f'sqrt({tokens[ixn]}^2 + {arg2}^2)'

    elif fn.upper() == 'SQR':
      tokens[ixn] = f'({tokens[ixn]})^2'

    return tokens, ixn

  elif tokens[ix0] in ['+', '-', '*', '/', '^']:
    op = tokens.pop(ix0)

    if ix0 == 0 or tokens[ix0-1] in ['+', '-', '*', '/', '^']:   # Unary minus
      tokens[ix0] = f'-{tokens[ix0]}'
      return tokens, ix0

    tokens, ixn = infixer(tokens, ix0, False)
    arg2 = tokens.pop(ixn)
    tokens, ixn = infixer(tokens, ixn, False)
    arg1 = tokens[ixn]

    if op == '+':
      op = ' + '

    elif op == '-':
      if '+' in arg2 or '-' in arg2: arg2 = f'({arg2})'
      op = ' - '

    elif op == '*':
      if '+' in arg1 or '-' in arg1: arg1 = f'({arg1})'
      if '+' in arg2 or '-' in arg2: arg2 = f'({arg2})'

    elif op == '/':
      if '+' in arg1 or '-' in arg1: arg1 = f'({arg1})'
      if '+' in arg2 or '-' in arg2 or '*' in arg2 or '/' in arg2: arg2 = f'({arg2})'
      
    elif op == '^':
      if '+' in arg1 or '-' in arg1 or '*' in arg1 or '/' in arg1: arg1 = f'({arg1})'
      if '+' in arg2 or '-' in arg2 or '*' in arg2 or '/' in arg2: arg2 = f'({arg2})'

    tokens[ixn] = f'{arg1}{op}{arg2}'
    return tokens, ixn

  else:
    if not toplevel: return tokens, ix0
    
    while ix0 > -1:
      tokens, ix0 = infixer(tokens, ix0) # Notice call is toplevel

    return tokens, ix0
      
# Rearrange expression from infix to postfix format

def postfix_to_infix(str, return_list = False):

  str = str.strip('\' "')

  # Split expression into tokens and recombine numbers like "3.4e-7" into a single token.
  # Also the minus sign in something like "a -1 *" is a unary minus.
  # To avoid splitting on a unary minus, substitute "@" for all unary minusus

  for i in range(len(str)-1):
    if str[i] == '-' and str[i+1] in '0123456789.': str = str[:i] + '@m' + str[i+1:]
    if str[i] == '+' and str[i+1] in '0123456789.': str = str[:i] + '@p' + str[i+1:]

  tokens = re.split('([/\-\*\+\^ ])', str)
  tokens = [val for val in tokens if val != '' and val != ' ']

  for i in range(len(tokens)):
    tokens[i] = tokens[i].replace('@m', '-').replace('@p', '+')

  #print (f'{tokens}')

  for ix in range(1,len(tokens)):
    if ix + 2 > len(tokens): break
    if tokens[ix] != '+' and tokens[ix] != '-': continue
    if tokens[ix-1][-1].upper() != 'E': continue
    if not tokens[ix-1][:-1].replace('.', '', 1).isdigit(): continue
    if not tokens[ix+1].isdigit(): continue
    tokens = tokens[:ix-1] + [''.join(tokens[ix-1:ix+2])] + tokens[ix+2:]

  # 

  tokens, ixn = infixer(tokens, len(tokens))
  if return_list:
    return tokens
  else:
    return tokens[0]

#------------------------------------------------------------------
#------------------------------------------------------------------
# Convert from elegant parameter name to bmad parameter name.

def bmad_param(param, ele_name):
  global common, bmad_param_name

  if len(param) == 2 and param[0] == 'c' and param[1] in '123456':
    return f'tt{param[1:]}'

  if len(param) == 3 and param[0] == 'r' and param[1] in '123456' and param[2] in '123456':
    return f'tt{param[1:]}'

  if len(param) == 4 and param[0] == 't' and param[1] in '123456' and param[2] in '123456' and param[3] in '123456':
    return f'tt{param[1:]}'

  #

  if param not in bmad_param_name: return '?'
  bparam = bmad_param_name[param]

  if ele_name in common.ele_dict:
    bmad_type = common.ele_dict[ele_name].bmad_type
    elegant_type = common.ele_dict[ele_name].elegant_type
    if bparam == 'tilt' and (bmad_type == 'sbend' or bmad_type == 'rbend'): return 'ref_tilt'
    if param == 'l' and bmad_type == 'patch': return '?'
    if param == 'phase' and bmad_type != 'rfcavity' and bmad_type != 'lcavity': return '?'
    if param == 'b' and bmad_type[1:] != 'bend': return '?'
    if param == 'fse' and bmad_type[1:] != 'bend': return '?'
    if param == 'angle' and bmad_type[1:] != 'bend': return '?'
  return bparam

#------------------------------------------------------------------
#------------------------------------------------------------------
# Return dictionary of "A = value" parameter definitions.
# Also see action_parameter_dictionary routine.
# The difference is that this routine will not split on commas.

def parameter_dictionary(word_lst):

  orig_word_lst = word_lst

  # replace "0." or "0.0" with "0"
  word_lst = ['0' if x == '0.0' or x == '0.' else x for x in word_lst]

  pdict = OrderedDict()
  while True:
    if len(word_lst) == 0: return pdict

    if word_lst[1] != '=':
      print ('PROBLEM PARSING PARAMETER LIST: ' + ''.join(orig_word_lst))
      return pdict

    if '=' in word_lst[2:]:
      ix = word_lst.index('=', 2)
      if word_lst[0] == 'eyaw':
        pdict['yaw'] = ''.join(word_lst[2:ix-2])
      if word_lst[0] == 'epitch':
        pdict['pitch'] = ''.join(word_lst[2:ix-2])
      else:
        pdict[word_lst[0]] = ''.join(word_lst[2:ix-2])
      word_lst = word_lst[ix-1:]

    else:
      pdict[word_lst[0]] = ''.join(word_lst[2:])
      return pdict

#-------------------------------------------------------------------
#------------------------------------------------------------------
# Construct the bmad lattice file name

def bmad_file_name(elegant_file):
  if elegant_file.lower()[-4:] == '.lte':
    return elegant_file[:-4] + '.bmad'

  else:
    return elegant_file + '.bmad'

#------------------------------------------------------------------
#------------------------------------------------------------------

def wrap_write(line, f_out):
  MAXLEN = 120
  tab = ''
  line = line.rstrip()

  while True:
    if len(line) <= MAXLEN+1:
      f_out.write(tab + line + '\n')
      return

    if ',' in line[:MAXLEN]:
      ix = line[:MAXLEN].rfind(',')
      f_out.write(tab + line[:ix+1] + '\n')  # Don't need '&' after a comma

    else:
      for char in ' -+/*':
        if char in line[:MAXLEN]:
          ix = line[:MAXLEN].rfind(char)
          f_out.write(tab + line[:ix+1] + ' &\n')
          break

    tab = '         '
    line = line[ix+1:]

#------------------------------------------------------------------
#------------------------------------------------------------------
# Adds parenteses around expressions with '+' or '-' operators.
# Otherwise just returns the expression.
# Eg: '-1.2'  -> '-1.2'
#      '7+3'  -> '(7+3)'
#      '7*3'  -> '7*3'

def add_parens (str):
  state = 'begin'
  for ch in str:
    if ch in '0123456789.':
      if state == 'out' or state == 'begin': state = 'r1'

    elif ch == 'e':
      if state == 'r1':  state = 'r2'
      else:              state = 'out'

    elif ch in '-+':
      if state == 'r2':
        state = 'r3'
      elif state == 'begin':
        state = 'out'
      else:
        return '(' + str + ')'

    else:
      state = 'out'

  return str

#------------------------------------------------------------------
#------------------------------------------------------------------

def float_val (str, default):
  try:
    return float(str)
  except:
    return default

#------------------------------------------------------------------
#------------------------------------------------------------------

def int_val (str, default):
  try:
    return int(str)
  except:
    return default

#------------------------------------------------------------------
#------------------------------------------------------------------

def negate(str):
  str = add_parens(str)
  if str[0] == '-':
    return str[1:]
  elif str[0] == '+':
    return '-' + add_parens(str[1:])
  else:
    return '-' + add_parens(str)

#------------------------------------------------------------------
#------------------------------------------------------------------
# Parse a lattice element

def parse_element(dlist):
  global common, ele_type_translate

  ele = ele_struct(dlist[0])
  f_out = common.f_out[-1]

  found = False
  for elegant_type in ele_type_translate:
    if elegant_type.startswith(dlist[2]):
      ele.elegant_type = elegant_type
      ele.bmad_type = ele_type_translate[elegant_type]
      found = True
      if ele.elegant_type in problematical_translation_list: 
        print (f'NOTE: {dlist[2].upper()} TYPE ELEMENT IN ELEGANT LATTICE. TRANSLATION IS POTENTIALLY PROBLEMATICAL!')
      break

  if not found:
    print (f'{dlist[2].upper()} TYPE ELEMENT NOT FOUND IN TRANSLATION TABLE. WILL BE TRANSLATED TO A DRIFT!')
    ele.elegant_type = dlist[2]
    ele.bmad_type = 'drift'
  

  params = parameter_dictionary(dlist[4:])

  # If the reference energy changes then a cavity should be an lcavity.

  if ele.bmad_type == 'rfcavity' and params.get('change_p0', '0') != '0': ele.bmad_type = 'lcavity'
  if elegant_type == 'rotate':
    if params.get('exclude_floor', '0') == '0' and params.get('exclude_optics', '0') == '0':
      ele.bmad_type = 'patch'
    elif params.get('exclude_floor', '0') != '0':
      ele.bmad_type = 'taylor'
    elif params.get('exclude_optics', '0') != '0':
      ele.bmad_type = 'floor_shift'

  #

  ele.param = params
  common.ele_dict[dlist[0]] = ele

  line = f'{ele.name}: {ele.bmad_type}, type = "{elegant_type}"'

  for eparam in ele.param:
    ## if eparam in ['dx', 'dy', 'dz'] and 'etilt' in params: continue   # Handled later

    # Malign -> Gkicker

    if elegant_type == 'malign':
      if eparam == 'dx': 
        bparam = 'x_kick'
      elif eparam == 'dy':
        bparam = 'y_kick'
      elif eparam == 'dz':
        bparam = 'z_kick'
      elif eparam == 'dxp':
        bparam = 'px_kick'
      elif eparam == 'dyp':
        bparam = 'py_kick'
      elif eparam == 'dp':
        bparam = 'pz_kick'
      else:
        bparam = '?'

    else:
      bparam = bmad_param(eparam, ele.name)
      if bparam == '?': continue
      if ele.bmad_type == 'drift' and bparam != 'l': continue

    #

    value = postfix_to_infix(params[eparam])

    if bparam == 'phi0':
      if ele.bmad_type == 'lcavity':
        value = f'({value} - 90)/360'
      else:
        value = f'({value})/360'

    if bparam == 'pitch': value = negate(str)   # Corresponds to Bmad y_pitch

    if float_val(value, 1) == 0: continue
    line += f', {bparam} = {value}'

  # Below for parameters that do not have a standard translation
  # Note that the  Elegant rotate element has a "confusing" sign convention for the tilt parameter.

  if elegant_type == 'rotate' and ele.bmad_type == 'taylor' and params.get('tilt', '0') != '0':
    t = postfix_to_infix(params['tilt'])
    line = f'''{ele.name}: {ele.bmad_type}, type = "{elegant_type}", tt11 = cos({t}), tt13 = sin({t}),
t31 = -sin({t}), t33 = cos({t}), tt22 = cos({t}), tt24 = sin({t}), t42 = -sin({t}), t44 = cos({t})'''

  # FSE

  if 'fse' in ele.param and 'bend' == ele.bmad_type[1:]: line += f', dg = {postfix_to_infix(params["fse"])} * {ele.name}[angle]/{ele.name}[L]'
  if 'fse_dipole' in ele.param and 'bend' == ele.bmad_type[1:]: line += f', dg = {postfix_to_infix(params["fse_dipole"])} * {ele.name}[angle]/{ele.name}[L]'
  if 'charge'     in ele.param:
    wrap_write('parameter[n_part] = 1.602176634e-19', f_out)
    line += f', charge = {postfix_to_infix(params["charge"])}'

  if 'knl' in ele.param: line += f', k{params.get("order", "1")}l = {value}'

  # Etilt

  if 'etilt' in params and 'bend' == ele.bmad_type[1:]:
    value = postfix_to_infix(params['etilt'])
    if 'etilt_sign' in params and int_val(params['etilt_sign'], 1) == -1: value = negate(value)
    ang2 = add_parens(params.get('angle', '0')) + '/2'
    line += f', roll = {add_parens(value)} * cos({ang2})'

  # edge effects

  ee1 = int(params.get('edge1_effects', '1'))
  ee2 = int(params.get('edge2_effects', '1'))

  if ee1 == 0 and ee2 == 0:
    line += f', fringe_at = no_end'
  elif ee1 != 0 and ee2 == 0:
    line += f', fringe_at = entrance_end'
  elif ee1 == 0 and ee2 != 0:
    line += f', fringe_at = exit_end'

  wrap_write(line, f_out)

  return ele

#------------------------------------------------------------------
#------------------------------------------------------------------

def parse_command(command, dlist):
  global common

  f_out = common.f_out[-1]
  if common.debug: print (str(dlist))

  if len(dlist) == 0:
    f_out.write('\n')
    return

  # &run_setup namelist

  if dlist[0] == '&run_setup':
    params = namelist_dict(dlist)
    if 'p_central'     in params: wrap_write(f'parameter[p0c] = {params["p_central"]}', f_out)
    if 'p_central_mev' in params: wrap_write(f'parameter[p0c] = 1e6*({params["p_central_mev"]})', f_out)
    if 'use_beamline'  in params: 
      name = params["use_beamline"].replace('"', '').replace("'", '')
      wrap_write(f'use, {name}', f_out)
      common.beam_line_name = '##'  # Prevent printing of use statement of last defined line
    return

  # &bunched_beam namelist

  if dlist[0] == '&bunched_beam':
    params = namelist_dict(dlist)
    if 'beta_x'     in params: wrap_write(f'beginning[beta_a] = {params["beta_x"]}', f_out)
    if 'beta_y'     in params: wrap_write(f'beginning[beta_b] = {params["beta_y"]}', f_out)
    if 'alpha_x'    in params: wrap_write(f'beginning[alpha_a] = {params["alpha_x"]}', f_out)
    if 'alpha_y'    in params: wrap_write(f'beginning[alpha_b] = {params["alpha_y"]}', f_out)
    if 'eta_x'      in params: wrap_write(f'beginning[eta_x] = {params["eta_x"]}', f_out)
    if 'eta_y'      in params: wrap_write(f'beginning[eta_y] = {params["eta_y"]}', f_out)
    if 'etap_x'     in params: wrap_write(f'beginning[etap_x] = {params["etap_x"]}', f_out)
    if 'etap_y'     in params: wrap_write(f'beginning[etap_y] = {params["etap_y"]}', f_out)
    if 'p0'         in params: wrap_write(f'parameter[p0c] = {params["p0"]}', f_out)
    if 'emit_x'     in params: wrap_write(f'particle_start[emittance_a] = {params["emit_x"]}', f_out)
    if 'emit_y'     in params: wrap_write(f'particle_start[emittance_b] = {params["emit_y"]}', f_out)
    if 'emit_nx'    in params: wrap_write(f'particle_start[emittance_a] = {params["emit_nx"]}/parameter[p0c]', f_out)
    if 'emit_ny'    in params: wrap_write(f'particle_start[emittance_b] = {params["emit_ny"]}/parameter[p0c]', f_out)
    return

  # &bunched_beam namelist

  if dlist[0] == '&floor_coordinates':
    params = namelist_dict(dlist)
    if 'x0'         in params: wrap_write(f'beginning[x_position] = {params["x0"]}', f_out)
    if 'y0'         in params: wrap_write(f'beginning[y_position] = {params["y0"]}', f_out)
    if 'z0'         in params: wrap_write(f'beginning[z_position] = {params["z0"]}', f_out)
    if 'theta0'     in params: wrap_write(f'beginning[theta_position] = {params["theta0"]}', f_out)
    if 'phi0'       in params: wrap_write(f'beginning[phi_position] = {params["phi0"]}', f_out)
    if 'psi0'       in params: wrap_write(f'beginning[psi_position] = {params["psi0"]}', f_out)
    return

  # Ignore other Namelists

  if dlist[0][0] == '&':
    return

  # "% <expression> sto <var>" construct

  if dlist[0] == '%':
    toks = postfix_to_infix(dlist[1], True)
    if len(toks) != 3 or toks[1] != 'sto':
      print (f'MALFORMED CONSTANT DEFINITION: {command}')
      return
    wrap_write (f'{toks[2]} = {toks[0]}', f_out)
    return

  # "#include"

  if dlist[0] == '#include:':

    file = dlist[1]
    if '"' in file or "'" in file:
      file = file.replace('"', '').replace("'", '')
    else:
      file = file.lower()    

    common.f_in.append(open(file, 'r'))  # Store file handle
    if common.one_file: 
      f_out.write(f'\n! In File: {common.f_in[-1].name}\n')
    else:
      f_out.write(f'call, file = {bmad_file_name(file)}\n')
      common.f_out.append(open(bmad_file_name(file), 'w'))
    return

  # Everything below has at least 3 words

  if len(dlist) < 3:
    print ('Unknown construct:\n  ' + command.strip())
    return

  # Is there a colon or equal sign?
  try:
    ix_colon = dlist.index(':')
  except:
    ix_colon = -1

  try:
    ix_equal = dlist.index('=')
  except:
    ix_equal = -1

  # Line

  if ix_colon > 0 and dlist[ix_colon+1] == 'line':
    if '-' in command: print('''
WARNING! LINE DEFINITION USES ELEMENT REVERSAL (NEGATIVE SIGN CHARACTER DETECTED).
IF REVERSAL INVOLVES A BEND WITH DIFFERING E1 AND E2 FACE ANGLES, THE TRANSLATION WILL BE OFF
SINCE WITH BMAD (AND MAD FOR THAT MATTER) REVERSAL DOES NOT FLIP E1 AND E2 BUT WITH ELEGANT IT DOES.
THAT IS, YOU WILL NEED TO EDIT THE BMAD LATTICE FILE TO FIX.''')
    wrap_write(command.replace(' ,', ','), f_out)
    if common.beam_line_name != '##': common.beam_line_name = dlist[0]
    return

  # Var definition.

  if dlist[1] == '%':
    # ...
    return

  # Ele parameter set

  if dlist[0] in common.ele_dict and dlist[1] == ',' and dlist[3] == '=':
    value = command.split('=')[1].strip()
    if dlist[0] in common.ele_dict:
      ele_name = common.ele_dict[dlist[0]].name
      name = f'{dlist[0]}[{bmad_param(dlist[2], ele_name)}]'
    else:  # In a complete valid lattice, parameter seets always happen after the element has been defined
      name = f'{dlist[0]}[{bmad_param(dlist[2], "???")}]'
    f_out.write(f'{name} = {value}\n')
    return

  # Element def

  if dlist[1] == ':':
    parse_element(dlist)
    return

  # Unknown

  print (f"Unknown construct:\n" + command + '\n')

#------------------------------------------------------------------
#------------------------------------------------------------------
# Get next Elegant command.
# Read in Elegant file line-by-line.  Assemble lines into commands.

def get_next_command ():
  global common

  quote_delim = ''  # Quote mark delimiting a string. Blank means not parsing a string yet.
  command = ''
  dlist = []

  # Loop until a command has been found

  while True:

    # Get a line

    if common.command == '':
      while True:
        f_in = common.f_in[-1]
        f_out = common.f_out[-1]

        line = f_in.readline()
        if len(line) > 0: break    # Check for end of file
        
        common.f_in[-1].close()
        common.f_in.pop()          # Remove last file handle
        if len(common.f_in) == 0: return ['', dlist]

        if not common.one_file:
          common.f_out[-1].close()
          common.f_out.pop()       # Remove last file handle

    else:
      f_in = common.f_in[-1]
      f_out = common.f_out[-1]
      line = common.command
      common.command = ''

    # Parse line

    if line.lstrip().startswith('!!verbatim'):
      f_out.write(line[10:].strip() + '\n')
      continue

    try:
      ix = line.index('!')
      f_out.write(line[ix:])
      line = line[:ix]
      if line.strip() == '': continue
    except:
      if line.strip() == '': 
        f_out.write('\n')
        continue

    try:
      ix = line.index(';')
      line = line[:ix]
      common.command = line[ix+1:].strip()
    except:
      if line.strip() == '': continue

    if line[0] == '%':
      command = line
      dlist = ['%', line[1:].strip()]
      return [command, dlist]

    if line.rstrip()[:9].lower() == '#include:':
      command = line
      dlist = ['#include:', line[9:].strip()]
      return [command, dlist]

    if line[0] == '&':    # Namelist
      while True:
        command += line
        dlist += [val for val in re.split('([=, ])', line.strip()) if val != '' and val != ' ']
        if '&end' in line: return [command, dlist]
        line = f_in.readline()

    while line != '':
      for ix in range(len(line)):
        #print (f'Ix: {ix}/{len(line)} "{line[ix]}" -{quote_delim}-|{line.rstrip()}')
        #print (f'  C: {command}')
        #print (f'  D: {dlist}')

        if line[ix] == '"' or line[ix] == "'":
          if line[ix] == quote_delim:      # Found end of string
            command += quote_delim + line[:ix+1]
            dlist.append(quote_delim + line[:ix+1])
            line = line[ix+1:]
            quote_delim = ''
            break

          elif quote_delim == '':          # Found start of string
            quote_delim = line[ix]
            command += line[:ix]
            if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
            line = line[ix+1:]
            break

        if quote_delim != '': continue     # Cycle if in quote string

        if line[ix] == '&' or (line[ix] == ',' and ix == len(line.rstrip())-1):
          line2 = f_in.readline().lstrip()
          while True:
            if line2[0] == '\n' or line2.lstrip()[0] == '!':
              f_out.write(line2)
              line2 = f_in.readline().lstrip()
            else:
              break

          if line[ix] == '&':
            line = line[:ix] + line2
          else:
            line = line[:ix+1] + line2

          break

        elif line[ix] in ':,=':
          command += line[:ix+1]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          dlist.append(line[ix])
          line = line[ix+1:]
          break

        elif line[ix] == '\n' or ix == len(line)-1:
          command += line
          if line.strip() != '': dlist.append(line.strip().lower())
          return [command, dlist]

        elif ix == len(line)-1:   # Happens at end of file
          command += line[:ix+1]
          if line.strip() != '': dlist.append(line.strip().lower())
          return [command, dlist]

#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
# Main program.

start_time = time.time()

# Read the parameter file specifying the Elegant lattice file, etc.

argp = argparse.ArgumentParser()
argp.add_argument('elegant_files', help = 'Name of input Elegant lattice file', nargs='+')
argp.add_argument('-d', '--debug', help = 'Print debug info (not of general interest).', action = 'store_true')
argp.add_argument('-f', '--many_files', help = 'Create a Bmad file for each Elegant input file.', action = 'store_true')
argp.add_argument('-c', '--constants', help = 'Add to lattice file a list of Elegant defined constants.', action = 'store_true')
arg = argp.parse_args()
## print(arg)

common = common_struct()
common.debug = arg.debug
common.one_file = not arg.many_files
common.add_constants = arg.constants

print ('*******Note: In beta testing! Please report any problems! **********')
print (f'Input lattice file(s) are: {arg.elegant_files}')

# Open files for reading and writing

# Loop over all input files

for ixf, elegant_lattice_file in enumerate(arg.elegant_files):
  common.f_in = [open(elegant_lattice_file, 'r')]

  if ixf == 0 or not common.one_file:
    bmad_lattice_file = bmad_file_name(elegant_lattice_file)
    print (f'Output lattice file: {bmad_lattice_file}')
    f_out = open(bmad_lattice_file, 'w')
    common.f_out = [f_out]

  common.command = ''  # init

  f_out.write (f'''
!+
! Translated by elegant_to_bmad.py from Elegant file(s): {arg.elegant_files}
!-

''')

  if common.add_constants and ixf == 0:
    f_out.write (f'''
c_cgs = 2.99792458e10
c_mks = 2.99792458e8 
e_cgs = 4.80325e-10
e_mks = 1.60217733e-19
me_cgs = 9.1093897e-28
me_mks = 9.1093897e-31
re_cgs = 2.81794092e-13
re_mks = 2.81794092e-15
kb_cgs = 1.380658e-16
kb_mks = 1.380658e-23
mev = 0.51099906
hbar_mks = 1.0545887e-34
hbar_MeVs = 6.582173e-22
mp_mks = 1.6726485e-27
mu_o = 4 * pi * 1e-7
eps_o = 1e7 / (4 * pi * c_mks^2)
Kas = 191.655e-2
Kaq = 75.0499e-2

''')

  # parse, convert and output elegant commands

  while True:
    [command, dlist] = get_next_command()
    if len(common.f_in) == 0: break
    parse_command(command, dlist)
    if len(common.f_in) == 0: break   # Hit Quit/Exit/Stop statement.

  #------------------------------------------------------------------
  f_out = common.f_out[0]  # Should be only one left
  if common.beam_line_name != '' and common.beam_line_name != '##': 
    f_out.write(f'\nuse, {common.beam_line_name}\n')
    common.beam_line_name = ''

#

print ('*******Note: In beta testing! Please report any problems! **********')

