#!/usr/bin/env python

#+
# Script to convert from MADX lattice format to Bmad lattice format.
# See the README file for more details
#-

import sys, re, math, argparse, time
from collections import OrderedDict

if sys.version_info[0] < 3 or sys.version_info[1] < 6:
  raise Exception("Must be using Python 3.6+")

#------------------------------------------------------------------

class ele_struct:
  def __init__(self, name):
    self.name = name
    self.madx_inherit = ''      # Is another element name or quadrupole, etc.
    self.madx_base_type = ''    # Is quadrupole, etc.
    self.bmad_inherit = ''      # Is another element name or quadrupole, etc.
    self.bmad_base_type = ''    # Is quadrupole, etc.
    self.at = '0'               # Used if element is in a sequence
    self.from_ref_ele = ''      # Used if element is in a sequence
    self.param = OrderedDict()
    self.count = 0

class seq_struct:
  def __init__(self, name = ''):
    self.name = name
    self.l = '0'
    self.refer = 'centre'
    self.refpos = ''
    self.seq_ele_dict = OrderedDict()
    self.last_ele_offset = ''
    self.line = ''                   # For when turning a sequence into a line
    self.drift_list = []

class common_struct:
  def __init__(self):
    self.debug = False               # Command line argument.
    self.prepend_vars = True         # Command line argument.
    self.superimpose_eles = False    # Command line argument.
    self.one_file = True             # Command line argument.
    self.in_seq = False              # Inside a sequence/endsequence construct?
    self.in_track = False            # Inside a track/endtrack construct?
    self.in_match = False            # Inside a match/endmatch construct?
    self.seqedit_name = ''           # Name of sequence in seqedit construct.
    self.last_seq = seq_struct()     # Current sequence being parsed.
    self.seq_dict = OrderedDict()    # List of all sequences.
    self.ele_dict = {}               # Dict of elements
    self.var_def_list = []           # List of "A = B" sets after translation to Bmad. Does not Include "A->P = B" parameter sets.
    self.var_name_list = []          # List of madx variable names.
    self.super_list = []             # List of superimpose statements to be prepended to the bmad file.
    self.f_in = []         # MADX input files
    self.f_out = []        # Bmad output files
    self.use = ''
    self.command = ''    # Scratch storage for read_madx_command routine.
    self.drift_count = 0

#------------------------------------------------------------------
#------------------------------------------------------------------

ele_param_factor = {
  'volt':     ' * 1e-6',
  'freq':     ' * 1e-6',
  'energy':   ' * 1e-9',
  'ex':       ' * 1e-6',
  'ey':       ' * 1e-6',
  'pc':       ' * 1e-9',
  'lag':      ' + 0.5',
}

ele_inv_param_factor = {
  'volt':     ' * 1e6',
  'freq':     ' * 1e6',
  'energy':   ' * 1e9',
  'ex':       ' * 1e6',
  'ey':       ' * 1e6',
  'pc':       ' * 1e9',
  'lag':      ' - 0.5',
}

negate_param = ['lag']

const_trans = {
  'e':       'e_log',
  'nmass':   'm_neutron * 1e9',
  'mumass':  'm_muon * 1e9',
  'clight':  'c_light',
  'qelect':  'e_charge',
  'hbar':    'h_bar * 1e6',
  'erad':    'r_e',
  'prad':    'r_p',
  'ceil':    'ceiling',
  'round':   'nint',
  'ranf':    'ran',
  'gauss':   'ran_gauss',
}

sequence_refer = {
  'entry':  'beginning',
  'centre': 'center',
  'exit':   'end'
}

ele_type_translate = {
    'tkicker':      'kicker', 
    'hacdipole':    'ac_kicker',
    'placeholder':  'instrument',
    'matrix':       'taylor', 
    'srotation':    'patch',
    'xrotation':    'patch',
    'yrotation':    'patch',
    'translation':  'patch',
    'changeref':    'patch',
    'monitor':      'monitor',
    'hmonitor':     'monitor',
    'vmonitor':     'monitor',
    'marker':       'marker',
    'drift':        'drift',
    'sbend':        'sbend',
    'rbend':        'rbend',
    'quadrupole':   'quadrupole',
    'sextupole':    'sextupole',
    'octupole':     'octupole',
    'multipole':    'multipole',
    'solenoid':     'solenoid',
    'hkicker':      'hkicker',
    'vkicker':      'vkicker',
    'tkiker':       'kicker',
    'kicker':       'kicker',
    'rfcavity':     'rfcavity',
    'twcavity':     'lcavity',
    'elseparator':  'elseparator',
    'instrument':   'instrument',
    'ecollimator':  'ecollimator',
    'rcollimator':  'rcollimator',
    'collimator':   'collimator',  # Will get sorted out afterwards
    'beambeam':     'beambeam',
    'crabcavity':   'crab_cavity',
    'hacdipole':    'ac_kicker',
    'vacdipole':    'ac_kicker',
    'rfmultipole':  '???',
    'nllens':       '???',
    'dipedge':      '???',
    'sequence':     '???',
    'twiss':        '???',
    'beam':         '???',
}

ignore_madx_param = ['lrad', 'slot_id', 'aper_tol', 'apertype', 'thick', 'add_angle', 'assembly_id', 
                     'mech_sep', 'betrf', 'tfill', 'shunt', 'pg', 'model']
ignore_madx_ele_param = {}
ignore_madx_ele_param['beambeam'] = ['particle', 'pc']

bmad_param_name = {
    'volt':   'voltage',
    'freq':   'rf_frequency',
    'lag':    'phi0',
    'ex':     'e_field',
    'ey':     'e_field',
    'lrad':   'l',  
    'xsize':  'x_limit',
    'ysize':  'y_limit',
    'dx':     'x_offset',
    'dy':     'y_offset',
    'ds':     'z_offset',
}

#------------------------------------------------------------------
#------------------------------------------------------------------
# Is character a valid character to be used in a label?

def is_label_char(char):
  return char.isalnum() or char in '._'

#------------------------------------------------------------------
#------------------------------------------------------------------
# Order var defs so that vars that depend upon other vars are come later.
# Also comment out first occurances if there are multiple defs of the same var.

def order_var_def_list():

  # Mark duplicates
  dependent_list = {}   # Stores names and dependencies
  new_def_list = []

  for vdef in reversed(common.var_def_list):
    if vdef[0] in dependent_list:
      new_def_list.insert(0, ['! Duplicate: ' + vdef[0], vdef[1]])
    else:
      new_def_list.insert(0, vdef)
      exp_list = re.split('\+|-|\*|/|\(|\)|\^|,', vdef[1])
      dependent_list[vdef[0]] = list(x.split() for x in exp_list)

  common.var_def_list = new_def_list

  # Move vars that are dependent upon vars defined further up the list.
  new_def_list = common.var_def_list[:]

  ix = 0
  while ix < len(new_def_list):
    vdef = new_def_list[ix]

    if vdef[0][0] == '!':
      ix += 1
      continue

    moved = False
    for ix2 in range(len(new_def_list)-1, ix, -1):
      sub = new_def_list[ix2][0]
      str = vdef[1]
      ns = len(str)
      found = False

      for match in re.finditer(sub, str):
        j = match.start()
        k = match.end()
        if (j == 0 or not is_label_char(str[j-1])) and (k > len(str)-1 or not is_label_char(str[k])): found = True

      if found:
        new_def_list.pop(ix)
        new_def_list.insert(ix2, vdef)
        moved = True
        break

    if not moved: ix += 1

  common.var_def_list = new_def_list

#------------------------------------------------------------------
#------------------------------------------------------------------
# Is an expression zero (to within 1e-11)?

def is_zero(input):
  if isinstance(input, str):
    try:
      v = eval(input)
      return abs(v) < 1e-11
    except:
      return False

  else:
    return abs(v) < 1e-11
  

#------------------------------------------------------------------
#------------------------------------------------------------------
# Convert from madx parameter name to bmad parameter name.

def bmad_param(param, ele_name):
  global bmad_param_name

  if ele_name in common.ele_dict:
    madx_type = common.ele_dict[ele_name].madx_base_type
  else:
    madx_type = 'xxxx'

  if param == 'tilt':
    if madx_type == 'sbend' or madx_type == 'rbend':
      return 'ref_tilt'
    else:
      return 'tilt'

  elif len(param) == 5 and param[0:4] == 'kick' and param[4] in '123456':
    return f'tt{param[4]}'

  elif len(param) == 4 and param[0:2] == 'rm' and param[2] in '123456' and param[3] in '123456':
    return f'tt{param[2:]}'

  elif len(param) == 5 and param[0:2] == 'tm' and param[2] in '123456' and param[3] in '123456' and param[4] in '123456':
    return f'tt{param[2:]}'

  # Translate something like "k1s" to "k1" in the hopes that k1 on the MADX side is zero so the 
  # translation is k1s -> k1, tilt -> tilt + pi/2
  elif len(param) == 3 and param[0] == 'k' and param[1].isdigit() and param[2] == 's':
    return f'{param[:-1]}'

  elif param in bmad_param_name:
    return bmad_param_name[param]

  else:
    return param

#------------------------------------------------------------------
#------------------------------------------------------------------
# Return dictionary of "A = value" parameter definitions.

def parameter_dictionary(word_lst):
  madx_logical = ['kill_ent_fringe', 'kill_exi_fringe', 'thick', 'no_cavity_totalpath', 'chrom']

  # Remove :, {, and } chars for something like "kn := {a, b, c}"
  word_lst = list(filter(lambda a: a not in ['{', '}', ':'], word_lst))

  # replace "0." or "0.0" with "0"
  word_lst = ['0' if x == '0.0' or x == '0.' else x for x in word_lst]

  # Logical can be of the form: "<logical_name>" or "-<logical_name>".
  # Put this into dict

  pdict = OrderedDict()

  for logical in madx_logical:
    if logical in word_lst:
      ix = word_lst.index(logical)
      if ix+1 < len(word_lst) and word_lst[ix+1] == '=': continue
      pdict[logical] = 'true'
    elif '-' + logical in word_lst:
      ix = word_lst.index('-' + logical)
      if ix+1 < len(word_lst) and word_lst[ix+1] == '=': continue
      pdict[logical] = 'false'
    else:
      continue

    word_lst.pop(ix)
    if ix < len(word_lst): word_lst.pop(ix)   # Must be comma


  # Fill dict
  while True:
    if len(word_lst) == 0: return pdict

    if word_lst[1] != '=':
      print ('PROBLEM PARSING PARAMETER LIST: ' + ''.join(word_lst))
      return pdict

    if '=' in word_lst[2:]:
      ix = word_lst.index('=', 2)
      pdict[word_lst[0]] = ''.join(word_lst[2:ix-2])
      word_lst = word_lst[ix-1:]

    else:
      pdict[word_lst[0]] = ''.join(word_lst[2:])
      return pdict

#------------------------------------------------------------------
#------------------------------------------------------------------
# Convert expression from MADX format to Bmad format
# To convert <expression> a construct that look like "<target_param> = <expression>".

def bmad_expression(line, target_param):
  global const_trans, ele_param_factor, negate_param, ele_inv_param_factor

  # Remove {, and } chars for something like "kn := {a, b, c}". Also remove leading and ending quote marks
  line = line.replace('{', '').replace('}', '').strip('"\'')

  #

  lst = re.split(r'(,|-|\+|\(|\)|\>|\*|/|\^)', line)
  lst = list(filter(lambda a: a != '', lst))    # Remove blank. EG: "->" => ["-", "", ">"] => ["-", ">"]

  out = ''

  while len(lst) != 0:
    if len(lst) >= 4 and lst[1] == '-' and lst[2] =='>':
      if lst[3] in ele_param_factor:
        if (len(lst) >= 5 and lst[4] == '^') or (len(out.strip()) > 0 and out.strip()[-1] == '/'):
          out += '(' + lst[0] + '[' + bmad_param(lst[3].strip(), lst[0]) + ']' + ele_param_factor[lst[3]]
        else:
          out += lst[0] + '[' + bmad_param(lst[3].strip(), lst[0]) + ']' + ele_param_factor[lst[3]]
      else:
        out += lst[0] + '[' + bmad_param(lst[3].strip(), lst[0]) + ']'
      lst = lst[4:]

    elif lst[0] in const_trans:
      out += const_trans[lst.pop(0)]

    else:
      out += lst.pop(0)

  # End while

  if target_param in ele_inv_param_factor: 
    if target_param in negate_param:
      out = '-' + add_parens(out, True) + ele_inv_param_factor[target_param]
    else:
      out = add_parens(out, True) + ele_inv_param_factor[target_param]


  return out

#-------------------------------------------------------------------
#------------------------------------------------------------------
# Construct the bmad lattice file name

def bmad_file_name(madx_file):

  if madx_file.lower()[-5:] == '.madx':
    return madx_file[:-5] + '.bmad'
  elif madx_file.lower()[-4:] == '.mad':
    return madx_file[:-4] + '.bmad'
  elif madx_file.lower()[-4:] == '.seq':
    return madx_file[:-4] + '.bmad'
  else:
    return madx_file + '.bmad'

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
# Eg: '-1.2'  -> '-1.2'    If ignore_leading_pm = True
# Eg: '-1.2'  -> '(-1.2)'  If ignore_leading_pm = False
#      '7+3'  -> '(7+3)'
#      '7*3'  -> '7*3'
# Note: Need to ignore +/- sybols in something like "3e-4"

def add_parens (str, ignore_leading_pm):
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
      elif state == 'begin' and ignore_leading_pm:
        state = 'out'
      else:
        return '(' + str + ')'

    else:
      state = 'out'

  return str

#------------------------------------------------------------------
#------------------------------------------------------------------

def negate(str):
  str = add_parens(str, True)
  if str[0] == '-':
    return str[1:]
  elif str[0] == '+':
    return '-' + str[1:]
  else:
    return '-' + str

#------------------------------------------------------------------
#------------------------------------------------------------------
# Parse a lattice element
# Assumed to be of the form dlist = ["name", ":", "type", ",", ...]

def parse_and_write_element(dlist, write_to_file, command):
  global common, ele_type_translate, ignore_madx_param

  if dlist[2] == 'dipedge':
    print ('DIPEDGE ELEMENT NOT TRANSLATED. SUGGESTION: MODIFY THE LATTICE FILE AND MERGE THE DIPEDGE ELEMENT WITH THE NEIGHBORING BEND.')
    return

  if dlist[2] in common.ele_dict:
    ele = ele_struct(dlist[0])
    ele.madx_inherit = dlist[2]
    ele.madx_base_type = common.ele_dict[dlist[2]].madx_base_type
    ele.bmad_inherit = dlist[2]
    ele.bmad_base_type = common.ele_dict[dlist[2]].bmad_base_type

  else:
    found = False
    for madx_type in ele_type_translate:
      if madx_type.startswith(dlist[2]): 
        ele = ele_struct(dlist[0])
        ele.madx_inherit = madx_type
        ele.madx_base_type = madx_type
        ele.bmad_inherit = ele_type_translate[madx_type]
        ele.bmad_base_type = ele.bmad_inherit
        found = True
        break

    if not found:
      print (dlist[2].upper() + ' TYPE ELEMENT IS UNKNOWN!')
      return
  #End if

  if ele.madx_base_type == '???':
    print (dlist[2].upper() + ' TYPE ELEMENT CANNOT BE TRANSLATED TO BMAD.')
    return

  params = parameter_dictionary(dlist[4:])

  if ele.madx_base_type == 'elseparator':
    if 'ex' in params:
      if 'ey' in params:
        if 'tilt' in params:
          params['tilt'] = params['tilt'] + ' - atan2(' + params['ex'] + ', ' + params['ey'] + ')'
        else:
          params['tilt'] = '-atan2(' + params['ex'] + ', ' + params['ey'] + ')'
        params['ey'] = 'sqrt((' + params['ex'] + ')^2 + (' + params['ey'] + ')^2)'
      else:
        if 'tilt' in params:
          params['tilt'] = params['tilt'] + ' - pi/2'
        else:
          params['tilt'] = '-pi/2'
        params['ey'] = params['ex']

  elif ele.madx_base_type == 'xrotation':
    if 'angle' in params: params['y_pitch'] = negate(params.pop('angle'))

  elif ele.madx_base_type == 'yrotation':
    if 'angle' in params: params['x_pitch'] = negate(params.pop('angle'))

  elif ele.madx_base_type == 'srotation':
    if 'angle' in params: params['tilt'] = params.pop('angle')

  elif ele.madx_base_type == 'changeref':
    if 'patch_ang' in params:
      angles = params.pop('patch_ang').split(',')
      params['y_pitch'] = angles[0]
      params['x_pitch'] = negate(angles[1])
      params['tilt']   = angles[2]
    if 'patch_trans' in params:
      trans = params.pop('patch_trans').split(',')
      params['x_offset'] = trans[0]
      params['y_offset'] = trans[1]
      params['z_offset'] = trans[2]

  elif ele.madx_base_type == 'rbend' or ele.madx_base_type == 'sbend':
    if 'tilt' in params: params['ref_tilt'] = params.pop('tilt')
    kill_ent = False; kill_exi = False
    if 'kill_ent_fringe' in params: kill_ent = (params['kill_ent_fringe'] == 'true')
    if 'kill_exi_fringe' in params: kill_exi = (params['kill_exi_fringe'] == 'true')
    if 'k0' in params: params['dg'] = params.pop('k0')
    if 'k0s' in params: params['a0'] = params.pop('k0s') + ' * ' + ele.name + '[angle]'
    if 'ktap' in params: params['dg'] = params.pop('ktap') + ' * ' + ele.name + '[angle] / ' + ele.name + '[l]'

    if kill_ent and kill_exi:
      params['fringe_at'] = 'no_end'
    elif kill_exi:
      params['fringe_at'] = 'entrance_end'
    elif kill_ent:
      params['fringe_at'] = 'exit_end'

  elif ele.madx_base_type == 'quadrupole':
    if 'k1' in params and 'k1s' in params:
      if 'tilt' in params:
        params['tilt'] = params['tilt'] + ' - atan2(' + params['k1s'] + ', ' + params['k1'] + ')/2'
      else:
        params['tilt'] = '-atan2(' + params['k1s'] + ', ' + params['k1'] + ')/2'
      params['k1'] = 'sqrt((' + params['k1'] + ')^2 + (' + params['k1s'] + ')^2)'
      params.pop('k1s')
    elif 'k1s' in params:
      if 'tilt' in params:
        params['tilt'] = params['tilt'] + ' - pi/4'
      else:
        params['tilt'] = '-pi/4'
      params.pop('k1s')

    # Need to differentiate:
    #   q: q0, ktap = 0.001;                 ! Here q inherits k1 from q0
    #   q: quad, k1 = XXX, ktap = 0.001;

    if 'ktap' in params: 
      if 'k1' in params:
        params['k1'] = params['k1'] + ' * (1 + ' + params.pop('ktap') + ')' 
      else:
        params['k1'] = ele.name + '[k1] * (1 + ' + params.pop('ktap') + ')' 

  elif ele.madx_base_type == 'sextupole':
    if 'k2' in params and 'k2s' in params:
      if 'tilt' in params:
        params['tilt'] = params['tilt'] + ' - atan2(' + params['k2s'] + ', ' + params['k2'] + ')/3'
      else:
        params['tilt'] = '-atan2(' + params['k2s'] + ', ' + params['k2'] + ')/3'
      params['k2'] = 'sqrt((' + params['k2'] + ')^2 + (' + params['k2s'] + ')^2)'
      params.pop('k2s')
    elif 'k2s' in params:
      if 'tilt' in params:
        params['tilt'] = params['tilt'] + ' - pi/6'
      else:
        params['tilt'] = '-pi/6'
      params.pop('k2s')

    if 'ktap' in params:
      if 'k2' in params:
        params['k2'] = params['k2'] + ' * (1 + ' + params.pop('ktap') + ')' 
      else:
        params['k2'] = ele.name + '[k2] * (1 + ' + params.pop('ktap') + ')' 


  elif ele.madx_base_type == 'octupole':
    if 'k3' in params and 'k3s' in params:
      if 'tilt' in params:
        params['tilt'] = params['tilt'] + ' - atan2(' + params['k3s'] + ', ' + params['k3'] + ')/4'
      else:
        params['tilt'] = '-atan2(' + params['k3s'] + ', ' + params['k3'] + ')/4'
      params['k3'] = 'sqrt((' + params['k3'] + ')^2 + (' + params['k3s'] + ')^2)'
      params.pop('k3s')
    elif 'k3s' in params:
      if 'tilt' in params:
        params['tilt'] = params['tilt'] + ' - pi/8'
      else:
        params['tilt'] = '-pi/8'
      params.pop('k3s')

  elif ele.madx_base_type == 'multipole':
    if 'knl' in params:
      for n, knl in enumerate(params.pop('knl').split(',')): 
        if knl == '0': continue
        params['k' + str(n) + 'l'] = bmad_expression(knl, '')
    if 'ksl' in params:
      for n, ksl in enumerate(params.pop('ksl').split(',')):  
        if ksl == '0': continue
        params['k' + str(n) + 'sl'] = bmad_expression(ksl, '')


  elif ele.madx_base_type == 'collimator':
    if params['apertype'] in ['ellipse', 'circle']:
      ele.bmad_inherit = 'ecollimator'
    else:
      ele.bmad_inherit = 'rcollimator'

  elif ele.madx_base_type == 'beambeam':
    if 'npart' in params:
      f_out = common.f_out[-1]
      f_out.write(f"parameter[n_part] = {params['npart']}\n")
      params.pop('npart')

  # collimator conversion

  if 'apertype' in params:
    aperture = bmad_expression(params.pop('aperture').replace('{', '').replace('}', ''), '')
    [params['x_limit'], params['y_limit']] = aperture.split(',')[2:4]

    if params['apertype'] in ['ellipse', 'circle']:
      params['aperture_type'] = 'elliptical'
    else:
      params['aperture_type'] = 'rectangular'

  if 'aper_offset' in params:
    params['x_offset'] = params['aper_offset'].split(',')[0]
    params['y_offset'] = params['aper_offset'].split(',')[1]

  # sequence params

  if 'at' in params: ele.at = params.pop('at')
  if 'from' in params: ele.from_ref_ele = params.pop('from')

  # Can happen in sequences that there are "name: name, at = X" constructs where name has already be defined.

  ele.param = params

  if write_to_file:
    line = ele.name + ': ' + ele.bmad_inherit
    for param in ele.param:
      if param in ignore_madx_param: continue
      if ele.madx_base_type in ignore_madx_ele_param and param in ignore_madx_ele_param[ele.madx_base_type]: continue
      line += ', ' + bmad_param(param, ele.name) + ' = ' + bmad_expression(params[param], param)
    f_out = common.f_out[-1]
    # Can have situation where an element is defined outside of a sequence ("this_name: that_class") and
    # inside of the sequence get the same definition.
    if dlist[0] in common.ele_dict:
      if ele.param != common.ele_dict[ele.name].param:
        print (f'''ERROR: ELEMENT WITH NAME {ele.name} IS BEING REDEFINED. THIS MAY LEAD TO PROBLEMS.''')
        wrap_write('!! Element redefined: ' + line, f_out)
    else:
      wrap_write(line, f_out)

  if dlist[0] not in common.ele_dict: common.ele_dict[dlist[0]] = ele

  return ele

#------------------------------------------------------------------
#------------------------------------------------------------------
# The "command" arg is the unsplit madx command.
# The "dlist" arg is the command split into pieces and converted to lower case.

def parse_command(command, dlist):
  global common, sequence_refer

  f_out = common.f_out[-1]

  if len(dlist) == 0: 
    f_out.write('\n')
    return

  # Get rid of "real", "int", "const" "const real", etc. prefixes

  if dlist[0].startswith('real ') or dlist[0].startswith('int ') or dlist[0].startswith('const '): dlist[0] = dlist[0].split(' ', 1)[1].strip()
  if dlist[0].startswith('real ') or dlist[0].startswith('int ') or dlist[0].startswith('const '): dlist[0] = dlist[0].split(' ', 1)[1].strip()
  if dlist[0].startswith('shared '): dlist[0] = dlist[0].split(' ', 1)[1].strip()

  # Transform: "a := 3" -> "a = 3"

  for ix in range(len(dlist)):
    if ix > len(dlist)-2: break
    if dlist[ix] == ':' and dlist[ix+1] == '=': dlist.pop(ix)

  # track and match constructs

  if dlist[0] == 'match': 
    common.in_match = True
    print ('Ignoring match construct: ' + command)
    return

  if dlist[0] == 'track': 
    common.in_track = True
    print ('Ignoring track construct: ' + command)
    return

  if dlist[0] == 'endmatch': 
    common.in_match = False
    return

  if dlist[0] == 'endtrack': 
    common.in_track = False
    return

  if common.in_match or common.in_track: return

  # The MADX parser takes a somewhat cavilier attitude towards commas and sometimes they can be omitted.
  # Examples: "call file" instead of  "call, file" and "q: quadrupole l = 7" instead of "q: quadrupole, l = 7".
  # So put the comma back in to make things uniform for easier parsing.
  # But do not do this with arithmetical expressions and quoted strings.

  if common.debug: print (str(dlist))

  i = 0
  while i < len(dlist):
    if ' ' in dlist[i] and not any(char in dlist[i] for char in '"\'+-*/^'):
      split = dlist[i].split()
      dlist = dlist[:i] + [split[0], ',', split[1]] + dlist[i+1:]
    i += 1

  # Ignore the following.

  if dlist[0] in ['exec' 'while', 'if']: 
    print (f'ERROR: "{dlist[0]}" COMMAND IGNORED: {command}\n' +
            '  THIS MEANS THAT IT IS LIKELY THAT THE BMAD LATTICE WILL BE DIFFERENT FROM THE MADX LATTICE!')
    return

  if dlist[0] in ['aperture', 'show', 'value', 'efcomp', 'print', 'select', 'optics', 'option', 'survey',
                  'emit', 'help', 'set', 'eoption', 'system', 'ealign', 'sixtrack', 'flatten', 
                  'elseif', 'else', 'savebeta', 'exec', 'makethin', 'save']:
    print ('Note! Ignoring command: ' + command)
    return

  if dlist[0] in ['twiss'] and len(dlist) == 1:
    print ('Note! Ignoring command: ' + command)
    return

  if 'macro' in dlist:
    return

  # Flag this

  if dlist[0] in ['cycle', 'reflect', 'move', 'remove', 'replace', 'extract']:
    print (f'WARNING! CANNOT TRANSLATE THE COMMAND: {dlist[0].upper()}')
    return

  # Seqedit

  if dlist[0] == 'seqedit':
    common.seqedit_name = dlist[4]
    return

  if dlist[0] == 'endedit':
    common.seqedit_name = ''
    return

  # Install

  if dlist[0] == 'install':
    params = parameter_dictionary(dlist[2:])
    if 'class' in params: f_out.write(f"{params['element']}: {params['class']}\n")   # Define new element

    if 'from' in params:
      f_out.write(f"superimpose, element = {params['element']}, ref = {params['from']}, offset = {params['at']}\n")
    else:
      f_out.write(f"superimpose, element = {params['element']}, ref = {common.seqedit_name}_mark, offset = {params['at']}\n")

    return


  # Return

  if dlist[0] == 'return':
    common.f_in[-1].close()
    common.f_in.pop()       # Remove last file handle
    if common.one_file:
      f_out.write(f'\n! Returned to File: {common.f_in[-1].name}\n')
    else:
      common.f_out[-1].close()
      common.f_out.pop()       # Remove last file handle

    return

  # Exit, Quit, Stop

  if dlist[0] == 'exit' or dlist[0] == 'quit' or dlist[0] == 'stop':
    common.f_in.pop()
    if not common.one_file: common.f_out.pop()
    return

  # Title
  # MADX will accept something like "title'abc'" which Bmad will not

  if dlist[0] == 'title':
    if len(dlist) > 1: 
      if dlist[1] == ',':
        f_out.write(command + '\n')
      else:
        f_out.write(f'title, {dlist[1]}\n')
    return

  # endsequence

  if dlist[0] == 'endsequence':
    common.in_seq = False
    seq = common.last_seq
    common.seq_dict[seq.name] = seq
    offset = f'{seq.l} - {add_parens(seq.last_ele_offset, False)}'

    # Replace "[[...]]" marker strings in offsets for elements that have been inserted when ref 
    # element has not yet been defined at the point the element was parsed.

    if not common.superimpose_eles and not is_zero(offset):
      drift_name = f'drift{common.drift_count}'
      seq.drift_list.append(f'{drift_name}: drift, l = {offset}')
      seq.line += drift_name + ', '
      common.drift_count += 1
    
    for ix, drift in enumerate(seq.drift_list):
      while '[[' in drift:
        ix1 = drift.find('[[')
        ix2 = drift.find(']]')
        ref_ele_name = drift[ix1+2:ix2]
        from_ref_ele = seq.seq_ele_dict[ref_ele_name]
        offset = from_ref_ele.at
        if 'l' in from_ref_ele.param:
          if seq.refer == 'entry': offset += f' + {add_parens(bmad_expression(from_ref_ele.param["l"], ""), False)/2}'
          if seq.refer == 'exit': offset += f' - {add_parens(bmad_expression(from_ref_ele.param["l"], ""), False)/2}'
        drift = f'{drift[:ix1]}({offset}){drift[ix2+2:]}'
        seq.drift_list[ix] = drift

    for drift in seq.drift_list:
      f_out.write(drift + '\n')

    #

    if not common.superimpose_eles:
      wrap_write (f'{seq.name}: line = ({seq.line[:-2]})', f_out)

    return

  # Everything below has at least 3 words

  if len(dlist) < 3:
    print ('Unrecognized construct:\n  ' + command.strip())
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

  # Sequence

  if dlist[1] == ':' and dlist[2] == 'sequence':
    common.in_seq = True
    common.last_seq = seq_struct(dlist[0])
    if len(dlist) > 4:
      param_dict = parameter_dictionary(dlist[4:])
      common.last_seq.l = param_dict.get('l', '0')
      common.last_seq.refer = param_dict.get('refer', 'centre')
      common.last_seq.refpos = param_dict.get('refpos', '')
      if 'add_pass' in param_dict: print ('Cannot handle "add_pass" construct in sequence.')
      if 'next_sequ' in param_dict: print ('Cannot handle "next_sequ" construct in sequence.')

    if common.superimpose_eles:
      f_out.write(f'{dlist[0]}_mark: null_ele\n')
      f_out.write(f'{dlist[0]}_drift: drift, l = {common.last_seq.l}\n')
      f_out.write(f'{dlist[0]}: line = ({dlist[0]}_mark, {dlist[0]}_drift)\n')

    return

  # In a sequence construct.

  if common.in_seq:
    seq = common.last_seq
    is_ele_here = True

    # This is an element in the sequence...
    # If "name: name, at = X" construct
    if dlist[0] == dlist[2] and dlist[1] == ':':
      ele = parse_and_write_element(dlist, False, command)
      offset = bmad_expression(ele.at, '')
      ele_name = ele.name

    # "name: type, ..." construct
    elif dlist[1] == ':':
      ele = parse_and_write_element(dlist, True, command)
      common.last_seq.seq_ele_dict[ele.name] = ele
      ele_name = ele.name
      offset = bmad_expression(ele.at, '')

    # If "name, at = X, ..." construct
    elif dlist[0] in common.ele_dict:
      ele = parse_and_write_element([dlist[0], ':']+dlist, False, command)
      offset = bmad_expression(ele.at, '')
      ele_name = ele.name
      # If element has modified parameters. Need to create a new element with a unique name with "__N" suffix.
      if len(ele.param) > 0:
        common.ele_dict[dlist[0]].count += 1
        ele_name = f'{dlist[0]}__{common.ele_dict[dlist[0]].count}'
        ele = parse_and_write_element([ele_name, ':']+dlist, True, command)
      seq.seq_ele_dict[ele_name] = ele    # In case this element is used as a positional reference

    else:   # Subsequence
      ele = ele_struct(dlist[0])
      ele.params = parameter_dictionary(dlist[2:])
      ele.at = ele.params['at']
      seq.seq_ele_dict[dlist[0]] = ele    # In case this element is used as a positional reference      
      is_ele_here = False

    # Finish ele in sequence.

    if is_ele_here:
      if ele.from_ref_ele != '':
        if ele.from_ref_ele in seq.seq_ele_dict:
          from_ref_ele = seq.seq_ele_dict[ele.from_ref_ele]
          offset += f' + {add_parens(bmad_expression(from_ref_ele.at, ""), False)}'
          if 'l' in from_ref_ele.param:
            if seq.refer == 'entry': offset += f' + {add_parens(bmad_expression(from_ref_ele.param["l"], ""), False)/2}'
            if seq.refer == 'exit': offset += f' - {add_parens(bmad_expression(from_ref_ele.param["l"], ""), False)/2}'
        else:
          # Ref element is not yet defined so put in marker string "[[...]]" that will be removed later to
          # be replaced by the actual offset.
          if offset == '':
            offset = f'[[{ele.from_ref_ele}]]'
          else:
            offset += f' + [[{ele.from_ref_ele}]]'

      if common.superimpose_eles:
        f_out.write(f'superimpose, element = {ele_name}, ref = {seq.name}_mark, ' + \
                    f'offset = {offset}, ele_origin = {sequence_refer[seq.refer]}\n')

      else:
        last_offset = f'{offset}'
        this_offset = f'{offset}'
        length = ''

        ele2 = ele
        while True:
          if 'l' in ele2.param: break
          if ele2.madx_inherit not in common.ele_dict: break
          ele2 = common.ele_dict[ele2.madx_inherit]

        if 'l' in ele2.param:
          if ele2.madx_base_type == 'rbend':
            length = f'{ele.name}[l]/sinc({ele2.name}[angle]/2)'
          else:
            length = ele2.param['l']

        if length != '': length = add_parens(bmad_expression(length, ''), False)

        if seq.refer == 'entry':
          if length != '': last_offset += f' + {length}'
        elif seq.refer == 'centre':
          if length != '': this_offset += f' - {length}/2'
          if length != '': last_offset += f' + {length}/2'
        else:
          if length != '': this_offset += f' - {length}'

        if seq.last_ele_offset != '': this_offset += f' - {add_parens(seq.last_ele_offset, False)}'

        if is_zero(this_offset):
          seq.line += f'{ele_name}, '
          seq.last_ele_offset = last_offset
        else:
          drift_name = f'drift{common.drift_count}'
          drift_line = f'{drift_name}: drift, l = {this_offset}'
          seq.drift_list.append(drift_line)
          seq.line += f'{drift_name}, {ele_name}, '
          seq.last_ele_offset = last_offset
          common.drift_count += 1

      return

    # Must be sequence within a sequence.

    ele = parse_and_write_element([dlist[0], ':', 'sequence']+dlist[1:], False, command)
    ele_name = ele.name

    try:
      seq2 = common.seq_dict[ele.name]
    except:
      print (f'CANNOT IDENTIFY THIS AS AN ELEMENT OR SEQUENCE: {dlist[0]}\n  IN LINE IN SEQUENCE: {command}')
      return

    offset = bmad_expression(ele.at, '')

    if ele.from_ref_ele != '':
      from_ref_ele = seq.ele_dict[ele.from_ref_ele]
      offset = f'{offset} - {add_parens(bmad_expression(from_ref_ele.at, ""), False)}'

    last_offset = offset
    length = add_parens(bmad_expression(seq2.l, ''), False)
    this_offset = f'{offset}'

    if seq2.refpos != '':
      refpos_ele = seq2.seq_ele_dict[seq2.refpos]
      offset += f' - {add_parens(refpos_ele.at, False)}'
      last_offset += f' + {refpos_ele.at} - {add_parens(seq2.l, False)}'
      print (f'A: {last_offset}')
    elif seq.refer == 'entry':
      if length != '': last_offset += f' + {length}'
      print (f'B: {last_offset}')
    elif seq.refer == 'centre':
      offset += f' - {add_parens(length, False)}/2'
      if length != '': this_offset += f' - {length}/2'
      if length != '': last_offset += f' + {length}/2'
      print (f'C: {last_offset}')
    else:
      offset += f' - {add_parens(length, False)}'
      if length != '': this_offset += f' - {length}'

    if common.superimpose_eles:
      common.super_list.append(f'superimpose, element = {ele.name}_mark, ref = {seq.name}_mark, offset = {offset}\n')
      f_out.write (f'!!** superimpose, element = {ele.name}_mark, ref = {seq.name}_mark, offset = {offset}\n')

    elif not is_zero(this_offset):
      drift_name = f'drift{common.drift_count}'
      drift_line = f'{drift_name}: drift, l = {this_offset}'
      common.drift_count += 1

      if seq.last_ele_offset != '': drift_line += f' - {add_parens(seq.last_ele_offset, False)}'
      seq.drift_list.append(drift_line)
      print (f'3: {seq.drift_list[-1]}')
      seq.line += f'{drift_name}, {ele_name}, '
      seq.last_ele_offset = last_offset

    return


  #-----------------------------------------------
  # Line.
  # Nonstandard "a: line = -b" must be converted to "a: line = (-b)"

  if ix_colon > 0 and dlist[ix_colon+1] == 'line':
    words = command.split('=')
    if words[1].strip()[0] != '(': command = words[0] + '= (' + words[1].strip() + ')'
    wrap_write(command, f_out)
    return

  # Var definition.
  # If a variable value is an expression that involves an element parameter, write
  # the line instead of adding the line to the var list. This is done to avoid moving
  # the def to before the point where the element is defined.

  if dlist[1] == '=' and not '->' in dlist[0]:
    if dlist[0] in common.var_name_list:
      print (f'Duplicate variable name: {dlist[0]}\n' + 
             f'  You may have to edit the Bmad lattice file by hand to resolve this.')
    common.var_name_list.append(dlist[0])
    name = dlist[0]
    value = bmad_expression(command.split('=')[1].strip(), '')
    if '[' in value or not common.prepend_vars:    # Involves an element parameter
      f_out.write(f'{name} = {value}\n')
    else:
      common.var_def_list.append([name, value])
    return

  # "qf, k1 = ..." parameter set

  if len(dlist) > 4 and dlist[0] in common.ele_dict and dlist[1] == ',' and dlist[3] == '=':
    f_out.write(dlist[0] + '[' + bmad_param(dlist[2]) + '] = ' + bmad_expression(''.join(dlist[4:]), dlist[2]) + '\n')
    return


  # Ele "qf->k1 = ..." parameter set

  if dlist[1] == '=' and '->' in dlist[0]:
    [ele_name, dummy, param] = dlist[0].partition('->')
    value = bmad_expression(command.split('=')[1].strip(), param)
    name = f'{ele_name}[{bmad_param(param, ele_name)}]'
    f_out.write(f'{name} = {value}\n')
    return

  # Call

  if dlist[0] == 'call':

    file = command.split('=')[1].strip()
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

  # Use

  if dlist[0] == 'use':
    if len(dlist) == 3:
      common.use = dlist[2]
    else:
      params = parameter_dictionary(dlist[2:])
      if 'sequence' in params: common.use = params.get('sequence')
      if 'period' in params:  common.use = params.get('period')

    f_out.write('use, ' + common.use + '\n')
    return

  # Beam

  if dlist[0] == 'beam' or dlist[2] == 'beam':
    if dlist[0] == 'beam': param = parameter_dictionary(dlist[2:])
    if dlist[2] == 'beam': param = parameter_dictionary(dlist[4:])
    if 'particle' in param:  f_out.write('parameter[particle] = ' + bmad_expression(param['particle'], '') + '\n')
    if 'energy'   in param:  f_out.write('parameter[E_tot] = ' + bmad_expression(param['energy'], 'energy') + '\n')
    if 'pc'       in param:  f_out.write('parameter[p0c] = ' + bmad_expression(param['pc'], 'pc') + '\n')
    if 'gamma'    in param:  f_out.write('parameter[E_tot] = mass_of(parameter[particle]) * ' + add_parens(bmad_expression(param['gamma'], ''), False) + '\n')
    if 'npart'    in param:  f_out.write('parameter[n_part] = ' + bmad_expression(param['npart'], '') + '\n')
    return

  # twiss

  if dlist[0] == 'twiss' or dlist[2] == 'beta0':
    if dlist[0] == 'twiss':
      param = parameter_dictionary(dlist[2:])
    else:
      param = parameter_dictionary(dlist[4:])
    if 'betx'   in param: f_out.write(f'beginning[beta_a] = {bmad_expression(param["betx"], "")}\n')
    if 'bety'   in param: f_out.write(f'beginning[beta_b] = {bmad_expression(param["bety"], "")}\n')
    if 'alfx'   in param: f_out.write(f'beginning[alpha_a] = {bmad_expression(param["alfx"], "")}\n')
    if 'alfy'   in param: f_out.write(f'beginning[alpha_a] = {bmad_expression(param["alfy"], "")}\n')
    if 'mux'    in param: f_out.write(f'beginning[phi_a] = twopi * {add_parens(bmad_expression(param["mux"], ""), False)}\n')
    if 'muy'    in param: f_out.write(f'beginning[phi_b] = twopi * {add_parens(bmad_expression(param["muy"], ""), False)}\n')
    if 'dx'     in param: f_out.write(f'beginning[eta_x] = {bmad_expression(param["dx"], "")}\n')
    if 'dy'     in param: f_out.write(f'beginning[eta_y] = {bmad_expression(param["dy"], "")}\n')
    if 'dpx'    in param: f_out.write(f'beginning[etap_x] = {bmad_expression(param["dpx"], "")}\n')
    if 'dpy'    in param: f_out.write(f'beginning[etap_y] = {bmad_expression(param["dpy"], "")}\n')
    if 'x'      in param: f_out.write(f'particle_start[x] = {bmad_expression(param["x"], "")}\n')
    if 'y'      in param: f_out.write(f'particle_start[y] = {bmad_expression(param["y"], "")}\n')
    if 'px'     in param: f_out.write(f'particle_start[px] = {bmad_expression(param["px"], "")}\n')
    if 'py'     in param: f_out.write(f'particle_start[py] = {bmad_expression(param["py"], "")}\n')
    return

  # Element def

  if dlist[1] == ':':
    parse_and_write_element(dlist, True, command)
    return

  # Unknown

  print (f"Unknown construct:\n    " + command.strip())

#------------------------------------------------------------------
#------------------------------------------------------------------
# Get next madx command.
# Read in MADX file line-by-line.  Assemble lines into commands, which are delimited by a ; (colon).

def read_madx_command ():
  global common

  quote_delim = ''  # Quote mark delimiting a string. Blank means not parsing a string yet.
  in_extended_comment = False
  command = ''
  dlist = []
  curly_brace_count = 0   # Count "{", "}" pairs

  # Loop until a command has been found.
  # Note: "macro" and "if" statements are strange since they are permitted to 
  # not end with a ';' but with a matching '}'

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
        if not common.one_file:
          common.f_out[-1].close()
          common.f_out.pop()       # Remove last file handle
        if len(common.f_in) == 0: return ['', dlist]  # If root file was closed

    else:
      f_in = common.f_in[-1]
      f_out = common.f_out[-1]
      line = common.command
      common.command = ''

    # Parse line

    line = line.strip()
    if line.strip() == '':
      f_out.write('\n')
      continue

    if len(line) > 1 and line[0:2] == '#!':   # "#!madx" line
      f_out.write('! ' + line + '\n')
      line = ''
      continue

    if in_extended_comment:
      if '*/' in line:
        ix = line.find('*/')
        f_out.write('! ' + line[:ix] + '\n')
        line = line[ix+2:].strip()
        in_extended_comment = False
      else:
        f_out.write ('! ' + line + '\n')
        line = ''
      continue

    while line != '':
      for ix in range(len(line)):
        #print (f'Ix: {ix} {line[ix]} -{quote_delim}- ' + line)
        #print (f'C: {command}')
        #print (f'D: {dlist}')

        if line[ix] == '{': curly_brace_count += 1
        if line[ix] == '}': 
          curly_brace_count -= 1
          if curly_brace_count == 0 and len(dlist) > 0 and (dlist[0] in ['if', 'elseif', 'else' 'while'] or 'macro' in dlist):
            command += line[:ix]
            if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
            common.command = line[ix+1:]
            return [command, dlist]

        if (line[ix] == '"' or line[ix] == "'"):
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

        if quote_delim != '': continue    # Cycle if in quote string

        if line[ix] == '!':
          if len(line) > ix+10 and line[ix:ix+10] == '!!verbatim':
            f_out.write(line[ix+10:].strip() + '\n')
          else:
            f_out.write(line[ix:] + '\n')
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          line = ''
          break

        # "if" or "macro" commands can have internal ";" characters that need to be ignored.
        elif line[ix] == ';' and not ((len(dlist) > 0 and dlist[0] in ['if', 'elseif', 'else', 'while']) or 'macro' in dlist):
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          common.command = line[ix+1:]
          return [command, dlist]

        # Need to split "if(" or "while(" constructs at "(". 
        # This only is necessary at the start of the command string.
        elif line[ix] in '{}:,=' or (len(dlist) == 0 and line[ix] == '('):
          command += line[:ix+1]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          dlist.append(line[ix])
          line = line[ix+1:]
          break

        elif ix == len(line) - 1:
          command += line
          dlist.append(line.strip())
          line = ''
          break

        if ix > len(line) - 2: continue

        if line[ix:ix+2] == '/*':
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          if '*/' in line[ix:]:
            ix2 = line.find('*/', ix)
            f_out.write('!' + line[ix+2:ix2] + '\n')
            line = line[ix+3:]
          else:
            f_out.write('!' + line[ix+2:] + '\n')
            line = ''
            in_extended_comment = True
          break

        elif line[ix:ix+2] == '//':
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          f_out.write('!' + line[ix+2:] + '\n')
          line = ''
          break

#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
# Main program.

start_time = time.time()

# Read the parameter file specifying the MADX lattice file, etc.

argp = argparse.ArgumentParser()
argp.add_argument('madx_file', help = 'Name of input MADX lattice file')
argp.add_argument('-d', '--debug', help = 'Print debug info (not of general interest).', action = 'store_true')
argp.add_argument('-f', '--many_files', help = 'Create a Bmad file for each MADX input file.', action = 'store_true')
argp.add_argument('-s', '--superimpose', help = 'Superimpose elements in a sequence.', action = 'store_true')
argp.add_argument('-v', '--no_prepend_vars', help = 'Do not move variables to the beginning of the Bmad file.', action = 'store_true')
arg = argp.parse_args()

common = common_struct()
common.debug = arg.debug
common.superimpose_eles = arg.superimpose
common.prepend_vars = not arg.no_prepend_vars
common.one_file = not arg.many_files

madx_lattice_file = arg.madx_file
bmad_lattice_file = bmad_file_name(madx_lattice_file)

print ('Input lattice file is:  ' + madx_lattice_file)
print ('Output lattice file is: ' + bmad_lattice_file)

# Open files for reading and writing

common.f_in.append(open(madx_lattice_file, 'r'))  # Store file handle
common.f_out.append(open(bmad_lattice_file, 'w'))

f_out = common.f_out[-1]

#------------------------------------------------------------------
# parse, convert and output madx commands

common.command = ''  # init

while True:
  [command, dlist] = read_madx_command()
  if len(common.f_in) == 0: break
  parse_command(command, dlist)
  if len(common.f_in) == 0: break   # Hit Quit/Exit/Stop statement.

f_out.close()

#------------------------------------------------------------------
# Prepend variables and superposition statements as needed

f_out = open(bmad_lattice_file, 'r')
lines = f_out.readlines()
f_out.close()

f_out = open(bmad_lattice_file, 'w')
f_out.write (f'!+\n! Translated from MADX to Bmad by madx_to_bmad.py\n! File: {madx_lattice_file}\n!-\n\n')

if common.prepend_vars:
  order_var_def_list()
  for vdef in common.var_def_list:
    wrap_write(f'{vdef[0]} = {vdef[1]}\n', f_out)
  f_out.write('\n')

if len(common.super_list) > 0:
  for line in common.super_list:
    f_out.write(line)
  f_out.write('\n')

for line in lines:
  f_out.write(line)

f_out.close()
