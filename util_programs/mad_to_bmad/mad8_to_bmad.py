#!/usr/bin/env python

#+
# Script to convert from MAD8 lattice format to Bmad lattice format.
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
    self.mad8_inherit_type = ''    # Is another element name or quadrupole, etc.
    self.mad8_base_type    = ''    # Is quadrupole, etc.
    self.bmad_inherit_type = ''    # Is another element name or quadrupole, etc.
    self.bmad_base_type    = ''    # Is quadrupole, etc.
    self.at = '0'                  # Used if element is in a sequence
    self.from_ref_ele = ''         # Used if element is in a sequence
    self.param = OrderedDict()
    self.count = 0

class seq_struct:
  def __init__(self, name = ''):
    self.name = name
    self.refer = 'centre'
    self.ele_dict = OrderedDict()

class common_struct:
  def __init__(self):
    self.debug = False
    self.prepend_vars = True
    self.one_file = True
    self.in_seq = False
    self.seqedit_name = ''           # Name of sequence in seqedit construct.
    self.last_seq = seq_struct()     # Current sequence being parsed.
    self.seq_dict = OrderedDict()    # List of all sequences.
    self.super_list = []             # List of superimpose statements to be prepended to the bmad file.
    self.ele_dict = {}               # Dict of elements
    self.var_def_list = []           # List of "A = B" sets after translation to Bmad. Does not Include "A,P = B" parameter sets.
    self.var_name_list = []          # List of mad8 variable names.
    self.f_in = []         # MAD8 input files
    self.f_out = []        # Bmad output files
    self.use = ''
    self.command = ''

#------------------------------------------------------------------
#------------------------------------------------------------------

ele_param_factor = {
  'deltae':   ' * 1e-6',
  'volt':     ' * 1e-6',
  'freq':     ' * 1e-6',
  'energy':   ' * 1e-9',
  'ex':       ' * 1e-6',
  'ey':       ' * 1e-6',
  'pc':       ' * 1e-9',
  'lag':      ' + 0.5',
}

ele_inv_param_factor = {
  'deltae':   ' * 1e6',
  'volt':     ' * 1e6',
  'freq':     ' * 1e6',
  'energy':   ' * 1e9',
  'ex':       ' * 1e6',
  'ey':       ' * 1e6',
  'pc':       ' * 1e9',
  'lag':      ' - 0.5',
}

const_trans = {
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
    'matrix':       'taylor', 
    'roll':         'patch',
    'srot':         'patch',
    'yrot':         'patch',
    'zrot':         'patch',
    'monitor':      'monitor',
    'hmonitor':     'monitor',
    'vmonitor':     'monitor',
    'blmonitor':    'monitor',
    'profile':      'monitor',
    'wire':         'monitor',
    'slmonitor':    'monitor',
    'imonitor':     'monitor',
    'marker':       'marker',
    'drift':        'drift',
    'sbend':        'sbend',
    'rbend':        'rbend',
    'quadrupole':   'quadrupole',
    'sextupole':    'sextupole',
    'octupole':     'octupole',
    'multipole':    'multipole',
    'dimultipole':  'multipole',
    'solenoid':     'solenoid',
    'hkicker':      'hkicker',
    'vkicker':      'vkicker',
    'kicker':       'kicker',
    'rfcavity':     'rfcavity',
    'elseparator':  'elseparator',
    'instrument':   'instrument',
    'ecollimator':  'ecollimator',
    'rcollimator':  'rcollimator',
    'beambeam':     'beambeam',
    'lcavity':      'lcavity',
    'lump':         '???',
    'arbitelm':     '???',
    'mtwiss':       '???'
}

ignore_mad8_param = ['lrad', 'slot_id', 'aper_tol', 'apertype', 'assembly_id', 'removed',
                     'mech_sep', 'betrf', 'tfill', 'shunt', 'pg', 'eloss', 'lfile', 'tfile', 'e0']

bmad_param_name = {
    'deltae':   'voltage',
    'volt':     'voltage',
    'freq':     'rf_frequency',
    'lag':      'phi0',
    'e ':       'e_field',
    'xsize':    'x_limit',
    'ysize':    'y_limit',
    'lrad':     'l',
    'swave':    'cavity_type',
}

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
      if new_def_list[ix2][0] in vdef[1]:
        new_def_list.pop(ix)
        new_def_list.insert(ix2, vdef)
        moved = True
        break

    if not moved: ix += 1

  common.var_def_list = new_def_list

#------------------------------------------------------------------
#------------------------------------------------------------------
# Convert from mad8 parameter name to bmad parameter name.

def bmad_param(param, ele_name):
  global common, bmad_param_name

  # For the SLAC version there are Rij and Tijk matrix elements

  if len(param) == 3 and param[0] == 'r' and param[1] in '123456' and param[2] in '123456':
    return f'tt{param[1:]}'

  if len(param) == 4 and param[0] == 't' and param[1] in '123456' and param[2] in '123456' and param[3] in '123456':
    return f'tt{param[1:]}'

  elif param in bmad_param_name:
    return bmad_param_name[param]

  if ele_name in common.ele_dict:
    mad8_type = common.ele_dict[ele_name].mad8_base_type
  else:
    mad8_type = 'xxxx'

  if mad8_type == 'dimultipole':
    d = '0123456789'
    if len(param) == 2 and param[0] == 'k' and param[1].isdigit(): return param_name + 'l'
    if len(param) == 3 and param[0] == 'k' and param[1].isdigit() and param[2].isdigit(): return param_name + 'l'

  if param == 'angle':
    if mad8_type == 'srot': return 'tilt'
    if mad8_type == 'roll': return 'tilt'
    if mad8_type == 'yrot': return 'x_pitch'
    if mad8_type == 'zrot': return 'x_pitch'

  if param == 'tilt':
    if mad8_type == 'sbend' or mad8_type == 'rbend': return 'ref_tilt'

  return param

#------------------------------------------------------------------
#------------------------------------------------------------------
# Return dictionary of "A = value" parameter definitions.
# Also see parameter_dictionary routine.
# The difference is that this routine will not split on commas.

def action_parameter_dictionary(word_lst):
  pdict = OrderedDict()
  while True:
    if len(word_lst) == 0:
      return pdict

    elif len(word_lst) == 1:
      pdict[word_lst[0]] = ''
      return pdict

    elif word_lst[1] == '=':
      try:
        ix = word_lst.index(',')
        pdict[word_lst[0]] = ' '.join(word_lst[2:ix])
        word_lst = word_lst[ix+1:]
      except:
        pdict[word_lst[0]] = ' '.join(word_lst[2:])
        return pdict

    elif word_lst[1] == ',':
      pdict[word_lst[0]] = ''
      word_lst = word_lst[2:]

    else:
      print (f'CANNOT PARSE: {word_lst}')
      return pdict

    
#------------------------------------------------------------------
#------------------------------------------------------------------
# Return dictionary of "A = value" parameter definitions.
# Also see action_parameter_dictionary routine.
# The difference is that this routine will not split on commas.

def parameter_dictionary(word_lst):

  orig_word_lst = word_lst

  # replace "0." or "0.0" with "0"
  word_lst = ['0' if x == '0.0' or x == '0.' else x for x in word_lst]

  # Replace "rm" and "tm" matrix constructions with Bmad equivalents

  for ix in range(len(word_lst)):
    if word_lst[ix].startswith('kick('):  # Something like: ['kick(3)']
      iz = min(ix+1, len(word_lst))
      word_lst = word_lst[:ix] + [f'tt{word_lst[ix][-2]}'] + word_lst[iz:]

    elif word_lst[ix].startswith('rm('):  # Something like: ['rm(3', ',', '4)']
      iz = min(ix+3, len(word_lst))
      word_lst = word_lst[:ix] + [f'tt{word_lst[ix][-1]}{word_lst[ix+2][0]}'] + word_lst[iz:]

    elif word_lst[ix].startswith('tm('):  # Something like: ['tm(3', ',', '4', ',', '2)']
      iz = min(ix+5, len(word_lst))
      word_lst = word_lst[:ix] + [f'tt{word_lst[ix][-1]}{word_lst[ix+2][0]}{word_lst[ix+4][0]}'] + word_lst[iz:]

    if ix > len(word_lst) - 2: break

  # Fill dict

  pdict = OrderedDict()
  while True:
    if len(word_lst) == 0: return pdict

    if 'tilt' in word_lst:
      word_lst = mad_defaulting_parameter_to_dict('tilt', word_lst, orig_word_lst, pdict)
      continue

    if 'swave' in word_lst:
      word_lst = mad_defaulting_parameter_to_dict('swave', word_lst, orig_word_lst, pdict)
      continue
      
    if word_lst[1] != '=':
      print ('PROBLEM PARSING PARAMETER LIST: ' + ''.join(orig_word_lst))
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
# A "defaulting parameter" is a parameter that may have a default value. For example: param = 'tilt'.

def mad_defaulting_parameter_to_dict(param, word_lst, orig_word_lst, pdict):
  ix = word_lst.index(param)
  if len(word_lst) > ix + 1 and word_lst[ix+1] == ',':
    pdict[param] = ''
    word_lst.pop(ix+1)
    word_lst.pop(ix)
  elif len(word_lst) == ix + 1 and word_lst[ix-1] == ',':
    pdict[param] = ''
    word_lst.pop(ix)
    word_lst.pop(ix-1)
  elif len(word_lst) > ix + 1 and word_lst[ix+1] == '=':
    if ',' in word_lst[ix+1:]:
      ixe = word_lst.index(',', ix+1)
      pdict[param] = ' '.join(word_lst[ix+2:ixe])
      word_lst = word_lst[:ix] + word_lst[ixe+1:]
    else:
      pdict[param] = ' '.join(word_lst[ix+2:])
      word_lst = word_lst[:ix]
      if ix > 0 and word_lst[-1] == ',': word_lst.pop()
  else:
    print ('PROBLEM PARSING "' + param + '" IN PARAMETER LIST: ' + ''.join(orig_word_lst))

  return word_lst

#------------------------------------------------------------------
#------------------------------------------------------------------
# Convert expression from MAD8 format to Bmad format
# To convert <expression> a construct that look like "<target_param> = <expression>".

def bmad_expression(line, target_param):
  global common, const_trans, ele_param_factor

  lst = re.split(r'(,|-|\+|\(|\)|\>|\*|/|\^|\[|\])', line)
  out = ''

  while len(lst) != 0:
    if len(lst) > 3 and lst[1] == '[' and lst[3] == ']':
      if lst[2] in ele_param_factor:
        if (len(lst) >= 5 and lst[4] == '^') or (len(out.strip()) > 0 and out.strip()[-1] == '/'):
          out += '(' + lst[0] + '[' + bmad_param(lst[2].strip(), lst[0]) + ']' + ele_param_factor[lst[2]]
        else:
          out += lst[0] + '[' + bmad_param(lst[2].strip(), lst[0]) + ']' + ele_param_factor[lst[2]]
      else:
        out += lst[0] + '[' + bmad_param(lst[2].strip(), lst[0]) + ']'

      lst = lst[4:]

    elif lst[0] in const_trans:
      out += const_trans[lst.pop(0)]

    else:
      out += lst.pop(0)

  # End while

  if target_param in ele_inv_param_factor: out = add_parens(out) + ele_inv_param_factor[target_param]
  return out

#-------------------------------------------------------------------
#------------------------------------------------------------------
# Construct the bmad lattice file name

def bmad_file_name(mad8_file):

  if mad8_file.lower()[-5:] == '.mad8':
    return mad8_file[:-5] + '.bmad'

  elif mad8_file.lower()[-5:] == '.xsif':
    return mad8_file[:-5] + '.bmad'

  elif mad8_file.lower()[-4:] == '.mad':
    return mad8_file[:-4] + '.bmad'

  elif mad8_file.lower()[-4:] == '.seq':
    return mad8_file[:-4] + '.bmad'

  else:
    return mad8_file + '.bmad'

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

def negate(str):
  str = add_parens(str)
  if str[0] == '-':
    return str[1:]
  elif str[0] == '+':
    return '-' + str[1:]
  else:
    return '-' + str

#------------------------------------------------------------------
#------------------------------------------------------------------
# Parse a lattice element

def parse_element(dlist, write_to_file):
  global common, ele_type_translate, ignore_mad8_param

  if dlist[2] in common.ele_dict:
    ele = ele_struct(dlist[0])
    ele.mad8_inherit_type = dlist[2]
    ele.mad8_base_type = common.ele_dict[dlist[2]].mad8_base_type
    ele.bmad_inherit_type = dlist[2]
    ele.bmad_base_type = common.ele_dict[dlist[2]].bmad_base_type

  else:
    found = False
    for mad8_type in ele_type_translate:
      if mad8_type.startswith(dlist[2]): 
        ele = ele_struct(dlist[0])
        ele.mad8_inherit_type = mad8_type
        ele.mad8_base_type = mad8_type
        ele.bmad_inherit_type = ele_type_translate[mad8_type]
        ele.bmad_base_type = ele.bmad_inherit_type
        found = True
        break

    if not found:
      print (dlist[2].upper() + ' TYPE ELEMENT IS UNKNOWN!')
      return

  if ele.mad8_base_type == '???':
    print (dlist[2].upper() + ' TYPE ELEMENT CANNOT BE TRANSLATED TO BMAD.')
    return

  params = parameter_dictionary(dlist[4:])

  if ele.mad8_base_type == 'yrot' or ele.mad8_base_type == 'zrot':
    if 'angle' in params: params['x_pitch'] = negate(params.pop('angle'))

  elif ele.mad8_base_type == 'srot' or ele.mad8_base_type == 'roll':
    if 'angle' in params: params['tilt'] = params.pop('angle')

  elif ele.mad8_base_type == 'rbend' or ele.mad8_base_type == 'sbend':
    if 'tilt' in params and params['tilt'] == '': params['tilt'] = 'pi/2'
    if 'tilt' in params: params['ref_tilt'] = params.pop('tilt')

  elif ele.mad8_base_type == 'quadrupole':
    if 'tilt' in params and params['tilt'] == '': params['tilt'] = 'pi/4'

  elif ele.mad8_base_type == 'sextupole':
    if 'tilt' in params and params['tilt'] == '': params['tilt'] = 'pi/6'

  elif ele.mad8_base_type == 'octupole':
    if 'tilt' in params and params['tilt'] == '': params['tilt'] = 'pi/8'

  elif ele.mad8_base_type == 'lcavity':
    if 'swave' not in params: params['swave'] = '.F.'

  #

  if 'swave' in params:
    if params['swave'] == '': params['swave'] = '.T.'
    if   params['swave'].upper() in ['.YES.', '.TRUE.', '.T.', '.ON.']:  params['swave'] = 'standing_wave' 
    elif params['swave'].upper() in ['.NO.', '.FALSE.', '.F.', '.OFF.']: params['swave'] = 'traveling_wave' 

  #

  ele.param = params
  common.ele_dict[dlist[0]] = ele

  if write_to_file:
    line = ele.name + ': ' + ele.bmad_inherit_type
    for param in ele.param:
      if param in ignore_mad8_param: continue
      line += ', ' + bmad_param(param, ele.name) + ' = ' + bmad_expression(params[param], param)
    f_out = common.f_out[-1]
    wrap_write(line, f_out)

  return ele

#------------------------------------------------------------------
#------------------------------------------------------------------

def parse_command(command, dlist):
  global common, sequence_refer

  f_out = common.f_out[-1]
  if common.debug: print (str(dlist))

  if len(dlist) == 0: 
    f_out.write('\n')
    return

  # Ignore this

  if dlist[0] in ['show', 'efcomp', 'print', 'select', 'optics', 'option', 'survey',
                  'emit', 'twiss', 'help', 'eoption', 'system', 'ealign', 'savebeta', 'assign']:
    return

  if dlist[0] in ['set']:
    print ('Note! Ignoring command: ' + command)
    return

  # Flag this

  if dlist[0] in ['cycle', 'move', 'remove', 'replace', 'extract']:
    print (f'WARNING! CANNOT TRANSLATE THE COMMAND: {dlist[0].upper()}')
    return

  # Return

  if dlist[0] == 'return':
    common.f_in[-1].close()
    common.f_in.pop()                  # Remove last file handle
    if len(common.f_in) == 0: return   # In case the return statement is in the master file.
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
  # MAD8 will accept something like "title'abc'" which Bmad will not

  if dlist[0] == 'title':
    if len(dlist) > 1: 
      if dlist[1] == ',':
        f_out.write(command + '\n')
      else:
        f_out.write(f'title, {dlist[1]}')
    return

  # Get rid of "real", "int", "const" "const real", etc. prefix

  if dlist[0].startswith('real ') or dlist[0].startswith('int ') or dlist[0].startswith('const '): dlist[0] = dlist[0].split(' ', 1)[1].strip()
  if dlist[0].startswith('real ') or dlist[0].startswith('int ') or dlist[0].startswith('const '): dlist[0] = dlist[0].split(' ', 1)[1].strip()
  if len(dlist) > 3 and dlist[1] == ':' and dlist[2] == 'constant': dlist.pop(2)

  # Get rid of constant in "name: constant = ..." expression

  if len(dlist) > 4 and dlist[1] == ':' and dlist[2] == 'constant' and dlist[3] == '=':
    dlist = [dlist[0], '='] + dlist[4:] 

  # Transform: "a := 3" -> "a = 3"

  for ix in range(len(dlist)):
    if ix > len(dlist)-2: break
    if dlist[ix] == ':' and dlist[ix+1] == '=': dlist.pop(ix)

  # The MAD8 parser takes a somewhat cavilier attitude towards commas and sometimes they can be omitted.
  # Examples: "call file" instead of  "call, file" and "q: quadrupole l = 7" instead of "q: quadrupole, l = 7
  # So put the comma back in to make things uniform for easier parsing.
  # But do not do this with arithmetical expressions and quoted strings

  i = 0
  while i < len(dlist):
    if ' ' in dlist[i] and not any(char in dlist[i] for char in '"\'+-*/^'):
      split = dlist[i].split()
      dlist = dlist[:i] + [split[0], ',', split[1]] + dlist[i+1:]
    i += 1

  #

  if dlist[0] == 'endsequence':
    common.in_seq = False
    seq = common.last_seq
    common.seq_dict[seq.name] = seq
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

    f_out.write(f'{dlist[0]}_mark: null_ele\n')
    f_out.write(f'{dlist[0]}_drift: drift, l = {common.last_seq.l}\n')
    f_out.write(f'{dlist[0]}: line = ({dlist[0]}_mark, {dlist[0]}_drift)\n')
    return

  # In a sequence construct.

  if common.in_seq:
    seq = common.last_seq

    if dlist[1] == ':':   # "name: type"  construct
      ele = parse_element(dlist, True)
      common.last_seq.ele_dict[ele.name] = ele
      f_out.write(f'superimpose, element = {ele.name}, ref = {seq.name}_mark, ' + \
                  f'offset = {ele.at}, ele_origin = {sequence_refer[seq.refer]}\n')

    elif dlist[0] in common.ele_dict:
      ele = parse_element([dlist[0], ':']+dlist, False)
      if len(ele.param) == 0:
        name = ele.name
      # element has modified parameters. Need to create a new element with a unique name with "__N" suffix.
      else:
        common.ele_dict[dlist[0]].count = common.ele_dict[dlist[0]].count + 1
        name = f'{dlist[0]}__{common.ele_dict[dlist[0]].count}'
        ele = parse_element([name, ':']+dlist, True)

      offset = ele.at
      if ele.from_ele != '':
        from_ele = seq.ele_dict[ele.from_ele]
        offset = f'{offset} - {add_parens(from_ele.at)}'

      f_out.write(f'superimpose, element = {name}, ref = {seq.name}_mark, ' + \
                  f'offset = {offset}, ele_origin = {sequence_refer[seq.refer]}\n')

    # Must be sequence name. So superimpose the corresponding marker.
    else:
      ele = parse_element([dlist[0], ':']+dlist, False)
      seq2 = common.seq_dict[ele.name]
      offset = ele.at

      if ele.from_ele != '':
        from_ele = seq.ele_dict[ele.from_ele]
        offset = f'{offset} - {add_parens(from_ele.at)}'

      if seq2.refpos != '':
        refpos_ele = seq.ele_dict[seq2.refpos]
        offset = f'{offset} - {add_parens(refpos_ele.at)}'
      elif seq2.refer == 'centre':
        offset = f'{offset} - {add_parens(seq2.l)} / 2'
      elif seq2.refer == 'exit':
        offset = f'{offset} - {add_parens(seq2.l)}'

      common.super_list.append(f'superimpose, element = {ele.name}_mark, ref = {seq.name}_mark, offset = {offset}\n')
      f_out.write (f'!!** superimpose, element = {ele.name}_mark, ref = {seq.name}_mark, offset = {offset}\n')

    return

  # Line

  if ix_colon > 0 and dlist[ix_colon+1] == 'line':
    wrap_write(command, f_out)
    return

  # Var definition.
  # If a variable value is an expression that involves an element parameter, write
  # the line instead of adding the line to the var list. This is done to avoid moving
  # the def to before the point where the element is defined.

  if dlist[1] == '=':
    if dlist[0] in common.var_name_list:
      print (f'Duplicate variable name: {dlist[0]}\n' + 
             f'  You may have to edit the Bmad lattice file by hand to resolve this problem.')

    common.var_name_list.append(dlist[0])
    name = dlist[0]
    value = bmad_expression(''.join(dlist[2:]), dlist[0])
    if '[' in value or not common.prepend_vars:    # Involves an element parameter
      f_out.write(f'{name} = {value}\n')
    else:
      common.var_def_list.append([name, value])

    return

  # Ele parameter set

  if dlist[0] in common.ele_dict and dlist[1] == ',' and dlist[3] == '=':
    value = bmad_expression(command.split('=')[1].strip(), dlist[0])
    if dlist[0] in common.ele_dict:
      ele_name = common.ele_dict[dlist[0]].name
      name = f'{dlist[0]}[{bmad_param(dlist[2], ele_name)}]'
    else:  # In a complete valid lattice, parameter seets always happen after the element has been defined
      name = f'{dlist[0]}[{bmad_param(dlist[2], "???")}]'
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
      params = action_parameter_dictionary(dlist[2:])
      if 'period' in params:  common.use = params.get('period')

    f_out.write('use, ' + common.use + '\n')
    return

  # Beam

  if dlist[0] == 'beam' or (dlist[1] == ':' and dlist[2] == 'beam'):
    if dlist[0] == 'beam':
      params = action_parameter_dictionary(dlist[2:])
    else:
      params = action_parameter_dictionary(dlist[4:])

    if 'particle' in params:  f_out.write('parameter[particle] = ' + bmad_expression(params['particle'], '') + '\n')
    if 'energy' in params:    f_out.write('parameter[E_tot] = ' + bmad_expression(params['energy'], 'energy') + '\n')
    if 'pc' in params:        f_out.write('parameter[p0c] = ' + bmad_expression(params['pc'], 'pc') + '\n')
    if 'gamma' in params:     f_out.write('parameter[E_tot] = mass_of(parameter[particle]) * ' + add_parens(bmad_expression(params['gamma'], '')) + '\n')
    if 'npart' in params:     f_out.write('parameter[n_part] = ' + bmad_expression(params['npart'], '') + '\n')
    return

  # Beta0

  if dlist[0] == 'beta0' or (dlist[1] == ':' and dlist[2] == 'beta0'):
    if dlist[0] == 'beta0':
      params = action_parameter_dictionary(dlist[2:])
    else:
      params = action_parameter_dictionary(dlist[4:])

    if 'betx' in params:      f_out.write(f'beginning[beta_a] = {bmad_expression(params["betx"], "")}\n')
    if 'bety' in params:      f_out.write(f'beginning[beta_b] = {bmad_expression(params["bety"], "")}\n')
    if 'alfx' in params:      f_out.write(f'beginning[alpha_a] = {bmad_expression(params["alfx"], "")}\n')
    if 'alfy' in params:      f_out.write(f'beginning[alpha_b] = {bmad_expression(params["alfy"], "")}\n')
    if 'mux' in params:       f_out.write(f'beginning[phi_a] = twopi * {add_parens(bmad_expression(params["mux"], ""))}\n')
    if 'muy' in params:       f_out.write(f'beginning[phi_b] = twopi * {add_parens(bmad_expression(params["muy"], ""))}\n')
    if 'dx' in params:        f_out.write(f'beginning[eta_x] = {bmad_expression(params["dx"], "")}\n')
    if 'dy' in params:        f_out.write(f'beginning[eta_y] = {bmad_expression(params["dy"], "")}\n')
    if 'dpx' in params:       f_out.write(f'beginning[etap_x] = {bmad_expression(params["dpx"], "")}\n')
    if 'dpy' in params:       f_out.write(f'beginning[etap_y] = {bmad_expression(params["dpy"], "")}\n')
    if 'x' in params:         f_out.write(f'particle_start[x] = {bmad_expression(params["x"], "")}\n')
    if 'y' in params:         f_out.write(f'particle_start[y] = {bmad_expression(params["y"], "")}\n')
    if 'px' in params:        f_out.write(f'particle_start[px] = {bmad_expression(params["px"], "")}\n')
    if 'py' in params:        f_out.write(f'particle_start[py] = {bmad_expression(params["py"], "")}\n')
    if 'energy' in params:    f_out.write(f'parameter[E_tot] = {bmad_expression(params["energy"], "energy")}\n')
    return

  # "qf, k1 = ..." parameter set

  if len(dlist) > 4 and dlist[0] in common.ele_dict and dlist[1] == ',' and dlist[3] == '=':
    f_out.write(dlist[0] + '[' + bmad_param(dlist[2], dlist[0]) + '] = ' + bmad_expression(''.join(dlist[4:]), dlist[2]))
    return

  # Element def

  if dlist[1] == ':':
    parse_element(dlist, True)
    return

  # Unknown

  print (f"Unknown construct:\n" + command + '\n')

#------------------------------------------------------------------
#------------------------------------------------------------------
# Get next MAD8 command.
# Read in MAD8 file line-by-line.  Assemble lines into commands.

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
        if not common.one_file:
          common.f_out[-1].close()
          common.f_out.pop()       # Remove last file handle
        if len(common.f_in) == 0: return ['', dlist]

    else:
      f_in = common.f_in[-1]
      f_out = common.f_out[-1]
      line = common.command
      common.command = ''

    # Parse line

    if line.strip() == '':
      f_out.write('\n')
      continue

    while line != '':
      for ix in range(len(line)):
        ##print (f'Ix: {ix:4} {line[ix]} -{quote_delim}--{line.rstrip()}')
        ##print (f'C: {command}')
        ##print (f'D: {dlist}')

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

        if line[ix] == '!':
          if len(line) > ix+10 and line[ix:ix+10] == '!!verbatim':
            f_out.write(line[ix+10:].strip() + '\n')
          else:
            f_out.write(line[ix:])
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          line = ''
          if len(dlist) != 0: return [command, dlist]
          break

        elif line[ix] == '&':
          line2 = f_in.readline().lstrip()
          while True:
            if line2[0] == '\n' or line2.lstrip()[0] == '!':
              f_out.write(line2)
              line2 = f_in.readline().lstrip()
            else:
              break
          line = line[:ix] + line2
          break

        elif line[ix] == ';':
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          if len(dlist) == 0:
            line = line[ix+1:]
            break
          else:
            common.command = line[ix+1:]
            return [command, dlist]

        elif line[ix] in ':,=':
          command += line[:ix+1]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          dlist.append(line[ix])
          line = line[ix+1:]
          break

        elif line[ix] == '\n' or ix == len(line)-1:
          command += line[:ix+1]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
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

# Read the parameter file specifying the MAD8 lattice file, etc.

argp = argparse.ArgumentParser()
argp.add_argument('mad8_file', help = 'Name of input MAD8 lattice file')
argp.add_argument('-d', '--debug', help = 'Print debug info (not of general interest).', action = 'store_true')
argp.add_argument('-f', '--many_files', help = 'Create a Bmad file for each MAD8 input file.', action = 'store_true')
argp.add_argument('-v', '--no_prepend_vars', help = 'Do not move variables to the beginning of the Bmad file.', action = 'store_true')
arg = argp.parse_args()

common = common_struct()
common.debug = arg.debug
common.prepend_vars = not arg.no_prepend_vars
common.one_file = not arg.many_files

mad8_lattice_file = arg.mad8_file
bmad_lattice_file = bmad_file_name(mad8_lattice_file)

print ('Input lattice file is:  ' + mad8_lattice_file)
print ('Output lattice file is: ' + bmad_lattice_file)

# Open files for reading and writing

common.f_in.append(open(mad8_lattice_file, 'r'))  # Store file handle
common.f_out.append(open(bmad_lattice_file, 'w'))

f_out = common.f_out[-1]

#------------------------------------------------------------------
# parse, convert and output mad8 commands

common.command = ''  # init

while True:
  [command, dlist] = get_next_command()
  if len(common.f_in) == 0: break
  parse_command(command, dlist)
  if len(common.f_in) == 0: break   # Hit Quit/Exit/Stop statement.

f_out.close()

#------------------------------------------------------------------
# Prepend variables and superposition statements as needed.

f_out = open(bmad_lattice_file, 'r')
lines = f_out.readlines()
f_out.close()

f_out = open(bmad_lattice_file, 'w')
f_out.write (f'!+\n! Translated from MAD8 file: {mad8_lattice_file}\n!-\n\n')

if common.prepend_vars :
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
