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
  def __init__(self, name = '', type = ''):
    self.name = name
    self.type = type
    self.at = '0'         # Used if element is in a sequence
    self.from_ele = ''    # Used if element is in a sequence
    self.param = OrderedDict()

class seq_struct:
  def __init__(self, name = ''):
    self.name = name
    self.l = '0'
    self.refer = 'centre'
    self.refpos = ''
    self.ele_dict = OrderedDict()

class common_struct:
  def __init__(self):
    self.debug = False
    self.prepend_vars = True
    self.one_file = True
    self.in_seq = False
    self.seqedit_name = ''           # Name of sequence in seqedit construct.
    self.macro_count = False         # Count "{}" braces. 0 -> Not in a "macro" or "if" statement
    self.last_seq = seq_struct()     # Current sequence being parsed.
    self.seq_dict = OrderedDict()    # List of all sequences.
    self.ele_dict = {}               # Dict of elements
    self.set_list = []               # List of "A = B" sets after translation to Bmad. Does not Include "A->P = B" parameter sets.
    self.var_name_list = []          # List of madx variable names.
    self.super_list = []             # List of superimpose statements to be prepended to the bmad file.
    self.f_in = []         # MADX input files
    self.f_out = []        # Bmad output files
    self.use = ''
    self.command = ''    # Scratch storage for get_next_command routine.

#------------------------------------------------------------------
#------------------------------------------------------------------

ele_param_factor = {
  'volt':     '1e-6',
  'freq':     '1e-6',
  'energy':   '1e-9',
  'ex':       '1e-6',
  'ey':       '1e-6',
  'pc':       '1e-9',
}

ele_inv_param_factor = {
  'volt':     '1e6',
  'freq':     '1e6',
  'energy':   '1e9',
  'ex':       '1e6',
  'ey':       '1e6',
  'pc':       '1e9',
}

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

  bmad_param_name = {
    'volt':   'voltage',
    'freq':   'rf_frequency',
    'lag':    'phi0',
    'ex':     'e_field',
    'ey':     'e_field',
    'lrad':   'l',  }

#------------------------------------------------------------------
#------------------------------------------------------------------
# Convert from madx parameter name to bmad parameter name.

def bmad_param(param, ele_name):
  global bmad_param_name

  if param == 'tilt':
    if ele_name in common.ele_dict:
      ele_type = common.ele_dict[ele_name].type
    else:
      ele_type = xxxx

    if 'sbend'.startswith(ele_type) or 'rbend'.startswith(ele_type):
      return 'ref_tilt'
    else:
      return 'tilt'

  elif param in bmad_param_name:
    return bmad_param_name[param]

  elif len(param) == 5 and param[0:4] == 'kick' and param[4] in '123456':
    return f'tt{param[4]}'

  elif len(param) == 4 and param[0:2] == 'rm' and param[2] in '123456' and param[3] in '123456':
    return f'tt{param[2:]}'

  elif len(param) == 5 and param[0:2] == 'tm' and param[2] in '123456' and param[3] in '123456' and param[4] in '123456':
    return f'tt{param[2:]}'

  else:
    return param

#------------------------------------------------------------------
#------------------------------------------------------------------
# Return dictionary of "A = value" parameter definitions.

def parameter_dictionary(word_lst):

  madx_logical = ['kill_ent_fringe', 'kill_exi_fringe', 'thick', 'no_cavity_totalpath']

  # Remove :, {, and } chars for something like "kn := {a, b, c}"
  word_lst = list(filter(lambda a: a not in ['{', '}', ':'], word_lst))

  # replace "0." or "0.0" with "0"
  word_lst = ['0' if x == '0.0' or x == '0.' else x for x in word_lst]

  # Logical can be of the form: "<logical_name>" or "-<logical_name>".
  # Put this into dict

  for logical in madx_logical:
    if logical not in word_lst: continue
    ix = word_lst.index(logical)
    if ix == len(word_lst) - 1: continue
    if word_lst[ix+1] == '=': continue
    pdict[logical] = 'true'
    word_lst.pop(ix)

  for logical in madx_logical:
    if '-' + logical not in word_lst: continue
    pdict[logical] = 'false'
    word_lst.pop(ix)

  # Fill dict
  pdict = OrderedDict()
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
  global const_trans, ele_param_factor

  # Remove {, and } chars for something like "kn := {a, b, c}"
  line = line.replace('{', '').replace('}', '')

  #

  lst = re.split(r'(,|-|\+|\(|\)|\>|\*|/|\^)', line)
  lst = list(filter(lambda a: a != '', lst))    # Remove blank. EG: "->" => ["-", "", ">"] => ["-", ">"]

  out = ''

  while len(lst) != 0:
    if len(lst) >= 4 and lst[1] == '-' and lst[2] =='>':
      if lst[3] in ele_param_factor:
        if (len(lst) >= 5 and lst[4] == '^') or (len(out.strip()) > 0 and out.strip()[-1] == '/'):
          out = out + '(' + lst[0] + '[' + bmad_param(lst[3].strip()) + '] * ' + ele_param_factor[lst[3]]
        else:
          out = out + lst[0] + '[' + bmad_param(lst[3].strip()) + '] * ' + ele_param_factor[lst[3]]
      else:
        out = out + lst[0] + '[' + bmad_param(lst[3].strip()) + ']'
      lst = lst[4:]

    elif lst[0] in const_trans:
      out += const_trans[lst.pop(0)]

    else:
      out += lst.pop(0)

  # End while

  if target_param in ele_inv_param_factor: out = add_parens(out) + ' * ' + ele_inv_param_factor[target_param]
  return out

#-------------------------------------------------------------------
#------------------------------------------------------------------
# Construct the bmad lattice file name

def bmad_file_name(madx_file):

  if madx_file.find('madx') != -1:
    return madx_file.replace('madx', 'bmad')
  elif madx_file.find('Madx') != -1:
    return madx_file.replace('Madx', 'bmad')
  elif madx_file.find('MADX') != -1:
    return madx_file.replace('MADX', 'bmad')
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

def parse_element(dlist, write_to_file, command):
  global common

  ele_type_trans = {
    'tkicker':      'kicker', 
    'hacdipole':    'ac_kicker',
    'hmonitor':     'monitor',
    'vmonitor':     'monitor',
    'placeholder':  'instrument',
    'matrix':       'taylor',
  }

  ignore_madx_param = ['lrad', 'slot_id', 'aper_tol', 'apertype', 'thick', 'add_angle', 'assembly_id', 'mech_sep']
  ele_type_ignore = ['nllens', 'rfmultipole']

  params = parameter_dictionary(dlist[4:])
  ele = ele_struct(dlist[0])

  if dlist[2] == 'dipedge':
    print ('DIPEDGE ELEMENT NOT TRANSLATED. SUGGESTION: MODIFY THE LATTICE FILE AND MERGE THE DIPEDGE ELEMENT WITH THE NEIGHBORING BEND.')
    return

  if dlist[2] in ele_type_ignore:
    print (dlist[2].upper() + ' TYPE ELEMENT CANNOT BE TRANSLATED TO BMAD.')
    return

  elif dlist[2] == 'elseparator':
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

  elif dlist[2] == 'xrotation':
    ele = ele_struct(dlist[0], 'patch')
    if 'angle' in params: params['ypitch'] = params.pop('angle')

  elif dlist[2] == 'yrotation':
    ele = ele_struct(dlist[0], 'patch')
    if 'angle' in params: params['xpitch'] = negate(params.pop('angle'))

  elif dlist[2] == 'srotation':
    ele = ele_struct(dlist[0], 'patch')
    if 'angle' in params: params['tilt'] = params.pop('angle')

  elif dlist[2] == 'changeref':
    ele = ele_struct(dlist[0], 'patch')
    if 'patch_ang' in params:
      angles = params.pop('patch_ang').split(',')
      params['ypitch'] = angles[0]
      params['xpitch'] = negate(angles[1])
      params['tilt']   = angles[2]
    if 'patch_trans' in params:
      trans = params.pop('patch_trans').split(',')
      params['x_offset'] = trans[0]
      params['y_offset'] = trans[1]
      params['z_offset'] = trans[2]

  elif dlist[2] == 'rbend' or dlist[2] == 'sbend':
    ele = ele_struct(dlist[0], dlist[2])
    if 'tilt' in params: params['tilt_ref'] = params.pop('tilt')
    kill_ent = False; kill_exi = False
    if 'kill_ent_fringe' in params: kill_ent = (params['kill_ent_fringe'] == 'true')
    if 'kill_exi_fringe' in params: kill_exi = (params['kill_exi_fringe'] == 'true')
    if kill_ent and kill_exi:
      params['fringe_at'] = 'no_end'
    elif kill_exi:
      params['fringe_at'] = 'entrance_end'
    elif kill_ent:
      params['fringe_at'] = 'exit_end'

  elif dlist[2] == 'quadrupole':
    ele = ele_struct(dlist[0], dlist[2])
    if 'k1' and 'k1s' in params:
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

  elif dlist[2] == 'sextupole':
    ele = ele_struct(dlist[0], dlist[2])
    if 'k2' and 'k2s' in params:
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


  elif dlist[2] == 'octupole':
    ele = ele_struct(dlist[0], dlist[2])
    if 'k3' and 'k3s' in params:
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

  elif dlist[2] == 'multipole':
    ele = ele_struct(dlist[0], dlist[2])
    if 'knl' in params:
      for n, knl in enumerate(params.pop('knl').split(',')): 
        if knl == '0': continue
        params['k' + str(n) + 'l'] = bmad_expression(knl, '')
    if 'ksl' in params:
      for n, ksl in enumerate(params.pop('ksl').split(',')):  
        if ksl == '0': continue
        params['k' + str(n) + 'sl'] = bmad_expression(ksl, '')

  elif dlist[2] in ele_type_trans:
    ele = ele_struct(dlist[0], ele_type_trans[dlist[2]])


  else:
    ele = ele_struct(dlist[0], dlist[2])

  # collimator conversion

  if 'apertype' in params:
    if dlist[2] == 'collimator':
      if params['apertype'] in ['ellipse', 'circle']:
        ele = ele_struct(dlist[0], 'ecollimator')
      else:
        ele = ele_struct(dlist[0], 'rcollimator')

    else:
      if params['apertype'] in ['ellipse', 'circle']:
        params['aperture_type'] = 'elliptical'
      else:
        params['aperture_type'] = 'rectangular'

  #

  if 'aperture' in params: 
    aperture = bmad_expression(params.pop('aperture').replace('{', '').replace('}', ''), '')
    [params['x_limit'], params['y_limit']] = aperture.split(',')[0:2]

  if 'aper_offset' in params:
    params['x_offset'] = params['aper_offset'].split(',')[0]
    params['y_offset'] = params['aper_offset'].split(',')[1]

  #

  if 'at' in params: ele.at = params.pop('at')
  if 'from' in params: ele.from_ele = params.pop('from')

  ele.param = params
  common.ele_dict[dlist[0]] = ele

  if write_to_file:
    line = ele.name + ': ' + ele.type
    for param in ele.param:
      if param in ignore_madx_param: continue
      line += ', ' + bmad_param(param, ele.name) + ' = ' + bmad_expression(params[param], param)
    f_out = common.f_out[-1]
    wrap_write(line, f_out)

  return ele

#------------------------------------------------------------------
#------------------------------------------------------------------

def parse_command(command, dlist):
  global common, sequence_refer

  f_out = common.f_out[-1]

  # Get rid of "real", "int", "const" "const real", etc. prefixes

  if dlist[0].startswith('real ') or dlist[0].startswith('int ') or dlist[0].startswith('const '): dlist[0] = dlist[0].split(' ', 1)[1].strip()
  if dlist[0].startswith('real ') or dlist[0].startswith('int ') or dlist[0].startswith('const '): dlist[0] = dlist[0].split(' ', 1)[1].strip()
  if dlist[0].startswith('shared '): dlist[0] = dlist[0].split(' ', 1)[1].strip()

  # Transform: "a := 3" -> "a = 3"

  for ix in range(len(dlist)):
    if ix > len(dlist)-2: break
    if dlist[ix] == ':' and dlist[ix+1] == '=': dlist.pop(ix)

  # The MADX parser takes a somewhat cavilier attitude towards commas and sometimes they can be omitted.
  # Examples: "call file" instead of  "call, file" and "q: quadrupole l = 7" instead of "q: quadrupole, l = 7
  # So put the comma back in to make things uniform for easier parsing.
  # But do not do this with arithmetical expressions and quoted strings

  if common.debug: print (str(dlist))

  i = 0
  while i < len(dlist):
    if ' ' in dlist[i] and not any(char in dlist[i] for char in '"\'+-*/^'):
      split = dlist[i].split()
      dlist = dlist[:i] + [split[0], ',', split[1]] + dlist[i+1:]
    i += 1

  # Ignore this

  if dlist[0].startswith('exec '): return
  if dlist[0] in ['aperture', 'show', 'value', 'efcomp', 'print', 'select', 'optics', 'option', 'survey',
                  'emit', 'twiss', 'help', 'set', 'eoption', 'system', 'ealign', 'sixtrack', 
                  'flatten']:
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


  # Macro and "if" statements are strange since they does not end with a ';' but with a matching '}'

  if (len(dlist) > 2 and dlist[1] == ':' and dlist[2] == 'macro') or dlist[0].startswith('if '): 
    common.macro_count = command.count('{') - command.count('}')
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

  #

  if dlist[0] == 'endsequence':
    common.in_seq = False
    seq = common.last_seq
    common.seq_dict[seq.name] = seq
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

    f_out.write(f'{dlist[0]}_mark: null_ele\n')
    f_out.write(f'{dlist[0]}_drift: drift, l = {common.last_seq.l}\n')
    f_out.write(f'{dlist[0]}: line = ({dlist[0]}_mark, {dlist[0]}_drift)\n')
    return

  # In a sequence construct.

  if common.in_seq:
    seq = common.last_seq

    if dlist[1] == ':':   # "name: type"  construct
      ele = parse_element(dlist, True, command)
      common.last_seq.ele_dict[ele.name] = ele
      f_out.write(f'superimpose, element = {ele.name}, ref = {seq.name}_mark, ' + \
                  f'offset = {ele.at}, ele_origin = {sequence_refer[seq.refer]}\n')

    elif dlist[0] in common.ele_dict:
      ele = parse_element([dlist[0], ':']+dlist, False, command)
      if len(ele.param) == 0:
        name = ele.name
      # element has modified parameters. Need to create a new element with a unique name with "__N" suffix.
      else:
        common.ele_dict[dlist[0]] = common.ele_dict[dlist[0]] + 1
        name = f'{dlist[0]}__{str(common.ele_dict[dlist[0]])}'
        ele = parse_element([name, ':']+dlist, True, command)

      offset = ele.at
      if ele.from_ele != '':
        from_ele = seq.ele_dict[ele.from_ele]
        offset = f'{offset} - {add_parens(from_ele.at)}'

      f_out.write(f'superimpose, element = {name}, ref = {seq.name}_mark, ' + \
                  f'offset = {offset}, ele_origin = {sequence_refer[seq.refer]}\n')

    # Must be sequence name. So superimpose the corresponding marker.
    else:
      ele = parse_element([dlist[0], ':']+dlist, False, command)

      try:
        seq2 = common.seq_dict[ele.name]
      except:
        print (f'CANNOT IDENTIFY THIS AS AN ELEMENT OR SEQUENCE: {ele.name}\n  IN LINE IN SEQUENCE: {command}')
        return

      offset = ele.at

      if ele.from_ele != '':
        from_ele = seq.ele_dict[ele.from_ele]
        offset = f'{offset} - {add_parens(from_ele.at)}'

      if seq2.refpos != '':
        refpos_ele = seq2.ele_dict[seq2.refpos]
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

  # Var set
  # Do not store variables whose value is an expression tht involves an element parameter

  if dlist[1] == '=' and not '->' in dlist[0]:
    if dlist[0] in common.var_name_list:
      print (f'Duplicate variable name: {dlist[0]}\n' + 
             f'  You will have to edit the lattice file by hand to resolve this problem.')
    common.var_name_list.append(dlist[0])
    name = dlist[0]
    value = bmad_expression(command.split('=')[1].strip(), dlist[0])
    if '[' in value or not common.prepend_vars:    # Involves an element parameter
      f_out.write(f'{name} = {value}\n')
    else:
      common.set_list.append([name, value])
    return

  # Ele param set

  if dlist[1] == '=' and '->' in dlist[0]:
    [ele_name, dummy, param] = dlist[0].partition('->')
    value = bmad_expression(command.split('=')[1].strip(), param)
    name = f'{ele_name}[{bmad_param(param)}]'
    f_out.write(f'{name} = {value}\n')
    return

  # Title

  if dlist[0] == 'title':
    f_out.write(command + '\n')
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
      common.f_out.append(open(bmad_file_name(file), 'r'))
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

  if dlist[0] == 'beam':
    params = parameter_dictionary(dlist[2:])
    if 'particle' in params:  f_out.write('parameter[particle] = ' + bmad_expression(params['particle'], '') + '\n')
    if 'energy' in params:    f_out.write('parameter[E_tot] = ' + bmad_expression(params['energy'], 'energy') + '\n')
    if 'pc' in params:        f_out.write('parameter[p0c] = ' + bmad_expression(params['pc'], 'pc') + '\n')
    if 'gamma' in params:     f_out.write('parameter[E_tot] = mass_of(parameter[particle]) * ' + add_parens(bmad_expression(params['gamma'], '')) + '\n')
    if 'npart' in params:     f_out.write('parameter[n_part] = ' + bmad_expression(params['npart'], '') + '\n')
    return

  # "qf, k1 = ..." parameter set

  if len(dlist) > 4 and dlist[0] in common.ele_dict and dlist[1] == ',' and dlist[3] == '=':
    f_out.write(dlist[0] + '[' + bmad_param(dlist[2]) + '] = ' + bmad_expression(''.join(dlist[4:]), dlist[2]))
    return

  # "qf->k1 = ..." parameter set

  if '->' in dlist[0] and dlist[1] == '=':
    p = dlist[0].split('->')
    f_out.write(p[0] + '[' + bmad_param(p[1]) + '] = ' + bmad_expression(''.join(dlist[2:]), p[1]))
    return

  # Element def

  if dlist[1] == ':':
    parse_element(dlist, True, command)
    return

  # Unknown

  print (f"Unknown construct:\n    " + command.strip())

#------------------------------------------------------------------
#------------------------------------------------------------------
# Get next madx command.
# Read in MADX file line-by-line.  Assemble lines into commands, which are delimited by a ; (colon).

def get_next_command ():
  global common

  quote = ''
  in_extended_comment = False
  command = ''
  f_out = common.f_out[-1]
  dlist = []

  # Loop until a command has been found

  while True:

    # Get a line

    if common.command == '':
      while True:
        f_in = common.f_in[-1]
        line = f_in.readline()
        if len(line) > 0: break    # Check for end of file
        
        common.f_in[-1].close()
        common.f_in.pop()          # Remove last file handle
        if len(common.f_in) == 0: return ['', dlist]

        if not common.one_file and write_to_bmad:
          common.f_out[-1].close()
          common.f_out.pop()       # Remove last file handle
          f_out = common.f_out[-1]
    else:
      line = common.command
      common.command = ''

    # Parse line

    line = line.strip()
    if line.strip() == '':
      f_out.write('\n')
      continue

    while line != '':
      for ix in range(len(line)):
        ##print (f'Ix: {ix} {line[ix]} -{quote}- ' + line)
        ##print (f'{dlist}')

        if common.macro_count != 0:   # Skip macros
          if line[ix] == '{': common.macro_count += 1
          if line[ix] == '}': common.macro_count -= 1
          if common.macro_count == 0:
            line = line[ix+1:]
            continue

        if (line[ix] == '"' or line[ix] == "'") and not in_extended_comment:
          if line[ix] == quote:
            command += quote + line[:ix+1]
            dlist.append(quote + line[:ix+1])
            line = line[ix+1:]
            quote = ''

          else:
            quote = line[ix]
            command += line[:ix]
            if line[:ix].strip() != '': dlist.append(line[:ix].strip())
            line = line[ix+1:]

          break

        if quote != '': continue

        if ix < len(line) - 1 and line[ix:ix+2] == '*/':
          f_out.write(line[:ix] + '\n')
          line = line[ix+2:]
          in_extended_comment = False
          break

        if in_extended_comment: continue

        if line[ix] == '!':
          f_out.write(line[ix:] + '\n')
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip())
          line = ''
          break

        elif line[ix] == ';':
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip())
          common.command = line[ix+1:]
          return [command, dlist]

        elif line[ix] in ':,=':
          command += line[:ix+1]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip())
          dlist.append(line[ix])
          line = line[ix+1:]
          break

        if ix > len(line) - 2: continue

        if line[ix:ix+2] == '/*':
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip())
          f_out.write('!' + line[ix+2:] + '\n')
          line = ''
          in_extended_comment = True
          break

        elif line[ix:ix+2] == '//':
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip())
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
argp.add_argument('-v', '--no_prepend_vars', help = 'Do not move variables to the beginning of the Bmad file.', action = 'store_true')
arg = argp.parse_args()

common = common_struct()
common.debug = arg.debug
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
  [command, dlist] = get_next_command()
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
f_out.write (f'!+\n! Translated from MADX file: {madx_lattice_file}\n!-\n')

if common.prepend_vars :
  for set in common.set_list:
    wrap_write(f'{set[0]} = {set[1]}\n', f_out)
  f_out.write('\n')

if len(common.super_list) > 0:
  for line in common.super_list:
    f_out.write(line)
  f_out.write('\n')

for line in lines:
  f_out.write(line)

f_out.close()
