#!/usr/bin/env python

#+
# Script to convert from Elegant lattice format to Bmad lattice format.
# See the README file for more details.
#-

import sys, re, math, argparse, time
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
    self.ele_dict = {}               # Dict of elements
    self.f_in = []         # Elegant input files
    self.f_out = []        # Bmad output files
    self.use = ''
    self.command = ''

#------------------------------------------------------------------
#------------------------------------------------------------------

ele_type_translate = {
  'beambeam':   'beambeam',
  'ccbend':     'sbend',
  'csbend':     'sbend',
  'csrcsbend':  'sbend',
  'nibend':     'sbend',
  'sben':       'sbend',
  'sbend':      'sbend',
  'rben':       'rbend',
  'rbend':      'rbend',
  'drift':      'drift',
  'csrdrift':   'drift',
  'drif':       'drift',
  'edrift':     'drift',
  'cwiggler':   'wiggler',
  'gfwiggler':  'wiggler',
  'wiggler':    'wiggler',
  'ehkick':     'hkicker',
  'hkick':      'hkicker',
  'vkick':      'vkicker',
  'evkick':     'vkicker',
  'ekicker':    'kicker',
  'kicker':     'kicker',
  'ematrix':    'taylor',
  'hmon':       'monitor',
  'vmon':       'monitor',
  'watch':      'monitor',
  'kquad':      'quadrupole',
  'kquse':      'quadrupole',
  'quad':       'quadrupole',
  'ksext':      'sextupole',
  'sext':       'sextupole',
  'octu':       'octupole',
  'koct':       'octupole',
  'malign':     'patch',
  'mark':       'marker',
  'mult':       'multipole',
  'ecol':       'ecollimator',
  'rcol':       'rcollimator',
  'scraper':    'rcollimator',
  'modrf':      'rfcavity',
  'rfca':       'rfcavity',
  'rfdf':       'rfcavity',
  'tmcf':       'lcavity',
  'shrfdf':     'crab_cavity',
  'sole':       'solenoid',
  'twiss':      'beginning',
  'floor':      'floor_shift',
  'ilmatrix':   'match',
  'rotate':     'marker', #  [exclude_floor != 0, exclude_optics != 0]
               # taylor      [exclude_floor != 0, exclude_optics == 0]
               # floor_shift [exclude_floor == 0, exclude_optics != 0]
               # patch       [exclude_floor == 0, exclude_optics == 0]
}

bmad_param_name = {
  'dx':           'x_offset',
  'dy':           'y_offset',
  'dz':           'z_offset',
  'pitch':        'y_pitch',
  'yaw':          'x_pitch',
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
  'h1':           'h_gap',
  'h2':           'h_gapx',
  'hgap':         'hgap',
  'hgapx':        'hgapx',
  'fint':         'fint',
  'fintx':        'fintx',
  'fint1':        'fint',
  'fint2':        'fintx',
  'xkick':        'hkick',
  'ykick':        'vkick',
  'kick':         'kick',
  'xcenter':      'x_offset',
  'ycenter':      'y_offset',
  'xsize':        'sig_x',
  'ysize':        'sig_y',
  'charge':       'charge',
  'b_max':        'b_max',
  'periods':      'n_period',
  'vertical':     'tilt = pi/2',
  'helical':      'field_calc = helical_model',
  'x_max':        'x_limit',
  'y_max':        'y_limit',
  'frequency':    'rf_frequency:1e6',
  'phase':        'phi0',
  'fse_dipole':   'dg:(%[angle]/%[L])',
  'fse':          'dg:(%[angle]/%[L])',
}

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

def postfix_to_infix(str):

  # Split expression into tokens and recombine numbers like "3.4e-7" into a single token.
  # Also the minus sign in something like "a -1 *" is a unary minus.
  # To avoid splitting on a unary minus, substitute "@" for all unary minusus

  for i in range(len(str)-1):
    if str[i] == '-' and str[i+1] in '0123456789.': str = str[:i] + '@' + str[i+1:]

  tokens = re.split('([/\-\*\+\^ ])', str)
  tokens = [val for val in tokens if val != '' and val != ' ']

  for i in range(len(tokens)):
    if tokens[i][0] == '@': tokens[i] = '-' + tokens[i][1:]

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
  return tokens

#------------------------------------------------------------------
#------------------------------------------------------------------
# Convert from elegant parameter name to bmad parameter name.

def bmad_param(param, ele_name):
  global common, bmad_param_name

  # For the SLAC version there are Rij and Tijk matrix elements

  if len(param) == 2 and param[0] == 'c' and param[1] in '123456':
    return f'tt{param[1:]}'

  if len(param) == 3 and param[0] == 'r' and param[1] in '123456' and param[2] in '123456':
    return f'tt{param[1:]}'

  if len(param) == 4 and param[0] == 't' and param[1] in '123456' and param[2] in '123456' and param[3] in '123456':
    return f'tt{param[1:]}'

  if param not in bmad_param_name: return '?'

  #

  bparam = bmad_param_name[param]

  if ele_name in common.ele_dict:
    bmad_type = common.ele_dict[ele_name].bmad_type
    if bparam == 'tilt' and (bmad_type == 'sbend' or bmad_type == 'rbend'): bparam = 'ref_tilt'

  return bparam

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

def parse_element(dlist):
  global common, ele_type_translate, bmad_param_name

  ele = ele_struct(dlist[0])

  found = False
  for elegant_type in ele_type_translate:
    if elegant_type.startswith(dlist[2]):
      ele.elegant_type = elegant_type
      ele.bmad_type = ele_type_translate[elegant_type]
      found = True
      break

  if not found:
    print (f'{dlist[2].upper()} TYPE ELEMENT WILL BE TRANSLATED TO A DRIFT!')
    ele.elegant_type = dlist[2]
    ele.bmad_type = 'drift'

  params = parameter_dictionary(dlist[4:])

  #

  ele.param = params
  common.ele_dict[dlist[0]] = ele

  line = ele.name + ': ' + ele.bmad_type
  for param in ele.param:
    bparam = bmad_param(param, ele.name)
    if bparam == '?': continue
    if ele.bmad_type == 'drift' and bparam != 'l': continue
    value = params[param]
    if value[0] == '"' or value[0] == "'": value = postfix_to_infix(value[1:-1])[0]
    line += f', {bparam} = {value}'

  ee1 = int(params.get('edge1_effects', '1'))
  ee2 = int(params.get('edge2_effects', '1'))

  if ee1 == 0 and ee2 == 0:
    line += f', fringe_at = no_end'
  elif ee1 != 0 and ee2 == 0:
    line += f', fringe_at = entrance_end'
  elif ee1 == 0 and ee2 != 0:
    line += f', fringe_at = exit_end'

  f_out = common.f_out[-1]
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

  # "% <expression> sto <var>" construct

  if dlist[0] == '%':
    toks = postfix_to_infix(dlist[1])
    if len(toks) != 3 or toks[1] != 'sto':
      print (f'MALFORMED CONSTANT DEFINITION: {command}')
      return
    wrap_write (f'{toks[2]} = {toks[0]}', f_out)
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
    wrap_write(command, f_out)
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

  # #include

  if dlist[0] == '#include':

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

    if line[0] == '%':
      command = line
      dlist = ['%', line[1:].strip()]
      return [command, dlist]

    if line.rstrip()[:9].lower() == '#include:':
      pass

    while line != '':
      for ix in range(len(line)):
        #print (f'Ix: {ix} {len(line)} {line[ix]} -{quote_delim}-|{line.rstrip()}')
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

        if line[ix] == '!':
          if line[ix:].startswith('!!verbatim'):
            f_out.write(line[ix+10:].strip() + '\n')
          else:
            f_out.write(line[ix:])
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          line = ''
          if len(dlist) != 0: return [command, dlist]
          break

        elif line[ix] == '&' or (line[ix] == ',' and ix == len(line.rstrip())-1):
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

        elif line[ix] == ';':
          command += line[:ix]
          if line[:ix].strip() != '': dlist.append(line[:ix].strip().lower())
          if len(dlist) == 0:
            line = line[ix+1:]
            break
          else:
            common.command = line[ix+1:].strip()
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

# Read the parameter file specifying the Elegant lattice file, etc.

argp = argparse.ArgumentParser()
argp.add_argument('elegant_file', help = 'Name of input Elegant lattice file')
argp.add_argument('-d', '--debug', help = 'Print debug info (not of general interest).', action = 'store_true')
argp.add_argument('-f', '--many_files', help = 'Create a Bmad file for each Elegant input file.', action = 'store_true')
arg = argp.parse_args()

common = common_struct()
common.debug = arg.debug
common.one_file = not arg.many_files

elegant_lattice_file = arg.elegant_file
bmad_lattice_file = bmad_file_name(elegant_lattice_file)

print ('*******Note: In beta testing! Please report any problems! **********')
print ('Input lattice file is:  ' + elegant_lattice_file)
print ('Output lattice file is: ' + bmad_lattice_file)

# Open files for reading and writing

common.f_in.append(open(elegant_lattice_file, 'r'))  # Store file handle
common.f_out.append(open(bmad_lattice_file, 'w'))

f_out = common.f_out[-1]

#------------------------------------------------------------------
# parse, convert and output elegant commands

common.command = ''  # init

while True:
  [command, dlist] = get_next_command()
  if len(common.f_in) == 0: break
  parse_command(command, dlist)
  if len(common.f_in) == 0: break   # Hit Quit/Exit/Stop statement.

f_out.close()

#------------------------------------------------------------------

f_out = open(bmad_lattice_file, 'r')
lines = f_out.readlines()
f_out.close()

f_out = open(bmad_lattice_file, 'w')
f_out.write (f'''
!+
! Translated by elegant_to_bmad.py from Elegant file: {elegant_lattice_file}
!-

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

for line in lines:
  f_out.write(line)

f_out.close()
print ('*******Note: In beta testing! Please report any problems! **********')
