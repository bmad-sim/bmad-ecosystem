#!/usr/bin/python

import sys, re

#------------------------------------------------------------------
#------------------------------------------------------------------
# Global variables

global_vars = []

#------------------------------------------------------------------
#-------------------------------------------------------------------

def wrap_write(directive, comment):
  MAXLEN = 120
  tab = ''

  if comment != '': directive = directive + '  !' + comment

  while True:
    if len(directive) <= MAXLEN:
      f_out.write(tab + directive + '\n')
      return

    ix = directive[:MAXLEN].rfind(',')

    if ix == -1: 
      ix = directive[:MAXLEN].rfind(' ')
      f_out.write(tab + directive[:ix+1] + ' ...\n')
    else:
      f_out.write(tab + directive[:ix+1] + '\n')  # Don't need '&' after a comma

    tab = '         '
    directive = directive[ix+1:]

#------------------------------------------------------------------
#------------------------------------------------------------------
# Each item of list_out is a list.

def get_arg_list(list_in, directive):

  if list_in[0] != '(' or list_in[-1] != ')':
    print ('ERROR. MISSING "(...)". CANNOT TRANSLATE: ' + directive + '\n')
    print (str(list_in) + '\n')
    wrap_write('???: ' + directive, comment)
    return []

  list_out = []
  this_out = []
  n_parens = 0

  for item in list_in[1:-1]:
    if item == ',' and n_parens == 0:
      list_out.append(this_out)
      this_out = []
    else:
      this_out.append(item)
      if item == '[' or item == '(': n_parens += 1
      if item == ']' or item == ')': n_parens -= 1

  list_out.append(this_out)
  return list_out

#------------------------------------------------------------------
#------------------------------------------------------------------

def print_help():
  print (''' \
Syntax:
  accelerator_toolkit_to_bmad <at_file_name> {<bmad_file_name>}
Default:
  <bmad_file_name> = <at_file_name>.bmad  # If <at_file_name> has ending ".at" this will be stripped.
''')

  sys.exit

#------------------------------------------------------------------
#------------------------------------------------------------------

def parse_directive (directive, comment):

  if directive.strip() == '' and comment.strip() == '':
    f_out.write('\n')
    return

  if directive.strip() == '':
    f_out.write('!' + comment + '\n')
    return

  wordlist = re.split("\s*([\[\]\(\)=, ])\s*", directive)
  wordlist = list(filter(lambda x: x != '' and x != ' ', wordlist))
  if wordlist[-1] == ';': wordlist = wordlist[:-1]
  word0 = wordlist[0]

  if word0 == 'buildlat':
    f_out.write('use, ' + wordlist[2] + '\n') 
    return

  if word0 == 'GLOBVAL.E0':
    f_out.write('parameter[E_tot] = ' + wordlist[2] + '\n')
    return

  if word0 in ['function', 'global', 'fprintf', 'findspos', 'FAMLIST', 'GLOBVAL.LatticeFile', 'clear', 'evalin']:
    wrap_write('! Ignored: ' + directive, comment)
    return

  # Ignore something like "[td,tune] = linopt(...)"
  if word0[0] == '[':    
    wrap_write('! Ignored: ' + directive, comment)
    return

  if wordlist[1] != '=':
    print ('ERROR. CANNOT TRANSLATE: ' + directive)
    wrap_write('???: ' + directive, comment)
    return

  word2 = wordlist[2]

  if word2 in ['version']:
    wrap_write('! Ignored: ' + directive, comment)
    return

  if word2 in ['setcellstruct']: 
    print ('ERROR. CANNOT TRANSLATE: ' + directive)
    wrap_write('???: ' + directive, comment)
    return

  if word2 == 'rfcavity':
    arg_list = get_arg_list(wordlist[3:], directive)
    if len(arg_list) == 0: return
    line = word0 + ': rfcavity, l = ' + ''.join(arg_list[1]) + ', voltage = ' + ''.join(arg_list[2]) + \
           ', rf_frequency = ' + ''.join(arg_list[3])
    wrap_write (line, comment)
    return 

  if word2 == 'quadrupole':
    arg_list = get_arg_list(wordlist[3:], directive)
    if len(arg_list) == 0: return
    line = word0 + ': quadrupole, l = ' + ''.join(arg_list[1]) + ', k1 = ' + ''.join(arg_list[2])
    wrap_write (line, comment)
    return 

  if word2 == 'sextupole':
    arg_list = get_arg_list(wordlist[3:], directive)
    if len(arg_list) == 0: return
    line = word0 + ': sextupole, l = ' + ''.join(arg_list[1]) + ', k2 = ' + ''.join(arg_list[2])
    wrap_write (line, comment)
    return 

  if word2 == 'rbend' or word2 == 'sbend':
    arg_list = get_arg_list(wordlist[3:], directive)
    if len(arg_list) == 0: return
    if word2 == 'rbend':
      line = word0 + ': rbend, l_arc = ' + ''.join(arg_list[1])
    else:
      line = word0 + ': rbend, l = ' + ''.join(arg_list[1])
    if len(arg_list) >= 4: line = line + ', angle = ' + ''.join(arg_list[2])
    if len(arg_list) >= 5: line = line + ', e1 = ' + ''.join(arg_list[3])
    if len(arg_list) >= 6: line = line + ', e2 = ' + ''.join(arg_list[4])
    if len(arg_list) >= 7: line = line + ', k1 = ' + ''.join(arg_list[5])
    wrap_write (line, comment)
    return

  if word2 == 'corrector':
    arg_list = get_arg_list(wordlist[3:], directive)
    if len(arg_list) == 0: return
    kick_list = list(filter(lambda x: x != ',', arg_list[2]))
    line = word0 + ': kicker, l = ' + ''.join(arg_list[1]) + ', hkick = ' + kick_list[1] + ', vkick = ' + kick_list[2]
    wrap_write (line, comment)
    return 

  if word2 == 'marker':
    line = word0 + ': marker'
    wrap_write (line, comment)
    return 

  if word2 == 'aperture':
    arg_list = get_arg_list(wordlist[3:], directive)
    if len(arg_list) == 0: return
    print (str(arg_list) + '\n')
    ap_list = list(filter(lambda x: x != ',', arg_list[1]))
    print (str(ap_list) + '\n')
    line = word0 + ': rcollimator, x1_limit = -(' + ''.join(ap_list[1]) + '), x2_limit = ' + ''.join(ap_list[2]) + \
                                ', y1_limit = -(' + ''.join(ap_list[3]) + '), y2_limit = ' + ''.join(ap_list[4])
    wrap_write (line, comment)
    return 

  if word2 == 'drift':
    arg_list = get_arg_list(wordlist[3:], directive)
    if len(arg_list) == 0: return
    line = word0 + ': drift, l = ' + ''.join(arg_list[1])
    wrap_write (line, comment)
    return 

  if word2 == '[':
    ap_list = list(filter(lambda x: x != ',', wordlist[2:]))
    line = word0 + ': line = (' + ', '.join(ap_list[1:-1]) + ')'
    wrap_write (line, comment)
    return 

  # Assume this is a simple variable definition
  wrap_write (directive, comment)
  return 




#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
# Main program.

# Read the parameter file specifying the SAD lattice file, etc.

if len(sys.argv) == 2:
  at_file = sys.argv[1]
  bmad_file = at_file
  if len(bmad_file) > 2 and bmad_file[-3:] == '.at': bmad_file = bmad_file[:-3]
  bmad_file = bmad_file + '.bmad'

elif len(sys.argv) == 3:
  at_file = sys.argv[1]
  bmad_file = sys.argv[2]

else:
  print_help()

print ('Input AT lattice file is:  ' + at_file)
print ('Output Bmad lattice file is: ' + bmad_file)

f_in = open(at_file, 'r')
f_out = open(bmad_file, 'w')

#------------------------------------------------------------------
# Read in AT file line-by-line.  Assemble lines into directives, which are delimited by a ; (colon).
# Call parse_directive whenever an entire directive has been obtained.

directive = ''
in_comment = False
comment = ''

for line in f_in:
  line = line.strip()

  ixc = line.find('%')
  if ixc != -1:
    comment = line[ixc+1:]
    line = line[:ixc]
    
  directive = directive + line

  ix = directive.find('...')
  if ix != -1: 
    directive = directive[:ix]
    continue

  parse_directive(directive, comment)
  directive = ''
  comment = ''

