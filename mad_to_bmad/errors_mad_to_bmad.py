#!/usr/bin/env python

#+
# Script to convert from MADX error data file to Bmad lattice format.
# See the README file for more details
#-

import sys, re, math, argparse, time

if sys.version_info[0] < 3 or sys.version_info[1] < 6:
  raise Exception("Must be using Python 3.6+")

#------------------------------------------------------------------
# Main program.

start_time = time.time()

# Read the parameter file specifying the MADX lattice file, etc.

argp = argparse.ArgumentParser()
argp.add_argument('madx_file', help = 'Name of input MADX lattice file')
argp.add_argument('-d', '--debug', help = 'Print debug info (not of general interest).', action = 'store_true')
arg = argp.parse_args()

debug = arg.debug

madx_lattice_file = arg.madx_file
bmad_lattice_file = madx_lattice_file + '.bmad'

print ('Input lattice file is:  ' + madx_lattice_file)
print ('Output lattice file is: ' + bmad_lattice_file)

# Open files for reading and writing

f_in  = open(madx_lattice_file, 'r')
f_out = open(bmad_lattice_file, 'w')
f_out.write (f'!+\n! Translated from MADX error file to Bmad lattice format by errors_madx_to_bmad.py\n! File: {madx_lattice_file}\n!-\n\n')
f_out.write ('*[scale_multipoles] = F\n')
f_out.write ('expand_lattice\n\n')

count = {}

for line in f_in.readlines():
  if line[0] == '@' or line[0] == '$': continue
  if line[0] == '*':
    param = line[1:].split()
    continue

  word = line.split()
  name = word[0].replace('"','')
  if name in count:
    count[name] += 1
  else:
    count[name] = 1

  for ix, val in enumerate(word[1:]):
    if float(val) == 0: continue
    pname = param[ix+1]
    if pname[0] == 'K' and pname[-2:] == 'SL':
      n = int(pname[1:-2])
      f_out.write(f'{name}##{count[name]}[a{n}] = {val}/{math.factorial(n)}\n')
    elif pname[0] == 'K' and pname[-1:] == 'L':
      n = int(pname[1:-1])
      f_out.write(f'{name}##{count[name]}[b{n}] = {val}/{math.factorial(n)}\n')
    elif pname == 'DX':
      f_out.write(f'{name}##{count[name]}[x_offset] = {val}\n')
    elif pname == 'DY':
      f_out.write(f'{name}##{count[name]}[y_offset] = {val}\n')
    elif pname == 'DS':
      f_out.write(f'{name}##{count[name]}[z_offset] = {val}\n')
    elif pname == 'DTHETA':
      f_out.write(f'{name}##{count[name]}[x_pitch] = {val}\n')
    elif pname == 'DPHI':
      f_out.write(f'{name}##{count[name]}[y_pitch] = {val}\n')
    elif pname == 'DPSI':
      f_out.write(f'{name}##{count[name]}[theta] = {val}\n')
    else:
      print (f'TRANSLATION OF "{pname}" FOR {name} NOT YET IMPLEMENTED. PLEASE CONTACT A BMAD MAINTAINER.')

