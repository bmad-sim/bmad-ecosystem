#!/usr/bin/env python

from itertools import tee
import sys

#

infile = sys.argv[1]
f = open(infile, 'r')
name_list = {}

for line in f:
  if line[0] == '-': continue
  if '1 END' in line: break
  words = line.split()
  if words[0] not in ['2', '3', '4', '5', '8', '9', '10', '11', '15']: continue
  if '_' in words[1]: continue
  if words[1] in ['RQ', 'CQ']: continue

  name_list[words[1]] = words[0]

  val2 = float(words[2])
  val3 = float(words[3])
  val4 = float(words[4])

  if words[0] == '2':
    print (f'{words[1]}: sbend, l = {val4}, angle = {val2}')
  elif words[0] == '3':
    print (f'{words[1]}: quadrupole, l = {val4}, k1 = {val2/val4:.8f}')
  elif words[0] == '4':
    print (f'{words[1]}: quadrupole, l = {val4}, k1 = {val2/val4:.8f}, tilt')
  elif words[0] == '5':
    print (f'{words[1]}: rfcavity, l = 0,  voltage = {val2} * 1e6')
  elif words[0] == '8':
    print (f'{words[1]}: sextupole, l = {val4}, b2 = {val2/2}')
  elif words[0] == '9':
    print (f'{words[1]}: sbend, l = {val4}, angle = {val2}, ref_tilt = -pi/2')
  elif words[0] == '10':
    print (f'{words[1]}: solenoid, l = {val4},  ks = {val2/val4}')
  elif words[0] == '15':
    print (f'{words[1]}: sbend, l = {val4}, angle = {val2}, k1 = {val3/val4}')
  elif words[0] == '16':
    print (f'{words[1]}: sbend, l = {val4}, angle = {val2}, k1 = {val3/val4}, ref_tilt = -pi/2')

#

for line in f:
  if line[0] == '-': continue
  if line.strip() == '': continue
  words = line.split()
  for name, s in zip(words[0::2], words[1::2]):
    if '_' in name: name = name[:name.index('_')]
    if name not in name_list: continue
    if name_list[name] == '10':  # Solenoid
      print (f'superimpose, element = {name}, ele_origin = beginning, offset = {1e-4*float(s):.4f}')
    else:
      print (f'superimpose, element = {name}, offset = {1e-4*float(s):.4f}')

#

print (f'''
parameter[e_tot] = 7.5e9
rfcavity::*[harmon] = 120
d: drift, l = {1e-4*float(s)}
ring: line = (d)
use, ring
''')
