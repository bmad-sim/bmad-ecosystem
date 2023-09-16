#!/usr/bin/env python3

def ix(n):
  return '{0:0>2}'.format(n)

def kicker(n, pre=''):
  ix1 = ix(n)
  ix2 = ix(2*n)
  ele = 'Pip'+ix2
  s = pre+'Cor'+ix1+': kicker, L = 0.05, field_master = T, superimpose, ref = '+pre+ele
  return s


def bpm(n, pre=''):
  ix1 = ix(n)
  ix2 = ix(2*n-1)
  ele = 'Pip'+ix2
  s = pre+'BPM'+ix1+': marker, superimpose, ref = '+pre+ele
  return s


def write_file(filename, ncells, prefix, func):
  with open(filename, 'w') as f:
    for n in range(1,ncells,1):
      f.write( func(n, prefix)+'\n')

write_file('fa.correctors.bmad', 38, 'FA.', kicker)
write_file('za.correctors.bmad', 20, 'ZA.', kicker)

write_file('fa.bpms.bmad', 38, 'FA.', bpm)
write_file('za.bpms.bmad', 20, 'ZA.', bpm)

