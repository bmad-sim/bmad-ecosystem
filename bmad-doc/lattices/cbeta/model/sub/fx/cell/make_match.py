#!/usr/bin/env python

import sys, re


def process(file, ele):
  with open(file, 'r') as f:
    print('! from', file)
    for x in f:
      y = re.sub(r'.*\[', ele+'[', x.split('\n')[0])
      y = re.sub(r'\]', '1]', y)
      print(y)
      

args = sys.argv
if len(args) < 3:
  print('usage: ./'+args[0], ' beam_start match_ele')
  # Hard coded
  for i in [1,2,3,4]:
    file = 'beam_start'+str(i)+'.bmad'
    ele  = 'FB.MATCH'+str(i)
    process(file, ele)
  sys.exit(0)


file = args[1]
ele = args[2]

process(file, ele)



