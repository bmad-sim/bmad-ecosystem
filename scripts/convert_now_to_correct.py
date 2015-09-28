#!/usr/bin/python

#+
# Script convert a "output.now" file to a format suitable to "output.correct".
#-

import re

p1 = re.compile(r'^(.*"[^"]+" +)(STR|\S+ +\S+) ')

now_file = open('output.now', 'r')  
temp_file = open('output.correct', 'w')

for now_line in now_file:
  if not now_line.strip():     # blank line
    temp_file.write(now_line)    

  elif now_line.strip()[0] == '!': 
    temp_file.write(now_line)    # comment line

  else:

    if p1.match(now_line):
      temp_file.write(p1.sub(r'\1 ', now_line))

    else:
       print('CANNOT SUBPARSE LINE: ' + now_line)

#

print('Created file: output.correct')
