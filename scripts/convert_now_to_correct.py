#!/usr/bin/python

#+
# Script convert a "output.now" file to a format suitable to "output.correct".
#-

now_file = open('output.now', 'r')  
temp_file = open('output.correct', 'w')

for now_line in now_file:
  if not now_line.strip():     # blank line
    temp_file.write(now_line)    

  elif now_line.strip()[0] == '!': 
    temp_file.write(now_line)    # comment line

  else:
    split = now_line.split('"', 4)

    if split[0].strip() != '':
      print('Cannot parse line: ' + now_line.strip())
      continue

    ## print split

    s2 = split[2].split()
    if len(s2) < 1:
      print('Cannot subparse line: ' + now_line)
      continue

    if (s2[0].strip() == 'STR'):
      temp_file.write('"' + split[1] + '"    "' + split[3] + '"\n')
    else:
      temp_file.write('"' + split[1] + '"   ' + '  '.join(s2[2:]) + '\n')

#

print('Created file: output.correct')
