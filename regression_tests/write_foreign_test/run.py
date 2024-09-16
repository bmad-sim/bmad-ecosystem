import subprocess
import os
import sys
import difflib

out_file = open('output.now', 'w')

#-----------

exe = sys.argv[1] + 'write_foreign_test'
results = subprocess.run([exe], stdout=subprocess.PIPE).stdout.decode('utf-8')
d = difflib.Differ()

##files = ['write_foreign_test.mad8', 'write_foreign_test.madx', 'write_foreign_test.sad', 
##             'write_foreign_test.lte', 'write_foreign_test.julia', 'write_foreign_test.opal']

files = ['write_foreign_test.mad8', 'write_foreign_test.madx', 
             'write_foreign_test.lte', 'write_foreign_test.julia']

for file in files:
  f1 = open(file + '.correct', 'r')
  lines1 = f1.readlines()

  f2 = open(file + '.now', 'r')
  lines2 = f2.readlines()

  differ = False
  for line in d.compare(lines1, lines2):
    if line[0] == ' ': continue
    if 'Bmad Lattice File' in line: continue
    if not differ:
      print('\n' + file)
    differ = True
    print(line)

  if differ:
    out_file.write ('"' + file + '" STR  "BAD"\n')
    print(results)
  else:
    out_file.write ('"' + file + '" STR  "GOOD"\n')

