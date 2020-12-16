import subprocess
import os

out_file = open('output.now', 'w')

#-----------

print ('Test1...')
os.chdir('test1')
results = subprocess.run(['../../../production/bin/tao', '-noplot', '-lat', 'lat.bmad'], stdout=subprocess.PIPE).stdout.decode('utf-8')
os.chdir('..')

if 'contact DCS' in results or 'FATAL' in results:
  out_file.write ('"Bookkeeper1" STR  "BAD"\n')
  print(results)
else:
  out_file.write ('"Bookkeeper1" STR  "GOOD"\n')

#-----------

print ('Test2...')
os.chdir('test2')
results = subprocess.run(['../../../production/bin/tao', '-noplot', '-lat', 'lat.bmad'], stdout=subprocess.PIPE).stdout.decode('utf-8')
os.chdir('..')

if 'contact DCS' in results or 'FATAL' in results:
  out_file.write ('"Bookkeeper2" STR  "BAD"\n')
  print(results)
else:
  out_file.write ('"Bookkeeper2" STR  "GOOD"\n')
