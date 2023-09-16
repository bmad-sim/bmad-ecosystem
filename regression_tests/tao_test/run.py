import subprocess
import os
import sys

out_file = open('output.now', 'w')

if os.path.isabs(sys.argv[1]):
  exe = sys.argv[1] + 'tao'
else:
  exe = '../' + sys.argv[1] + 'tao'

#-----------

print ('Test1...')
os.chdir('test1')
results = subprocess.run([exe, '-noplot', '-lat', 'lat.bmad'], stdout=subprocess.PIPE).stdout.decode('utf-8')
os.chdir('..')

if 'contact DCS' in results or 'FATAL' in results:
  out_file.write ('"Bookkeeper1" STR  "BAD"\n')
  print(results)
else:
  out_file.write ('"Bookkeeper1" STR  "GOOD"\n')

#-----------

print ('Test2...')
os.chdir('test2')
results = subprocess.run([exe, '-noplot', '-lat', 'lat.bmad'], stdout=subprocess.PIPE).stdout.decode('utf-8')
os.chdir('..')

if 'contact DCS' in results or 'FATAL' in results:
  out_file.write ('"Bookkeeper2" STR  "BAD"\n')
  print(results)
else:
  out_file.write ('"Bookkeeper2" STR  "GOOD"\n')
