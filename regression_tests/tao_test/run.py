import subprocess
import os
import sys

exe1 = sys.argv[1] + 'tao -noplot -lat tao_test.bmad -startup tao_1.startup'
exe2 = sys.argv[1] + 'tao -noplot -lat small.bmad -startup tao_2.startup'
exe3 = sys.argv[1] + 'tao -noplot -lat small.bmad --rf_on -startup tao_2.startup'

#-----------

print ('Test #1')
results = subprocess.run(exe1, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
print ('Test #2')
results = subprocess.run(exe2, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
print ('Test #3')
results = subprocess.run(exe3, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
