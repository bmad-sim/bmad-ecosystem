import subprocess
import os
import sys

exe1 = sys.argv[1] + 'tao -noplot -lat tao_test.bmad -startup tao_1.startup'
exe2 = sys.argv[1] + 'tao -noplot -lat small.bmad -startup tao_2.startup'
exe3 = sys.argv[1] + 'tao -noplot -lat small.bmad --rf_on -startup tao_2.startup'

#-----------

print ('Test #1')
pp = subprocess.Popen(exe1, shell=True)
pp.wait()

print ('Test #2')
pp = subprocess.Popen(exe2, shell=True)
pp.wait()

print ('Test #3')
pp = subprocess.Popen(exe3, shell=True)
pp.wait()
