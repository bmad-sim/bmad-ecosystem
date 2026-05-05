import subprocess
import sys

exe1 = sys.argv[1] + 'tao -noplot -lat tao_test.bmad -startup tao_1.startup'

#-----------

print ('Test #1: Symbolic numbers in datum meas fields')
pp = subprocess.Popen(exe1, shell=True)
pp.wait()
