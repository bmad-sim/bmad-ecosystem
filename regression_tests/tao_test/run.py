import subprocess
import os
import sys

exe = sys.argv[1] + 'tao'

#-----------

print ('Test...')
#os.chdir('test1')
results = subprocess.run([exe, '-noplot'], stdout=subprocess.PIPE).stdout.decode('utf-8')
