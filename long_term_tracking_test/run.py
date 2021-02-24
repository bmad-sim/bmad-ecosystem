import os
import sys

exe = sys.argv[1] + 'long_term_tracking'
os.system(exe + ' sim1.init')
os.system(exe + ' sim2.init')
os.system('cat sim1.dat sim2.dat > output.now')
os.system('rm sim1.dat sim2.dat')
