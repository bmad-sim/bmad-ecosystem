import os
import sys

os.environ["Test_EV"] = "ZZZ"
exe = sys.argv[1] + 'sim_utils_test'
os.system(exe)
