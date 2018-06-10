import os
import sys

BASE = os.environ['ACC_ROOT_DIR']
sys.path.insert (0, BASE + '/tao/python/ctypes')   # directory containing pytao.py

import pytao

tao = pytao.pytao()     # Assumes Tao library is in the standard directory.
bmadfile = BASE + '/tao/examples/cesr/bmad_L9A18A000-_MOVEREC.lat'
tao.init('-lat ' + bmadfile)
tao.cmd('show ele 1')
