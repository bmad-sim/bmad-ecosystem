import pytao
import os

tao = pytao.pytao()
BASE = os.environ['ACC_ROOT_DIR']
bmadfile = BASE + '/tao/examples/cesr/bmad_L9A18A000-_MOVEREC.lat'
tao.init('-lat ' + bmadfile)
tao.cmd('show ele 1')
