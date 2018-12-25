import os
import sys
TAO_PYTHON_DIR=os.environ['ACC_ROOT_DIR'] + '/tao/python'
sys.path.insert(0, TAO_PYTHON_DIR)
import pytao

tao = pytao.Tao()

LAT = os.environ['ACC_ROOT_DIR'] + '/tao/examples/cesr/bmad_L9A18A000-_MOVEREC.lat'
tao.init('-noinit -lat ' + LAT)    # Should pop up a plotting window.

out = tao.cmd('show ele 10')   # Output is a list
print('\n'.join(out))          # To look nice, print each list item on its own line

out = tao.cmd_real('python plot_line r13.g.a x')  # r13 = beta plot in this example.
print(out[0:5])                                   # Just print a few x values.
