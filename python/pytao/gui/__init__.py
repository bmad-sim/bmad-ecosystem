'''
This package implements a GUI for tao using tkinter, including interactive
plots with matplotlib.  Matplotlib is required for the GUI.
The GUI is organized into one root window from which various other
toplevel windows can be spawned, including plots, element windows, and tables
for data, variables, and lattice elements.  The GUI also supports defining
new data, variables, and plot templates on the fly.
To start the gui, just run

python -m pytao.gui <args>

on the command line.  <args> can be any switches that tao supports.
Additionally, the GUI supports the -gui_init switch, which takes a
gui.init file as its argument (see GUI documentation for more info)
'''
from .main import tao_root_window
