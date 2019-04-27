import tkinter as tk
from tkinter import messagebox
from tkinter import Menu

import sys

#----------------------------------------------------------------

def read_cmd():
  print ('Read called')

def write_cmd():
  print ('Write called')

def reinit_cmd():
  print ('Reinit called')

def quit_cmd(event = ''):
  result = messagebox.askquestion("Quit", "Are You Sure?", icon='warning')
  if result == 'yes':
    sys.exit(0)
  else:
    return

def optimizer_cmd():
  print ('Optimizer called')
  win = tk.Toplevel()
  win.title('Optimizer')

def plotting_cmd():
  print ('Plotting called')

def wave_cmd():
  print ('Wave called')

def xx_cmd():
  print ('Xx called')

#---------------------------------------------------------------

root = tk.Tk()
root.title("Tao")
root.geometry('350x200')
root.protocol("WM_DELETE_WINDOW", quit_cmd)

menubar = Menu(root)

file_menu = Menu(menubar)
file_menu.add_command(label = 'Read...', command = read_cmd)
file_menu.add_command(label = 'Write...', command = write_cmd)
file_menu.add_command(label = 'Reinit...', command = reinit_cmd)
file_menu.add_separator()
file_menu.add_command(label = 'Quit', command = quit_cmd, accelerator = 'Ctrl+Q')
root.bind_all("<Control-q>", quit_cmd)
menubar.add_cascade(label = 'File', menu = file_menu)

window_menu = Menu(menubar)
window_menu.add_command(label = 'Optimizer...', command = optimizer_cmd)
window_menu.add_command(label = 'Plotting...', command = plotting_cmd)
window_menu.add_command(label = 'Wave...', command = wave_cmd)
menubar.add_cascade(label = 'Window', menu = window_menu)



root.config(menu=menubar)
root.mainloop()
