import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
import string

class tk_tao_parameter():
  '''
  Takes a tao_parameter (defined in parameters.py) and a 
  tk frame, and creates an object containing the parameter and appropriate tk widget(s) for displaying and modifying the parameter and value
  '''

  def __init__(self, tao_parameter, frame, pipe=0):
    self.param = tao_parameter
    self.tk_label = tk.Label(frame, text=self.param.name)

    if self.param.type in ['STR', 'INT', 'REAL']:
      self.tk_var = tk.StringVar()
      self.tk_var.set(str(self.param.value))
      self.tk_wid = tk.Entry(frame, textvariable=self.tk_var)
    elif self.param.type == 'ENUM':
      self.tk_var = tk.StringVar()
      self.tk_var.set(self.param.value)
      options = enum_fetch(self.param.name,pipe)
      self.tk_wid = tk.OptionMenu(frame, self.tk_var, *options) 
    elif self.param.type == 'FILE':
      self.tk_var = tk.StringVar()
      self.tk_var.set(self.param.value)
      if self.tk_var.get() == "":
        self.tk_var.set("Browse...")
      self.tk_wid = tk.Button(frame, textvariable=self.tk_var, command=self.open_file)
    elif self.param.type == 'LOGIC':
      self.tk_var = tk.BooleanVar()
      self.tk_var.set(self.param.value)
      self.tk_wid = tk.Checkbutton(frame, variable=self.tk_var)

    if not self.param.can_vary:
      self.tk_wid.config(state="disabled")

  def open_file(self):
    filename = filedialog.askopenfilename(title = "Select " + self.param.name)
    self.tk_var.set(filename)

def enum_fetch(enum,pipe):
  '''
  Takes the name of an enum variable and returns a list of its allowed values using the given pipe
  '''
  if pipe != 0:
    list_string = pipe.cmd_in("python enum " + enum)
    option_list = list_string.splitlines()
    for i in range(len(option_list)):
      sc = option_list[i].find(';')
      option_list[i] = option_list[i][sc+1:]
  else:
    option_list = ["TAO NOT STARTED"]
  return option_list
