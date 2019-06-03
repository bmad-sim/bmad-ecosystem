import tkinter as tk
from tkinter import messagebox
from collections import OrderedDict 
import sys
import os
import string

startup_list = [
  'beam_file;STR;T;xx',
  'beam_all_file;STR;T;',
  'beam_init_position_file;STR;T;',
  'building_wall_file;STR;T;',
  'data_file;STR;T;',
  'hook_init_file;STR;T;',
  'init_file;STR;T;',
  'lattice_file;STR;T;',
  'plot_file;STR;T;',
  'startup_file;STR;T;',
  'var_file;STR;T;',
  'slice_lattice;STR;T;',
  'disable_smooth_line_calc;LOGIC;T;F',
  'log_startup;LOGIC;T;F',
  'no_stopping;LOGIC;T;F',
  'rf_on;LOGIC;T;F',
  'tracking_method;ENUM;T;',
]

#-------------------------------------------------

class tao_parameter():

  def __init__(self, param_name, param_type, can_vary, param_value):
    self.name = param_name
    self.type = param_type
    self.can_vary = (can_vary == 'T')

    if param_type == 'STR':
      self.value = param_value
    elif param_type == 'INT':
      self.value = int(param_value)
    elif param_type == 'REAL':
      self.value = float(param_value)
    elif param_type == 'LOGIC':
      self.value = (param_value == 'T')
    elif param_type == 'ENUM':
      self.value = param_value
    else:
      print ('UNKNOWN PARAMETER TYPE: ' + param_type)

  def __str__(self):
    return str(self.value)

  def __repr__(self):
    return self.type + ';' + str(self.can_vary) + ';' + str(self.value)

# An item in the parameter list is a string that looks like:
#        'lattice_file;STR;T;bmad.lat'   

def tao_parameter_dict(param_list):
    this_dict = OrderedDict()
    for param in param_list:
      v = param.split(';')
      this_dict[v[0]] = tao_parameter(v[0], v[1], v[2], v[3])
    return this_dict

#-------------------------------------------------

class vertical_param_frame(tk.Frame):

  def __init__(self, param):
    self.label = tk.Label(self, text = param.name)
    self.label.grid(sticky="W", row=0, column=0)

    if param.type == 'STR':
      pass
    elif param.type == 'INT':
      self.value = int(param_value)
    elif param.type == 'REAL':
      self.value = float(param_value)
    elif param.type == 'LOGIC':
      self.value = (param_value == 'T')
    elif param.type == 'ENUM':
      self.value = param_value


#-------------------------------------------------

param_dict = tao_parameter_dict(startup_list)

#Temporary setup
class tao_root_window(tk.Tk):

  def __init__(self, *args, **kwargs):
    tk.Tk.__init__(self)

    self.title("Tao")
    self.geometry('350x400')
    self.protocol("WM_DELETE_WINDOW", self.quit_cmd)

    # Tao startup

    sys.path.append(os.environ['ACC_ROOT_DIR'] + '/tao/python/tao_pexpect')
    #SHOULD BE $DIST_BASE_DIR/tao/python/tao_pexpect
    #sys.path.append('/home/dcs16/linux_lib/tao/python/tao_pexpect')

    import tao_pipe
    self.pipe = tao_pipe.tao_io()


    # Init GUI

    init_frame = tk.Frame(self, width = 20, height = 30)
    init_frame.pack()
    tk_dict = {} #Items: [tk variable, label, tk object]
    k = 0 #row number counter
    for param, tao_param in param_dict.items():
      #create new entry in tk_dict and display it
      if tao_param.type in ['STR', 'INT', 'REAL']:
        tk_dict[param] = [tk.StringVar()]

        tk_dict[param].append(tk.Label(init_frame,text=param))
        
        tk_dict[param].append(tk.Entry(init_frame, textvariable=tk_dict[param][0]))
      elif tao_param.type == 'LOGIC':
        tk_dict[param] = [tk.BooleanVar()]

        tk_dict[param].append(tk.Label(init_frame,text=param))

        tk_dict[param].append(tk.Checkbutton(init_frame,variable=tk_dict[param][0]))
      elif tao_param.type == 'ENUM':
        tk_dict[param] = [tk.StringVar()]

        tk_dict[param].append(tk.Label(init_frame,text=param))

        options = self.pipe.cmd_in("python enum " + param)
        options = options.split('\r\n')
        for i in range(len(options)):
          options[i] = options[i].strip('0123456789;') #removes number at beginning

        tk_dict[param].append(tk.OptionMenu(init_frame,tk_dict[param][0],*options))

      tk_dict[param][0] = param_dict[param].value
      tk_dict[param][1].grid(row=k,sticky='E')
      tk_dict[param][2].grid(row=k, column=1,sticky='W')
      k = k+1

    def tao_load():
      for param, tao_param in param_dict.items():
        if tao_param.type in ['STR','ENUM']:
          tao_param.value = tk_dict[param][0].get()
        elif tao_param.type == 'INT':
          tao_param.value = int(tk_dict[param][0].get())
        elif tao_param.type == 'REAL':
          tao_param.value = float(tk_dict[param][0].get())
        elif tao_param.type == 'LOGIC':
          tao_param.value = tk_dict[param][0].get()
      print(param_dict)

    b = tk.Button(init_frame, text="Load", command=tao_load)
    b.grid(row=k, columnspan=2)


  def quit_cmd(self, event = ''):
    result = messagebox.askquestion("Quit", "Are You Sure?", icon='warning')
    if result == 'yes':
      sys.exit(0)
    else:
      return
     
root = tao_root_window(sys.argv)
root.mainloop()
