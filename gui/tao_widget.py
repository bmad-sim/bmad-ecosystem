import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
import sys
import os
sys.path.append(os.environ['ACC_ROOT_DIR'] + '/tao/gui')
from parameters import tao_parameter_dict
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

    self.tk_wid.config(disabledforeground="black")
    if not self.param.can_vary:
      self.tk_wid.config(state="disabled")

  def open_file(self):
    filename = filedialog.askopenfilename(title = "Select " + self.param.name)
    self.tk_var.set(filename)

#-----------------------------------------------------------------

class d2_data_frame():
  '''
  Sets up an object with the following attributes:
  self.frame = a tk frame that will hold all relevant widgets and labels
  self.name = holds the name of the d2_data this frame represents
  self.d1_data_list = list of all d1_data items contained by this d2_datum
  '''

  def __init__(self, master, pipe, d2_data_name, u_ix):
    self.frame = tk.Frame(master)
    self.name = d2_data_name
    self.d1_data_list = []
    self.d1_ix_lb_list = [] #holds index lower bounds
    self.d1_ix_ub_list = [] #holds index upper bounds
    d2_info = pipe.cmd_in("python data_d1_array " + u_ix + "@" + self.name)
    d2_info = d2_info.splitlines()
    for item in d2_info:
      item = item.split(';')
      self.d1_data_list.append(item[3])
      self.d1_ix_lb_list.append(int(item[5]))
      self.d1_ix_ub_list.append(int(item[6]))
    tk.Label(self.frame, text=self.name, font=('Helvetica', 12, 'bold')).grid(row=0)
    for i in range(len(self.d1_data_list)):
      tk.Label(self.frame, text=self.d1_data_list[i], font=('Helvetica',12)).grid(row=i+1,column=0)
      tk.Button(self.frame, text="View...", font=('Helvetica',12), command=self.open_d1_callback(self.name, self.d1_data_list[i], pipe, self.d1_ix_lb_list[i], self.d1_ix_ub_list[i], u_ix)).grid(row=i+1,column=1)

  def open_d1_callback(self, d2_data_name, d1_data_name, pipe, ix_lb, ix_ub, u_ix):
    return lambda : self.open_d1(d2_data_name, d1_data_name, pipe, ix_lb, ix_ub, u_ix)

  def open_d1(self, d2_data_name, d1_data_name, pipe, ix_lb, ix_ub, u_ix):
    '''
    Opens a window with detailed information for d2_data_name.d1_data_name
    '''
    from main import tao_d1_data_window
    win = tao_d1_data_window(None, pipe, d2_data_name + '.' + d1_data_name, u_ix, ix_lb, ix_ub)

#-----------------------------------------------------------------
class d1_data_list_entry():
  '''
  Creates various tk widgets to display attributes of a single datum
  '''
  def __init__(self, master, d1_data_name, d1_data_ix, u_ix, pipe):
    self.tk_list = [] # holds all the tk widgets for convenience
    self.index = d1_data_ix
    self.ix_label = tk.Label(master, text=str(d1_data_ix))
    self.tk_list.append(self.ix_label)
    self.name = d1_data_name
    self.name_label = tk.Label(master, text=d1_data_name)
    self.tk_list.append(self.name_label)

    param_list = pipe.cmd_in("python data1 " + str(u_ix) + "@" + d1_data_name + "[" + str(self.index) + "]")
    self.param_dict = tao_parameter_dict(param_list.splitlines())

    self.merit_type = tk.Label(master, text=self.param_dict["merit_type"].value)
    self.tk_list.append(self.merit_type)
    self.ref_ele = tk.Label(master, text=self.param_dict["ele_ref_name"].value)
    self.tk_list.append(self.ref_ele)
    self.start_ele = tk.Label(master, text=self.param_dict["ele_start_name"].value)
    self.tk_list.append(self.start_ele)
    self.ele_name = tk.Label(master, text=self.param_dict["ele_name"].value)
    self.tk_list.append(self.ele_name)

    self.meas = tk_tao_parameter(self.param_dict["meas_value"], master, pipe)
    self.tk_list.append(self.meas.tk_wid)
    self.model = tk_tao_parameter(self.param_dict["model_value"], master, pipe)
    self.tk_list.append(self.model.tk_wid)
    self.design = tk_tao_parameter(self.param_dict["design_value"], master, pipe)
    self.tk_list.append(self.design.tk_wid)

    self.useit_opt = tk_tao_parameter(self.param_dict["useit_opt"], master, pipe)
    self.tk_list.append(self.useit_opt.tk_wid)
    self.useit_plot = tk_tao_parameter(self.param_dict["useit_plot"], master, pipe)
    self.tk_list.append(self.useit_plot.tk_wid)

    self.good_user = tk_tao_parameter(self.param_dict["good_user"],master, pipe)
    self.tk_list.append(self.good_user.tk_wid)
    self.weight = tk_tao_parameter(self.param_dict["weight"],master, pipe)
    self.tk_list.append(self.weight.tk_wid)

#-----------------------------------------------------------------

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
