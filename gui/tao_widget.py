import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
import sys
import os
sys.path.append(os.environ['ACC_ROOT_DIR'] + '/tao/gui')
from parameters import tao_parameter_dict
from parameters import tao_parameter
import string

class tk_tao_parameter():
  '''
  Takes a tao_parameter (defined in parameters.py) and a tk frame, and creates an object containing the parameter and appropriate tk widget(s) for displaying and modifying the parameter and value
  '''

  def __init__(self, tao_parameter, frame, pipe=0):
    self.param = tao_parameter
    self.tk_label = tk.Label(frame, text=self.param.name)

    if self.param.type in ['STR', 'INT', 'REAL']:
      self.tk_var = tk.StringVar()
      if self.param.value == None:
        self.tk_var.set("")
      else:
        self.tk_var.set(str(self.param.value))
      self.tk_wid = tk.Entry(frame, textvariable=self.tk_var)
    elif self.param.type == 'ENUM':
      self.tk_var = tk.StringVar()
      self.tk_var.set(self.param.value)
      options = enum_fetch(self.param.name,pipe)
      if options == [""]:
        options = [self.param.value]
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
    #Only set if the user has selected an actual file
    if isinstance(filename, str) & (filename != ""):
      self.tk_var.set(filename)
    else:
      self.tk_var.set("Browse...")

#-----------------------------------------------------------------

class d2_data_frame():
  '''
  Sets up an object with the following attributes:
  self.frame = a tk frame that will hold all relevant widgets and labels
  self.name = holds the name of the d2_data this frame represents
  self.d1_data_list = list of all d1_data items contained by this d2_datum
  NOTE: master should be the frame that this frame is gridded into
  (usually a list_frame), while root should be the application root window
  '''

  def __init__(self, master, root, pipe, d2_data_name, u_ix):
    self.frame = tk.Frame(master)
    self.root = root
    self.name = d2_data_name
    self.d1_data_list = []
    self.d1_using_list = [] #holds using information
    self.d1_ix_lb_list = [] #holds index lower bounds
    self.d1_ix_ub_list = [] #holds index upper bounds
    d2_info = pipe.cmd_in("python data_d1_array " + u_ix + "@" + self.name)
    d2_info = d2_info.splitlines()
    for item in d2_info:
      item = item.split(';')
      self.d1_data_list.append(item[3])
      self.d1_using_list.append(item[4])
      self.d1_ix_lb_list.append(int(item[5]))
      self.d1_ix_ub_list.append(int(item[6]))
    tk.Label(self.frame, text=self.name).grid(row=0, column=0, columnspan=2)
    tk.Label(self.frame, text="Indices").grid(row=0, column=2)
    tk.Label(self.frame, text="Using").grid(row=0, column=3)
    for i in range(4):
      self.frame.grid_columnconfigure(i, pad=10)
    for i in range(len(self.d1_data_list)):
      tk.Label(self.frame, text=self.d1_data_list[i]).grid(row=i+1,column=0)
      tk.Button(self.frame, text="View...", command=self.open_d1_callback(self.name, self.d1_data_list[i], pipe, self.d1_ix_lb_list[i], self.d1_ix_ub_list[i], u_ix)).grid(row=i+1,column=1)
      mytext = str(self.d1_ix_lb_list[i]) + ":" + str(self.d1_ix_ub_list[i])
      tk.Label(self.frame, text=mytext).grid(row=i+1, column=2)
      tk.Label(self.frame, text=self.d1_using_list[i]).grid(row=i+1, column=3)

  def open_d1_callback(self, d2_data_name, d1_data_name, pipe, ix_lb, ix_ub, u_ix):
    return lambda : self.open_d1(d2_data_name, d1_data_name, pipe, ix_lb, ix_ub, u_ix)

  def open_d1(self, d2_data_name, d1_data_name, pipe, ix_lb, ix_ub, u_ix):
    '''
    Opens a window with detailed information for d2_data_name.d1_data_name
    '''
    from tao_windows import tao_d1_data_window
    win = tao_d1_data_window(self.root, pipe, d2_data_name + '.' + d1_data_name, u_ix, ix_lb, ix_ub)

#-----------------------------------------------------------------
class d1_data_list_entry():
  '''Creates various tk widgets to display attributes of a single datum.  Takes one line of output from python data_d_array as input.
  '''
  def __init__(self, master, d_list):
    d_list = d_list.split(';')
    self.tk_wids = [] #Holds all the widgets
    self.tk_tao_params = {} #Holds the tk_tao_parameters for items that are not labels
    self.index = d_list[0]
    for i in range(6):
      self.tk_wids.append(tk.Label(master, text=d_list[i]))

    #Meas value
    meas_value = tao_parameter("meas_value", "REAL", "T", d_list[6])
    self.tk_tao_params["meas_value"] = tk_tao_parameter(meas_value, master)
    self.tk_wids.append(self.tk_tao_params["meas_value"].tk_wid)

    #Model value
    model_value = tao_parameter("model_value", "REAL", "F", d_list[7])
    self.tk_tao_params["model_value"] = tk_tao_parameter(model_value, master)
    self.tk_wids.append(self.tk_tao_params["model_value"].tk_wid)

    #Design value
    design_value = tao_parameter("design_value", "REAL", "F", d_list[8])
    self.tk_tao_params["design_value"] = tk_tao_parameter(design_value, master)
    self.tk_wids.append(self.tk_tao_params["design_value"].tk_wid)

    #Useit_opt and Useit_plot
    useit_opt = tao_parameter("useit_opt", "LOGIC", "F", d_list[9])
    self.tk_tao_params["useit_opt"] = tk_tao_parameter(useit_opt, master)
    self.tk_wids.append(self.tk_tao_params["useit_opt"].tk_wid)

    useit_plot = tao_parameter("useit_plot", "LOGIC", "F", d_list[10])
    self.tk_tao_params["useit_plot"] = tk_tao_parameter(useit_plot, master)
    self.tk_wids.append(self.tk_tao_params["useit_plot"].tk_wid)

    #Good User
    good_user = tao_parameter("good_user", "LOGIC", "T", d_list[11])
    self.tk_tao_params["good_user"] = tk_tao_parameter(good_user, master)
    self.tk_wids.append(self.tk_tao_params["good_user"].tk_wid)

    #Weight
    weight = tao_parameter("weight", "REAL", "T", d_list[12])
    self.tk_tao_params["weight"] = tk_tao_parameter(weight, master)
    self.tk_wids.append(self.tk_tao_params["weight"].tk_wid)

class v1_var_list_entry():
  '''
  Creates the widgets needed to display a single variable held in a v1_var array
  v_list should be one row of output from
  python var_v_array [variable name]
  '''
  def __init__(self, master, v_list):
    v_list = v_list.split(';')
    self.tk_wids = [] #Holds all the widgets
    self.tk_tao_params = {} #Holds the tk_tao_parameters for items that are not labels

    self.index = v_list[0]
    for i in range(2):
      self.tk_wids.append(tk.Label(master, text=v_list[i]))

    # Meas value
    meas_value = tao_parameter("meas_value", "REAL", "T", v_list[2])
    self.tk_tao_params["meas_value"] = tk_tao_parameter(meas_value, master)
    self.tk_wids.append(self.tk_tao_params["meas_value"].tk_wid)

    # Model value
    model_value = tao_parameter("model_value", "REAL", "F", v_list[3])
    self.tk_tao_params["model_value"] = tk_tao_parameter(model_value, master)
    self.tk_wids.append(self.tk_tao_params["model_value"].tk_wid)

    # Design value
    design_value = tao_parameter("design_value", "REAL", "F", v_list[4])
    self.tk_tao_params["design_value"] = tk_tao_parameter(design_value, master)
    self.tk_wids.append(self.tk_tao_params["design_value"].tk_wid)

    # Useit_opt
    useit_opt = tao_parameter("useit_opt", "LOGIC", "F", v_list[5])
    self.tk_tao_params["useit_opt"] = tk_tao_parameter(useit_opt, master)
    self.tk_wids.append(self.tk_tao_params["useit_opt"].tk_wid)

    # good_user
    good_user = tao_parameter("good_user", "LOGIC", "T", v_list[6])
    self.tk_tao_params["good_user"] = tk_tao_parameter(good_user, master)
    self.tk_wids.append(self.tk_tao_params["good_user"].tk_wid)

    # Weight
    weight = tao_parameter("weight", "REAL", "T", v_list[7])
    self.tk_tao_params["weight"] = tk_tao_parameter(weight, master)
    self.tk_wids.append(self.tk_tao_params["weight"].tk_wid)

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
  if option_list == []:
    option_list = [""]
  return option_list
