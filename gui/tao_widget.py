import tkinter as tk
import ttk
from tkinter import messagebox
from tkinter import filedialog
import sys
import os
from parameters import tao_parameter_dict
from parameters import tao_parameter
from data_type_list import data_type_list
import string

class tk_tao_parameter():
  '''
  Takes a tao_parameter (defined in parameters.py) and a tk frame,
  and creates an object containing the parameter and appropriate tk widget(s)
  for displaying and modifying the parameter and value
  pipe: the tao_interface object
  data_source: for DAT_TYPE and DAT_TYPE_Z, filters allowed data types
  plot: for DAT_TYPE_Z, the plot where x_axis_type should be checked
  '''

  def __init__(self, tao_parameter, frame, pipe=0, data_source='', plot=''):
    self.param = tao_parameter
    self.pipe = pipe

    if self.param.type == 'DAT_TYPE_Z':
      # Check if operating as ENUM or as DAT_TYPE
      if plot == '': # should never occur
        self.param.type = 'STR' # most fail safe option
      else:
        plot_info = self.pipe.cmd_in('python plot1 ' + plot)
        plot_info = plot_info.splitlines()
        for item in plot_info:
          if item.find('x_axis_type;') == 0:
            x_axis_type = item.split(';')[3]
            break
          else:
            x_axis_type = None
        if x_axis_type == 'data':
          self.param.type = 'DAT_TYPE'
        else:
          self.param.type = 'ENUM_Z'

    if self.param.type in ['STR', 'INT', 'REAL']:
      self.tk_var = tk.StringVar()
      if self.param.value == None:
        self.tk_var.set("")
      else:
        self.tk_var.set(str(self.param.value))
      self.tk_wid = tk.Entry(frame, textvariable=self.tk_var)
    elif self.param.type in ['ENUM', 'ENUM_Z']:
      self.tk_var = tk.StringVar()
      self.tk_var.set(self.param.value)
      if self.param.type == 'ENUM':
        options = enum_fetch(self.param.name,pipe)
      elif self.param.type == 'ENUM_Z': #added to handle DAT_TYPE_Z
        options = enum_fetch('data_type_z', pipe)
        self.param.type = 'ENUM'
      if options == [""]:
        options = [self.param.value]
      if self.param.value == "":
        self.tk_var.set(options[0])
      self.tk_wid = tk.OptionMenu(frame, self.tk_var, *options)
      # Check for and remove num^ from self.param.name
      self.param.name = self.param.name.split('^')[-1]
    elif self.param.type == 'INUM':
      self.tk_var = tk.StringVar()
      self.tk_var.set(self.param.value)
      options = inum_fetch(self.param.name,pipe)
      if options == [""]:
        options = [self.param.value]
      self.tk_wid = tk.OptionMenu(frame, self.tk_var, *options)
      # Check for and remove num^ from self.param.name
      self.param.name = self.param.name.split('^')[-1]
    elif self.param.type == 'FILE':
      self.tk_var = tk.StringVar()
      self.tk_var.set(self.param.value)
      if self.tk_var.get() == "":
        self.tk_var.set("Browse...")
      self.tk_wid = tk.Button(
          frame, textvariable=self.tk_var, command=self.open_file)
    elif self.param.type == 'LOGIC':
      self.tk_var = tk.BooleanVar()
      self.tk_var.set(self.param.value)
      if self.param.can_vary:
        self.tk_wid = tk.Checkbutton(frame, variable=self.tk_var)
      else:
        self.tk_wid = tk.Label(frame, text=str(self.tk_var.get()))
    elif self.param.type == 'REAL_ARR':
      self.tk_var = tk.StringVar()
      val = ""
      for i in range(len(self.param.value)):
        val = val+str(self.param.value[i])
        val = val+';'
      self.tk_var.set(val) #; delimited list of values
      self.tk_wid = tk.Frame(frame) #holds the entries
      self._svar = [] #slave variables
      self._s = [] #slave widgets
      for i in range(len(self.param.value)):
        self._svar.append(tk.StringVar())
        self._svar[i].set(str(self.param.value[i]))
        self._s.append(tk.Entry(self.tk_wid, textvariable=self._svar[i]))
        self._s[i].pack(side="left", expand=1)
        self._svar[i].trace("w", self._update_real_arr)
    elif self.param.type == 'DAT_TYPE':
      self._no_s_refresh = False
      self.tk_var = tk.StringVar() #This is for the entire value
      self.tk_var.set(self.param.value)
      self.tk_wid = tk.Frame(frame) #The "widget" that should be placed by the gui
      if data_source != '':
        self._data_source = data_source
      self._mvar = tk.StringVar() # The master variable
      self._mvar.set((self.tk_var.get()).split('.')[0])
      self._mvar_old = self._mvar.get() # Tracks changes inn self._mvar
      self._m = ttk.Combobox(self.tk_wid, textvariable=self._mvar,
          values = self._get_dat_types(True), state="readonly")
      self._m.bind("<<ComboboxSelected>>", self._s_refresh)
      self._m.pack(side='left', fill='both', expand=1)
      self._svar = [] #list of slave variables
      self._stype = [] #types of slave widgets (needed for input validation)
      self._s = [] # list of slave widgets
      self._traced = False # tracks whether self._trace_cb exists
      self._s_refresh()
      self._trace_cb = self.tk_var.trace('w', self._fill_widgets)
      self._traced = True

    if self.param.type not in ['DAT_TYPE', 'REAL_ARR']:
      self.tk_wid.config(disabledforeground="black")
    else:
      for widget in self._s:
        widget.config(disabledforeground="black")
    self.tk_label = tk.Label(frame, text=self.param.name)
    if not self.param.can_vary:
      if self.param.type not in ['DAT_TYPE', 'REAL_ARR']:
        self.tk_wid.config(state="disabled")
      else:
        for widget in self._s:
          widget.config(state="disabled")

  def _update_real_arr(self, event=None, *args, **kwargs):
    '''
    Takes a real array's slave variables and puts their data
    into the master tk_var
    '''
    new_val = ""
    for i in range(len(self._svar)):
      # Check that the variable contains a float
      try:
        x = float(self._svar[i].get())
        new_val = new_val + str(float(self._svar[i].get()))
      except:
        # Use zero if not a float
        new_val = new_val + str(0.0)
      new_val = new_val + ';'
    self.tk_var.set(new_val)

  def _get_dat_types(self, filter_data_source=False):
    '''
    Parses the data_type_list (from data_type_list.py) and
    returns a list of allowed master data types
    If filter_data_source is set True, data types whose
    data_source does not match self._data_source are removed
    (no difference if self._data_source is not defined)
    '''
    master_list = []
    for item in data_type_list:
      # Filter out data_types not allowed for this data source
      if filter_data_source:
        try:
          if self._data_source not in item['data_source']:
            continue
        except:
          pass
      item = item['name'].split('<')[0]
      if item != "":
        if item[-1] == '.':
          item = item[:-1]
      master_list.append(item)
    return master_list

  def _s_refresh(self, event=None, *args):
    '''
    Clears the existing slave widgets and variables,
    makes new slave widgets and variables, and populates them
    if necessary
    '''
    # Clear the existing slaves
    for item in self._s:
      item.destroy()
    self._svar = []
    self._stype = []
    self._s = []

    try:
      m_ix = (self._get_dat_types()).index(self._mvar.get())
    except ValueError:
      m_ix = 0

    # Make new slave widgets
    dat_dict = data_type_list[m_ix]
    p_start = dat_dict['name'].find('<') - 1
    if p_start != -2:
      slave_params = (dat_dict['name'][p_start:]).split('.')[1:]
    else:
      slave_params = []
    k = 0 #loop counter
    for p in slave_params:
      self._svar.append(tk.StringVar())
      self._stype.append(p)
      if p.find('<enum') != -1: #Enums
        self._s.append(tk.OptionMenu(self.tk_wid, self._svar[k],
          *dat_dict[p]))
        self._svar[k].set(dat_dict[p][0]) # Default to first list item
      elif p.find('<digit:') != -1: #Digit dropdown box
        p = p.split(':')[1]
        p = p.split('>')[0] #p is now "low-high"
        low, high = p.split('-')
        low = int(low)
        high = int(high)
        digits = list(range(low, high+1))
        self._s.append(tk.OptionMenu(self.tk_wid, self._svar[k],
          *digits))
        self._svar[k].set(digits[0])
      elif p.find('<str>') != -1: #Strings
        self._s.append(tk.Entry(self.tk_wid, textvariable=self._svar[k]))
      elif p.find('<digits') != -1: #Fixed length int
        p = p.split('s')[1]
        p = p.split('>')[0]
        length = int(p)
        self._s.append(tk.Entry(self.tk_wid, textvariable=self._svar[k]))
      elif p.find('<int>') != -1: #Integer
        self._s.append(tk.Entry(self.tk_wid, textvariable=self._svar[k]))
      elif p.find('<real>') != -1: #Float
        self._s.append(tk.Entry(self.tk_wid, textvariable=self._svar[k]))
      self._svar[k].trace("w", self._update_tk_var)
      k = k+1

    # Set the slave variables appropriately
    #current_mvar = (self.tk_var.get()).split('<')[0]
    #if current_mvar != "":
    #  if current_mvar[-1] == '.':
    #    current_mvar = current_mvar[:-1]
    #current_mvar = current_mvar.split('.')[0]
    if self._mvar.get() == self._mvar_old:
      for k in range(len(self._svar)):
        # Special case: velocity -> velocity.
        if self.tk_var.get() == 'velocity':
          self._svar[k].set("velocity.".split('.')[k+1])
        else:
          self._svar[k].set((self.tk_var.get()).split('.')[k+1])
    else: # Update self.tk_var if self._mvar has changed
      self._update_tk_var()
      self._mvar_old = self._mvar.get()

    # Pack the new slaves
    for k in range(len(self._s)):
      self._s[k].pack(side="left", fill="both")

  def _update_tk_var(self, event=None, *args, **kwargs):
    '''
    Updates self.tk_var with the current contents of
    self._mvar and self._svar
    '''
    new_tk_var = self._mvar.get()
    for k in range(len(self._svar)):
      # Input validation (TODO)
      p = self._stype[k]
      if p.find('<digits') != -1:
        p = p.split('s')[1]
        p = p.split('>')[0]
        length = int(p)
        #Check if the last character typed was a digit
        if len(self._svar[k].get()) != 0:
          try:
            x = int(self._svar[k].get()[-1])
          except:
            self._svar[k].set(self._svar[k].get()[:-1])
        #Prevent length from being greater than 6 digits
        self._svar[k].set((self._svar[k].get())[-1*length:])
        #Write the corrent number of digits (pad with zeros)
        try: #see if self._svar[k] is an int
          x = int(self._svar[k].get())
          x = len(self._svar[k].get())
          if x < 6:
            new_tk_var = new_tk_var + '.' + (6-x) * '0' + self._svar[k].get()
          else:
            new_tk_var = new_tk_var + '.' + self._svar[k].get()
        except ValueError:
          new_tk_var = new_tk_var + '.' + length * '0'
      elif p == '<int>':
        try: #see if self._svar[k] is an int
          x = int(self._svar[k].get())
          new_tk_var = new_tk_var + '.' + self._svar[k].get()
        except ValueError:
          new_tk_var = new_tk_var + '.' + '0'
      elif p == '<real>':
        try: #see if self._svar[k] is an float
          x = float(self._svar[k].get())
          new_tk_var = new_tk_var + '.' + self._svar[k].get()
        except ValueError:
          new_tk_var = new_tk_var + '.' + '0'
      else:
        new_tk_var = new_tk_var + '.' + self._svar[k].get()
    # Special case: velocity. -> velocity
    if new_tk_var == "velocity.":
      new_tk_var = "velocity"
    # Un-trace tk_var to prevent repeatedly running this method
    self._no_s_refresh = True
    self.tk_var.set(new_tk_var)
    # Re-trace tk_var
    self._no_s_refresh = False
    #self._trace_cb = self.tk_var.trace('w', self._fill_widgets)
    #self._traced = True

  def _is_valid_dat_type(self, x):
    '''
    Returns True if the string x is a valid data_type and False otherwise
    '''
    if not isinstance(x, str):
      return False
    # Special case: normal.h
    if x.find('normal.h') == 0:
      x = ['normal.h'] + x[8:].split('.')
    else:
      x = x.split('.')
    try:
      m_ix = (self._get_dat_types()).index(x[0])
    except ValueError:
      #print('Failed at m_ix')
      return False

    # Check what subvalues are allowed
    dat_dict = data_type_list[m_ix]
    p_start = dat_dict['name'].find('<') - 1
    if p_start != -2:
      slave_params = (dat_dict['name'][p_start:]).split('.')[1:]
    else:
      slave_params = []
    # Reject if length of slave_params and x don't match up
    if len(slave_params) != len(x)-1:
      #print('Failed due to length mismatch')
      #print(slave_params)
      #print(len(x)-1)
      return False

    # Check each subvalue
    k = 0 #loop counter
    for p in slave_params:
      if p.find('<enum') != -1: #Enums
        if x[k+1] not in dat_dict[p]:
          #print('failed because enum not found')
          return False
      elif p.find('<digit:') != -1: #Digit dropdown box
        try:
          x[k+1] = int(x[k+1])
        except ValueError:
          #print('failed bc not int')
          return False
        p = p.split(':')[1]
        p = p.split('>')[0] #p is now "low-high"
        low, high = p.split('-')
        low = int(low)
        high = int(high)
        if (x[k+1] < low) | (x[k+1] > high):
          #print('failed bc out of range')
          return False
      elif p.find('<str>') != -1: #Strings
        if len(x[k+1]) == 0:
          #print('failed bc 0 length')
          return False
      elif p.find('<digits') != -1: #Fixed length int
        try:
          junk_var = int(x[k+1]) #don't need the value
        except ValueError:
          #print('failed bc not an int')
          return False
        p = p.split('s')[1]
        p = p.split('>')[0]
        length = int(p)
        if len(x[k+1]) != length:
          #print('failed bc wrong length')
          return False
      elif p.find('<int>') != -1: #Integer
        try:
          junk_var = int(x[k+1]) #don't need the value
        except ValueError:
          #print('failed because not an int')
          return False
      elif p.find('<real>') != -1: #Float
        try:
          junk_var = float(x[k+1]) #don't need the value
        except ValueError:
          #print('failed bc not a float')
          return False
      k = k+1
    return True # All tests passed

  def _fill_widgets(self, event=None, *args):
    '''
    Runs self._s_refresh() and then fills the slave widgets
    with data from self.tk_var, but only if self.is_valid_dat_type(self.tk_var.get()) is True
    '''
    if not self._is_valid_dat_type(self.tk_var.get()):
      return

    # Refresh slave widgets
    self._mvar.set((self.tk_var.get()).split('.')[0])
    self._mvar_old = self._mvar.get() # gets the slave variables to fill
    if self._no_s_refresh:
      return
    self._s_refresh()

  def _has_ele(self, *args):
    '''
    Returns True if self.tk_var.get() can take associated ele names
    If pipe == 0, returns True by default
    If self.param.type != DAT_TYPE, returns false
    '''
    if self.param.type != 'DAT_TYPE':
      return False
    if self.pipe == 0:
      return True

    cmd_string = 'python datum_has_ele '
    cmd_string += self.tk_var.get()
    if self.pipe.cmd_in(cmd_string) == 'no':
      return False
    else:
      return True

  def _has_s_offset(self, *args):
    '''
    Returns True if self.tk_var.get() can take an s_offset
    Returns False if self.param.type != DAT_TYPE
    '''
    if self.param.type != 'DAT_TYPE':
      return False

    # Locate the correct dictionary in data_type_list
    for data_dict in data_type_list:
      if data_dict['name'].find(self._mvar.get()) != -1:
        found=True
        break
      else:
        found=False

    if not found:
      return True #default

    return data_dict['s_offset']

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
      tk.Button(self.frame, text="View...",
          command=self.open_d1_callback(self.name, self.d1_data_list[i], pipe,
            self.d1_ix_lb_list[i], self.d1_ix_ub_list[i], u_ix)).grid(
                row=i+1,column=1)
      mytext = str(self.d1_ix_lb_list[i]) + ":" + str(self.d1_ix_ub_list[i])
      tk.Label(self.frame, text=mytext).grid(row=i+1, column=2)
      tk.Label(self.frame, text=self.d1_using_list[i]).grid(row=i+1, column=3)

  def open_d1_callback(self, d2_data_name, d1_data_name,
      pipe, ix_lb, ix_ub, u_ix):
    return lambda : self.open_d1(
        d2_data_name, d1_data_name, pipe, ix_lb, ix_ub, u_ix)

  def open_d1(self, d2_data_name, d1_data_name, pipe, ix_lb, ix_ub, u_ix):
    '''
    Opens a window with detailed information for d2_data_name.d1_data_name
    '''
    from tao_windows import tao_d1_data_window
    win = tao_d1_data_window(
        self.root, pipe, d2_data_name + '.' + d1_data_name, u_ix, ix_lb, ix_ub)

#-----------------------------------------------------------------
class d1_data_list_entry():
  '''Creates various tk widgets to display attributes of a single datum.
  Takes one line of output from python data_d_array as input.
  '''
  def __init__(self, master, d_list):
    d_list = d_list.split(';')
    self.tk_wids = [] #Holds all the widgets
    self.tk_tao_params = {} #Holds the tk_tao_parameters for non-label items
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
    self.tk_tao_params = {} #Holds the tk_tao_parameters for non-label items

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
  Takes the name of an enum variable and returns a list of its allowed values
  using the given pipe
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


def inum_fetch(inum,pipe):
  '''
  Takes the name of an inum variable and returns a list of its allowed values
  using the given pipe
  '''
  if pipe != 0:
    list_string = pipe.cmd_in("python inum " + inum)
    option_list = list_string.splitlines()
  else:
    option_list = ["TAO NOT STARTED"]
  if option_list == []:
    option_list = [""]
  return option_list
