import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
import sys
import os
sys.path.append(os.environ['ACC_ROOT_DIR'] + '/tao/gui')
from tao_widget import *
from parameters import str_to_tao_param
import string


#---------------------------------------------------------------
# List window

class tao_list_window(tk.Toplevel):

  def __init__(self, root, title, *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
    self.title(title)

    self.geometry('400x600')

    self.outer_frame=tk.Frame(self)

    canvas=tk.Canvas(self.outer_frame)
    self.list_frame=tk.Frame(canvas)
    scrollbar=tk.Scrollbar(self.outer_frame,orient="vertical",command=canvas.yview)
    canvas.configure(yscrollcommand=scrollbar.set)

    def scrollhelper(event):
      canvas.configure(scrollregion=canvas.bbox("all"),width=200,height=200)
    self.list_frame.bind("<Configure>",scrollhelper)

    def mouse_scroll(event):
      #canvas.yview_scroll(direction,"units")
      if event.num == 4:
        canvas.yview_scroll(-1,"units")
      elif event.num == 5:
        canvas.yview_scroll(1,"units")
      #canvas.update_idletasks()

    def bind_mouse(event):
      self.outer_frame.bind_all("<Button-4>", mouse_scroll)
      self.outer_frame.bind_all("<Button-5>", mouse_scroll)

    def unbind_mouse(event):
      self.outer_frame.unbind_all("<Button-4>")
      self.outer_frame.unbind_all("<Button-5>")

    self.outer_frame.bind("<Enter>", bind_mouse)
    self.outer_frame.bind("<Leave>", unbind_mouse)

    self.outer_frame.pack(side="top",fill="both",expand=1)
    scrollbar.pack(side="right",fill="y")
    canvas.pack(side="left",fill="both",expand=1)
    canvas.create_window((0,0),window=self.list_frame,anchor='nw')


#---------------------------------------------------------------
# Parameter window
class tao_parameter_window(tao_list_window):

  def __init__(self, root, title, tao_list, pipe, *args, **kwargs):
    tao_list_window.__init__(self, root, title, *args, **kwargs)
    self.button_frame = tk.Frame(self)
    self.button_frame.pack(side="top", fill="both", expand=0)
    self.tao_list = tao_list
    for k in range(len(self.tao_list)):
      self.tao_list[k] = tk_tao_parameter(self.tao_list[k], self.list_frame, pipe)
      tk.Label(self.list_frame,text=self.tao_list[k].param.name).grid(row=k,column=0,sticky="E")
      self.tao_list[k].tk_wid.grid(row=k,column=1,sticky="W")
      k = k+1


#---------------------------------------------------------------
# d2_data window

class tao_d2_data_window(tao_list_window):

  def __init__(self, root, pipe, *args, **kwargs):
    tao_list_window.__init__(self, root, "Data", *args, **kwargs)
    self.pipe = pipe
    univ_list = self.pipe.cmd_in("python super_universe")
    n_universe = str_to_tao_param(univ_list.splitlines()[0])
    self.univ_frame = tk.Frame(self)
    tk.Label(self.univ_frame, text="Universe: ", font=('Helvetica',12)).grid(row=0,column=0,sticky="E")
    self.u_ix = tk.StringVar()
    u_ix_list = []
    for i in range(n_universe.value):
      u_ix_list.append(str(i+1))
    self.u_ix.set(u_ix_list[0])
    u_ix_box = tk.OptionMenu(self.univ_frame, self.u_ix, *u_ix_list, command=self.refresh)
    u_ix_box.grid(row=0,column=1,sticky="W")
    self.univ_frame.pack(fill="both", expand=0)

    # Populate self.list_frame
    self.refresh(self.u_ix.get())


  def refresh(self,u_ix):
    '''
    Clears self.list_frame and fills it with the current universe's
    d2/d1 data
    '''
    # Clear self.list_frame:
    for child in self.list_frame.winfo_children():
      child.destroy()
    # Get this universe's d2_data
    d2_data_list = self.pipe.cmd_in("python data_d2_array " + u_ix)
    d2_data_list = d2_data_list.splitlines()
    for d2_data_item in d2_data_list:
      new_frame = d2_data_frame(self.list_frame, self.pipe, d2_data_item, u_ix)
      new_frame.frame.pack()


#---------------------------------------------------------------
# d1_data window

class tao_d1_data_window(tao_list_window):

  def __init__(self, root, pipe, d1_data_name, u_ix, ix_lb, ix_ub, *args, **kwargs):
    tao_list_window.__init__(self, root, d1_data_name, *args, **kwargs)
    self.geometry('1550x600')
    self.pipe = pipe
    self.d1_data_name = d1_data_name
    self.u_ix = u_ix
    self.ix_lb = ix_lb
    self.ix_ub = ix_ub
    self.refresh()

    self.button_frame = tk.Frame(self)
    self.button_frame.pack(side="bottom", fill="both", expand=0)
    b1 = tk.Button(self.button_frame, text="Apply Changes", command=self.apply)
    b2 = tk.Button(self.button_frame, text="Discard Changes and Refresh", command=self.refresh)

    b1.pack(side="left")
    b2.pack(side="right")

  def refresh(self):
    # Clear self.list_frame:
    for child in self.list_frame.winfo_children():
      child.destroy()
    # Row titles
    title_list = ["Index",
        "d1_data_name",
        "Merit type",
        "Ref Element",
        "Start Element",
        "Element Name",
        "Meas value",
        "Model value",
        "Design value",
        "Useit_opt",
        "Useit_plot",
        "good_user",
        "Weight"]
    j = 0
    for item in title_list:
      tk.Label(self.list_frame, text=item).grid(row=0, column=j)
      self.list_frame.grid_columnconfigure(j, pad=10)
      j=j+1

    #Bulk editing
    #Only meas_value, good user, and weight can vary
    tk.Label(self.list_frame, text="Bulk editing:").grid(row=1, column=0, columnspan=6)
    tk.Label(self.list_frame, text="Click to fill:").grid(row=2, column=0, columnspan=6)

    self.bulk_params = [] #Holds the current bulk edit widgets
    self.bulk_filled = [] #Holds whether or not bulk filling has been used
    self.bulk_value = [] #Holds the last value that was filled to the cells
    self.bulk_apply = [] #Holds fill buttons

    self.bulk_params.append(tk_tao_parameter(str_to_tao_param("meas_value;REAL;T;"), self.list_frame))
    self.bulk_params[0].tk_wid.grid(row=1, column=6)
    self.bulk_filled.append(False)
    self.bulk_value.append(self.bulk_params[0].tk_var.get())
    self.bulk_apply.append(tk.Button(self.list_frame, text="Fill...", command= lambda : self.fill(0)))
    self.bulk_apply[0].grid(row=2, column=6)

    self.bulk_params.append(tk_tao_parameter(str_to_tao_param("good_user;LOGIC;T;"), self.list_frame))
    self.bulk_params[1].tk_wid.grid(row=1, column=11)
    self.bulk_filled.append(False)
    self.bulk_value.append(self.bulk_params[1].tk_var.get())
    self.bulk_apply.append(tk.Button(self.list_frame, text="Fill...", command= lambda : self.fill(1)))
    self.bulk_apply[1].grid(row=2, column=11)

    self.bulk_params.append(tk_tao_parameter(str_to_tao_param("weight;REAL;T;"), self.list_frame))
    self.bulk_params[2].tk_wid.grid(row=1, column=12)
    self.bulk_filled.append(False)
    self.bulk_value.append(self.bulk_params[2].tk_var.get())
    self.bulk_apply.append(tk.Button(self.list_frame, text="Fill...", command= lambda : self.fill(2)))
    self.bulk_apply[2].grid(row=2, column=12)

    #Fetch and fill in the data
    self.list_rows = []
    d_list = self.pipe.cmd_in("python data_d_array " + self.u_ix + '@' + self.d1_data_name)
    d_list = d_list.splitlines()
    #i = row counter, j = column counter
    #grid to row i+3 because row 0 is titles, row 1 is bulk editing widgets, row 2 is apply checkboxes
    for i in range(self.ix_ub - self.ix_lb):
      self.list_rows.append(d1_data_list_entry(self.list_frame, d_list[i]))
      for j in range(len(self.list_rows[i].tk_wids)):
        self.list_rows[i].tk_wids[j].grid(row=i+3, column=j)
      tk.Button(self.list_frame, text="View More...", command=self.open_datum_window_callback(i+self.ix_lb)).grid(row=i+3, column=j+1)

  def fill(self, index):
    '''
    Fills the meas_value, good_user, or weight column (determined by index) to the bulk edit value, saves the bulk edit value for efficient calls to tao_set, and clears the bulk edit box
    '''
    # Save the bulk parameter state
    self.bulk_value[index] = self.bulk_params[index].tk_var.get()
    self.bulk_filled[index] = True
    # Clear the bulk parameter widget
    if index != 1:
      self.bulk_params[index].tk_var.set("")
    else:
      self.bulk_params[index].tk_var.set(False)
    # Fill the appropriate variable
    p_names = ["meas_value", "good_user", "weight"]
    for i in range(self.ix_ub - self.ix_lb):
      self.list_rows[i].tk_tao_params[p_names[index]].tk_var.set(self.bulk_value[index])


  def apply(self):
    #Apply bulk changes
    for i in range(len(self.bulk_params)):
      print(self.bulk_filled[i])
      if self.bulk_filled[i]:
        set_str = "set data " + self.u_ix + '@' + self.d1_data_name + '|'
        self.bulk_params[i].tk_var.set(self.bulk_value[i])
        tao_set([self.bulk_params[i]], set_str, self.pipe, overide=(i==1)) #overide is necessary for good_user

    #Apply individual changes that are different from bulk changes
    for i in range(self.ix_ub - self.ix_lb):
      set_list = []
      set_str = "set data " + self.u_ix + '@' + self.d1_data_name + '[' + str(i+self.ix_lb) + ']|'
      #Meas value
      c1 = (self.list_rows[i].tk_tao_params["meas_value"].tk_var.get() != self.bulk_value[0])
      c2 = not self.bulk_filled[0]
      try:
        c3 = (float(self.list_rows[i].tk_tao_params["meas_value"].tk_var.get()) != self.list_rows[i].tk_tao_params["meas_value"].param.value)
      except ValueError:
        c3 = False
      if (c1 | c2) & c3:
        set_list.append(self.list_rows[i].tk_tao_params["meas_value"])

      #Good user
      c1 = (self.list_rows[i].tk_tao_params["good_user"].tk_var.get() != self.bulk_value[1])
      c2 = not self.bulk_filled[1]
      try:
        c3 = (bool(self.list_rows[i].tk_tao_params["good_user"].tk_var.get()) != self.list_rows[i].tk_tao_params["good_user"].param.value)
      except ValueError:
        c3 = False
      if (c1 | c2) & c3:
        set_list.append(self.list_rows[i].tk_tao_params["good_user"])

      #Weight
      c1 = (self.list_rows[i].tk_tao_params["weight"].tk_var.get() != self.bulk_value[2])
      c2 = not self.bulk_filled[2]
      try:
        c3 = (float(self.list_rows[i].tk_tao_params["weight"].tk_var.get()) != self.list_rows[i].tk_tao_params["weight"].param.value)
      except ValueError:
        c3 = False
      if (c1 | c2) & c3:
        set_list.append(self.list_rows[i].tk_tao_params["weight"])

      if set_list != []:
        tao_set(set_list, set_str, self.pipe)

    #Refresh
    self.refresh()

  def open_datum_window_callback(self, d1_data_ix):
    return lambda : self.open_datum_window(d1_data_ix)

  def open_datum_window(self, d1_data_ix):
    param_list = self.pipe.cmd_in("python data " + str(self.u_ix) + '@' + self.d1_data_name + '[' + str(d1_data_ix) + ']')
    param_list = param_list.splitlines()
    for i in range(len(param_list)):
      param_list[i]=str_to_tao_param(param_list[i])
    win = tao_parameter_window(None, self.d1_data_name + '[' + str(d1_data_ix) + ']', param_list, self.pipe)

    set_str = "set data " + str(self.u_ix) + '@' + self.d1_data_name + '[' + str(d1_data_ix) + ']|'
    b = tk.Button(win.button_frame, text="Apply changes", command=lambda : self.datum_set_callback(win.tao_list,set_str))
    b.pack()

  def datum_set_callback(self, tao_list, set_str):
    tao_set(tao_list, set_str, self.pipe)
    self.refresh()



#---------------------------------------------------------------
# Variable Window

#class tao_var_general_window(tk.Toplevel):

#---------------------------------------------------------------
# Root window

class tao_root_window(tk.Tk):

  def __init__(self, *args, **kwargs):
    tk.Tk.__init__(self)

    self.title("Tao")
    self.protocol("WM_DELETE_WINDOW", self.quit_cmd)
    self.tk.call('tk', 'scaling', 1.0)

    # Menu bar

    self.menubar = tk.Menu(self)

    file_menu = tk.Menu(self.menubar)
    file_menu.add_command(label = 'Read...', command = self.read_cmd)
    file_menu.add_command(label = 'Write...', command = self.write_cmd)
    file_menu.add_command(label = 'Reinit...', command = self.reinit_cmd)
    file_menu.add_separator()
    file_menu.add_command(label = 'Quit', command = self.quit_cmd, accelerator = 'Ctrl+Q')
    self.menubar.add_cascade(label = 'File', menu = file_menu)

    window_menu = tk.Menu(self.menubar)
    window_menu.add_command(label = 'Optimizer...', command = self.optimizer_cmd)
    window_menu.add_command(label = 'Plotting...', command = self.plotting_cmd)
    window_menu.add_command(label = 'Wave...', command = self.wave_cmd)
    window_menu.add_command(label = 'Global Variables...', command = self.set_global_vars_cmd)
    window_menu.add_command(label = 'Data...', command = self.view_data_cmd)
    self.menubar.add_cascade(label = 'Window', menu = window_menu)

    self.config(menu=self.menubar)

    # Init GUI

    init_frame = tk.Frame(self, width = 20, height = 30)
    init_frame.pack()
    self.tao_load(init_frame)

    # Key bindings

    self.bind_all("<Control-q>", self.quit_cmd)

  #-------------------------
  # Tao startup

  def start_main(self):
    self.menubar.entryconfig("File", state="normal")
    self.menubar.entryconfig("Window", state="normal")
    self.main_frame = tk.Frame(self, width = 20, height = 30)
    self.main_frame.pack()
    self.bind_all('<Control-g>', self.global_vars_event)

    # Call
    tk.Label(self.main_frame, text="Call command file:").grid(row=0, column=0)
    self.call_file = tk_tao_parameter(str_to_tao_param("call_file;FILE;T;"), self.main_frame, self.pipe)
    self.call_file.tk_wid.grid(row=0, column=1)
    tk.Button(self.main_frame, text="Run", command=self.tao_call).grid(row=0, column=2)

    # Optimization
    tk.Label(self.main_frame, text="Optimization:").grid(row=1, column=0)
    tk.Button(self.main_frame, text="Setup").grid(row=1, column=1)
    tk.Button(self.main_frame, text="Run").grid(row=1, column=2)

    # Command line
    self.cmd_frame = tk.Frame(self)
    self.cmd_frame.pack(side="bottom", fill="x")
    self.history = [] #holds the history
    self.history.append([]) #Tao history
    self.history.append([]) #Shell history
    self.history.append([]) #Call history
    tk.Button(self.cmd_frame, text="View History...").pack(side="top", fill="x")

    self.command = tk_tao_parameter(str_to_tao_param("command;STR;T;"), self.cmd_frame, self.pipe)
    self.command.tk_wid.bind("<Return>", self.tao_command)
    self.command.tk_wid.pack(side="left", fill="x", expand=1)
    tk.Button(self.cmd_frame, text="Run (System Shell)", command=self.tao_spawn).pack(side="right")
    tk.Button(self.cmd_frame, text="Run (Tao)", command=self.tao_command).pack(side="right")

  def tao_command(self, event=None):
    '''
    Runs the text in self.command at the Tao command line, appends it to the history, and clears self.command
    '''
    if self.command.tk_var.get() != "":
      self.pipe.cmd(self.command.tk_var.get())
      self.history[0].append(self.command.tk_var.get())
      self.command.tk_var.set("")

  def tao_spawn(self):
    '''
    Runs the text in self.command at the system command line, appends it to the history, and clears self.command
    '''
    cmd_txt = self.command.tk_var.get()
    self.pipe.cmd("spawn " + cmd_txt)
    self.history[1].append(cmd_txt)
    self.command.tk_var.set("")

  def tao_call(self):
    '''
    Runs the command file in self.call_file, appends it to the history, and clears self.call_file
    '''
    if self.call_file.tk_var.get() != "Browse...":
      self.pipe.cmd("call " + self.call_file.tk_var.get())
      self.history[2].append(self.call_file.tk_var.get())
      self.call_file.tk_var.set("Browse...")

  def tao_load(self,init_frame):
    self.menubar.entryconfig("File", state="disabled")
    self.menubar.entryconfig("Window", state="disabled")
    from parameters import param_dict
    tk_list = [] #Items: tk_tao_parameter()'s (see tao_widget.py)
    init_frame.grid_columnconfigure(0, weight=1, uniform="test", pad=10)
    init_frame.grid_columnconfigure(1, weight=1, uniform="test")

    #Look for and read gui.init
    #gui.init should be in the same directory that
    #main.py is run from
    #each line should have the form
    #parameter:value
    #E.g. init_file:tao.init
    #E.g. rf_on:T
    try:
      init_file = open('gui.init')
      init_list = init_file.read()
      init_list = init_list.splitlines()
    except:
      init_list = []
    init_dict = {}
    for entry in init_list:
      entry = entry.strip()
      #Check for proper formatting and not a comment
      if entry.find(':') != -1 & entry.find('#') != 0:
        name = entry[:entry.find(':')]
        name = name.strip()
        value = entry[entry.find(':')+1:]
        value = value.strip()
        init_dict[name] = value
    k = 0 #row number counter
    for param, tao_param in param_dict.items():
      tk_list.append(tk_tao_parameter(tao_param,init_frame))
      #Possibly set value from init file
      if tk_list[k].param.name in init_dict:
        if tk_list[k].param.type == 'FILE':
          #Check that a good filename has been given
          try:
            filename = init_dict[tk_list[k].param.name]
            #Expand environment variables and ~
            filename = os.path.expanduser(filename)
            filename = os.path.expandvars(filename)
            f = open(filename)
            f.close()
            tk_list[k].tk_var.set(filename)
          except:
            messagebox.showwarning(tk_list[k].param.name, "File not found: " + filename)
        elif tk_list[k].param.type == 'LOGIC':
          state = (init_dict[tk_list[k].param.name] == 'T') | (init_dict[tk_list[k].param.name] == 'True')
          tk_list[k].tk_var.set(state)
      tk.Label(init_frame,text=param).grid(row=k,sticky="E")
      tk_list[k].tk_wid.grid(row=k, column=1, sticky="W")
      k = k+1

    def param_load(event=None):
      init_args = ""
      tao_exe = ""
      for tk_param in tk_list:
        if tk_param.param.type in ['STR','ENUM']:
          tk_param.param.value = tk_param.tk_var.get()
          if tk_param.param.value != "":
            init_args = init_args + "-" + tk_param.param.name + " " + tk_param.param.value + " "
        elif tk_param.param.type == 'FILE':
          tk_param.param.value = tk_param.tk_var.get()
          if tk_param.param.value == "Browse...":
            tk_param.param.value = ""
          if tk_param.param.value != "":
            if tk_param.param.name == "tao_exe":
              tao_exe = tk_param.param.value
            else:
              init_args = init_args + "-" + tk_param.param.name + " " + tk_param.param.value + " "
        elif tk_param.param.type == 'LOGIC':
          tk_param.param.value = tk_param.tk_var.get()
          if tk_param.param.value:
            init_args = init_args + "-" + tk_param.param.name + " "
      # Run Tao, clear the init_frame, and draw the main frame

      from tao_interface import tao_interface
      self.pipe = tao_interface(init_args, tao_exe)
      #sys.path.append(os.environ['ACC_ROOT_DIR'] + '/tao/python/tao_pexpect')
      #import tao_pipe
      #self.pipe = tao_pipe.tao_io(init_args)

      init_frame.destroy()
      self.start_main()

    load_b = tk.Button(init_frame, text="Start Tao", command=param_load)
    load_b.grid(row=k, columnspan=2)
    self.bind_all("<Return>", param_load)

  #-------------------------
  # Menu bar callbacks

  def read_cmd(self):
    print ('Read called')

  def write_cmd(self):
    print ('Write called')

  def reinit_cmd(self):
    print ('Reinit called')

  def quit_cmd(self, event = ''):
    result = messagebox.askquestion("Quit", "Are You Sure?", icon='warning')
    if result == 'yes':
      sys.exit(0)
    else:
      return

  def optimizer_cmd(self):
    print ('Optimizer called')
    win = tk.Toplevel(self)
    win.title('Optimizer')

  def plotting_cmd(self):
    print ('Plotting called')

  def wave_cmd(self):
    print ('Wave called')

  def set_global_vars_cmd(self):
    global_list = root.pipe.cmd_in("python global")
    global_list = global_list.splitlines()
    for i in range(len(global_list)):
      global_list[i]=str_to_tao_param(global_list[i])
    win = tao_parameter_window(None, "Global Variables", global_list, self.pipe)

    b = tk.Button(win.button_frame, text="Set Global Variables", command=lambda : tao_set(win.tao_list, "set global ", self.pipe))
    b.pack()

  def view_data_cmd(self):
    win = tao_d2_data_window(None, self.pipe)

  # Other callbacks

  def global_vars_event(self, event):
    self.set_global_vars_cmd()

#---------------------------------------------------

def tao_set(tao_list,set_str,pipe, overide=False):
  '''
  Takes a list of tk_tao_parameters and makes a call to tao
  to set the parameters to the values input by the user
  set_str should be "set global ", "set data orbit.x[10]|", or whatever is appropriate
  Use the overide option to run set commands even if no change has been made.
  '''
  # Exit imediately if tao_list is empty
  if tao_list == []:
    return 0
  # STOP lattice calculation here
  pipe.cmd_in("set global lattice_calc_on = F")
  pipe.cmd_in("set global plot_on = F")
  #Freeze input fields:
  for item in tao_list:
    item.tk_wid.config(state="disabled")
  update_dict = {} #Record of which variables were changed
  set_lattice_calc = False
  set_plot = False
  for item in tao_list:
    #Type casting and validation
    if item.tk_var.get() != "": #Skip empty boxes
      if item.param.type == 'INT':
        try:
          new_val = int(item.tk_var.get())
        except ValueError:
          messagebox.showwarning("Error",item.param.name + " must be an integer ")
          new_val = item.param.value
      elif item.param.type == 'REAL':
        try:
          new_val = float(item.tk_var.get())
        except ValueError:
          messagebox.showwarning("Error",item.param.name + " must be a real number")
          new_val = item.param.value
      else:
        new_val = item.tk_var.get()
    else:
      new_val = item.param.value
    #Check for any change
    if new_val != item.param.value:
      item.param.value = new_val
      update_dict[item.param.name] = True
    else:
      update_dict[item.param.name] = overide

    #Wait till the end to set lattice_calc_on and plot_on
    if item.param.name == 'lattice_calc_on':
      lattice_calc_on = str(item.param.value)
      set_lattice_calc = True
    elif item.param.name == 'plot_on':
      plot_on = str(item.param.value)
      set_plot = True
    elif update_dict[item.param.name]:
      msg = pipe.cmd_in(set_str + item.param.name + " = " + str(item.param.value))
      if msg.find("ERROR") != -1:
        messagebox.showwarning(item.param.name,msg)
  #Now set lattice_calc_on and plot_on
  if set_plot:
    pipe.cmd_in("set global plot_on = " + plot_on)
  else:
    pipe.cmd_in("set global plot_on = True")
  if set_lattice_calc:
    pipe.cmd_in("set global plot_on = " + lattice_calc_on)
  else:
    pipe.cmd_in("set global lattice_calc_on = True")
  #Re-enable input for parameters that can vary
  for item in tao_list:
    if item.param.can_vary:
      item.tk_wid.configure(state="normal")

#---------------------------------------------------------------

if __name__ == "__main__":
  root = tao_root_window(sys.argv)
  root.mainloop()
