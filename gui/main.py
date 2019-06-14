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

    self.wm_geometry(newGeometry='400x600')

    self.outer_frame=tk.Frame(self)

    canvas=tk.Canvas(self.outer_frame)
    self.list_frame=tk.Frame(canvas)
    scrollbar=tk.Scrollbar(self.outer_frame,orient="vertical",command=canvas.yview)
    canvas.configure(yscrollcommand=scrollbar.set)

    def scrollhelper(event):
      canvas.configure(scrollregion=canvas.bbox("all")) #,width=200,height=200)
      new_width = event.width + 15
      old_geo = self.wm_geometry()
      if old_geo[2] != '1':
        new_geo = old_geo[old_geo.find('x'):old_geo.find('+')]
        new_geo = str(new_width) + new_geo
      else:
        new_geo = str(new_width) + 'x500'
      self.geometry(new_geo)
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
    #self.geometry('1550x600')
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

    b2.pack(side="right")
    b1.pack(side="right")

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
    self.bulk_params[0].tk_wid.bind("<Return>", lambda event : self.fill(0))
    self.bulk_filled.append(False)
    self.bulk_value.append(self.bulk_params[0].tk_var.get())
    self.bulk_apply.append(tk.Button(self.list_frame, text="Fill...", command= lambda : self.fill(0)))
    self.bulk_apply[0].grid(row=2, column=6)

    self.bulk_params.append(tk_tao_parameter(str_to_tao_param("good_user;LOGIC;T;"), self.list_frame))
    self.bulk_params[1].tk_wid.grid(row=1, column=11)
    self.bulk_params[1].tk_wid.bind("<Return>", lambda event : self.fill(1))
    self.bulk_filled.append(False)
    self.bulk_value.append(self.bulk_params[1].tk_var.get())
    self.bulk_apply.append(tk.Button(self.list_frame, text="Fill...", command= lambda : self.fill(1)))
    self.bulk_apply[1].grid(row=2, column=11)

    self.bulk_params.append(tk_tao_parameter(str_to_tao_param("weight;REAL;T;"), self.list_frame))
    self.bulk_params[2].tk_wid.grid(row=1, column=12)
    self.bulk_params[1].tk_wid.bind("<Return>", lambda event : self.fill(1))
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

  def fill(self, index, event=None):
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

class tao_var_general_window(tao_list_window):

  def __init__(self, root, pipe, *args, **kwargs):
    tao_list_window.__init__(self, root, "v1 Variables", *args, **kwargs)
    self.pipe = pipe
    self.list_frame.grid_columnconfigure(0, pad=10)
    self.list_frame.grid_columnconfigure(1, pad=10)
    self.list_frame.grid_columnconfigure(2, pad=10)
    self.refresh()

  def refresh(self):
    for child in self.list_frame.winfo_children():
      child.destroy()
    v1_var_list = self.pipe.cmd_in("python var_general")
    v1_var_list = v1_var_list.splitlines()
    for i in range(len(v1_var_list)):
      v1_var_list[i] = v1_var_list[i].split(';')

    tk.Label(self.list_frame, text="Variable").grid(row=0, column=0, columnspan=2)
    tk.Label(self.list_frame, text="Indices").grid(row=0, column=2)
    tk.Label(self.list_frame, text="Using").grid(row=0, column=3)

    i=1
    for item in v1_var_list:
      tk.Label(self.list_frame, text=item[0]).grid(row=i, column=0)
      tk.Button(self.list_frame, text="View...").grid(row=i, column=1)
      tk.Label(self.list_frame, text=item[1] + ':' + item[2]).grid(row=i, column=2)
      i = i+1

#---------------------------------------------------------------
# History Window

class tao_history_window(tao_list_window):

  def __init__(self, *args, **kwargs):
    tao_list_window.__init__(self, None, "History", *args, **kwargs)
    #self.geometry('400x600')
    self.refresh()

  def refresh(self):
    for child in self.list_frame.winfo_children():
      child.destroy()

    tk.Label(self.list_frame, text="Commands").grid(row=0, column=0)
    #tk.Label(self.list_frame, text="System shell").grid(row=0, column=1)
    tk.Label(self.list_frame, text="Command files").grid(row=0, column=1)

    for j in range(len(root.history)):
      ii = len(root.history[j]) #Actual row counter
      for i in range(len(root.history[j])):
        b = tk.Button(self.list_frame, text=root.history[j][i])
        b.configure(command=self.re_run_callback(root.history[j][i], j))
        b.bind("<Button-3>", self.re_run_callback(root.history[j][i], j+2))
        b.grid(row=ii, column=j)
        ii -= 1

  def re_run(self, cmd_string, mode, event=None):
    '''
    Re-runs the given command or command file (cmd_string),
    using the specified mode.
    mode == 0 -> Run in Tao/system shell
    mode == 1 -> Run in system shell
    mode == 2 -> Re-run command file
    Using modes 4, 5, and 6 simply respawns cmd_string in
    the command line or call_file box, and does not run it
    '''
    if mode ==0:
      root.command.tk_var.set(cmd_string)
      root.tao_command()
    elif mode == 1:
      root.call_file.tk_var.set(cmd_string)
      root.tao_call()
    elif mode ==2:
      root.command.tk_var.set(cmd_string)
    elif mode == 3:
      root.call_file.tk_var.set(cmd_string)
    #if mode == 0:
    #  root.command.tk_var.set(cmd_string)
    #  root.tao_command()
    #elif mode == 1:
    #  root.command.tk_var.set(cmd_string)
    #  root.tao_spawn()
    #elif mode == 2:
    #  root.call_file.tk_var.set(cmd_string)
    #  root.tao_call()
    #elif mode == 3:
    #  root.command.tk_var.set(cmd_string)
    #elif mode == 4:
    #  root.command.tk_var.set(cmd_string)
    #elif mode == 5:
    #  root.call_file.tk_var.set(cmd_string)

  def re_run_callback(self, cmd_string, mode, event=None):
    '''
    Formats a callback to self.re_run
    '''
    return lambda event=None : self.re_run(cmd_string, mode)


#---------------------------------------------------------------
# Root window

class tao_root_window(tk.Tk):

  def __init__(self, *args, **kwargs):
    tk.Tk.__init__(self)

    self.title("Tao")
    self.protocol("WM_DELETE_WINDOW", self.quit_cmd)
    self.tk.call('tk', 'scaling', 1.0)

    # Menu bar
    self.menubar_init()

    #self.menubar = tk.Menu(self)

    #file_menu = tk.Menu(self.menubar)
    #file_menu.add_command(label = 'Read...', command = self.read_cmd)
    #file_menu.add_command(label = 'Write...', command = self.write_cmd)
    #file_menu.add_command(label = 'Reinit...', command = self.reinit_cmd)
    #file_menu.add_separator()
    #file_menu.add_command(label = 'Quit', command = self.quit_cmd, accelerator = 'Ctrl+Q')
    #self.menubar.add_cascade(label = 'File', menu = file_menu)

    #window_menu = tk.Menu(self.menubar)
    #window_menu.add_command(label = 'Optimizer...', command = self.optimizer_cmd)
    #window_menu.add_command(label = 'Plotting...', command = self.plotting_cmd)
    #window_menu.add_command(label = 'Wave...', command = self.wave_cmd)
    #window_menu.add_command(label = 'Global Variables...', command = self.set_global_vars_cmd)
    #window_menu.add_command(label = 'Data...', command = self.view_data_cmd)
    #self.menubar.add_cascade(label = 'Window', menu = window_menu)

    #self.config(menu=self.menubar)

    # Init GUI

    init_frame = tk.Frame(self)
    init_frame.pack()
    self.tao_load(init_frame)

    # Key bindings

    self.bind_all("<Control-q>", self.quit_cmd)

  #-------------------------
  # Tao startup

  def start_main(self):
    self.unbind_all("<Return>")
    self.lift()
    self.focus_force()
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
    self.history.append([]) #Tao and shell history
    self.history_pos = 0 #Used for scrolling in history on command line
    #self.history.append([]) #Shell history
    self.history.append([]) #Call history
    tk.Button(self.cmd_frame, text="View History...", command=self.view_history_cmd).pack(side="top", fill="x")

    self.command = tk_tao_parameter(str_to_tao_param("command;STR;T;"), self.cmd_frame, self.pipe)
    self.command.tk_wid.bind("<Return>", self.tao_command)
    self.command.tk_wid.bind("<Shift-Return>", self.tao_spawn)
    self.command.tk_wid.bind("<Up>", self.hist_scroll_up)
    self.command.tk_wid.bind("<Down>", self.hist_scroll_down)
    self.command.tk_wid.pack(side="left", fill="x", expand=1)
    tk.Button(self.cmd_frame, text="Run (System Shell)", command=self.tao_spawn).pack(side="right")
    tk.Button(self.cmd_frame, text="Run (Tao)", command=self.tao_command).pack(side="right")

  def hist_scroll_up(self, event=None):
    if len(self.history[0]) > self.history_pos: #if there's room left to scroll up in history
      self.history_pos += 1
      self.command.tk_var.set(self.history[0][-1*self.history_pos])
      self.command.tk_wid.icursor(len(self.command.tk_var.get()))

  def hist_scroll_down(self, event=None):
    if self.history_pos > 0:
      self.history_pos -= 1
      if self.history_pos > 0:
        self.command.tk_var.set(self.history[0][-1*self.history_pos])
      else:
        self.command.tk_var.set("")
      self.command.tk_wid.icursor(len(self.command.tk_var.get()))

  def tao_command(self, event=None):
    '''
    Runs the text in self.command at the Tao command line, appends it to the history, and clears self.command
    '''
    if self.command.tk_var.get() != "":
      self.pipe.cmd(self.command.tk_var.get())
      self.history[0].append(self.command.tk_var.get())
      self.command.tk_var.set("")
    self.history_pos = 0
    #Try to refresh history window
    try:
      self.history_window.refresh()
    except:
      pass

  def tao_spawn(self, event=None):
    '''
    Runs the text in self.command at the system command line, appends it to the history, and clears self.command
    '''
    if self.command.tk_var.get() != "":
      cmd_txt = self.command.tk_var.get()
      self.command.tk_var.set("spawn " + cmd_txt)
      self.tao_command()
    #if self.command.tk_var.get() != "":
    #  cmd_txt = self.command.tk_var.get()
    #  self.pipe.cmd("spawn " + cmd_txt)
    #  self.history[1].append(cmd_txt)
    #  self.command.tk_var.set("")
    #self.history_pos = 0
    ##Try to refresh history window
    #try:
    #  self.history_window.refresh()
    #except:
    #  pass

  def tao_call(self):
    '''
    Runs the command file in self.call_file, appends it to the history, and clears self.call_file
    '''
    if self.call_file.tk_var.get() != "Browse...":
      self.pipe.cmd("call " + self.call_file.tk_var.get())
      self.history[1].append(self.call_file.tk_var.get())
      self.call_file.tk_var.set("Browse...")
    #Try to refresh history window
    try:
      self.history_window.refresh()
    except:
      pass

  def menubar_init(self):
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
    window_menu.add_command(label = 'Variables...', command = self.view_vars_cmd)
    window_menu.add_command(label = 'Global Variables...', command = self.set_global_vars_cmd)
    window_menu.add_command(label = 'Data...', command = self.view_data_cmd)
    self.menubar.add_cascade(label = 'Window', menu = window_menu)

    self.config(menu=self.menubar)

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
      init_file.close()
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
      c1 = (name in ["tao_exe", "tao_lib"])
      try:
        c2 = (param_dict[name].type == 'FILE')
      except KeyError:
        c2 = False
      if c1 | c2:
        #Check that a good filename has been given
        try:
          filename = init_dict[name]
          #Expand environment variables and ~
          filename = os.path.expanduser(filename)
          filename = os.path.expandvars(filename)
          f = open(filename)
          f.close()
          init_dict[name] = filename
        except:
          messagebox.showwarning(tk_list[k].param.name, "File not found: " + filename)
          init_dict.pop(name)
          #Remove bad files from init dict
    k = 0 #row number counter
    for param, tao_param in param_dict.items():
      tk_list.append(tk_tao_parameter(tao_param,init_frame))
      #Possibly set value from init file
      if tk_list[k].param.name in init_dict:
        if tk_list[k].param.type == 'FILE':
          tk_list[k].tk_var.set(init_dict[tk_list[k].param.name])
        elif tk_list[k].param.type == 'LOGIC':
          state = (init_dict[tk_list[k].param.name] == 'T') | (init_dict[tk_list[k].param.name] == 'True')
          tk_list[k].tk_var.set(state)
      tk.Label(init_frame,text=param).grid(row=k,sticky="E")
      tk_list[k].tk_wid.grid(row=k, column=1, sticky="W")
      k = k+1
    # Choosing whether to use pexpect or ctypes must be handled separately
    tk.Label(init_frame, text="Choose interface").grid(row=k, sticky="E")
    def swap_box(event=None):
      if chosen_interface.get() == "pexpect":
        ctype_label.grid_forget()
        tao_lib.tk_wid.grid_forget()
        pexp_label.grid(row=k+1, sticky="E")
        tao_exe.tk_wid.grid(row=k+1, column=1, sticky="W")
      elif chosen_interface.get() == "ctypes":
        pexp_label.grid_forget()
        tao_exe.tk_wid.grid_forget()
        ctype_label.grid(row=k+1, sticky="E")
        tao_lib.tk_wid.grid(row=k+1, column=1, sticky="W")

    chosen_interface = tk.StringVar()
    chosen_interface.set("pexpect")
    if "chosen_interface" in init_dict:
      if init_dict["chosen_interface"] in ["pexpect", "ctypes"]:
        chosen_interface.set(init_dict["chosen_interface"])
    tk.OptionMenu(init_frame, chosen_interface, "pexpect", "ctypes", command=swap_box).grid(row=k, column=1, sticky="W")
    pexp_label = tk.Label(init_frame, text="Tao executable")
    ctype_label = tk.Label(init_frame, text="Shared Library")
    tao_exe = tk_tao_parameter(str_to_tao_param("tao_exe;FILE;T;"), init_frame)
    if "tao_exe" in init_dict:
      tao_exe.tk_var.set(init_dict["tao_exe"])
    tao_lib = tk_tao_parameter(str_to_tao_param("tao_lib;FILE;T;"), init_frame)
    if "tao_lib" in init_dict:
      tao_exe.tk_var.set(init_dict["tao_exe"])
    swap_box()

    def param_load(event=None):
      if chosen_interface.get() == "ctypes":
        messagebox.showwarning("Error", "ctypes is not currently supported.  Please use pexpect.")
        return 0
      init_args = ""
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
            init_args = init_args + "-" + tk_param.param.name + " " + tk_param.param.value + " "
        elif tk_param.param.type == 'LOGIC':
          tk_param.param.value = tk_param.tk_var.get()
          if tk_param.param.value:
            init_args = init_args + "-" + tk_param.param.name + " "
      # Run Tao, clear the init_frame, and draw the main frame
      from tao_interface import tao_interface
      if chosen_interface.get() == "pexpect":
        mode = "pexpect"
        if tao_exe.tk_var.get() == "Browse...":
          tao_exe.tk_var.set("")
        self.pipe = tao_interface(mode, init_args, tao_exe.tk_var.get())
      elif chosen_interface.get() == "ctypes":
        if tao_lib.tk_var.get() == "Browse...":
          tao_lib.tk_var.set("")
        mode = "ctypes"
        self.pipe = tao_interface(mode, init_args, tao_lib.tk_var.get())

      init_frame.destroy()
      self.start_main()

    load_b = tk.Button(init_frame, text="Start Tao", command=param_load)
    load_b.grid(row=k+2, columnspan=2)
    self.bind_all("<Return>", param_load)

    #Start Tao immediately if skip_setup is set true
    if "skip_setup" in init_dict:
      if init_dict["skip_setup"] in ['T', 'True']:
        param_load()

  #-------------------------
  # Menu bar callbacks

  def read_cmd(self):
    print ('Read called')

  def write_cmd(self):
    print ('Write called')

  def reinit_cmd(self):
    '''
    Quit Tao, destroy the main frame, and respawn the init frame
    '''
    result = messagebox.askquestion("Reinit", "This will close all windows and restart Tao.  Are you sure?", icon="warning")
    if result != 'yes':
      return

    for child in root.winfo_children():
      child.destroy()

    self.pipe.cmd("quit")
    self.menubar_init()
    init_frame = tk.Frame(self)
    init_frame.pack()
    self.tao_load(init_frame)

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

  def view_vars_cmd(self):
    win = tao_var_general_window(None, self.pipe)

  def view_data_cmd(self):
    win = tao_d2_data_window(None, self.pipe)

  # Other callbacks

  def global_vars_event(self, event):
    self.set_global_vars_cmd()

  def view_history_cmd(self):
    #Just focus the existing history window, if it exists
    try:
      self.history_window.lift()
      #self.history_window.force_focus()
    except:
      self.history_window = tao_history_window()

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
    pass
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
