import tkinter as tk
import ttk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import font
import sys
import os
import copy
sys.path.append(os.environ['ACC_ROOT_DIR'] + '/tao/gui')
from tao_widget import *
from taoplot import taoplot
from parameters import str_to_tao_param
from elements import *
from main import tao_set
import string


#-----------------------------------------------------
# List window

class tao_list_window(tk.Toplevel):

  def __init__(self, root, title, use_upper=False, min_width=0, *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
    self.title(title)
    self.root = root

    #self.wm_geometry(newGeometry='400x600')

    self.upper_frame=tk.Frame(self)
    self.outer_frame=tk.Frame(self)

    canvas=tk.Canvas(self.outer_frame)
    self.list_frame=tk.Frame(canvas)
    scrollbar=tk.Scrollbar(self.outer_frame,orient="vertical",command=canvas.yview)
    canvas.configure(yscrollcommand=scrollbar.set)

    def scrollhelper(event):
      canvas.configure(scrollregion=canvas.bbox("all")) #,width=200,height=200)
      new_width = event.width + 15
      # Don't resize to a smaller size
      old_geo = self.wm_geometry()
      old_width = int(old_geo[:old_geo.find('x')])
      if old_width > new_width:
        new_width = old_width
      # Don't resize to be smaller than the
      # upper frame
      self.upper_frame.update()
      upper_width = self.upper_frame.winfo_width()
      if upper_width > new_width:
        new_width = upper_width
      # Don't resize to less than the min width
      if min_width > new_width:
        new_width = min_width
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

    if use_upper:
      self.upper_frame.pack(side="top",fill="both", expand=1)
    self.outer_frame.pack(side="top",fill="both",expand=1)
    scrollbar.pack(side="right",fill="y")
    canvas.pack(side="left",fill="both",expand=1)
    canvas.create_window((0,0),window=self.list_frame,anchor='nw')

#-----------------------------------------------------
# Parameter frame
class tao_parameter_frame(tk.Frame):
  '''
  Meant to display a list of parameters in a given
  number of columns
  tao_list should be a list of tao_parameters
  '''

  def __init__(self, parent, tao_list, n_col, pipe, *args, **kwargs):
    tk.Frame.__init__(self, parent, *args, **kwargs)
    self.tao_list = [] #List for tk_tao_parameters
    for p in tao_list:
      self.tao_list.append(tk_tao_parameter(p, self, pipe))
    cols = [] # A list of the parts of tao_list,
    # after being divided into n_col columns
    len_col = int(len(self.tao_list)/n_col) # How many elements per column
    if len_col < len(self.tao_list)/n_col:
      #if the result was rounded down
      len_col = len_col+1
    for i in range(n_col):
      cols.append(self.tao_list[i*len_col : (i+1)*len_col])

    # Grid the widgets in each column
    j = 0 #controls which column the widgets get gridded to
    for c in cols:
      for i in range(len(c)):
        c[i].tk_label.grid(row=i, column=2*j, sticky = 'E')
        c[i].tk_wid.grid(row=i, column=2*j+1, sticky = 'EW')
      j = j+1

#-----------------------------------------------------
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
      self.tao_list[k].tk_wid.grid(row=k,column=1,sticky="EW")
      k = k+1

#-----------------------------------------------------
# Table window
class table_window(tao_list_window):
  '''
  Meant for showing large amounts of information in a table (e.g. d1_data and v1_vars).  Comes with bulk editing, detailed view of individual rows, and editing of individual parameters in the table.

  Input parameters:
  root: parent widget
  pipe: the tao_interface object that allows interface with Tao
  array_name: the name of the object this table displays (e.g. the name of the d1_datum or v1_variable)
  title_list: the column titles, in order from left to right
  bulk_template: a list with elements [tao_parameter, column], where tao_parameter is a generic copy of the parameter that will be bulk filled (i.e. initialized to "blank"), and column is the column number it should be gridded to (start counting from 0 on the left)
  bulk_set_format: format string for the bulk set string, to be used with the str.format() method
  set_format: format string for individual rows' set strings, to be used with the str.format() method
    example format strings: "set data 1@{}|", "set data 1@{}[{}]|"
  '''

  def __init__(self, root, pipe, array_name, title_list, bulk_template, bulk_set_format, set_format, *args, **kwargs):
    tao_list_window.__init__(self, root, array_name, *args, **kwargs)
    self.pipe = pipe
    #Make the font a bit smaller
    #default_font = font.nametofont("TkDefaultFont")
    #default_font.configure(size=14)
    #self.option_add("*Font", default_font)
    self.array_name = array_name
    self.title_list = title_list
    self.bulk_template = bulk_template
    self.bulk_set_format = bulk_set_format
    self.set_format = set_format
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

    # Grid the column titles
    j = 0
    for item in self.title_list:
      tk.Label(self.list_frame, text=item).grid(row=0, column=j)
      self.list_frame.grid_columnconfigure(j, pad=10)
      j=j+1

    # Bulk editing
    tk.Label(self.list_frame, text="Bulk editing:").grid(row=1, column=0, columnspan=self.bulk_template[0][1])
    tk.Label(self.list_frame, text="Click to fill:").grid(row=2, column=0, columnspan=self.bulk_template[0][1])
    #self.bulk_template[0][1] is the appropriate
    #columnspan because these labels can fill all
    #the space until the first bulk fill box

    self.bulk_params = [] #Holds the current bulk edit widgets
    self.bulk_filled = [] #Holds whether or not bulk filling has been used
    self.bulk_value = [] #Holds the last value that was filled to the cells
    self.bulk_apply = [] #Holds fill buttons

    j=0
    for item in self.bulk_template:
      self.bulk_params.append(tk_tao_parameter(copy.copy(item[0]), self.list_frame, self.pipe))
      self.bulk_params[j].tk_wid.grid(row=1, column=item[1])
      self.bulk_params[j].tk_wid.bind("<Return>", self.fill_callback(j))
      self.bulk_filled.append(False)
      self.bulk_value.append(self.bulk_params[j].tk_var.get())
      self.bulk_apply.append(tk.Button(self.list_frame, text="Fill...", command = self.fill_callback(j) ))
      self.bulk_apply[j].grid(row=2, column=item[1])
      j=j+1

    #Fetch and fill in the data
    # IT IS EXPECTED THAT SUBCLASSES WILL DEFINE
    # self.row_num AND A METHOD FOR POPULATING
    # self.list_rows, AND THEN CALL THIS METHOD
    #i = row counter, j = column counter
    #grid to row i+3 because row 0 is titles, row 1 is bulk editing widgets, row 2 is fill buttons
    self.list_rows = []
    for i in range(self.row_num):
      self.list_rows.append(self.list_row_fetch(i))
      for j in range(len(self.list_rows[i].tk_wids)):
        self.list_rows[i].tk_wids[j].grid(row=i+3, column=j)
      tk.Button(self.list_frame, text="View More...", command=self.open_detail_window_callback(self.list_rows[i].index)).grid(row=i+3, column=len(self.list_rows[i].tk_wids))

  def fill(self, index, event=None):
    '''
    Fills the column specified by bulk_template[index][1] to the bulk edit value, saves the bulk edit value for efficient calls to tao_set, and clears the bulk edit box
    '''
    # Save the bulk parameter state
    self.bulk_value[index] = self.bulk_params[index].tk_var.get()
    self.bulk_filled[index] = True
    # Clear the bulk parameter widget
    if self.bulk_params[index].param.type in ['STR', 'INT', 'REAL']:
      self.bulk_params[index].tk_var.set("")
    elif self.bulk_params[index].param.type == 'LOGIC':
      self.bulk_params[index].tk_var.set(False)
    # Fill the appropriate variable
    for i in range(len(self.list_rows)):
      self.list_rows[i].tk_tao_params[self.bulk_template[index][0].name].tk_var.set(self.bulk_value[index])

  def fill_callback(self, index, event=None):
    return lambda event=None : self.fill(index)


  def apply(self):
    #Apply bulk changes
    for i in range(len(self.bulk_params)):
      if self.bulk_filled[i]:
        set_str = self.bulk_set_format.format(self.array_name)
        self.bulk_params[i].tk_var.set(self.bulk_value[i])
        tao_set([self.bulk_params[i]], set_str, self.pipe, overide=(self.bulk_params[i].param.type=='LOGIC')) #overide is necessary for LOGIC parameters

    #Apply individual changes that are different from bulk changes
    for i in range(len(self.list_rows)):
      set_list = []
      #NOTE: it is expected that self.list_rows[i].index exist, and be equal to that row's index
      set_str = self.set_format.format(self.array_name, self.list_rows[i].index)

      #Find elements in row that need setting
      for j in range(len(self.bulk_template)):
        name = self.bulk_template[j][0].name
        c1 = (self.list_rows[i].tk_tao_params[name].tk_var.get() != self.bulk_value[j])
        c2 = not self.bulk_filled[j]
        try:
          if self.bulk_template[j][0].type == 'REAL':
            c3 = (float(self.list_rows[i].tk_tao_params[name].tk_var.get()) != self.list_rows[i].tk_tao_params[name].param.value)
          elif self.bulk_template[j][0].type == 'INT':
            c3 = (int(self.list_rows[i].tk_tao_params[name].tk_var.get()) != self.list_rows[i].tk_tao_params[name].param.value)
          elif self.bulk_template[j][0].type == 'STR':
            c3 = (str(self.list_rows[i].tk_tao_params[name].tk_var.get()) != self.list_rows[i].tk_tao_params[name].param.value)
          elif self.bulk_template[j][0].type == 'LOGIC':
            c3 = (bool(self.list_rows[i].tk_tao_params[name].tk_var.get()) != self.list_rows[i].tk_tao_params[name].param.value)
        except ValueError:
          c3 = False
        if (c1 | c2) & c3:
          set_list.append(self.list_rows[i].tk_tao_params[name])

      if set_list != []:
        tao_set(set_list, set_str, self.pipe)

    #Refresh
    self.refresh()

  def open_detail_window_callback(self, index):
    return lambda : self.open_detail_window(index)

  def open_detail_window(self, index):
    # SUBCLASSES ARE EXPECTED TO DEFINE self.param_list
    # FOR USE IN A tao_parameter_window, THEN CALL
    # THIS METHOD
    detail_title = self.array_name + '[' + str(index) + ']'
    win = tao_parameter_window(self, detail_title, self.param_list, self.pipe)

    set_str = self.set_format.format(self.array_name, index)
    b = tk.Button(win.button_frame, text="Apply changes", command=lambda : self.detail_set_callback(win.tao_list,set_str))
    b.pack()

  def detail_set_callback(self, tao_list, set_str):
    tao_set(tao_list, set_str, self.pipe)
    self.refresh()


#-----------------------------------------------------
# d2_data window

class tao_d2_data_window(tao_list_window):

  def __init__(self, root, pipe, *args, **kwargs):
    tao_list_window.__init__(self, root, "Data", *args, **kwargs)
    self.pipe = pipe
    univ_list = self.pipe.cmd_in("python super_universe")
    n_universe = str_to_tao_param(univ_list.splitlines()[0])
    self.univ_frame = tk.Frame(self)
    tk.Label(self.univ_frame, text="Universe: ").grid(row=0,column=0,sticky="E")
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
      new_frame = d2_data_frame(self.list_frame, self.root, self.pipe, d2_data_item, u_ix)
      new_frame.frame.pack()


#-----------------------------------------------------
# d1_data window

class tao_d1_data_window(table_window):

  def __init__(self, root, pipe, d1_data_name, u_ix, ix_lb, ix_ub, *args, **kwargs):
    self.u_ix = u_ix
    self.ix_lb = ix_lb
    self.ix_ub = ix_ub
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
    bulk_template = []
    bulk_template.append([str_to_tao_param("meas_value;REAL;T;"), 6])
    bulk_template.append([str_to_tao_param("good_user;LOGIC;T;"), 11])
    bulk_template.append([str_to_tao_param("weight;REAL;T;"), 12])

    bulk_set_format = "set data " + str(self.u_ix) + '@{}|'
    set_format = "set data " + str(self.u_ix) + '@{}[{}]|'
    table_window.__init__(self, root, pipe, d1_data_name, title_list, bulk_template, bulk_set_format, set_format, *args, **kwargs)

  def refresh(self):
    self.d_list = self.pipe.cmd_in("python data_d_array " + self.u_ix + '@' + self.array_name)
    self.d_list = self.d_list.splitlines()
    self.row_num = len(self.d_list)
    table_window.refresh(self)

  def list_row_fetch(self, index):
    '''
    Returns row number index to be appended to self.list_rows
    '''
    return d1_data_list_entry(self.list_frame, self.d_list[index])

  def open_detail_window(self, index):
    self.param_list = self.pipe.cmd_in("python data " + str(self.u_ix) + '@' + self.array_name + '[' + str(index) + ']')
    self.param_list = self.param_list.splitlines()
    for i in range(len(self.param_list)):
      self.param_list[i]=str_to_tao_param(self.param_list[i])
    table_window.open_detail_window(self, index)

#-----------------------------------------------------
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
      tk.Button(self.list_frame, text="View...", command=self.open_v1_callback(item[0])).grid(row=i, column=1)
      tk.Label(self.list_frame, text=item[2] + ':' + item[3]).grid(row=i, column=2)
      tk.Label(self.list_frame, text=item[1]).grid(row=i, column=3)
      i = i+1

  def open_v1_callback(self, v1_var_name):
    return lambda : self.open_v1(v1_var_name)

  def open_v1(self, v1_var_name):
    win = tao_v1_var_window(self.root, self.pipe, v1_var_name)

#-----------------------------------------------------
# v1_var window

class tao_v1_var_window(table_window):

  def __init__(self, root, pipe, v1_var_name, *args, **kwargs):
    title_list = ["Index",
        "Name",
        "Meas value",
        "Model value",
        "Design value",
        "Useit_opt",
        "good_user",
        "Weight"]
    bulk_template = []
    bulk_template.append([str_to_tao_param("meas_value;REAL;T;"), 2])
    bulk_template.append([str_to_tao_param("good_user;LOGIC;T;"), 6])
    bulk_template.append([str_to_tao_param("weight;REAL;T;"), 7])

    bulk_set_format = "set var {}|"
    set_format = "set var {}[{}]|"
    table_window.__init__(self, root, pipe, v1_var_name, title_list, bulk_template, bulk_set_format, set_format, *args, **kwargs)

  def refresh(self):
    self.v_list = self.pipe.cmd_in("python var_v_array " + self.array_name)
    self.v_list = self.v_list.splitlines()
    self.row_num = len(self.v_list)
    table_window.refresh(self)

  def list_row_fetch(self, index):
    '''
    Returns row number index to be appended to self.list_rows
    '''
    return v1_var_list_entry(self.list_frame, self.v_list[index])

  def open_detail_window(self, index):
    self.param_list = self.pipe.cmd_in("python var " + self.array_name + '[' + str(index) + ']')
    self.param_list = self.param_list.splitlines()
    for i in range(len(self.param_list)):
      self.param_list[i]=str_to_tao_param(self.param_list[i])
    table_window.open_detail_window(self, index)


#-----------------------------------------------------
# History Window

class tao_history_window(tao_list_window):

  def __init__(self, root, *args, **kwargs):
    tao_list_window.__init__(self, root, "History", *args, **kwargs)
    self.refresh()

  def refresh(self):
    for child in self.list_frame.winfo_children():
      child.destroy()

    tk.Label(self.list_frame, text="Commands").grid(row=0, column=0)
    tk.Label(self.list_frame, text="Command files").grid(row=0, column=1)

    for j in range(len(self.root.history)):
      ii = len(self.root.history[j]) #Actual row counter
      for i in range(len(self.root.history[j])):
        b = tk.Button(self.list_frame, text=self.root.history[j][i])
        b.configure(command=self.re_run_callback(self.root.history[j][i], j))
        b.bind("<Button-3>", self.re_run_callback(self.root.history[j][i], j+2))
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
      self.root.command.tk_var.set(cmd_string)
      self.root.tao_command()
    elif mode == 1:
      self.root.call_file.tk_var.set(cmd_string)
      self.root.tao_call()
    elif mode ==2:
      self.root.command.tk_var.set(cmd_string)
    elif mode == 3:
      self.root.call_file.tk_var.set(cmd_string)

  def re_run_callback(self, cmd_string, mode, event=None):
    '''
    Formats a callback to self.re_run
    '''
    return lambda event=None : self.re_run(cmd_string, mode)
#-----------------------------------------------------
# Plot template window
#TODO: source input validation (base model design etc)

class tao_plot_tr_window(tao_list_window):
  '''
  Displays information about existing plot templates and regions
  Use a drop-down list to select template to view.
  '''
  def __init__(self, root, pipe, mode, *args, **kwargs):
    tao_list_window.__init__(self, root, "Plot Templates", *args, **kwargs)
    self.minsize(325, 100)
    self.pipe = pipe
    self.mode = mode
    if self.mode == "R":
      self.title("Plot Regions")

    self.temp_frame = tk.Frame(self)
    if self.mode == "T":
      tk.Label(self.temp_frame, text="Choose template:").grid(row=0, column=0)
    elif self.mode == "R":
      tk.Label(self.temp_frame, text="Choose plot:").grid(row=0, column=0)
    self.temp_frame.columnconfigure(0, pad=10)
    self.temp_frame.pack(fill="both", expand=0)

    if self.mode == "T":
      plot_list = self.pipe.cmd_in("python plot_list t")
    elif self.mode == "R":
      plot_list = self.pipe.cmd_in("python plot_list r")
    # Populate self.plot_list and self.index_list
    if self.mode == "T":
      self.plot_list = plot_list.splitlines()
      self.index_list = len(self.plot_list)*[0] #get correct length
      for i in range(len(self.plot_list)):
        self.index_list[i], self.plot_list[i] = self.plot_list[i].split(';')
    elif self.mode == "R":
      #Only use regions that contain a plot
      self.plot_list = []
      self.index_list = []
      self.region_list = [] #needed to open the correct graph/curve windows
      plot_list = plot_list.splitlines()
      for i in range(len(plot_list)):
        if plot_list[i].split(';')[2] != "":
          self.plot_list.append(plot_list[i].split(';')[2])
          self.index_list.append(plot_list[i].split(';')[0])
          self.region_list.append(plot_list[i].split(';')[1])

    if self.plot_list == []:
      tk.Label(self.list_frame, text="NO PLOTS FOUND").pack()
      return

    self.temp = tk.StringVar()
    self.temp.set(self.plot_list[0])
    self.temp_select = ttk.Combobox(self.temp_frame, textvariable=self.temp, values=self.plot_list)
    self.temp_select.bind("<<ComboboxSelected>>", self.refresh)
    self.temp_select.bind("<Return>", self.refresh)
    self.temp_select.grid(row=0, column=1)

    self.refresh()

    self.button_frame = tk.Frame(self)
    self.button_frame.pack(fill="both", expand=0)
    b = tk.Button(self.button_frame, text="Apply changes", command=self.plot_apply)
    b.pack(side="left", fill="both", expand=0)
    if (self.root.plot_mode == "matplotlib") & (self.mode == "T"):
      plot_b = tk.Button(self.button_frame, text="Plot!", command=self.mpl_plot)
      plot_b.pack(side="right", fill="both", expand=0)

  def mpl_plot(self, event=None):
    '''
    Creates a new matplotlib window with the currently selected plot
    '''
    # Check if the template has been placed in a region already, and place it if necessary
    if self.plot in self.root.placed.keys():
      pass
    else:
      #Place the plot in the next available region
      #and make visible
      r_index = 1
      while (("r" + str(r_index)) in self.root.placed.values()):
        r_index = r_index + 1
      self.pipe.cmd_in("place r" + str(r_index) + " " + self.plot)
      msg = self.pipe.cmd_in("set plot r" + str(r_index) + ' visible = T')
      self.root.placed[self.plot] = 'r' + str(r_index)
    # Plot with matplotlib
    x = taoplot(self.pipe, self.root.placed[self.plot])
    x.plot()


  def refresh(self, event=None):
    '''
    Clears self.list_frame and populates it with information relevant to self.temp
    '''
    # Don't run if self.temp is not a valid template
    if self.temp.get() not in self.plot_list:
      return

    # Clear list_frame
    for child in self.list_frame.winfo_children():
      child.destroy()

    # Get info on self.temp
    self.plot = self.temp.get()
    data_list = self.pipe.cmd_in("python plot1 " + self.plot)
    data_list = data_list.splitlines()
    num_graphs = data_list.pop(0)
    num_graphs = int(num_graphs.split(';')[3])
    self.graph_list = []
    for i in range(num_graphs):
      self.graph_list.append(data_list.pop(0))
      self.graph_list[i] = self.graph_list[i].split(';')[3]

    # Grid the name and the graphs
    name = tk_tao_parameter(str_to_tao_param(data_list.pop(0)), self.list_frame, self.pipe)
    name.tk_label.grid(row=0, column=0)
    name.tk_wid.grid(row=0, column=1, sticky='EW')
    tk.Label(self.list_frame, text="Graphs").grid(row=1, column=0, rowspan=num_graphs)
    i=1
    for graph in self.graph_list:
      tk.Button(self.list_frame, text=graph, command=self.open_graph_callback(name.param.value, graph)).grid(row=i, column=1, sticky='EW')
      i = i+1

    # Grid the rest of the information
    self.list_frame.columnconfigure(1, weight=1)
    for i in range(len(data_list)):
      data_list[i] = tk_tao_parameter(str_to_tao_param(data_list[i]), self.list_frame, self.pipe)
      # Make sure the Entry boxes are big enough to
      # display their contents
      if data_list[i].param.type in ['STR', 'INT', 'REAL']:
        data_list[i].tk_wid.configure(width=len(str(data_list[i].param.value))+1)
      # Grid to row i+1+num_graphs because row 0
      # is for the name and there are num_graphs rows
      # of graph buttons
      data_list[i].tk_label.grid(row=i+1+num_graphs, column=0)
      data_list[i].tk_wid.grid(row=i+1+num_graphs, column=1, sticky='EW')
    self.tao_list = data_list

  def open_graph_callback(self, plot, graph):
    return lambda : self.open_graph(plot, graph)

  def open_graph(self, plot, graph):
    '''
    Opens a window to display information about plot.graph
    '''
    if self.mode == "T":
      graph = plot+'.'+graph
    elif self.mode == "R":
      region = self.region_list[self.plot_list.index(self.plot)]
      graph = region+'.'+graph
    index = self.index_list[self.plot_list.index(self.plot)]
    win = tao_plot_graph_window(self.root, graph, self.pipe, self.mode, index)

  def plot_apply(self):
    index = self.index_list[self.plot_list.index(self.plot)]
    set_str = "set plot @" + self.mode + str(index) + ' '
    tao_set(self.tao_list, set_str, self.pipe)

#-----------------------------------------------------
# plot_graph window
class tao_plot_graph_window(tao_list_window):
  '''
  Displays information about a given graph
  mode: passed from the plot template/region window,
  should be "T" or "R"
  index: passed from the plot template/region window,
  should be the index of this graph's plot template/region
  '''
  def __init__(self, root, graph, pipe, mode, index, *args, **kwargs):
    tao_list_window.__init__(self, root, graph, *args, **kwargs)
    self.pipe = pipe
    self.graph = graph
    self.mode = mode
    self.index = index

    self.refresh()

    self.button_frame = tk.Frame(self)
    self.button_frame.pack(fill="both", expand=0)
    b = tk.Button(self.button_frame, text="Apply changes", command=self.graph_apply)
    b.pack()

  def refresh(self):
    # Clear self.list_frame
    for child in self.list_frame.winfo_children():
      child.destroy

    # Fetch data to display
    data_list = self.pipe.cmd_in("python plot_graph " + self.graph)
    data_list = data_list.splitlines()
    num_curves = data_list.pop(0)
    num_curves = int(num_curves.split(';')[3])
    curve_list = []
    for i in range(num_curves):
      curve_list.append(data_list.pop(0))
      curve_list[i] = curve_list[i].split(';')[3]

    # Display the name
    name = tk_tao_parameter(str_to_tao_param(data_list.pop(0)), self.list_frame, self.pipe)
    name.tk_label.grid(row=0, column=0)
    name.tk_wid.grid(row=0, column=1, sticky='EW')

    # Curve buttons
    if num_curves > 0:
      tk.Label(self.list_frame, text="Curves").grid(row=1, column=0, rowspan=num_curves)
      i=1
      for curve in curve_list:
        tk.Button(self.list_frame, text=curve, command=self.open_curve_callback(self.graph, curve)).grid(row=i, column=1, sticky='EW')
        i = i+1

    # Grid everything else
    self.list_frame.columnconfigure(1, weight=1)
    for i in range(len(data_list)):
      data_list[i] = tk_tao_parameter(str_to_tao_param(data_list[i]), self.list_frame, self.pipe)
    self.tao_list = data_list
    for i in range(len(data_list)):
      # Make sure the Entry boxes are big enough to
      # display their contents
      if data_list[i].param.type in ['STR', 'INT', 'REAL']:
        data_list[i].tk_wid.configure(width=len(str(data_list[i].param.value))+1)
      data_list[i].tk_label.grid(row=i+1+num_curves, column=0)
      data_list[i].tk_wid.grid(row=i+1+num_curves, column=1, sticky='EW')

  def open_curve_callback(self, graph, curve):
    return lambda : self.open_curve(graph, curve)

  def open_curve(self, graph, curve):
    '''
    Opens a window to display info for a given curve
    '''
    curve = graph + '.' + curve
    win = tao_plot_curve_window(self.root, curve, self.pipe)

    b = tk.Button(win.button_frame, text="Apply changes", command=lambda : self.curve_apply(win))
    b.pack()

  def curve_apply(self, win):
    #Convert from plot.graph.curve to index.graph.curve
    curve_name = str(self.index) + win.curve[win.curve.index('.'):]
    set_str = "set curve @" + self.mode + curve_name + " "
    tao_set(win.tao_list, set_str, win.pipe)

  def graph_apply(self):
    #Convert from plot.graph to index.graph
    graph_name = str(self.index) + self.graph[self.graph.index('.'):]
    set_str = "set graph @" + self.mode + graph_name + ' '
    tao_set(self.tao_list, set_str, self.pipe)

#-----------------------------------------------------
# plot_curve window

class tao_plot_curve_window(tao_parameter_window):
  '''
  Displays info for a given curve
  '''
  def __init__(self, root, curve, pipe, *args, **kwargs):
    # Get the parameters
    self.curve = curve
    self.pipe = pipe
    data_list = self.pipe.cmd_in("python plot_curve " + curve)
    data_list = data_list.splitlines()
    for i in range(len(data_list)):
      data_list[i] = str_to_tao_param(data_list[i])

    tao_parameter_window.__init__(self, root, curve, data_list, self.pipe, *args, **kwargs)

#-----------------------------------------------------
# Element window

class tao_ele_window(tao_list_window):
  '''
  Window for displaying and modifying element info.
  default specifies the element that should be
  displayed when the window is created
  Format for default: [universe, branch, element, base/model/design]
  '''
  def __init__(self, root, pipe, default=None, *args, **kwargs):
    tao_list_window.__init__(self, root, "Lattice Elements", use_upper=True, min_width=600, *args, **kwargs)
    self.pipe = pipe
    self.default = default

    # Get list of universes, branches, and
    # number of elements in each branch
    n_uni = self.pipe.cmd_in("python super_universe")
    n_uni = n_uni.splitlines()
    n_uni = int(n_uni[0].split(';')[3])
    self.u_list = [] # List of universes
    self.b_list = {} # b_list[i] is a list of branches in universe i
    self.b_name_list = {} # b_name_list[i] is a list of branch names in universe i
    self.e_list = {} #e_list[i][j] is a range from 0 to the max element number
    self.e_name_list = {} #e_name_list[i][j] is a list of ele names in branch j of uni i
    self.e_display_list = {} # indexed the same as e_list and e_name_list, but for displaying the number and name
    # There are three separate lists for the elements because
    # users can specify elements by name, by number,
    # or by selecting from a dropdown, which displays
    # the name and number
    for i in range(n_uni):
      self.u_list.append(str(i+1))

    for u in self.u_list:
      self.b_list[u] = self.pipe.cmd_in("python lat_general " + str(u))
      self.b_list[u] = self.b_list[u].splitlines()
      self.b_name_list[u] = []
      self.e_list[u] = {}
      self.e_name_list[u] = {}
      self.e_display_list[u] = {}
      for i in range(len(self.b_list[u])):
        branch_num = self.b_list[u][i].split(';')[0]
        branch_name = self.b_list[u][i].split(';')[1]
        ele_num = int(self.b_list[u][i].split(';')[3])
        self.b_list[u][i] = branch_num
        self.b_name_list[u].append('(' + branch_num + ') ' + branch_name)
        self.e_list[u][branch_num] = range(ele_num+1)
        self.e_display_list[u][branch_num] = []
        ele_names = self.pipe.cmd_in("python lat_ele_list " + u + '@' + branch_num)
        ele_names = ele_names.splitlines()
        for j in range(len(ele_names)):
          ele_names[j] = ele_names[j].split(';')[1]
          display_name = '(' + str(j) + ') ' + ele_names[j]
          self.e_display_list[u][branch_num].append(display_name)
        self.e_name_list[u][branch_num] = ele_names

    # Set up frames to structure the window
    self.top_frame = tk.Frame(self.upper_frame)  # Holds the selection widgets
    self.head_frame = tk.Frame(self.upper_frame) # Holds the general info
    #self.body_frame = tk.Frame(self) # Holds everything else
    self.top_frame.pack(fill="both", expand=1)
    for i in range(5):
      self.top_frame.grid_columnconfigure(i, weight=1)
    self.head_frame.pack(fill="both", expand=1)
    for i in range(4):
      self.head_frame.grid_columnconfigure(i, weight=1, pad=5)
    #self.body_frame.pack(fill="both", expand=1)

    # The Element selection fields
    tk.Label(self.top_frame, text="Universe").grid(row=0, column=0)
    tk.Label(self.top_frame, text="Branch").grid(row=0, column=1)
    self.ele_label = tk.StringVar()
    tk.Label(self.top_frame, textvariable=self.ele_label).grid(row=0, column=2)
    tk.Label(self.top_frame, text="Base/Model/Design").grid(row=0, column=3)

    self.uni = tk.StringVar()
    self.branch = tk.StringVar()
    self.branch_name = tk.StringVar()
    self.ele = tk.StringVar()
    self.bmd = tk.StringVar() # Base/Model/Design
    if (default != None):
      self.uni.set(str(default[0]))
      self.branch.set(str(default[1]))
      self.branch.set(str(default[2]))
      self.bmd.set(str(default[3]))
    else:
      self.uni.set(str(self.u_list[0]))
      self.branch.set(str(self.b_list[self.uni.get()][0]))
      self.ele.set('0')
      self.bmd.set("Base")
    self.branch_name.set(self.b_name_list[self.uni.get()][0])
    try: # in case self.ele was set by name and not by index
      ele_num = int(self.ele.get())
      self.ele.set(self.e_display_list[self.uni.get()][self.branch.get()][ele_num])
    except ValueError:
      pass
    self.ele_label.set("Element (0 to " + str(self.e_list[self.uni.get()][self.branch.get()][-1]) + ")")

    self.uni_chooser = tk.OptionMenu(self.top_frame, self.uni, *self.u_list, command=self.make_branch)
    self.branch_chooser = tk.OptionMenu(self.top_frame, self.branch_name, *self.b_name_list[self.uni.get()], command=self.make_ele_chooser)
    #self.ele_chooser = tk.Entry(self.top_frame, textvariable=self.ele)
    self.ele_chooser = ttk.Combobox(self.top_frame, textvariable=self.ele, values=self.e_display_list[self.uni.get()][self.branch.get()])
    self.ele_chooser.bind("<<ComboboxSelected>>", self.refresh)
    self.ele_chooser.bind("<Return>", self.refresh)
    self.bmd_chooser = tk.OptionMenu(self.top_frame, self.bmd, "Base", "Model", "Design")

    self.uni_chooser.grid(row=1, column=0, sticky='EW')
    self.branch_chooser.grid(row=1, column=1, sticky='EW')
    self.ele_chooser.grid(row=1, column=2, sticky='EW')
    self.bmd_chooser.grid(row=1, column=3, sticky='EW')
    tk.Button(self.top_frame, text="Load", command=self.refresh).grid(row=0, column=4, rowspan=2, sticky="NSEW")

    self.refresh()

  def make_branch(self, event=None):
    '''
    This is necessary because the list of options
    to show in the branch chooser depends on what
    universe you are in
    '''
    self.branch_chooser.destroy()
    self.branch_chooser = tk.OptionMenu(self.top_frame, self.branch_name, *self.b_name_list[self.uni.get()], command=self.make_ele_chooser)
    self.branch.set(self.b_list[self.uni.get()][0])
    self.branch_name.set(self.b_name_list[self.uni.get()][0])
    self.branch_chooser.grid(row=1, column=1, sticky='EW')
    self.make_ele_chooser()

  def make_ele_chooser(self, event=None):
    '''
    This is necessary because different branches
    have different lists of elements to display
    '''
    self.ele_chooser.destroy()
    self.ele_chooser = ttk.Combobox(self.top_frame, textvariable=self.ele, values=self.e_display_list[self.uni.get()][self.branch.get()])
    self.ele_chooser.bind("<<ComboboxSelected>>", self.refresh)
    self.ele_chooser.bind("<Return>", self.refresh)
    self.ele_label.set("Element (0 to " + str(self.e_list[self.uni.get()][self.branch.get()][-1]) + ")")
    # set self.ele to element 0 in the first branch
    self.ele.set(self.e_display_list[self.uni.get()][self.branch.get()][0])
    self.ele_chooser.grid(row=1, column=2, sticky='EW')

  def refresh(self, event=None):
    '''
    This is where most of the element information is actually created
    '''
    # Make sure the element field has an actual element in it
    if (self.ele.get() not in self.e_name_list[self.uni.get()][self.branch.get()]) \
        & (self.ele.get() not in self.e_display_list[self.uni.get()][self.branch.get()]):
      try:
        if int(self.ele.get()) not in self.e_list[self.uni.get()][self.branch.get()]:
          messagebox.showwarning("Error", "Element not found")
          return
      except ValueError:
        messagebox.showwarning("Error", "Element not found")
        return

    # Set self.ele, in case the element
    # was specified by name or display name
    if self.ele.get() in self.e_name_list[self.uni.get()][self.branch.get()]:
      ele_num = self.e_name_list[self.uni.get()][self.branch.get()].index(self.ele.get())
      self.ele.set(str(ele_num))
    elif self.ele.get() in self.e_display_list[self.uni.get()][self.branch.get()]:
      ele_num = self.e_display_list[self.uni.get()][self.branch.get()].index(self.ele.get())
      self.ele.set(str(ele_num))

    # Clear existing window contents
    for child in self.head_frame.winfo_children():
      child.destroy()
    for child in self.list_frame.winfo_children():
      child.destroy()

    # Create an element object
    self.element = lat_element(self.uni.get(), self.branch.get(), self.ele.get(), (self.bmd.get()).lower(), self.pipe)

    # Populate self.head_frame
    self.head_tk_tao_params = []
    for p in self.element.params.keys():
      self.element.params[p] = tk_tao_parameter(self.element.params[p], self.head_frame, self.pipe)
      if self.element.params[p].param.can_vary:
        self.head_tk_tao_params.append(self.element.params[p])
    #name = tk.Label(self.head_frame, text=self.element.params["name"].param.value)
    name = tk.Label(self.head_frame, text=self.element.params["name"].param.value)
    name.configure(font=("Sans", 16, "bold"))
    name.grid(row=0, column = 0, columnspan = 4)

    #Fixed parameters
    tk.Label(self.head_frame, text="Key").grid(row=1, column=0, sticky='E')
    self.element.params["key"].tk_wid.grid(row=1, column=1, sticky='EW')
    tk.Label(self.head_frame, text="s").grid(row=2, column=0, sticky='E')
    self.element.params["s"].tk_wid.grid(row=2, column=1, sticky='EW')
    tk.Label(self.head_frame, text="s_start").grid(row=3, column=0, sticky='E')
    self.element.params["s_start"].tk_wid.grid(row=3, column=1, sticky='EW')
    tk.Label(self.head_frame, text="Ref time").grid(row=4, column=0, sticky='E')
    self.element.params["ref_time"].tk_wid.grid(row=4, column=1, sticky='EW')

    #Variable parameters
    tk.Label(self.head_frame, text="Type").grid(row=1, column=2, sticky='E')
    self.element.params["type"].tk_wid.grid(row=1, column=3, sticky='EW')
    tk.Label(self.head_frame, text="Alias").grid(row=2, column=2, sticky='E')
    self.element.params["alias"].tk_wid.grid(row=2, column=3, sticky='EW')
    tk.Label(self.head_frame, text="Description").grid(row=3, column=2, sticky='E')
    self.element.params["descrip"].tk_wid.grid(row=3, column=3, sticky='EW')
    tk.Label(self.head_frame, text="is_on").grid(row=4, column=2, sticky='E')
    self.element.params["is_on"].tk_wid.grid(row=4, column=3, sticky='EW')

    # Body frame
    self.sh_b_list = [] # Show/hide button list
    self.tao_lists = [] # tk_tao_parameters for each frame
    self.p_frames = [] # Parameter frames
    i = 0 #counter for above list indices
    for key in self.element.has.keys():
      if self.element.has[key].value: #If the element has the given property
        if key == "multipoles_elec": # multipoles_elec is called
          key = "elec_multipoles"    # elec_multipoles in tao_python_cmd.f90
        if key == "mat6": # mat6 not yet implemented`
          continue
        if key == "floor": # floor currently broken
          continue
        if key == "lord_slave": # also currently broken
          continue
        # Make a button
        self.sh_b_list.append(tk.Button(self.list_frame, text=key))
        tao_list = self.pipe.cmd_in("python ele:"
            + key + ' ' + self.uni.get() + '@'
            + self.branch.get() + '>>'
            + self.ele.get() + '|'
            + (self.bmd.get()).lower())
        self.tao_lists.append(tao_list.splitlines())
        for j in range(len(self.tao_lists[i])):
          self.tao_lists[i][j] = str_to_tao_param(self.tao_lists[i][j])
        # Configure the buttons with commands
        self.sh_b_list[i].configure(command=self.s_callback(i))
        self.p_frames.append(tao_parameter_frame(self.list_frame, self.tao_lists[i], 2, self.pipe))
        self.sh_b_list[i].grid(row=2*i, column=0, sticky='EW')
        #self.p_frames[i].pack()
        i = i+1

    # Reset self.ele to be the display name
    self.ele.set(self.e_display_list[self.uni.get()][self.branch.get()][int(self.ele.get())])
    self.list_frame.columnconfigure(0, weight=1)

  def s_callback(self, index):
    return lambda : self.show(index)

  def show(self, index):
    '''
    Grids the parameter frame with given index to the body frame
    '''
    self.p_frames[index].grid(row=2*index+1, column=0, sticky='EW')
    # Button should now hide instead of show
    self.sh_b_list[index].configure(command=self.h_callback(index))

  def h_callback(self, index):
    return lambda : self.hide(index)

  def hide(self, index):
    '''
    Un-grids the paramter frame with the given index
    '''
    self.p_frames[index].grid_forget()
    # Button should now show instead of hide
    self.sh_b_list[index].configure(command=self.s_callback(index))














