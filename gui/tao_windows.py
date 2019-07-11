import tkinter as tk
import ttk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import font
import sys
import os
import copy
if 'ACC_LOCAL_DIR' in os.environ.keys():
  sys.path.append(os.environ['ACC_LOCAL_DIR']+'/tao/python/tao_pexpect')
else:
  sys.path.append(os.environ['ACC_ROOT_DIR']+'/tao/python/tao_pexpect')
from tao_widget import *
from taoplot import taoplot
from parameters import str_to_tao_param
from elements import *
from tao_set import tao_set
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backends._backend_tk import FigureManagerTk
from matplotlib.backend_bases import key_press_handler
from tao_ele_location import in_element


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

    self.canvas=tk.Canvas(self.outer_frame)
    self.list_frame=tk.Frame(self.canvas)
    scrollbar=tk.Scrollbar(self.outer_frame,orient="vertical",command=self.canvas.yview)
    self.canvas.configure(yscrollcommand=scrollbar.set)

    def scrollhelper(event):
      self.canvas.configure(scrollregion=self.canvas.bbox("all")) #,width=200,height=200)
      new_width = event.width + 15
      # Don't resize to a smaller size
      old_geo = self.wm_geometry()
      old_width = int(old_geo[:old_geo.find('x')])
      old_height = int(old_geo[old_geo.find('x')+1:old_geo.find('+')])
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
      if old_height != 1:
        #new_geo = old_geo[old_geo.find('x'):old_geo.find('+')]
        new_geo = old_geo[old_geo.find('x'):]
        new_geo = str(new_width) + new_geo
      else:
        new_geo = str(new_width) + 'x500'
      # Place the new window near the root window
      if new_geo.find('+') == -1:
        new_pos = self.root.wm_geometry()
        new_pos = new_pos[new_pos.find('+'):]
        new_geo = new_geo + new_pos
      self.geometry(new_geo)
      self.canvas.configure(width=new_width)
    self.list_frame.bind("<Configure>",scrollhelper)

    self.outer_frame.bind("<Enter>", self.bind_mouse)
    self.outer_frame.bind("<Leave>", self.unbind_mouse)

    if use_upper:
      self.upper_frame.pack(side="top",fill="both", expand=0)
    self.outer_frame.pack(side="top",fill="both",expand=1)
    scrollbar.pack(side="right",fill="y")
    self.canvas.pack(side="left",fill="both",expand=1)

    self.canvas_list_window = self.canvas.create_window((0,0),window=self.list_frame,anchor='nw')

  def bind_mouse(self, event):
    self.outer_frame.bind_all("<Button-4>", self.mouse_scroll)
    self.outer_frame.bind_all("<Button-5>", self.mouse_scroll)

  def unbind_mouse(self, event):
    self.outer_frame.unbind_all("<Button-4>")
    self.outer_frame.unbind_all("<Button-5>")

  def mouse_scroll(self, event):
    #self.canvas.yview_scroll(direction,"units")
    if event.num == 4:
      self.canvas.yview_scroll(-1,"units")
    elif event.num == 5:
      self.canvas.yview_scroll(1,"units")


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
    # Don't count units# in list length
    real_len = len(self.tao_list)
    for item in self.tao_list:
      if item.param.name.find('units#') == 0:
        real_len -= 1
    # real_len now actually represents how many items there are
    len_col = int(real_len/n_col) # How many elements per column
    if len_col < len(self.tao_list)/n_col:
      #if the result was rounded down
      len_col = len_col+1
    #for i in range(n_col):
    #  cols.append(self.tao_list[i*len_col : (i+1)*len_col])

    # Grid widgets (and units)
    i = 0
    for item in self.tao_list:
      if item.param.name.find('units#') != 0:
        r = (i % len_col) + 1 # row 0 reserved for titles
        c = int(i / len_col) * 3
        item.tk_label.grid(row=r, column=c, sticky='E')
        item.tk_wid.grid(row=r, column=c+1, sticky='EW')
        i = i+1
      else:
        tk.Label(self, text=item.param.value).grid(row=r, column=c+2, sticky='W')

    # Grid the widgets in each column
    #j = 0 #controls which column the widgets get gridded to
    #for c in cols:
    #  for i in range(len(c)):
    #    c[i].tk_label.grid(row=i+1, column=2*j, sticky = 'E')
    #    c[i].tk_wid.grid(row=i+1, column=2*j+1, sticky = 'EW')
    #  j = j+1

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
class lw_table_window(tk.Toplevel):
  '''
  Light-weight version of table_window, intended to be less
  graphically/computationally intensive (for use with networked
  machines, for example).  Table is read only but individual
  entries can still be opened for editing with double click.
  '''
  def __init__(self, root, pipe, array_name, title_list, set_format, *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
    self.title(array_name)
    self.pipe = pipe
    self.array_name = array_name
    self.title_list = title_list
    self.set_format = set_format
    self.table_frame = tk.Frame(self) #holds the table
    self.table_frame.pack(fill='both', expand=1)
    self.refresh()

  def refresh(self):
    '''
    Creates a treeview to display the table information,
    and binds double click to open one item
    '''
    # Clear the existing table
    for child in self.table_frame.winfo_children():
      child.destroy()

    widths = [0]*len(self.title_list) # tracks column widths

    # Create table
    self.tree = ttk.Treeview(self.table_frame, columns=self.title_list, show='headings')
    # Column titles
    for title in self.title_list:
      self.tree.heading(title, text=title)
      self.tree.column(title, stretch=True, anchor='center')

    # Fill rows
    #Fetch and fill in the data
    # IT IS EXPECTED THAT SUBCLASSES WILL DEFINE
    # self.row_num AND A METHOD FOR POPULATING
    # self.list_rows, AND THEN CALL THIS METHOD
    #i = row counter, j = column counter
    self.list_rows = []
    for i in range(self.row_num):
      row = self.lw_list_row_fetch(i)
      self.tree.insert("", "end", values=row)
      for j in range(len(row)):
        if len(row[j])*12 > widths[j]:
          widths[j] = len(row[j])*12

    # Set column widths appropriately
    for j in range(len(self.title_list)):
      if len(self.title_list[j])*12 > widths[j]:
        widths[j] = len(self.title_list[j])*12
      self.tree.column(self.title_list[j], width=widths[j], minwidth=widths[j])

    # Scrollbars
    hbar = ttk.Scrollbar(self.table_frame, orient="horizontal", command=self.tree.xview)
    vbar = ttk.Scrollbar(self.table_frame, orient="vertical", command=self.tree.yview)
    self.tree.configure(xscrollcommand=hbar.set, yscrollcommand=vbar.set)

    vbar.pack(side="right", fill="y", expand=0)
    hbar.pack(side="bottom", fill='x', expand=0)
    self.tree.pack(side="left", fill="both", expand=1)
    self.widths = widths

    # Double click to open element
    self.tree.bind('<Double-Button-1>', self.lw_detail_callback)

    tot = 0
    for w in widths:
      tot = tot+w
    self.maxsize(1800, 1000)
    self.minsize(1300, 100)

  def open_detail_window_callback(self, event=None):
    '''
    Checks the table for the currently selected row and
    makes a call to open_detail_window_callback
    '''
    x = self.tree.focus()
    row = self.tree.item(x)
    self.open_detail_window(int(row['values'][0]))

  def open_detail_window(self, index):
    '''
    Opens up a detail window for the given index.
    '''
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
  '''
  With lw set to True, opens a lw_table_window instead
  '''

  def __init__(self, root, pipe, d1_data_name, u_ix, ix_lb, ix_ub, *args, **kwargs):
    self.u_ix = u_ix
    self.ix_lb = ix_lb
    self.ix_ub = ix_ub
    self.lw = root.lw
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
    if self.lw:
      lw_table_window.__init__(self, root, pipe, d1_data_name, title_list, set_format, *args, **kwargs)
    else:
      table_window.__init__(self, root, pipe, d1_data_name, title_list, bulk_template, bulk_set_format, set_format, *args, **kwargs)

  def refresh(self):
    self.d_list = self.pipe.cmd_in("python data_d_array " + self.u_ix + '@' + self.array_name)
    self.d_list = self.d_list.splitlines()
    self.row_num = len(self.d_list)
    if self.lw:
      lw_table_window.refresh(self)
    else:
      table_window.refresh(self)

  def list_row_fetch(self, index):
    '''
    Returns row number index to be appended to self.list_rows
    '''
    return d1_data_list_entry(self.list_frame, self.d_list[index])

  def lw_list_row_fetch(self, index):
    '''
    For use with lw_table_window.  Returns a list of the values
    to be displayed in each column of row #index.
    '''
    return self.d_list[index].split(';')

  def lw_detail_callback(self, event=None):
    '''
    Callback for the light-weight table to open a detail window
    '''
    lw_table_window.open_detail_window_callback(self, event)

  def open_detail_window(self, index):
    self.param_list = self.pipe.cmd_in("python data " + str(self.u_ix) + '@' + self.array_name + '[' + str(index) + ']')
    self.param_list = self.param_list.splitlines()
    for i in range(len(self.param_list)):
      self.param_list[i]=str_to_tao_param(self.param_list[i])
    if self.lw:
      lw_table_window.open_detail_window(self, index)
    else:
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
    self.lw = root.lw #set light-weight flag
    bulk_template.append([str_to_tao_param("meas_value;REAL;T;"), 2])
    bulk_template.append([str_to_tao_param("good_user;LOGIC;T;"), 6])
    bulk_template.append([str_to_tao_param("weight;REAL;T;"), 7])

    bulk_set_format = "set var {}|"
    set_format = "set var {}[{}]|"
    if self.lw:
      lw_table_window.__init__(self, root, pipe, v1_var_name, title_list, set_format, *args, **kwargs)
    else:
      table_window.__init__(self, root, pipe, v1_var_name, title_list, bulk_template, bulk_set_format, set_format, *args, **kwargs)

  def refresh(self):
    self.v_list = self.pipe.cmd_in("python var_v_array " + self.array_name)
    self.v_list = self.v_list.splitlines()
    self.row_num = len(self.v_list)
    if self.lw:
      lw_table_window.refresh(self)
    else:
      table_window.refresh(self)

  def list_row_fetch(self, index):
    '''
    Returns row number index to be appended to self.list_rows
    '''
    return v1_var_list_entry(self.list_frame, self.v_list[index])

  def lw_list_row_fetch(self, index):
    '''
    For use with lw_table_window.  Returns a list of the values
    to be displayed in each column of row #index.
    '''
    return self.v_list[index].split(';')

  def lw_detail_callback(self, event=None):
    '''
    Callback for the light-weight table to open a detail window
    '''
    lw_table_window.open_detail_window_callback(self, event)

  def open_detail_window(self, index):
    self.param_list = self.pipe.cmd_in("python var " + self.array_name + '[' + str(index) + ']')
    self.param_list = self.param_list.splitlines()
    for i in range(len(self.param_list)):
      self.param_list[i]=str_to_tao_param(self.param_list[i])
    if self.lw:
      lw_table_window.open_detail_window(self, index)
    else:
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
      self.title("Active Plots")

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
    win = tao_plot_window(self.root, self.plot, self.pipe)
    self.root.plot_windows.append(win)

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
    if self.mode == "T":
      data_list = self.pipe.cmd_in("python plot1 " + self.plot)
    elif self.mode == "R":
      data_list = self.pipe.cmd_in("python plot1 "
          + self.region_list[self.plot_list.index(self.plot)])
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
    # In matplotlib mode, apply for the appropriate region too
    c1 = self.root.plot_mode == "matplotlib"
    c2 = self.mode == "T"
    c3 = self.plot in self.root.placed.keys()
    if c1 & c2 & c3:
      set_str = "set plot " + self.root.placed[self.plot] + ' '
      tao_set(self.tao_list, set_str, self.pipe)
    # Refresh any existing plot windows
    for win in self.root.plot_windows:
      win.refresh()

#----------------------------------------------------
# Matplotlib plotting window

class tao_plot_window(tk.Toplevel):
  '''
  Displays one (perhaps multiple) matplotlib plots
  that the user has specified from the plotting
  template window that they want to plot. Creating a
  window in tkinter is necessary rather than using
  matplotlib's built in system for creating windows
  because using that system will halt the tkinter
  mainloop until the plots are closed.
  '''
  def __init__(self, root, template, pipe, *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
    self.template = template #The template plot being plotted
    self.title(template)
    self.root = root
    self.pipe = pipe
    self.fig = False #Default value

    # Check if the template has been placed in a region already, and place it if necessary
    if self.template in self.root.placed.keys():
      print("found in root.placed")
      pass
    else:
      #Place the plot in the next available region
      #and make visible
      r_index = 1
      while (("r" + str(r_index)) in self.root.placed.values()):
        r_index = r_index + 1
      self.pipe.cmd_in("place r" + str(r_index) + " " + self.template)
      self.pipe.cmd_in("set plot r" + str(r_index) + ' visible = T')
      self.root.placed[self.template] = 'r' + str(r_index)

    self.mpl = taoplot(pipe, self.root.placed[self.template])
    self.refresh()

  def refresh(self, event=None):
    '''
    Makes the call to matplotlib to draw the plot to the window
    '''
    #Clear the window
    for child in self.winfo_children():
      child.destroy()

    #Get plotting results
    self.plot_output = self.mpl.plot()

    #Get the figure
    self.fig = self.plot_output[0]

    #Get figure information
    self.fig_info = self.plot_output[1]

    #Create widgets to display the figure
    canvas = FigureCanvasTkAgg(self.fig, master=self)
    canvas.draw()
    canvas.get_tk_widget().pack(side="top", fill="both", expand=1)
    # DO NOT TOUCH
    canvas.manager = FigureManagerTk(canvas, self.fig.number, tk.Toplevel(self.root))

    toolbar = NavigationToolbar2Tk(canvas, self)
    toolbar.update()
    canvas._tkcanvas.pack(side="top", fill="both", expand=1)

    def on_key_press(event):
      key_press_handler(event, canvas, toolbar)

    canvas.mpl_connect("key_press_event", on_key_press)

    def on_click(event):
      if event.dblclick:
        eleList = in_element(event.xdata,event.ydata,self.fig_info)
        for i in eleList:
          tao_ele_window(self.root,self.pipe,default=[self.fig_info[1],self.fig_info[2],i,self.fig_info[3]])

    canvas.mpl_connect("button_press_event", on_click)

  def destroy(self):
    # Note: lat_layout should not be automatically removed from r1
    if self.template != "lat_layout":
      # Unplace the template from its region
      self.pipe.cmd_in("place " + self.root.placed[self.template] + " none")
      # Remove self from root.plot_windows
      try:
        self.root.plot_windows.pop(self.root.plot_windows.index(self))
      except ValueError: #incase the window never got added to the list
        pass
    tk.Toplevel.destroy(self)


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
    # In matplotlib mode, apply for the appropriate region too
    c1 = self.root.plot_mode == "matplotlib"
    c2 = self.mode == "T"
    plot = win.curve[:win.curve.index('.')]
    c3 = plot in self.root.placed.keys()
    if c1 & c2 & c3:
      curve_name = self.root.placed[win.curve[:win.curve.index('.')]] \
          + win.curve[win.curve.index('.'):]
      set_str = "set curve " + curve_name + ' '
      tao_set(win.tao_list, set_str, win.pipe)
    # Refresh any existing plot windows
    for win in self.root.plot_windows:
      win.refresh()

  def graph_apply(self):
    #Convert from plot.graph to index.graph
    graph_name = str(self.index) + self.graph[self.graph.index('.'):]
    set_str = "set graph @" + self.mode + graph_name + ' '
    tao_set(self.tao_list, set_str, self.pipe)
    # In matplotlib mode, apply for the appropriate region too
    c1 = self.root.plot_mode == "matplotlib"
    c2 = self.mode == "T"
    plot = self.graph[:self.graph.index('.')]
    c3 = plot in self.root.placed.keys()
    if c1 & c2 & c3:
      graph_name = self.root.placed[self.graph[:self.graph.index('.')]] \
          + self.graph[self.graph.index('.'):]
      set_str = "set graph " + graph_name + ' '
      tao_set(self.tao_list, set_str, self.pipe)
    # Refresh any existing plot windows
    for win in self.root.plot_windows:
      win.refresh()

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
# Branch/element choosing widgets
class tao_branch_widgets:
  '''
  Provides several widgets for slecting universe,
  branch, and elements
  Available widgets:
  self.uni_chooser: OptionMenu to pick the universe index
  self.branch_chooser: OptionMenu to pick the branch (displays name and index)
  self.ele_chooser: ttk Combobox for selecting an element.  Maybe be specified by name or index
  self.bmd_chooser: OptionMenu for choosing base, model, or design
  self.bme_chooser: OptionMenu for choosing beginning, middle, or end
  In addition, the class provides tk variables for each of these widgets, named the same but without _chooser at the end
  '''
  def __init__(self, parent, pipe, default=None):
    '''
    parent: the parent widget where these widgets will be placed
    pipe: tao_interface object
    default: specify the default state of the widgets
    '''
    self.parent = parent
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

    # The variables and widgets
    self.uni = tk.StringVar()
    self.branch = tk.StringVar()
    self.branch_name = tk.StringVar()
    self.ele = tk.StringVar()
    self.ele_label = tk.StringVar() # For showing index range
    self.bmd = tk.StringVar() # Base/Model/Design
    self.bme = tk.StringVar() # Beginning/Middle/End
    if (default != None):
      self.uni.set(str(default[0]))
      self.branch.set(str(default[1]))
      self.ele.set(str(default[2]))
      self.bmd.set(str(default[3]))
    else:
      self.uni.set(str(self.u_list[0]))
      self.branch.set(str(self.b_list[self.uni.get()][0]))
      self.ele.set('0')
      self.bmd.set("Model")
    self.branch_name.set(self.b_name_list[self.uni.get()][0])
    try: # in case self.ele was set by name and not by index
      ele_num = int(self.ele.get())
      self.ele.set(self.e_display_list[self.uni.get()][self.branch.get()][ele_num])
    except ValueError:
      pass
    self.ele_label.set("Element (0 to " + str(self.e_list[self.uni.get()][self.branch.get()][-1]) + ")")

    self.uni_chooser = tk.OptionMenu(self.parent, self.uni, *self.u_list, command=self.make_branch)
    self.branch_chooser = tk.OptionMenu(self.parent, self.branch_name, *self.b_name_list[self.uni.get()], command=self.make_ele_chooser)
    self.ele_chooser = ttk.Combobox(self.parent, textvariable=self.ele, values=self.e_display_list[self.uni.get()][self.branch.get()])
    #self.ele_chooser.bind("<<ComboboxSelected>>", self.refresh)
    #self.ele_chooser.bind("<Return>", self.refresh)
    self.bmd_chooser = tk.OptionMenu(self.parent, self.bmd, "Base", "Model", "Design")
    self.bme_chooser = tk.OptionMenu(self.parent, self.bme, "Beginning", "Middle", "End")
    self.bme.set("End")

  def make_branch(self, event=None):
    '''
    This is necessary because the list of options
    to show in the branch chooser depends on what
    universe you are in
    '''
    # Get the current geometry manager and properties
    manager = self.branch_chooser.winfo_manager()
    if manager == 'pack':
      props = self.branch_chooser.pack_info()
    elif manager == 'grid':
      props = self.branch_chooser.grid_info()
    elif manager == 'place':
      props = self.branch_chooser.place_info()
    else: #No need to remake the branch_chooser
      return

    # Remake the branch_chooser
    self.branch_chooser.destroy()
    self.branch_chooser = tk.OptionMenu(self.parent, self.branch_name, *self.b_name_list[self.uni.get()], command=self.make_ele_chooser)
    self.branch.set(self.b_list[self.uni.get()][0])
    self.branch_name.set(self.b_name_list[self.uni.get()][0])
    self.make_ele_chooser()

    # Put the branch_chooser back correctly
    if manager == 'pack':
      self.branch_chooser.pack(**props)
    if manager == 'grid':
      self.branch_chooser.grid(**props)
    if manager == 'place':
      self.branch_chooser.place(**props)

  def make_ele_chooser(self, event=None):
    '''
    This is necessary because different branches
    have different lists of elements to display
    '''
    # Get the current geometry manager and properties
    manager = self.ele_chooser.winfo_manager()
    if manager == 'pack':
      props = self.ele_chooser.pack_info()
    elif manager == 'grid':
      props = self.ele_chooser.grid_info()
    elif manager == 'place':
      props = self.ele_chooser.place_info()
    else: #No need to remake the branch_chooser
      return

    self.ele_chooser.destroy()
    self.ele_chooser = ttk.Combobox(self.parent, textvariable=self.ele, values=self.e_display_list[self.uni.get()][self.branch.get()])
    self.ele_chooser.bind("<<ComboboxSelected>>", self.update)
    self.ele_chooser.bind("<Return>", self.update)
    self.ele_label.set("Element (0 to " + str(self.e_list[self.uni.get()][self.branch.get()][-1]) + ")")
    # set self.ele to element 0 in the first branch
    self.ele.set(self.e_display_list[self.uni.get()][self.branch.get()][0])

    # Put the branch_chooser back correctly
    if manager == 'pack':
      self.ele_chooser.pack(**props)
    if manager == 'grid':
      self.ele_chooser.grid(**props)
    if manager == 'place':
      self.ele_chooser.place(**props)

  def update(self, event=None):
    '''
    Sets self.ele to its index value (as a string)
    Returns 0 if setting self.ele failed, 1 if succeeded
    '''
    # Make sure the element field has an actual element in it
    if (self.ele.get() not in self.e_name_list[self.uni.get()][self.branch.get()]) \
        & (self.ele.get() not in self.e_display_list[self.uni.get()][self.branch.get()]):
      try:
        if int(self.ele.get()) not in self.e_list[self.uni.get()][self.branch.get()]:
          messagebox.showwarning("Error", "Element not found")
          return 0
      except ValueError:
        messagebox.showwarning("Error", "Element not found")
        return 0

    # Set self.ele, in case the element
    # was specified by name or display name
    if self.ele.get() in self.e_name_list[self.uni.get()][self.branch.get()]:
      ele_num = self.e_name_list[self.uni.get()][self.branch.get()].index(self.ele.get())
      self.ele.set(str(ele_num))
    elif self.ele.get() in self.e_display_list[self.uni.get()][self.branch.get()]:
      ele_num = self.e_display_list[self.uni.get()][self.branch.get()].index(self.ele.get())
      self.ele.set(str(ele_num))
    return 1

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
    self.open_frames = [] #keeps track of which p_frames are open

    # Set up frames to structure the window
    self.top_frame = tk.Frame(self.upper_frame)  # Holds the selection widgets
    self.head_frame = tk.Frame(self.upper_frame) # Holds the general info
    #self.body_frame = tk.Frame(self) # Holds everything else
    self.top_frame.pack(fill="both", expand=1)
    for i in range(5):
      self.top_frame.grid_columnconfigure(i, weight=1)
    self.head_frame.pack(fill="both", expand=1)
    for i in range(4):
      self.head_frame.grid_columnconfigure(i, weight=i%2, pad=5)
    #self.body_frame.pack(fill="both", expand=1)

    # The Element selection fields
    self.ele_wids = tao_branch_widgets(self.top_frame, self.pipe, self.default)
    tk.Label(self.top_frame, text="Universe").grid(row=0, column=0)
    tk.Label(self.top_frame, text="Branch").grid(row=0, column=1)
    self.ele_label = self.ele_wids.ele_label #Might not work
    tk.Label(self.top_frame, textvariable=self.ele_label).grid(row=0, column=2)
    tk.Label(self.top_frame, text="Base/Model/Design").grid(row=0, column=3)

    # Configure and place widgets
    self.ele_wids.ele_chooser.bind("<<ComboboxSelected>>", self.refresh)
    self.ele_wids.ele_chooser.bind("<Return>", self.refresh)
    self.ele_wids.uni_chooser.grid(row=1, column=0, sticky='EW')
    self.ele_wids.branch_chooser.grid(row=1, column=1, sticky='EW')
    self.ele_wids.ele_chooser.grid(row=1, column=2, sticky='EW')
    self.ele_wids.bmd_chooser.grid(row=1, column=3, sticky='EW')
    tk.Button(self.top_frame, text="Load", command=self.refresh).grid(row=0, column=4, rowspan=2, sticky="NSEW")

    self.refresh()

  def refresh(self, event=None):
    '''
    This is where most of the element information is actually created
    '''
    # Update the ele_wids and check for successful update
    if not self.ele_wids.update():
      messagebox.showwarning("Error", "Element not found")
      return

    # Clear existing window contents
    for child in self.head_frame.winfo_children():
      child.destroy()
    for child in self.list_frame.winfo_children():
      child.destroy()

    # Create an element object
    self.element = lat_element(self.ele_wids.uni.get(), self.ele_wids.branch.get(), self.ele_wids.ele.get(), (self.ele_wids.bmd.get()).lower(), self.pipe)

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
    self.p_names = [] # Parameter frame names
    i = 0 #counter for above list indices
    for key in self.element.has.keys():
      if self.element.has[key].value: #If the element has the given property
        if key == "multipoles_elec": # multipoles_elec is called
          key = "elec_multipoles"    # elec_multipoles in tao_python_cmd.f90
        if key == "mat6": # mat6 not yet implemented`
          continue
        if key in ["multipoles", "elec_multipoles"]:
          continue #not yet implemented
        if key == "lord_slave": # extremely special case
          self.sh_b_list.append(tk.Button(self.list_frame, text=key))
          ls_frame = tk.Frame(self.list_frame)
          ls_list = self.pipe.cmd_in("python ele:lord_slave "
              + self.ele_wids.uni.get() + '@' + self.ele_wids.branch.get()
              + '>>' + self.ele_wids.ele.get() + '|' + (self.ele_wids.bmd.get()).lower())
          ls_list = ls_list.splitlines()
          self.tao_lists.append(ls_list) # Don't want to misalign indices
          ls_cols = ['Index', 'Name', 'Key', 'Slave/Lord Status']
          # Table:
          ls_tree = ttk.Treeview(ls_frame, columns=ls_cols)
          # Column configuration
          ls_tree.heading('#0', text='Lord/Slave')
          ls_tree.column('#0', stretch=True, anchor='center')
          for title in ls_cols:
            ls_tree.heading(title, text=title)
            ls_tree.column(title, stretch=True, anchor='center')
          # Fill rows with data
          for line in ls_list:
            line = line.split(';')
            if line[0] == 'Element':
              current_level = ls_tree.insert("", "end", text=line[0], values=line[1:])
            else:
              ls_tree.insert(current_level, "end", text=line[0], values=line[1:])
          # Fix scrolling
          ls_tree.bind('<Enter>', self.unbind_mouse)
          ls_tree.bind('<Leave>', self.bind_mouse)
          # Double click to open new element window
          def open_ele_window(event=None):
            '''
            Opens an element window for the currently selected row
            '''
            x = ls_tree.focus()
            row = ls_tree.item(x)
            bele = row['values'][0] #branch>>ele_ix
            settings = [self.ele_wids.uni.get(),
                bele[:bele.find('>>')],
                bele[bele.find('>>')+2:],
                self.ele_wids.bmd.get()]
            win = tao_ele_window(self.root, self.pipe, settings)
            return('break')
          ls_tree.bind('<Double-Button-1>', open_ele_window)
          ls_tree.pack(fill='both', expand=1)
          self.p_frames.append(ls_frame)
          self.p_names.append(key)
          self.sh_b_list[i].configure(command=self.s_callback(i))
          self.sh_b_list[i].grid(row=2*i, column=0, sticky='W')
          i = i+1
          continue

        # GENERIC CASE

        # Make a button
        self.sh_b_list.append(tk.Button(self.list_frame, text=key))
        tao_list = self.pipe.cmd_in("python ele:"
            + key + ' ' + self.ele_wids.uni.get() + '@'
            + self.ele_wids.branch.get() + '>>'
            + self.ele_wids.ele.get() + '|'
            + (self.ele_wids.bmd.get()).lower())
        self.tao_lists.append(tao_list.splitlines())
        for j in range(len(self.tao_lists[i])):
          self.tao_lists[i][j] = str_to_tao_param(self.tao_lists[i][j])
        # Configure the buttons with commands
        self.sh_b_list[i].configure(command=self.s_callback(i))
        if key == 'floor': #should just have one column
          n_cols = 1
        else:
          n_cols = 2
        self.p_frames.append(tao_parameter_frame(self.list_frame, self.tao_lists[i], n_cols, self.pipe))
        self.p_names.append(key)
        if key == 'floor': #column titles
          title_frame = tk.Frame(self.p_frames[i])
          tk.Label(title_frame, text='X').grid(row=0, column=0)
          tk.Label(title_frame, text='Y').grid(row=0, column=1)
          tk.Label(title_frame, text='Z').grid(row=0, column=2)
          tk.Label(title_frame, text='Theta').grid(row=0, column=3)
          tk.Label(title_frame, text='Phi').grid(row=0, column=4)
          tk.Label(title_frame, text='Psi').grid(row=0, column=5)
          for x in range(6):
            title_frame.columnconfigure(x, weight=1)
          title_frame.grid(row=0, column=1, sticky='EW')
        self.sh_b_list[i].grid(row=2*i, column=0, sticky='W')
        #self.p_frames[i].pack()
        i = i+1

    # Reset self.ele to be the display name
    self.ele_wids.ele.set(self.ele_wids.e_display_list[self.ele_wids.uni.get()][self.ele_wids.branch.get()][int(self.ele_wids.ele.get())])
    self.list_frame.columnconfigure(0, weight=1)

    # Open the frames listed in self.open_frames
    for i in range(len(self.p_names)):
      if self.p_names[i] in self.open_frames:
        self.show(i)

  def s_callback(self, index):
    return lambda : self.show(index)

  def show(self, index):
    '''
    Grids the parameter frame with given index to the body frame
    '''
    self.p_frames[index].grid(row=2*index+1, column=0, sticky='EW')
    # Allow lord_slave to stretch to bottom of window
    if self.p_names[index] == 'lord_slave':
      self.list_frame.rowconfigure(index, weight=1)
    # Add to self.open_frames, if not already there
    if self.p_names[index] not in self.open_frames:
      self.open_frames.append(self.p_names[index])
    # Button should now hide instead of show
    self.sh_b_list[index].configure(command=self.h_callback(index))

  def h_callback(self, index):
    return lambda : self.hide(index)

  def hide(self, index):
    '''
    Un-grids the paramter frame with the given index
    '''
    self.p_frames[index].grid_forget()
    # Remove from self.open_frames (if possible)
    if self.p_names[index] in self.open_frames:
      self.open_frames.pop(self.open_frames.index(self.p_names[index]))
    # Button should now show instead of hide
    self.sh_b_list[index].configure(command=self.s_callback(index))

#---------------------------------------------------
# Lattice Window
class tao_lattice_window(tk.Toplevel):
  '''
  Shows lattice elements in a read-only table view
  with an interface to select which rows/columns
  are displayed
  '''
  # TODO: replace blank
  def __init__(self, root, pipe, switches=""):
    tk.Toplevel.__init__(self, root)
    self.title("Lattice")
    self.root = root
    self.pipe = pipe
    self.switches = switches
    self.bind("<Return>", self.refresh)

    # Set up top_frame for switches and table_frame for the table
    self.top_frame = tk.Frame(self)
    self.top_frame.pack(fill="both", expand=0)
    self.top_frame.columnconfigure(6, weight=1)
    self.table_frame = tk.Frame(self)
    self.table_frame.pack(fill="both", expand=1)

    # Switches

    # Load/Save table template
    tk.Label(self.top_frame, text="Template File: ").grid(row=4, column=0, sticky='W')
    self.template_file = tk_tao_parameter(str_to_tao_param("template_file;FILE;T;"), self.top_frame, self.pipe)
    self.template_file.tk_wid.configure(command=self.temp_file_load)
    self.template_file.tk_wid.grid(row=4, column=1, sticky='EW')
    tk.Label(self.top_frame, text="Template:").grid(row=4, column=2, sticky='W')
    self.temp_var = tk.StringVar() # Holds the chosen template
    self.temp_var.set("NO TEMPLATE FILE SELECTED")
    self.temp_chooser = tk.OptionMenu(self.top_frame, self.temp_var, [])
    self.temp_chooser.grid(row=4, column=3, sticky='EW')
    tk.Button(self.top_frame, text="Save", command=self.save_template).grid(row=4, column=5, sticky='W')
    self.temp_save = tk_tao_parameter(str_to_tao_param("name;STR;T;"), self.top_frame, self.pipe)
    self.temp_save.tk_wid.bind('<Return>', self.save_template)
    self.temp_save.tk_wid.grid(row=4, column=4, sticky='EW')

    # Branch/General
    self.branch_wids = tao_branch_widgets(self.top_frame, self.pipe)
    #tk.Label(self.top_frame, text="PLACEHOLDER").grid(row=0, column=0, columnspan=5, sticky='EW')
    tk.Label(self.top_frame, text="Universe:").grid(row=0, column=0, sticky='W')
    #self.branch_wids.uni_chooser.configure(command=self.set_uni)
    self.branch_wids.uni_chooser.grid(row=0, column=1, sticky='EW')
    tk.Label(self.top_frame, text="Branch:").grid(row=0, column=2, sticky='W')
    self.branch_wids.branch_chooser.grid(row=0, column=3, sticky='EW')
    self.branch_wids.bmd_chooser.grid(row=0, column=4, sticky='EW')
    self.branch_wids.bme_chooser.grid(row=0, column=5, sticky='EW')


    # Columns
    tk.Label(self.top_frame, text="Columns:").grid(row=1, column=0, sticky='W')
    self.col_filter = tk.StringVar()
    col_opts = ["Default", "Floor Coordinates",
        "Orbit", "Spin", "Orbit + Spin",
        "Radiation Integrals",
        "Cumulative Radiation Integrals",
        "List Attributes",
        "Use Template File"]
    self.col_filter.set("Default")
    self.cf_box = tk.OptionMenu(self.top_frame, self.col_filter,
        *col_opts, command=self.col_filter_callback)
    self.cf_box.grid(row=1, column=1, sticky='EW')
    tk.Label(self.top_frame, text="Attributes:").grid(row=1, column=2, sticky='EW')
    self.col_atts = tk.StringVar()
    self.att_box = tk.Entry(self.top_frame, textvariable=self.col_atts)
    self.att_box.configure(state="disabled")
    self.att_box.grid(row=1, column=3, sticky='EW')
    tk.Label(self.top_frame, text="Template file:").grid(row=1, column=4, sticky='EW')
    #tk.Button(self.top_frame, text="Browse... (placeholder)").grid(row=1, column=5, sticky='EW')
    self.col_file = tk_tao_parameter(str_to_tao_param("col_file;FILE;T;"), self.top_frame, self.pipe)
    self.col_file.tk_wid.configure(state="disabled")
    self.col_file.tk_wid.configure(disabledforeground="grey")
    self.col_file.tk_wid.grid(row=1, column=5, sticky='EW')

    # Rows
    tk.Label(self.top_frame, text="Rows:").grid(row=2, column=0, sticky='W')
    self.f_button = tk.Button(self.top_frame, text="Filters", command=self.open_filter_menu)
    self.f_button.grid(row=2, column=1, sticky='EW')

    self.filter_vars = []
    for i in range(5):
      self.filter_vars.append(tk.BooleanVar())
      self.filter_vars[i].set(False)
    self.remove_if_zero = tk.StringVar()
    self.s_range = tk.StringVar()

    tk.Label(self.top_frame, text="Element List:").grid(row=2, column=2, sticky='E')
    self.ele_list = tk.StringVar()
    self.ele_list_opt = tk.StringVar()
    ele_list_opts = ["Default (first 200)", "All", "Tracking elements", "Custom"]
    self.ele_list_chooser = tk.OptionMenu(self.top_frame, self.ele_list_opt, *ele_list_opts, command=self.ele_list_callback)
    self.ele_list_opt.set(ele_list_opts[0])
    self.ele_list_chooser.grid(row=2, column=3, sticky='EW')
    self.ele_list_box = tk.Entry(self.top_frame, textvariable=self.ele_list)
    self.ele_list_box.configure(state="disabled")
    self.ele_list_box.grid(row=2, column=4, columnspan=3, sticky='EW')

    # Advanced Options
    self.use_advanced = False
    tk.Button(self.top_frame, text="Advanced\nOn/Off", command=self.toggle_advanced).grid(row=2, column=7, rowspan=2, sticky='NSEW')
    tk.Label(self.top_frame, text="Advanced:").grid(row=3, column=0, sticky='W')
    self.advanced_var = tk.StringVar()
    self.advanced_var.set(self.switches)
    self.advanced_box = tk.Entry(self.top_frame, textvariable=self.advanced_var)
    self.advanced_box.configure(state="disabled")
    self.advanced_box.grid(row=3, column=1, columnspan=6, sticky='EW') #pack(side="left", fill="both", expand=1)
    #self.advanced_box.bind("<Return>", self.refresh)
    #self.tree.focus() and self.tree.selection() for current items

    b = tk.Button(self.top_frame, text="Refresh", command=self.refresh)
    b.grid(row=0, column=7, rowspan=2, sticky='NSEW') #pack(side="left", fill="both", expand=0)
    self.refresh()

  def save_template(self, event=None):
    '''
    Writes the current switches to the file in self.template_file
    using the name given in self.temp_save.tk_var (if any)
    '''
    if self.template_file.tk_var.get() == "Browse...":
      return
    #Don't save if no switches have been set
    if self.switches == '':
      return
    t_file = open(self.template_file.tk_var.get(), mode='a')
    # Write to file and update self.temp_dict
    if self.temp_save.tk_var.get() != '':
      t_file.write("name:" + self.temp_save.tk_var.get() + '\n')
      self.temp_dict[self.temp_save.tk_var.get()] = self.switches
    else:
      self.temp_dict[self.switches] = self.switches
    t_file.write(self.switches + '\n')
    t_file.close()
    # Remake the template chooser to list the new template
    self.temp_chooser.destroy()
    temp_opts = list(self.temp_dict.keys())
    self.temp_var.set(temp_opts[-1])
    self.temp_chooser = tk.OptionMenu(self.top_frame, self.temp_var, *temp_opts, command=self.temp_chooser_callback)
    self.temp_chooser.grid(row=4, column=3, sticky='EW')

  def temp_file_load(self, event=None):
    '''
    Tries to load the specified template file, and creates an OptionMenu to pick from defined templates if successful.
    Also creates a save template button and entry for the name
    self.temp_dict is a dictionary whose keys are the names of templates and whose values are the switches for those templates
    '''
    # First run the tk_tao_paramter open file method that is being overloaded here
    tk_tao_parameter.open_file(self.template_file)
    if self.template_file.tk_var.get() == "Browse...":
      return
    # Attempt to parse the file for templates
    t_file = open(self.template_file.tk_var.get())
    templates = t_file.read().splitlines()
    t_file.close()
    self.temp_dict = {}
    given_name = False
    for item in templates:
      #Skip the line if it starts with # (comment)
      if item.find('#') == 0:
        continue
      #Set the name flag if line starts with name:
      if item.find('name:') == 0:
        given_name = True
        name = item[5:] #strip off the name: part
        name = name.strip() # remove trailing whitespace
        continue
      #Add the current item to self.temp_dict
      if given_name:
        given_name = False
        self.temp_dict[name] = item
      else:
        self.temp_dict[item] = item
    # Stop here if temp_dict is empty
    if self.temp_dict == {}:
      return
    #Make the template chooser widget properly
    self.temp_chooser.destroy()
    temp_opts = list(self.temp_dict.keys())
    self.temp_var.set(temp_opts[-1])
    self.temp_chooser = tk.OptionMenu(self.top_frame, self.temp_var, *temp_opts, command=self.temp_chooser_callback)
    self.temp_chooser.grid(row=4, column=3, sticky='EW')
    #Load the last template
    self.temp_chooser_callback()

  def temp_chooser_callback(self, event=None):
    '''
    Fills in the switches specified by the template and refreshes
    '''
    self.switches = self.temp_dict[self.temp_var.get()]
    self.fill_switches()
    self.refresh()

  def col_filter_callback(self, event=None):
    if self.col_filter.get() == "List Attributes":
      self.att_box.configure(state="normal")
      self.col_file.tk_wid.configure(state="disabled")
    elif self.col_filter.get() == "Use Template File":
      self.col_file.tk_wid.configure(state="normal")
      self.att_box.configure(state="disabled")
    else:
      self.att_box.configure(state="disabled")
      self.col_file.tk_wid.configure(state="disabled")

  def ele_list_callback(self, event=None):
    if self.ele_list_opt.get() == "Custom":
      self.ele_list_box.configure(state="normal")
    else:
      self.ele_list_box.configure(state="disabled")

  def open_filter_menu(self, event=None):
    '''
    Opens a menu to set the various row filters
    '''
    win = tk.Toplevel(self)
    win.title("Lattice Row Filters")
    filters = ["Lords only", "No Slaves", "No Super Slaves", "Remove Line if zero", "s"]
    buttons = []
    for i in range(len(self.filter_vars)):
      buttons.append(tk.Checkbutton(win, text=filters[i], variable=self.filter_vars[i]))
      buttons[i].grid(row=i, column=0, sticky='W')
    tk.Entry(win, textvariable=self.remove_if_zero).grid(row=3, column=1, sticky='EW')
    tk.Entry(win, textvariable=self.s_range).grid(row=4, column=1, sticky='EW')

  def toggle_advanced(self, event=None):
    '''
    Toggles self.use_advanced and enables/disables the appropriate widgets
    '''
    widgets = [self.cf_box, self.att_box,
        self.f_button, self.ele_list_chooser,
        self.ele_list_box]
    if self.use_advanced: #Turn off
      for w in widgets:
        w.configure(state="normal")
        # Turn off att_box and ele_list_box, maybe
        self.ele_list_callback()
        self.col_filter_callback()
        self.advanced_box.configure(state="disabled")
    else:
      for w in widgets:
        w.configure(state="disabled")
        self.advanced_box.configure(state="normal")
    self.use_advanced = not self.use_advanced

  def fill_switches(self, event=None):
    '''
    Fills the switch widgets to reflect the values in self.switches
    '''
    # Format of switches is [-option [value]] ... [ele list]
    # Switches with arguments:
    # -att -blank_replacement -branch -custom -remove_line_if_zero -s
    # -universe
    # Switches without arguments:
    # -0undef -all -base -design -floor_coords -lords -middle -no_slaves
    # -no_super_slaves -orbit -radiation_integrals -spin
    # -sum_radiation_integrals -tracking_elements -undef0
    def single_switch(switch, value=''):
      '''
      Fills a single switch widget (switch) with value if appropriate
      '''
      if switch == '-universe':
        self.branch_wids.uni.set(value)
      elif switch == '-branch':
        self.branch_wids.branch.set(value)
      elif switch == '-base':
        self.branch_wids.bmd.set('Base')
      elif switch == '-design':
        self.branch_wids.bmd.set('Design')
      elif switch == '-middle':
        self.branch_wids.bme.set('Middle')
      elif switch == '-floor_coords':
        self.col_filter.set('Floor Coordinates')
      elif switch == '-orbit':
        if self.col_filter.get() == 'Spin':
          self.col_filter.set('Orbit + Spin')
        else:
          self.col_filter.set('Orbit')
      elif switch == '-spin':
        if self.col_filter.get() == 'Orbit':
          self.col_filter.set('Orbit + Spin')
        else:
          self.col_filter.set('Spin')
      elif switch == '-radiation_integrals':
        self.col_filter.set('Radiation Integrals')
      elif switch == '-sum_radiation_integrals':
        self.col_filter.set('Cumulative Radiation Integrals')
      elif switch == '-att':
        self.col_filter.set('List Attributes')
        self.col_atts.set(self.col_atts.get() + value + ' ')
      elif switch == '-custom':
        self.col_filter.set('Use Template File')
        self.col_file.tk_var.set(value)
      elif switch == '-lords':
        self.filter_vars[0].set(True)
      elif switch == '-no_slaves':
        self.filter_vars[1].set(True)
      elif switch == '-no_super_slaves':
        self.filter_vars[2].set(True)
      elif switch == '-remove_line_if_zero':
        self.filter_vars[3].set(True)
        self.remove_if_zero.set(value)
      elif switch == '-s':
        self.filter_vars[4].set(True)
        self.s_range.set(value)
      elif switch == '-all':
        self.ele_list_opt.set('All')
      elif switch == '-tracking_elements':
        self.ele_list_opt.set('Tracking elements')

    #Set all widgets to default state:
    self.branch_wids.uni.set('1')
    self.branch_wids.branch.set('0')
    self.branch_wids.bmd.set('Model')
    self.branch_wids.bme.set('End')
    self.col_filter.set('Default')
    self.col_atts.set('')
    self.col_file.tk_var.set('Browse...')
    for i in range(len(self.filter_vars)):
      self.filter_vars[i].set(False)
    self.remove_if_zero.set('')
    self.s_range.set('')
    self.ele_list_opt.set('Default (first 200)')
    self.ele_list.set('')

    switch_list = self.switches.split(' ')
    arg_switch = False #tracks whether we're looking for a switch argument
    ele_list = False #tracks if we've hit the ele list portion
    for item in switch_list:
      #skip empty items
      if item == '':
        continue
      # Scan for switches/arguments
      if not ele_list:
        if (item[0] == '-') & (not arg_switch): #is a switch, last item was not an arg switch
          if item in ['-att', '-blank_replacement', '-branch', '-custom', '-remove_line_if_zero', '-s', '-universe']:
            arg_switch = item
            continue
          else: #item is a switch, not an arg switch
            single_switch(item)
        elif bool(arg_switch) & (not (item[0] == '-')):
          #item is the argument for arg_switch
          single_switch(arg_switch, item)
          arg_switch = False
        elif (not arg_switch) & (not (item[0] == '-')):
          #item is not a switch or an argument --> ele list started
          ele_list = True
        else: #item is a switch, but we were looking for an argument
          messagebox.showwarning('Warning', 'No argument given for switch "' + arg_switch + '" (one required).  Skipping this switch...')
          arg_switch = False
          continue
      if ele_list:
        self.ele_list.set(self.ele_list.get() + item + ' ')
    #Enable the ele list if it is non_empty
    if self.ele_list.get().strip() != '':
      self.ele_list_opt.set('Custom')
    #Toggle advanced on and off to set widget active/inactive states correctly
    self.toggle_advanced()
    self.toggle_advanced()

  def get_switches(self, event=None):
    '''
    Scans all widgets to get the switches that will be used
    '''
    if self.use_advanced:
      self.switches = self.advanced_var.get()
    else:
      switches = ""
      # Parse uni/branch/base-model-design/beg-middle-end
      switches += "-universe " + self.branch_wids.uni.get() + ' '
      switches += "-branch " + self.branch_wids.branch.get() + ' '
      if self.branch_wids.bmd.get() == "Base":
        switches += '-base '
      elif self.branch_wids.bmd.get() == "Design":
        switches += '-design '
      if self.branch_wids.bme.get() == "Middle":
        switches += '-middle '
      #TODO: beginning switch
      #elif self.branch_wids.bme.get() == "Beginning":
      #  switches += '-beginning '
      # Parse column switches
      if self.col_filter.get() == "Floor Coordinates":
        switches += '-floor_coords '
      elif self.col_filter.get() == "Orbit":
        switches += '-orbit '
      elif self.col_filter.get() == "Spin":
        switches += '-spin '
      elif self.col_filter.get() == "Orbit + Spin":
        switches += '-orbit -spin '
      elif self.col_filter.get() == "Radiation Integrals":
        switches += '-radiation_integrals '
      elif self.col_filter.get() == "Cumulative Radiation Integrals":
        switches += '-sum_radiation_integrals '
      elif self.col_filter.get() == "List Attributes":
        for att in (self.col_atts.get()).split(' '):
          if att != '': #handles more than one space between attributes
            switches += '-att ' + att + ' '
      elif self.col_filter.get() == "Use Template File":
        if self.col_file.tk_var.get() != "Browse...":
          switches += '-custom ' + self.col_file.tk_var.get() + ' '

      # Row switches
      if self.filter_vars[0].get():
        switches += '-lords '
      if self.filter_vars[1].get():
        switches += '-no_slaves '
      if self.filter_vars[2].get():
        switches += '-no_super_slaves '
      if self.filter_vars[3].get():
        switches += '-remove_line_if_zero ' + self.remove_if_zero.get() + ' '
      if self.filter_vars[4].get():
        switches += '-s ' + self.s_range.get() + ' '

      if self.ele_list_opt.get() == "All":
        switches += '-all '
      elif self.ele_list_opt.get() == "Tracking Elements":
        switches += '-tracking_elements '
      elif self.ele_list_opt.get() == "Custom":
        switches += self.ele_list.get()

      self.switches = switches
      self.advanced_var.set(switches)

  def refresh(self, event=None):
    '''
    Fetches the lattice table with the new table
    parameters and updates the window
    '''
    # Clear the existing table
    for child in self.table_frame.winfo_children():
      child.destroy()

    # Get the new table data using the given switches
    self.get_switches()
    lattice = self.pipe.cmd_in("python show lattice -python " + self.switches)
    lattice = lattice.splitlines()
    #Remove error messages
    while lattice[0][:6] in ['[ERROR', '[FATAL']:
      lattice = lattice[1:] #remove line with [ERROR or [FATAL
      while lattice[0].find('    ') == 0:
        lattice = lattice[1:] #remove error description lines
        if len(lattice) == 0:
          break
      if len(lattice) == 0:
        break
    if len(lattice) == 0:
      print("No lattice found")
      tk.Label(self.table_frame, text="NO LATTICE FOUND").pack()
      return
    for i in range(len(lattice)):
      lattice[i] = lattice[i].split(';')
    #lattice[i][j] --> row i, column j
    widths = [0]*len(lattice[0]) # tracks column widths

    # Create table
    self.tree = ttk.Treeview(self.table_frame, columns=lattice[0], show='headings')
    # Column titles
    for title in lattice[0]:
      self.tree.heading(title, text=title)
      self.tree.column(title, stretch=True, anchor='center')

    # Fill rows
    for row in lattice[1:]:
      self.tree.insert("", "end", values=row)
      for j in range(len(row)):
        if len(row[j])*15 > widths[j]:
          widths[j] = len(row[j])*15

    # Set column widths appropriately
    for j in range(len(lattice[0])):
      if len(lattice[0][j])*15 > widths[j]:
        widths[j] = len(lattice[0][j])*15
      self.tree.column(lattice[0][j], width=widths[j], minwidth=widths[j])

    # Scrollbars
    hbar = ttk.Scrollbar(self.table_frame, orient="horizontal", command=self.tree.xview)
    vbar = ttk.Scrollbar(self.table_frame, orient="vertical", command=self.tree.yview)
    self.tree.configure(xscrollcommand=hbar.set, yscrollcommand=vbar.set)

    vbar.pack(side="right", fill="y", expand=0)
    hbar.pack(side="bottom", fill='x', expand=0)
    self.tree.pack(side="left", fill="both", expand=1)
    self.widths = widths

    # Double click to open element
    self.tree.bind('<Double-Button-1>', self.open_ele_window)

    tot = 0
    for w in widths:
      tot = tot+w
    self.maxsize(1800, 1000)
    self.minsize(1300, 100)

  def open_ele_window(self, event=None):
    '''
    Opens an element window for the currently selected row
    '''
    x = self.tree.focus()
    row = self.tree.item(x)
    if isinstance(row['values'][0], int): #checks that the first value is an element index
      settings = [self.branch_wids.uni.get(),
          self.branch_wids.branch.get(),
          str(row['values'][0]),
          self.branch_wids.bmd.get()]
      win = tao_ele_window(self.root, self.pipe, settings)













