import tkinter as tk
import ttk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import font
import sys
import os
import copy
from tao_widget import *
from taoplot import taoplot
from parameters import str_to_tao_param
from elements import *
from tao_set import *
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backends._backend_tk import FigureManagerTk
from matplotlib.backend_bases import key_press_handler
from tao_ele_location import in_element
from tao_mpl_toolbar import taotoolbar
from matplotlib.widgets import Slider

#-----------------------------------------------------
# List window

class tao_list_window(tk.Toplevel):

  def __init__(self, root, title, use_upper=False,
      min_width=0, *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
    self.title(title)
    self.root = root

    #self.wm_geometry(newGeometry='400x600')

    self.upper_frame=tk.Frame(self)
    self.outer_frame=tk.Frame(self)

    self.canvas=tk.Canvas(self.outer_frame)
    self.list_frame=tk.Frame(self.canvas)
    scrollbar=tk.Scrollbar(self.outer_frame,orient="vertical",
        command=self.canvas.yview)
    self.canvas.configure(yscrollcommand=scrollbar.set)

    def scrollhelper(event):
      self.canvas.configure(scrollregion=self.canvas.bbox("all"))
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

    self.canvas_list_window = self.canvas.create_window(
        (0,0),window=self.list_frame,anchor='nw')

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
    self.pipe = pipe
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

    # Grid widgets (and units)
    i = 0
    for item in self.tao_list:
      if item.param.name.find('units#') != 0:
        r = (i % len_col) + 2 # rows 0,1 reserved
        c = int(i / len_col) * 3
        item.tk_label.grid(row=r, column=c, sticky='E')
        item.tk_wid.grid(row=r, column=c+1, sticky='EW')
        i = i+1
      else:
        tk.Label(self,text=item.param.value).grid(row=r, column=c+2, sticky='W')

  def set_params(self, set_str, event=None):
    '''
    Runs tao_set on self.tao_list
    '''
    tao_set(self.tao_list, set_str, self.pipe)

  def check_for_changes(self):
    '''
    Runs check_for_changes (from tao_set) on self.tao_list
    '''
    return check_for_changes(self.tao_list)


#-----------------------------------------------------
# Parameter window
class tao_parameter_window(tao_list_window):

  def __init__(self, root, title, tao_list, pipe, plot="", *args, **kwargs):
    tao_list_window.__init__(self, root, title, *args, **kwargs)
    self.button_frame = tk.Frame(self)
    self.button_frame.pack(side="top", fill="both", expand=0)
    self.tao_list = tao_list
    for k in range(len(self.tao_list)):
      self.tao_list[k] = tk_tao_parameter(
          self.tao_list[k], self.list_frame, pipe, plot=plot)
      tk.Label(self.list_frame,text=self.tao_list[k].param.name).grid(
          row=k,column=0,sticky="E")
      self.tao_list[k].tk_wid.grid(row=k,column=1,sticky="EW")
      k = k+1

#-----------------------------------------------------
# Table window
class table_window(tao_list_window):
  '''
  Meant for showing large amounts of information in a table
  (e.g. d1_data and v1_vars).  Comes with bulk editing, detailed view of
  individual rows, and editing of individual parameters in the table.

  Input parameters:
  root: parent widget
  pipe: the tao_interface object that allows interface with Tao
  array_name: the name of the object this table displays
      (e.g. the name of the d1_datum or v1_variable)
  title_list: the column titles, in order from left to right
  bulk_template: a list with elements [tao_parameter, column], where
      tao_parameter is a generic copy of the parameter that will be bulk filled
      (i.e. initialized to "blank"), and column is the column number it should
      be gridded to (start counting from 0 on the left)
  bulk_set_format: format string for the bulk set string, to be used with the
      str.format() method
  set_format: format string for individual rows' set strings, to be used with
      the str.format() method
      example format strings: "set data 1@{}|", "set data 1@{}[{}]|"
  '''

  def __init__(self, root, pipe, array_name, title_list, bulk_template,
      bulk_set_format, set_format, *args, **kwargs):
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
    b2 = tk.Button(self.button_frame, text="Discard Changes and Refresh",
        command=self.refresh)

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
    tk.Label(self.list_frame, text="Bulk editing:").grid(
        row=1, column=0, columnspan=self.bulk_template[0][1])
    tk.Label(self.list_frame, text="Click to fill:").grid(
        row=2, column=0, columnspan=self.bulk_template[0][1])
    #self.bulk_template[0][1] is the appropriate
    #columnspan because these labels can fill all
    #the space until the first bulk fill box

    self.bulk_params = [] #Holds the current bulk edit widgets
    self.bulk_filled = [] #Holds whether or not bulk filling has been used
    self.bulk_value = [] #Holds the last value that was filled to the cells
    self.bulk_apply = [] #Holds fill buttons

    j=0
    for item in self.bulk_template:
      self.bulk_params.append(tk_tao_parameter(copy.copy(item[0]),
        self.list_frame, self.pipe))
      self.bulk_params[j].tk_wid.grid(row=1, column=item[1])
      self.bulk_params[j].tk_wid.bind("<Return>", self.fill_callback(j))
      self.bulk_filled.append(False)
      self.bulk_value.append(self.bulk_params[j].tk_var.get())
      self.bulk_apply.append(tk.Button(self.list_frame, text="Fill...",
        command = self.fill_callback(j) ))
      self.bulk_apply[j].grid(row=2, column=item[1])
      j=j+1

    #Fetch and fill in the data
    # IT IS EXPECTED THAT SUBCLASSES WILL DEFINE
    # self.row_num AND A METHOD FOR POPULATING
    # self.list_rows, AND THEN CALL THIS METHOD
    #i = row counter, j = column counter
    #grid to row i+3 because row 0 is titles, row 1 is bulk editing widgets,
    #row 2 is fill buttons
    self.list_rows = []
    for i in range(self.row_num):
      self.list_rows.append(self.list_row_fetch(i))
      for j in range(len(self.list_rows[i].tk_wids)):
        self.list_rows[i].tk_wids[j].grid(row=i+3, column=j)
      tk.Button(self.list_frame, text="View More...",
          command=self.open_detail_window_callback(
            self.list_rows[i].index)).grid(
                row=i+3, column=len(self.list_rows[i].tk_wids))

  def fill(self, index, event=None):
    '''
    Fills the column specified by bulk_template to the bulk edit value,
    saves the bulk edit value for efficient calls to tao_set,
    and clears the bulk edit box
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
      self.list_rows[i].tk_tao_params[
          self.bulk_template[index][0].name].tk_var.set(self.bulk_value[index])

  def fill_callback(self, index, event=None):
    return lambda event=None : self.fill(index)


  def apply(self):
    #Apply bulk changes
    for i in range(len(self.bulk_params)):
      if self.bulk_filled[i]:
        set_str = self.bulk_set_format.format(self.array_name)
        self.bulk_params[i].tk_var.set(self.bulk_value[i])
        #overide is necessary for LOGIC parameters
        tao_set([self.bulk_params[i]], set_str, self.pipe, overide=(
          self.bulk_params[i].param.type=='LOGIC'))

    #Apply individual changes that are different from bulk changes
    for i in range(len(self.list_rows)):
      set_list = []
      #NOTE: it is expected that self.list_rows[i].index exist,
      #and be equal to that row's index
      set_str = self.set_format.format(self.array_name, self.list_rows[i].index)

      #Find elements in row that need setting
      for j in range(len(self.bulk_template)):
        name = self.bulk_template[j][0].name
        c1 = (self.list_rows[i].tk_tao_params[name].tk_var.get()
            != self.bulk_value[j])
        c2 = not self.bulk_filled[j]
        try:
          if self.bulk_template[j][0].type == 'REAL':
            c3 = (float(self.list_rows[i].tk_tao_params[name].tk_var.get())
                != self.list_rows[i].tk_tao_params[name].param.value)
          elif self.bulk_template[j][0].type == 'INT':
            c3 = (int(self.list_rows[i].tk_tao_params[name].tk_var.get())
                != self.list_rows[i].tk_tao_params[name].param.value)
          elif self.bulk_template[j][0].type == 'STR':
            c3 = (str(self.list_rows[i].tk_tao_params[name].tk_var.get())
                != self.list_rows[i].tk_tao_params[name].param.value)
          elif self.bulk_template[j][0].type == 'LOGIC':
            c3 = (bool(self.list_rows[i].tk_tao_params[name].tk_var.get())
                != self.list_rows[i].tk_tao_params[name].param.value)
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
    b = tk.Button(win.button_frame, text="Apply changes",
        command=lambda : self.detail_set_callback(win.tao_list,set_str))
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
  def __init__(self, root, pipe, array_name, title_list, bulk_format,
      bulk_set_format, set_format, *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
    self.title(array_name)
    self.pipe = pipe
    self.array_name = array_name
    self.title_list = title_list
    self.bulk_format = bulk_format
    self.bulk_set_format = bulk_set_format
    self.set_format = set_format
    self.button_frame = tk.Frame(self) #holds the buttons
    self.button_frame.pack(fill='x', expand=0)
    self.table_frame = tk.Frame(self) #holds the table
    self.table_frame.pack(fill='both', expand=1)

    if len(bulk_format) > 0:
      tk.Button(self.button_frame, text="Bulk fill...",
          command=self.open_bulk_window).pack(side='left')

    tk.Button(self.button_frame, text="Show all",
        command=self.show_all).pack(side='right')
    self.hide_nonexist = True
    self.refresh()

  def show_all(self, event=None):
    self.hide_nonexist = False
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
    self.tree = ttk.Treeview(
        self.table_frame, columns=self.title_list, show='headings')
    # Column titles
    for title in self.title_list:
      self.tree.heading(title, text=title)
      self.tree.column(title, stretch=True, anchor='center')

    # Fill rows
    # Fetch and fill in the data
    # IT IS EXPECTED THAT SUBCLASSES WILL DEFINE
    # self.row_num AND A METHOD FOR POPULATING
    # self.list_rows, AND THEN CALL THIS METHOD
    #i = row counter, j = column counter
    self.list_rows = []
    for i in range(self.row_num):
      row = self.lw_list_row_fetch(i)
      if self.hide_nonexist:
        # Filter out if exists is False
        if row[-1] == 'F':
          continue
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
    hbar = ttk.Scrollbar(
        self.table_frame, orient="horizontal", command=self.tree.xview)
    vbar = ttk.Scrollbar(
        self.table_frame, orient="vertical", command=self.tree.yview)
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

  def open_bulk_window(self, event=None):
    '''
    Opens a window with bulk settings for meas_value,
    ref_value, weight, and good_user.
    '''
    win = tk.Toplevel(self)
    win.title("Bulk settings for " + self.array_name)

    j = 0 #column counter
    fill_choices = [] #what is being used to fill (tk.StringVar()'s)
    for item in self.bulk_format:
      fill_choices.append(j)
      fill_frame, where_frame, fill_choices[j] = self.make_bulk_frame(
          win, item[0], usebmd=(item[0].name=='meas_value'))
      fill_frame.grid(row=0, column=2*j, sticky='NSEW')
      where_frame.grid(row=2, column=2*j, sticky='NSEW')
      if j != 0:
        ttk.Separator(win, orient='vertical').grid(
            row=0, column=2*j-1, rowspan=3, sticky='NS')
      win.columnconfigure(j, weight=1)
      j = j+1
    ttk.Separator(win, orient='horizontal').grid(
        row=1, column=0, columnspan=2*j-1, sticky='EW')
    ttk.Separator(win, orient='horizontal').grid(
        row=3, column=0, columnspan=2*j-1, sticky='EW')

    def fill_cmd(event=None):
      self.bulk_set(fill_choices, win)

    tk.Button(win, text="Fill and apply", command=fill_cmd).grid(
        row=4, column=0, columnspan=2*j-1)

  def bulk_set(self, fill_choices, parent):
    '''
    Runs set commands to set variables specified in a bulk window
    appropriately.
    '''
    for i in range(len(self.bulk_format)):
      if fill_choices[i]['choice'].get() == 'none':
        continue
      else:
        fill_val = fill_choices[i][fill_choices[i]['choice'].get()].get()

      # Check that the formula gives the correct length array
      if fill_choices[i]['choice'].get() == 'formula':
        if fill_choices[i]['where'].get() == 'all':
          len1 = len(self.pipe.cmd_in('python evaluate '
            + self.array_name + '|'
            + self.bulk_format[i][0].name).splitlines())
        elif fill_choices[i]['where'].get() == 'range':
          len1 = len(self.pipe.cmd_in('python evaluate '
            + self.array_name + '[' + fill_choices[i]['range'].get()
            + ']|' + self.bulk_format[i][0].name).splitlines())
        len2 = len(self.pipe.cmd_in('python evaluate '
          + fill_val).splitlines())
        if len1 != len2:
          messagebox.showwarning('Error',
              'Entered formula yields an array of incorrect length.',
              parent=parent)
          continue

      if fill_choices[i]['where'].get() == 'all':
        set_str = self.bulk_set_format.format(self.array_name)
      elif fill_choices[i]['where'].get() == 'range':
        set_str = self.bulk_set_format.format(self.array_name + '['
            + fill_choices[i]['range'].get() + ']')

      if fill_choices[i]['choice'].get() == 'const':
        # Make sure fill_val is True/False for LOGIC parameters
        if self.bulk_format[i][0].type == 'LOGIC':
          fill_val = str(bool(fill_val))
        self.pipe.cmd_in(set_str + self.bulk_format[i][0].name
            + ' = ' + fill_val)
      elif fill_choices[i]['choice'].get() == 'bmd':
        fill_val = fill_val.lower()
        self.pipe.cmd_in(set_str + self.bulk_format[i][0].name
            + ' = ' + set_str.split(' ')[-1] + fill_val)
      elif fill_choices[i]['choice'].get() == 'formula':
        self.pipe.cmd_in(set_str + self.bulk_format[i][0].name
            + ' = ' + fill_val)
    self.refresh()


  def make_bulk_frame(self, parent, bulk_item, usebmd=False):
    '''
    Creates a frame with all the widgets needed for one parameter in the
    bulk settings window.
    parent: parent widget for this frame
    bulk_item: tao_parameter instance for this frame's parameter
    usebmd: if set true, and option will be given to set to the
        base/model/design/meas/ref value
    Returns a tuple (fill_frame, where_frame, fill_vars)
    '''
    fill_frame = tk.Frame(parent)
    fill_frame.columnconfigure(2, weight=1)
    # Title the column
    tk.Label(fill_frame, text=bulk_item.name).grid(row=0, column=0,
        columnspan=3, sticky='EW')
    # What to fill with
    fill_choice = tk.StringVar()
    fill_choice.set('none')
    fill_vars = {}
    fill_vars['choice'] = fill_choice

    # No change
    none_button = tk.Radiobutton(fill_frame, text="",
        variable=fill_choice, value='none')
    none_button.grid(row=1, column=0, sticky='W')
    tk.Label(fill_frame, text="No Change").grid(
        row=1, column=1, sticky='W')
    # Constant
    const_button = tk.Radiobutton(fill_frame, text="",
        variable=fill_choice, value='const')
    const_button.grid(row=2, column=0, sticky='W')
    tk.Label(fill_frame, text="Constant:").grid(
        row=2, column=1, sticky='W')
    if bulk_item.type == 'REAL':
      const_var = tk.StringVar()
      const_wid = tk.Entry(fill_frame, textvariable=const_var)
      const_wid.grid(row=2, column=2, sticky='EW')
    elif bulk_item.type == 'LOGIC':
      const_var = tk.BooleanVar()
      const_wid = tk.Checkbutton(fill_frame, variable=const_var)
      const_wid.grid(row=2, column=2, sticky='W')
    fill_vars['const'] = const_var
    # Base/Model/Design/Meas/Ref
    if usebmd:
      bmd_var = tk.StringVar()
      bmd_var.set('Base')
      bmd_button = tk.Radiobutton(fill_frame, text="",
          variable=fill_choice, value='bmd')
      bmd_button.grid(row=3, column=0, sticky='W')
      tk.Label(fill_frame, text="From:").grid(
          row=3, column=1, sticky='W')
      bmd_wid = tk.OptionMenu(fill_frame, bmd_var, "Base", "Model",
          "Design", "Measure", "Reference")
      bmd_wid.grid(row=3, column=2, sticky='EW')
      fill_vars['bmd'] = bmd_var
    # Formula
    if bulk_item.type == 'REAL':
      formula_button = tk.Radiobutton(fill_frame, text="",
          variable=fill_choice, value='formula')
      formula_button.grid(row=4, column=0, sticky='W')
      tk.Label(fill_frame, text="Formula:").grid(
          row=4, column=1, sticky='W')
      formula_var = tk.StringVar()
      formula_wid = tk.Entry(fill_frame, textvariable=formula_var)
      formula_wid.grid(row=4, column=2)
      fill_vars['formula'] = formula_var

    # Where to fill
    where_frame = tk.Frame(parent)
    fill_where = tk.StringVar()
    fill_where.set('all')
    fill_vars['where'] = fill_where
    fill_range = tk.StringVar()
    # All
    all_button = tk.Radiobutton(where_frame, text="",
        variable=fill_where, value='all')
    all_button.grid(row=0, column=0, sticky='W')
    tk.Label(where_frame, text="All").grid(row=0, column=1, sticky='W')
    # Range
    range_button = tk.Radiobutton(where_frame, text="",
        variable=fill_where, value='range')
    range_button.grid(row=1, column=0, sticky='W')
    tk.Label(where_frame, text="Range:").grid(
        row=1, column=1, sticky='W')
    range_wid = tk.Entry(where_frame, textvariable=fill_range)
    range_wid.grid(row=1, column=2, sticky='EW')
    fill_vars['range'] = fill_range

    return (fill_frame, where_frame, fill_vars)


  def open_detail_window(self, index):
    '''
    Opens up a detail window for the given index.
    '''
    detail_title = self.array_name + '[' + str(index) + ']'
    win = tao_parameter_window(self, detail_title, self.param_list, self.pipe)

    set_str = self.set_format.format(self.array_name, index)
    b = tk.Button(win.button_frame, text="Apply changes",
        command=lambda : self.detail_set_callback(win.tao_list,set_str))
    b.pack()

  def detail_set_callback(self, tao_list, set_str):
    tao_set(tao_list, set_str, self.pipe)
    # set exists for the array
    # SUBCLASSES MUST DEFINE self.set_exists
    self.pipe.cmd_in(self.set_exists)
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
    u_ix_box = tk.OptionMenu(
        self.univ_frame, self.u_ix, *u_ix_list, command=self.refresh)
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
      new_frame = d2_data_frame(
          self.list_frame, self.root, self.pipe, d2_data_item, u_ix)
      new_frame.frame.pack()


#-----------------------------------------------------
# d1_data window

class tao_d1_data_window(lw_table_window):
  '''
  With lw set to True, opens a lw_table_window instead
  '''
  def __init__(self, root, pipe, d1_data_name,
      u_ix, ix_lb, ix_ub, *args, **kwargs):
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
        "Weight",
        "Exists"]
    bulk_template = []
    bulk_template.append([str_to_tao_param("meas_value;REAL;T;"), 6])
    bulk_template.append([str_to_tao_param("good_user;LOGIC;T;"), 11])
    bulk_template.append([str_to_tao_param("weight;REAL;T;"), 12])

    bulk_set_format = "set data " + str(self.u_ix) + '@{}|'
    set_format = "set data " + str(self.u_ix) + '@{}[{}]|'
    self.set_exists = "set data " + str(self.u_ix) + '@' + d1_data_name + "|exists = T"
    if self.lw:
      lw_table_window.__init__(self, root, pipe, d1_data_name, title_list,
          bulk_template, bulk_set_format, set_format, *args, **kwargs)
    else:
      table_window.__init__(self, root, pipe, d1_data_name, title_list,
          bulk_template, bulk_set_format, set_format, *args, **kwargs)

  def refresh(self):
    self.d_list = self.pipe.cmd_in(
        "python data_d_array " + self.u_ix + '@' + self.array_name)
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
    self.param_list = self.pipe.cmd_in("python data " + str(self.u_ix) + '@'
        + self.array_name + '[' + str(index) + ']')
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

    tk.Label(self.list_frame, text="Variable").grid(
        row=0, column=0, columnspan=2)
    tk.Label(self.list_frame, text="Indices").grid(row=0, column=2)
    tk.Label(self.list_frame, text="Using").grid(row=0, column=3)

    i=1
    for item in v1_var_list:
      tk.Label(self.list_frame, text=item[0]).grid(row=i, column=0)
      tk.Button(self.list_frame, text="View...",
          command=self.open_v1_callback(item[0])).grid(row=i, column=1)
      tk.Label(self.list_frame,text=item[2] +':'+ item[3]).grid(row=i,column=2)
      tk.Label(self.list_frame, text=item[1]).grid(row=i, column=3)
      i = i+1

  def open_v1_callback(self, v1_var_name):
    return lambda : self.open_v1(v1_var_name)

  def open_v1(self, v1_var_name):
    win = tao_v1_var_window(self.root, self.pipe, v1_var_name)

#-----------------------------------------------------
# v1_var window

class tao_v1_var_window(lw_table_window):

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
    self.set_exists = "set var " + v1_var_name + "|exists = T"

    bulk_set_format = "set var {}|"
    set_format = "set var {}[{}]|"
    if self.lw:
      lw_table_window.__init__(self, root, pipe, v1_var_name, title_list,
          bulk_template, bulk_set_format, set_format, *args, **kwargs)
    else:
      table_window.__init__(self, root, pipe, v1_var_name, title_list,
          bulk_template, bulk_set_format, set_format, *args, **kwargs)

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
    self.param_list = self.pipe.cmd_in(
        "python var " + self.array_name + '[' + str(index) + ']')
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
      self.root.console.set_command(cmd_string)
      self.root.console.run_command()
      #self.root.command.tk_var.set(cmd_string)
      #self.root.tao_command()
    elif mode == 1:
      self.root.call_file.tk_var.set(cmd_string)
      self.root.tao_call()
    elif mode ==2:
      self.root.console.set_command(cmd_string)
      #self.root.command.tk_var.set(cmd_string)
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
    self.temp_select = ttk.Combobox(
        self.temp_frame, textvariable=self.temp, values=self.plot_list)
    self.temp_select.bind("<<ComboboxSelected>>", self.refresh)
    self.temp_select.bind("<Return>", self.refresh)
    self.temp_select.grid(row=0, column=1)

    self.refresh()

    self.button_frame = tk.Frame(self)
    self.button_frame.pack(fill="both", expand=0)
    b = tk.Button(
        self.button_frame, text="Apply changes", command=self.plot_apply)
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
    Clears self.list_frame and populates it with information relevant to
    self.temp
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
    name = tk_tao_parameter(
        str_to_tao_param(data_list.pop(0)), self.list_frame, self.pipe)
    name.tk_label.grid(row=0, column=0)
    name.tk_wid.grid(row=0, column=1, sticky='EW')
    tk.Label(self.list_frame, text="Graphs").grid(
        row=1, column=0, rowspan=num_graphs)
    i=1
    for graph in self.graph_list:
      tk.Button(self.list_frame, text=graph,
          command=self.open_graph_callback(name.param.value, graph)).grid(
              row=i, column=1, sticky='EW')
      i = i+1

    # Grid the rest of the information
    self.list_frame.columnconfigure(1, weight=1)
    for i in range(len(data_list)):
      data_list[i] = tk_tao_parameter(
          str_to_tao_param(data_list[i]), self.list_frame, self.pipe)
      # Make sure the Entry boxes are big enough to
      # display their contents
      if data_list[i].param.type in ['STR', 'INT', 'REAL']:
        data_list[i].tk_wid.configure(
            width=len(str(data_list[i].param.value))+1)
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

    # Check if the template has been placed in a region already,
    # and place it if necessary
    if self.template in self.root.placed.keys():
      print("found in root.placed")
      pass
    else:
      #Place the plot in the next available region
      #and make visible
      r_index = 1
      while (("r" + str(r_index)) in self.root.placed.values()):
        r_index = r_index + 1
      self.pipe.cmd_in("place -no_buffer r" + str(r_index) + " " + self.template)
      self.pipe.cmd_in("set plot r" + str(r_index) + ' visible = T')
      self.root.placed[self.template] = 'r' + str(r_index)

    self.mpl = taoplot(pipe, self.root.placed[self.template])
    self.refresh()

  def refresh(self, event=None, width=1):
    '''
    Makes the call to matplotlib to draw the plot to the window
    '''
    #Clear the window
    for child in self.winfo_children():
      child.destroy()

    #Get plotting results
    self.plot_output = self.mpl.plot(width)

    #Get the figure
    self.fig = self.plot_output[0]

    #Get figure information
    self.fig_info = self.plot_output[1]

    #Create widgets to display the figure
    canvas = FigureCanvasTkAgg(self.fig, master=self)
    canvas.draw()
    canvas.get_tk_widget().pack(side="top", fill="both", expand=1)
    # DO NOT TOUCH
    canvas.manager = FigureManagerTk(
        canvas, self.fig.number, tk.Toplevel(self.root))

    toolbar = taotoolbar(canvas, self)
    toolbar.update()
    canvas._tkcanvas.pack(side="top", fill="both", expand=1)

    def on_key_press(event):
      key_press_handler(event, canvas, toolbar)

    canvas.mpl_connect("key_press_event", on_key_press)

    def on_click(event):
      if event.dblclick:
        eleList = in_element(event.xdata,event.ydata,self.fig_info)
        for i in eleList:
          tao_ele_window(self.root,self.pipe,
              default=[self.fig_info[1],i[0],i[1],self.fig_info[3]])

    canvas.mpl_connect("button_press_event", on_click)

    if self.fig_info[0] == 'floor_plan':
      self.fig.subplots_adjust(bottom=0.2) #adds room below graph for slider
      width_slider = Slider(self.fig.add_axes([.1,.05,.8,.05]), 'width', 0, 2, width) #element width slider

      def update_slider(width):
        self.refresh(width=width_slider.val)

      width_slider.on_changed(update_slider) #call update when slider moves

  def destroy(self):
    # Note: lat_layout should not be automatically removed from r1
    if self.template != "lat_layout":
      # Unplace the template from its region
      self.pipe.cmd_in("place -no_buffer " + self.root.placed[self.template] + " none")
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
    b = tk.Button(
        self.button_frame, text="Apply changes", command=self.graph_apply)
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
    name = tk_tao_parameter(
        str_to_tao_param(data_list.pop(0)), self.list_frame, self.pipe)
    name.tk_label.grid(row=0, column=0)
    name.tk_wid.grid(row=0, column=1, sticky='EW')

    # Curve buttons
    if num_curves > 0:
      tk.Label(self.list_frame, text="Curves").grid(
          row=1, column=0, rowspan=num_curves)
      i=1
      for curve in curve_list:
        tk.Button(self.list_frame, text=curve,
            command=self.open_curve_callback(self.graph, curve)).grid(
                row=i, column=1, sticky='EW')
        i = i+1

    # Grid everything else
    self.list_frame.columnconfigure(1, weight=1)
    for i in range(len(data_list)):
      data_list[i] = tk_tao_parameter(
          str_to_tao_param(data_list[i]), self.list_frame, self.pipe)
    self.tao_list = data_list
    for i in range(len(data_list)):
      # Make sure the Entry boxes are big enough to
      # display their contents
      if data_list[i].param.type in ['STR', 'INT', 'REAL']:
        data_list[i].tk_wid.configure(
            width=len(str(data_list[i].param.value))+1)
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

    b = tk.Button(win.button_frame, text="Apply changes",
        command=lambda : self.curve_apply(win))
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

    tao_parameter_window.__init__(
        self, root, curve, data_list, self.pipe, plot=curve.split('.')[0], *args, **kwargs)

#-----------------------------------------------------
# Branch/element choosing widgets
class tao_branch_widgets:
  '''
  Provides several widgets for slecting universe,
  branch, and elements
  Available widgets:
  self.uni_chooser: OptionMenu to pick the universe index
  self.branch_chooser: OptionMenu to pick the branch (displays name and index)
  self.ele_chooser: ttk Combobox for selecting an element.
      May be specified by name or index
  self.bmd_chooser: OptionMenu for choosing base, model, or design
  self.bme_chooser: OptionMenu for choosing beginning, middle, or end
  In addition, the class provides tk variables for each of these widgets,
      named the same but without _chooser at the end
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
    self.b_list = {} # b_list[i] = branches in universe i
    self.b_name_list = {} # b_name_list[i] = branch names in universe i
    self.e_list = {} #e_list[i][j] = range from 0 to the max element number
    self.e_name_list = {} #e_name_list[i][j] = ele names in branch j of uni i
    self.e_display_list = {} # indexed the same as e_list and e_name_list,
                             # but for displaying the number and name
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
        ele_names = self.pipe.cmd_in(
            "python lat_ele_list " + u + '@' + branch_num)
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
    self.branch_name.trace('w', self.update_branch)
    try: # in case self.ele was set by name and not by index
      ele_num = int(self.ele.get())
      self.ele.set(
          self.e_display_list[self.uni.get()][self.branch.get()][ele_num])
    except ValueError:
      pass
    self.ele_label.set("Element (0 to "
        + str(self.e_list[self.uni.get()][self.branch.get()][-1]) + ")")

    self.uni_chooser = tk.OptionMenu(
        self.parent, self.uni, *self.u_list, command=self.make_branch)
    self.branch_chooser = tk.OptionMenu(self.parent, self.branch_name,
        *self.b_name_list[self.uni.get()], command=self.make_ele_chooser)
    self.ele_chooser = ttk.Combobox(self.parent, textvariable=self.ele,
        values=self.e_display_list[self.uni.get()][self.branch.get()])
    #self.ele_chooser.bind("<<ComboboxSelected>>", self.refresh)
    #self.ele_chooser.bind("<Return>", self.refresh)
    self.bmd_chooser = tk.OptionMenu(
        self.parent, self.bmd, "Base", "Model", "Design")
    self.bme_chooser = tk.OptionMenu(
        self.parent, self.bme, "Beginning", "Middle", "End")
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
    self.branch_chooser = tk.OptionMenu(self.parent, self.branch_name,
        *self.b_name_list[self.uni.get()], command=self.make_ele_chooser)
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

    # Update self.branch to match self.branch_name
    ix = self.b_name_list[self.uni.get()].index(self.branch_name.get())
    self.branch.set(self.b_list[self.uni.get()][ix])

    #self.ele_chooser.destroy()
    #self.ele_chooser = ttk.Combobox(self.parent, textvariable=self.ele,
    self.ele_chooser.configure(values=self.e_display_list[self.uni.get()][self.branch.get()])
    #self.ele_chooser.bind("<<ComboboxSelected>>", self.update)
    #self.ele_chooser.bind("<Return>", self.update)
    self.ele_label.set("Element (0 to "
        + str(self.e_list[self.uni.get()][self.branch.get()][-1]) + ")")
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
    if (self.ele.get() not in
        self.e_name_list[self.uni.get()][self.branch.get()]) \
        & (self.ele.get() not in
            self.e_display_list[self.uni.get()][self.branch.get()]):
      try:
        if (int(self.ele.get()) not in
            self.e_list[self.uni.get()][self.branch.get()]):
          messagebox.showwarning("Error", "Element not found", parent=self.parent)
          return 0
      except ValueError:
        messagebox.showwarning("Error", "Element not found", parent=self.parent)
        return 0

    # Set self.ele, in case the element
    # was specified by name or display name
    if self.ele.get() in self.e_name_list[self.uni.get()][self.branch.get()]:
      ele_num = self.e_name_list[self.uni.get()][self.branch.get()].index(
          self.ele.get())
      self.ele.set(str(ele_num))
    elif self.ele.get() in self.e_display_list[self.uni.get()][self.branch.get()]:
      ele_num = self.e_display_list[self.uni.get()][self.branch.get()].index(
          self.ele.get())
      self.ele.set(str(ele_num))
    return 1

  def update_branch(self, *args):
    '''
    Trace callback for self.branch_name to update self.branch
    '''
    ix = self.b_name_list[self.uni.get()].index(self.branch_name.get())
    self.branch.set(self.b_list[self.uni.get()][ix])


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
    tao_list_window.__init__(self, root, "Lattice Elements", use_upper=True,
        min_width=600, *args, **kwargs)
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
    tk.Label(self.top_frame, textvariable=self.ele_wids.ele_label).grid(row=0, column=2)
    tk.Label(self.top_frame, text="Base/Model/Design").grid(row=0, column=3)

    # Configure and place widgets
    self.ele_wids.ele_chooser.bind("<<ComboboxSelected>>", self.refresh)
    self.ele_wids.ele_chooser.bind("<Return>", self.refresh)
    self.ele_wids.bmd.trace('w', self.refresh)
    self.ele_wids.uni_chooser.grid(row=1, column=0, sticky='EW')
    self.ele_wids.branch_chooser.grid(row=1, column=1, sticky='EW')
    self.ele_wids.ele_chooser.grid(row=1, column=2, sticky='EW')
    self.ele_wids.bmd_chooser.grid(row=1, column=3, sticky='EW')

    self.refresh()

  def refresh(self, event=None, *args):
    '''
    This is where most of the element information is actually created
    '''
    # Ask to save changes
    something_changed = False
    try:
      something_changed = check_for_changes(self.head_tk_tao_params)
    except:
      pass
    try:
      for i in range(len(self.p_frames)):
        if self.p_names[i] not in ['lord_slave', 'mat6', 'floor']:
          something_changed = something_changed | \
              self.p_frames[i].check_for_changes()
    except:
      pass
    if something_changed:
      x = messagebox.askyesnocancel(title="Unsaved Changes",
          message="Apply changes before switching elements?", parent=self)
      if x:
        self.ele_set()
      if x == None:
        return #don't refresh if "Cancel is picked"
    # Update the ele_wids and check for successful update
    if not self.ele_wids.update():
      messagebox.showwarning("Error", "Element not found", parent=self)
      return

    # Clear existing window contents
    for child in self.head_frame.winfo_children():
      child.destroy()
    for child in self.list_frame.winfo_children():
      child.destroy()

    # Create an element object
    self.element = lat_element(self.ele_wids.uni.get(),
        self.ele_wids.branch.get(), self.ele_wids.ele.get(),
        (self.ele_wids.bmd.get()).lower(), self.pipe)

    # Populate self.head_frame
    self.head_tk_tao_params = []
    for p in self.element.params.keys():
      self.element.params[p] = tk_tao_parameter(
          self.element.params[p], self.head_frame, self.pipe)
      if self.element.params[p].param.can_vary:
        self.head_tk_tao_params.append(self.element.params[p])
    name = tk.Label(
        self.head_frame, text=self.element.params["name"].param.value)
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
    tk.Label(self.head_frame, text="Description").grid(
        row=3, column=2, sticky='E')
    self.element.params["descrip"].tk_wid.grid(row=3, column=3, sticky='EW')
    tk.Label(self.head_frame, text="is_on").grid(row=4, column=2, sticky='E')
    self.element.params["is_on"].tk_wid.grid(row=4, column=3, sticky='EW')
    # Set button
    tk.Button(self.top_frame, text='Apply all\nchanges',
        command=self.ele_set).grid(row=0, column=4, rowspan=2, sticky='EW')

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
        if key in ['ab_multipoles', 'kt_multipoles']:
          key = "multipoles"
        if key in ["multipoles", "elec_multipoles"]:
          self.sh_b_list.append(tk.Button(self.list_frame, text=key))
          self.tao_lists.append([]) #FIX THIS
          # set up shared variables for multipoles_on and scale_multipoles
          try:
            x = self.list_frame.mp_on_var.get()
          except:
            self.list_frame.mp_on_var = tk.BooleanVar()
          try:
            x = self.list_frame.scale_mp_var.get()
          except:
            self.list_frame.scale_mp_var = tk.BooleanVar()
          tao_output = self.pipe.cmd_in(
              "python ele:" + key + ' ' + self.element.id)
          self.p_frames.append(
              tao_multipole_frame(self.list_frame, tao_output, self.pipe))
          self.p_names.append(key)
          self.sh_b_list[i].configure(command=self.s_callback(i))
          self.sh_b_list[i].grid(row=2*i, column=0, sticky='W')
          i = i+1
          continue
        if key == "lord_slave": # extremely special case
          self.sh_b_list.append(tk.Button(self.list_frame, text=key))
          ls_frame = tk.Frame(self.list_frame)
          ls_list = self.pipe.cmd_in("python ele:lord_slave "
              + self.ele_wids.uni.get() + '@' + self.ele_wids.branch.get()
              + '>>' + self.ele_wids.ele.get() + '|'
              + (self.ele_wids.bmd.get()).lower())
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
              current_level = ls_tree.insert(
                  "", "end", text=line[0], values=line[1:])
            else:
              ls_tree.insert(
                  current_level, "end", text=line[0], values=line[1:])
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
        if key == 'mat6':
          #tao_list = '\n' + self.pipe.cmd_in("python ele:mat6 "
          #    + self.element.id + ' err')
          tao_list = self.pipe.cmd_in("python ele:mat6 "
              + self.element.id + ' mat6')
          #tao_list += '\n' + self.pipe.cmd_in("python ele:mat6 "
          #    + self.element.id + ' vec0')
        else:
          tao_list = self.pipe.cmd_in("python ele:" + key
              + ' ' + self.element.id)
        self.tao_lists.append(tao_list.splitlines())
        for j in range(len(self.tao_lists[i])):
          self.tao_lists[i][j] = str_to_tao_param(self.tao_lists[i][j])
        # Configure the buttons with commands
        self.sh_b_list[i].configure(command=self.s_callback(i))
        if key in ['mat6', 'floor']: #should just have one column
          n_cols = 1
        else:
          n_cols = 2
        self.p_frames.append(tao_parameter_frame(
            self.list_frame, self.tao_lists[i], n_cols, self.pipe))
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
        if key == 'mat6':
          # add symplectic error
          sym_err = self.pipe.cmd_in('python ele:mat6 '
              + self.element.id + ' err')
          sym_err = tk_tao_parameter(str_to_tao_param(sym_err),
              self.p_frames[i], self.pipe)
          sym_err.tk_label.grid(row=0, column=0, sticky='E')
          sym_err.tk_wid.grid(row=0, column=1, sticky='W')
          # add vec0
          vec0 = self.pipe.cmd_in('python ele:mat6 '
              + self.element.id + ' vec0')
          vec0 = tk_tao_parameter(str_to_tao_param(vec0),
              self.p_frames[i], self.pipe)
          separator = tk.Frame(self.p_frames[i], height=10, bd=1).grid(
              row=8, column=1, columnspan=6, sticky='EW')
          vec0.tk_label.grid(row=9, column=0, sticky='E')
          vec0.tk_wid.grid(row=9, column=1, sticky='W')
          # labels
          mat6_list = ['x', 'px', 'y', 'py', 'z', 'pz']
          title_frame = tk.Frame(self.p_frames[i])
          for j in range(1,7):
            self.p_frames[i].tao_list[j-1].tk_label.grid_forget()
            tk.Label(title_frame, text=mat6_list[j-1]).grid(
                row=1, column=j, sticky='EW')
            title_frame.columnconfigure(j, weight=1)
            tk.Label(self.p_frames[i], text=mat6_list[j-1]).grid(
                row=j+1, column=0, sticky='E')
          tk.Label(self.p_frames[i], text='mat6').grid(
              row=1, column=0, sticky='EW')
          title_frame.grid(row=1, column=1, sticky='EW')
        self.sh_b_list[i].grid(row=2*i, column=0, sticky='W')
        #self.p_frames[i].pack()
        i = i+1

    # Reset self.ele to be the display name
    u = self.ele_wids.uni.get()
    b = self.ele_wids.branch.get()
    e = self.ele_wids.ele.get()
    self.ele_wids.ele.set(self.ele_wids.e_display_list[u][b][int(e)])
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

  def ele_set(self, event=None):
    '''
    Runs set commands for all the parameter frames and refreshes the window
    '''
    # Don't need |base/model/design when setting
    set_str = "set element " + self.element.id[:self.element.id.find('|')] + ' '
    # Set the head parameters
    tao_set(self.head_tk_tao_params, set_str, self.pipe)
    # Set the parameters in self.p_frames
    for i in range(len(self.p_frames)):
      if self.p_names[i] not in ['lord_slave', 'mat6', 'floor']:
        self.p_frames[i].set_params(set_str)
      elif self.p_names[i] == 'floor':
        for p in self.p_frames[i].tao_list:
          if p.param.can_vary:
            #need to set x,y,z,theta,phi,psi_position
            #to the values in p
            floor_list = []
            names = ['x_position',
                'y_position',
                'z_position',
                'theta_position',
                'phi_position',
                'psi_position']
            for i in range(len(p._svar)):
              floor_list.append(
                  str_to_tao_param(names[i]+';REAL;T;' + p._svar[i].get()))
              floor_list[i] = tk_tao_parameter(
                  floor_list[i], self.head_frame, self.pipe)
            # Run the set command
            tao_set(floor_list, set_str, self.pipe, overide=True)
    # Refresh the element window
    self.refresh()

#---------------------------------------------------
class tao_multipole_frame(tk.Frame):
  '''
  Displays multipole information (output of a
  python ele:multipoles or ele:elec_multipoles
  command in a table with certain elements editable
  parent: the parent widget of this frame
  tao_output: the raw, unfiltered output of the tao command that gave
      the multipole information
  pipe: tao_interface object
  '''
  def __init__(self, parent, tao_output, pipe, *args, **kwargs):
    tk.Frame.__init__(self, parent, *args, **kwargs)
    self.pipe = pipe
    self.top_frame = tk.Frame(self) #holds multipoles_on and scale_multipoles
    self.top_frame.pack(fill="both", expand=1)
    self.table_frame = tk.Frame(self) #holds multipole information
    self.table_frame.pack(fill="both", expand=1)
    self.button_frame = tk.Frame(self) #holds the button
    self.button_frame.pack(fill="both", expand=0)
    self.tk_tao_list = []
    self.shown_orders = [] #stores which multipole orders are shown
    self.titles = [] #stores column titles

    # Rip off first two lines of tao_output
    tao_output = tao_output.splitlines()
    self.top_info = tao_output[:2]
    # Display in top_frame
    for i in range(len(self.top_info)):
      self.top_info[i] = tk_tao_parameter(
          str_to_tao_param(self.top_info[i]), self.top_frame, self.pipe)
      # Link the multipoles and elec_multipoles variables
      self.top_info[i].tk_var = (parent.mp_on_var if i == 0
          else parent.scale_mp_var)
      self.top_info[i].tk_var.set(self.top_info[i].param.value)
      self.top_info[i].tk_wid = tk.Checkbutton(
          self.top_frame, variable=self.top_info[i].tk_var)
      self.top_info[i].tk_label.grid(row=0, column=2*i, sticky='E')
      self.top_info[i].tk_wid.grid(row=0, column=2*i+1, sticky='W')
    tk.Button(self.top_frame, text="Show all orders",
        command=self.show_all_orders).grid(row=0, column=4, sticky='W')

    # Display the table (rest of output)
    if len(tao_output) < 3:
      return
    self.titles = tao_output[2]
    self.titles = ['Order'] + self.titles.split(';')
    for i in range(len(self.titles)):
      title = self.titles[i]
      tk.Label(self.table_frame, text=title).grid(row=0, column=i, sticky='EW')
    i = 1 #row counter
    for line in tao_output[3:]:
      line = line.split(';')
      tk.Label(self.table_frame, text=line[0]).grid(
          row=i, column=0, sticky='EW')
      self.shown_orders.append(line[0])
      j = 1 #column counter
      for item in line[1:]:
        name = self.titles[j]
        name = name[:name.find('n')] + str(j) + name[name.find('n')+1:]
        can_vary = 'T' if j<3 else 'F'
        self.tk_tao_list.append(tk_tao_parameter(str_to_tao_param(
          name + ';REAL;' + can_vary + ';' + item),
          self.table_frame, self.pipe))
        self.tk_tao_list[-1].tk_wid.grid(row=i, column=j, sticky='EW')
        j = j+1
      i = i+1

  def show_all_orders(self, event=None):
    '''
    Expands the table to include all orders of multipoles
    '''
    # Do nothing if all orders already shown
    if self.shown_orders == 'all':
      return

    # Clear self.table_frame
    for child in self.table_frame.winfo_children():
      child.grid_forget()
    # Re-grid the titles
    for j in range(len(self.titles)):
      tk.Label(self.table_frame, text=self.titles[j]).grid(
          row=0, column=j, sticky='EW')
    # the multipoles:
    x = 0 # counts position in self.tk_tao_list
    for i in range(22): #supports orders 0 through 21
      tk.Label(self.table_frame, text=str(i)).grid(
          row=i+1, column=0, sticky='EW')
      # Check if widgets exist for this order
      if str(i) in self.shown_orders:
        for j in range(len(self.titles)-1):
          self.tk_tao_list[x].tk_wid.grid(row=i+1, column=j+1, sticky='EW')
          x += 1
      else:
        #Create new widgets, grid them, and append to self.tk_tao_list
        for j in range(1, len(self.titles)):
          name = self.titles[j]
          name = name[:name.find('n')] + str(i) + name[name.find('n')+1:]
          can_vary = 'T' if j<3 else 'F'
          self.tk_tao_list.append(tk_tao_parameter(str_to_tao_param(
            name + ';REAL;' + can_vary + ';' + str(0)),
            self.table_frame, self.pipe))
          self.tk_tao_list[-1].tk_wid.grid(row=i+1, column=j, sticky='EW')
    self.shown_orders = 'all'

  def set_params(self, set_str, event=None):
    tao_set(self.top_info, set_str, self.pipe)
    tao_set(self.tk_tao_list, set_str, self.pipe)

  def check_for_changes(self):
    return (check_for_changes(self.tk_tao_list)
        | check_for_changes(self.top_info))

#---------------------------------------------------
# Lattice Window
class tao_lattice_window(tk.Toplevel):
  '''
  Shows lattice elements in a read-only table view
  with an interface to select which rows/columns
  are displayed
  '''
  # TODO: replace blank
  def __init__(self, root, pipe, switches="", *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
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
    self.temp_label_1 = tk.Label(self.top_frame, text="Template File: ")
    self.temp_label_1.grid(row=4, column=0, sticky='W')
    self.template_file = tk_tao_parameter(
        str_to_tao_param("template_file;FILE;T;"), self.top_frame, self.pipe)
    self.template_file.tk_wid.configure(command=self.temp_file_load)
    self.template_file.tk_wid.grid(row=4, column=1, sticky='EW')
    self.temp_label_2 = tk.Label(self.top_frame, text="Template:")
    self.temp_label_2.grid(row=4, column=2, sticky='W')
    self.temp_var = tk.StringVar() # Holds the chosen template
    self.temp_var.set("NO TEMPLATE FILE SELECTED")
    self.temp_chooser = tk.OptionMenu(self.top_frame, self.temp_var, [])
    self.temp_chooser.grid(row=4, column=3, sticky='EW')
    self.save_button = tk.Button(self.top_frame, text="Save", command=self.save_template)
    self.save_button.grid(row=4, column=5, sticky='W')
    self.temp_save = tk_tao_parameter(
        str_to_tao_param("name;STR;T;"), self.top_frame, self.pipe)
    self.temp_save.tk_wid.bind('<Return>', self.save_template)
    self.temp_save.tk_wid.grid(row=4, column=4, sticky='EW')

    # Branch/General
    self.branch_wids = tao_branch_widgets(self.top_frame, self.pipe)
    tk.Label(self.top_frame, text="Universe:").grid(row=0, column=0, sticky='W')
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
    tk.Label(self.top_frame, text="Attributes:").grid(
        row=1, column=2, sticky='EW')
    self.col_atts = tk.StringVar()
    self.att_box = tk.Entry(self.top_frame, textvariable=self.col_atts)
    self.att_box.configure(state="disabled")
    self.att_box.grid(row=1, column=3, sticky='EW')
    tk.Label(self.top_frame, text="Template file:").grid(
        row=1, column=4, sticky='EW')
    self.col_file = tk_tao_parameter(
        str_to_tao_param("col_file;FILE;T;"), self.top_frame, self.pipe)
    self.col_file.tk_wid.configure(state="disabled")
    self.col_file.tk_wid.configure(disabledforeground="grey")
    self.col_file.tk_wid.grid(row=1, column=5, sticky='EW')

    # Rows
    tk.Label(self.top_frame, text="Rows:").grid(row=2, column=0, sticky='W')
    self.f_button = tk.Button(
        self.top_frame, text="Filters", command=self.open_filter_menu)
    self.f_button.grid(row=2, column=1, sticky='EW')

    self.filter_vars = []
    for i in range(5):
      self.filter_vars.append(tk.BooleanVar())
      self.filter_vars[i].set(False)
    self.remove_if_zero = tk.StringVar()
    self.s_range = tk.StringVar()

    tk.Label(self.top_frame, text="Element List:").grid(
        row=2, column=2, sticky='E')
    self.ele_list = tk.StringVar()
    self.ele_list_opt = tk.StringVar()
    ele_list_opts = ["All", "Tracking elements", "Custom"]
    self.ele_list_chooser = tk.OptionMenu(self.top_frame, self.ele_list_opt,
        *ele_list_opts)
    self.ele_list_opt.set(ele_list_opts[0])
    self.ele_list_opt.trace('w', self.ele_list_callback)
    self.ele_list_chooser.grid(row=2, column=3, sticky='EW')
    self.ele_list_box = tk.Entry(self.top_frame, textvariable=self.ele_list)
    self.ele_list_box.configure(state="disabled")
    self.ele_list_box.grid(row=2, column=4, columnspan=3, sticky='EW')

    # Advanced Options
    self.use_advanced = False
    tk.Button(self.top_frame, text="Advanced\nOn/Off",
        command=self.toggle_advanced).grid(
            row=2, column=7, rowspan=2, sticky='NSEW')
    tk.Label(self.top_frame, text="Advanced:").grid(row=3, column=0, sticky='W')
    self.advanced_var = tk.StringVar()
    self.advanced_var.set(self.switches)
    self.advanced_box = tk.Entry(self.top_frame, textvariable=self.advanced_var)
    self.advanced_box.configure(state="disabled")
    self.advanced_box.grid(row=3, column=1, columnspan=6, sticky='EW')
    #self.advanced_box.bind("<Return>", self.refresh)
    #self.tree.focus() and self.tree.selection() for current items

    b = tk.Button(self.top_frame, text="Refresh", command=self.refresh)
    b.grid(row=0, column=7, rowspan=2, sticky='NSEW')
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
    self.temp_chooser = tk.OptionMenu(self.top_frame, self.temp_var, *temp_opts,
        command=self.temp_chooser_callback)
    self.temp_chooser.grid(row=4, column=3, sticky='EW')

  def temp_file_load(self, event=None):
    '''
    Tries to load the specified template file, and creates an OptionMenu to pick
    from defined templates if successful.  Also creates a save template button
    and entry for the name
    self.temp_dict is a dictionary whose keys are the names of templates and
    whose values are the switches for those templates
    '''
    # First run the tk_tao_paramter open file method that is being overloaded
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
    self.temp_chooser = tk.OptionMenu(self.top_frame, self.temp_var, *temp_opts,
        command=self.temp_chooser_callback)
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

  def ele_list_callback(self, event=None, *args):
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
    filters = ["Lords only", "No Slaves", "No Super Slaves",
        "Remove Line if zero", "s"]
    buttons = []
    for i in range(len(self.filter_vars)):
      buttons.append(
          tk.Checkbutton(win, text=filters[i], variable=self.filter_vars[i]))
      buttons[i].grid(row=i, column=0, sticky='W')
    tk.Entry(win, textvariable=self.remove_if_zero).grid(
        row=3, column=1, sticky='EW')
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
        #is a switch, last item was not an arg switch
        if (item[0] == '-') & (not arg_switch):
          if item in ['-att', '-blank_replacement', '-branch', '-custom',
              '-remove_line_if_zero', '-s', '-universe']:
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
          messagebox.showwarning('Warning', 'No argument given for switch "'
              + arg_switch + '" (one required).  Skipping this switch...', parent=self)
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
      tk.Label(self.table_frame, text="NO LATTICE FOUND").pack()
      return
    for i in range(len(lattice)):
      lattice[i] = lattice[i].split(';')
    #lattice[i][j] --> row i, column j
    widths = [0]*len(lattice[0]) # tracks column widths

    # Create table
    self.tree = ttk.Treeview(
        self.table_frame, columns=lattice[0], show='headings')
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
      widths[0] = 120 #prevent giant index column
      self.tree.column(lattice[0][j], width=widths[j], minwidth=widths[j])

    # Scrollbars
    hbar = ttk.Scrollbar(
        self.table_frame, orient="horizontal", command=self.tree.xview)
    vbar = ttk.Scrollbar(
        self.table_frame, orient="vertical", command=self.tree.yview)
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
    #checks that the first value is an element index
    if isinstance(row['values'][0], int):
      settings = [self.branch_wids.uni.get(),
          self.branch_wids.branch.get(),
          str(row['values'][0]),
          self.branch_wids.bmd.get()]
      win = tao_ele_window(self.root, self.pipe, settings)

class tao_new_data_window(tk.Toplevel):
  '''
  Provides a window for creating new d2 data arrays (and their associated
  d1 arrays)
  '''
  def __init__(self, root, pipe, *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
    self.root = root
    self.pipe = pipe
    self.rowconfigure(0, weight=1)
    self.title('New Data')
    self.name = ""

    #Frame for inputting d2 parameters
    self.d2_frame = tk.Frame(self)
    self.d2_frame.grid(row=0, column=0, sticky='NSEW')

    #Frames for inputting d1 parameters
    self.d1_frame = tk.Frame(self)
    self.nb_exists = False #used to track if the ttk Notebook has been set up
    self.d1_frame_list = []

    self.fill_d2_frame()

  def fill_d2_frame(self):
    tk.Label(self.d2_frame, text="New d2 Data",
        font=('Sans', 16, 'bold')).grid(
            row=0, column=0, columnspan=2, sticky='EW')

    # Small bit of setup
    self.uni = tk.StringVar()
    self.uni.set("All")
    uni_max = self.pipe.cmd_in("python super universe")
    uni_max = str_to_tao_param(uni_max.splitlines()[0]).value
    uni_list = ["All"] + list(range(1, uni_max+1))
    def my_ttp(x):
      ''' shorcut for the following commonly used construct '''
      return tk_tao_parameter(str_to_tao_param(x), self.d2_frame, self.pipe)

    # Clone existing d2
    existing_d2_arrays = self.pipe.cmd_in("python data_d2_array 1")
    existing_d2_arrays = ['None'] + existing_d2_arrays.splitlines()
    self.clone_d2 = tk.StringVar()
    self.clone_d2.set('None')

    # Widgets
    self.d2_param_list = [my_ttp("name;STR;T;"),
        tk.OptionMenu(self.d2_frame, self.uni, *uni_list),
        my_ttp("data_source;ENUM;T;"),
        my_ttp("data^merit_type;ENUM;T;"),
        my_ttp("weight;REAL;T;"),
        my_ttp("good_user;LOGIC;T;T"),
        tk.OptionMenu(self.d2_frame, self.clone_d2, *existing_d2_arrays)]

    # Labels
    self.d2_label_list = [tk.Label(self.d2_frame, text="d2_data Name:"),
        tk.Label(self.d2_frame, text="Which universe:"),
        tk.Label(self.d2_frame, text="Default data source:"),
        tk.Label(self.d2_frame, text="Default merit type:"),
        tk.Label(self.d2_frame, text="Default weight:"),
        tk.Label(self.d2_frame, text="Default good_user:"),
        tk.Label(self.d2_frame, text="Clone existing d2_array:")]

    # Grid widgets and labels
    for i in range(len(self.d2_param_list)):
      self.d2_label_list[i].grid(row=i+1, column=0, sticky='W')
      if i in [1, 6]: #non tk_tao_parameters
        self.d2_param_list[i].grid(row=i+1, column=1, sticky='EW')
      else:
        self.d2_param_list[i].tk_wid.grid(row=i+1, column=1, sticky='EW')

    # Warning labels
    self.name_warning_1 = tk.Label(self.d2_frame, text="Cannot be empty")
    # Next button
    self.next_b = tk.Button(self.d2_frame, text="Next", command=self.load_d1_frame)
    self.next_b.grid(row=i+2, column=1, sticky='W')

    # Focus the name entry
    self.d2_param_list[0].tk_wid.focus_set()

  def load_d2_frame(self):
    self.d1_frame.pack_forget()
    self.d2_frame.grid(row=0, column=0, sticky='NSEW')

  def load_d1_frame(self):
    '''
    Ungrids self.d2_frame, grids self.d1_frame, and sets up a notebook
    for the d1_frames if necessary.
    '''
    clone_dict = {} # keys=d2_array names, values=lists of d1s
    # Check if d2_name is nonempty
    name = self.d2_param_list[0].tk_var.get().strip()
    if name == "":
      self.name_warning_1.grid(row=1, column=2, sticky='W')
      return

    # Conditions
    c1 = (name in self.pipe.cmd_in('python data_d2_array 1').splitlines()) #name in use
    c2 = (self.clone_d2.get() != 'None') #clone has been specified
    c3 = (self.name != name) #name has changed
    c4 = False #discard existing work
    ans_var = tk.StringVar()
    # Ask if user wants to keep existing data
    if c3 and self.name:
      self.name = name
      tao_message_box(self.root, self, ans_var, title='Warning', message='Would you like to keep or discard the d1_arrays you defined for ' + self.name + '?', choices=['Keep', 'Discard'])
      if ans_var.get() == 'Keep':
        c4 = False
      elif ans_var.get() == 'Discard':
        c4 = True
      else:
        return
    # Ask if user wants to load existing data
    if c1:
      ans = messagebox.askyesno('Warning', name + " already exists as a d2 data array.  Would you like to clone its existing d1_arrays?", parent=self)
      if ans:
        # Will need to read in data for existing d2
        clone_dict[name] = self.pipe.cmd_in('python data_d1_array ' + name).splitlines()
    # Clone existing data
    if c2 and (self.clone_d2.get() != name):
      clone_dict[self.clone_d2.get()] = self.pipe.cmd_in('python data_d1_array ' + self.clone_d2.get()).splitlines()
    # Remake d1_frames and notebook if necessary
    if c4 and self.nb_exists:
      for frame in self.d1_frame_list:
        frame.destroy()
      self.d1_frame_list = []
      self.notebook.destroy()
      self.new_tab_frame.destroy()
      self.back_b.destroy()
      self.create_b.destroy()
      self.nb_exists = False
    else:
      self.name = name
      self.name_warning_1.grid_forget()

    # Possibly create self.notebook
    if not self.nb_exists:
      self.notebook = ttk.Notebook(self.d1_frame)
      self.notebook.pack(side='top', fill='both', expand=1)
      self.nb_exists = True

      self.d1_index = 0 #marks current tab index

      # New tab button
      self.new_tab_frame = tk.Frame(self.notebook)
      self.notebook.insert('end', self.new_tab_frame)
      self.notebook.tab(len(self.d1_frame_list), text='+')
      self.notebook.bind('<<NotebookTabChanged>>', self.tab_handler)

      # Back button
      self.back_b = tk.Button(self.d1_frame, text="Back", command=self.load_d2_frame)
      self.back_b.pack(side='left')

      # Create button
      self.create_b = tk.Button(self.d1_frame, text="Create!", command=self.create_data)
      self.create_b.pack(side='right')

    self.d2_frame.grid_forget()
    self.d1_frame.pack(fill='both', expand=1)
    self.title("New data: " + self.d2_param_list[0].tk_var.get())

    # Clone the requested d2 array(s)
    if clone_dict != {}:
      # Progress window
      num_bars = 0
      for k in clone_dict.keys():
        num_bars += len(clone_dict[k])
      self.pw = tao_progress_window(self.root, self, num_bars)
      self.pw.title("Loading...")
      tk.Label(self.pw, text='Loading existing d1_arrays').grid(row=0, column=0, columnspan=2, sticky='EW')
      # Find correct universe number
      if self.uni.get() == "All":
        u = "1"
      else:
        u = self.uni.get()
      if name in clone_dict.keys():
        for d1 in clone_dict[name]:
          self.add_d1_frame(d1, u+'@'+name)
      if (self.clone_d2.get() in clone_dict.keys()) and (name != self.clone_d2.get()):
        for d1 in clone_dict[self.clone_d2.get()]:
          self.add_d1_frame(d1, u+'@'+self.clone_d2.get())
      self.notebook.select(0)
      self.pw.destroy()

    # Copy d2 defaults into d1_arrays
    for i in range(len(self.d1_frame_list)):
      val = self.d2_param_list[2].tk_var.get()
      self.d1_frame_list[i].d1_array_wids[1].tk_var.set(val)
      val = self.d2_param_list[3].tk_var.get()
      self.d1_frame_list[i].d1_array_wids[3].tk_var.set(val)
      val = self.d2_param_list[4].tk_var.get()
      self.d1_frame_list[i].d1_array_wids[4].tk_var.set(val)

  def create_data(self, event=None):
    '''
    Takes the information from the d2_frame and d1_frames and runs
    the necessary commands to create the data in tao, then closes
    the create data window
    '''
    # Input validation
    messages = []
    for d1_frame in self.d1_frame_list:
      # Check names
      if d1_frame.name_handler():
        name_m = "Please check d1_array names."
        if name_m not in messages:
          messages.append(name_m)
      if d1_frame.ix_min_handler():
        messages.append("Please check the start index for " + d1_frame.name)
      if d1_frame.ix_max_handler():
        messages.append("Please check the end index for " + d1_frame.name)
      if d1_frame.d1_array_wids[2].tk_var.get() == "":
        messages.append("Please choose a data type for " + d1_frame.name)
    for m in messages:
      messagebox.showwarning("Error", m, parent=self)
    if messages != []:
      return
    # Book-keeping
    datum_params = ['data_type', 'ele_ref_name', 'ele_start_name', 'ele_name',
        'data^merit_type', 'meas_value', 'ref_value', 'weight', 'good_user',
        'data_source', 'eval_point', 's_offset', '1^ix_bunch',
        'invalid_value', 'spin_n0_x', 'spin_n0_y', 'spin_n0_z']
    d1_params = ['name', 'data_source', 'data_type', 'data^merit_type',
        'weight', 'good_user']
    d2_params = ['name', 'uni', 'data_source', 'data^merit_type', 'weight', 'good_user']
    # Create the data array
    if self.uni.get() == 'All':
      uni_max = self.pipe.cmd_in("python super universe")
      uni_max = str_to_tao_param(uni_max.splitlines()[0]).value
      uni_list = list(range(1, uni_max+1))
    else:
      uni_list = [int(self.uni.get())]
    d1_count = 0 # used to make the corrent number of progress bars
    for u in uni_list:
      cmd_str = 'python data_d2_create ' + str(u) + '@' + self.name
      cmd_str += ' ' + str(len(self.d1_frame_list)) + ' '
      for d1_frame in self.d1_frame_list:
        d1_count += 1
        # min/max indices for each d1_array
        cmd_str += str(d1_frame.ix_min) + ' '
        cmd_str += str(d1_frame.ix_max) + ' '
      # Create the d2/d1_arrays
      self.pipe.cmd_in(cmd_str)
    # Progress bars
    self.pw = tao_progress_window(self.root, self, d1_count)
    for u in uni_list:
      # Set parameters at the d2 level
      #set_str = 'set data ' + str(u) + '@' + self.name + '|'
      #tao_set(self.d2_param_list[2:5], set_str, self.pipe)
      # Set parameters at the d1 level
      # Set the names last for convenience
      i = 1 #temp names set by tao
      for d1_frame in self.d1_frame_list:
        self.pw.label_vars[self.pw.ix].set(
            'Creating' + str(u) + '@' + self.name + '.' + d1_frame.name)
        self.pw.set_max(self.pw.ix, d1_frame.ix_max-d1_frame.ix_min+1)
        #self.update_idletasks()
        #set_str = 'set data ' + str(u) + '@' + self.name + '.' + str(i) + '|'
        #tao_set(d1_frame.d1_array_wids[1:5], set_str, self.pipe)
        # set individual data parameters
        #set_format = set_str[:-1] + '[{}]|'
        #tao_dict_set(d1_frame.data_dict, set_format, self.pipe)
        # set name now
        #tao_set(d1_frame.d1_array_wids[0:1], set_str, self.pipe)
        for j in range(d1_frame.ix_min, d1_frame.ix_max+1):
          self.pw.set_val(self.pw.ix, j-d1_frame.ix_min)
          #self.update_idletasks()
          cmd_str = 'python datum_create '
          cmd_str += str(u) + '@' + self.name + '.' + str(i) + '[' + str(j) + ']'
          for p in datum_params:
            #look in d1_frame.data_dict
            if (j in d1_frame.data_dict.keys()) and (p in d1_frame.data_dict[j].keys()):
              cmd_str += '^' + d1_frame.data_dict[j][p]
            elif p in d1_params:
              cmd_str += '^' + d1_frame.d1_array_wids[d1_params.find(p)].tk_var.get()
            elif p in d2_params:
              cmd_str += '^' + self.d2_param_list[d1_params.find(p)].tk_var.get()
            else:
              cmd_str += '^'
          self.pipe.cmd_in(cmd_str)
        i = i+1
        self.pw.ix += 1
      # set data|exists = T
      #self.pipe.cmd_in('set data ' + str(u) + '@' + self.name + '|exists = T')
    # Close the window
    self.destroy()

  def tab_handler(self, event=None):
    '''
    Handles new tab creation and updates self.d1_index as necessary
    '''
    # Check if the new tab frame has been selected
    if self.notebook.select() == self.new_tab_frame._w:
      # Add new tab
      self.d1_frame_list.append(new_d1_frame(self))
      self.d1_index = len(self.d1_frame_list)-1
      self.notebook.insert(self.d1_index, self.d1_frame_list[-1])
      self.notebook.tab(self.d1_index, text='New d1_array')
      self.notebook.select(self.d1_index)
    else:
      # Update self.d1_index
      for i in range(len(self.d1_frame_list)):
        frame = self.d1_frame_list[i]
        if self.notebook.select() == frame._w:
          self.d1_index = i
          # Unblock frame's handlers
          frame.handler_block = False
          break

  def add_d1_frame(self, d1_name, d2_name, event=None):
    '''
    Creates a new d1_frame for d1_name and reads existing data in from Tao if present
    Takes a line of output from python data_d1_array for d1_name
    d2_name should be "uni@d2_name"
    '''
    d1 = d1_name.split(';')[3]
    self.d1_frame_list.append(
        new_d1_frame(self, name=d1, full_name=d2_name + '.' + d1))
    ix = len(self.d1_frame_list) - 1
    self.notebook.insert(ix, self.d1_frame_list[-1])
    self.notebook.tab(ix, text=d1)


class new_d1_frame(tk.Frame):
  '''
  Provides a frame for inputting properties of a d1_data array.
  To load an existing d1_array, pass the d1_array name as name
  (e.g. x,y, NOT orbit.x, orbit.y)
  Also pass the full name that should be used for data lookup
  (e.g. 1@orbit.x)
  '''
  def __init__(self, d2_array, name="", full_name=""):
    tk.Frame.__init__(self, d2_array.notebook)
    self.d2_array = d2_array
    self.pipe = self.d2_array.pipe
    self.handler_block = False
    if name == "":
      self.name = "New d1_array" #Default
    else:
      self.name = name

    # d1 Widgets
    def d1_ttp(x):
      ''' Shortcut for commonly used construct '''
      return tk_tao_parameter(str_to_tao_param(x), self, self.pipe)

    self.d1_array_wids = [d1_ttp('name;STR;T;'),
        d1_ttp('data_source;ENUM;T;'),
        d1_ttp('data_type;DAT_TYPE;T;'),
        d1_ttp('data^merit_type;ENUM;T;'),
        d1_ttp('weight;REAL;T;'),
        d1_ttp('good_user;LOGIC;T;T'),
        d1_ttp('ix_min;INT;T;'),
        d1_ttp('ix_max;INT;T;')]
    # d1 labels (NAMES AS STRINGS ONLY)
    self.d1_array_labels = ["d1_array Name:", "Default data source:",
        "Default data type:", "Default merit type:", "Default weight",
        "Default good_user:", "Start index", "End index"]
    # Read in defaults from d2 level
    val = self.d2_array.d2_param_list[2].tk_var.get()
    self.d1_array_wids[1].tk_var.set(val)
    val = self.d2_array.d2_param_list[3].tk_var.get()
    self.d1_array_wids[3].tk_var.set(val)
    val = self.d2_array.d2_param_list[4].tk_var.get()
    self.d1_array_wids[4].tk_var.set(val)
    val = self.d2_array.d2_param_list[5].tk_var.get()
    self.d1_array_wids[5].tk_var.set(val)
    # Grid widgets and labels:
    for i in range(len(self.d1_array_wids)):
      tk.Label(self, text=self.d1_array_labels[i]).grid(row=i+2, column=0, sticky='W')
      self.d1_array_wids[i].tk_wid.grid(row=i+2, column=1, sticky='EW')
    i = i+2
    # Set name
    if self.name != "New d1_array":
      self.d1_array_wids[0].tk_var.set(self.name)

    # Warning labels
    # (defined here to be gridded/ungridded as necessary)
    self.name_warning_1 = tk.Label(self, text="Must not be empty")
    self.name_warning_2 = tk.Label(self, text="d1 name already in use")
    self.ix_min_warning_1 = tk.Label(self, text="Must be a non-negative integer")
    self.ix_min_warning_2 = tk.Label(self, text="Cannot be larger than the maximum index")
    self.ix_max_warning_1 = tk.Label(self, text="Must be a non-negative integer")
    self.ix_max_warning_2 = tk.Label(self, text="Cannot be smaller than the minimum index")

    # Responses to edits
    self.d1_array_wids[0].tk_wid.bind('<FocusOut>', self.name_handler)
    self.d1_array_wids[6].tk_wid.bind('<FocusOut>', self.ix_min_handler)
    self.d1_array_wids[7].tk_wid.bind('<FocusOut>', self.ix_max_handler)
    self.d1_array_wids[1].tk_var.trace('w', self.data_source_handler)
    self.d1_array_wids[2].tk_var.trace('w', self.data_type_handler)

    ttk.Separator(self, orient='horizontal').grid(row=i+1, column=0, columnspan=3, sticky='EW')
    i = i+1

    # Element browsers
    tk.Label(self, text="Choose elements:").grid(row=i+1, column=0, sticky='W')
    self.ele_name_button = tk.Button(self, text="Browse...",
        command=lambda : self.lat_browser('name'))
    self.ele_name_button.grid(row=i+1, column=1, sticky='EW')

    tk.Label(self, text="Start elements:").grid(row=i+2, column=0, sticky='W')
    self.ele_start_name_button = tk.Button(self, text="Browse...",
        command=lambda : self.lat_browser('start_name'))
    self.ele_start_name_button.grid(row=i+2, column=1, sticky='EW')

    tk.Label(self, text="Ref elements:").grid(row=i+3, column=0, sticky='W')
    self.ele_ref_name_button = tk.Button(self, text="Browse...",
        command=lambda : self.lat_browser('ref_name'))
    self.ele_ref_name_button.grid(row=i+3, column=1, sticky='EW')

    ttk.Separator(self, orient='horizontal').grid(row=i+4, column=0, columnspan=3, sticky='EW')
    i = i+1

    # Individual data
    self.data_dict = {}
    self.datum_frame = tk.Frame(self)
    tk.Label(self, text="Datum:").grid(row=i+4, column=0, sticky='W')
    self.data_ix = tk.StringVar()
    self.ix_min = -1
    self.ix_max = -1
    self.data_chooser = ttk.Combobox(self, textvariable=self.data_ix,
        values=['PLEASE SPECIFY MIN/MAX INDICES'], state='readonly')
    self.data_chooser.bind('<<ComboboxSelected>>', self.make_datum_frame)
    self.data_chooser.grid(row=i+4, column=1, sticky='EW')
    # Set ix_min and ix_max for existing d1_arrays
    if self.name != "New d1_array":
      ix_data = self.pipe.cmd_in('python data_d1_array ' + full_name)
      ix_data = ix_data.splitlines()
      for line in ix_data:
        if line.split(';')[3] == self.name:
          break
      #line is now set to the relevant line
      self.d1_array_wids[-2].tk_var.set(line.split(';')[5])
      self.d1_array_wids[-1].tk_var.set(line.split(';')[6])
      self.ix_min_handler()
      self.ix_max_handler()
    # Fill self.data_dict for existing d1 arrays
    data_dict_params = ['data_source', 'data_type', 'ele_name', 'ele_start_name', 'ele_ref_name',
        'data^merit_type', 'meas_value', 'ref_value', 'weight', 'good_user', '1^ix_bunch',
        'eval_point', 's_offset']
    if self.name != "New d1_array":
      self.d2_array.pw.label_vars[self.d2_array.pw.ix].set(
          'Loading ' + self.d2_array.name + '.' + self.name + '...')
      self.d2_array.pw.set_max(self.d2_array.pw.ix, self.ix_max-self.ix_min+1)
      #self.d2_array.update_idletasks()
      for i in range(self.ix_min, self.ix_max+1):
        self.d2_array.pw.set_val(self.d2_array.pw.ix, i-self.ix_min)
        #self.d2_array.update_idletasks()
        self.data_dict[i] = {}
        existing_data = self.pipe.cmd_in(
            'python data ' + full_name + '[' + str(i) + ']')
        existing_data = existing_data.splitlines()
        for p in data_dict_params:
          # Find correct line
          for line in existing_data:
            if line.find(p + ';') == 0:
              break
          # Write to self.data_dict
          self.data_dict[i][p] = line.split(';')[3].strip()
      self.d2_array.pw.ix += 1
      # Set default data_type to data_type of first datum
      self.d1_array_wids[1].tk_var.set(self.data_dict[self.ix_min]['data_source'])
      self.d1_array_wids[2].tk_var.set(self.data_dict[self.ix_min]['data_type'])
      # Load first datum
      self.data_ix.set(self.ix_min)
      self.make_datum_frame()

    # Delete button
    tk.Button(self, text="DELETE THIS D1_ARRAY", fg='red', command=self.delete).grid(
        row=0, column=0, columnspan=3, sticky='EW')

    # Duplicate button
    tk.Button(self, text="Duplicate this d1_array", command=self.duplicate).grid(
        row=1, column=0, columnspan=3, sticky='EW')

    # Focus the d1 name widget
    self.d1_array_wids[0].tk_wid.focus_set()

  def delete(self, ask=True, event=None):
    '''
    Deletes this d1_array frame
    Call with ask = False to skip confirmation
    '''
    # Ask for confirmation
    if ask:
      ans = messagebox.askokcancel("Delete " + self.name, "Delete this d1_array and its associated data?", parent=self.d2_array)
      if not ans:
        return

    # Remove from parent d1_frame_list
    ix = self.d2_array.d1_index
    if ix > 0:
      # Prefer to move to the tab to the left
      self.d2_array.notebook.select(ix-1)
    self.d2_array.d1_frame_list.pop(ix)

    # Destroy self
    self.destroy()

  def duplicate(self, event=None):
    '''
    Adds a new d1_frame to self.d2_array.d1_frame_list that is a copy of
    this frame, and changes focus to that frame
    '''
    # Don't run any handlers for this d1_frame
    self.handler_block = True
    self.d2_array.d1_frame_list.append(new_d1_frame(self.d2_array))
    ix = len(self.d2_array.d1_frame_list) - 1
    #self.d2_array.d1_index = ix
    # Copy properties into new frame
    self.d2_array.d1_frame_list[-1].name = self.name + '_copy'
    self.d2_array.d1_frame_list[-1].ix_min = self.ix_min
    self.d2_array.d1_frame_list[-1].ix_max = self.ix_max
    self.d2_array.d1_frame_list[-1].data_dict = self.data_dict
    for i in range(len(self.d1_array_wids)):
      if i == 0:
        self.d2_array.d1_frame_list[-1].d1_array_wids[i].tk_var.set(self.d1_array_wids[i].tk_var.get() + '_copy')
      else:
        self.d2_array.d1_frame_list[-1].d1_array_wids[i].tk_var.set(self.d1_array_wids[i].tk_var.get())
    # Run all input validation handlers
    self.d2_array.notebook.insert(ix, self.d2_array.d1_frame_list[-1])
    self.d2_array.notebook.select(ix)
    self.d2_array.tab_handler()
    self.update_idletasks()
    self.d2_array.d1_frame_list[-1].name_handler()
    self.d2_array.d1_frame_list[-1].ix_min_handler()
    self.d2_array.d1_frame_list[-1].ix_max_handler()
    self.d2_array.d1_frame_list[-1].data_source_handler()
    self.d2_array.d1_frame_list[-1].data_type_handler()

  def make_datum_frame(self, event=None):
    '''
    Adds an entry to self.data_dict if necessary, and updates
    self.datum_frame to show the values for the current datum
    '''
    ix = self.data_ix.get()
    if ix == 'PLEASE SPECIFY MIN/MAX INDICES':
      return
    ix = int(ix)
    if ix not in self.data_dict.keys():
      self.data_dict[ix] = {}
    self.datum_frame.destroy()
    self.datum_frame = tk.Frame(self)
    def datum_ttp(x):
      ''' Shortcut for commonly used construct '''
      return tk_tao_parameter(x, self.datum_frame, self.pipe)
    # Parameters
    param_list = ['data_source;ENUM;T;', 'data_type;DAT_TYPE;T;',
        'ele_name;STR;T;', 'ele_start_name;STR;T;', 'ele_ref_name;STR;T;',
        'data^merit_type;ENUM;T;', 'meas_value;REAL;T;', 'ref_value;REAL;T;',
        'weight;REAL;T;', 'good_user;LOGIC;T;T', '1^ix_bunch;INUM;T;',
        'eval_point;ENUM;T;', 's_offset;REAL;T;']
    param_list = list(map(str_to_tao_param, param_list))
    # Fill in defaults set at d1 level
    param_list[0].value = self.d1_array_wids[1].tk_var.get()
    param_list[1].value = self.d1_array_wids[2].tk_var.get()
    param_list[5].value = self.d1_array_wids[3].tk_var.get()
    param_list[8].value = self.d1_array_wids[4].tk_var.get()
    param_list[9].value = bool(self.d1_array_wids[5].tk_var.get())
    # Set parameter values if specified
    for i in range(len(param_list)):
      if param_list[i].name in self.data_dict[ix].keys():
        param_list[i].value = self.data_dict[ix][param_list[i].name]
    # Create widgets
    self.datum_wid_list = list(map(datum_ttp, param_list))
    def data_writer(i):
      ''' Writes the contents of self.datum_wid_list[i] into self.data_dict '''
      if self.datum_wid_list[i].param.type == 'LOGIC':
        val = bool(self.datum_wid_list[i].tk_var.get())
      else:
        val = self.datum_wid_list[i].tk_var.get()
      self.data_dict[ix][self.datum_wid_list[i].param.name] = val
    def data_writer_callback(i):
      ''' Callback for data_writer() '''
      return lambda *args : data_writer(i)

    # Set state of ele_names and s_offset for selected data_type
    def data_type_callback(*args):
      # ele_names
      if self.datum_wid_list[1]._has_ele():
        for i in range(2, 5):
          self.datum_wid_list[i].tk_wid.configure(state='normal')
      else:
        for i in range(2, 5):
          self.datum_wid_list[i].tk_wid.configure(state='disabled')
      # s_offset
      if self.datum_wid_list[1]._has_s_offset():
        self.datum_wid_list[11].tk_wid.configure(state='normal')
      else:
        self.datum_wid_list[11].tk_wid.configure(state='disabled')
    # Grid widgets
    for i in range(len(self.datum_wid_list)):
      self.datum_wid_list[i].tk_label.grid(row=i, column=0, sticky='E')
      self.datum_wid_list[i].tk_wid.grid(row=i, column=1, sticky='EW')
      # Write changes into self.data_dict
      self.datum_wid_list[i].tk_var.trace('w', data_writer_callback(i))
    self.datum_wid_list[1].tk_var.trace('w', data_type_callback)
    self.datum_frame.grid(row=20, column=0, columnspan=3, sticky='EW')


  def lat_browser(self, which):
    '''
    Opens a modified lattice table window to select elements en masse
    which should be either 'name', 'start_name', or 'ref_name'
    When the elements are chosen, ele_which will be set for each of this
    d1_array's data sequentially
    '''
    # Make sure the name, ix_min, and ix_max are set
    if (bool(self.name_handler()) | bool(self.ix_min_handler())
        | bool(self.ix_max_handler())):
      return
    #
    name = self.d2_array.name + '.' + self.name
    win = tao_ele_browser(self.d2_array.root, self.pipe, name, self,
        'data', which, self.d2_array.uni.get())

  def name_handler(self, event=None):
    '''
    Changes the tab name to match the d1_name
    Returns 1 if unsuccessful
    '''
    if self.handler_block:
      return
    name = self.d1_array_wids[0].tk_var.get().strip()
    if name != "":
      # Make sure the name isn't already in use
      i = 0
      for d1 in self.d2_array.d1_frame_list:
        if (d1.name == name) & (self != self.d2_array.d1_frame_list[i]):
          self.name_warning_1.grid_forget()
          self.name_warning_2.grid(row=2, column=2, sticky='W')
          return 1
        i = i+1

      self.d2_array.notebook.tab(self.d2_array.d1_index, text=self.name)
      self.name_warning_1.grid_forget()
      self.name_warning_2.grid_forget()
      if self.name != name:
        if self.d1_array_wids[2]._is_valid_dat_type(self.d2_array.name + '.' + name):
          # set data type
          self.d1_array_wids[2].tk_var.set(self.d2_array.name + '.' + name)
      self.name = name
    else:
      self.name_warning_2.grid_forget()
      self.name_warning_1.grid(row=2, column=2, sticky='W')
      self.name = "New d1_array"
      self.d2_array.notebook.tab(self.d2_array.d1_index, text="New d1_array")
      return 1

  def ix_min_handler(self, event=None):
    if self.handler_block:
      return
    try:
      ix_min = int(self.d1_array_wids[6].tk_var.get())
    except ValueError:
      ix_min = -1
    if ix_min > -1:
      self.ix_min = ix_min
      if (ix_min <= self.ix_max) | (self.ix_max == -1):
        self.ix_min_warning_1.grid_forget()
        self.ix_min_warning_2.grid_forget()
      else:
        self.ix_min_warning_1.grid_forget()
        self.ix_min_warning_2.grid(row=8, column=2, sticky='W')
        return 1
    else:
      self.ix_min_warning_2.grid_forget()
      self.ix_min_warning_1.grid(row=8, column=2, sticky='W')
      return 1
    # Update the data index range if possible
    if (self.ix_min > -1) & (self.ix_max >= self.ix_min):
      self.data_chooser.configure(values=list(range(self.ix_min, self.ix_max+1)))

  def ix_max_handler(self, event=None):
    try:
      ix_max = int(self.d1_array_wids[7].tk_var.get())
    except ValueError:
      ix_max = -1
    if ix_max > -1:
      self.ix_max = ix_max
      if ix_max >= self.ix_min:
        self.ix_max_warning_1.grid_forget()
        self.ix_max_warning_2.grid_forget()
      else:
        self.ix_max_warning_1.grid_forget()
        self.ix_max_warning_2.grid(row=9, column=2, sticky='W')
        return 1
    else:
      self.ix_max_warning_2.grid_forget()
      self.ix_max_warning_1.grid(row=9, column=2, sticky='W')
      return 1
    # Update the data index range if possible
    if (self.ix_min > -1) & (self.ix_max >= self.ix_min):
      self.data_chooser.configure(values=list(range(self.ix_min, self.ix_max+1)))

  def data_source_handler(self, *args):
    '''
    Updates the data_source widget to only show valid data types
    for the currently selected data source
    '''
    if self.handler_block:
      return
    source = self.d1_array_wids[1].tk_var.get()
    if source in ['beam', 'lat']:
      # Remake self.d1_array_wids[2] (data type)
      old_val = self.d1_array_wids[2].tk_var.get()
      self.d1_array_wids[2].tk_wid.grid_forget()
      self.d1_array_wids[2] = tk_tao_parameter(
          str_to_tao_param('data_type;DAT_TYPE;T;'),
          self, self.pipe, data_source=source)
      self.d1_array_wids[2].tk_var.set(old_val)
      self.d1_array_wids[2].tk_var.trace('w', self.data_type_handler)
      self.d1_array_wids[2].tk_wid.grid(row=4, column=1, sticky='EW')

  def data_type_handler(self, *args):
    '''
    Sets the state of the ele buttons appropriately for the currently selected data type
    '''
    if self.handler_block:
      return
    if not self.d1_array_wids[2]._has_ele():
      self.ele_name_button.configure(state='disabled')
      self.ele_start_name_button.configure(state='disabled')
      self.ele_ref_name_button.configure(state='disabled')
    else:
      self.ele_name_button.configure(state='normal')
      self.ele_start_name_button.configure(state='normal')
      self.ele_ref_name_button.configure(state='normal')

class tao_new_var_window(tk.Toplevel):
  '''
  Provides a window for creating new v1_variable arrays
  '''
  def __init__(self, root, pipe, *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
    self.root = root
    self.pipe = pipe
    self.title('New Variables')
    self.v1_frame_list = []

    # Possibly create self.notebook
    self.notebook = ttk.Notebook(self)
    self.notebook.pack(side='top', fill='both', expand=1)

    # Set up a v1 tab
    self.v1_frame_list.append(new_v1_frame(self))
    self.notebook.insert('end', self.v1_frame_list[0])
    self.notebook.tab(0, text='New v1_array')
    self.v1_index = 0 #marks current tab index

    # New tab button
    self.new_tab_frame = tk.Frame(self.notebook)
    self.notebook.insert('end', self.new_tab_frame)
    self.notebook.tab(1, text='+')
    self.notebook.bind('<<NotebookTabChanged>>', self.tab_handler)

    # Create button
    self.create_b = tk.Button(self, text="Create!", command=self.create_variables)
    self.create_b.pack(side='right')

  def create_variables(self, event=None):
    '''
    Takes the information from the variable frames and runs
    the necessary commands to create the variables in tao, then closes
    the create variable window
    '''
    # Input validation
    messages = []
    missing_atts = [] #variables with element must have an attribute
    for v1_frame in self.v1_frame_list:
      # Check names
      if v1_frame.name_handler():
        name_m = "Please check v1_array names."
        if name_m not in messages:
          messages.append(name_m)
      # Check min indices
      if v1_frame.ix_min_handler():
        messages.append("Please check the start index for " + v1_frame.name)
      # Check max indices
      if v1_frame.ix_max_handler():
        messages.append("Please check the end index for " + v1_frame.name)
      # Check low/high limits
      if v1_frame.low_high_handler():
        messages.append("Please check the low and high limits for " + v1_frame.name)
      # Check for semicolons in any fields
      semi_message = "Semicolons not allowed in any input field"
      caret_message = "Carets not allowed in any input field"
      broken = False #Used to break out of the below for loops
      # Check for semicolons/carets
      for ttp in v1_frame.v1_array_wids:
        if str(ttp.tk_var.get()).find(';') != -1:
          messages.append(semi_message)
          broken = True
          break
        if str(ttp.tk_var.get()).find('^') != -1:
          messages.append(caret_message)
          broken = True
          break
      for var_dict in v1_frame.var_dict.values():
        if broken:
          break
        for v in var_dict.values():
          if str(v).find(';') != -1:
            messages.append(semi_message)
            broken = True
            break
          if str(v).find('^') != -1:
            messages.append(caret_message)
            broken = True
            break
        if broken:
          break
      # Check that variables with elements have attributes
      for ix, var_dict in v1_frame.var_dict.items():
        if 'ele_name' in var_dict.keys():
          if var_dict['ele_name'] != "":
            if 'attribute' in var_dict.keys():
              if var_dict['attribute'] == "":
                missing_atts.append(ix)
            elif v1_frame.v1_array_wids[2].tk_var.get() == "":
              missing_atts.append(ix)
    if missing_atts != []:
      missing_atts.sort()
      att_msg = "The following variables were assigned elements but not attributes:\n"
      for ix in missing_atts:
        att_msg += v1_frame.name + '[' + str(ix) + ']' + '\n'
      messages.append(att_msg)
    for m in messages:
      messagebox.showwarning("Error", m, parent=self)
    if messages != []:
      return
    for v1_frame in self.v1_frame_list:
      # Create the variable array
      cmd_str = 'python var_v1_create ' + v1_frame.name + ' '
      cmd_str += str(v1_frame.ix_min) + ' '
      cmd_str += str(v1_frame.ix_max)
      self.pipe.cmd_in(cmd_str)
      # Create the individual variables
      for i in range(v1_frame.ix_min, v1_frame.ix_max+1):
        cmd_str = 'python var_create ' + v1_frame.name
        cmd_str += '[' + str(i) + ']^'
        if i in v1_frame.var_dict.keys():
          var_dict = v1_frame.var_dict[i]
        else:
          var_dict = {}
        params = ['ele_name', 'attribute', 'universes', 'weight', 'step',
            'low_lim', 'high_lim', 'merit_type', 'good_user', 'key_bound',
            'key_delta']
        # Used to look up v1-level default parameters
        v1_params = ['name', 'universes', 'attribute', 'weight', 'step', 'merit_type',
            'low_lim', 'high_lim', 'good_user', 'key_bound', 'key_delta']
        for j in range(len(params)):
          p = params[j]
          if p in var_dict.keys():
            if p == 'universes': #cannot be empty
              u = v1_frame.v1_array_wids[v1_ix].tk_var.get()
              if u == "":
                pass
            else:
              cmd_str += str(var_dict[p]) + '^'
          elif j == 0:
            cmd_str += '^'
          else:
            v1_ix = v1_params.index(p)
            if v1_frame.v1_array_wids[v1_ix].param.type == 'LOGIC':
              cmd_str += ('T^' if v1_frame.v1_array_wids[v1_ix].tk_var.get() else 'F^')
            elif p == 'universes': #cannot be empty
              u = v1_frame.v1_array_wids[v1_ix].tk_var.get()
              if u == "":
                u = '*'
              cmd_str += u + '^'
            else:
              cmd_str += v1_frame.v1_array_wids[v1_ix].tk_var.get() + '^'
        cmd_str = cmd_str[:-1] # remove ending ^
        self.pipe.cmd_in(cmd_str)
      # Set parameters at the v1 level
      #set_str = 'set var ' + v1_frame.name + '|'
      #tao_set(v1_frame.v1_array_wids[1:10], set_str, self.pipe)
      ## set individual data parameters
      #set_format = set_str[:-1] + '[{}]|'
      #tao_dict_set(v1_frame.var_dict, set_format, self.pipe)
      ## set var|exists = T
      #self.pipe.cmd_in('set var ' + v1_frame.name + '|exists = T')
    # Close the window
    self.destroy()

  def tab_handler(self, event=None):
    '''
    Handles new tab creation and updates self.d1_index as necessary
    '''
    # Check if the new tab frame has been selected
    if self.notebook.select() == self.new_tab_frame._w:
      # Add new tab
      self.v1_frame_list.append(new_v1_frame(self))
      self.v1_index = len(self.v1_frame_list)-1
      self.notebook.insert(self.v1_index, self.v1_frame_list[-1])
      self.notebook.tab(self.d1_index, text='New v1_array')
      self.notebook.select(self.v1_index)
    else:
      # Update self.v1_index
      for i in range(len(self.v1_frame_list)):
        frame = self.v1_frame_list[i]
        if self.notebook.select() == frame._w:
          self.v1_index = i
          break


class new_v1_frame(tk.Frame):
  '''
  Provides a frame for inputting properties of a v1_variable array.
  '''
  def __init__(self, parent):
    tk.Frame.__init__(self, parent.notebook)
    self.parent = parent
    self.pipe = self.parent.pipe
    self.name = "New v1_array" #Default

    # d1 Widgets
    def v1_ttp(x):
      ''' Shortcut for commonly used construct '''
      return tk_tao_parameter(str_to_tao_param(x), self, self.pipe)

    self.v1_array_wids = [v1_ttp('name;STR;T;'),
        v1_ttp('universes;STR;T;'),
        v1_ttp('attribute;STR;T;'),
        v1_ttp('weight;REAL;T;'),
        v1_ttp('step;REAL;T;'),
        v1_ttp('var^merit_type;ENUM;T;'),
        v1_ttp('low_lim;REAL;T;'),
        v1_ttp('high_lim;REAL;T;'),
        v1_ttp('good_user;LOGIC;T;T'),
        v1_ttp('key_bound;LOGIC;T;'),
        v1_ttp('key_delta;REAL;T;'),
        v1_ttp('ix_min;INT;T;'),
        v1_ttp('ix_max;INT;T;')]
    # v1 labels (NAMES AS STRINGS ONLY)
    self.v1_array_labels = ["v1_array Name:", "Default universe:", "Default attribute:",
        "Default weight:", "Default step:", "Default merit type",
        "Default low_lim:", "Default high_lim", "Default good_user",
        "Default key_bound", "Default key_delta", "Start index", "End index"]
    # Grid widgets and labels:
    for i in range(len(self.v1_array_wids)):
      tk.Label(self, text=self.v1_array_labels[i]).grid(row=i+1, column=0, sticky='W')
      self.v1_array_wids[i].tk_wid.grid(row=i+1, column=1, sticky='EW')
    i = i+1

    # Warning labels
    # (defined here to be gridded/ungridded as necessary)
    self.name_warning_1 = tk.Label(self, text="Must not be empty")
    self.name_warning_2 = tk.Label(self, text="v1 name already in use")
    self.ix_min_warning_1 = tk.Label(self, text="Must be a non-negative integer")
    self.ix_min_warning_2 = tk.Label(self, text="Cannot be larger than the maximum index")
    self.ix_max_warning_1 = tk.Label(self, text="Must be a non-negative integer")
    self.ix_max_warning_2 = tk.Label(self, text="Cannot be smaller than the minimum index")
    self.low_high_warning_1 = tk.Label(self, text="Must be a real number")
    self.low_high_warning_2 = tk.Label(self, text="Must be a real number")
    self.low_high_warning_3 = tk.Label(self, text="Min value must be smaller than max value")

    # Responses to edits
    self.v1_array_wids[0].tk_wid.bind('<FocusOut>', self.name_handler)
    self.v1_array_wids[11].tk_wid.bind('<FocusOut>', self.ix_min_handler)
    self.v1_array_wids[12].tk_wid.bind('<FocusOut>', self.ix_max_handler)
    self.v1_array_wids[6].tk_wid.bind('<FocusOut>', self.low_high_handler)
    self.v1_array_wids[7].tk_wid.bind('<FocusOut>', self.low_high_handler)

    ttk.Separator(self, orient='horizontal').grid(row=i+1, column=0, columnspan=3, sticky='EW')
    i = i+1

    # Element browsers
    tk.Label(self, text="Choose elements:").grid(row=i+1, column=0, sticky='W')
    self.ele_name_button = tk.Button(self, text="Browse...",
        command=lambda : self.lat_browser())
    self.ele_name_button.grid(row=i+1, column=1, sticky='EW')

    ttk.Separator(self, orient='horizontal').grid(row=i+4, column=0, columnspan=3, sticky='EW')
    i = i+1

    # Individual variables
    self.var_dict = {}
    self.var_frame = tk.Frame(self)
    tk.Label(self, text="Variable:").grid(row=i+4, column=0, sticky='W')
    self.var_ix = tk.StringVar()
    self.ix_min = -1
    self.ix_max = -1
    self.var_chooser = ttk.Combobox(self, textvariable=self.var_ix,
        values=['PLEASE SPECIFY MIN/MAX INDICES'], state='readonly')
    self.var_chooser.bind('<<ComboboxSelected>>', self.make_var_frame)
    self.var_chooser.grid(row=i+4, column=1, sticky='EW')

    # Delete button
    tk.Button(self, text="DELETE THIS V1_ARRAY", fg='red', command=self.delete).grid(
        row=0, column=0, columnspan=3, sticky='EW')

    # Focus the v1 name widget
    self.v1_array_wids[0].tk_wid.focus_set()

  def delete(self, ask=True, event=None):
    '''
    Deletes this v1_array frame
    '''
    # Ask for confirmation
    ans = messagebox.askokcancel("Delete " + self.name, "Delete this v1_array and its associated data?", parent=self.parent)
    if not ans:
      return

    # Remove from parent
    ix = self.parent.v1_index
    if ix > 0:
      # Prefer to move to the tab to the left
      self.parent.notebook.select(ix-1)
    self.parent.v1_frame_list.pop(ix)

    # Destroy self
    self.destroy()

  def make_var_frame(self, event=None):
    '''
    Adds an entry to self.var_dict if necessary, and updates
    self.var_frame to show the values for the current variable
    '''
    ix = self.var_ix.get()
    if ix == 'PLEASE SPECIFY MIN/MAX INDICES':
      return
    ix = int(ix)
    if ix not in self.var_dict.keys():
      self.var_dict[ix] = {}
    self.var_frame.destroy()
    self.var_frame = tk.Frame(self)
    def var_ttp(x):
      ''' Shortcut for commonly used construct '''
      return tk_tao_parameter(x, self.var_frame, self.pipe)
    # Parameters
    param_list = ['ele_name;STR;T;', 'universes;STR;T;', 'attribute;STR;T;',
        'weight;REAL;T;', 'step;REAL;T;', 'var^merit_type;ENUM;T;',
        'low_lim;REAL;T;', 'high_lim;REAL;T;', 'good_user;LOGIC;T;T',
        'key_bound;LOGIC;T;', 'key_delta;REAL;T;']
    param_list = list(map(str_to_tao_param, param_list))
    # Fill in defaults set at v1 level
    for i in range(1, len(param_list)):
      param_list[i].value = self.v1_array_wids[i].tk_var.get()
    # Set parameter values if specified
    for i in range(len(param_list)):
      if param_list[i].name in self.var_dict[ix].keys():
        param_list[i].value = self.var_dict[ix][param_list[i].name]
    # Create widgets
    self.var_wid_list = list(map(var_ttp, param_list))
    def var_writer(i):
      ''' Writes the contents of self.var_wid_list[i] into self.data_dict '''
      if i in [6,7]:
        if low_high_handler(grid=False):
          return # Don't write to low/high limit if contents are invalid
      if self.var_wid_list[i].param.type == 'LOGIC':
        val = bool(self.var_wid_list[i].tk_var.get())
      else:
        val = self.var_wid_list[i].tk_var.get()
      self.var_dict[ix][self.var_wid_list[i].param.name] = val
    def var_writer_callback(i):
      ''' Callback for data_writer() '''
      return lambda *args : var_writer(i)

    # low/high validation
    def low_high_handler(self, grid=True, event=None):
      '''
      Verifies that the min and max values are set to good values
      returns 1 if they are not
      if grid is set False, messages will not be gridded/ungridded
      '''
      # Low limit
      fail = False
      if self.v1_array_wids[6].tk_var.get() != "":
        try:
          low_lim = float(self.v1_array_wids[6].tk_var.get())
        except ValueError:
          low_lim = None
        if low_lim == None:
          self.low_high_warning_1.grid(row=7, column=2, sticky='W')
          fail = True
        else:
          self.low_high_warning_1.grid_forget()
      else:
        low_lim = None
      # High limit
      if self.v1_array_wids[7].tk_var.get() != "":
        try:
          high_lim = float(self.v1_array_wids[7].tk_var.get())
        except ValueError:
          high_lim = None
        if high_lim == None:
          self.low_high_warning_2.grid(row=8, column=2, sticky='W')
          fail = True
        else:
          self.low_high_warning_2.grid_forget()
      else:
        high_lim = None
      # Make sure low < high
      if (low_lim != None) & (high_lim != None):
        if low_lim > high_lim:
          self.low_high_warning_3.grid(row=7, column=2, sticky='W')
          fail = True
        else:
          self.low_high_warning_3.grid_forget()
      else:
        self.low_high_warning_3.grid_forget()

      if fail:
        return 1

    # Grid widgets
    low_high_warning_1 = tk.Label(self.var_frame, text="Must be a real number")
    low_high_warning_2 = tk.Label(self.var_frame, text="Must be a real number")
    low_high_warning_3 = tk.Label(self.var_frame, text="Low limit must be less than high limit")
    for i in range(len(self.var_wid_list)):
      self.var_wid_list[i].tk_label.grid(row=i, column=0, sticky='E')
      self.var_wid_list[i].tk_wid.grid(row=i, column=1, sticky='EW')
      # Write changes into self.var_dict
      self.var_wid_list[i].tk_var.trace('w', var_writer_callback(i))
    self.var_wid_list[6].tk_wid.bind('<FocusOut>', low_high_handler)
    self.var_wid_list[7].tk_wid.bind('<FocusOut>', low_high_handler)
    self.var_frame.grid(row=20, column=0, columnspan=3, sticky='EW')


  def lat_browser(self):
    '''
    Opens a modified lattice table window to select elements en masse
    which should be either 'name', 'start_name', or 'ref_name'
    When the elements are chosen, ele_which will be set for each of this
    v1_array's variables sequentially
    '''
    # Make sure the name, ix_min, and ix_max are set
    if (bool(self.name_handler()) | bool(self.ix_min_handler())
        | bool(self.ix_max_handler()) | bool(self.low_high_handler())):
      return
    #
    win = tao_ele_browser(self.parent.root, self.pipe, self.name, self, 'var')

  def name_handler(self, event=None):
    '''
    Changes the tab name to match the v1_name
    Returns 1 if unsuccessful
    '''
    name = self.v1_array_wids[0].tk_var.get().strip()
    if name != "":
      # Make sure the name isn't already in use
      i = 0
      for v1 in self.parent.v1_frame_list:
        if (v1.name == name) & (self != self.parent.v1_frame_list[i]):
          self.name_warning_1.grid_forget()
          self.name_warning_2.grid(row=1, column=2, sticky='W')
          return 1
        i = i+1

      self.parent.notebook.tab(self.parent.v1_index, text=name)
      self.name_warning_1.grid_forget()
      self.name_warning_2.grid_forget()
      self.name = name
    else:
      self.name_warning_2.grid_forget()
      self.name_warning_1.grid(row=1, column=2, sticky='W')
      self.name = "New v1_array"
      self.parent.notebook.tab(self.parent.v1_index, text="New v1_array")
      return 1

  def ix_min_handler(self, event=None):
    try:
      ix_min = int(self.v1_array_wids[11].tk_var.get())
    except ValueError:
      ix_min = -1
    if ix_min > -1:
      self.ix_min = ix_min
      if (ix_min <= self.ix_max) | (self.ix_max == -1):
        self.ix_min_warning_1.grid_forget()
        self.ix_min_warning_2.grid_forget()
      else:
        self.ix_min_warning_1.grid_forget()
        self.ix_min_warning_2.grid(row=12, column=2, sticky='W')
        return 1
    else:
      self.ix_min_warning_2.grid_forget()
      self.ix_min_warning_1.grid(row=12, column=2, sticky='W')
      return 1
    # Update the data index range if possible
    if (self.ix_min > -1) & (self.ix_max >= self.ix_min):
      self.var_chooser.configure(values=list(range(self.ix_min, self.ix_max+1)))

  def ix_max_handler(self, event=None):
    if self.handler_block:
      return
    try:
      ix_max = int(self.v1_array_wids[12].tk_var.get())
    except ValueError:
      ix_max = -1
    if ix_max > -1:
      self.ix_max = ix_max
      if ix_max >= self.ix_min:
        self.ix_max_warning_1.grid_forget()
        self.ix_max_warning_2.grid_forget()
      else:
        self.ix_max_warning_1.grid_forget()
        self.ix_max_warning_2.grid(row=13, column=2, sticky='W')
        return 1
    else:
      self.ix_max_warning_2.grid_forget()
      self.ix_max_warning_1.grid(row=13, column=2, sticky='W')
      return 1
    # Update the variable index range if possible
    if (self.ix_min > -1) & (self.ix_max >= self.ix_min):
      self.var_chooser.configure(values=list(range(self.ix_min, self.ix_max+1)))

  def low_high_handler(self, event=None):
    '''
    Verifies that the min and max values are set to good values
    returns 1 if they are not
    '''
    # Low limit
    fail = False
    if self.v1_array_wids[6].tk_var.get() != "":
      try:
        low_lim = float(self.v1_array_wids[6].tk_var.get())
      except ValueError:
        low_lim = None
      if low_lim == None:
        self.low_high_warning_1.grid(row=7, column=2, sticky='W')
        fail = True
      else:
        self.low_high_warning_1.grid_forget()
    else:
      low_lim = None

    # High limit
    if self.v1_array_wids[7].tk_var.get() != "":
      try:
        high_lim = float(self.v1_array_wids[7].tk_var.get())
      except ValueError:
        high_lim = None
      if high_lim == None:
        self.low_high_warning_2.grid(row=8, column=2, sticky='W')
        fail = True
      else:
        self.low_high_warning_2.grid_forget()
    else:
      high_lim = None

    # Make sure low < high
    if (low_lim != None) & (high_lim != None):
      if low_lim > high_lim:
        self.low_high_warning_3.grid(row=7, column=2, sticky='W')
        fail = True
      else:
        self.low_high_warning_3.grid_forget()
    else:
      self.low_high_warning_3.grid_forget()

    if fail:
      return 1



class tao_ele_browser(tao_lattice_window):
  '''
  Very similar to the tao_lattice_window, but is
  intendend for selecting elements to be assigned
  to data as ele_name/ele_start_name/ele_ref_name
  root: the tk root window
  pipe: tao_interface instance
  name: the name of the d1/v1 array this is browsing for
  parent: the new_d1/v1_frame this is browsing for
  parent_type: must be either "data" or "var"
  which: must be "name", "start_name", or "ref_name"
  uni: a universe index as a string, or "All"
  '''
  def __init__(self, root, pipe, name, parent, parent_type,
      which='name', uni="All", *args, **kwargs):
    if uni == "All":
      uni = "1"
    self._full_init = False # used to delay self.refresh until init finishes
    tao_lattice_window.__init__(
        self, root, pipe, switches="-universe " + uni + " -all", *args, **kwargs)
    self.name = name
    self.title(self.name + ': ele_' + which + ' Browser')
    self.parent = parent
    self.parent_type = parent_type
    self.which = which

    # Remove template widgets
    self.temp_label_1.grid_forget()
    self.template_file.tk_wid.grid_forget()
    self.temp_label_2.grid_forget()
    self.temp_chooser.grid_forget()
    self.save_button.grid_forget()
    self.temp_save.tk_wid.grid_forget()

    # Default to custom ele list
    self.ele_list_opt.set('Custom')
    self.ele_list_box.focus_set()

    # Data index widgets
    self.data_frame = tk.Frame(self.top_frame)
    self.data_frame.grid(row=0, column=6)
    tk.Label(self.data_frame, text="Apply ele names to " + self.name + "[").pack(side="left")
    self.ix_min_var = tk.StringVar()
    tk.Entry(self.data_frame, width=5, textvariable=self.ix_min_var).pack(side="left")
    tk.Label(self.data_frame, text="] through " + self.name + "[").pack(side="left")
    self.ix_max_var = tk.StringVar()
    tk.Entry(self.data_frame, width=5, textvariable=self.ix_max_var).pack(side="left")
    tk.Label(self.data_frame, text="]").pack(side="left")
    # placeholder:
    self.ele_count = tk.Label(self.data_frame, text="")
    self.ix_min_var.set(self.parent.ix_min)
    self.ix_max_var.set(self.parent.ix_max)

    self.apply_button = tk.Button(self.top_frame, text="Apply Element Names", command=self.apply_callback)
    self.apply_button.grid(row=1, column=6)
    # Refresh now
    self._full_init = True
    self.refresh()

  def apply_callback(self, event=None):
    '''
    Checks that good indices have been input in self.ix_min_var
    and self.ix_max_var, and then copies ele names into ele_which
    for the selected data/var indices
    '''
    # Checks
    messages = []
    try:
      ix_min = int(self.ix_min_var.get())
    except ValueError:
      ix_min = -1
    if ix_min < 0:
      messages.append("Minimum index must be a non-negative integer")

    try:
      ix_max = int(self.ix_max_var.get())
    except ValueError:
      ix_max = -1
    if ix_max < 0:
      messages.append("Maximum index must be a non-negative integer")

    if ix_min > ix_max:
      messages.append("Minimum index must not be greater than maximum index")

    for message in messages:
      if self.parent_type == 'data':
        title = "Data index error"
      elif self.parent_type == 'var':
        title = "Variable index error"
      messagebox.showwarning(title, message, parent=self)
    if messages != []:
      return

    # Get ele names
    names = list(self.tree.get_children())
    for i in range(len(names)):
      if len(self.tree.item(names[i])['values']) < 2:
        # Remove "Lord Elements:" row
        names.pop(i)
        break
    for i in range(len(names)):
      names[i] = self.tree.item(names[i])['values'][1]
    # Stop if no eles found
    if len(names) == 0:
      return

    # Check that the number of names matches the number of indices
    if len(names) > len(range(ix_min, ix_max+1)):
      # Too many names
      message = "You have selected more elements than "
      if self.parent_type == 'data':
        message += "data indices.  "
      elif self.parent_type == 'var':
        message += 'variable indices.  '
      message += "Would you like to apply the first element names in the list to your "
      if self.parent_type == 'data':
        message += "data and discard the rest?"
      elif self.parent_type == 'var':
        message += "variables and discard the rest?"
      ans = messagebox.askokcancel("Too many elements", message, parent=self)
      if not ans:
        return
    elif len(names) < len(range(ix_min, ix_max+1)):
      # Not enough names
      message = "You have not selected enough element names for "
      message += self.name + '[{}]'.format(str(ix_min)) + ' through '
      message += self.name + '[{}]'.format(str(ix_max)) + '.  '
      message += "Would you like to repeat the names to fill "
      message += self.name + '[{}]'.format(str(ix_min)) + ' through '
      message += self.name + '[{}]'.format(str(ix_max)) + '?'
      ans = messagebox.askokcancel("Too few elements", message, parent=self)
      if not ans:
        return

    # Write the ele names into ele_which
    for i in range(ix_min, ix_max+1):
      if self.parent_type == 'data':
        if i not in self.parent.data_dict.keys():
          self.parent.data_dict[i] = {}
        self.parent.data_dict[i]['ele_'+self.which] = names[(i-ix_min)%len(names)]
      elif self.parent_type == 'var':
        if i not in self.parent.var_dict.keys():
          self.parent.var_dict[i] = {}
        self.parent.var_dict[i]['ele_'+self.which] = names[(i-ix_min)%len(names)]

    # Close this window
    self.destroy()

  def refresh(self, event=None):
    '''
    Overload of tao_lattice_window.refresh() that prints the number
    of elements found
    '''
    if not self._full_init:
      return
    # Call original version
    tao_lattice_window.refresh(self, event)
    # Count elements found
    try:
      names = list(self.tree.get_children())
      for i in range(len(names)):
        if len(self.tree.item(names[i])['values']) < 2:
          # Remove "Lord Elements:" row
          names.pop(i)
          break
    except tk._tkinter.TclError:
      names = [] # no lattice found
    # Display ele count
    self.ele_count.destroy()
    self.ele_count = tk.Label(self.data_frame,
        text=" (" + str(len(names)) + " elements found)")
    self.ele_count.pack(side='left')

class tao_progress_window(tk.Toplevel):
  '''
  Provides a window for displaying progress bars for long running processes
  Use the grid geometry manager for placing extra widgets in this window
  Row 0 is reserved for external manipulation (e.g. for a label)
  Properties:
  self.label_vars: list of variables used to set labels
  self.ix: counter used for external bookkeeping
  Init args:
  root: the Tk root window
  parent: the parent window for these progress bars
  num_bars: the number of progress bars
  Other init arguments are passed to the tk.Toplevel.__init__ method
  '''
  def __init__(self, root, parent, num_bars, *args, **kwargs):
    tk.Toplevel.__init__(self, parent, *args, **kwargs)
    self.parent = parent
    self._labels = [] # not intended to be externally manipulated
    self.label_vars = []
    self._bars = []
    self.ix = 0 #used for external bookkeeping
    for i in range(num_bars):
      self.label_vars.append(tk.StringVar())
      self._labels.append(tk.Label(self, textvariable=self.label_vars[-1]))
      self._labels[-1].grid(row=len(self._labels), column=0)
      self._bars.append(ttk.Progressbar(self, mode='determinate'))
      self._bars[-1].grid(row=len(self._bars), column=1)
  def set_max(self, ix, val):
    '''
    Sets the max value for self.bars[ix] to val
    '''
    self._bars[ix].configure(maximum=val)
  def set_val(self, ix, val, update=True):
    '''
    Sets to value of self.bars[ix] to val
    Note that self.set_max should be used first to set the progress bar max
    Use update=False to prevent running self.parent.update_idletaks()
    (usually necessary to make progress bar reflect actual progress)
    '''
    self._bars[ix].configure(value=val)
    if update:
      self.parent.update_idletasks()

class tao_message_box(tk.Toplevel):
  '''
  Custom messageboxes for when tkinter.messagebox does not suffice
  root: the Tk root window
  parent: the parent window for this message box
  tk_var: the tk.StringVar to write the results to
  title: title for this window
  message: the text to be displayed
  choices: list of strings for the user to pick from
  orient: grid the buttons horizontally or vertically (may be 'horizontal' or 'vertical')
  Other arguments are passed to the tk.Toplevel.__init__ method
  '''
  def __init__(self, root, parent, tk_var, title='', message='', choices=[],
      orient='horizontal', *args, **kwargs):
    tk.Toplevel.__init__(self, parent, *args, **kwargs)
    self.transient(parent)
    if title:
      self.title(title)
    self.tk_var = tk_var
    # grid buttons first to get appropriate window width
    i = 0
    for c in choices:
      if orient == 'horizontal':
        tk.Button(self, text=c, command=self.button_callback(c)).grid(
            row=1, column=i, sticky='EW', padx=5)
      elif orient == 'vertical':
        tk.Button(self, text=c, command=self.button_callback(c)).grid(
            row=i+1, column=0, sticky='EW', padx=5)
      i = i+1
    self.update_idletasks()
    tk.Label(self, text=message, wraplength=self.winfo_width()).grid(row=0, column=0, columnspan=i)
    self.protocol("WM_DELETE_WINDOW", self.cancel)
    self.wait_window(self)

  def button_callback(self, choice):
    return lambda: self.button_cmd(choice)

  def button_cmd(self, choice):
    self.tk_var.set(choice)
    self.destroy()

  def cancel(self, *args):
    self.tk_var.set("")













