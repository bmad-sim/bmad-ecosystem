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

    # Poopulate self.list_frame
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
    self.geometry('1500x600')
    self.pipe = pipe
    list_rows = []
    # row counter: i-ix_lb
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
      j=j+1
    for i in range(ix_lb, ix_ub+1):
      list_rows.append(d1_data_list_entry(self.list_frame, d1_data_name, i, u_ix, self.pipe))
      j = 0 # column counter
      for widget in list_rows[i-ix_lb].tk_list:
        widget.grid(row=i-ix_lb+1, column=j)
        j = j+1
      tk.Button(self.list_frame, text="Edit...", command=self.open_datum_window_callback(d1_data_name, i, u_ix)).grid(row=i-ix_lb+1, column=j)

  def open_datum_window_callback(self, d1_data_name, d1_data_ix, u_ix):
    return lambda : self.open_datum_window(d1_data_name, d1_data_ix, u_ix)

  def open_datum_window(self, d1_data_name, d1_data_ix, u_ix):
    param_list = self.pipe.cmd_in("python data1 " + str(u_ix) + '@' + d1_data_name + '[' + str(d1_data_ix) + ']')
    param_list = param_list.splitlines()
    for i in range(len(param_list)):
      print(param_list[i])
      param_list[i]=str_to_tao_param(param_list[i])
    win = tao_parameter_window(None, d1_data_name + '[' + str(d1_data_ix) + ']', param_list, self.pipe)


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

    # Menu bar

    menubar = tk.Menu(self)

    file_menu = tk.Menu(menubar)
    file_menu.add_command(label = 'Read...', command = self.read_cmd)
    file_menu.add_command(label = 'Write...', command = self.write_cmd)
    file_menu.add_command(label = 'Reinit...', command = self.reinit_cmd)
    file_menu.add_separator()
    file_menu.add_command(label = 'Quit', command = self.quit_cmd, accelerator = 'Ctrl+Q')
    menubar.add_cascade(label = 'File', menu = file_menu)

    window_menu = tk.Menu(menubar)
    window_menu.add_command(label = 'Optimizer...', command = self.optimizer_cmd)
    window_menu.add_command(label = 'Plotting...', command = self.plotting_cmd)
    window_menu.add_command(label = 'Wave...', command = self.wave_cmd)
    window_menu.add_command(label = 'Global Variables...', command = self.set_global_vars_cmd)
    window_menu.add_command(label = 'Data...', command = self.view_data_cmd)
    menubar.add_cascade(label = 'Window', menu = window_menu)

    self.config(menu=menubar)

    # Init GUI

    init_frame = tk.Frame(self, width = 20, height = 30)
    init_frame.pack()
    self.tao_load(init_frame)

    # Key bindings

    self.bind_all("<Control-q>", self.quit_cmd)

  #-------------------------
  # Tao startup

  def start_main(self):
    main_frame = tk.Frame(self, width = 20, height = 30)
    main_frame.pack()
    self.bind_all('<Return>', self.return_event)
    self.bind_all('<Control-g>', self.global_vars_event)
    tao_command = tk.StringVar()
    tk.Entry(main_frame, textvariable=tao_command).pack()
    def run_command():
      print(tao_command.get())
      self.pipe.cmd(tao_command.get())

    cmd_button = tk.Button(main_frame, text = "Run Command", command=run_command)
    cmd_button.pack()

  def tao_load(self,init_frame):
    from parameters import param_dict
    tk_list = [] #Items: tk_tao_parameter()'s (see tao_widget.py)
    init_frame.grid_columnconfigure(0, weight=1, uniform="test")
    init_frame.grid_columnconfigure(1, weight=1, uniform="test")
    k = 0 #row number counter
    for param, tao_param in param_dict.items():
      tk_list.append(tk_tao_parameter(tao_param,init_frame))
      tk.Label(init_frame,text=param).grid(row=k,sticky="E")
      tk_list[k].tk_wid.grid(row=k, column=1, sticky="W")
      k = k+1

#    for param, tao_param in param_dict.items():
#      #create new entry in tk_dict and display it
#      if tao_param.type == 'STR':
#        tk_dict[param] = [tk.StringVar()]
#
#        tk_dict[param].append(tk.Label(init_frame,text=param))
#
#        tk_dict[param].append(tk.Entry(init_frame, textvariable=tk_dict[param][0]))
#      elif tao_param.type == 'FILE':
#        tk_dict[param] = [tk.StringVar()]
#        tk_dict[param][0].set("Browse...")
#
#        tk_dict[param].append(tk.Label(init_frame,text=param))
#
#        tk_dict[param].append(tk.Button(init_frame,width=20,justify='left',wraplength=150,textvariable=tk_dict[param][0],command=self.open_file_callback(tk_dict,param) ))
#      elif tao_param.type == 'LOGIC':
#        tk_dict[param] = [tk.BooleanVar()]
#
#        tk_dict[param].append(tk.Label(init_frame,text=param))
#
#        tk_dict[param].append(tk.Checkbutton(init_frame,variable=tk_dict[param][0]))
#        tk_dict[param][0].set(param_dict[param].value)

    def param_load(event=None):
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
      self.pipe = tao_interface(init_args)
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
    global_list.pop(-1) # last line is blank
    global_list.pop(-1)
    for i in range(len(global_list)):
      global_list[i]=str_to_tao_param(global_list[i])
    win = tao_parameter_window(None, "Global Variables", global_list, self.pipe)

    b = tk.Button(win.button_frame, text="Set Global Variables", command=lambda : self.tao_set(win.tao_list))
    b.pack()

  def view_data_cmd(self):
    win = tao_d2_data_window(None, self.pipe)

  # Other callbacks

  def return_event(self, event):
    print("You hit return.")

  def global_vars_event(self, event):
    self.set_global_vars_cmd()

  def tao_set(self,tao_list):
    '''
    Takes a list of tk_tao_parameters and makes a call to tao
    to set the parameters to the values input by the user
    '''
    # STOP lattice calculation here
    self.pipe.cmd("set global lattice_calc_on = F")
    self.pipe.cmd("set global plot_on = F")
    #Freeze input fields:
    for item in tao_list:
      item.tk_wid.config(state="disabled")
    for item in tao_list:
      #Type casting and validation
      if item.param.type == 'INT':
        try:
          item.param.value = int(item.tk_var.get())
        except ValueError:
          messagebox.showwarning("Error",item.param.name + " must be an integer ")
      elif item.param.type == 'REAL':
        try:
          item.param.value = float(item.tk_var.get())
        except ValueError:
          messagebox.showwarning("Error",item.param.name + " must be a real number")
      else:
        item.param.value = item.tk_var.get()

      #Wait till the end to set lattice_calc_on and plot_on
      if item.param.name == 'lattice_calc_on':
        lattice_calc_on = str(item.param.value)
      elif item.param.name == 'plot_on':
        plot_on = str(item.param.value)
      else:
        #self.pipe.cmd("set global " + item.param.name + " = " + str(item.param.value))
        msg = self.pipe.cmd_in("set global " + item.param.name + " = " + str(item.param.value))
        if msg.find("ERROR") != -1:
          messagebox.showwarning(item.param.name,msg)
    #Now set lattice_calc_on and plot_on
    self.pipe.cmd("set global plot_on = " + plot_on)
    self.pipe.cmd("set global lattice_calc_on = " + lattice_calc_on)
    #Re-enable input
    for item in tao_list:
      item.tk_wid.configure(state="normal")

#---------------------------------------------------------------

if __name__ == "__main__":
  root = tao_root_window(sys.argv)
  root.mainloop()
