# Check for required modules:
import sys
import os
from .module_check import module_check
module_check()

import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import font

from .tao_widget import *
from .tao_set import tao_set
from pytao.util.parameters import str_to_tao_param, tao_startup_param_dict
from pytao.util.parameters import tao_parameter_dict
from .tao_console import tao_console
from .tao_plot_dict import *
from .tao_interface import tao_interface

from .tao_data_windows import *
from .tao_lat_windows import *
from .tao_misc_windows import *
from .tao_plot_windows import *
from .tao_var_windows import *
from .tao_beam_windows import *

#---------------------------------------------------------------
# Root window

class tao_root_window(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, className='Tao')

        self.title("Tao")
        self.protocol("WM_DELETE_WINDOW", self.quit_cmd)
        self.tk.call('tk', 'scaling', 1.0)
        #default_font = font.nametofont("TkDefaultFont")
        #default_font.configure(size=14)
        #self.option_add("*Font", default_font)

        self.GUI_DIR = (os.environ['ACC_LOCAL_ROOT'] if
                'ACC_LOCAL_ROOT' in os.environ.keys() else os.environ['ACC_ROOT_DIR'])
        self.GUI_DIR += '/tao/python/pytao/gui'
        if sys.platform == "linux": # setting icon only works on linux for now
            self.icon = '@' + self.GUI_DIR + "/assets/tao-logo.xbm"
            self.iconbitmap(self.icon)

        # Window lists (accessible for refreshing)
        self.refresh_windows = {}
        self.tao_id_list = ['plot', 'data', 'var', 'ele', 'lat']
        for w in self.tao_id_list:
            self.refresh_windows[w] = []

        # Key bindings

        self.bind_all("<Control-q>", self.quit_cmd)

        # Init GUI

        self.old_init_args = ""
        init_frame = tk.Frame(self)
        init_frame.pack()
        self.tao_load(init_frame)

    #-----------------------------------------------------------------------------------
    # Tao startup

    def start_main(self):
        self.unbind("<Return>")
        self.lift()
        self.focus_force()
        #self.menubar.entryconfig("File", state="normal")
        #self.menubar.entryconfig("Window", state="normal")
        self.menubar_init()
        self.main_frame = tk.Frame(self, width = 20, height = 30)
        self.main_frame.pack()
        #Key bindings
        self.bind_all('<Control-n>', self.new_data_event)
        self.bind_all('<Alt-n>', self.new_var_event)
        self.bind_all('<Alt-p>', self.new_template_event)
        self.bind_all('<Control-g>', self.global_vars_event)
        self.bind_all('<Alt-g>', self.bmad_com_event)
        self.bind_all('<Control-d>', self.view_data_event)
        self.bind_all('<Control-e>', self.view_ele_event)
        self.bind_all('<Control-l>', self.view_lattice_event)
        self.bind_all('<Control-t>', self.plot_template_event)
        self.bind_all('<Control-r>', self.plot_region_event)
        self.bind_all('<Alt-v>', self.view_vars_event)
        self.bind_all('<Alt-q>', self.reinit_event)
        self.bind_all('<Control-b>', self.beam_init_event)
        # Clear out python lat_calc_done
        self.pipe.cmd_in('python lat_calc_done')

        # Call
        tk.Label(self.main_frame, text="Call command file:").grid(row=0, column=0)
        self.call_file = tk_tao_parameter(str_to_tao_param("call_file;FILE;T;"), self.main_frame, self.pipe)
        self.call_file.tk_wid.grid(row=0, column=1, sticky='EW')
        tk.Button(self.main_frame, text="Run", command=self.tao_call).grid(row=0, column=2, sticky='EW')
        # Arguments
        tk.Label(self.main_frame, text="Arguments:").grid(row=1, column=0)
        self.cf_args = tk_tao_parameter(str_to_tao_param("cf_args;STR;T;"), self.main_frame, self.pipe)
        self.cf_args.tk_wid.configure(width=30)
        self.cf_args.tk_wid.bind("<Return>", self.tao_call)
        self.cf_args.tk_wid.grid(row=1, column=1, columnspan=2, sticky='EW')

        # Optimization
        tk.Label(self.main_frame, text="Optimization:").grid(row=2, column=0)
        tk.Button(self.main_frame, text="Setup").grid(row=2, column=1, sticky='EW')
        tk.Button(self.main_frame, text="Run").grid(row=2, column=2, sticky='EW')

        # Debug
        tk.Button(self.main_frame, text="Debug", command=self.debug_cmd).grid(
                row=3, column=0, columnspan=3, sticky='EW')

        # Command line
        self.cmd_frame = tk.Frame(self)
        self.cmd_frame.pack(side="bottom", fill='both', expand=1)
        self.history = [] #holds the history
        self.history.append([]) #Tao and shell history
        self.history_pos = 0 #Used for scrolling in history on command line
        self.history.append([]) #Call history
        tk.Button(self.cmd_frame, text="View History...",
                command=self.view_history_cmd).grid(row=0, column=1, sticky='EW')
        tk.Button(self.cmd_frame, text="Show/Hide Console",
                command=self.sh_console).grid(row=0, column=0, sticky='EW')
        self.console = tao_console(self.cmd_frame, self, self.pipe)
        self.cmd_frame.columnconfigure(0, weight=1)
        self.cmd_frame.columnconfigure(1, weight=1)
        self.console.grid(row=1, column=0, columnspan=2, sticky = 'NSEW')
        self.cmd_frame.rowconfigure(1, weight=1)
        self.console_packed = True

    #-----------------------------------------------------------------------------------

    def sh_console(self, event=None):
        '''
        Packs or unpacks self.console as appropriate
        '''
        if self.console_packed:
            self.console.grid_forget()
        else:
            self.console.grid(row=1, column=0, columnspan=2, sticky = 'NSEW')
        self.console_packed = not self.console_packed

    #-----------------------------------------------------------------------------------

    def tao_call(self, event=None):
        '''
        Runs the command file in self.call_file, appends it to the history,
        and clears self.call_file
        '''
        if self.call_file.tk_var.get() != "Browse...":
            self.pipe.cmd("call " + self.call_file.tk_var.get() + ' '
                    + self.cf_args.tk_var.get())
            self.history[1].append(self.call_file.tk_var.get())
            self.call_file.tk_var.set("Browse...")
        #Try to refresh history window
        try:
            self.history_window.refresh()
        except:
            pass
        # Check the place buffer and place plots
        self.default_plots(include_init=False)
        # Update pgplot settings
        if self.plot_mode == "pgplot/plplot":
            self.placed.pgplot_update()

    #-----------------------------------------------------------------------------------

    def menubar_init(self):
        self.menubar = tk.Menu(self)

        file_menu = tk.Menu(self.menubar, tearoff=0)
        file_menu.add_command(label = 'Write Data Namelist...', command = self.write_data_cmd)
        file_menu.add_command(label = 'Write Var Namelist...', command = self.write_var_cmd)
        file_menu.add_separator()
        file_menu.add_command(label = 'Reinit...',
                command = self.reinit_cmd,accelerator = 'Alt+Q')
        file_menu.add_command(label = 'Quit', command = self.quit_cmd,
                accelerator = 'Ctrl+Q')
        self.menubar.add_cascade(label = 'File', menu = file_menu)

        create_menu = tk.Menu(self.menubar, tearoff=0)
        create_menu.add_command(label = "New Data...",
                command = self.new_data_cmd, accelerator = 'Ctrl+N')
        create_menu.add_command(label = "New Variables...",
                command = self.new_var_cmd, accelerator = 'Alt+N')
        create_menu.add_command(label = "New Plot Template...",
                command = self.new_template_cmd, accelerator = 'Alt+P')
        self.menubar.add_cascade(label = 'Create', menu = create_menu)

        view_menu = tk.Menu(self.menubar, tearoff=0)
        view_menu.add_command(label = 'Data...',
                command = self.view_data_cmd, accelerator = 'Ctrl+D')
        view_menu.add_command(label = 'Variables...',
                command = self.view_vars_cmd, accelerator = 'Alt+V')
        view_menu.add_command(label = 'Lattice Elements...',
                command = self.view_ele_cmd, accelerator = 'Ctrl+E')
        view_menu.add_command(label = 'Lattice...',
                command = self.view_lattice_cmd, accelerator = 'Ctrl+L')
        view_menu.add_command(label = 'Lat_layout Shapes...',
                command = self.lat_layout_shape_cmd)
        view_menu.add_command(label = 'Floor_plan Shapes...',
                command = self.floor_plan_shape_cmd)
        view_menu.add_command(label = 'Global Variables...',
                command = self.set_global_vars_cmd, accelerator = 'Ctrl+G')
        view_menu.add_command(label = 'Bmad Parameters...',
                command = self.set_bmad_com_cmd, accelerator = 'Alt+G')
        self.menubar.add_cascade(label = 'View', menu = view_menu)

        plot_menu = tk.Menu(self.menubar, tearoff=0)
        plot_menu.add_command(label = 'New Plot...',
                command = self.plot_template_cmd, accelerator = 'Ctrl+T')
        plot_menu.add_command(label = 'Edit Plot...',
                command = self.plot_region_cmd, accelerator = 'Ctrl+R')
        plot_menu.add_command(label = 'Building Wall Settings...',
                command = self.building_wall_cmd)
        self.menubar.add_cascade(label = 'Plot', menu = plot_menu)

        beam_menu = tk.Menu(self.menubar, tearoff=0)
        beam_menu.add_command(label = 'Beam settings...',
                command = self.beam_init_cmd, accelerator = 'Ctrl+B')
        self.menubar.add_cascade(label = 'Beam', menu = beam_menu)

        #window_menu = tk.Menu(self.menubar)
        #window_menu.add_command(label = 'Optimizer...',
        #        command = self.optimizer_cmd)
        #window_menu.add_command(label = 'Wave...', command = self.wave_cmd)
        #self.menubar.add_cascade(label = 'Window', menu = window_menu)

        self.config(menu=self.menubar)

    #-----------------------------------------------------------------------------------

    def default_plots(self, include_init=True):
        '''
        If self.plot_mode == "matplotlib", opens all the
        plot templates listed in plot.gui.init in
        separate matplotlib windows
        Also opens matplotlib windows for any plots placed by the tao init file
        Will only place plots from the place buffer if include_init is False
        '''
        if include_init:
            # Open plot.gui.init
            try:
                plot_file = open('plot.gui.init')
                plot_list = plot_file.read()
                plot_list = plot_list.splitlines()
                plot_file.close()
            except:
                plot_list = []
        else:
            plot_list = []
        # Read the place buffer
        init_plots = self.pipe.cmd_in('python place_buffer')
        init_plots = init_plots.splitlines()
        init_regions = [None]*len(init_plots)
        for i in range(len(init_plots)):
            init_regions[i] = init_plots[i].split(';')[0]
            init_plots[i] = init_plots[i].split(';')[1]
        plot_list = init_plots + plot_list
        plot_regions = init_regions + [None]*len(init_plots)

        # Get list of plot templates to check input against
        plot_templates = self.pipe.cmd_in("python plot_list t")
        plot_templates = plot_templates.splitlines()
        for i in range(len(plot_templates)):
            plot_templates[i] = plot_templates[i].split(';')[1]

        # Place the plots
        # Validate entries in plot.gui.init
        for i in range(len(plot_list)):
            plot = plot_list[i]
            if plot in plot_templates:
                if self.plot_mode == "matplotlib":
                    # Make a window for the plot
                    win = tao_plot_window(self, plot, self.pipe, region=plot_regions[i])
                elif self.plot_mode == "pgplot/plplot":
                    # self.placed (tao_plot_dict) will keep track of how many
                    # rows/columns of the pgplot window are available
                    self.placed.place_template(plot, plot_regions[i])

        # Determine pgplot layout from currently placed plots
        if self.plot_mode == "pgplot/plplot":
            self.placed.pgplot_update()


    #-----------------------------------------------------------------------------------

    def tao_load(self, init_frame):
        '''
        Handles the startup of tao, including parsing gui.init and commandline
        arguments.
        '''

        def toggle_init(*args):
            '''
            Enables/disables the init file box if noinit is selected
            '''
            ix_noinit = tk_list.index('noinit')
            ix_init_file = tk_list.index('init_file')
            if tk_list[ix_noinit].tk_var.get():
                tk_list[ix_init_file].tk_wid.configure(state='disabled')
            else:
                tk_list[ix_init_file].tk_wid.configure(state='normal')
        # End: def toggle_init

        def swap_box(event=None):
            if interface_to_tao.get() == "pexpect":
                ctype_label.grid_forget()
                tao_lib.tk_wid.grid_forget()
                pexp_label.grid(row=k+1, sticky="E")
                tao_exe.tk_wid.grid(row=k+1, column=1, sticky="W")
            elif interface_to_tao.get() == "ctypes":
                pexp_label.grid_forget()
                tao_exe.tk_wid.grid_forget()
                ctype_label.grid(row=k+1, sticky="E")
                tao_lib.tk_wid.grid(row=k+1, column=1, sticky="W")
        # End: def swap_box

        def param_load(event=None):
            # set the font
            try:
                new_size = int(font_var.get())
                new_font = font.nametofont("TkDefaultFont")
                new_font.configure(size=new_size)
                self.option_add("*Font", new_font)
                self.font_size = new_size
            except:
                self.font_size = None
            # set self.lw
            self.lw = lw_var.get()
            self.plot_mode = plot_mode.get()
            # With noinit, make sure lattice file is chosen
            ix_noinit = tk_list.index('noinit')
            ix_lat_file = tk_list.index('lattice_file')
            if tk_list[ix_noinit].tk_var.get() and tk_list[ix_lat_file].tk_var.get() in ["", "Browse..."]:
                messagebox.showwarning("Warning", "A lattice file is required when starting Tao without an init file.")
                return
            init_args = ""
            for tk_param in tk_list:
                if tk_param.param.type in ['STR','ENUM']:
                    tk_param.param.value = tk_param.tk_var.get()
                    if tk_param.param.value != "":
                        init_args = (init_args + "-" + tk_param.param.name + " "
                                + tk_param.param.value + " ")
                    elif tk_param.param.name in self.old_init_args:
                        init_args += '--' + tk_param.param.name + ' '
                elif tk_param.param.type == 'FILE':
                    tk_param.param.value = tk_param.tk_var.get()
                    if tk_param.param.value == "Browse...":
                        tk_param.param.value = ""
                    if tk_param.param.value != "":
                        init_args = (init_args + "-" + tk_param.param.name + " "
                                + tk_param.param.value + " ")
                    elif tk_param.param.name in self.old_init_args:
                        init_args += '--' + tk_param.param.name + ' '
                elif tk_param.param.type == 'LOGIC':
                    tk_param.param.value = tk_param.tk_var.get()
                    if tk_param.param.value:
                        init_args = init_args + "-" + tk_param.param.name + " "
                    elif tk_param.param.name in self.old_init_args:
                        init_args += '--' + tk_param.param.name + ' '
            #
            if plot_mode.get() != "pgplot/plplot":
                init_args = init_args + "-noplot -external_plotting"
            # Run Tao, clear the init_frame, and draw the main frame
            if interface_to_tao.get() == "pexpect":
                mode = "pexpect"
                if tao_exe.tk_var.get() == "Browse...":
                    tao_exe.tk_var.set("")
                self.pipe = tao_interface(mode, init_args, tao_exe.tk_var.get())
            elif interface_to_tao.get() == "ctypes":
                mode = "ctypes"
                if tao_lib.tk_var.get() == "Browse...":
                    tao_lib.tk_var.set("")
                #Ctypes needs a more sophisticated reinit using reinit tao
                #Check if Tao is open
                if 'pipe' in self.__dict__ and self.pipe.mode == 'ctypes':
                    self.pipe.cmd_in('reinit tao ' + init_args)
                else:
                    self.pipe = tao_interface(mode, init_args, so_lib=tao_lib.tk_var.get())
            self.old_init_args = init_args

            init_frame.destroy()
            self.start_main()
            # Clear the plotting lists
            self.placed = tao_plot_dict(self.plot_mode, self.pipe)
            if plot_mode.get() == "matplotlib":
                # place the lattice layout in layout by default
                self.placed.place_template('lat_layout', 'layout')
            self.default_plots()
        # End: def param_load

        def update_font(event=None):
            try:
                new_size = int(font_var.get())
                new_font = font.nametofont("TkDefaultFont")
                new_font.configure(size=new_size)
                self.option_add("*Font", new_font)
            except:
                pass
            return 'break' # Do not load into Tao yet

        # End of local defs
        #-----------------------

        tk_list = [] #Items: tk_tao_parameter() (see tao_widget.py)
        init_frame.grid_columnconfigure(0, weight=1, pad=10)
        init_frame.grid_columnconfigure(1, weight=1)
        warning_messages = []
        gui_init_file = None

        # Parse command line arguments
        clargs = {} # keys: command line switches, values: switch value
        looking_for = 'switch'
        for i in range(len(sys.argv)):
            if looking_for == 'switch':
                arg = sys.argv[i]
                if arg[0] == '-': # -switch
                    arg = arg[1:]

                    normal = True
                    if arg[0] == '-':   # "--switch" negates the effect of "-switch"
                      normal = False
                      arg = arg[1:]

                    # Determine if arg is a valid switch name
                    matches = [] # list of startup params that arg might refer to
                    # Special handling for gui_init
                    if arg == 'gui_init':
                        gui_init_file = ""
                        looking_for = 'file'
                        continue
                    # Generic case
                    for param in tao_startup_param_dict.keys():
                        if param.find(arg) == 0:
                            matches.append(param)
                    if len(matches) == 1: #exactly one match -> good switch
                        arg = matches[0]
                        if (tao_startup_param_dict[matches[0]].type == 'FILE') & (i != len(sys.argv)-1):
                            clargs[arg] = "" #add arg to clargs
                            if normal: looking_for = 'file' # file switches need files
                        elif tao_startup_param_dict[matches[0]].type == 'LOGIC':
                            if normal:
                              clargs[arg] = "T" #add arg to clargs
                            else:
                              clargs[arg] = "F" #add arg to clargs
            elif looking_for == 'file':
                arg_file = sys.argv[i]
                if arg_file.find('-') == 0: #not a file
                    if arg == 'gui_init':
                        pass
                    else:
                        # Remove previous argument from list
                        clargs.pop(arg)
                else:
                    if arg == 'gui_init':
                        gui_init_file = arg_file
                    else:
                        clargs[arg] = arg_file
                looking_for = 'switch'

        #Look for and read gui.init
        #each line should have the form
        #parameter:value
        #E.g. init_file:tao.init
        #E.g. rf_on:T
        if gui_init_file == None:
            # Try to use gui.init if init file not specified
            gui_init_file = os.path.expanduser('gui.init')
        else:
            gui_init_file = os.path.expanduser(gui_init_file)
            gui_init_file = os.path.expandvars(gui_init_file)

        try:
            init_file = open(gui_init_file)
            init_list = init_file.read()
            init_list = init_list.splitlines()
            init_file.close()
        except:
            init_list = []

        for key, value in clargs.items(): # add clargs to init_list
            init_list.append(key + ':' + value)

        init_dict = {}
        for entry0 in init_list:
            entry = entry0.strip()
            entry = entry.split('#', 1)[0]   # Remove comment
            if entry == '': continue
            #Check for proper formatting and not a comment
            if entry.find(':') != -1:
                name = entry[:entry.find(':')]
                name = name.strip()
                value = entry[entry.find(':')+1:]
                value = value.strip()
            else:
                messagebox.showwarning('Error!', 'Bad line in gui.init file:\n' + entry0)
                continue

            c1 = (name in ["tao_executable", "tao_shared_library"])
            try:
                c2 = (tao_startup_param_dict[name].type == 'FILE')
            except KeyError:
                c2 = False
            if c1 | c2:
                #Check that a good filename has been given
                try:
                    filename = value
                    #Expand environment variables and ~
                    filename = os.path.expanduser(filename)
                    filename = os.path.expandvars(filename)
                    f = open(filename)
                    f.close()
                    init_dict[name] = filename
                except FileNotFoundError:
                    #Don't put bad files in init_dict
                    #messagebox.showwarning(name, "File not found: " + filename)
                    warning_messages.append([name, filename])
            else:
                init_dict[name] = value

        k = 0 #row number counter
        for param, tao_param in tao_startup_param_dict.items():
            if param == 'noplot': continue   # The setting of noplot affects the "Plot Mode"
            tk_list.append(tk_tao_parameter(tao_param,init_frame))
            #Possibly set value from init file
            if tk_list[k].param.name in init_dict:
                if tk_list[k].param.type == 'FILE':
                    tk_list[k].tk_var.set(init_dict[tk_list[k].param.name])
                elif tk_list[k].param.type == 'LOGIC':
                    state = (init_dict[tk_list[k].param.name] == 'T') | (init_dict[tk_list[k].param.name] == 'True')
                    tk_list[k].tk_var.set(state)

            if param == 'noinit': #gets a special label
                tk.Label(init_frame, text="Don't use an init file (noinit)").grid(row=k,sticky="E")
                tk_list[k].tk_var.trace('w', toggle_init)
            else:
                tk.Label(init_frame, text=param).grid(row=k,sticky="E")

            tk_list[k].tk_wid.grid(row=k, column=1, sticky="W")
            k = k+1

        # Add separator
        tk.ttk.Separator(init_frame, orient="horizontal").grid(row=k, columnspan=2, sticky='ew')
        k = k+1

        # Choosing whether to use pexpect or ctypes must be handled separately
        tk.Label(init_frame, text="Interface to Tao:").grid(row=k, sticky="E")

        interface_to_tao = tk.StringVar()
        interface_to_tao.set("ctypes")
        if "interface_to_tao" in init_dict:
            if init_dict["interface_to_tao"] in ["pexpect", "ctypes"]:
                interface_to_tao.set(init_dict["interface_to_tao"])
        tk.OptionMenu(init_frame, interface_to_tao, "pexpect", "ctypes", command=swap_box).grid(row=k, column=1, sticky="W")
        pexp_label = tk.Label(init_frame, text="Tao Executable")
        ctype_label = tk.Label(init_frame, text="Tao Shared Library")
        tao_exe = tk_tao_parameter(str_to_tao_param("tao_executable;FILE;T;"), init_frame)
        if "tao_executable" in init_dict:
            tao_exe.tk_var.set(init_dict["tao_executable"])
        elif 'ACC_LOCAL_ROOT' in os.environ.keys():
            tao_exe.tk_var.set(os.environ['ACC_LOCAL_ROOT']+'/production/bin/tao')
        elif 'ACC_EXE' in os.environ.keys():
            tao_exe.tk_var.set(os.environ['ACC_EXE']+'/tao')
        tao_lib = tk_tao_parameter(str_to_tao_param("tao_shared_library;FILE;T;"), init_frame)
        if "tao_shared_library" in init_dict:
            tao_lib.tk_var.set(init_dict["tao_shared_library"])
        elif 'ACC_LOCAL_ROOT' in os.environ.keys() and \
                               os.path.isfile(os.environ['ACC_LOCAL_ROOT']+'/production/lib/libtao.so'):
            tao_lib.tk_var.set(os.environ['ACC_LOCAL_ROOT']+'/production/lib/libtao.so')
        elif 'ACC_ROOT_DIR' in os.environ.keys() and \
                                  os.path.isfile(os.environ['ACC_ROOT_DIR']+'/production/lib/libtao.so'):
            tao_lib.tk_var.set(os.environ['ACC_ROOT_DIR']+'/production/lib/libtao.so')
        swap_box()

        #Choosing plot mode must also be handled separately
        tk.Label(init_frame, text="Plot Mode").grid(row=k+2, sticky='E')
        plot_mode = tk.StringVar()
        plot_mode.set("matplotlib")
        plot_options = ["matplotlib", "pgplot/plplot", "none"]
        plot_chooser = tk.OptionMenu(init_frame, plot_mode, *plot_options)
        plot_chooser.grid(row=k+2, column=1, sticky='W')
        # Set plot_mode from init_dict if specified
        if "plot_mode" in init_dict:
            if init_dict["plot_mode"] == "pgplot":
                plot_mode.set("pgplot/plplot")
            elif init_dict["plot_mode"] in plot_options:
                plot_mode.set(init_dict["plot_mode"])
            else:
                messagebox.showwarning('Error', 'Bad setting of "plot_mode" in gui.init: ' + init_dict["plot_mode"])

        if 'noplot' in init_dict and init_dict['noplot'] == 'T': plot_mode.set('none')

        # Font size
        tk.Label(init_frame, text="Font size").grid(row=k+3, sticky='E')
        font_var = tk.StringVar()
        if 'font_size' in init_dict.keys():
            font_var.set(str(init_dict['font_size']))
        else:
            font_var.set('14')
        font_box = tk.Entry(init_frame, textvariable=font_var)

        update_font() # Set the font for the first time
        font_box.bind('<Return>', update_font)
        font_box.grid(row=k+3, column=1, sticky='W')

        # Light-weight switch
        lw_var = tk.BooleanVar()
        if 'lw' in init_dict:
            lw_var.set((init_dict['lw']=='T') | (init_dict['lw'] == 'True'))
        else:
            lw_var.set(True)
        # Start button
        load_b = tk.Button(init_frame, text="Start Tao", command=param_load)
        load_b.grid(row=k+5, columnspan=2)
        self.bind("<Return>", param_load)

        # Show warning messages about bad filenames from startup
        for warning in warning_messages:
            messagebox.showwarning(warning[0], "File not found: " + warning[1])

        #Start Tao immediately if start_tao is set true
        c1 = "start_tao" in init_dict
        c2 = len(warning_messages) == 0
        c3 = 'pipe' not in self.__dict__
        c4 = ('start_tao' in init_dict.keys()) and (init_dict["start_tao"] in ['T', 'True'])
        if c1 & c2 & c3 & c4: param_load()

        if ('do_mainloop' in init_dict.keys()) and (init_dict["do_mainloop"] in ['F', 'False']):
          self.do_mainloop = False
        else:
          self.do_mainloop = True

    #-----------------------------------------------------------------------------------
    # Menu bar callbacks

    def new_data_cmd(self):
        win = tao_new_data_window(self, self.pipe)

    def new_var_cmd(self):
        win = tao_new_var_window(self, self.pipe)

    def new_template_cmd(self):
        win = tao_new_plot_template_window(self, self.pipe)

    def write_data_cmd(self):
        '''
        Writes a fortran namelist file for all d2_arrays
        '''
        out_file = filedialog.asksaveasfilename(title="Save data namelist as...")
        if out_file not in [(), ""]:
            self.pipe.cmd_in('write namelist -data ' + out_file)

    def write_var_cmd(self):
        '''
        Writes a fortran namelist file for all v1_arrays
        '''
        out_file = filedialog.asksaveasfilename(title="Save var namelist as...")
        if out_file not in [(), ""]:
            self.pipe.cmd_in('write namelist -variable ' + out_file)

    def reinit_cmd(self):
        '''
        Quit Tao, destroy the main frame, and respawn the init frame
        '''
        result = messagebox.askquestion("Reinit",
                "This will close all windows and restart Tao.  Are you sure?")
        if result != 'yes':
            return

        for child in self.winfo_children():
            child.destroy()

        if self.pipe.mode == 'pexpect':
            # Hack to prevent pexpect from hanging:
            self.pipe.cmd_in(" spawn echo Tao\\>; quit", no_warn=True)
        init_frame = tk.Frame(self)
        init_frame.pack()
        self.tao_load(init_frame)

    def quit_cmd(self, event = ''):
        result = messagebox.askquestion("Quit", "Are You Sure?", icon='warning')
        if result == 'yes':
            for child in self.winfo_children():
                # better safe than sorry
                try:
                    child.destroy()
                except:
                    pass
            sys.exit(0)
        else:
            return

    def optimizer_cmd(self):
        print ('Optimizer called')
        win = tk.Toplevel(self)
        win.title('Optimizer')

    def beam_init_cmd(self):
        win = tao_beam_init_window(self, self.pipe)

    def plot_template_cmd(self):
        win = tao_place_plot_window(self, self.pipe)

    def plot_region_cmd(self):
        #win = tao_plot_tr_window(self, self.pipe, "R")
        win = tao_new_plot_template_window(self, self.pipe, None, 'T')

    def building_wall_cmd(self):
        win = tao_building_wall_window(self, self.pipe)

    def wave_cmd(self):
        print ('Wave called')

    def set_global_vars_cmd(self):
        win = tao_global_vars_window(self)

    def set_bmad_com_cmd(self):
        win = tao_bmad_com_window(self)

    def view_vars_cmd(self):
        win = tao_var_general_window(self, self.pipe)

    def view_data_cmd(self):
        win = tao_d2_data_window(self, self.pipe)

    def view_ele_cmd(self):
        win = tao_ele_window(self, self.pipe)

    def view_lattice_cmd(self):
        win = tao_lattice_window(self, self.pipe)

    def lat_layout_shape_cmd(self):
        win = tao_ele_shape_window(self, self.pipe, "lat_layout")

    def floor_plan_shape_cmd(self):
        win = tao_ele_shape_window(self, self.pipe, "floor_plan")

    # Other callbacks

    def new_data_event(self, event):
        self.new_data_cmd()

    def new_var_event(self, event):
        self.new_var_cmd()

    def new_template_event(self, event):
        self.new_template_cmd()

    def global_vars_event(self, event):
        self.set_global_vars_cmd()

    def bmad_com_event(self, event):
        self.set_bmad_com_cmd()

    def view_data_event(self, event):
        self.view_data_cmd()

    def view_ele_event(self, event):
        self.view_ele_cmd()

    def view_lattice_event(self, event):
        self.view_lattice_cmd()

    def plot_template_event(self, event):
        self.plot_template_cmd()

    def plot_region_event(self, event):
        self.plot_region_cmd()

    def view_vars_event(self, event):
        self.view_vars_cmd()

    def beam_init_event(self, event):
        self.beam_init_cmd()

    def reinit_event(self, event):
        self.reinit_cmd()

    def view_history_cmd(self):
        #Just focus the existing history window, if it exists
        try:
            self.history_window.lift()
            #self.history_window.force_focus()
        except:
            self.history_window = tao_history_window(self)

    def debug_cmd(self):
        win = tk.Toplevel(self)
        win.title('Debug')

        tk.Label(win, text='Commands to print to stdout:').grid(row=0, column=0, sticky='EW')
        var = tk.StringVar()
        opts = ['None', 'Omit querries (show, python var_general, etc)', 'All']
        i=1
        for o in opts:
            tk.Radiobutton(win, text=o, variable=var, value=str(i)).grid(row=i, column=0, sticky='EW')
            i += 1
        def set_debug():
            '''Command for the button'''
            if var.get() == '1':
                print('1')
                self.pipe.debug = False
                self.pipe.debug_patterns = []
            elif var.get() == '2':
                print(2)
                self.pipe.debug = True
                self.pipe.debug_patterns = ['set', 'place',
                        'python data_d2_create', 'python datum_create',
                        'python var_create', 'python var_v1_create',
                        'python plot_manage_plot',
                        'python plot_manage_curve',
                        'python plot_manage_graph',
                        'python lat_calc_done', 'python place_buffer'
                        'write', 'call']
            elif var.get() == '3':
                print(3)
                self.pipe.debug = True
                self.pipe.debug_patterns = 'All'
        tk.Button(win, text='Apply', command=set_debug).grid(row=i, column=0, sticky='EW')

#-----------------------------------------------------

if __name__ == "__main__":
    root = tao_root_window(sys.argv)
    root.mainloop()
