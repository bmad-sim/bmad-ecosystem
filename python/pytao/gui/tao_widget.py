'''
This module defines some basic constructs used by the GUI, including the
tk_tao_parameter class, which handles displaying and editing a single
parameter in tao
'''
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import sys
import os

from pytao.util.parameters import tao_parameter, str_to_tao_param, tao_parameter_dict
from .data_type_list import data_type_list
#from .tao_data_windows import tao_new_data_window, tao_d1_data_window

class tk_tao_parameter():
    '''
    Takes a tao_parameter (defined in parameters.py) and a tk frame,
    and creates an object containing the parameter and appropriate tk widget(s)
    for displaying and modifying the parameter and value
    pipe: the tao_interface object
    data_source: for DAT_TYPE and DAT_TYPE_Z, filters allowed data types
    plot: for DAT_TYPE_Z, the plot where x_axis_type should be checked
    struct_name: for the components of a STRUCT, pass the struct name here
    '''

    def __init__(self, tao_parameter, frame, pipe=0, data_source='', plot='', struct_name=''):
        self.param = tao_parameter
        self.pipe = pipe
        self.sub_wid = None # to be set externally, e.g. by tk_tao_linker
        if struct_name != '':
            struct_name = struct_name + '.'

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

        if self.param.is_ignored:
            self.tk_var = tk.StringVar()
            self.tk_var.set(str(self.param.value))
            self.tk_wid = tk.Label(frame, textvariable=self.tk_var)
        elif self.param.type in ['STR', 'INT', 'REAL', 'SPECIES']:
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
                if self.param.prefix == None:
                    options = enum_fetch(self.param.name,pipe)
                else:
                    options = enum_fetch(self.param.prefix + '^' + self.param.name,pipe)
            elif self.param.type == 'ENUM_Z': #added to handle DAT_TYPE_Z
                options = enum_fetch('data_type_z', pipe)
                self.param.type = 'ENUM'
            if options == [""]:
                options = [self.param.value]
            if options[0].find('[FATAL') == 0:
                options = ['[ERROR]']
            if self.param.value == "":
                self.tk_var.set(options[0])
            if self.param.can_vary:
                self.tk_wid = tk.OptionMenu(frame, self.tk_var, *options)
            else:
                self.tk_wid = tk.Label(frame, text=self.param.value)
            # Check for and remove num^ from self.param.name
            self.param.name = self.param.name.split('^')[-1]
        elif self.param.type == 'INUM':
            self.tk_var = tk.StringVar()
            self.tk_var.set(self.param.value)
            if self.param.prefix == None:
                options = inum_fetch(self.param.name,pipe)
            else:
                options = inum_fetch(self.param.prefix + '^' + self.param.name,pipe)
            if options == [""]:
                options = [self.param.value]
            if self.param.value == None:
                self.tk_var.set(options[0])
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
        elif self.param.type == 'DAT_TYPE_E':
            self.tk_var = tk.StringVar()
            self.tk_var.set(self.param.value)
            self.tk_wid = tk.Frame(frame)
            # Sub-widget for the DAT_TYPE portion
            # Parse out the extra info around the DAT_TYPE
            # dat_type_str == part of self.param.value relevant to the dat_type
            dat_type_str = self.param.value
            if self.param.value.find('[') == 0:
                if self.param.value.find('@') != -1:
                    dat_type_str = self.param.value[self.param.value.find('@')+1:]
            if dat_type_str.find('[') != -1:
                dat_type_str = dat_type_str[:dat_type_str.find('[')]
            dat_type_str = self.param.name + ';DAT_TYPE;' + ('T;' if self.param.can_vary else 'F;') + dat_type_str
            self._dat_type = tk_tao_parameter(
                    str_to_tao_param(dat_type_str),
                    self.tk_wid, self.pipe, data_source=data_source)
            # DAT_TYPE_E gets a few extra input fields
            self._uvar = tk.StringVar() # Universe range
            self._u = tk.Entry(self.tk_wid, textvariable=self._uvar)
            self._u.configure(width=10)
            self._evar = tk.StringVar() # Element index
            self._e = tk.Entry(self.tk_wid, textvariable=self._evar)
            self._e.configure(width=10)
            self._cvar = tk.StringVar() # Component
            self._c = tk.Entry(self.tk_wid, textvariable=self._cvar)
            self._c.configure(width=10)
            # Extra labels to put around uni, ele, component
            self._labels = [ tk.Label(self.tk_wid, text='['),
                    tk.Label(self.tk_wid, text=']@'),
                    tk.Label(self.tk_wid, text='['),
                    tk.Label(self.tk_wid, text=']|')]
            # Pack everything
            self._labels[0].pack(side='left')
            self._u.pack(side='left')
            self._labels[1].pack(side='left')
            self._dat_type.tk_wid.pack(side='left', fill='both', expand=1)
            self._labels[2].pack(side='left')
            self._e.pack(side='left')
            self._labels[3].pack(side='left')
            self._c.pack(side='left')
        elif self.param.type == 'STRUCT':
            self.tk_var = tk.StringVar()
            self.tk_var.set('STRUCT')
            self.tk_wid = tk.Frame(frame)
            self.tk_wid.grid_columnconfigure(1, weight=1)
            self._shown = False
            self._m_label = tk.StringVar()
            self._m_label.set("Configure..." if self.param.can_vary else "View...")
            self._m = tk.Button(self.tk_wid, textvariable=self._m_label, command=self._show_hide_struct)
            self._m.grid(row=0, column=0, columnspan=2, sticky='EW')
            self._s = [] # list of tk_tao_parameters
            for component in self.param.value:
                # SPECIAL CASE: x, x2, y, y2 structs
                # do not transfer their name as the prefix
                if self.param.name not in ['x', 'x2', 'y', 'y2']:
                    self._s.append(tk_tao_parameter(component, self.tk_wid, pipe, struct_name=self.param.name))
                else:
                    self._s.append(tk_tao_parameter(component, self.tk_wid, pipe))
        elif self.param.type == 'COMPONENT':
            self.tk_var = tk.StringVar()
            self.tk_var.set(self.param.value)
            if not self.param.can_vary:
                self.tk_wid = tk.Entry(frame, textvariable=self.tk_var)
            else:
                self.tk_wid = tk.Frame(frame)
                self._copts = ['model', 'base', 'design', 'meas', 'ref']
                self._aopts = ['', '+', '-']
                self._handle_block = False
                self.tk_var.trace('w', self._write_comps)
                self._c_refresh()

        if self.param.type not in ['DAT_TYPE', 'DAT_TYPE_E', 'REAL_ARR', 'STRUCT', 'COMPONENT']:
            if (self.param.type != 'FILE') and (not self.param.is_ignored):
                self.tk_wid.config(disabledforeground="black")
        elif self.param.type not in ['DAT_TYPE_E', 'STRUCT']:
            for widget in self._s:
                widget.config(disabledforeground="black")
        elif self.param.type == 'DAT_TYPE_E':
            self._u.config(disabledforeground="black")
            self._e.config(disabledforeground="black")
            self._c.config(disabledforeground="black")
        self.tk_label = tk.Label(frame, text=self.param.name)
        if not self.param.can_vary and not self.param.is_ignored:
            if self.param.type not in ['DAT_TYPE', 'DAT_TYPE_E', 'REAL_ARR', 'STRUCT']:
                self.tk_wid.config(state="disabled")
            elif self.param.type not in ['STRUCT', 'DAT_TYPE_E']:
                for widget in self._s:
                    widget.config(state="disabled")
            elif self.param.type == 'DAT_TYPE_E':
                self._u.config(state='disabled')
                self._e.config(state='disabled')
                self._c.config(state='disabled')
        #Bind info printing
        if pipe != 0:
            if pipe.debug == True:
                self.tk_wid.bind("<Button-3>", self.print_info)

    def print_info(self, *args):
        '''
        Prints diagnostic info about this widget for testing purposes
        '''
        print("Param name: " + self.param.name)
        print("Param type: " + self.param.type)
        print("Param can_vary: " + str(self.param.can_vary))
        print("Param value: " + str(self.param.value))
        print("Widget contents: " + str(self.tk_var.get()))
        print("---------------------")


    def value(self):
        '''
        Returns the value in the input field(s) of self, appropriately typed
        If an invalid value is input, returns None
        '''
        if self.param.type in ['STR', 'ENUM', 'FILE', 'SPECIES']:
            return self.tk_var.get()
        elif self.param.type in ['INT', 'INUM']:
            try:
                return int(self.tk_var.get())
            except:
                return None
        elif self.param.type == 'REAL':
            try:
                return float(self.tk_var.get())
            except:
                return None
        elif self.param.type == 'LOGIC':
            return bool(self.tk_var.get())
        elif self.param.type == 'DAT_TYPE':
            if self._is_valid_dat_type(self.tk_var.get()):
                return self.tk_var.get()
            else:
                return None
        elif self.param.type == 'REAL_ARR':
            output = []
            for i in range(len(self._svar)):
                try:
                    output.append(float(self._svar[i].get()))
                except:
                    output.append(None)
            return output
        elif self.param.type == 'STRUCT':
            d = {}
            for ttp in self._s:
                if ttp.value() != None:
                    d[ttp.param.name] = ttp.value()
            return d

    def copy(self, ttp):
        '''
        Copy the the value(s) of the tk_var(s) of ttp into self
        Does nothing if ttp.param.type != self.param.type
        '''
        if not isinstance(ttp, tk_tao_parameter):
            return
        if self.param.type != ttp.param.type:
            return
        if self.param.type == 'REAL_ARR':
            if len(ttp._svar) != len(self._svar):
                return
            for i in range(len(self._svar)):
                if ttp.value()[i] != None:
                    self._svar[i].set(ttp.value()[i])
        elif self.param.type == 'STRUCT':
            for i in range(len(self._s)):
                if self._s[i].param.name in ttp.value().keys():
                    self._s[i].tk_var.set(ttp.value()[self._s[i].param.name])

    def param_copy(self, tao_param):
        '''
        Copies tao_param.value into self.tk_var (and slaves if appropriate)
        '''
        if not isinstance(tao_param, tao_parameter):
            return
        if self.param.type in ['STR', 'INT', 'REAL', 'SPECIES']:
            if tao_param.value == None:
                self.tk_var.set("")
            else:
                self.tk_var.set(str(tao_param.value))
        elif self.param.type in ['ENUM', 'ENUM_Z']:
            self.tk_var.set(tao_param.value)
        elif self.param.type == 'INUM':
            self.tk_var.set(tao_param.value)
        elif self.param.type == 'FILE':
            self.tk_var.set(tao_param.value)
            if self.tk_var.get() == "":
                self.tk_var.set("Browse...")
        elif self.param.type == 'LOGIC':
            self.tk_var.set(tao_param.value)
            if not self.param.can_vary:
                self.tk_wid.configure(text=str(self.tk_var.get()))
        elif self.param.type == 'REAL_ARR':
            self.tk_var = tk.StringVar()
            val = ""
            for i in range(len(tao_param.value)):
                val = val+str(tao_param.value[i])
                val = val+';'
            self.tk_var.set(val) #; delimited list of values
            # Hack to prevent _update_real_arr from running:
            self.param.type = '~REAL_ARR'
            for i in range(len(tao_param.value)):
                self._svar[i].set(str(tao_param.value[i]))
            # Fix self.param.type
            self.param.type = 'REAL_ARR'
        elif self.param.type == 'DAT_TYPE':
            self.tk_var.set(tao_param.value)
            if self._data_source in ["data", "var"]:
                self._mvar.set((self.tk_var.get()))
                self._mvar_old = self._mvar.get() # Tracks changes in self._mvar
                self._s_refresh()
            else:
                self._mvar.set((self.tk_var.get()).split('.')[0])
                self._mvar_old = self._mvar.get() # Tracks changes in self._mvar
                self._s_refresh()
        elif self.param.type == 'STRUCT':
            for i in range(len(tao_param.value)):
                for j in range(len(self.param.value)):
                    # tao_param.value and self.param.value may not be in the same order
                    if self.param.value[j].name == tao_param.value[i].name:
                        self._s[j].param_copy(tao_param.value[i])
                        break

    def _show_hide_struct(self, event=None, *args):
        '''
        Shows or hides a struct's components as appropriate
        '''
        # Safeguard
        if self.param.type != 'STRUCT':
            return
        if self._shown:
            for ttp in self._s:
                ttp.tk_wid.grid_forget()
                ttp.tk_label.grid_forget()
            self._shown = False
            self._m_label.set("Configure..." if self.param.can_vary else "View...")
        else:
            for i in range(len(self._s)):
                ttp = self._s[i]
                ttp.tk_label.grid(row=i+1, column=0, sticky='W')
                ttp.tk_wid.grid(row=i+1, column=1, sticky='EW')
            self._shown = True
            self._m_label.set("Collapse...")

    def _update_real_arr(self, event=None, *args, **kwargs):
        '''
        Takes a real array's slave variables and puts their data
        into the master tk_var
        '''
        # Safeguard
        if self.param.type != 'REAL_ARR':
            return
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
        Returns existing d1_data arrays if self._data_source == 'data'
        and existing v1_var arrays if self._data_source == 'var'
        '''
        master_list = []
        # Don't filter if self._data_source = ""
        filter_data_source = False
        if self._data_source not in ['data', 'var']:
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
        elif self._data_source == 'data':
            # Fetch d2 arrays
            d2_arr = self.pipe.cmd_in('python data_d2_array 1').splitlines()
            for d2 in d2_arr:
                d1_list = self.pipe.cmd_in('python data_d1_array ' + d2).splitlines()
                for line in d1_list:
                    master_list.append(d2 + '.' + line.split(';')[3])
        else:
            # Fetch v1 arrays
            var_gen = self.pipe.cmd_in('python var_general').splitlines()
            for line in var_gen:
                master_list.append(line.split(';')[0])
        return master_list

    def _s_refresh(self, event=None, dat_source_swap=False, *args):
        '''
        Clears the existing slave widgets and variables,
        makes new slave widgets and variables, and populates them
        if necessary
        '''
        # Safeguard in case this method gets called for a normal tk_tao_param
        if self.param.type not in ['DAT_TYPE', 'DAT_TYPE_E']:
            return
        if self.param.type == 'DAT_TYPE_E':
            return self._dat_type._s_refresh()
        # Clear the existing slaves
        for item in self._s:
            item.destroy()
        self._svar = []
        self._stype = []
        self._s = []
        self._m.configure(values=self._get_dat_types(True))
        if dat_source_swap:
            if self._data_source in ['data', 'var']:
                self._mvar.set(self.tk_var.get())
            else:
                if self._mvar.get().find('normal.h') == 0:
                    self._mvar.set('normal.h')
                else:
                    self._mvar.set(self._mvar.get().split('.')[0])
            self._mvar_old = self._mvar.get()

        try:
            if self._data_source in ['data' , 'var']:
                m_ix = 0
            else:
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
        #        current_mvar = current_mvar[:-1]
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
        #print("_update_tk_var called for " + self.param.name)
        #print("param.type = " + self.param.type)
        #print("self._data_soure = " + self._data_source)
        if self.param.type == 'DAT_TYPE_E':
            self._dat_type._update_tk_var()
            self.tk_var.set('['+self._uvar.get()+']@'
                    +self._dat_type.tk_var.get()
                    +'['+self._evar.get()+']|'
                    +self._cvar.get())
        if self.param.type != 'DAT_TYPE':
            return
        #if self._data_source not in ['lat', 'beam']:
        #    return
        new_tk_var = self._mvar.get()
        #print("new_tk_var = " + new_tk_var)
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
            #print("new_tk_var = " + new_tk_var)
        # Special case: velocity. -> velocity
        if new_tk_var == "velocity.":
            new_tk_var = "velocity"
        # Un-trace tk_var to prevent repeatedly running this method
        self._no_s_refresh = True
        #print("new_tk_var = " + new_tk_var)
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
            return False

        # Check each subvalue
        k = 0 #loop counter
        for p in slave_params:
            if p.find('<enum') != -1: #Enums
                if x[k+1] not in dat_dict[p]:
                    return False
            elif p.find('<digit:') != -1: #Digit dropdown box
                try:
                    x[k+1] = int(x[k+1])
                except ValueError:
                    return False
                p = p.split(':')[1]
                p = p.split('>')[0] #p is now "low-high"
                low, high = p.split('-')
                low = int(low)
                high = int(high)
                if (x[k+1] < low) | (x[k+1] > high):
                    return False
            elif p.find('<str>') != -1: #Strings
                if len(x[k+1]) == 0:
                    return False
            elif p.find('<digits') != -1: #Fixed length int
                try:
                    junk_var = int(x[k+1]) #don't need the value
                except ValueError:
                    return False
                p = p.split('s')[1]
                p = p.split('>')[0]
                length = int(p)
                if len(x[k+1]) != length:
                    return False
            elif p.find('<int>') != -1: #Integer
                try:
                    junk_var = int(x[k+1]) #don't need the value
                except ValueError:
                    return False
            elif p.find('<real>') != -1: #Float
                try:
                    junk_var = float(x[k+1]) #don't need the value
                except ValueError:
                    return False
            k = k+1
        return True # All tests passed

    def _fill_widgets(self, event=None, *args):
        '''
        Runs self._s_refresh() and then fills the slave widgets
        with data from self.tk_var, but only if self.is_valid_dat_type(self.tk_var.get()) is True
        '''
        if self.param.type != 'DAT_TYPE':
            return
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
        if self.param.type == 'DAT_TYPE_E':
            return self._dat_type._has_ele()
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
        if self.param.type == 'DAT_TYPE_E':
            return self._dat_type._has_s_offset()
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

    def _c_refresh(self):
        '''
        Reads the contents of self.tk_var and creates the proper number of
        OptionMenu widgets
        '''
        if self.param.type != 'COMPONENT':
            return
        for child in self.tk_wid.winfo_children():
            child.destroy()
        # parse self.tk_var
        val = self.tk_var.get()
        def parse_comps(x):
            '''
            Recursive function to read off one term in the component string.
            Returns a list in the form [comp1, +/-, comp2, +/-, ..., compN]
            '''
            x = x.strip()
            pix = x.find('+')
            if pix == -1:
                pix = len(x)
            mix = x.find('-')
            if mix == -1:
                mix = len(x)
            if pix < mix:
                a = '+'
            elif mix < pix:
                a = '-'
            else:
                a = None
            if a == None:
                return [x]
            else:
                return [x[:x.find(a)]] + [a] + parse_comps(x[x.find(a)+1:])
        val = parse_comps(val) #len(val) is odd
        # create widgets
        self._s = [] # widgets
        self._svar = [] # tk variables
        for i in range(len(val)+1):
            if i%2 == 0: #component
                opts = self._copts
            else: #+/-
                opts = self._aopts
            self._svar.append(tk.StringVar())
            self._svar[-1].set(val[i] if i<len(val) else "")
            self._s.append(tk.OptionMenu(self.tk_wid, self._svar[-1], *opts))
            self._s[-1].pack(side='left')
            self._svar[-1].trace('w', self._write_comps)

    def _write_comps(self, *args):
        '''
        Writes the contents of self._svar into self.tk_var (for COMPONENT
        parameters only), then calls self._c_refresh()
        '''
        if self._handle_block == True: # prevent this function from triggering itself
            return
        if self.param.type != 'COMPONENT':
            return
        val = ""
        for v in self._svar:
            if v.get() != "":
                val += v.get()
            else:
                break
        self._handle_block = True
        self.tk_var.set(val)
        #print(self.param.name + 'set to ' + val)
        self._handle_block = False
        self._c_refresh()

    def open_file(self):
        '''
        Provides a dialog window for selecting a file
        '''
        filename = filedialog.askopenfilename(title = "Select " + self.param.name)
        #Only set if the user has selected an actual file
        if isinstance(filename, str) & (filename != ""):
            self.tk_var.set(filename)
        else:
            self.tk_var.set("Browse...")

    # This makes it easier to search lists of tk_tao_parameters
    def __eq__(self, name):
        return self.param.name == name

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
        # Special case: ele shapes
        # Have to add in custom shapes
        if enum == "shape^shape":
            custom_shapes = pipe.cmd_in("python shape_pattern_list")
            custom_shapes = custom_shapes.splitlines()
            for i in range(len(custom_shapes)):
                custom_shapes[i] = custom_shapes[i].split(';')[0].strip()
            # Insert in custom pattern slot if possible
            if "Pattern:<pattern-name>" in option_list:
                ix = option_list.index("Pattern:<pattern-name>")
                option_list = option_list[:ix] + custom_shapes + option_list[ix+1:]
            else:
                option_list = option_list + custom_shapes
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

def tk_tao_linker(tk_tao_list):
    '''
    Takes a list of tk_tao_parameters as input and returns a list of tk_tao_parameters
    The output list will have parameters marked with I either removed or
    added to their appropriate parent parameter's tk_wid
    '''
    # determined from ttp.is_ignored:
    sub_widget_names = [] # the names
    sub_widgets = [] # the widgets
    sub_widget_indices = [] # the indices in tk_tao_list
    # First pass: collect sub_widgets
    for i in range(len(tk_tao_list)):
        ttp = tk_tao_list[i]
        if ttp.param.is_ignored:
            # add parentheses around value
            ttp.tk_var.set('(' + ttp.tk_var.get() + ')')
            sub_widget_names.append(ttp.param.name)
            sub_widgets.append(ttp.tk_wid)
            sub_widget_indices.append(i)
    # Second pass: link sub_widgets to their masters
    for ttp in tk_tao_list:
        sub = ttp.param.sub_param
        if sub in sub_widget_names:
            ttp.sub_wid = sub_widgets[sub_widget_names.index(sub)]
    # Third pass: remove ignored widgets
    output = []
    sub_widget_indices.append(len(tk_tao_list)) #fixes the loop below
    for i in range(len(sub_widget_indices)):
        ix = sub_widget_indices[i]
        if i==0:
            output += tk_tao_list[:ix]
        else:
            last_ix = sub_widget_indices[i-1]
            output += tk_tao_list[last_ix+1:ix]
    return output
