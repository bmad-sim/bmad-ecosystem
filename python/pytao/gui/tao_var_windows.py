'''
Provides windows for viewing and editing variable arrays in tao
'''
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog

from pytao.util.parameters import str_to_tao_param, tao_parameter_dict
from .tao_widget import *
from .tao_set import *
from .tao_base_windows import *
from .tao_lat_windows import tao_ele_browser
#-----------------------------------------------------
# Variable Window

class tao_var_general_window(tao_list_window):

    def __init__(self, root, pipe, *args, **kwargs):
        self.tao_id = 'var'
        tao_list_window.__init__(self, root, "v1 Variables", *args, **kwargs)
        self.pipe = pipe
        for i in [0,4,5]:
            self.list_frame.grid_columnconfigure(i, pad=10)
        self.refresh()

    def refresh(self):
        for child in self.list_frame.winfo_children():
            child.destroy()
        v1_var_list = self.pipe.cmd_in("python var_general")
        v1_var_list = v1_var_list.splitlines()
        # close if there are no variables
        if len(v1_var_list) == 0:
            messagebox.showwarning('Warning', 'No variables defined', parent=self)
            self.destroy()
            return
        for i in range(len(v1_var_list)):
            v1_var_list[i] = v1_var_list[i].split(';')

        tk.Label(self.list_frame, text="Variable").grid(
                row=0, column=0, columnspan=4, sticky='W')
        tk.Label(self.list_frame, text="Indices").grid(row=0, column=4)
        tk.Label(self.list_frame, text="Using").grid(row=0, column=5)

        i=1
        for item in v1_var_list:
            tk.Label(self.list_frame, text=item[0]).grid(row=i, column=0, sticky='W')
            tk.Button(self.list_frame, text="View...",
                    command=self.open_v1_callback(item[0])).grid(row=i, column=1)
            tk.Button(self.list_frame, text="Edit...",
                    command=self.edit_v1_callback(item[0])).grid(row=i, column=2)
            tk.Button(self.list_frame, text="Delete...",
                    command=self.delete_v1_callback(item[0])).grid(row=i, column=3)
            #tk.Button(self.list_frame, text="Write...",
            #        command=self.write_v1).grid(row=i, column=3)
            tk.Label(self.list_frame,text=item[2] +':'+ item[3]).grid(row=i,column=4)
            tk.Label(self.list_frame, text=item[1]).grid(row=i, column=5)
            i = i+1

    def open_v1_callback(self, v1_var_name):
        return lambda : self.open_v1(v1_var_name)

    def open_v1(self, v1_var_name):
        win = tao_v1_var_window(self.root, self.pipe, v1_var_name)

    def edit_v1_callback(self, v1_var_name):
        return lambda : self.edit_v1(v1_var_name)

    def edit_v1(self, v1_var_name):
        win = tao_new_var_window(self.root, self.pipe, default=v1_var_name)

    def delete_v1_callback(self, v1_var_name):
        return lambda : self.delete_v1(v1_var_name)

    def delete_v1(self, v1_var_name):
        ans = messagebox.askokcancel(title="Delete " + v1_var_name + "?",
                message="Are you sure you want to delete " + v1_var_name
                + "?")
        if ans:
            self.pipe.cmd_in("python var_v1_destroy " + v1_var_name)
            for win in self.root.refresh_windows['var']:
                win.refresh()
            for win in self.root.refresh_windows['plot']:
                win.refresh()


#-----------------------------------------------------
# v1_var window

class tao_v1_var_window(lw_table_window):

    def __init__(self, root, pipe, v1_var_name, *args, **kwargs):
        self.tao_id = 'var'
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

    def bulk_set(self, fill_choices, parent):
        '''
        Overload of lw_table_window.bulk_set that calls refresh data and plot windows
        '''
        lw_table_window.bulk_set(self, fill_choices, parent)
        # Refresh data-related windows
        for win in self.root.refresh_windows['var']:
            win.refresh()
        for win in self.root.refresh_windows['plot']:
            win.refresh()

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

    def detail_set_callback(self, tao_list, set_str):
        '''
        Overload of lw_table_window.detail_set_callback that refreshes
        data-related windows
        '''
        lw_table_window.detail_set_callback(self, tao_list, set_str)
        # Refresh data-related windows
        for win in self.root.refresh_windows['var']:
            win.refresh()
        for win in self.root.refresh_windows['plot']:
            win.refresh()

class tao_new_var_window(Tao_Toplevel):
    '''
    Provides a window for creating new v1_variable arrays
    '''
    def __init__(self, root, pipe, default=None, *args, **kwargs):
        self.root = root
        Tao_Toplevel.__init__(self, root, *args, **kwargs)
        self.pipe = pipe
        self.title('New Variables')
        self.v1_frame_list = []

        # Possibly create self.notebook
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(side='top', fill='both', expand=1)

        # Set up a v1 tab
        self.v1_frame_list.append(new_v1_frame(self))
        self.notebook.insert('end', self.v1_frame_list[0])
        self.notebook.tab(0, text=self.v1_frame_list[0].BLANK_TITLE)
        self.v1_index = 0 #marks current tab index

        # New tab button
        self.new_tab_frame = tk.Frame(self.notebook)
        self.notebook.insert('end', self.new_tab_frame)
        self.notebook.tab(1, text='+')
        self.notebook.bind('<<NotebookTabChanged>>', self.tab_handler)

        # Create button
        self.create_b = tk.Button(self, text="Create!", command=self.create_variables)
        self.create_b.pack(side='right')

        # Clone default if it exists
        if default != None:
            self.v1_frame_list[0].clone(default)

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
            if not v1_frame.name_handler(strict=True):
                name_m = "Please check v1_array names."
                if name_m not in messages:
                    messages.append(name_m)
            # Check min indices
            if not v1_frame.ix_min_handler(strict=True):
                messages.append("Please check the start index for " + v1_frame.name)
            # Check lengths
            if not v1_frame.length_handler(strict=True):
                messages.append("Please check the end index for " + v1_frame.name)
            # Check low/high limits
            if not v1_frame.low_high_handler(strict=True):
                messages.append("Please check the low and high limits for " + v1_frame.name)
            # Check for semicolons in any fields
            semi_message = "Semicolons not allowed in any input field"
            caret_message = "Carets not allowed in any input field"
            space_message = "Spaces not allowed in any input field"
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
                if str(ttp.tk_var.get()).find(' ') != -1:
                    messages.append(space_message)
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
                    if str(ttp.tk_var.get()).find(' ') != -1:
                        messages.append(space_message)
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
                cmd_str += '[' + str(i) + ']^^'
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
                    if j != 0:
                        v1_ix = v1_params.index(p)
                    else:
                        v1_ix = 0
                    # Highest priority: individual var settings
                    if p in var_dict.keys():
                        if p == 'universes' and var_dict[p]=="": #cannot be empty
                            u = v1_frame.v1_array_wids[v1_ix].tk_var.get()
                            if u == "":
                                u = '*'
                            cmd_str += u + '^^'
                        else:
                            if v1_frame.v1_array_wids[v1_ix].param.type == 'LOGIC':
                                cmd_str += ('T^^' if var_dict[p] else 'F^^')
                            else:
                                cmd_str += str(var_dict[p]) + '^^'
                    else: #Use defaults if var parameter not set
                        if j == 0:
                            cmd_str += '^^'
                            continue
                        if v1_frame.v1_array_wids[v1_ix].param.type == 'LOGIC':
                            cmd_str += ('T^^' if v1_frame.v1_array_wids[v1_ix].tk_var.get() else 'F^^')
                        elif p == 'universes': #cannot be empty
                            u = v1_frame.v1_array_wids[v1_ix].tk_var.get()
                            if u == "":
                                u = '*'
                            cmd_str += u + '^^'
                        else:
                            cmd_str += v1_frame.v1_array_wids[v1_ix].tk_var.get() + '^^'
                cmd_str = cmd_str[:-2] # remove ending ^^
                self.pipe.cmd_in(cmd_str)
        # Close the window
        self.destroy()
        # Refresh variable-related windows
        for win in self.root.refresh_windows['var']:
            win.refresh()
        for win in self.root.refresh_windows['plot']:
            win.refresh()

    def tab_handler(self, event=None):
        '''
        Handles new tab creation and updates self.v1_index as necessary
        '''
        # Check if the new tab frame has been selected
        if self.notebook.select() == self.new_tab_frame._w:
            # Add new tab
            self.v1_frame_list.append(new_v1_frame(self))
            self.v1_index = len(self.v1_frame_list)-1
            self.notebook.insert(self.v1_index, self.v1_frame_list[-1])
            self.notebook.tab(self.v1_index, text=self.v1_frame_list[-1].BLANK_TITLE)
            self.notebook.select(self.v1_index)
        else:
            # Update self.v1_index
            for i in range(len(self.v1_frame_list)):
                frame = self.v1_frame_list[i]
                if self.notebook.select() == frame._w:
                    self.v1_index = i
                    # Unblock frame's handlers
                    frame.handler_block = False
                    break

    def refresh(self):
        '''
        Only here in case something tries to refresh this window
        '''
        pass


class new_v1_frame(tk.Frame):
    '''
    Provides a frame for inputting properties of a v1_variable array.
    '''
    def __init__(self, parent):
        tk.Frame.__init__(self, parent.notebook)
        self.parent = parent
        self.pipe = self.parent.pipe
        self.BLANK_TITLE = "[UNTITLED]"
        self.name = self.BLANK_TITLE #Default
        self.clone_of = "" # Used to track clones of existing v1_arrays
        self.handler_block = False

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
                "Default weight:", "Default step:", "Default merit type:",
                "Default low_lim:", "Default high_lim:", "Default good_user:",
                "Default key_bound:", "Default key_delta:", "Base index:", "Array length:"]
        # Grid widgets and labels:
        for i in range(len(self.v1_array_wids)):
            tk.Label(self, text=self.v1_array_labels[i]).grid(row=i+3, column=0, sticky='W')
            self.v1_array_wids[i].tk_wid.grid(row=i+3, column=1, sticky='EW')
        i = i+3

        # Warning labels
        # (defined here to be gridded/ungridded as necessary)
        self.name_warning_empty = tk.Label(self, text="Must not be empty")
        self.name_warning_invalid = tk.Label(self, text="v1 name already in use")
        self.ix_min_warning = tk.Label(self, text="Must be a non-negative integer")
        self.length_warning = tk.Label(self, text="Must be a positive integer")
        self.low_warning = tk.Label(self, text="Must be a real number")
        self.high_warning = tk.Label(self, text="Must be a real number")
        self.low_high_warning = tk.Label(self, text="Min value must be smaller than max value")

        # Grid settings for the above warnings
        self.name_warning_gs = {"row":3, "column":2, "sticky":'EW'}
        self.ix_min_warning_gs = {"row":14, "column":2, "sticky":'EW'}
        self.length_warning_gs = {"row":15, "column":2, "sticky":'EW'}
        self.low_warning_gs = {"row":9, "column":3, "sticky":'EW'}
        self.high_warning_gs = {"row":10, "column":3, "sticky":'EW'}
        self.low_high_warning_gs = {"row":9, "column":3, "rowspan":2, "sticky":'EW'}

        # Responses to edits
        self.v1_array_wids[0].tk_wid.bind('<FocusOut>', self.name_handler)
        self.v1_array_wids[11].tk_wid.bind('<FocusOut>', self.ix_min_handler)
        self.v1_array_wids[12].tk_wid.bind('<FocusOut>', self.length_handler)
        self.v1_array_wids[6].tk_wid.bind('<FocusOut>', self.low_high_handler)
        self.v1_array_wids[7].tk_wid.bind('<FocusOut>', self.low_high_handler)

        # "Fill" to individual vars
        def var_fill_callback(ix):
            '''Helper function for buttons below'''
            return lambda : self.var_fill(ix)
        for ix in range(1, 11):
            tk.Button(self, text="Fill to vars", command=var_fill_callback(ix)).grid(row=ix+3, column=2)

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

        # Duplicate button
        tk.Button(self, text="Duplicate this v1_array", command=self.duplicate).grid(
                row=1, column=0, columnspan=3, sticky='EW')

        # Clone existing v1_array
        tk.Label(self, text="Clone existing v1:").grid(row=2, column=0, sticky='E')
        v1_list = self.pipe.cmd_in('python var_general').splitlines()
        for x in range(len(v1_list)):
            v1_list[x] = v1_list[x].split(';')[0]
        self.v1_clone = tk.StringVar()
        self.v1_clone.set(v1_list[0])
        tk.OptionMenu(self, self.v1_clone, *v1_list).grid(row=2, column=1, sticky='EW')
        tk.Button(self, text="Clone", command= lambda: self.clone(self.v1_clone.get())).grid(row=2, column=2, sticky='W')

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

    def duplicate(self, event=None):
        '''
        Adds a new v1_frame to self.parent.v1_frame_list that is a copy of
        this frame, and changes focus to that frame
        '''
        # Don't run any handlers for this d1_frame
        self.handler_block = True
        self.parent.v1_frame_list.append(new_v1_frame(self.parent))
        ix = len(self.parent.v1_frame_list) - 1
        #self.d2_array.d1_index = ix
        # Copy properties into new frame
        self.parent.v1_frame_list[-1].name = self.name + '_copy'
        self.parent.v1_frame_list[-1].ix_min = self.ix_min
        self.parent.v1_frame_list[-1].ix_max = self.ix_max
        self.parent.v1_frame_list[-1].var_dict = self.var_dict
        for i in range(len(self.v1_array_wids)):
            if i == 0:
                self.parent.v1_frame_list[-1].v1_array_wids[i].tk_var.set(self.v1_array_wids[i].tk_var.get() + '_copy')
            else:
                self.parent.v1_frame_list[-1].v1_array_wids[i].tk_var.set(self.v1_array_wids[i].tk_var.get())
        # Run all input validation handlers
        self.parent.notebook.insert(ix, self.parent.v1_frame_list[-1])
        self.parent.notebook.select(ix)
        self.parent.tab_handler()
        self.update_idletasks()
        self.parent.v1_frame_list[-1].name_handler()
        self.parent.v1_frame_list[-1].ix_min_handler()
        self.parent.v1_frame_list[-1].length_handler()

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
        self.var_frame.grid(row=50, column=0, columnspan=3, sticky='EW')

    def var_fill(self, ix):
        '''
        "Fills" the value in self.v1_array_wids[ix] to all vars in self.var_dict
        Actually just erases the values in self.var_dict, which produces the same effect
        '''
        for i in self.var_dict.keys():
            param = self.v1_array_wids[ix].param.name
            if param in self.var_dict[i].keys():
                self.var_dict[i].pop(param)
            elif 'var^'+param in self.var_dict[i].keys():
                self.var_dict[i].pop('var^'+param)
        # Redraw datum frame to reflect changes
        try:
            int(self.var_ix.get())
            self.make_var_frame()
        except ValueError:
            pass

    def lat_browser(self):
        '''
        Opens a modified lattice table window to select elements en masse
        which should be either 'name', 'start_name', or 'ref_name'
        When the elements are chosen, ele_which will be set for each of this
        v1_array's variables sequentially
        '''
        # Make sure ix_min and length are not invalidly set
        if (self.ix_min_handler() == False) or (self.length_handler() == False):
            return
        # Autosize the array if self.length == -1
        if self.length == -1:
            autosize = True
        else:
            autosize = False
        win = tao_ele_browser(self.parent.root, self.pipe, self.name, self, 'var', autosize=autosize)

    def name_handler(self, event=None, strict=False):
        '''
        Checks whether the input name is valid, and
        sets the tab name to it if it is
        Returns True if the name is valid,
        False if the name is already in use by another v1,
        and None if the name is empty
        If strict is set True, a warning will be displayed
        for an empty name
        '''
        if self.handler_block:
            return
        name = self.v1_array_wids[0].tk_var.get().strip()
        if name == "":
            self.name_warning_invalid.grid_forget()
            self.name_warning_empty.grid_forget()
            # Warning in strict mode
            if strict:
                self.name_warning_empty.grid(**self.name_warning_gs)
            self.name = self.BLANK_TITLE
            self.parent.notebook.tab(self.parent.v1_index, text=self.BLANK_TITLE)
            return None
        # Name is nonempty
        # -> check if it's already in use
        i = 0
        for v1 in self.parent.v1_frame_list:
            if (v1.name == name) & (self != self.parent.v1_frame_list[i]):
                self.name_warning_empty.grid_forget()
                self.name_warning_invalid.grid(**self.name_warning_gs)
                return False
            i = i+1

        # Name is nonempty and valid
        self.parent.notebook.tab(self.parent.v1_index, text=name)
        self.name_warning_empty.grid_forget()
        self.name_warning_invalid.grid_forget()
        self.name = name
        if self.name == self.clone_of:
            return True # prevent repeatedly calling this method
        # Check if v1_array already exists in tao
        v1_list = self.pipe.cmd_in('python var_general').splitlines()
        for line in v1_list:
            if line.split(';')[0] == name:
                ans_var = tk.StringVar()
                tao_message_box(self.parent.root, self.parent, ans_var,
                        title=name+ " already exists",
                        message=name + " is already defined.    Would you like to edit its existing properties or overwrite them?",
                        choices=['Edit', 'Overwrite'])
                if ans_var.get() == 'Edit':
                    self.clone(self.name)
        return True

    def ix_min_handler(self, event=None, strict=False):
        '''
        Checks if the starting index is valid,
        and updates self.ix_min if it is.
        Returns True if a good index is set,
        False if a bad index is set
        (non-integer or greater than ix_max),
        and None if the ix_min field is empty
        If strict is set true, prints a warning message
        for an empty ix_min
        '''
        if self.handler_block:
            return
        # Check if min index is empty
        ix_min = self.v1_array_wids[11].tk_var.get().strip()
        if ix_min == "":
            self.ix_min = -1
            self.ix_min_warning.grid_forget()
            if strict:
                self.ix_min_warning.grid(**self.ix_min_warning_gs)
            #print("ix_min empty")
            return None
        # Check if min index is a non-negative integer
        try:
            ix_min = int(ix_min)
        except ValueError:
            ix_min = -1
        if ix_min < 0:
            self.ix_min_warning.grid(**self.ix_min_warning_gs)
            self.ix_min = -1
            #print("ix_min invalid")
            return False
        # ix_min is a good value -> update self.ix_min
        #print("valid ix_min")
        self.ix_min = ix_min
        self.ix_min_warning.grid_forget()
        # Update the data index range if self.length is set
        if self.length != -1:
            #print("updating ix_max and data_chooser")
            self.ix_max = self.ix_min + self.length - 1
            self.var_chooser.configure(values=list(range(self.ix_min, self.ix_max+1)))
        return True

    def length_handler(self, event=None, strict=False):
        '''
        Checks that the v1 array length has been set
        to a valid value, and updates self.length if it has.
        Returns True if the length has been set to a valid
        (positive integer) value, False if set to an invalid
        value, and None if left blank.
        '''
        if self.handler_block:
            #print("handler blocked")
            return
        # Check if length is empty
        length = self.v1_array_wids[12].tk_var.get().strip()
        if length == "":
            self.length = -1
            self.length_warning.grid_forget()
            if strict:
                self.length_warning.grid(**self.length_warning_gs)
            #print("length empty")
            return None
        # Check if length is a positive integer
        try:
            length = int(length)
        except ValueError:
            length = -1
        if length < 1:
            self.length_warning.grid(**self.length_warning_gs)
            self.length = -1
            #print("length invalid")
            return False
        # length is a good value -> update self.length
        self.length = length
        self.length_warning.grid_forget()
        #print("length valid")
        #print("self.ix_min:" + str(self.ix_min))
        # Update the data index range if self.ix_min is set
        if self.ix_min != -1:
            #print("updating self.ix_max and data_chooser")
            self.ix_max = self.ix_min + self.length - 1
            self.var_chooser.configure(values=list(range(self.ix_min, self.ix_max+1)))
        return True

    def low_high_handler(self, event=None, strict=False):
        '''
        Verifies that the min and max values are set to good values
        Returns True if they are,
        False if at least one is set invalidly,
        or None otherwise (at least one is blank)
        Prints an error for an empty low or high limit
        if strict is set True
        '''
        if self.handler_block:
            return
        self.low_warning.grid_forget()
        self.high_warning.grid_forget()
        self.low_high_warning.grid_forget()
        # Check if both limits are empty
        low_lim = self.v1_array_wids[6].tk_var.get()
        high_lim = self.v1_array_wids[7].tk_var.get()
        if (low_lim == "") and (high_lim == ""): #both are empty
            if strict:
                self.low_warning.grid(**self.low_warning_gs)
                self.high_warning.grid(**self.high_warning_gs)
            return None
        elif (low_lim == ""): #low is empty
            try:
                high_lim = float(high_lim)
            except ValueError:
                high_lim = None
            if strict:
                self.low_warning.grid(**self.low_warning_gs)
            if high_lim == None:
                self.high_warning.grid(**self.high_warning_gs)
                return False
            return None
        elif (high_lim == ""): #high is empty
            try:
                low_lim = float(low_lim)
            except ValueError:
                low_lim = None
            if strict:
                self.high_warning.grid(**self.high_warning_gs)
            if low_lim == None:
                self.low_warning.grid(**self.low_warning_gs)
                return False
            return None
        else: #neither are empty
            try:
                low_lim = float(low_lim)
            except ValueError:
                low_lim = None
                self.low_warning.grid(**self.low_warning_gs)
            try:
                high_lim = float(high_lim)
            except ValueError:
                high_lim = None
                self.high_warning.grid(**self.high_warning_gs)
            if None in [low_lim, high_lim]:
                return False
            # Check that low_lim < high_lim
            if not (low_lim < high_lim):
                self.low_high_warning.grid(**self.low_high_warning_gs)
                return False
            return True

    def clone(self, v1_array):
        '''
        Clones the existing v1_array into this v1_frame
        '''
        # Read existing var info from tao
        var_general = self.pipe.cmd_in('python var_general').splitlines()
        for line in var_general:
            if line.find(v1_array) == 0:
                break
        line = line.split(';')
        # set name
        self.v1_array_wids[0].tk_var.set(v1_array)
        self.clone_of = v1_array
        self.name_handler()
        # set min/max indices
        self.v1_array_wids[11].tk_var.set(line[2])
        self.v1_array_wids[12].tk_var.set(line[3])
        self.ix_min_handler()
        self.length_handler()
        # copy inidividual variable info
        var_params = ['ele_name', 'attribute', 'weight', 'step', 'var^merit_type', 'low_lim', 'high_lim', 'good_user', 'key_bound', 'key_delta']
        self.var_dict = {}
        for i in range(self.ix_min, self.ix_max+1):
            self.var_dict[i] = {}
            var_info = self.pipe.cmd_in('python var ' + v1_array + '[' + str(i) + ']')
            var_info = var_info.splitlines()
            for p in var_params:
                found=False
                for line in var_info:
                    if line.find(p) == 0:
                        found=True
                        break
                if found:
                    self.var_dict[i][p] = line.split(';')[3].strip()
            # Special case: attribute (vs attrib_name)
            found=False
            for line in var_info:
                if line.find('attrib_name') == 0:
                    found=True
                    break
            if found:
                self.var_dict[i]['attribute'] = line.split(';')[3].strip()
            # Special case: universes
            unis = []
            for line in self.pipe.cmd_in('python var ' + v1_array + '[' + str(i) + '] slaves').splitlines():
                u = int(line.split(';')[0])
                if u not in unis:
                    unis.append(u)
            unis.sort()
            universes = ""
            for u in unis:
                universes += str(u) + ','
            universes = universes[:-1] # remove extra comma
            self.var_dict[i]['universes'] = universes
        # Set v1 defaults to var_dict[self.ix_min] values
        for i in range(1, len(var_params)):
            self.v1_array_wids[i+1].tk_var.set(self.var_dict[self.ix_min][var_params[i]])

