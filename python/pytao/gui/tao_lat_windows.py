'''
Provides windows for viewing the lattice and editting element attributes
'''
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

from pytao.util.parameters import str_to_tao_param, tao_parameter_dict
from pytao.util.lattice_element import *
from .tao_widget import *
from .tao_set import *
from .tao_base_windows import *
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
        self.tao_id = 'ele'
        tao_list_window.__init__(self, root, "Lattice Elements", use_upper=True,
                min_width=600, *args, **kwargs)
        # Make sure the user can't close the window without being asked
        # to save changes to the current element
        self.protocol("WM_DELETE_WINDOW", self.close_cmd)
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

    def refresh(self, event=None, ask=True, *args):
        '''
        This is where most of the element information is actually created
        '''
        # Ask to save changes
        if ask and self.check_for_changes():
            x = messagebox.askyesnocancel(title="Unsaved Changes",
                    message="Apply changes before switching elements?", parent=self)
            if x:
                self.ele_set()
            if x == None:
                return #don't refresh if "Cancel" is picked
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
        tk.Label(self.head_frame, text="Key").grid(row=1, column=0, sticky='W')
        self.element.params["key"].tk_wid.grid(row=1, column=1, sticky='EW')
        tk.Label(self.head_frame, text="s").grid(row=2, column=0, sticky='W')
        self.element.params["s"].tk_wid.grid(row=2, column=1, sticky='EW')
        tk.Label(self.head_frame, text="s_start").grid(row=3, column=0, sticky='W')
        self.element.params["s_start"].tk_wid.grid(row=3, column=1, sticky='EW')
        tk.Label(self.head_frame, text="Ref time").grid(row=4, column=0, sticky='W')
        self.element.params["ref_time"].tk_wid.grid(row=4, column=1, sticky='EW')

        #Variable parameters
        tk.Label(self.head_frame, text="Type").grid(row=1, column=2, sticky='W')
        self.element.params["type"].tk_wid.grid(row=1, column=3, sticky='EW')
        tk.Label(self.head_frame, text="Alias").grid(row=2, column=2, sticky='W')
        self.element.params["alias"].tk_wid.grid(row=2, column=3, sticky='EW')
        tk.Label(self.head_frame, text="Description").grid(
                row=3, column=2, sticky='W')
        self.element.params["descrip"].tk_wid.grid(row=3, column=3, sticky='EW')
        tk.Label(self.head_frame, text="is_on").grid(row=4, column=2, sticky='W')
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
                    key = "elec_multipoles"      # elec_multipoles in tao_python_cmd.f90
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
                    #        + self.element.id + ' err')
                    tao_list = self.pipe.cmd_in("python ele:mat6 "
                            + self.element.id + ' mat6')
                    #tao_list += '\n' + self.pipe.cmd_in("python ele:mat6 "
                    #        + self.element.id + ' vec0')
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
                    sym_err.tk_label.grid(row=0, column=0, sticky='W')
                    sym_err.tk_wid.grid(row=0, column=1, sticky='W')
                    # add vec0
                    vec0 = self.pipe.cmd_in('python ele:mat6 '
                            + self.element.id + ' vec0')
                    vec0 = tk_tao_parameter(str_to_tao_param(vec0),
                            self.p_frames[i], self.pipe)
                    separator = tk.Frame(self.p_frames[i], height=10, bd=1).grid(
                            row=8, column=1, columnspan=6, sticky='EW')
                    vec0.tk_label.grid(row=9, column=0, sticky='EW')
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

    def close_cmd(self, event=None):
        '''Ask the user to save any unsaved changes before closing'''
        if self.check_for_changes():
            name = self.element.params["name"].param.value
            x = messagebox.askyesnocancel(title="Unsaved Changes",
                    message="Save changes to " + name + " before closing?", parent=self)
            if x:
                self.ele_set(refresh_self=False)
            if x == None:
                return #don't close if "Cancel" is picked
        # Actually close the window
        self.destroy()

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

    def ele_set(self, event=None, refresh_self=True):
        '''
        Runs set commands for all the parameter frames and refreshes the window
        This window will only be refreshed if refresh_self=True
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
        # Refresh ele-dependent windows (including self)
        for win in self.root.refresh_windows['ele']:
            if (self!=win):
                win.refresh()
            elif refresh_self:
                win.refresh(ask=False)
        for win in self.root.refresh_windows['data']:
            win.refresh()
        for win in self.root.refresh_windows['var']:
            win.refresh()
        for win in self.root.refresh_windows['lat']:
            win.refresh()
        if self.pipe.cmd_in('python lat_calc_done').find('T') != -1:
            for win in self.root.refresh_windows['plot']:
                win.refresh()

    def check_for_changes(self):
        '''Returns True if the current element has unsaved changes'''
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
        return something_changed




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
                name = name[:name.find('n')] + str(line[0]) + name[name.find('n')+1:]
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
class tao_lattice_window(Tao_Toplevel):
    '''
    Shows lattice elements in a read-only table view
    with an interface to select which rows/columns
    are displayed
    '''
    # TODO: replace blank
    def __init__(self, root, pipe, switches="", *args, **kwargs):
        self.root = root
        self.tao_id = 'lat'
        Tao_Toplevel.__init__(self, root, *args, **kwargs)
        self.title("Lattice")
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
        win = Tao_Toplevel.__new__(Tao_Toplevel)
        win.root = self.root
        Tao_Toplevel.__init__(win, self)
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
        # Apply button
        tk.Button(win, text="Apply Filters", command=self.refresh).grid(
                row=5, column=0, columnspan=2)

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
            NO_COLS=False # Hack to allow no attribute columns
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
                if self.col_atts.get().strip() == "":
                    NO_COLS=True
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
                if self.ele_list.get().strip() != "":
                    switches += self.ele_list.get().strip()
                else:
                    switches += '-all '

            if NO_COLS:
                self.switches = switches + " -att"
            else:
                self.switches = switches
            self.advanced_var.set(switches)

    def refresh(self, event=None):
        '''
        Fetches the lattice table with the new table
        parameters and updates the window
        Returns the number of elements found
        '''
        LENGTH_SCALE = 10 #pixels per character that should be allotted for each column
        # Clear the existing table
        for child in self.table_frame.winfo_children():
            child.destroy()

        # Get the new table data using the given switches
        self.get_switches()
        lattice = self.pipe.cmd_in("python show lattice -python " + self.switches)
        lattice = lattice.splitlines()
        #Remove error messages
        #while lattice[0][:6] in ['[ERROR', '[FATAL']:
        #    lattice = lattice[1:] #remove line with [ERROR or [FATAL
        #    while lattice[0].find('    ') == 0:
        #        lattice = lattice[1:] #remove error description lines
        #        if len(lattice) == 0:
        #            break
        #    if len(lattice) == 0:
        #        break
        if lattice[0][:6] in ['[ERROR', '[FATAL']:
            tk.Label(self.table_frame, text="NO LATTICE FOUND").pack()
            return 0
        if len(lattice) == 0:
            tk.Label(self.table_frame, text="NO LATTICE FOUND").pack()
            return 0
        if lattice[0].find("ELEMENT(S) NOT FOUND") != -1:
            tk.Label(self.table_frame, text="NO LATTICE FOUND").pack()
            return 0
        if lattice[0].find("ERROR READING RANGE SELECTION") != -1:
            tk.Label(self.table_frame, text="NO LATTICE FOUND").pack()
            return 0
        for i in range(len(lattice)):
            lattice[i] = lattice[i].split(';')
        #lattice[i][j] --> row i, column j
        dat_types = lattice[1] # marks columns as STR, INT, etc
        lattice = lattice[:1] + lattice[2:]
        # Trim off last column if empty
        if lattice[0][-1] == "":
            for i in range(len(lattice)):
                lattice[i] = lattice[i][:-1]
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
                if len(row[j])*LENGTH_SCALE > widths[j]:
                    widths[j] = len(row[j])*LENGTH_SCALE

        # Set column widths appropriately
        for j in range(len(lattice[0])):
            if len(lattice[0][j])*LENGTH_SCALE > widths[j]:
                widths[j] = len(lattice[0][j])*LENGTH_SCALE
            widths[0] = 100 #prevent giant index column
            self.tree.column(lattice[0][j], width=widths[j], minwidth=widths[j])
            # Text alignment
            if dat_types[j] == 'STR':
                self.tree.column(lattice[0][j], anchor='w')
            else:
                self.tree.column(lattice[0][j], anchor='e')

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
        self.minsize(100, 100)
        # Determine number of elements
        for row in self.tree.get_children():
            if len(self.tree.item(row)['values']) < 2:
                # Don't count the "Lord Elements:" row
                return len(lattice[1:]) - 1
        return len(lattice[1:])

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
    autosize: if set True, the parent d1/v1 array's length
    will be set to match the number of elements found,
    and the base index for the array will be set to 1 if
    it is not already set
    '''
    def __init__(self, root, pipe, name, parent, parent_type,
            which='name', uni="All", autosize=False, *args, **kwargs):
        if uni == "All":
            uni = "-1"
        self._full_init = False # used to delay self.refresh until init finishes
        tao_lattice_window.__init__(
                self, root, pipe, switches="-universe " + uni + " -all", *args, **kwargs)
        self.name = name
        self.title(self.name + ': ele_' + which + ' Browser')
        self.parent = parent
        self.parent_type = parent_type
        self.which = which
        self.autosize = autosize

        # Ungrid all widgets and re-grid the window
        # with only the widgets needed for the ele_browser
        for child in self.top_frame.grid_slaves():
            child.grid_forget()

        # First row: branch, row filters
        tk.Label(self.top_frame, text="Branch:").grid(row=0, column=0, sticky='W')
        self.branch_wids.branch_chooser.grid(row=0, column=1, sticky='EW')
        self.f_button.grid(row=0, column=2, sticky='W')
        self.f_button.configure(text="Row Filters...")

        # Only show name, key, s, l columns
        self.col_filter.set("List Attributes")
        self.col_atts.set("")

        # Second row: ele search
        tk.Label(self.top_frame, text="Search:").grid(row=1, column=0, sticky='W')
        self.ele_list_box.grid(row=1, column=1, columnspan=2, sticky='EW')
        self.ele_list_opt.set('Custom')
        self.ele_list_box.focus_set()
        self.search_button = tk.Button(self.top_frame, text="Search!", command=self.refresh)
        self.search_button.grid(row=0, column=3, rowspan=2, sticky='NSW')


        # Remove template widgets
        #self.temp_label_1.grid_forget()
        #self.template_file.tk_wid.grid_forget()
        #self.temp_label_2.grid_forget()
        #self.temp_chooser.grid_forget()
        #self.save_button.grid_forget()
        #self.temp_save.tk_wid.grid_forget()


        # Data index widgets
        self.data_frame = tk.Frame(self.top_frame)
        self.data_frame.grid(row=0, column=6, rowspan=2, sticky='NS')
        if not autosize:
            tk.Label(self.data_frame, text="Apply ele names to " + self.name + "[").pack(side="left")
            self.ix_min_var = tk.StringVar()
            tk.Entry(self.data_frame, width=5, textvariable=self.ix_min_var).pack(side="left")
            tk.Label(self.data_frame, text="] through " + self.name + "[").pack(side="left")
            self.ix_max_var = tk.StringVar()
            tk.Entry(self.data_frame, width=5, textvariable=self.ix_max_var).pack(side="left")
            tk.Label(self.data_frame, text="]").pack(side="left")
            self.ix_min_var.set(self.parent.ix_min)
            self.ix_max_var.set(self.parent.ix_max)
        # placeholder:
        self.ele_count = tk.Label(self.data_frame, text="")

        self.apply_button = tk.Button(self.top_frame, text="Apply Element Names", command=self.apply_callback)
        self.apply_button.grid(row=0, column=7, rowspan=2, sticky='NS')
        # Refresh now
        self._full_init = True
        self.refresh()

    def apply_callback(self, event=None):
        '''
        Checks that good indices have been input in self.ix_min_var
        and self.ix_max_var, and then copies ele names into ele_which
        for the selected data/var indices
        '''
        if self.autosize:
            self._autosize_apply_callback()
        else:
            self._fixed_size_apply_callback()
        # Close this window
        self.destroy()

    def _autosize_apply_callback(self):
        '''
        Apply callback implementation when operating with autosize == True
        '''
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
        # Set ix_min and ix_max
        if self.parent.ix_min == -1: # default to 1
            ix_min = 1
            if self.parent_type == 'data':
                #print("setting ix_min to default 1")
                self.parent.d1_array_wids[7].tk_var.set(str(ix_min))
            elif self.parent_type == 'var':
                self.parent.v1_array_wids[11].tk_var.set(str(ix_min))
            self.parent.ix_min_handler()
        else:
            ix_min = self.parent.ix_min
        # Set parent array length and ix_max
        if self.parent_type == 'data':
            self.parent.d1_array_wids[8].tk_var.set(str(len(names)))
        elif self.parent_type == 'var':
            self.parent.v1_array_wids[12].tk_var.set(str(len(names)))
        self.parent.length_handler()
        ix_max = self.parent.ix_max
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


    def _fixed_size_apply_callback(self):
        '''
        Apply callback implementation when operating with autosize == False
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


    def refresh(self, event=None):
        '''
        Overload of tao_lattice_window.refresh() that prints the number
        of elements found
        '''
        if not self._full_init:
            return
        # Call original version
        num_eles = tao_lattice_window.refresh(self, event)
        if num_eles != 0:
            # Disable clicking on table
            self.tree.configure(selectmode="none")
            self.tree.unbind('<Double-Button-1>')
        # Resize columns
        #self.tree.column("Index", width=50, minwidth=50)
        #self.minsize(300,100)
        #self.maxsize(1000,1000)
        # Count elements found
        #try:
        #    names = list(self.tree.get_children())
        #    for i in range(len(names)):
        #        if len(self.tree.item(names[i])['values']) < 2:
        #            # Remove "Lord Elements:" row
        #            names.pop(i)
        #            break
        #except tk._tkinter.TclError:
        #    names = [] # no lattice found
        # Display ele count
        self.ele_count.destroy()
        count_text = str(num_eles) + " elements found"
        if not self.autosize:
            count_text = "(" + count_text + ")"
        self.ele_count = tk.Label(self.data_frame, text=count_text)
        self.ele_count.pack(side='left')
