'''
Contains the base classes used by many (if not all) of the windows
in the GUI for tao.
'''
import sys
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import copy
from .tao_widget import *
from .tao_set import *

#-----------------------------------------------------
# Tao_Toplevel window
class Tao_Toplevel(tk.Toplevel):
    '''
    Skeleton class that handles adding/removing windows
    from the root refresh windows dictionary
    Subclasses should set self.tao_id before running this
    __init__ method.
    Subclasses MUST set self.root before running this
    __init__ method.
    '''
    def __init__(self, parent, *args, **kwargs):
        tk.Toplevel.__init__(self, parent, class_='Tao', *args, **kwargs)
        if sys.platform == "linux":
            self.iconbitmap(self.root.icon)
        # Handle root window list placement
        if 'tao_id' not in self.__dict__:
            self.tao_id = None
        for tao_id in self.root.tao_id_list:
            if self.tao_id == tao_id:
                self.root.refresh_windows[tao_id].append(self)
                break

    def destroy(self):
        '''
        Overloaded here to remove self from root.refresh_windows
        '''
        # Try to remove self from root.refresh_windows
        for tao_id in self.root.tao_id_list:
            if self.tao_id == tao_id:
                try:
                    self.root.refresh_windows[tao_id].pop(self.root.refresh_windows[tao_id].index(self))
                except ValueError:
                    pass
        # Call Toplevel.destroy()
        tk.Toplevel.destroy(self)


class Tao_Popup(Tao_Toplevel):
    '''
    Provides a more convenient way to initialize a toplevel
    window that doesn't need to be its own subclass, used for
    things like ele_shape windows, dialog boxes, etc

    parent: the parent window
    root: the tao root window
    tao_id: which set of refresh windows this window should belong to (optional)
    '''
    def __init__(self, parent, root, tao_id=None, *args, **kwargs):
        self.root = root
        self.tao_id = tao_id
        Tao_Toplevel.__init__(self, parent, *args, **kwargs)


#-----------------------------------------------------
# List window

class tao_list_window(Tao_Toplevel):
    '''
    Skeleton class with several commonly used features.
    Comes with a scrollable region and supports adding other widgets
    Subclasses set tao_id before running this __init__ to handle
    adding/removing themselves from root window lists
    '''
    def __init__(self, root, title, use_upper=False,
            min_width=0, parent=None, *args, **kwargs):
        self.root = root
        if parent != None:
            Tao_Toplevel.__init__(self, parent, *args, **kwargs)
        else:
            Tao_Toplevel.__init__(self, root, *args, **kwargs)
        self.title(title)

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
        try:
            if event.num == 4:
                self.canvas.yview_scroll(-1,"units")
            elif event.num == 5:
                self.canvas.yview_scroll(1,"units")
        except: #in case this gets called when it shouldn't
            pass

#-----------------------------------------------------
# Scrolling frame
class tao_scroll_frame(tk.Frame):
    '''
    Provides a frame with a vertical scrollbar attached
    parent: the parent of this frame
    min_width: the minimum width that this frame must have
    self.frame should be used for placing widgets into
    this frame
    '''
    #TODO: upgrade tao_list_window to use this for its scrolling
    def __init__(self, parent, min_width=0, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.min_width = min_width
        self._old_width=0
        self._canvas=tk.Canvas(self)
        self.frame=tk.Frame(self._canvas)
        scrollbar=tk.Scrollbar(self,orient="vertical",
                command=self._canvas.yview)
        self._canvas.configure(yscrollcommand=scrollbar.set)

        self.frame.bind("<Configure>", self.frame_handler)
        self._canvas.bind("<Configure>", self.canvas_handler)

        self.bind("<Enter>", self.bind_mouse)
        self.bind("<Leave>", self.unbind_mouse)

        self._canvas.grid(row=0, column=0, sticky='NSEW')
        scrollbar.grid(row=0, column=1, sticky='NS')
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self._canvas_window = self._canvas.create_window(
                (0,0),window=self.frame,anchor='nw')

    def frame_handler(self, event):
        self._canvas.configure(scrollregion=self._canvas.bbox("all"))
        frame_width = self.frame.winfo_reqwidth()
        #frame_width = event.width
        # try using winfo_reqwidth() TODO
        self.grid_columnconfigure(0, minsize=frame_width)
        self.update_idletasks()

    def canvas_handler(self, event):
        #canvas_width = self._canvas.winfo_reqwidth()#event.width
        canvas_width = event.width
        self._canvas.itemconfig(self._canvas_window, width=canvas_width)
        self.update_idletasks()

    def scrollhelper(self, event):
        self._canvas.configure(scrollregion=self._canvas.bbox("all"))
        new_width = self.frame.winfo_geometry()
        new_width = int(new_width[:new_width.find('x')])
        # Don't resize to a smaller size
        #if self._old_width > new_width:
        #    return
        # Don't resize to less than the min width
        if self.min_width > new_width:
            new_width = self.min_width
        self.grid_columnconfigure(0, minsize=new_width)
        #self._canvas.configure(width=new_width)
        # Set the width of the frame to match the canvas
        new_frame_width = self._canvas.winfo_width()
        if new_frame_width > 1:
            self._canvas.itemconfigure("self.frame", width=new_frame_width)
        self._old_width = new_width

    def bind_mouse(self, event):
        self.bind_all("<Button-4>", self.mouse_scroll)
        self.bind_all("<Button-5>", self.mouse_scroll)

    def unbind_mouse(self, event):
        self.unbind_all("<Button-4>")
        self.unbind_all("<Button-5>")

    def mouse_scroll(self, event):
        #self.canvas.yview_scroll(direction,"units")
        try:
            if event.num == 4:
                self._canvas.yview_scroll(-1,"units")
            elif event.num == 5:
                self._canvas.yview_scroll(1,"units")
        except: #in case this gets called when it shouldn't
            pass

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
    '''
    Generic window for displaying and editing a list of parameters
    '''
    def __init__(self, root, title, tao_list, pipe, plot="", parent=None, *args, **kwargs):
        tao_list_window.__init__(self, root, title, parent=parent, *args, **kwargs)
        self.button_frame = tk.Frame(self)
        self.button_frame.pack(side="top", fill="both", expand=0)
        self.tao_list = tao_list
        for k in range(len(self.tao_list)):
            self.tao_list[k] = tk_tao_parameter(
                    self.tao_list[k], self.list_frame, pipe, plot=plot)
        # Link/filter ignored parameters
        self.tao_list = tk_tao_linker(self.tao_list)
        for k in range(len(self.tao_list)):
            self.tao_list[k].tk_label.grid(row=k,column=0,sticky="W")
            self.tao_list[k].tk_wid.grid(row=k,column=1,sticky="EW")
            if self.tao_list[k].sub_wid != None:
                self.tao_list[k].sub_wid.grid(row=k, column=2, sticky='W')
        # Attempt to bind data_source and data_types together
        data_source_ix = None
        for k in range(len(self.tao_list)):
            if self.tao_list[k].param.name == "data_source":
                data_source_ix = k
                break
        if data_source_ix != None:
            # Define _data_source_handler
            def _dat_source_handler(*args):
                '''Updates _data_source for DAT_TYPE widgets in this window'''
                for w in self.tao_list:
                    if w.param.type == 'DAT_TYPE':
                        w._data_source = self.tao_list[data_source_ix].tk_var.get()
                        w._s_refresh(dat_source_swap=True)
            # Trace tk_var for the data_source widget
            self.tao_list[data_source_ix].tk_var.trace('w', _dat_source_handler)
            _dat_source_handler()

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
        win = tao_parameter_window(self.root, detail_title, self.param_list, self.pipe, parent=self)

        set_str = self.set_format.format(self.array_name, index)
        b = tk.Button(win.button_frame, text="Apply changes",
                command=lambda : self.detail_set_callback(win.tao_list,set_str))
        b.pack()

    def detail_set_callback(self, tao_list, set_str):
        tao_set(tao_list, set_str, self.pipe)
        self.refresh()

#-----------------------------------------------------
class lw_table_window(Tao_Toplevel):
    '''
    Light-weight version of table_window, intended to be less
    graphically/computationally intensive (for use with networked
    machines, for example).  Table is read only but individual
    entries can still be opened for editing with double click.
    '''
    def __init__(self, root, pipe, array_name, title_list, bulk_format,
            bulk_set_format, set_format, *args, **kwargs):
        self.root = root
        Tao_Toplevel.__init__(self, root, *args, **kwargs)
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

        self.sh_button = tk.Button(self.button_frame, text="Show all",
                command=self.toggle_show_all)
        self.sh_button.pack(side='left')
        self.hide_nonexist = True
        self.refresh()

    def toggle_show_all(self, event=None):
        self.hide_nonexist = not self.hide_nonexist
        newtext = ("Show all" if self.hide_nonexist
                else "Hide nonexistent rows")
        self.sh_button.configure(text=newtext)
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
        win = Tao_Toplevel.__new__(Tao_Toplevel)
        win.root = self.root
        Tao_Toplevel.__init__(win, self)
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
        win = tao_parameter_window(self.root, detail_title, self.param_list, self.pipe, parent=self)

        set_str = self.set_format.format(self.array_name, index)
        b = tk.Button(win.button_frame, text="Apply changes",
                command=lambda : self.detail_set_callback(win.tao_list,set_str))
        b.pack()

    def detail_set_callback(self, tao_list, set_str):
        tao_set(tao_list, set_str, self.pipe)
        self.refresh()

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
        tk.Toplevel.__init__(self, parent, class_='Tao', *args, **kwargs)
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
        tk.Toplevel.__init__(self, parent, class_='Tao', *args, **kwargs)
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


class tabbed_frame(tk.Frame):
    '''
    Provides a frame for showing multiple frames one at a time
    in a tabbed view.
    parent: The parent widget where this frame will be placed
    new_tab_func: returns a new frame to add to a new tab
    self will be passed as the first positional argument and should
    be the only required positional arguments
    new_tab_func should return a tk.Frame with an instance variable
    self.name, a string that should be used for the initial tab name
    new_tab_func will likely be the init method for a tk.Frame subclass
    pipe: the tao_interface in use
    '''
    def __init__(self, parent, new_tab_func):
        tk.Frame.__init__(self, parent)
        self.new_tab_func = new_tab_func
        # Create self.notebook
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(side='top', fill='both', expand=1)

        self.tab_index = 0 # Marks the currently selected tab
        self.tab_list = []

        # New tab button
        self.new_tab_frame = tk.Frame(self.notebook)
        self.notebook.insert('end', self.new_tab_frame)
        self.notebook.tab(len(self.tab_list), text='+')
        self.notebook.bind('<<NotebookTabChanged>>', self._tab_handler)

    def _tab_handler(self, event=None):
        '''
        Handles new tab creation and updates self.tab_index as necessary
        '''
        # Check if the new tab frame has been selected
        if self.notebook.select() == self.new_tab_frame._w:
            # Add new tab at the end
            self.add_tab(len(self.tab_list), self)
            self.notebook.select(len(self.tab_list)-1)
        else:
            # Update self.tab_index
            for i in range(len(self.tab_list)):
                frame = self.tab_list[i]
                if self.notebook.select() == frame._w:
                    self.tab_index = i
                    # Unblock frame's handlers
                    frame.handler_block = False
                    break

    def add_tab(self, ix=None, *args, **kwargs):
        '''
        Adds a new tab at the specified index
        If ix is not specified, the tab is added at the end of the list
        *args and **kwargs are passed to self.new_tab_func after self
        '''
        if ix == None:
            ix = len(self.tab_list)
        if ix < len(self.tab_list):
            self.tab_list = self.tab_list[:ix] + [self.new_tab_func(*args, **kwargs)] + self.tab_list[ix:]
        else:
            ix = len(self.tab_list)
            self.tab_list.append(self.new_tab_func(*args, **kwargs))
        self.notebook.insert(ix, self.tab_list[ix])
        self.notebook.tab(ix, text=self.tab_list[ix].name)

    def remove_tab(self, ix, destroy=False):
        '''
        Removes the tab at the specified position
        Does not destroy the removed frame unless destroy==True
        '''
        # Make sure ix is not too large
        if ix >= len(self.tab_list):
            ix = len(self.tab_list) - 1
        # Prefer the previous index tab
        if ix > 0:
            self.notebook.select(ix-1)
        else:
            self.notebook.select(0)
        # Remove the tab
        frame = self.tab_list.pop(ix)
        self.notebook.forget(ix)
        if destroy:
            frame.destroy()

    def update_name(self, ix=None):
        '''
        Updates the tab title to match the frame name for the
        tab in the specified position
        If ix==None, runs this method for all tabs
        '''
        if ix not in range(len(self.tab_list)):
            for i in range(len(self.tab_list)):
                self.update_name(i)
        else:
            self.notebook.tab(ix, text=self.tab_list[ix].name)

class ele_shape_frame(tk.Frame):
    '''
    Provides a frame for viewing and editing the shapes that are used to
    display lat_layout and floor_plan graphs in tao

    parent: the parent frame/Toplevel where this frame will be placed
    pipe: the tao_interface used to querry/set ele shapes
    which: either "lat_layout" or "floor_plan"
    '''
    def __init__(self, parent, root, pipe, which):
        self.parent = parent
        self.pipe = pipe
        self.root = root
        tk.Frame.__init__(self, parent)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.type = which

        # Format string for the shape_set command
        self.shape_set_format = "python shape_set " + self.type + "^^" \
                + "{shape_ix}^^{ele_id}^^{shape}^^{color}^^{size}" \
                + "^^{label}^^{draw}^^{multi}^^{line_width}"

        self.title_text = tk.StringVar()
        self.title = tk.Label(self, textvariable=self.title_text, font=('Sans', 16, 'bold'))
        self.title.grid(row=0, column=0, sticky='EW')

        # Shape table
        self.table_frame = tk.Frame(self)
        self.table_frame.grid(row=1, column=0, sticky='NSEW')
        self.title_list = ["Index", "Ele ID", "Shape", "Color",
                "Size", "Label", "Draw", "Multi", "Line Width"]
        self.keys = ["shape_ix", "ele_id", "shape", "color", "size", "label", "draw", "multi", "line_width"]
        self.shape_table = ttk.Treeview(self.table_frame, columns=self.title_list, show='headings')
        self.shape_list = []
        self.widths = []
        # Column titles
        for title in self.title_list:
            self.shape_table.heading(title, text=title)
            self.shape_table.column(title, stretch=True, anchor='center')
        # Scrollbars
        hbar = ttk.Scrollbar(
                self.table_frame, orient="horizontal", command=self.shape_table.xview)
        vbar = ttk.Scrollbar(
                self.table_frame, orient="vertical", command=self.shape_table.yview)
        self.shape_table.configure(xscrollcommand=hbar.set, yscrollcommand=vbar.set)

        vbar.pack(side="right", fill="y", expand=0)
        hbar.pack(side="bottom", fill='x', expand=0)
        self.shape_table.pack(side="left", fill="both", expand=1)

        # Button Frame
        self.button_frame = tk.Frame(self)
        self.button_frame.grid(row=2, column=0, sticky='NSEW')
        self.move_up_b = tk.Button(self.button_frame, text="Move up",
                command = self.move_up)
        self.move_down_b = tk.Button(self.button_frame, text="Move down",
                command = self.move_down)
        self.dup_b = tk.Button(self.button_frame, text="Duplicate selected row",
                command = self.duplicate_item)
        self.del_b = tk.Button(self.button_frame, text="Delete selected row",
                command = self.delete_item)
        self.new_top_b = tk.Button(self.button_frame, text="New at top",
                command = self.new_at_top)
        self.new_bottom_b = tk.Button(self.button_frame, text="New at bottom",
                command = self.new_at_bottom)
        self.move_up_b.grid(row=0, column=0, sticky='EW')
        self.move_down_b.grid(row=0, column=1, sticky='EW')
        self.dup_b.grid(row=0, column=2, sticky='EW')
        self.del_b.grid(row=1, column=2, sticky='EW')
        self.new_top_b.grid(row=1, column=0, sticky='EW')
        self.new_bottom_b.grid(row=1, column=1, sticky='EW')
        for i in range(3):
            self.button_frame.grid_columnconfigure(i, weight=1)

        self.refresh()

    def refresh(self):
        '''
        Sets the table contents equal to the output of
        python shape_list lat_layout/floor_plan as requested
        '''
        if self.type not in ["lat_layout", "floor_plan"]:
            return
        # Clear existing rows
        for item in self.shape_table.get_children():
            self.shape_table.delete(item)
        self.shape_list = []

        # Fill rows
        ele_shapes = self.pipe.cmd_in('python shape_list ' + self.type)
        ele_shapes = ele_shapes.splitlines()
        for row in ele_shapes:
            row = row.split(';')
            self.shape_list.append(row)
            self.shape_table.insert("", "end", values=row)

        # Set column widths appropriately
        self.width_adjust()

        # Double click to edit shape
        self.shape_table.bind('<Double-Button-1>', self.edit_shape)

        # Set title
        self.title_text.set("Lat Layout Shapes" if self.type=="lat_layout" else "Floor Plan Shapes")

    def width_adjust(self, event=None):
        '''
        Adjusts column widths and total table width to show all contents
        '''
        # Determine column widths
        widths = [0]*len(self.title_list)
        for i in range(len(self.shape_list)):
            for j in range(len(self.shape_list[i])):
                widths[j] = max([widths[j], 15*len(str(self.shape_list[i][j]))])
        # Update column widths and table width
        for j in range(len(widths)):
            widths[j] = max([len(self.title_list[j])*12, widths[j]])
            self.shape_table.column(self.title_list[j], width=widths[j], minwidth=widths[j])
        self.grid_columnconfigure(0, minsize=sum(widths))
        self.widths = widths

    def get_focus_ix(self):
        '''
        Returns the index of the focus row (lowest index = 0, not 1)
        '''
        x = self.shape_table.focus()
        current_row = self.shape_table.item(x)
        current_row = current_row['values']
        if current_row == "":
            return None
        return int(current_row[0]) - 1

    def edit_shape(self, event=None):
        '''
        Opens a window for editing the selected ele_shape
        '''
        # Open new Toplevel window for shape editing
        win = tk.Toplevel(self, class_='Tao')
        win.grid_columnconfigure(0, weight=1)
        widget_frame = tk.Frame(win)
        widget_frame.grid(row=0, column=0, sticky='NSEW')
        button_frame = tk.Frame(win)
        button_frame.grid(row=1, column=0, sticky='NSEW')
        # Get currently selected row
        current_ix = self.get_focus_ix()
        current_row = self.shape_list[current_ix]

        # Create widgets
        title = "ele_shape(" + str(current_ix+1) + ")"
        win.title(title)
        params = []
        def wid_maker(i):
            '''Helper function'''
            return tao_parameter(self.keys[i], types[i-1], "T", current_row[i])
        types = ["STR", "ENUM", "ENUM", "REAL", "ENUM", "LOGIC", "LOGIC", "INT"]
        for i in range(1, len(current_row)):
            params.append(tk_tao_parameter(wid_maker(i), widget_frame, self.pipe,
                    prefix="shape" if types[i-1]=="ENUM" else ""))
            params[i-1].tk_label.grid(row=i, column=0, sticky='W')
            params[i-1].tk_wid.grid(row=i, column=1, sticky='EW')

        # Apply button
        def write_to_table():
            '''
            Writes the contents of this window into the table,
            makes the necessary changes in tao, and closes this window
            '''
            tao_set(params, "set " + self.type + " ele_shape(" + str(current_ix+1) + ")%", self.pipe)
            win.destroy()
            self.refresh()

        b = tk.Button(button_frame, text="Apply", command=write_to_table)
        b.pack()

    def swap_shapes(self, ix1, ix2):
        '''
        Swaps the shapes with indices ix1 and ix2
        (0 based) and refreshes the window
        '''
        # Input validation
        if not (isinstance(ix1, int) and isinstance(ix2, int)):
            return
        if ix1 == ix2:
            return
        if ix1 >= len(self.shape_list) or ix2 >= len(self.shape_list):
            return
        # Swap the shapes
        self.shape_list[ix1][0] = str(ix2+1)
        self.shape_list[ix2][0] = str(ix1+1)
        cmd_str = self.shape_set_format.format(**dict(zip(self.keys, self.shape_list[ix1])))
        self.pipe.cmd_in(cmd_str)
        cmd_str = self.shape_set_format.format(**dict(zip(self.keys, self.shape_list[ix2])))
        self.pipe.cmd_in(cmd_str)
        self.refresh()

    def move_up(self, event=None):
        '''
        Moves the selected item up one position, and shifts
        the others down by one
        '''
        # Get currently selected row
        current_ix = self.get_focus_ix()
        if current_ix == None:
            return
        self.swap_shapes(current_ix, current_ix-1)

    def move_down(self, event=None):
        '''
        Moves the selected item down one position, and shifts
        the others down by one
        '''
        # Get currently selected row
        current_ix = self.get_focus_ix()
        if current_ix == None:
            return
        self.swap_shapes(current_ix, current_ix+1)

    def duplicate_item(self, event=None):
        '''
        Copies the selected shape and adds a duplicate of it to
        the shape table at the next index
        '''
        ix = self.get_focus_ix()
        if ix == None:
            return
        # Add a new row in tao
        self.pipe.cmd_in("python shape_manage " + self.type + " "
                + str(ix+2) + " add")
        # Set the new row to be a copy of the selected row
        self.shape_list[ix][0] = str(ix+2)
        self.pipe.cmd_in(self.shape_set_format.format(**dict(zip(self.keys, self.shape_list[ix]))))
        self.refresh()

    def delete_item(self, event=None):
        '''
        Removes the selected item from the ele_shape table
        and shifts other ele_shapes accordingly
        '''
        ix = self.get_focus_ix()
        if ix==None:
            return
        self.pipe.cmd_in("python shape_manage " + self.type + " " + str(ix+1) + " delete")
        self.refresh()


    def new_at_top(self, event=None):
        '''
        Adds a new empty shape to the top of the ele_list
        and shifts other elements down
        '''
        # Add the shape
        self.pipe.cmd_in("python shape_manage " + self.type + " 1 add")
        # Set the ele_id non-empty to make it show up in the table
        cmd_str = self.shape_set_format.format(shape_ix=1, ele_id="None", shape="box",
                color="black", shape_size=0, label="none", draw="F", multi="F", line_width=0)
        self.pipe.cmd_in(cmd_str)
        # Refresh table
        self.refresh()

    def new_at_bottom(self, event=None):
        '''
        Adds a new empty shape to the bottom of the ele_list
        and shifts other elements down
        '''
        ix = len(self.shape_list)+1
        # Add the shape
        self.pipe.cmd_in("python shape_manage " + self.type
                + " " + str(ix) + " add")
        # Set the ele_id non-empty to make it show up in the table
        cmd_str = self.shape_set_format.format(shape_ix=ix, ele_id="None", shape="box",
                color="black", size=0, label="none", draw="F", multi="F", line_width=0)
        self.pipe.cmd_in(cmd_str)
        # Refresh table
        self.refresh()
