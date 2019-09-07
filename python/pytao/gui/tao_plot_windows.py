'''
Provides windows for viewing and editing plots in tao
'''
import tkinter as tk
import ttk
from tkinter import messagebox

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backends._backend_tk import FigureManagerTk
from matplotlib.backend_bases import key_press_handler
from matplotlib.widgets import Slider

from pytao.util.parameters import str_to_tao_param, tao_parameter_dict
from .tao_widget import *
from .taoplot import taoplot
from .tao_set import *
from .tao_base_windows import *
from .tao_lat_windows import tao_ele_window
from .tao_ele_location import in_element
from .tao_mpl_toolbar import taotoolbar
#-----------------------------------------------------
# Plot placing window

class tao_place_plot_window(Tao_Toplevel):
    '''
    Allows the user to choose from defined template
    plots and plot them in matplotlib windows
    Currently only supported in matplotlib mode
    '''
    def __init__(self, root, pipe, *args, **kwargs):
        self.root = root
        self.tao_id = 'plot'
        Tao_Toplevel.__init__(self, root, *args, **kwargs)
        self.title('Choose Plot')
        self.pipe = pipe
        self.grid_rowconfigure(0, weight=1)
        self.list_frame = tk.Frame(self)
        #self.list_frame.pack(side='top', fill='both', expand=1)
        self.list_frame.grid(row=0, column=0, sticky='NSEW')
        self.button_frame = tk.Frame(self)
        self.button_frame.grid(row=1, column=0, sticky='EW')
        self.button_frame.grid_columnconfigure(0, weight=1)
        self.button_frame.grid_columnconfigure(1, weight=1)
        #self.button_frame.pack(side='bottom', fill='x', expand=1)
        self.refresh()
        tk.Button(self.button_frame, text="Edit Template", command=self.edit_template).grid(row=0, column=0, sticky='EW')
        tk.Button(self.button_frame, text="Plot!", command=self.mpl_plot).grid(row=0, column=1, sticky='EW')

    def refresh(self):
        '''
        Responsible for creating widgets and filling them with plots
        '''
        for child in self.list_frame.winfo_children():
            child.destroy()

        # List of plots w/descriptions
        plots = self.pipe.cmd_in("python plot_list t")
        plots = plots.splitlines()
        for i in range(len(plots)):
            # get the description
            plot = plots[i].split(';')[1]
            plot_info = self.pipe.cmd_in('python plot1 ' + plot)
            plot_info = plot_info.splitlines()
            for line in plot_info:
                if line.find('description;') == 0:
                    d = line.split(';')[3]
                    break
                else:
                    d = ""
            plots[i] = [plot, d]
        widths = [0,0] # track column widths

        # Create list
        titles = ['Name', 'Description']
        self.tree = ttk.Treeview(
                self.list_frame, columns=titles, show='headings')
        # Column titles
        for title in titles:
            self.tree.heading(title, text=title)
            self.tree.column(title, stretch=True, anchor='w')

        # Fill rows
        for plot in plots:
            self.tree.insert("", "end", values=plot)
            for j in range(len(plot)):
                if len(plot[j])*10 > widths[j]:
                    widths[j] = len(plot[j])*10

        # Set column widths appropriately
        for j in range(len(titles)):
            if len(titles[j])*10 > widths[j]:
                widths[j] = len(titles[j])*10
            self.tree.column(titles[j], width=widths[j], minwidth=widths[j])

        # Scrollbars
        vbar = ttk.Scrollbar(
                self.list_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vbar.set)

        vbar.pack(side="right", fill="y", expand=0)
        self.tree.pack(side="left", fill="both", expand=1)
        self.widths = widths

        # Double click to open plot
        self.tree.bind('<Double-Button-1>', self.mpl_plot)
        # Single click to set self.plot
        self.tree.bind('<Button-1>', self.set_plot)

        tot = 0
        for w in widths:
            tot = tot+w
        self.maxsize(1800, 1000)
        self.minsize(tot, 100)

    def mpl_plot(self, event=None):
        '''
        Creates a new matplotlib window with the currently selected plot
        '''
        if self.set_plot():
            messagebox.showwarning('Error', "Please choose a template to plot.", parent=self)
            return
        self.set_plot()
        win = tao_plot_window(self.root, self.plot, self.pipe)

    def set_plot(self, event=None):
        '''
        Sets self.plot to the selected plot in the list
        '''
        x = self.tree.focus()
        if x == "":
            return 1
        row = self.tree.item(x)
        self.plot = row['values'][0]
        return 0

    def edit_template(self, event=None):
        '''
        Opens up a plot editting window and loads the selected template
        '''
        if self.set_plot():
            messagebox.showwarning('Error', "Please choose a template to edit.", parent=self)
            return
        win = tao_plot_tr_window(self.root, self.pipe, 'T')
        win.temp.set(self.plot)
        win.refresh()


#-----------------------------------------------------
# Plot editting window
#TODO: source input validation (base model design etc)

class tao_plot_tr_window(tao_list_window):
    '''
    Displays information about existing plot templates and regions
    Use a drop-down list to select template to view.
    '''
    def __init__(self, root, pipe, mode, *args, **kwargs):
        self.tao_id = 'plot'
        tao_list_window.__init__(self, root, "Edit Plots", *args, **kwargs)
        self.minsize(325, 100)
        self.pipe = pipe
        self.mode = mode

        # Widgets for swapping between templates and plots
        self.temp_frame = tk.Frame(self)
        self.temp_frame.pack(fill="both", expand=0)
        self.temp_frame.columnconfigure(0, pad=10)

        # Template/Active Plots
        tk.Label(self.temp_frame, text="Show:").grid(row=0, column=0, sticky='W')
        self.mode_var = tk.StringVar()
        if self.mode == "T":
            self.mode_var.set('Templates')
        elif self.mode == 'R':
            self.mode_var.set('Active Plots')
        tk.OptionMenu(self.temp_frame, self.mode_var, 'Templates', 'Active Plots',
                command = self.swap_mode).grid(row=0, column=1, sticky='EW')

        # Plot chooser
        tk.Label(self.temp_frame, text="Choose plot:").grid(row=1, column=0, sticky='W')
        self.temp = tk.StringVar()
        self.temp_select = ttk.Combobox(
                self.temp_frame, textvariable=self.temp, values=[])
        self.temp_select.bind("<<ComboboxSelected>>", self.refresh)
        self.temp_select.bind("<Return>", self.refresh)
        self.temp_select.grid(row=1, column=1)

        # Apply/transfer buttons
        self.button_frame = tk.Frame(self)
        self.button_frame.grid_columnconfigure(0, weight=1)
        self.button_frame.grid_columnconfigure(1, weight=1)
        self.button_frame.pack(fill="both", expand=0)
        apply_button = tk.Button(
                self.button_frame, text="Apply changes", command=self.plot_apply)
        apply_button.grid(row=0, column=0, sticky='EW')

        self.transfer_text = tk.StringVar()
        self.transfer_button = tk.Button(self.button_frame,
                textvariable=self.transfer_text, command=self.plot_transfer)
        self.transfer_button.grid(row=0, column=1, sticky='EW')

        self.swap_mode(overide=True) # also runs self.refresh()

    def swap_mode(self, overide=False, event=None):
        '''
        Swaps between showing templates and active plots
        Will always run if overide is set to True
        '''
        # Switch self.mode
        mode_dict = {'Templates':'T', 'Active Plots':'R'}
        new_mode = mode_dict[self.mode_var.get()]
        if (self.mode == new_mode) and not overide:
            return #no action necessary
        self.mode = new_mode

        # Template plots
        t_plot_list = self.pipe.cmd_in("python plot_list t")
        t_plot_list = t_plot_list.splitlines()
        t_index_list = len(t_plot_list)*[0] #get correct length
        for i in range(len(t_plot_list)):
            t_index_list[i], t_plot_list[i] = t_plot_list[i].split(';')

        # Active plots
        r_plot_list = self.pipe.cmd_in("python plot_list r")
        new_r_plot_list = []
        r_index_list = []
        r_region_list = [] #needed to open the correct graph/curve windows
        r_plot_list = r_plot_list.splitlines()
        for i in range(len(r_plot_list)):
            if r_plot_list[i].split(';')[2] != "": #region contains a plot
                new_r_plot_list.append(r_plot_list[i].split(';')[2])
                r_index_list.append(r_plot_list[i].split(';')[0])
                r_region_list.append(r_plot_list[i].split(';')[1])
        r_plot_list = new_r_plot_list

        # Populate self.plot_list and self.index_list
        if self.mode == "T":
            self.plot_list = t_plot_list
            self.plot_display_list = t_plot_list
            self.index_list = t_index_list
            self.transfer_text.set('Transfer properties to active plots')
        elif self.mode == "R":
            self.plot_list = r_plot_list
            self.plot_display_list = []
            self.index_list = r_index_list
            self.region_list = r_region_list
            for i in range(len(self.plot_list)):
                self.plot_display_list.append(
                        self.plot_list[i] + " (" + self.region_list[i] + ")")
            self.transfer_text.set('Transfer properties to template')
        # Save t_plot_list and r_plot_list for use in self.refresh
        self.t_plot_list = t_plot_list
        self.r_plot_list = r_plot_list

        self.temp_select.configure(values=self.plot_display_list)
        if self.plot_list == []:
            for child in self.list_frame.winfo_children():
                child.destroy()
            tk.Label(self.list_frame, text="NO PLOTS FOUND").pack()
            return
        self.temp.set(self.plot_display_list[0])
        self.refresh()

    def refresh(self, event=None):
        '''
        Clears self.list_frame and populates it with information relevant to
        self.temp
        '''
        # Don't run if self.temp is not a valid template
        if self.temp.get() not in self.plot_display_list:
            return

        # Clear list_frame
        for child in self.list_frame.winfo_children():
            child.destroy()

        # Get info on self.temp
        self.plot = self.temp.get()
        self.plot = self.plot_list[self.plot_display_list.index(self.plot)]
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
        name.tk_label.grid(row=0, column=0, sticky='W')
        name.tk_wid.grid(row=0, column=1, sticky='EW')
        if num_graphs > 0:
            tk.Label(self.list_frame, text="Graphs").grid(
                    row=1, column=0, rowspan=num_graphs, sticky='W')
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
            data_list[i].tk_label.grid(row=i+1+num_graphs, column=0, sticky='W')
            data_list[i].tk_wid.grid(row=i+1+num_graphs, column=1, sticky='EW')
        self.tao_list = data_list

        # Activate/deactivate transfer button
        if ((self.plot in self.t_plot_list) and (self.mode == 'R')
                or (self.plot in self.r_plot_list) and (self.mode == 'T')):
            self.transfer_button.configure(state='normal')
        else:
            self.transfer_button.configure(state='disabled')


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
        win = tao_plot_graph_window(self.root, graph, self, self.pipe, self.mode, index)

    def plot_apply(self):
        index = self.index_list[self.plot_list.index(self.plot)]
        set_str = "set plot @" + self.mode + str(index) + ' '
        tao_set(self.tao_list, set_str, self.pipe)
        ## In matplotlib mode, apply for the appropriate region too
        #c1 = self.root.plot_mode == "matplotlib"
        #c2 = self.mode == "T"
        #c3 = self.plot in self.root.placed.keys()
        #if c1 & c2 & c3:
        #  set_str = "set plot " + self.root.placed[self.plot] + ' '
        #  tao_set(self.tao_list, set_str, self.pipe)
        # Refresh any existing plot windows
        for win in self.root.refresh_windows['plot']:
            #if isinstance(win, tao_plot_window):
            #  # only refresh relevant plot windows
            #  if self.plot == win.region:
            #        win.refresh()
            #else:
            #  win.refresh()
            win.refresh()

    def plot_transfer(self):
        '''
        Transfer properties between templates and active plots
        '''
        ix = self.plot_list.index(self.plot)
        ix = self.index_list[ix]
        if self.mode == 'T':
            self.pipe.cmd_in('python plot_transfer @T' + str(ix))
        elif self.mode == 'R':
            self.pipe.cmd_in('python plot_transfer @R' + str(ix))
        self.refresh()

#----------------------------------------------------
# Matplotlib plotting window

class tao_plot_window(Tao_Toplevel):
    '''
    Displays one (perhaps multiple) matplotlib plots
    that the user has specified from the plotting
    template window that they want to plot. Creating a
    window in tkinter is necessary rather than using
    matplotlib's built in system for creating windows
    because using that system will halt the tkinter
    mainloop until the plots are closed.
    If the region to place the graph is not specified, one will be selected automatically
    '''
    def __init__(self, root, template, pipe, region=None, *args, **kwargs):
        if region == 'layout': # do not place plots in the layout region
            return
        if template == "key_table":
            messagebox.showwarning("Warning", "Key table not available in the GUI")
            return
        # verify that the graphs for this template are valid
        # must place the template first to get accurate info
        tmp_reg = root.placed.place_template(template)
        plot1 = pipe.cmd_in('python plot1 ' + tmp_reg).splitlines()
        valid = True
        for i in range(str_to_tao_param(plot1[0]).value):
            plot_graph = pipe.cmd_in("python plot_graph " + tmp_reg + '.'
                    + str_to_tao_param(plot1[i+1]).value).splitlines()
            if not tao_parameter_dict(plot_graph)['valid'].value:
                valid = False
                break
        root.placed.unplace_region(tmp_reg)
        if not valid:
            messagebox.showwarning("Warning", "The plot you have selected ("
                    + template + ") has one or more invalid graphs.")
            return
        self.root = root
        self.tao_id = 'plot'
        Tao_Toplevel.__init__(self, root, *args, **kwargs)
        self.template = template #The template plot being plotted
        self.pipe = pipe
        self.fig = False #Default value

        self.region = self.root.placed.place_template(self.template, region)

        self.mpl = taoplot(pipe, self.region)
        self.title(template + ' (' + self.region + ')')
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
        #canvas.get_tk_widget().pack(side="top", fill="both", expand=1)
        toolbar = taotoolbar(canvas, self, width, self.root.GUI_DIR)
        toolbar.update()
        # DO NOT TOUCH
        canvas.manager = FigureManagerTk(
                canvas, self.fig.number, tk.Toplevel(self.root))

        #toolbar = taotoolbar(canvas, self)
        #toolbar.update()
        #canvas._tkcanvas.pack(side="top", fill="both", expand=1)

        def on_key_press(event):
            key_press_handler(event, canvas, toolbar)

        canvas.mpl_connect("key_press_event", on_key_press)

        def on_click(event):
            '''opens element information window if an element is double clicked on'''
            if event.dblclick:
                eleList = in_element(event.xdata,event.ydata,self.fig_info)
                for i in eleList:
                    tao_ele_window(self.root,self.pipe,
                            default=[self.fig_info[1],i[0],i[1],self.fig_info[3]])

        canvas.mpl_connect("button_press_event", on_click)
        '''
        if self.fig_info[0] == 'floor_plan':
            self.fig.subplots_adjust(bottom=0.2) #adds room below graph for slider
            width_slider = Slider(self.fig.add_axes([.1,.05,.8,.05]), 'width', 0, 2, width) #element width slider

            def update_slider(width):
                self.refresh(width=width_slider.val)

            width_slider.on_changed(update_slider) #call update when slider moves
        '''
        self.update_idletasks()
        self.pack_propagate(False)

    def destroy(self):
        # Clear self.region
        self.root.placed.unplace_region(self.region)
        Tao_Toplevel.destroy(self)

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
    def __init__(self, root, graph, parent, pipe, mode, index, *args, **kwargs):
        tao_list_window.__init__(self, root, graph, parent=parent, *args, **kwargs)
        # Set the title properly
        if graph.split('.')[0] in self.root.placed.keys():
            title = self.root.placed[graph.split('.')[0]]
            title += '.' + graph.split('.')[1]
            title += " (" + graph.split('.')[0] + ")"
            self.title(title)
        self.transient(parent)
        #self.grab_set()
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

        self.wait_window()

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
        name.tk_label.grid(row=0, column=0, sticky='W')
        name.tk_wid.grid(row=0, column=1, sticky='EW')

        # Curve buttons
        if num_curves > 0:
            tk.Label(self.list_frame, text="Curves").grid(
                    row=1, column=0, rowspan=num_curves, sticky='W')
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
            data_list[i].tk_label.grid(row=i+1+num_curves, column=0, sticky='W')
            data_list[i].tk_wid.grid(row=i+1+num_curves, column=1, sticky='EW')

    def open_curve_callback(self, graph, curve):
        return lambda : self.open_curve(graph, curve)

    def open_curve(self, graph, curve):
        '''
        Opens a window to display info for a given curve
        '''
        full_curve = graph + '.' + curve
        win = tao_plot_curve_window(self.root, full_curve, self, self.pipe)
        # Set the title properly
        if graph.split('.')[0] in self.root.placed.keys():
            title = self.root.placed[graph.split('.')[0]]
            title += '.' + graph.split('.')[1] + '.' + curve
            title += " (" + graph.split('.')[0] + ")"
            win.title(title)

    def curve_apply(self, win):
        #Convert from plot.graph.curve to index.graph.curve
        curve_name = str(self.index) + win.curve[win.curve.index('.'):]
        set_str = "set curve @" + self.mode + curve_name + " "
        tao_set(win.tao_list, set_str, win.pipe)
        ## In matplotlib mode, apply for the appropriate region too
        #c1 = self.root.plot_mode == "matplotlib"
        #c2 = self.mode == "T"
        #plot = win.curve[:win.curve.index('.')]
        #c3 = plot in self.root.placed.keys()
        #if c1 & c2 & c3:
        #  curve_name = self.root.placed[win.curve[:win.curve.index('.')]] \
        #            + win.curve[win.curve.index('.'):]
        #  set_str = "set curve " + curve_name + ' '
        #  tao_set(win.tao_list, set_str, win.pipe)
        # Refresh any existing plot windows
        for win in self.root.refresh_windows['plot']:
            win.refresh()

    def graph_apply(self):
        #Convert from plot.graph to index.graph
        graph_name = str(self.index) + self.graph[self.graph.index('.'):]
        set_str = "set graph @" + self.mode + graph_name + ' '
        tao_set(self.tao_list, set_str, self.pipe)
        ## In matplotlib mode, apply for the appropriate region too
        #c1 = self.root.plot_mode == "matplotlib"
        #c2 = self.mode == "T"
        #plot = self.graph[:self.graph.index('.')]
        #c3 = plot in self.root.placed.keys()
        #if c1 & c2 & c3:
        #  graph_name = self.root.placed[self.graph[:self.graph.index('.')]] \
        #            + self.graph[self.graph.index('.'):]
        #  set_str = "set graph " + graph_name + ' '
        #  tao_set(self.tao_list, set_str, self.pipe)
        # Refresh any existing plot windows
        for win in self.root.refresh_windows['plot']:
            #if isinstance(win, tao_plot_window):
            #  # only refresh plot windows for this graph
            #  if self.graph[self.graph.index('.'):] == win.region:
            #        win.refresh()
            #else:
            #  win.refresh()
            win.refresh()

#-----------------------------------------------------
# plot_curve window

class tao_plot_curve_window(tao_parameter_window):
    '''
    Displays info for a given curve
    '''
    def __init__(self, root, curve, parent, pipe, *args, **kwargs):
        # Get the parameters
        self.curve = curve
        self.pipe = pipe
        data_list = self.pipe.cmd_in("python plot_curve " + curve)
        data_list = data_list.splitlines()
        for i in range(len(data_list)):
            data_list[i] = str_to_tao_param(data_list[i])

        tao_parameter_window.__init__(
                self, root, curve, data_list, self.pipe, plot=curve.split('.')[0], parent=parent, *args, **kwargs)

        b = tk.Button(self.button_frame, text="Apply changes",
                command=lambda : parent.curve_apply(self))
        b.pack()
        self.transient(parent)
        #self.grab_set()
        #self.wait_window(self)


class tao_new_plot_template_window(Tao_Toplevel):
    '''
    Provides a window for creating new plot templates (and their
    associated graphs and curves)
    Pass the name of an existing plot to open that plot and start
    editing its graphs and curves
    '''
    def __init__(self, root, pipe, default=None, *args, **kwargs):
        self.root = root
        Tao_Toplevel.__init__(self, root, *args, **kwargs)
        self.pipe = pipe
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=2)
        self.rowconfigure(1, weight=1)
        self.title('New Plot Template')
        self.name = ""
        self.x_axis_type = 'index'

        # Frame for inputting plot parameters
        self.plot_frame = tk.Frame(self)
        self.plot_frame.grid(row=1, column=0, sticky='NSEW')

        self.fill_plot_frame()
        # X-axis/graph type compatibility dict
        self.compat_dict = {}
        self.compat_dict['index'] = self.compat_dict['ele_index'] \
                = self.compat_dict['lat'] = self.compat_dict['var'] \
                = self.compat_dict['data'] = ['data']
        self.compat_dict['s'] = ['data', 'lat_layout']
        self.compat_dict['floor'] = ['floor_plan']
        self.compat_dict['phase_space'] = [
                'dynamic_aperture', 'histogram', 'phase_space']
        self.compat_dict['histogram'] = ['histogram']
        self.compat_dict['none'] = ['floor_plan'] #FIXME
        # Graph frame
        self.graph_frame = tabbed_frame(self, lambda arg : new_graph_frame(arg, self))
        self.graph_frame.grid(row=1, column=1, sticky='NSEW')
        if default != None:
            self.plot_param_list[0].tk_var.set(default)
            self.clone_plot(ask=False)

    def fill_plot_frame(self):
        tk.Label(self, text="New Plot Template",
                font=('Sans', 16, 'bold')).grid(
                        row=0, column=0, columnspan=2, sticky='EW')

        # Small bit of setup
        def my_ttp(x):
            ''' shorcut for the following commonly used construct '''
            return tk_tao_parameter(str_to_tao_param(x), self.plot_frame, self.pipe)

        # Widgets
        params = ["name;STR;T;",
                "description;STR;T;",
                "x_axis_type;ENUM;T;",
                "autoscale_gang_x;LOGIC;T;T",
                "autoscale_gang_y;LOGIC;T;T",
                "autoscale_x;LOGIC;T;F",
                "autoscale_y;LOGIC;T;F",
                "n_curve_pts;INT;T;"]
        self.plot_param_list = list(map(my_ttp, params))

        # Parameter index lookup (for convenience)
        self.ixd = {} #index dictionary
        for i in range(len(self.plot_param_list)):
            self.ixd[self.plot_param_list[i].param.name] = i

        # Labels
        def plot_label_maker(x):
            ''' Helper function '''
            return tk.Label(self.plot_frame, text=x)
        labels = ["Plot Name:",
                "Description:",
                "X-axis Type:",
                "Autoscale gang X:",
                "Autoscale gang Y:",
                "Autoscale X:",
                "Autoscale Y:",
                "Number of curve points:"]
        self.plot_label_list = list(map(plot_label_maker, labels))

        # Grid widgets and labels
        for i in range(len(self.plot_param_list)):
            self.plot_label_list[i].grid(row=i+1, column=0, sticky='W')
            self.plot_param_list[i].tk_wid.grid(row=i+1, column=1, sticky='EW')

        # Warning labels
        self.name_warning_1 = tk.Label(self.plot_frame, text="Cannot be empty")
        self.name_warning_2 = tk.Label(self.plot_frame, text="Cannot contain whitespace")

        # Responses to edits
        self.plot_param_list[0].tk_wid.bind('<FocusOut>', self.plot_name_handler)
        self.plot_param_list[self.ixd['x_axis_type']].tk_var.trace('w', self.x_axis_type_handler)

        # Clone existing plot
        existing_plot_templates = self.pipe.cmd_in("python plot_list t").splitlines()
        for i in range(len(existing_plot_templates)):
            existing_plot_templates[i] = existing_plot_templates[i].split(';')[1]
        existing_plot_templates = ['None'] + existing_plot_templates
        self.clone_plot = tk.StringVar()
        self.clone_plot.set('None')
        self.clone_chooser = ttk.Combobox(self.plot_frame,
                textvariable=self.clone_plot, values=existing_plot_templates, state='readonly')
        self.clone_b = tk.Button(self.plot_frame, text="Clone", command=self.clone_plot)

        tk.Label(self.plot_frame, text="Clone existing plot:").grid(row=i+2, column=0, sticky='W')
        self.clone_chooser.grid(row=i+2, column=1, sticky='EW')
        self.clone_b.grid(row=i+2, column=2, sticky='W')

        # Focus the name entry
        self.plot_param_list[0].tk_wid.focus_set()

    def x_axis_type_handler(self, *args):
        '''
        Updates self.x_axis_type to match the chosen x_axis_type, then calls
        graph_type_handler for each graph
        '''
        self.x_axis_type = self.plot_param_list[2].tk_var.get()
        for graph in self.graph_frame.tab_list:
            graph.graph_type_handler()

    def clone_plot(self, plot_name, ask=True):
        '''
        Clone the plot specified by plot_name
        '''
        clone_graphs = []
        if plot_name == 'None':
            return
        ans_var = tk.StringVar()
        # Ask if user wants to keep existing graphs
        if ask:
            tao_message_box(self.root, self, ans_var, title='Warning', message='Would you like to keep or discard the graphs you defined for ' + self.name + '?', choices=['Keep', 'Discard'])
        else:
            ans_var.set('Discard')
        if ans_var.get() == 'Keep':
            c4 = False
        elif ans_var.get() == 'Discard':
            c4 = True
        else:
            return
        # Specified plot to clone
        plot1 = tao_parameter_dict(self.pipe.cmd_in('python plot1 ' + plot_name).splitlines())
        num_graphs = plot1['num_graphs'].value
        for i in range(1, num_graphs+1):
            clone_graphs.append(plot1['graph['+str(i)+']'].value)
        # Ask what to do about plot properties
        use_plot_props = None
        msg = "Would you like to use the plot-level properties of "
        msg += plot_name + ", or use the properties you have specified?"
        choices = []
        choices.append("Use properties of " + plot_name)
        choices.append("Use the properties I have specified")
        tao_message_box(self.root, self, ans_var, title='Plot parameters', message=msg, choices=choices)
        if ans_var.get() == choices[0]:
            use_plot_props = plot_name
        # Copy in plot properties if necessary
        if use_plot_props != None:
            plot1 = tao_parameter_dict(
                    self.pipe.cmd_in("python plot1 " + use_plot_props).splitlines())
            for w in self.plot_param_list:
                if w.param.name in plot1.keys():
                    w.param_copy(plot1[w.param.name])
        # Copy the graphs
        for i in range(len(clone_graphs)):
            self.graph_frame.add_tab(i)
            graph = self.clone_plot.get() + '.' + graph
            self.graph_frame.tab_list[i].clone(graph)
        # Delete old graphs if requested
        if c4:
            for j in range(i+1, len(self.graph_frame.tab_list)):
                self.graph_frame.remove_tab(i+1, destroy=True)
        # Switch to first tab
        self.graph_frame.notebook.select(0)

    def create_template(self, event=None):
        '''
        Takes the information from the d2_frame and d1_frames and runs
        the necessary commands to create the data in tao, then closes
        the create data window
        '''
        # Input validation (more TODO)
        messages = []
        for graph_frame in self.graph_frame_list:
            # Check names
            if graph_frame.name_handler():
                name_m = "Please check graph names."
                if name_m not in messages:
                    messages.append(name_m)
            # Check for semicolons in any fields
            semi_message = "Semicolons not allowed in any input field"
            caret_message = "Carets not allowed in any input field"
            curve_name_m = "Curve names cannot contain whitespace"
            broken = False #Used to break out of the below for loops
            # Check for semicolons/carets
            for ttp in graph_frame.graph_wids:
                if str(ttp.tk_var.get()).find(';') != -1:
                    messages.append(semi_message)
                    broken = True
                    break
                if str(ttp.tk_var.get()).find('^') != -1:
                    messages.append(caret_message)
                    broken = True
                    break
            for curve_list in graph_frame.curve_dict.values():
                if broken:
                    break
                for ttp in curve_list:
                    if str(ttp.tk_var.get()).find(';') != -1:
                        messages.append(semi_message)
                        broken = True
                        break
                    if str(ttp.tk_var.get()).find('^') != -1:
                        messages.append(caret_message)
                        broken = True
                        break
                    if (str(ttp.tk_var.get()).find(' ') != -1) and (ttp.param.name == 'name'):
                        messages.append(curve_name_m)
                        broken = True
                        break
                if broken:
                    break
        for m in messages:
            messagebox.showwarning("Error", m, parent=self)
        if messages != []:
            return
        # Get the appropriate template index
        plot_list_t = self.pipe.cmd_in('python plot_list t').splitlines()
        indices = [""]*len(plot_list_t)
        names = [""]*len(plot_list_t)
        for i in range(len(plot_list_t)):
            indices[i] = plot_list_t[i].split(';')[0]
            names[i] = plot_list_t[i].split(';')[1]
        if self.name in names:
            template_ix = indices[names.index(self.name)]
        else:
            template_ix = str(int(indices[-1]) + 1)
        # Create the template
        n_graph = str(len(self.graph_frame_list))
        cmd_str = 'python plot_manage_plot @T' + template_ix + '^^'
        cmd_str += self.name + '^^' + n_graph
        for graph_frame in self.graph_frame_list:
            cmd_str += '^^' + graph_frame.name
        self.pipe.cmd_in(cmd_str)
        # Set the plot properties (but not name)
        set_list = []
        for ttp in self.plot_param_list:
            if ttp.param.name != 'name':
                set_list.append(ttp)
        tao_set(set_list, "set plot @T" + template_ix + ' ', self.pipe)
        # Set graph properties
        for gf in self.graph_frame_list:
            graph_name = "@T" + template_ix + '.' + gf.name
            # Set graph properties (but not name or n_curve)
            set_list = []
            for ttp in gf.graph_wids:
                if ttp.param.name not in ['name', 'n_curve']:
                    set_list.append(ttp)
            tao_set(set_list, "set graph " + graph_name + ' ', self.pipe)
            # Create curves
            curve_nums = range(1, gf.n_curve+1)
            for c in curve_nums:
                curve_name = gf.curve_dict[c][0].tk_var.get().strip()
                if curve_name == "": #default to c1, c2, etc
                    curve_name = "c" + str(c)
                self.pipe.cmd_in('python plot_manage_curve ' + graph_name
                        + '^^' + str(c) + '^^' + curve_name)
                curve_name = graph_name + '.' + curve_name
                # Set curve properties (but not the name)
                set_list = []
                for ttp in gf.curve_dict[c]:
                    if ttp.param.name != 'name':
                        set_list.append(ttp)
                tao_set(set_list, 'set curve ' + curve_name + ' ', self.pipe)
        # Close the window
        self.destroy()
        # Refresh plot-related windows
        for win in self.root.refresh_windows['plot']:
            win.refresh()

    def plot_name_handler(self, event=None):
        '''
        Reads the plot name into self.name, and warns the user if left blank
        '''
        name = self.plot_param_list[self.ixd['name']].tk_var.get().strip()
        self.name_warning_1.grid_forget()
        self.name_warning_2.grid_forget()
        if name == "":
            self.name_warning_1.grid(row=self.ixd['name']+1, column=2, sticky='W')
        elif name.find(' ') != -1:
            self.name_warning_2.grid(row=self.ixd['name']+1, column=2, sticky='W')
        else:
            self.title(name)
            self.name = name

    def x_axis_type_handler(self, *args):
        '''
        Sets self.x_axis_type to the selected x_axis_type
        '''
        self.x_axis_type = self.plot_param_list[self.ixd['x_axis_type']].tk_var.get()
        # Run graph_type_handler for each child graph
        for graph in self.graph_frame.tab_list:
            graph.graph_type_handler()


    def refresh(self):
        '''
        Only here in case something tries to refresh this window
        '''
        pass

class new_graph_frame(tk.Frame):
    '''
    New and improved frame for editing graph properties
    '''
    def __init__(self, parent, plot):
        tk.Frame.__init__(self, parent.notebook)
        self.parent = parent
        self.plot = plot
        self.pipe = plot.pipe
        self.handler_block = False
        self._uf = tk.Frame(self) # Used to make graph frame and curve notebook uniform width
        self._uf.grid(row=0, column=0, sticky='NSEW')
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)#, uniform="uf")
        self.grid_columnconfigure(1, weight=1)#, uniform="uf")
        self._scroll_frame = tao_scroll_frame(self._uf)
        self.name = "New_graph"
        # Default to first type compatible with self.plot.x_axis_type
        self.type = self.plot.compat_dict[self.plot.x_axis_type][0]

        #TEMPORARY
        def f(e):
            self.handler_block = False
        self.bind_all('<Control-b>', f)


        # Delete button
        tk.Button(self._uf, text="DELETE THIS GRAPH", fg='red', command=self.delete).grid(
                row=0, column=0, columnspan=3, sticky='EW')

        # Duplicate button
        tk.Button(self._uf, text="Duplicate this graph", command=self.duplicate).grid(
                row=1, column=0, columnspan=3, sticky='EW')

        self._uf.grid_columnconfigure(2, weight=1)

        # Setup
        self.graph_frame = self._scroll_frame.frame
        self.graph_frame.grid_columnconfigure(1, weight=1)
        # Helper functions
        def graph_ttp(x):
            ''' Shortcut for commonly used construct '''
            p = str_to_tao_param(x)
            if x in self.head_wids:
                return tk_tao_parameter(p, self._uf, self.pipe)
            else:
                return tk_tao_parameter(p, self.graph_frame, self.pipe)
        def qp_axis_props(x):
            '''Adds the props of a qp-axis-struct to the parameter string x'''
            if self.plot.root.plot_mode == 'matplotlib':
                x += ';label;STR;;min;REAL;;max;REAL;;draw_label;LOGIC;T;draw_numbers;LOGIC;T'
            else:
                x += ';label;STR;;min;REAL;;max;REAL;;number_offset;REAL;'
                x += ';label_offset;REAL;;label_color;ENUM;;major_tick_len;REAL;'
                x += ';minor_tick_len;REAL;;major_div;INT;;major_div_nominal;INT;'
                x += ';minor_div;INT;;minor_div_max;INT;;places;INT;;type;ENUM;'
                x += ';bounds;ENUM;;tick_side;ENUM;;number_side;ENUM;'
                x += ';draw_label;LOGIC;T;draw_numbers;LOGIC;T'
            return x
        def graph_label_maker(x):
            '''Helper function'''
            if x in self.head_labels:
                return tk.Label(self._uf, text=x)
            else:
                return tk.Label(self.graph_frame, text=x)

        # Graph widgets

        # Basics
        self.head_wids = ["name;STR;T;",
                "graph^type;ENUM;T;"+self.type,
                "title;STR;T;"]
        self.head_labels = ["Graph name:", "Graph type:", "Graph title"]
        self.head_wids = list(map(graph_ttp, self.head_wids))
        self.head_labels = list(map(graph_label_maker, self.head_labels))
        # Type-specific options
        self.wids = {}
        self.labels = {}
        # Data graphs
        self.wids['data'] = ["component;COMPONENT;T;model",
                "ix_universe;INT;T;",
                "draw_only_good_user_data_or_vars;LOGIC;T;T"]
        self.labels['data'] = ["Component:", "Universe:",
                "Draw only good_user data/variables:"]
        self.wids['data'] = list(map(graph_ttp, self.wids['data']))
        self.labels['data'] = list(map(graph_label_maker, self.labels['data']))
        # Floor plans
        self.wids['floor_plan'] = [
                qp_axis_props("x2;STRUCT;T"),
                "ix_universe;INT;T;",
                "floor_plan_rotation;REAL;T;",
                "floor_plan_view;STR;T;zx",
                "floor_plan_orbit_scale;REAL;T;",
                "floor_plan_orbit_color;STR;T;",
                "floor_plan_flip_label_side;LOGIC;T;F",
                "floor_plan_size_is_absolute;LOGIC;T;F",
                "floor_plan_draw_only_first_pass;LOGIC;T;F"]
        self.labels['floor_plan'] = ["X2-axis:", "Universe:",
                "Rotation:", "View:", "Orbit scale:", "Orbit color:",
                "Flip label side:", "Absolute size:", "Draw only first pass:"]
        self.wids['floor_plan'] = list(map(graph_ttp, self.wids['floor_plan']))
        self.labels['floor_plan'] = list(map(graph_label_maker, self.labels['floor_plan']))
        # Lat layouts
        self.wids['lat_layout'] = [
                "ix_universe;INT;T;",
                "ix_branch;INT;T;"]
        self.labels['lat_layout'] = ["Universe:", "Branch:"]
        self.wids['lat_layout'] = list(map(graph_ttp, self.wids['lat_layout']))
        self.labels['lat_layout'] = list(map(graph_label_maker, self.labels['lat_layout']))
        # Dynamic aperture
        self.wids['dynamic_aperture'] = ["ix_universe;INT;T;"]
        self.labels['dynamic_aperture'] = ["Universe:"]
        self.wids['dynamic_aperture'] = list(map(graph_ttp, self.wids['dynamic_aperture']))
        self.labels['dynamic_aperture'] = list(map(graph_label_maker, self.labels['dynamic_aperture']))
        # Histograms
        self.wids['histogram'] = ["ix_universe;INT;T;"]
        self.labels['histogram'] = ["Universe:"]
        self.wids['histogram'] = list(map(graph_ttp, self.wids['histogram']))
        self.labels['histogram'] = list(map(graph_label_maker, self.labels['histogram']))
        # Phase space plots
        self.wids['phase_space'] = ["ix_universe;INT;T;"]
        self.labels['phase_space'] = ["Universe:"]
        self.wids['phase_space'] = list(map(graph_ttp, self.wids['phase_space']))
        self.labels['phase_space'] = list(map(graph_label_maker, self.labels['phase_space']))
        # Graph styling
        self.style_wids = [
                qp_axis_props("x;STRUCT;T"),
                qp_axis_props("y;STRUCT;T"),
                qp_axis_props("y2;STRUCT;T"),
                "clip;LOGIC;T;T",
                "draw_axes;LOGIC;T;T",
                "draw_grid;LOGIC;T;T",
                "allow_wrap_around;LOGIC;T;F",
                "symbol_size_scale;REAL;T;",
                "correct_xy_distortion;LOGIC;T;F",
                "x_axis_scale_factor;REAL;T;"]
        self.style_labels = [
                "X-axis:", "Y-axis:", "Y2-axis",
                "Clip at boundary:", "Draw axes:", "Draw grid:",
                "Allow wrap-around:", "Symbol size scale:",
                "Correct xy distortion:", "X-axis scale factor:"]
        if self.plot.root.plot_mode != 'matplotlib':
            self.style_wids = [
                    "box;STR;T;",
                    "margin;REAL;T;",
                    "scale_margin;REAL;T;"] + self.style_wids
            self.style_labels = ["Box:", "Margin:", "Scale Margin:"] + self.style_labels
        self.style_wids = list(map(graph_ttp, self.style_wids))
        self.style_labels = list(map(graph_label_maker, self.style_labels))

        # Grid head widgets
        for i in range(len(self.head_wids)):
            self.head_labels[i].grid(row=i+2, column=0, sticky='W')
            self.head_wids[i].tk_wid.grid(row=i+2, column=1, sticky='EW')

        # Warning labels
        # (defined here to be gridded/ungridded as necessary)
        self.name_warning_1 = tk.Label(self._uf, text="Must not be empty")
        self.name_warning_2 = tk.Label(self._uf, text="Graph name already in use")
        self.name_warning_3 = tk.Label(self._uf, text="Cannot contain whitespace")

        # Responses to edits
        self.head_wids[0].tk_wid.bind('<FocusOut>', self.name_handler)
        self.head_wids[1].tk_var.trace('w', self.graph_type_handler)

        # Curves
        self.curve_frame = tabbed_frame(self, lambda arg : new_curve_frame(arg, self))

        # Element shapes (for lat_layouts and floor_plans)
        self.ele_frame = tk.Frame(self)
        tk.Label(self.ele_frame, text="COMING SOON").grid(row=0, column=0, sticky='NSEW')

        # Grid everything else
        self._scroll_frame.grid(row=2+len(self.head_wids), column=0, columnspan=3, sticky='NSEW')
        self._uf.grid_rowconfigure(2+len(self.head_wids), weight=1)
        self.refresh()
        self.update_idletasks()
        self.refresh() #called again to get the head widgets sized correctly

    def refresh(self):
        '''
        Grids the appropriate widgets for the current graph type,
        grids self.curve_frame or self.ele_frame as appropriate,
        then calls refresh for each of the curves if appropriate
        (i.e. if self.type not in ['lat_layout', 'floor_plan'])
        '''
        # Ungrid the non-style widgets
        for child in self.graph_frame.winfo_children():
            child.grid_forget()
        # Grid the appropriate widgets
        for i in range(len(self.wids[self.type])):
            self.labels[self.type][i].grid(row=i, column=0, sticky='W')
            self.wids[self.type][i].tk_wid.grid(row=i, column=1, sticky='EW')
        offset = i+1
        # Grid the style widgets
        for i in range(len(self.style_wids)):
            ix = i + offset
            self.style_labels[i].grid(row=ix, column=0, sticky='W')
            self.style_wids[i].tk_wid.grid(row=ix, column=1, sticky='EW')
        #self.update_idletasks() #let widgets obtain their sizes
        # Set label widths properly for head widgets
        label_width=0
        for child in self.graph_frame.grid_slaves():
            if isinstance(child, tk.Label):
                label_width = max(label_width, child.winfo_width())
        self._uf.grid_columnconfigure(0, minsize=label_width)
        # Swap between self.curve_frame or self.ele_frame as necessary
        if self.type in ['lat_layout', 'floor_plan']:
            self.curve_frame.grid_forget()
            self.ele_frame.grid(row=0, column=1, sticky='NSEW')
        else:
            self.ele_frame.grid_forget()
            self.curve_frame.grid(row=0, column=1, sticky='NSEW')
            for curve in self.curve_frame.tab_list:
                curve.refresh()

    def delete(self, ask=True, event=None):
        '''
        Deletes this graph frame
        Call with ask = False to skip confirmation
        '''
        # Ask for confirmation
        if ask:
            ans = messagebox.askokcancel("Delete " + self.name, "Delete this graph and its associated curves?", parent=self.parent)
            if not ans:
                return

        # Remove from tabbed frame
        self.parent.remove_tab(self.parent.tab_list.index(self))

        # Destroy self
        self.destroy()

    def duplicate(self, event=None):
        '''
        Adds a new graph_frame to self.parent.graph_frame_list that is a copy of
        this frame, and changes focus to that frame
        '''
        # Don't run any handlers for this graph_frame
        self.handler_block = True
        self.parent.graph_frame_list.append(new_graph_frame(self.parent))
        new_frame = self.parent.graph_frame_list[-1]
        ix = len(self.parent.graph_frame_list) - 1
        # Copy properties into new frame
        #new_frame.name = self.name + '_copy'
        #new_frame.n_curve = self.n_curve
        # Copy graph widgets
        for i in range(len(self.graph_wids)):
            w = self.graph_wids[i]
            if i == self.ixd['name']:
                new_frame.graph_wids[i].tk_var.set(w.tk_var.get() + '_copy')
            else:
                new_frame.graph_wids[i].copy(w)
        # Run all input validation handlers
        self.parent.notebook.insert(ix, new_frame)
        self.parent.notebook.select(ix)
        self.parent.tab_handler()
        self.update_idletasks()
        new_frame.name_handler()
        new_frame.n_curve_handler()
        # Copy curve widgets
        for i in range(1, self.n_curve + 1):
            new_frame.curve_ix.set(i)
            new_frame.fill_curve_frame() # create the widgets
            for j in range(len(self.curve_dict[i])):
                w = self.curve_dict[i][j]
                new_frame.curve_dict[i][j].copy(w)

    def clone(self, graph, event=None):
        '''
        Clone an existing graph (already defined in tao, use self.duplicate() to
        make copies of graphs already defined in the new plot template window
        Does not clone plot properties
        '''
        #TODO
        # Grab the graph info
        plot_graph = self.pipe.cmd_in("python plot_graph " + graph)
        graph_dict = tao_parameter_dict(plot_graph.splitlines())
        # Copy in the graph info
        for key in graph_dict.keys():
            gkey = key.split('^')[-1] #key as it appears in self.ixd
            if gkey in self.ixd.keys():
                self.graph_wids[self.ixd[gkey]].param_copy(graph_dict[key])
        self.name_handler()
        # Copy the curves
        if ('num_curves' in graph_dict.keys()) and (graph_dict['num_curves'] != None):
            num_curves = graph_dict['num_curves'].value
            self.graph_wids[self.ixd['n_curve']].tk_var.set(num_curves)
            self.n_curve_handler()
            for i in range(1, num_curves+1):
                self.curve_ix.set(i)
                self.fill_curve_frame() # create the widgets if necessary
                curve = graph_dict['curve[' + str(i) + ']'].value
                plot_curve = self.pipe.cmd_in("python plot_curve " + graph + '.' + curve)
                curve_dict = tao_parameter_dict(plot_curve.splitlines())
                for key in curve_dict.keys():
                    ckey = key.split('^')[-1] #key as it appears in self.curve_ixd
                    if ckey in self.curve_ixd.keys():
                        self.curve_dict[i][self.curve_ixd[ckey]].param_copy(curve_dict[key])

    def graph_type_handler(self, *args):
        '''
        Checks that plot.x_axis_type is compatible with the selected graph type.
        If it is, updates self.type and calls self.refresh()
        If it is not, warns the user of the incompatibility.
        '''
        if self.handler_block:
            return
        self.handler_block = True
        new_type = self.head_wids[1].tk_var.get()
        cdict = self.plot.compat_dict
        if new_type in cdict[self.plot.x_axis_type]:
            self.type = new_type
            self.refresh()
        else:
            title = "X-axis/graph type mismatch"
            msg = "The X-axis type you have selected (" + self.plot.x_axis_type
            msg += ") is not compatible with the graph type you have selected for "
            msg += self.name + " (" + new_type + ").\n"
            msg += "The " + self.plot.x_axis_type + " X-axis type is compatible with "
            msg += "the following graph types:"
            for t in cdict[self.plot.x_axis_type]:
                msg += "\n" + t
            messagebox.showwarning(title, msg, parent=self.plot)
        self.handler_block = False

    def name_handler(self, event=None):
        '''
        Checks that a good name has been input for the graph and
        updates self.name as necessary.  Warns the user if the name
        is taken or empty
        '''
        if self.handler_block:
            return
        new_name = self.head_wids[0].tk_var.get().strip()
        self.name_warning_1.grid_forget()
        self.name_warning_2.grid_forget()
        self.name_warning_3.grid_forget()
        # Check what names are taken for this plot
        taken_names = []
        for graph in self.parent.tab_list:
            taken_names.append(graph.name)
        if new_name == "":
            self.name_warning_1.grid(row=2, column=2, sticky='W')
            self.name = "New_graph"
        elif new_name.find(' ') != -1:
            self.name_warning_3.grid(row=2, column=2, sticky='W')
        elif (taken_names.count(new_name) > 1) or (
                (taken_names.count(new_name) == 1) and (new_name != self.name)):
            self.name_warning_2.grid(row=2, column=2, sticky='W')
        else:
            self.name = new_name
        # Update the tab text for this graph
        self.parent.update_name(self.parent.tab_list.index(self))

class new_curve_frame(tk.Frame):
    '''
    Provides a frame to configure curve properties as necessary
    '''
    def __init__(self, parent, graph):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.graph = graph
        self.pipe = graph.pipe
        self._scroll_frame = tao_scroll_frame(self)
        self.name = "New_curve"
        self.grid_columnconfigure(1, weight=1)

        # Delete button
        self.delete_b = tk.Button(
                self, text="DELETE THIS CURVE", fg='red', command=self.delete)

        # Duplicate button
        self.dup_b = tk.Button(
                self, text="Duplicate this curve", command=self.duplicate)

        self.delete_b.grid(row=0, column=0, columnspan=3, sticky='EW')
        self.dup_b.grid(row=1, column=0, columnspan=3, sticky='EW')
        self.grid_columnconfigure(2, weight=1)

        # Curve configuration widgets
        self.curve_frame = self._scroll_frame.frame
        self.curve_frame.grid_columnconfigure(1, weight=1)
        # Helper functions
        def curve_ttp(x):
            '''Returns a tk_tao_parameter for the string x'''
            p = str_to_tao_param(x)
            if x in self.head_wids:
                return tk_tao_parameter(p, self, self.pipe)
            else:
                return tk_tao_parameter(p, self.curve_frame, self.pipe)
        def curve_label(x):
            '''Returns a tk.Label for the string x'''
            if x in self.head_labels:
                return tk.Label(self, text=x)
            else:
                return tk.Label(self.curve_frame, text=x)
        # Widget definitions
        # General
        self.head_wids = ["name;STR;T;"]
        self.head_labels = ["Name:"]
        self.head_wids = list(map(curve_ttp, self.head_wids))
        self.head_labels = list(map(curve_label, self.head_labels))
        # X-axis/graph type specific settings
        # Dictionary of dictionary of lists (self.wids[key1][key2][ix])
        # key1: x_axis_type; key2: graph_type
        self.wids = {}
        self.labels = {}
        for key in ['index', 'ele_index', 'lat', 'var', 's', 'floor',
                'phase_space', 'histogram', 'data', 'none']:
            self.wids[key] = {}
            self.labels[key] = {}
        # Normal data plot
        self.wids['index']['data'] = [
                    'data_source;ENUM;T;',
                    'data_type_x;DAT_TYPE_Z;T;',
                    'data_type;DAT_TYPE;T;',
                    'component;COMPONENT;T;model']
        self.labels['index']['data'] = [ "Data source:",
                "X data type:", "Y data type:",
                "Component:", "Data index:"]
        # Data slice
        self.wids['data']['data'] = [
                'data_source;ENUM;F;data',
                #'data_type_x;DAT_TYPE_E;T;',
                #'data_type;DAT_TYPE_E;T;',
                #'data_index;DAT_TYPE_E;T;',
                'component;COMPONENT;T;model',
                'ele_ref_name;STR;T;']
        self.labels['data']['data'] = [ "Data source:",
                "X data type:", "Y data type:",
                "Data index:", "Component:",
                "Reference Element:"]
        # Data vs lat/var
        self.wids['lat']['data'] = [
                'data_source;ENUM;F;data',
                'data_type_x;STR;T;',
                'data_type;STR;T;',
                'component;COMPONENT;T;model',
                'ele_ref_name;STR;T;']
        self.labels['lat']['data'] = [ "Data source:",
                "X data type:", "Y data type:",
                "Data index:", "Component:",
                "Reference Element:"]
        # Histograms
        self.wids['histogram']['histogram'] = [
                'data_source;ENUM;T;',
                'data_type;ENUM_Z;T;',
                'ele_ref_name;STR;T;',
                'hist;STRUCT;T;density_normalized;LOGIC;T;weight_by_charge;LOGIC;T;number;INT;100;width;INT;']
        self.labels['histogram']['histogram'] = [
                "Data source:", "Data type:", "Reference element:", "Bin settings:"]
        # Phase space plots
        self.wids['phase_space']['phase_space'] = [
                'data_source;ENUM;T;',
                'data_type_x;ENUM_Z;T;',
                'data_type;ENUM_Z;T;',
                'ele_ref_name;STR;T;',
                'use_z_color;LOGIC;T;F',
                'data_type_z;ENUM_Z;T;',
                'z_color0;REAL;T;0',
                'z_color1;REAL;T;0',
                'autoscale_z_color;LOGIC;T;T']
        self.labels['phase_space']['phase_space'] = [
                "Data source:", "X data type:", "Y data type:", "Reference element:",
                "Use z color:", "Color data type:", "Color min:", "Color max:",
                "Autoscale color:"]
        # Map strings in self.wids and self.labels to tk widgets
        for key1 in self.wids.keys():
            for key2 in self.wids[key1].keys():
                self.wids[key1][key2] = list(map(curve_ttp, self.wids[key1][key2]))
                self.labels[key1][key2] = list(map(curve_label, self.labels[key1][key2]))
        # Copy lists as necessary
        self.wids['ele_index']['data'] = self.wids['s']['data'] = self.wids['index']['data']
        self.labels['ele_index']['data'] = self.labels['s']['data'] = self.labels['index']['data']
        self.wids['var']['data'] = self.wids['lat']['data']
        self.labels['var']['data'] = self.labels['lat']['data']
        self.wids['phase_space']['histogram'] = self.wids['histogram']['histogram']
        self.labels['phase_space']['histogram'] = self.labels['histogram']['histogram']
        # Set up empty lists for other axis/graph type compbos to prevent key errors
        x_types = ['index', 'ele_index', 's', 'data', 'lat', 'var', 'floor', 'phase_space', 'histogram', 'none']
        g_types = ['data', 'lay_layout', 'floor_plan', 'dynamic_aperture', 'histogram', 'phase_space']
        for x in x_types:
            if x not in self.wids.keys():
                self.wids[x] = {}
                self.labels[x] = {}
            for g in g_types:
                if g not in self.wids[x].keys():
                    self.wids[x][g] = []
                    self.labels[x][g] = []
        # Style
        self.style_wids = [
                'units;STR;T;',
                'legend_text;STR;T;',
                'y_axis_scale_factor;REAL;T;',
                'use_y2;LOGIC;T;F',
                'draw_line;LOGIC;T;T',
                'draw_symbols;LOGIC;T;T;',
                'draw_symbol_index;LOGIC;T;F',
                'line;STRUCT;T;width;INT;1;color;ENUM;;pattern;ENUM;',
                'symbol;STRUCT;T;type;ENUM;;height;REAL;6.0;color;ENUM;;fill_pattern;ENUM;;line_width;INT;1',
                'symbol_every;INT;T;',
                'smooth_line_calc;LOGIC;T;']
        self.style_labels = [ "Units:", "Legend text:",
                "Y axis scale factor:",
                "Use Y2 axis:",
                "Draw line:",
                "Draw symbols:",
                "Draw symbol index:",
                "Line:",
                "Symbols:",
                "Symbol frequency:",
                "Smooth line calc:"]
        self.style_wids = list(map(curve_ttp, self.style_wids))
        self.style_labels = list(map(curve_label, self.style_labels))

        # Grid head widgets
        for i in range(len(self.head_wids)):
            self.head_labels[i].grid(row=i+2, column=0, sticky='W')
            self.head_wids[i].tk_wid.grid(row=i+2, column=1, sticky='EW')

        # Warning labels
        self.name_warning_1 = tk.Label(self, text="Must not be empty")
        self.name_warning_2 = tk.Label(self, text="Curve name already in use")
        self.name_warning_3 = tk.Label(self, text="Cannot contain whitespace")

        # Responses to edits
        self.head_wids[0].tk_wid.bind('<FocusOut>', self.name_handler)

        # Grid everything else
        self._scroll_frame.grid(row=2+len(self.head_wids), column=0, columnspan=3, sticky='NSEW')
        self.grid_rowconfigure(2+len(self.head_wids), weight=1)
        self.refresh()
        self.update_idletasks()
        self.refresh() #called again to set label widths properly

    def refresh(self):
        '''
        Draws appropriate widgets to the frame depending on graph.type and
        plot.x_axis_type.  The parent graph is expected to verify that
        the graph and x_axis type are compatible before calling this method
        '''
        # Remove existing widgets
        for child in self.curve_frame.winfo_children():
            child.grid_forget()

        # Grid the x-axis/graph type specific widgets
        graph_type = self.graph.type
        x_axis_type = self.graph.plot.x_axis_type
        i=0 # In case self.wids[x_axis_type][graph_type] is empty
        for i in range(len(self.wids[x_axis_type][graph_type])):
            self.labels[x_axis_type][graph_type][i].grid(row=i, column=0, sticky='W')
            self.wids[x_axis_type][graph_type][i].tk_wid.grid(row=i, column=1, sticky='EW')
        offset = i+1

        # Grid the style widgets
        for i in range(len(self.style_wids)):
            ix = i + offset
            self.style_labels[i].grid(row=ix, column=0, sticky='W')
            self.style_wids[i].tk_wid.grid(row=ix, column=1, sticky='EW')

        # Set head label widths properly
        label_width=0
        for child in self.curve_frame.grid_slaves():
            if child.grid_info()['column'] == 0:
                label_width = max(label_width, child.winfo_width())
        self.grid_columnconfigure(0, minsize=label_width)

    def delete(self):
        '''
        Deletes this graph frame
        Call with ask = False to skip confirmation
        '''
        # Ask for confirmation
        if ask:
            ans = messagebox.askokcancel("Delete " + self.name, "Delete this curve?", parent=self.graph.plot)
            if not ans:
                return

        # Remove from tabbed frame
        self.parent.remove_tab(self.parent.tab_list.index(self))

        # Destroy self
        self.destroy()

    def duplicate(self):
        pass

    def name_handler(self, event=None):
        '''
        Checks that a good name has been input for the curve and
        updates self.name as necessary.  Warns the user if the name
        is taken or empty
        '''
        if self.handler_block:
            return
        new_name = self.head_wids[0].tk_var.get().strip()
        self.name_warning_1.grid_forget()
        self.name_warning_2.grid_forget()
        self.name_warning_3.grid_forget()
        # Check what names are taken for this plot
        taken_names = []
        for curve in self.parent.tab_list:
            taken_names.append(curve.name)
        if new_name == "":
            self.name_warning_1.grid(row=2, column=2, sticky='W')
            self.name = "New_curve"
        elif new_name.find(' ') != -1:
            self.name_warning_3.grid(row=2, column=2, sticky='W')
        elif (taken_names.count(new_name) > 1) or (
                (taken_names.count(new_name) == 1) and (new_name != self.name)):
            self.name_warning_2.grid(row=2, column=2, sticky='W')
        else:
            self.name = new_name
        # Update the tab text for this graph
        self.parent.update_name(self.parent.tab_list.index(self))

