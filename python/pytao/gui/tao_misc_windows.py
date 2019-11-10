'''
Provides miscellaneous windows needed in the GUI for tao
'''
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

from .tao_widget import *
from pytao.util.parameters import str_to_tao_param, tao_parameter_dict
from .tao_set import *
from .tao_base_windows import *

#----------------------------------------------------
# Global Variables window

class tao_global_vars_window(tao_parameter_window):
    def __init__(self, root):
        global_list = root.pipe.cmd_in("python global")
        global_list = global_list.splitlines()
        for i in range(len(global_list)):
            global_list[i]=str_to_tao_param(global_list[i])
        tao_parameter_window.__init__(self, root, "Global Variables", global_list, root.pipe)
        b = tk.Button(self.button_frame, text="Set Global Variables",
                command=self.set_callback)
        b.pack()

    def set_callback(self):
        tao_set(self.tao_list, "set global ", self.root.pipe)
        # Refresh windows
        if self.root.pipe.cmd_in('python lat_calc_done') == 'T':
            for k in self.root.refresh_windows.keys():
                for win in self.root.refresh_windows[k]:
                    win.refresh()

#----------------------------------------------------
# bmad_com window

class tao_bmad_com_window(tao_parameter_window):
    def __init__(self, root):
        bmad_list = root.pipe.cmd_in("python bmad_com")
        bmad_list = bmad_list.splitlines()
        for i in range(len(bmad_list)):
            bmad_list[i]=str_to_tao_param(bmad_list[i])
        tao_parameter_window.__init__(self, root, "Bmad Parameters", bmad_list, root.pipe)
        b = tk.Button(self.button_frame, text="Set Bmad Parameters",
                command=self.set_callback)
        b.pack()

    def set_callback(self):
        tao_set(self.tao_list, "set bmad_com ", self.root.pipe)
        # Refresh windows
        if self.root.pipe.cmd_in('python lat_calc_done') == 'T':
            for k in self.root.refresh_windows.keys():
                for win in self.root.refresh_windows[k]:
                    win.refresh()

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



