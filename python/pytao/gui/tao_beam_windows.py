'''
Provides windows for viewing and editing beam properties in tao
'''
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

from pytao.util.parameters import str_to_tao_param, tao_parameter_dict
from .tao_widget import *
from .taoplot import taoplot
from .tao_set import *
from .tao_base_windows import *
from .tao_lat_windows import tao_ele_window
from .tao_ele_location import in_element
from .tao_mpl_toolbar import taotoolbar

#-------------------------------------------------------
# Beam_init_window
class tao_beam_init_window(tao_parameter_window):
    '''
    Window for viewing and editing the beam_init struct
    '''
    def __init__(self, root, pipe):
        self.tao_id = "beam"
        self.pipe = pipe
        tao_list_window.__init__(self, root, "Beam Settings")
        self.button_frame = tk.Frame(self)
        self.button_frame.pack(side="bottom", fill="both", expand=0)

        self.apply_b = tk.Button(self.button_frame, text="Apply",
                command=self.apply)
        self.apply_b.pack()

        self.beam_head_list = []
        self.beam_body_list = []
        self.refresh()

    def refresh(self, event=None):
        '''
        Pulls beam settings from tao and updates the window contents
        '''
        for child in self.list_frame.winfo_children():
            child.destroy()

        beam_head_list = self.pipe.cmd_in("python beam 1")
        beam_body_list = self.pipe.cmd_in("python beam_init 1")
        #TODO: commands above should not be hardcoded to universe 1
        beam_head_list = beam_head_list.splitlines()
        beam_body_list = beam_body_list.splitlines()
        for k in range(len(beam_head_list)):
            beam_head_list[k] = tk_tao_parameter(str_to_tao_param(beam_head_list[k]),
                    self.list_frame, self.pipe)
            beam_head_list[k].tk_label.grid(row=k, column=0, sticky='W')
            beam_head_list[k].tk_wid.grid(row=k, column=1, sticky='EW')
        b_offset = len(beam_head_list)
        for k in range(len(beam_body_list)):
            beam_body_list[k] = tk_tao_parameter(str_to_tao_param(beam_body_list[k]),
                    self.list_frame, self.pipe)
            beam_body_list[k].tk_label.grid(row=k+b_offset, column=0, sticky='W')
            beam_body_list[k].tk_wid.grid(row=k+b_offset, column=1, sticky='EW')
        self.beam_head_list = beam_head_list
        self.beam_body_list = beam_body_list

    def apply(self, event=None):
        '''
        Runs set commands to apply changes from this window to tao
        '''
        tao_set(self.beam_head_list, "set beam ", self.pipe)
        tao_set(self.beam_body_list, "set beam_init ", self.pipe)
        self.refresh()

