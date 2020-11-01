import tkinter as tk
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_tools import *
from matplotlib import rcParams
from matplotlib.widgets import Slider
from matplotlib.backend_bases import ToolContainerBase
from matplotlib.backends._backend_tk import FigureManagerTk, ToolbarTk
from matplotlib.backend_managers import ToolManager

from pytao.util.parameters import *

#---------------------------------------------
# Note: Help and recalculate buttons are commented out since later version of matplotlib/tkinter will not
# work with absolute path names. It would be good to fix this at some point...

class taotoolbar(NavigationToolbar2Tk):
    def __init__(self,canvas_,parent_,width_,GUI_DIR_):
        self.basedir = os.path.join(rcParams['datapath'], 'images')
        self.GUI_DIR = GUI_DIR_
        self.toolitems = (
            ('Home', 'Reset original view', 'home', 'home'),
            ('Back', 'Back to previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            ('Pan', 'Pan axes with left mouse, zoom with right mouse or scroll wheel', 'move', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
##            ('Redraw','Recalculate points',self.GUI_DIR+'/assets/recalculate','redraw'),
            ('Slider','Width slider','subplots','slider'),
            (None,None,None,None),
            ('Save', 'Save image of figure', 'filesave', 'save_figure'),
            (None,None,None,None),
##            ('Help','Graph help',self.GUI_DIR+'/assets/help','help'),
            )
        self.parent = parent_
        self.canvas = canvas_
        self.width = width_ #width used with slider tool
        self.axes_list = self.canvas.figure.get_axes()
        self.templateDict = self.parent.root.placed
        self.graph_list = self.parent.fig_info[18]
        self.template = self.templateDict[self.parent.region] #corresponding template for a region in tao
        self.templateGraphList = []
        for i in self.graph_list:
            self.templateGraphList.append(self.template + '.' + i.split('.')[1])
        self.gRangeList = [] #list of original x min, x max, y min, and y max for each graph
        try:
            for i in self.templateGraphList:
                gInfo=self.parent.pipe.cmd_in('python plot_graph '+i,no_warn = True).splitlines()
                gRange = []
                gInfoDict = {}
                for i in range(len(gInfo)):
                    gInfoDict[gInfo[i].split(';')[0]]=str_to_tao_param(gInfo[i])
                gRange.append(gInfoDict['x'].get_component('min'))
                gRange.append(gInfoDict['x'].get_component('max'))
                gRange.append(gInfoDict['y'].get_component('min'))
                gRange.append(gInfoDict['y'].get_component('max'))
                self.gRangeList.append(gRange)
        except IndexError: #incase number of graphs changes, like in wave analysis
            pass
        NavigationToolbar2Tk.__init__(self,canvas_,parent_)

    cid = 'none' #connection id for scroll wheel
    cidKeyP = 'none' #connection id for key press
    cidKeyR = 'none' #connection id for key release

    cur_ax = 'none' #axes instance that the mouse is currently over or 'none'
    xzoom = True #allow for x axis scaling
    yzoom = True #allow for y axis scaling



    def enter_axes(self,event):
        '''handles the mouse entering axes'''
        if event.inaxes is not None:
            self.cur_ax = event.inaxes


    def leave_axes(self,event):
        '''handles the mouse leaving axes'''
        if event.inaxes is None:
            self.cur_ax = 'none'

    def onKeyPress(self,event):
        '''handles key presses'''
        if event.key == 'x':
            self.xzoom = True
            self.yzoom = False
        if event.key == 'y':
            self.xzoom = False
            self.yzoom = True

    def onKeyRelease(self,event):
        '''handles releasing key presses'''
        self.xzoom = True
        self.yzoom = True



    def zoom_factory(self,axes,canv,event,on = False,cid='none',base_scale = 1.5):
        '''controls connections for scroll wheel zooming'''

        fig = canv[0].get_figure() # get the figure of interest

        enter = fig.canvas.mpl_connect('axes_enter_event', self.enter_axes)
        leave = fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)

        self.cidKeyP = fig.canvas.mpl_connect('key_press_event',self.onKeyPress)
        self.cidKeyR = fig.canvas.mpl_connect('key_release_event',self.onKeyRelease)

        def zoom_fun(event):
            '''changes graph axes if scroll wheel is used'''
            if self.cur_ax != 'none':
                try:
                    ax = self.cur_ax

                    cur_xlim = ax.get_xlim()
                    cur_ylim = ax.get_ylim()
                    #get the current x and y limits

                    xdata = event.xdata #get event x location
                    ydata = event.ydata #get event y location
                    x_left = xdata - cur_xlim[0]
                    x_right = cur_xlim[1] - xdata
                    y_top = ydata - cur_ylim[0]
                    y_bottom = cur_ylim[1] - ydata
                    if event.button == 'up': # dealwith zoom in
                        scale_factor = 1/base_scale
                    elif event.button == 'down': #deal with zoom out
                        scale_factor = base_scale
                    else:
                        scale_factor = 1
                    if self.xzoom == True: #new x limits
                        ax.set_xlim([xdata - x_left*scale_factor,xdata + x_right*scale_factor])
                    if self.yzoom == True: #new y limits
                        ax.set_ylim([ydata - y_top*scale_factor,ydata + y_bottom*scale_factor])
                    self.canvas.draw_idle()
                except TypeError:
                    pass #scrolling outside of graph
        # attach the call backs
        if on == True and cid == 'none':
            self.cid=fig.canvas.mpl_connect('scroll_event',zoom_fun)

        if on != True:
            fig.canvas.mpl_disconnect(cid)
            self.cid='none'

        if self._nav_stack() is None: #make home button work with scroll wheel zoom
            self.push_current()

        #return the function
        return zoom_fun



    def home(self, *args):
        """Restore the original view."""
        try:
            for i in range(len(self.graph_list)):
                self.parent.pipe.cmd_in('set graph '+self.graph_list[i]+' x.min = '+str(self.gRangeList[i][0]),no_warn = True)
                self.parent.pipe.cmd_in('set graph '+self.graph_list[i]+' x.max = '+str(self.gRangeList[i][1]),no_warn = True)
                self.parent.pipe.cmd_in('set graph '+self.graph_list[i]+' y.min = '+str(self.gRangeList[i][2]),no_warn = True)
                self.parent.pipe.cmd_in('set graph '+self.graph_list[i]+' y.max = '+str(self.gRangeList[i][3]),no_warn = True)
        except IndexError:
            pass
        self.parent.refresh()
        self._nav_stack.home()
        self.set_history_buttons()
        self._update_view()



    def pan(self, *args):
        """Activate the pan/zoom tool. pan with left button, zoom with right or with with scroll wheel"""
        #set the pointer icon and button press funcs to the appropriate callbacks
        #default matplotlib pan tool, but with added scroll wheel zooming in pan/zoom mode
        if self._active == 'PAN':
            self._active = None
        else:
            self._active = 'PAN'
        if self._idPress is not None:
            self._idPress = self.canvas.mpl_disconnect(self._idPress)
            self.mode = ''

        if self._idRelease is not None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)
            self.mode = ''

        if self._active:

            self._idPress = self.canvas.mpl_connect(
                'button_press_event', self.press_pan)
            self._idRelease = self.canvas.mpl_connect(
                'button_release_event', self.release_pan)
            self.mode = 'pan/zoom'
            scale_canv = self.canvas
            scale = self.zoom_factory(self,scale_canv.figure.get_axes(),scale_canv,True,self.cid)
            self.canvas.widgetlock(self)
        else:
            self.canvas.widgetlock.release(self)
            scale_canv = self.canvas
            scale = self.zoom_factory(self,scale_canv.figure.get_axes(),scale_canv,False,self.cid)
        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self._active)

        self.set_message(self.mode)



    def help(self):
        '''help tool which opens a tk window with instructions on graph interaction'''
        image = None

        win = tk.Toplevel(self.parent)

        # Title string
        title = 'Graph Help'

        #Help items and descriptions
        help_items = ['Home',  'Back', 'Forward', 'Pan/Zoom', 'Zoom Rectangle', 'Recalculate Points', 'Width Slider', 'Save', '<Double Click>']
        help_descrips = [
            "Returns to original view with original points, shortcuts are 'h' or 'r'.",
            "Returns to previous view, shortcuts are 'c' or 'left arrow'.",
            "Undoes the last back command, shortcuts are 'v' or 'right arrow'.",
            """Toggles panning by left clicking and dragging, and zooming by using the scroll wheel or by right clicking and dragging. 
Holding 'x' restricts panning and zooming to the x axis, holding 'y' restricts panning and zooming to the y axis, and holding control maintains the aspect ratio.""",
            """Toggles zooming by left clicking and dragging to select the new window. 
Holding 'x' restricts zooming to the x axis, holding 'y' restricts zooming to the y axis, and holding control maintains the aspect ratio.""",
            "Recalculates points to better fit the current window size",
            "Makes a slider that changes the size of a lat layout if one is drawn below a graph, or changes width of elements in a floor plan.",
            "Saves the current figure as an image file, shortcut is 'ctrl + s'.",
            "With Lat_Layout or Floor_Plan: Double click on Shape to open element parameter window."]

        wl = 600
        tk.Label(win, text=title).grid(row=0, column=0, columnspan=2, sticky='EW')
        for i in range(len(help_items)):
            tk.Label(win, text=help_items[i]).grid(row=2*i+1, column=0, sticky='NW', padx = 10)
            tk.Label(win, text='').grid(row=2*i+2, column=0, sticky='W')
            tk.Label(win, text=help_descrips[i], wraplength=wl, justify='left').grid(row=2*i+1, column=1, sticky='W')
            tk.Label(win, text='', wraplength=wl, justify='left').grid(row=2*i+2, column=1, sticky='W')
        win.grid_columnconfigure(1, weight=1)



    def redraw(self):
        '''forces recalculation of points in window'''

        for i in range(len(self.graph_list)):
            xRange = self.axes_list[i].get_xlim()
            yRange = self.axes_list[i].get_ylim()
            self.parent.pipe.cmd_in('set graph '+self.graph_list[i]+' x.min = '+str(xRange[0]),no_warn = True)
            self.parent.pipe.cmd_in('set graph '+self.graph_list[i]+' x.max = '+str(xRange[1]),no_warn = True)
            self.parent.pipe.cmd_in('set graph '+self.graph_list[i]+' y.min = '+str(yRange[0]),no_warn = True)
            self.parent.pipe.cmd_in('set graph '+self.graph_list[i]+' y.max = '+str(yRange[1]),no_warn = True)

        self.parent.refresh()



    def slider(self):
        '''opens a new matplotlib figure with a slider that adjusts the size of a lat layout if one is present or adjusts the width of floor plan elements if a floor plan is present'''
        rcParams['toolbar'] = 'None' #hides toolbar on figure with slider
        slidefig = plt.figure(figsize=[4,.7])
        slideax = slidefig.add_subplot(1,1,1)
        plt.title('Element Width Slider')
        width_slider = Slider(slideax, 'width', 0, 2, self.width) #element width slider
        slidefig.tight_layout(pad=.5)
        slidefig.show()
        def update_slider(width):
            for i in range(len(self.graph_list)):
                xRange = self.axes_list[i].get_xlim()
                yRange = self.axes_list[i].get_ylim()
            self.parent.refresh(width=width_slider.val)
            for i in range(len(self.graph_list)):
                self.axes_list[i].set_xlim(xRange)
                self.axes_list[i].set_ylim(yRange)
        width_slider.on_changed(update_slider) #call update when slider moves



