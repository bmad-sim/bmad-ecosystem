import numpy as np
import tkinter as tk
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_tools import *
from matplotlib import rcParams

class taotoolbar(NavigationToolbar2Tk):
	def __init__(self,canvas_,parent_):
		self.toolitems = (
			('Home', 'Reset original view', 'home', 'home'),
			('Back', 'Back to previous view', 'back', 'back'),
			('Forward', 'Forward to next view', 'forward', 'forward'), 
			('Pan', 'Pan axes with left mouse, zoom with right mouse or scroll wheel', 'move', 'pan'),
			('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
			(None,None,None,None),
			('Save', 'Save image of figure', 'filesave', 'save_figure'),
			)
		NavigationToolbar2Tk.__init__(self,canvas_,parent_)


	cid = 'none' #connection id for scroll wheel
	cur_ax = 'none' #axes instance that the mouse is currently over or 'none'
	cidKeyP = 'none'
	cidKeyR = 'none'
	xzoom = True
	yzoom = True

	def enter_axes(self,event):
		if event.inaxes is not None:
			self.cur_ax = event.inaxes


	def leave_axes(self,event):
		if event.inaxes is None:
			self.cur_ax = 'none'

	def onKeyPress(self,event):
		if event.key == 'x':
			self.xzoom = True
			self.yzoom = False
		if event.key == 'y':
			self.xzoom = False
			self.yzoom = True

	def onKeyRelease(self,event):
		self.xzoom = True
		self.yzoom = True



	def zoom_factory(self,axes,canv,event,on = False,cid='none',base_scale = 1.5):
		'''controls connections for scroll wheel zooming'''

		fig = canv[0].get_figure() # get the figure of interest

		enter=fig.canvas.mpl_connect('axes_enter_event', self.enter_axes)
		leave=fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)

		self.cidKeyP = fig.canvas.mpl_connect('key_press_event',self.onKeyPress)
		self.cidKeyR = fig.canvas.mpl_connect('key_release_event',self.onKeyRelease)

		def zoom_fun(event):
			'''changes graph axes if scroll wheel is used'''
			if self.cur_ax != 'none':
				ax = self.cur_ax
				# get the current x and y limits
				cur_xlim = ax.get_xlim()
				cur_ylim = ax.get_ylim()
				xdata = event.xdata # get event x location
				ydata = event.ydata # get event y location
				x_left = xdata - cur_xlim[0]
				x_right = cur_xlim[1] - xdata 
				y_top = ydata - cur_ylim[0]
				y_bottom = cur_ylim[1] - ydata
				if event.button == 'up':
					# deal with zoom in
					scale_factor = 1/base_scale
				elif event.button == 'down':
					# deal with zoom out
					scale_factor = base_scale
				else:
					# deal with something that should never happen
					scale_factor = 1
				# set new limits
				if self.xzoom == True:
					ax.set_xlim([xdata - x_left*scale_factor,xdata + x_right*scale_factor])
				if self.yzoom == True:
					ax.set_ylim([ydata - y_top*scale_factor,ydata + y_bottom*scale_factor])
				self.canvas.draw_idle()

		# attach the call back
		if on == True and cid == 'none':
			self.cid=fig.canvas.mpl_connect('scroll_event',zoom_fun)

		if on != True:
			fig.canvas.mpl_disconnect(cid)
			self.cid='none'

		if self._nav_stack() is None:
			# set the home button to this view
			self.push_current()

		#return the function
		return zoom_fun

	

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

		





