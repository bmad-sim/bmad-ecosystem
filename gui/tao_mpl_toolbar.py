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
			('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
			('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
			(None,None,None,None),
			('Save', 'Save image of figure', 'filesave', 'save_figure'),
			)
		NavigationToolbar2Tk.__init__(self,canvas_,parent_)





	def zoom_factory(self,axes,canv,event,on = False,cid='none',base_scale = 2):
		ax = canv[0].get_figure().get_axes()[0]
		def zoom_fun(event):
			# get the current x and y limits
			cur_xlim = ax.get_xlim()
			cur_ylim = ax.get_ylim()
			cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
			cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
			xdata = event.xdata # get event x location
			ydata = event.ydata # get event y location
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
			ax.set_xlim([xdata - cur_xrange*scale_factor,
				xdata + cur_xrange*scale_factor])
			ax.set_ylim([ydata - cur_yrange*scale_factor,
				ydata + cur_yrange*scale_factor])
			self.draw()

		fig = canv[0].get_figure() # get the figure of interest

		# attach the call back
		if on == True and cid == 'none':
			cid=fig.canvas.mpl_connect('scroll_event',zoom_fun)

		if on != True:
			fig.canvas.mpl_disconnect(cid)
			cid='none'

		#return the function
		return zoom_fun,cid

	

	def pan(self, *args):
		"""Activate the pan/zoom tool. pan with left button, zoom with right"""
		# set the pointer icon and button press funcs to the
		# appropriate callbacks
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
			scale = self.zoom_factory(self,scale_canv.figure.get_axes(),scale_canv,True)
			cid=scale[1]
			self.canvas.widgetlock(self)
		else:
			self.canvas.widgetlock.release(self)
			scale_canv = self.canvas
			scale = self.zoom_factory(self,scale_canv.figure.get_axes(),scale_canv,False,cid)
			cid=scale[1]
		for a in self.canvas.figure.get_axes():
			a.set_navigate_mode(self._active)

		self.set_message(self.mode)

		





