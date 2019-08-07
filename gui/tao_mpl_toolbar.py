import numpy as np
import tkinter as tk
import matplotlib as mpl
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
					
			self.canvas.widgetlock(self)
		else:
			self.canvas.widgetlock.release(self)

		for a in self.canvas.figure.get_axes():
			a.set_navigate_mode(self._active)

		self.set_message(self.mode)

		





