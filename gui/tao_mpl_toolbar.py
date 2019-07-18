import numpy as np
import tkinter as tk
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

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


