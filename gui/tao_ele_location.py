from taoplot import taoplot
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math


def in_rectangle(a,b,c,d,x,y):
	'''checks if the point (x,y) is in a rectangle with corners a, b, c, and d'''
	ABC = .5*abs(a[0]*(b[1]-c[1])+b[0]*(c[1]-a[1])+c[0]*(a[1]-b[1]))
	BCD = .5*abs(b[0]*(c[1]-d[1])+c[0]*(d[1]-b[1])+d[0]*(b[1]-c[1]))
	areaRectangle = ABC + BCD

	ABP = .5*abs(a[0]*(b[1]-y)+b[0]*(y-a[1])+x*(a[1]-b[1]))
	BDP = .5*abs(b[0]*(d[1]-y)+d[0]*(y-b[1])+x*(b[1]-d[1]))
	CDP = .5*abs(c[0]*(d[1]-y)+d[0]*(y-c[1])+x*(c[1]-d[1]))
	CAP = .5*abs(c[0]*(a[1]-y)+a[0]*(y-c[1])+x*(c[1]-a[1]))

	if areaRectangle != 0:
		return math.isclose(areaRectangle,(ABP+BDP+CDP+CAP))
	else:
		return False


def in_element(x,y,fig_info):
	'''takes coordinates and taoplot.plot() output, returns list of element indices of elements that contain the specified coordinates'''
	inIndexList = []
	returnList = fig_info
	try:
		if returnList[0] == 'lat_layout' or returnList[0] == 'data':
			for i in returnList[4]:
				if returnList[5][str(i)] < returnList[6][str(i)]:
					if returnList[5][str(i)] < x < returnList[6][str(i)]:
						inIndexList.append(i)
				if returnList[5][str(i)] > returnList[6][str(i)]:
					if x > returnList[5][str(i)] or x < returnList[6][str(i)]:
						inIndexList.append(i)
		#find lat_layout elements containing the specified point
		#checks x ranges of each element to locate the specified point
	
		if returnList[0] == 'floor_plan':
			for i in returnList[7]:
				if returnList[8][str(i)] == 'circle':
					if ((x-returnList[9][str(i)][0])**2 + (y-returnList[9][str(i)][1])**2) <= ((returnList[10][str(i)])**2):
						inIndexList.append(i)
				
				else:
					if in_rectangle(returnList[11][str(i)],returnList[12][str(i)],returnList[13][str(i)],returnList[14][str(i)],x,y) == True:
						
						inIndexList.append(i)
		#find floor_plan elements containing the specified point
		#checks a circle for the specified point for circle elements, checks a box otherwise
		#checks a box made from the element start and end coordinates with y1 and y2 for the height
				
	except TypeError:
		pass
	except IndexError:
		pass

	return inIndexList
