from taoplot import taoplot
import matplotlib
import matplotlib.pyplot as plt

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
		#find lat_layout elements containing the point
	
		if returnList[0] == 'floor_plan':
			for i in returnList[7]:
				if ((x-returnList[8][str(i)][0])**2 + (y-returnList[8][str(i)][1])**2) <= ((returnList[9][i])**2):
					inIndexList.append(i)
	except TypeError:
		pass
	except IndexError:
		pass

	return inIndexList
