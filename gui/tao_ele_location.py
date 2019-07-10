from taoplot import taoplot
import matplotlib
import matplotlib.pyplot as plt

def in_element(x,y,taoplot):
	'''takes coordinates and taoplot instance, returns list of element indices of elements that contain the specified coordinates'''
	inIndexList = []
	returnList = taoplot.plot()[1]
	try:
		if returnList[0] == 'lat_layout' or returnList[0] == 'data':
			for i in returnList[4]:
				if returnList[5][str(i)] < returnList[6][str(i)]:
					if returnList[5][str(i)] < x < returnList[6][str(i)]:
						inIndexList.append(i)
				if returnList[5][str(i)] > returnList[6][str(i)]:
					if x > returnList[5][str(i)] or x < returnList[6][str(i)]:
						inIndexList.append(i)

	except TypeError:
		pass

	return inIndexList
