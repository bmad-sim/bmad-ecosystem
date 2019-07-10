from taoplot import taoplot
import matplotlib
import matplotlib.pyplot as plt

def in_element(x,y,taoplot):
	'''takes coordinates and taoplot instance, returns list of element indices of elements that contain the specified coordinates'''
	returnList = taoplot.plot()[1]
	inIndexList = []
	if returnList[0] == 'lat_layout' or returnList[0] == 'data':
		for i in returnList[1]:
			if returnList[2][str(i)] < returnList[3][str(i)]:
				if returnList[2][str(i)] < x < returnList[3][str(i)]:
					inIndexList.append(i)
			if returnList[2][str(i)] > returnList[3][str(i)]:
				if x > returnList[2][str(i)] or x < returnList[3][str(i)]:
					inIndexList.append(i)

	return inIndexList
