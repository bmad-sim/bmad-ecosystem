#Kevin Kowalski kjk226 6/17/2019

import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from tao_interface import *
from parameters import *
import time


#pipe=tao_interface('pexpect','-init_file ../examples/cesr/tao.init')
#determines tao settings and lattice to be used eg: '-init_file ../examples/cesr/tao.init' for CESR

#GraphRegion='top'
#region in tao of plot

class taoplot:
	def __init__(self,pipe,GraphRegion):
		self.pipe= pipe
		#tao_interface object
		self.GraphRegion= GraphRegion
		#string describing region in tao of the desired plot

	def plot(self):
		fig = plt.figure()
		#creates plotting figure
		pipe = self.pipe
		GraphRegion = self.GraphRegion

		def pgp_to_mpl(x):
			'''takes string with pgplot characters and returns string with characters replaced with matplotlib equivalent, raises NotImplementedError if an unknown pgplot character is used'''
			x=x.replace('\\','\\\\')
			if '\\\\' in x:
				lx = '$'+x+'$'
				while lx.find('\\\\d') != -1 and lx.find('\\\\u') != -1:
					if lx.find('\\\\d') <	lx.find('\\\\u'):
						if lx.find('\\\\d') != -1:
							start=lx.find('\\\\d')
							end=lx.find('\\\\u')
							sx=lx[start:end+3]
							lx=lx.replace(sx,'_'+sx[3:-3])
							sx=''
							start=''
					else:
						if lx.find('\\\\u') != -1:
							start=lx.find('\\\\u')
							end=lx.find('\\\\d')
							sx=lx[start:end+3]
							lx=lx.replace(sx,'^'+sx[3:-3])
							sx=''
							start=''
							end=''

				lx=lx.replace(' ','\\ ')

				lx=lx.replace('\\\\ga','\\alpha')
				lx=lx.replace('\\\\gb','\\beta')
				lx=lx.replace('\\\\gg','\\gamma')
				lx=lx.replace('\\\\gd','\\delta')
				lx=lx.replace('\\\\ge','\\epsilon')
				lx=lx.replace('\\\\gz','\\zeta')
				lx=lx.replace('\\\\gy','\\eta')
				lx=lx.replace('\\\\gh','\\theta')
				lx=lx.replace('\\\\gi','\\iota')
				lx=lx.replace('\\\\gk','\\kappa')
				lx=lx.replace('\\\\gl','\\lambda')
				lx=lx.replace('\\\\gm','\\mu')
				lx=lx.replace('\\\\gn','\\nu')
				lx=lx.replace('\\\\gc','\\xi')
				lx=lx.replace('\\\\go','\\omicron')
				lx=lx.replace('\\\\gp','\\pi')
				lx=lx.replace('\\\\gr','\\rho')
				lx=lx.replace('\\\\gs','\\sigma')
				lx=lx.replace('\\\\gt','\\tau')
				lx=lx.replace('\\\\gu','\\upsilon')
				lx=lx.replace('\\\\gf','\\phi')
				lx=lx.replace('\\\\gx','\\chi')
				lx=lx.replace('\\\\gq','\\psi')
				lx=lx.replace('\\\\gw','\\omega')

				lx=lx.replace('\\\\gA','A')
				lx=lx.replace('\\\\gB','B')
				lx=lx.replace('\\\\gG','\\Gamma')
				lx=lx.replace('\\\\gD','\\Delta')
				lx=lx.replace('\\\\gE','E')
				lx=lx.replace('\\\\gZ','Z')
				lx=lx.replace('\\\\gY','H')
				lx=lx.replace('\\\\gH','\\Theta')
				lx=lx.replace('\\\\gI','I')
				lx=lx.replace('\\\\gK','\\Kappa')
				lx=lx.replace('\\\\gL','\\Lambda')
				lx=lx.replace('\\\\gM','M')
				lx=lx.replace('\\\\gN','N')
				lx=lx.replace('\\\\gC','\\Xi')
				lx=lx.replace('\\\\gO','O')
				lx=lx.replace('\\\\gP','\\Pi')
				lx=lx.replace('\\\\gR','P')
				lx=lx.replace('\\\\gS','\\Sigma')
				lx=lx.replace('\\\\gT','T')
				lx=lx.replace('\\\\gU','\\Upsilon')
				lx=lx.replace('\\\\gF','\\Phi')
				lx=lx.replace('\\\\gX','X')
				lx=lx.replace('\\\\gQ','\\Psi')
				lx=lx.replace('\\\\gW','\\Omega')
				if '\\\\' in lx:
					raise NotImplementedError('unknown character in string, character not yet added to pgp_to_mpl')
				return lx
			else:
				return x

		SymbolsDict={
			'square':'s', #no fill
			'dot':'.',
			'plus':'+',
			'times':(6,2,0),
			'circle':'$\\circ$', 
			'x':'x',
			'triangle':'^', #no fill
			'circle_plus':'$\\oplus$',
			'circle_dot':'$\\odot$',
			'square_concave':(4,1,45),
			'diamond':'d', #no fill
			'star5':'*', #no fill
			'triangle_filled':'^', #fill
			'red_cross':'P',
			'star_of_david':(6,1,0),
			'square_filled':'s', #fill
			'circle_filled':'o', #fill
			'star5_filled':'*', #fill
			'0':'s', #no fill
			'1':'.',
			'2':'+',
			'3':(6,2,0),
			'4':'$\\circ$', 
			'5':'x',
			'6':'s',
			'7':'^', #no fill
			'8':'$\\oplus$',
			'9':'$\\odot$',
			'10':(4,1,45),
			'11':'d', #no fill
			'12':'*', #no fill
			'13':'^', #fill
			'14':'P',
			'15':(6,1,0),
			'16':'s', #fill
			'17':'o', #fill
			'18':'*', #fill
			'-1':',',
			'-2':',',
			'-3':'(3,0,0)',
			'-4':'(4,0,0)',
			'-5':'(5,0,0)',
			'-6':'(6,0,0)',
			'-7':'(7,0,0)',
			'-8':'(8,0,0)',
			'-9':'(9,0,0)',
			'-10':'(10,0,0)',
			'-11':'(11,0,0)',
			'-12':'(12,0,0)'
		}
		#Dictionary with pgplot symbol strings as keys and corresponding matplotlib symbol strings as values




		'''''''''Graphing'''''''''

		'''Region Data'''

		rInfo=pipe.cmd_in('python plot1 '+GraphRegion).splitlines()
		#list of plotting parameter strings from tao command python plot1

		rInfoList = []
		rInfoDict = {}
		for i in range(len(rInfo)):
			rInfoDict[rInfo[i].split(';')[0]]=str_to_tao_param(rInfo[i])
			rInfoList.append(rInfo[i].split(';')[0])
		#list of tao_parameter object names from python plot_1
		#dictionary of tao_parameter name string keys to the corresponding tao_parameter object

		gList = []
		for i in range(rInfoDict['num_graphs'].value):
			gList.append(rInfoDict[('graph['+str(i+1)+']')].value)
		#list of string names of graphs

		GraphDict = {}
			

		for gNumber in range(len(gList)):

			if gNumber == 0:
				GraphDict['graph'+str(gNumber+1)]=fig.add_subplot(len(gList)+1,1,gNumber+1)
			elif gNumber > 0:
				GraphDict['graph'+str(gNumber+1)]=fig.add_subplot(len(gList)+1,1,gNumber+1,sharex=GraphDict['graph1'])
			#create plots in figure, second line also makes x axes scale together



			'''Graph Data'''

			gType = GraphRegion+'.'+gList[gNumber]
			#graph type, like r13.g or top.x


			gInfo=pipe.cmd_in('python plot_graph '+gType).splitlines()
			#list of plotting parameter strings from tao command python plot_graph


			gInfoList = []
			gInfoDict = {}
			for i in range(len(gInfo)):
				gInfoDict[gInfo[i].split(';')[0]]=str_to_tao_param(gInfo[i])
				gInfoList.append(gInfo[i].split(';')[0])
			#list of tao_parameter object names from python plot_graph
			#dictionary of tao_parameter name string keys to the corresponding tao_parameter object




			'''Curve Data'''

			cList=[]
			for i in range(gInfoDict['num_curves'].value):
				cList.append(gInfoDict[('curve['+str(i+1)+']')].value)
			#list of string names of curves

			cInfo=[]
			for i in cList:
				cInfo.append(pipe.cmd_in('python plot_curve '+gType+'.'+i).splitlines())
			#list of lists of plotting parameter strings from tao command python plot_curve for each curve

			cInfoSuperList=[]
			cInfoDictList=[]
			for i in range(len(cList)):
				cInfoList = []
				cInfoDict = {}
				for j in range(len(cInfo[i])):
					cInfoDict[cInfo[i][j].split(';')[0]]=str_to_tao_param(cInfo[i][j])
					cInfoList.append(cInfo[i][j].split(';')[0])
				cInfoSuperList.append(cInfoList)
				cInfoDictList.append(cInfoDict)
				cInfoList = []
				cInfoDict = {}
			#list of lists of tao_parameter object names from python plot_graph for each curve
			#list of dictionaries of tao_parameter name string keys to the corresponding tao_parameter object for each curve

			#print(cInfoDictList[1]['name'].value) #gives name of 2nd curve which is y for orbit




			'''Line Data'''

			lInfo=[]
			try:
				for i in cList:
					lInfo.append(pipe.cmd_in('python plot_line '+gType+'.'+i).splitlines())
				#list of points from tao command python plot_line for each curve
				if len(lInfo) != 0:
					PointsSuperList=[]
					for i in range(len(cList)):
						pList = []
						for j in range(len(lInfo[i])):
							pList.append(lInfo[i][j].split(';'))
						LineCoords = []
						for j in range(len(pList)):
							LineCoords.append([float(pList[j][1]),float(pList[j][2])])
						PointsSuperList.append(LineCoords)
						LineCoords=[]
					#list of lists of points used to draw each curve
			except:
				lInfo=[]
				PointsSuperList = [[[0,0]]]
				for i in cList:
					PointsSuperList.append([[0,0]])


			'''Symbol Data'''
			try:
				sInfo=[]
				for i in cList:
					sInfo.append(pipe.cmd_in('python plot_symbol '+gType+'.'+i).splitlines())
				#list of points from tao command python plot_symbol for each curve

				SymbolSuperList=[]
				for i in range(len(cList)):
					sList = []
					for j in range(len(sInfo[i])):
						sList.append(sInfo[i][j].split(';'))
					SymCoords = []
					for j in range(len(sList)):
						SymCoords.append([float(sList[j][2]),float(sList[j][3])])
					SymbolSuperList.append(SymCoords)
					SymCoords=[]
				#list of lists of points used to draw symbols on each curve
			except:
				sInfo = []
				SymbolSuperList = []
				for i in cList:
					SymbolSuperList.append([[0,0]])



			'''Plot Data'''

			CurvesList = []
			for i in range(len(cList)):
				CurveData = []
				CurveData.append(PointsSuperList[i])
				CurveData.append(SymbolSuperList[i])
				CurveData.append(cInfoDictList[i]['line.color'].value)
				CurveData.append(cInfoDictList[i]['line.pattern'].value.lower())
				CurveData.append(cInfoDictList[i]['line.width'].value)
				CurveData.append(cInfoDictList[i]['symbol.color'].value)
				if (cInfoDictList[i]['symbol.fill_pattern'] == 'solid_fill'):
					CurveData.append(cInfoDictList[i]['symbol.color'].value)
				else:
					CurveData.append('none')
				if (cInfoDictList[i]['draw_symbols'].value == True):
					CurveData.append(cInfoDictList[i]['symbol.height'].value)
				else:
					CurveData.append(0)
				CurveData.append(SymbolsDict[cInfoDictList[i]['symbol.type'].value])
				#if fill and sizing are not automatic for certain symbols, override here 
				CurveData.append(cInfoDictList[i]['symbol.line_width'].value)
				CurvesList.append(CurveData)
				CurveData = []
			#list of [PointsList,SymbolsList,line color,linetype,linewidth,marker color,marker fill type,markertype,markersize] for each curve
			#example: [[[(4,2),(2,5),(0,1),(3,5),(1,2)],[(4,2),(2,5),(0,1)],'Blue','dashed',2,'Blue','s',5],[[(4,3),(2,4),(0,7),(1,6),(3,3)],[(4,3),(2,4),(0,7)],'Orange','solid',2,'Orange','o',5]]




			'''Plotting'''

			#x=[[(4,2),(2,5),(0,1),(3,5),(1,2)],[(4,2),(2,5),(0,1)],'tab:blue','dashed',2,'tab:blue','s',5]
			#y=[[(4,3),(2,4),(0,7),(1,6),(3,3)],[(4,3),(2,4),(0,7)],'tab:orange','solid',2,'tab:orange','o',5]

			#infofromtao [CurvesList,title,xlabel,ylabel,grid]

			LineList = []
			for i in CurvesList:
				PointsList = i[0]
				SymbolsList = i[1]
				PointsList.sort()
				SymbolsList.sort()
				xpList = []
				ypList = []
				for j in PointsList:
					xpList.append(j[0])
					ypList.append(j[1])
				xsList = []
				ysList = []
				for k in SymbolsList:
					xsList.append(k[0])
					ysList.append(k[1])

				if gInfoDict['graph^type'].value == 'data':
					LineList.append(GraphDict['graph'+str(gNumber+1)].plot(xpList,ypList,color=i[2],linestyle=i[3],linewidth=i[4]))
					GraphDict['graph'+str(gNumber+1)].plot(xsList,ysList,color=i[5],linewidth=0,markerfacecolor=i[6],markersize=i[7]/2,marker=i[8],mew=i[9]/2)
				#line and symbol graphs

				elif gInfoDict['graph^type'].value == 'phase_space':
					LineList.append(GraphDict['graph'+str(gNumber+1)].plot(xsList,ysList,color=i[5],linewidth=0,markerfacecolor=i[6],markersize=i[7]/2,marker=i[8],mew=i[9]/2))
				#phase space graphs

				elif gInfoDict['graph^type'].value == 'histogram':
					LineList.append(GraphDict['graph'+str(gNumber+1)].hist(xpList,bins=100,weights=ypList,histtype='step',color=i[5]))
				#histogram	

			#plot curves
			#LineList gives names of curves




			plt.title(pgp_to_mpl(gInfoDict['title'].value)+' '+gInfoDict['title_suffix'].value)
			#plot title

			LegendList = []
			LabelList = []
			for i in range(len(CurvesList)):
				LegendList.append(LineList[i][0])
				LabelList.append(pgp_to_mpl(cInfoDictList[i]['legend_text'].value))

			if (gInfoDict['draw_curve_legend'].value == True and LabelList != ['']):
				GraphDict['graph'+str(gNumber+1)].legend(LegendList,LabelList)
			#plot legend

			plt.xlabel(pgp_to_mpl(gInfoDict['x.label'].value))
			plt.ylabel(pgp_to_mpl(gInfoDict['y.label'].value))
			#plot axis labels

			xmajorLocator=MultipleLocator((gInfoDict['x.max'].value-gInfoDict['x.min'].value)/gInfoDict['x.major_div'].value)
			ymajorLocator=MultipleLocator((gInfoDict['y.max'].value-gInfoDict['y.min'].value)/gInfoDict['y.major_div'].value)
			GraphDict['graph'+str(gNumber+1)].xaxis.set_major_locator(xmajorLocator)
			GraphDict['graph'+str(gNumber+1)].yaxis.set_major_locator(ymajorLocator)
			GraphDict['graph'+str(gNumber+1)].grid(gInfoDict['draw_grid'].value,which='major',axis='both')
			plt.xlim(gInfoDict['x.min'].value,gInfoDict['x.max'].value)
			plt.ylim(gInfoDict['y.min'].value,gInfoDict['y.max'].value)
			GraphDict['graph'+str(gNumber+1)].set_axisbelow(True)
			#plot grid

		fig.tight_layout()

		return fig
