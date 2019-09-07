import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.backend_tools import ToolBase, ToolToggleBase
from matplotlib.widgets import Slider, Button, RadioButtons

from .tao_interface import *
from pytao.util.parameters import *



class taoplot:
        def __init__(self,pipe,GraphRegion):
                '''initializer, takes a tao interface and a graph region'''
                self.pipe= pipe #tao_interface object
                self.GraphRegion= GraphRegion #string describing region in tao of the desired plot

        def plot(self, width):
                '''returns a figure containing graphs using the data in the region GraphRegion of the tao instance in pipe, and plots a lat_layout below if applicable, also returns information about the indices and locations of elements, width is used to modify the width of floor plan elements and the size of the lat layout'''
                fig = plt.figure()
                #creates plotting figure

                pipe = self.pipe #tao interface instance to plot in
                GraphRegion = self.GraphRegion #graph region to plot in

                LatLayout=False #drawn below data plots if x label is not indices
                FloorPlan=False
                #determines if shapes are drawn on graph, changed to True later if needed

                eleIndexList=[]
                eleStartDict={}
                eleEndDict={}
                fpeIndexList=[]
                fpeCenterDict={}
                fpeRadiusDict={}
                pathDict = {}
                #records information about element locations to be returned with the figure

                def color(x):
                        '''takes string containing pgplot color and returns corresponding matplotlib color'''
                        if x == 'White' or x == 'white':
                                return 'white'
                        elif x == 'Black' or x == 'black':
                                return 'black'
                        elif x == 'Red' or x == 'red':
                                return 'red'
                        elif x == 'Green' or x == 'green':
                                return 'green'
                        elif x == 'Blue' or x == 'blue':
                                return 'blue'
                        elif x == 'Cyan' or x == 'cyan':
                                return 'cyan'
                        elif x == 'Magenta' or x == 'magenta':
                                return 'magenta'
                        elif x == 'Yellow' or x == 'yellow':
                                return 'yellow'
                        elif x == 'Orange' or x == 'orange':
                                return 'orange'
                        elif x == 'Yellow_Green' or x == 'yellow_green':
                                return 'greenyellow'
                        elif x == 'Light_Green' or x == 'light_green':
                                return 'limegreen'
                        elif x == 'Navy_Blue' or x == 'navy_blue':
                                return 'navy'
                        elif x == 'Purple' or x == 'purple':
                                return 'purple'
                        elif x == 'Reddish_Purple' or x == 'reddish_purple' or x == 'Redish_Purple' or x == 'redish_purple':
                                return 'mediumvioletred'
                        elif x == 'Dark_Grey' or x == 'dark_grey':
                                return 'gray'
                        elif x == 'Light_Grey' or x == 'light_grey':
                                return 'lightgray'
                        elif x == 'Transparent' or x == 'transparent':
                                return 'none'
                        else:
                                return x

                def pgp_to_mpl(x):
                        '''takes string with pgplot characters and returns string with characters replaced with matplotlib equivalent, raises NotImplementedError if an unknown pgplot character is used'''
                        x=x.replace('\\','\\\\')
                        if '\\\\' in x:
                                lx = '$'+x+'$'
                                while lx.find('\\\\d') != -1 and lx.find('\\\\u') != -1:
                                        if lx.find('\\\\d') <       lx.find('\\\\u'):
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
                                lx=lx.replace('%','\\%')

                                lx=lx.replace('\\\\(2265)','\\partial')
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



                StylesDict = {
                        'solid':'solid',
                        'dashed':'dashed',
                        'dash_dot':'dashdot',
                        'dotted':'dotted',
                        'dash_dot3':'dashdot', #currently the same as dashdot
                        '1':'solid',
                        '2':'dashed',
                        '3':'dashdot',
                        '4':'dotted',
                        '5':'dashdot'
                        }

                FillDict = {
                        'solid_fill':'full',
                        'no_fill':'none',
                        'hatched':'full',
                        'cross_hatched':'full',
                        '1':'full',
                        '2':'none',
                        '3':'full',
                        '4':'full'
                        }

                SymbolsDict = {
                        'do_not_draw':'',
                        'square':'s', #no fill
                        'dot':'.',
                        'plus':'+',
                        'times':(6,2,0),
                        'circle':'$\\circ$',
                        'x':'x',
                        'x_symbol':'x',
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



                '''''''''Graphing Data'''''''''
                #obtains information about the number and data of graphs from tao

                '''Region Data'''

                rInfo=pipe.cmd_in('python plot1 '+GraphRegion,no_warn = True).splitlines()
                #list of plotting parameter strings from tao command python plot1

                rInfoDict = {}
                for i in range(len(rInfo)):
                        rInfoDict[rInfo[i].split(';')[0]]=str_to_tao_param(rInfo[i])
                #dictionary of tao_parameter name string keys to the corresponding
                #tao_parameter object from python plot_1

                gList = []
                heightsList = []
                for i in range(rInfoDict['num_graphs'].value):
                        gList.append(rInfoDict[('graph['+str(i+1)+']')].value)
                        heightsList.append(1) #makes graphs the same size
                #list of string names of graphs


                number_graphs = len(gList) + 1
                layout_height = .2*width
                #relative height of lat layout to graph if there is one graph

                heightsList.append(layout_height*len(gList) + layout_height*(len(gList)-1))
                #makes layout approximately a fixed size as number of graphs changes

                gs = fig.add_gridspec(nrows=number_graphs,ncols=1, height_ratios=heightsList)

                GraphDict = {} #dictionary of graph name keys with the graph's subplot as it's value
                graph_list = [] #list of graphs, eg: r1.g or r3.x




                for gNumber in range(len(gList)):

                        if gNumber == 0:
                                GraphDict['graph'+str(gNumber+1)]=fig.add_subplot(gs[gNumber,0])
                        elif gNumber > 0:
                                GraphDict['graph'+str(gNumber+1)]=fig.add_subplot(gs[gNumber,0],sharex=GraphDict['graph1'])
                        #create plots in figure, second line also makes x axes scale together



                        '''Graph Data'''

                        gType = GraphRegion+'.'+gList[gNumber]
                        graph_list.append(gType)
                        #graph type, eg: r13.g or top.x


                        gInfo=pipe.cmd_in('python plot_graph '+gType,no_warn = True).splitlines()
                        #list of plotting parameter strings from tao command python plot_graph


                        gInfoDict = {}
                        for i in range(len(gInfo)):
                                gInfoDict[gInfo[i].split(';')[0]]=str_to_tao_param(gInfo[i])
                        #tao_parameter object names from python plot_graph
                        #dictionary of tao_parameter name string keys to the corresponding tao_parameter object




                        '''Curve Data'''

                        cList=[]
                        for i in range(gInfoDict['num_curves'].value):
                                cList.append(gInfoDict[('curve['+str(i+1)+']')].value)
                        #list of string names of curves

                        cInfo=[]
                        for i in cList:
                                cInfo.append(pipe.cmd_in('python plot_curve '+gType+'.'+i,no_warn = True).splitlines())
                        #list of lists of plotting parameter strings from tao command python plot_curve for each curve


                        cInfoDictList=[]
                        for i in range(len(cList)):
                                cInfoDict = {}
                                for j in range(len(cInfo[i])):
                                        cInfoDict[cInfo[i][j].split(';')[0]]=str_to_tao_param(cInfo[i][j])
                                cInfoDictList.append(cInfoDict)
                                cInfoDict = {}
                        #list of dictionaries of tao_parameter name string keys to the corresponding tao_parameter object for each curve





                        '''Line Data'''

                        lInfo=[]
                        try:
                                for i in cList:
                                        lInfo.append(pipe.cmd_in('python plot_line '+gType+'.'+i,no_warn = True).splitlines())
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
                        except: #handle graph with no lines
                                lInfo=[]
                                PointsSuperList = []
                                for i in cList:
                                        PointsSuperList.append([[0,0]])


                        '''Symbol Data'''

                        try:
                                sInfo=[]
                                for i in cList:
                                        sInfo.append(pipe.cmd_in('python plot_symbol '+gType+'.'+i,no_warn = True).splitlines())
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
                        except: #handle graph with no symbols
                                sInfo = []
                                SymbolSuperList = []
                                for i in cList:
                                        SymbolSuperList.append([[0,0]])




                        '''Histogram Data'''

                        try:
                                hInfo=[]
                                for i in cList:
                                        hInfo.append(pipe.cmd_in('python plot_histogram '+gType+'.'+i,no_warn = True).splitlines())
                                hInfoDictList=[]
                                for i in range(len(cList)):
                                        hInfoDict = {}
                                        for j in range(len(hInfo[i])):
                                                hInfoDict[hInfo[i][j].split(';')[0]]=str_to_tao_param(hInfo[i][j])
                                        hInfoDictList.append(hInfoDict)
                                        hInfoDict = {}
                                #list of lists of dictionaries of plot_histogram data for each curve

                        except:
                                hInfo=[]




                        '''Plot Data'''

                        CurvesList = []
                        for i in range(len(cList)):
                                CurveData = []
                                CurveData.append(PointsSuperList[i]) #points for each curve
                                CurveData.append(SymbolSuperList[i]) #symbols for each curve
                                CurveData.append(color(cInfoDictList[i]['line'].get_component('color'))) #line color
                                CurveData.append(StylesDict[cInfoDictList[i]['line'].get_component('pattern').lower()]) #line style
                                if (cInfoDictList[i]['draw_line'].value == True): #line width if drawn
                                        CurveData.append(cInfoDictList[i]['line'].get_component('width'))
                                else:
                                        CurveData.append(0)
                                CurveData.append(color(cInfoDictList[i]['symbol'].get_component('color'))) #symbol color

                                if cInfoDictList[i]['symbol'].get_component('type') == 'dot' or cInfoDictList[i]['symbol'].get_component('type') == '1': #determine if symbol should be filled
                                        CurveData.append(cInfoDictList[i]['symbol'].get_component('color'))
                                elif cInfoDictList[i]['symbol'].get_component('type')[-6:] == 'filled':
                                        CurveData.append(cInfoDictList[i]['symbol'].get_component('color'))
                                elif cInfoDictList[i]['symbol'].get_component('type')[:1] == '-':
                                        CurveData.append(cInfoDictList[i]['symbol'].get_component('color'))
                                elif (FillDict[cInfoDictList[i]['symbol'].get_component('fill_pattern')] == 'solid'):
                                        CurveData.append(cInfoDictList[i]['symbol'].get_component('color'))
                                else:
                                        CurveData.append('none')

                                if (cInfoDictList[i]['draw_symbols'].value == True) and SymbolsDict[cInfoDictList[i]['symbol'].get_component('type')] != '': #symbol size if drawn
                                        CurveData.append(cInfoDictList[i]['symbol'].get_component('height'))
                                else:
                                        CurveData.append(0)

                                if SymbolsDict[cInfoDictList[i]['symbol'].get_component('type')] != '':
                                        CurveData.append(SymbolsDict[cInfoDictList[i]['symbol'].get_component('type')]) #symbol type
                                else:
                                        CurveData.append('.')
                                CurveData.append(cInfoDictList[i]['symbol'].get_component('line_width')) #symbol line width
                                CurvesList.append(CurveData)
                                CurveData = []
                        #list of data needed to plot line and symbol graphs
                        #includes points, and line and symbol information for each curve




                        '''''''''Plotting'''''''''
                        #plots line graphs, histograms, phase space plots, dynamic aperture graphs, and plots for wave analysis

                        LineList = []
                        for i in CurvesList:
                                PointsList = i[0]
                                SymbolsList = i[1]
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
                                try:
                                        yMax = max(.5*max(max(ypList),max(ysList)),2*max(max(ypList),max(ysList)))
                                        yMin = min(.5*min(min(ypList),min(ysList)),2*min(min(ypList),min(ysList)))
                                except ValueError:
                                        raise ValueError('no points found, make sure data is properly initialized')
                                #boundaries for wave analysis rectangles

                                if gInfoDict['graph^type'].value == 'data' or gInfoDict['graph^type'].value == 'wave.0' or gInfoDict['graph^type'].value == 'wave.a' or gInfoDict['graph^type'].value == 'wave.b':
                                        LineList.append(GraphDict['graph'+str(gNumber+1)].plot(xpList,ypList,color=i[2],linestyle=i[3],linewidth=i[4]/2))
                                        GraphDict['graph'+str(gNumber+1)].plot(xsList,ysList,color=i[5],linewidth=0,markerfacecolor=i[6],markersize=i[7]/2,marker=i[8],mew=i[9]/2)

                                        if rInfoDict['x_axis_type'].value == 's':
                                                LatLayout = True
                                        #add lat layout if x axis matches layout axis

                                        if gInfoDict['graph^type'].value != 'data': #wave analysis rectangles
                                                wInfo = pipe.cmd_in('python wave params',no_warn = True).splitlines()
                                                a1 = float(wInfo[1].split(';')[3])
                                                a2 = float(wInfo[2].split(';')[3])
                                                b1 = float(wInfo[3].split(';')[3])
                                                b2 = float(wInfo[4].split(';')[3])
                                                #wave region boundaries

                                        if i[5] == 'blue' or i[5] == 'navy' or i[5] == 'cyan' or i[5] == 'green' or i[5] == 'purple':
                                                waveColor = 'orange'
                                        else:
                                                waveColor = 'blue'
                                        #wave analysis rectangle color

                                        if gInfoDict['graph^type'].value == 'wave.0' or gInfoDict['graph^type'].value == 'wave.a':
                                                GraphDict['graph'+str(gNumber+1)].add_patch(patches.Rectangle((a1,yMin),a2-a1,yMax-yMin,fill=False,color=waveColor))
                                        if gInfoDict['graph^type'].value == 'wave.0' or gInfoDict['graph^type'].value == 'wave.b':
                                                GraphDict['graph'+str(gNumber+1)].add_patch(patches.Rectangle((b1,yMin),b2-b1,yMax-yMin,fill=False,color=waveColor))
                                #line and symbol graphs

                                elif gInfoDict['graph^type'].value == 'dynamic_aperture':
                                        LineList.append(GraphDict['graph'+str(gNumber+1)].plot(xpList,ypList,color=i[2],linestyle=i[3],linewidth=i[4]/2))
                                        GraphDict['graph'+str(gNumber+1)].plot(xsList,ysList,color=i[5],linewidth=0,markerfacecolor=i[6],markersize=i[7]/2,marker=i[8],mew=i[9]/2)
                                #dynamic aperture graphs

                                elif gInfoDict['graph^type'].value == 'phase_space':
                                        if lInfo != []:
                                                LineList.append(GraphDict['graph'+str(gNumber+1)].plot(xpList,ypList,color=i[2],linestyle=i[3],linewidth=i[4]/2))
                                                GraphDict['graph'+str(gNumber+1)].plot(xsList,ysList,color=i[5],linewidth=0,markerfacecolor=i[6],markersize=i[7]/2,marker=i[8],mew=i[9]/2)
                                        else:
                                                LineList.append(GraphDict['graph'+str(gNumber+1)].plot(xsList,ysList,color=i[5],linewidth=0,markerfacecolor=i[6],markersize=i[7]/2,marker=i[8],mew=i[9]/2))
                                #phase space graphs

                                elif gInfoDict['graph^type'].value == 'histogram':
                                        LineList.append(GraphDict['graph'+str(gNumber+1)].hist(xpList,bins=int(hInfoDictList[CurvesList.index(i)]['number'].value),weights=ypList,histtype='step',color=i[5]))
                                #histogram

                        if gInfoDict['graph^type'].value == 'key_table':
                                raise NotImplementedError('key table is not available in the GUI')

                        if gInfoDict['graph^type'].value == 'lat_layout':
                                LatLayout = True
                                gList = []
                                plt.axis('off')
                        #sets up lat layout plot

                        if gInfoDict['graph^type'].value == 'floor_plan':
                                FloorPlan = True
                                gList = []
                                plt.axis('off')
                        #sets up floor plan plot

                        if gInfoDict['draw_axes'].value == False:
                                plt.axis('off')
                        #hides axes if draw_axes is turned off

                        if gInfoDict['why_invalid'].value != '':
                                raise ValueError(gInfoDict['why_invalid'].value + ', make sure graph is properly initialized')

                        #plots line and symbol graphs and histograms, lat layouts and floor plans are drawn later
                        #LineList gives names of curves




                        plt.title(pgp_to_mpl(gInfoDict['title'].value)+' '+gInfoDict['title_suffix'].value)
                        #plot title

                        LegendList = [] #legends for each graph
                        LabelList = [] #labels in each legend
                        try:
                                for i in range(len(CurvesList)):
                                        LegendList.append(LineList[i][0])
                                        if pgp_to_mpl(cInfoDictList[i]['legend_text'].value) != '':
                                                LabelList.append(pgp_to_mpl(cInfoDictList[i]['legend_text'].value))
                                        elif pgp_to_mpl(cInfoDictList[i]['data_type'].value) == 'physical_aperture':
                                                LabelList.append(pgp_to_mpl(cInfoDictList[i]['data_type'].value))
                                        else:
                                                LabelList.append('')
                                #list of curves to be added to a legend and list of labels for each curve in the legend
                        except IndexError:
                                raise NotImplementedError('unknown graph type')


                        if (gInfoDict['draw_curve_legend'].value == True and LabelList != ['']) and gInfoDict['graph^type'].value != 'lat_layout' and gInfoDict['graph^type'].value != 'floor_plan':
                                GraphDict['graph'+str(gNumber+1)].legend(LegendList,LabelList)
                        #plot legend

                        plt.xlabel(pgp_to_mpl(gInfoDict['x'].get_component('label')))
                        plt.ylabel(pgp_to_mpl(gInfoDict['y'].get_component('label')))
                        #plot axis labels

                        GraphDict['graph'+str(gNumber+1)].grid(gInfoDict['draw_grid'].value,which='major',axis='both')
                        #plot grid

                        plt.xlim(gInfoDict['x'].get_component('min'),gInfoDict['x'].get_component('max'))
                        plt.ylim(gInfoDict['y'].get_component('min'),gInfoDict['y'].get_component('max'))
                        #set axis limits

                        GraphDict['graph'+str(gNumber+1)].set_axisbelow(True)
                        #place graphs over grid lines





                '''''''''Lattice Layout'''''''''
                #plots lat layouts

                if LatLayout == True:
                        if gInfoDict['graph^type'].value != 'lat_layout': #add space for lat layout below graph
                                GraphDict['LatLayout']=fig.add_subplot(gs[len(gList),0],sharex=GraphDict['graph1'])

                        else: #standalone lat layout graph
                                GraphDict['LatLayout']=fig.add_subplot(len(gList)+1,1,len(gList)+1,sharex=GraphDict['graph1'])

                        layInfo=pipe.cmd_in('python plot_graph layout.g',no_warn = True).splitlines()
                        #list of plotting parameter strings from tao command python plot_graph

                        layInfoDict = {}
                        for i in range(len(layInfo)):
                                layInfoDict[layInfo[i].split(';')[0]]=str_to_tao_param(layInfo[i])
                        #tao_parameter object names from python plot_graph
                        #dictionary of tao_parameter name string keys to the corresponding tao_parameter object



                        twinAxes=GraphDict['LatLayout'].axes.twinx()
                        plt.xlim(gInfoDict['x'].get_component('min'),gInfoDict['x'].get_component('max'))
                        plt.ylim(layInfoDict['y'].get_component('min'),layInfoDict['y'].get_component('max'))
                        twinAxes.set_navigate(True)
                        GraphDict['LatLayout'].axis('off')
                        twinAxes.axis('off')
                        #makes lat layout only have horizontal axis for panning and zooming

                        GraphDict['LatLayout'].axes.set_navigate(False)
                        GraphDict['LatLayout'].axhline(y=0,xmin=1.1*layInfoDict['x'].get_component('min'),xmax=1.1*layInfoDict['x'].get_component('max'),color='Black')
                        #sets axis limits and creates second axis to allow x panning and zooming


                        if layInfoDict['ix_universe'].value != -1:
                                universe = layInfoDict[ix_universe].value
                        else:
                                universe = 1
                        branch = layInfoDict['-1^ix_branch'].value
                        #lat layout branch and universe information


                        eleInfo=pipe.cmd_in('python plot_lat_layout '+str(universe)+'@'+str(branch),no_warn = True).splitlines()
                        #list of strings containing information about each element

                        eleIndexList = []
                        eleStartDict = {} #element start coordinate
                        eleEndDict = {} #element end coordinate
                        eleLwDict = {} #element line width
                        eleShapeDict = {} #element shape (eg: box)
                        eleY1Dict = {} #height above 0
                        eleY2Dict = {} #height below 0
                        eleColorDict = {}
                        eleNameDict = {}
                        for i in range(len(eleInfo)):
                                eleIndexList.append(int(eleInfo[i].split(';')[0]))
                                eleStartDict[eleInfo[i].split(';')[0]]= float(eleInfo[i].split(';')[1])
                                eleEndDict[eleInfo[i].split(';')[0]]= float(eleInfo[i].split(';')[2])
                                eleLwDict[eleInfo[i].split(';')[0]]= float(eleInfo[i].split(';')[3])
                                eleShapeDict[eleInfo[i].split(';')[0]]= eleInfo[i].split(';')[4].lower()
                                eleY1Dict[eleInfo[i].split(';')[0]]= float(eleInfo[i].split(';')[5])
                                eleY2Dict[eleInfo[i].split(';')[0]]= float(eleInfo[i].split(';')[6])
                                eleColorDict[eleInfo[i].split(';')[0]]= color(eleInfo[i].split(';')[7].lower())
                                eleNameDict[eleInfo[i].split(';')[0]]= eleInfo[i].split(';')[8]
                        #all dict keys and entries are strings which match a lattice layout element index (eg: '1') string to the corresponding information


                        for i in eleIndexList:

                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[eleY1Dict[str(i)],-1.7*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)],alpha=0)
                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[-1.7*eleY2Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)],alpha=0)
                                #invisible elements to give enough space for actual lat layout elements

                                try:
                                        if eleShapeDict[str(i)] == 'box' and eleEndDict[str(i)]-eleStartDict[str(i)] > 0:
                                                GraphDict['LatLayout'].add_patch(patches.Rectangle((eleStartDict[str(i)],-1*eleY2Dict[str(i)]),eleEndDict[str(i)]-eleStartDict[str(i)],eleY1Dict[str(i)]+eleY2Dict[str(i)],lw=eleLwDict[str(i)],color=eleColorDict[str(i)],fill=False))
                                        #draw box element




                                        elif eleShapeDict[str(i)] == 'xbox' and eleEndDict[str(i)]-eleStartDict[str(i)] > 0:
                                                GraphDict['LatLayout'].add_patch(patches.Rectangle((eleStartDict[str(i)],-1*eleY2Dict[str(i)]),eleEndDict[str(i)]-eleStartDict[str(i)],eleY1Dict[str(i)]+eleY2Dict[str(i)],lw=eleLwDict[str(i)],color=eleColorDict[str(i)],fill=False))
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[eleY1Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[-1*eleY2Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                        #draw xbox element


                                        elif eleShapeDict[str(i)] == 'x' and eleEndDict[str(i)]-eleStartDict[str(i)] > 0:
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[eleY1Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[-1*eleY2Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                        #draw x element


                                        elif eleShapeDict[str(i)] == 'bow_tie' and eleEndDict[str(i)]-eleStartDict[str(i)] > 0:
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[eleY1Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[-1*eleY2Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[eleY1Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleEndDict[str(i)]],[-1*eleY2Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                        #draw bow_tie element


                                        elif eleShapeDict[str(i)] == 'diamond' and eleEndDict[str(i)]-eleStartDict[str(i)] > 0:
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleStartDict[str(i)]+(eleEndDict[str(i)]-eleStartDict[str(i)])/2],[0,-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleStartDict[str(i)]+(eleEndDict[str(i)]-eleStartDict[str(i)])/2],[0,eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)]+(eleEndDict[str(i)]-eleStartDict[str(i)])/2,eleEndDict[str(i)]],[-1*eleY2Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)]+(eleEndDict[str(i)]-eleStartDict[str(i)])/2,eleEndDict[str(i)]],[eleY1Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                        #draw diamond element


                                        elif eleShapeDict[str(i)] == 'circle':
                                                GraphDict['LatLayout'].add_patch(patches.Ellipse((eleStartDict[str(i)]+(eleEndDict[str(i)]-eleStartDict[str(i)])/2,0),eleY1Dict[str(i)]+eleY2Dict[str(i)],eleY1Dict[str(i)]+eleY2Dict[str(i)],lw=eleLwDict[str(i)],color=eleColorDict[str(i)],fill=False))
                                        #draw circle element


                                        elif eleShapeDict[str(i)] == 'box' and eleEndDict[str(i)]-eleStartDict[str(i)] < 0:
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[eleY1Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[eleY1Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[-1*eleY2Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[-1*eleY2Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleStartDict[str(i)]],[eleY1Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleEndDict[str(i)],eleEndDict[str(i)]],[eleY1Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                        #draw wrapped box element


                                        elif eleShapeDict[str(i)] == 'xbox' and eleEndDict[str(i)]-eleStartDict[str(i)] < 0:
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[eleY1Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[eleY1Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[-1*eleY2Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[-1*eleY2Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],eleStartDict[str(i)]],[eleY1Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleEndDict[str(i)],eleEndDict[str(i)]],[eleY1Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[eleY1Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[0,eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[-1*eleY2Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[0,-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                        #draw wrapped xbox element


                                        elif eleShapeDict[str(i)] == 'x' and eleEndDict[str(i)]-eleStartDict[str(i)] < 0:
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[eleY1Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[0,eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[-1*eleY2Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[0,-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                        #draw wrapped x element


                                        elif eleShapeDict[str(i)] == 'bow_tie' and eleEndDict[str(i)]-eleStartDict[str(i)] < 0:
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[eleY1Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[eleY1Dict[str(i)],eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[-1*eleY2Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[-1*eleY2Dict[str(i)],-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[eleY1Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[0,eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[-1*eleY2Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[0,-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                        #draw wrapped bow tie element


                                        elif eleShapeDict[str(i)] == 'diamond' and eleEndDict[str(i)]-eleStartDict[str(i)] < 0:
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[0,eleY1Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[eleY1Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([eleStartDict[str(i)],layInfoDict['x'].get_component('max')],[0,-1*eleY2Dict[str(i)]],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].plot([layInfoDict['x'].get_component('min'),eleEndDict[str(i)]],[-1*eleY2Dict[str(i)],0],lw=eleLwDict[str(i)],color=eleColorDict[str(i)])
                                        #draw wrapped diamond element


                                        if eleEndDict[str(i)]-eleStartDict[str(i)] > 0:
                                                GraphDict['LatLayout'].text(eleStartDict[str(i)]+(eleEndDict[str(i)]-eleStartDict[str(i)])/2,-1.1*eleY2Dict[str(i)],eleNameDict[str(i)],ha='center',va='top',clip_on=True,color=eleColorDict[str(i)])
                                        #draw element name


                                        elif eleEndDict[str(i)]-eleStartDict[str(i)] < 0:
                                                GraphDict['LatLayout'].text(layInfoDict['x'].get_component('max'),-1.1*eleY2Dict[str(i)],eleNameDict[str(i)],ha='right',va='top',clip_on=True,color=eleColorDict[str(i)])
                                                GraphDict['LatLayout'].text(layInfoDict['x'].get_component('min'),-1.1*eleY2Dict[str(i)],eleNameDict[str(i)],ha='left',va='top',clip_on=True,color=eleColorDict[str(i)])
                                        #draw wrapped element name



                                except KeyError:
                                        pass


                else:
                        GraphDict['LatLayout']=fig.add_subplot(gs[len(gList),0],sharex=GraphDict['graph1'])
                        GraphDict['LatLayout'].remove()



                '''''''''Floor Plan'''''''''
                #plots floor plans

                if FloorPlan == True:
                        GraphDict['FloorPlan']=fig.add_subplot(len(gList)+1,1,len(gList)+1,sharex=GraphDict['graph1'])

                        floInfo=pipe.cmd_in('python plot_graph '+gType,no_warn = True).splitlines()
                        #list of plotting parameter strings from tao command python plot_graph


                        floInfoDict = {}
                        for i in range(len(floInfo)):
                                floInfoDict[floInfo[i].split(';')[0]]=str_to_tao_param(floInfo[i])
                        #tao_parameter object names from python plot_graph for a floor plan
                        #dictionary of tao_parameter name string keys to the corresponding tao_parameter object

                        if floInfoDict['ix_universe'].value != -1:
                                universe = floInfoDict[ix_universe].value

                        else:
                                universe = 1

                        fpeInfo=pipe.cmd_in('python floor_plan '+gType,no_warn = True).splitlines()
                        #list of plotting parameter strings from tao command python floor_plan


                        fpeIndexList = [] #contains lists of branch index then element index for each point
                        fpeTypeDict = {}
                        fpeSxDict = {} #start x coordinate
                        fpeSyDict = {} #start y coordinate
                        fpeSaDict = {} #start angle
                        fpeExDict = {} #end x coordinate
                        fpeEyDict = {} #end y coordinate
                        fpeEaDict = {} #end angle
                        fpeLwDict = {} #line width
                        fpeShapeDict = {}
                        fpeY1Dict = {} #distance above
                        fpeY2Dict = {} #distance below
                        fpeColorDict = {}
                        fpeNameDict = {}
                        fpeAlDict = {} #arc length
                        fpeBaDict = {} #bend angle
                        fpeSfaDict = {} #relative angle of starting face to incoming line
                        fpeEfaDict = {} #relative angle of ending face to incoming line
                        corner1 = {} #corner coordinates of objects
                        corner2 = {}
                        corner3 = {}
                        corner4 = {}
                        for i in range(len(fpeInfo)):
                                fpeIndexList.append([int(fpeInfo[i].split(';')[0]),int(fpeInfo[i].split(';')[1])])
                                fpeTypeDict[str(fpeIndexList[i])]= fpeInfo[i].split(';')[2].lower()
                                fpeSxDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[3])
                                fpeSyDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[4])
                                fpeSaDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[5])
                                fpeExDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[6])
                                fpeEyDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[7])
                                fpeEaDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[8])
                                fpeLwDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[9])
                                fpeShapeDict[str(fpeIndexList[i])]= fpeInfo[i].split(';')[10].lower()
                                fpeY1Dict[str(fpeIndexList[i])]= width*float(fpeInfo[i].split(';')[11])
                                fpeY2Dict[str(fpeIndexList[i])]= width*float(fpeInfo[i].split(';')[12])
                                fpeColorDict[str(fpeIndexList[i])]= color(fpeInfo[i].split(';')[13].lower())
                                fpeNameDict[str(fpeIndexList[i])]= fpeInfo[i].split(';')[14]
                                try:
                                        fpeAlDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[15])
                                        fpeBaDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[16])
                                        fpeSfaDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[17])
                                        fpeEfaDict[str(fpeIndexList[i])]= float(fpeInfo[i].split(';')[18])
                                except IndexError:
                                        pass
                        #dict keys and entries are strings which match a floor plan element index (eg: '1') to the corresponding information

                        def line(p1, p2):
                                '''returns lines based on given points to be used with intersect'''
                                A = (p1[1] - p2[1])
                                B = (p2[0] - p1[0])
                                C = (p1[0]*p2[1] - p2[0]*p1[1])
                                return A, B, -C

                        def intersect(L1, L2):
                                '''returns intersection point of 2 lines from the line funciton, or false if the lines don't intersect'''
                                D = L1[0]*L2[1] - L1[1]*L2[0]
                                Dx = L1[2]*L2[1] - L1[1]*L2[2]
                                Dy = L1[0]*L2[2] - L1[2]*L2[0]

                                if D != 0:
                                        x = Dx / D
                                        y = Dy / D
                                        return x,y
                                else:
                                        return False



                        conv = (180)/(np.pi) #radian to degree conversion
                        for i in fpeIndexList:
                                fpeCenterDict[str(i)]=([fpeSxDict[str(i)] + (fpeExDict[str(i)]-fpeSxDict[str(i)])/2,fpeSyDict[str(i)] + (fpeEyDict[str(i)]-fpeSyDict[str(i)])/2])
                                fpeRadiusDict[str(i)]=fpeY1Dict[str(i)] #for click detection

                                try:
                                        if fpeTypeDict[str(i)] == 'drift' or fpeTypeDict[str(i)] == 'kicker':
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)],fpeExDict[str(i)]],[fpeSyDict[str(i)],fpeEyDict[str(i)]],color='black')
                                        #draw drift element



                                        if fpeY1Dict[str(i)] == 0 and fpeY2Dict[str(i)] == 0 and fpeTypeDict[str(i)] != 'sbend' and fpeColorDict[str(i)] != '':
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)],fpeExDict[str(i)]],[fpeSyDict[str(i)],fpeEyDict[str(i)]],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                        #draw line element



                                        elif fpeShapeDict[str(i)] == 'box' and fpeTypeDict[str(i)] != 'sbend' and fpeColorDict[str(i)] != '':
                                                GraphDict['FloorPlan'].add_patch(patches.Rectangle((fpeSxDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeSyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])),np.sqrt((fpeExDict[str(i)]-fpeSxDict[str(i)])**2 + (fpeEyDict[str(i)]-fpeSyDict[str(i)])**2),fpeY1Dict[str(i)]+fpeY2Dict[str(i)],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)],fill=False,angle=fpeSaDict[str(i)]*conv))
                                        #draw box element


                                        elif fpeShapeDict[str(i)] == 'xbox' and fpeTypeDict[str(i)] != 'sbend' and fpeColorDict[str(i)] != '':
                                                GraphDict['FloorPlan'].add_patch(patches.Rectangle((fpeSxDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeSyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])),np.sqrt((fpeExDict[str(i)]-fpeSxDict[str(i)])**2 + (fpeEyDict[str(i)]-fpeSyDict[str(i)])**2),fpeY1Dict[str(i)]+fpeY2Dict[str(i)],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)],fill=False,angle=fpeSaDict[str(i)]*conv))
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                        #draw xbox element


                                        elif fpeShapeDict[str(i)] == 'x' and fpeTypeDict[str(i)] != 'sbend' and fpeColorDict[str(i)] != '':
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                        #draw x element


                                        elif fpeShapeDict[str(i)] == 'bow_tie' and fpeTypeDict[str(i)] != 'sbend' and fpeColorDict[str(i)] != '':
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                        #draw bow_tie element


                                        elif fpeShapeDict[str(i)] == 'diamond' and fpeTypeDict[str(i)] != 'sbend' and fpeColorDict[str(i)] != '':
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)],fpeSxDict[str(i)] + (fpeExDict[str(i)]-fpeSxDict[str(i)])/2 - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)],fpeSyDict[str(i)] + (fpeEyDict[str(i)]-fpeSyDict[str(i)])/2 + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] + (fpeExDict[str(i)]-fpeSxDict[str(i)])/2 - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)]],[fpeSyDict[str(i)] + (fpeEyDict[str(i)]-fpeSyDict[str(i)])/2 + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)]],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)],fpeSxDict[str(i)] + (fpeExDict[str(i)]-fpeSxDict[str(i)])/2 + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)])],[fpeSyDict[str(i)],fpeSyDict[str(i)] + (fpeEyDict[str(i)]-fpeSyDict[str(i)])/2 - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)] + (fpeExDict[str(i)]-fpeSxDict[str(i)])/2 + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeExDict[str(i)]],[fpeSyDict[str(i)] + (fpeEyDict[str(i)]-fpeSyDict[str(i)])/2 - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeEyDict[str(i)]],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                        #draw diamond element


                                        elif fpeShapeDict[str(i)] == 'circle' and fpeTypeDict[str(i)] != 'sbend' and fpeColorDict[str(i)] != '':
                                                GraphDict['FloorPlan'].add_patch(patches.Circle((fpeSxDict[str(i)] + (fpeExDict[str(i)]-fpeSxDict[str(i)])/2,fpeSyDict[str(i)] + (fpeEyDict[str(i)]-fpeSyDict[str(i)])/2),fpeY1Dict[str(i)],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)],fill=False))
                                        #draw circle element


                                        elif fpeShapeDict[str(i)] == 'box' and fpeTypeDict[str(i)] == 'sbend' and fpeColorDict[str(i)] != '':
                                                GraphDict['FloorPlan'].plot([fpeSxDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)]),fpeSxDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])],[fpeSyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)]),fpeSyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                GraphDict['FloorPlan'].plot([fpeExDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeEaDict[str(i)]+fpeEfaDict[str(i)]),fpeExDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeEaDict[str(i)]+fpeEfaDict[str(i)])],[fpeEyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeEaDict[str(i)]+fpeEfaDict[str(i)]),fpeEyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeEaDict[str(i)]+fpeEfaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                #draws straight sbend edges

                                                intersection = intersect(line([fpeSxDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeSyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)])],[fpeSxDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeSyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])]), line([fpeExDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeEaDict[str(i)]),fpeEyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeEaDict[str(i)])],[fpeExDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeEaDict[str(i)]),fpeEyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeEaDict[str(i)]+fpeEfaDict[str(i)])]))
                                                #center of circle used to draw arc edges of sbends

                                                if intersection == False:
                                                        GraphDict['FloorPlan'].plot([fpeSxDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)]),fpeExDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeEaDict[str(i)]+fpeEfaDict[str(i)])],[fpeSyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)]),fpeEyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeEaDict[str(i)]+fpeEfaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                        GraphDict['FloorPlan'].plot([fpeSxDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)]),fpeExDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeEaDict[str(i)]+fpeEfaDict[str(i)])],[fpeSyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)]),fpeEyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeEaDict[str(i)]+fpeEfaDict[str(i)])],lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)])
                                                #draw sbend edges if bend angle is 0


                                                elif intersection != False:
                                                        angle1 = 360 + conv*np.arctan2(fpeSyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1],fpeSxDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])
                                                        angle2 = 360 + conv*np.arctan2(fpeEyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeEaDict[str(i)]+fpeEfaDict[str(i)])-intersection[1],fpeExDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeEaDict[str(i)]+fpeEfaDict[str(i)])-intersection[0])
                                                        #angles of further curve endpoints relative to center of circle
                                                        angle3 = 360 + conv*np.arctan2(fpeSyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1],fpeSxDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])
                                                        angle4 = 360 + conv*np.arctan2(fpeEyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeEaDict[str(i)]+fpeEfaDict[str(i)])-intersection[1],fpeExDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeEaDict[str(i)]+fpeEfaDict[str(i)])-intersection[0])
                                                        #angles of closer curve endpoints relative to center of circle

                                                        if abs(angle1-angle2)<180 and abs(angle3-angle4)<180:
                                                                if angle1 > angle2 and angle3 > angle4:
                                                                        a1=angle2
                                                                        a2=angle1
                                                                        a3=angle4
                                                                        a4=angle3
                                                                elif angle1 < angle2 and angle3 > angle4:
                                                                        a1=angle1
                                                                        a2=angle2
                                                                        a3=angle4
                                                                        a4=angle3
                                                                elif angle1 > angle2 and angle3 < angle4:
                                                                        a1=angle2
                                                                        a2=angle1
                                                                        a3=angle3
                                                                        a4=angle4
                                                                else:
                                                                        a1=angle1
                                                                        a2=angle2
                                                                        a3=angle3
                                                                        a4=angle4

                                                        elif abs(angle1-angle2)>180 and abs(angle3-angle4)<180:
                                                                if angle1 > angle2 and angle3 > angle4:
                                                                        a1=angle1
                                                                        a2=angle2
                                                                        a3=angle4
                                                                        a4=angle3
                                                                elif angle1 < angle2 and angle3 > angle4:
                                                                        a1=angle2
                                                                        a2=angle1
                                                                        a3=angle4
                                                                        a4=angle3
                                                                elif angle1 > angle2 and angle3 < angle4:
                                                                        a1=angle1
                                                                        a2=angle2
                                                                        a3=angle3
                                                                        a4=angle4
                                                                else:
                                                                        a1=angle2
                                                                        a2=angle1
                                                                        a3=angle3
                                                                        a4=angle4

                                                        elif abs(angle1-angle2)<180 and abs(angle3-angle4)>180:
                                                                if angle1 > angle2 and angle3 > angle4:
                                                                        a1=angle2
                                                                        a2=angle1
                                                                        a3=angle3
                                                                        a4=angle4
                                                                elif angle1 < angle2 and angle3 > angle4:
                                                                        a1=angle1
                                                                        a2=angle2
                                                                        a3=angle3
                                                                        a4=angle4
                                                                elif angle1 > angle2 and angle3 < angle4:
                                                                        a1=angle2
                                                                        a2=angle1
                                                                        a3=angle4
                                                                        a4=angle3
                                                                else:
                                                                        a1=angle1
                                                                        a2=angle2
                                                                        a3=angle4
                                                                        a4=angle3

                                                        else:
                                                                if angle1 > angle2 and angle3 > angle4:
                                                                        a1=angle1
                                                                        a2=angle2
                                                                        a3=angle3
                                                                        a4=angle4
                                                                elif angle1 < angle2 and angle3 > angle4:
                                                                        a1=angle2
                                                                        a2=angle1
                                                                        a3=angle3
                                                                        a4=angle4
                                                                elif angle1 > angle2 and angle3 < angle4:
                                                                        a1=angle1
                                                                        a2=angle2
                                                                        a3=angle4
                                                                        a4=angle3
                                                                else:
                                                                        a1=angle2
                                                                        a2=angle1
                                                                        a3=angle4
                                                                        a4=angle3
                                                        #determines correct start and end angles for arcs

                                                        GraphDict['FloorPlan'].add_patch(patches.Arc((intersection[0],intersection[1]),np.sqrt((fpeSxDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])**2 + (fpeSyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1])**2)*2,np.sqrt((fpeSxDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])**2 + (fpeSyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1])**2)*2,theta1=a1,theta2=a2,lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)]))
                                                        GraphDict['FloorPlan'].add_patch(patches.Arc((intersection[0],intersection[1]),np.sqrt((fpeSxDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])**2 + (fpeSyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1])**2)*2,np.sqrt((fpeSxDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])**2 + (fpeSyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1])**2)*2,theta1=a3,theta2=a4,lw=fpeLwDict[str(i)],color=fpeColorDict[str(i)]))
                                                        #draw sbend edges if bend angle is nonzero
                                        #draw sbend element


                                        if fpeNameDict[str(i)] != '' and fpeColorDict[str(i)] != '' and np.sin(((fpeEaDict[str(i)]+fpeSaDict[str(i)])/2)) > 0:
                                                GraphDict['FloorPlan'].text(fpeSxDict[str(i)]+(fpeExDict[str(i)]-fpeSxDict[str(i)])/2 - 1.3*fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeSyDict[str(i)]+(fpeEyDict[str(i)]-fpeSyDict[str(i)])/2 + 1.3*fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeNameDict[str(i)],ha='right',va='center',color='black',rotation=-90+((fpeEaDict[str(i)]+fpeSaDict[str(i)])/2)*conv,clip_on=True,rotation_mode='anchor')

                                        elif fpeNameDict[str(i)] != '' and fpeColorDict[str(i)] != '' and np.sin(((fpeEaDict[str(i)]+fpeSaDict[str(i)])/2)) <= 0:
                                                GraphDict['FloorPlan'].text(fpeSxDict[str(i)]+(fpeExDict[str(i)]-fpeSxDict[str(i)])/2 - 1.3*fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeSyDict[str(i)]+(fpeEyDict[str(i)]-fpeSyDict[str(i)])/2 + 1.3*fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]),fpeNameDict[str(i)],ha='left',va='center',color='black',rotation=90+((fpeEaDict[str(i)]+fpeSaDict[str(i)])/2)*conv,clip_on=True,rotation_mode='anchor')
                                        #draw element name




                                        '''floor plan click detection'''

                                        if fpeTypeDict[str(i)] == 'sbend' and intersection != False: #for sbend click detection
                                                c1 = [fpeSxDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)]),fpeSyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])]
                                                c2 = [fpeExDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeEaDict[str(i)]+fpeEfaDict[str(i)]),fpeEyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeEaDict[str(i)]+fpeEfaDict[str(i)])]
                                                c3 = [fpeSxDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)]),fpeSyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])]
                                                c4 = [fpeExDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeEaDict[str(i)]+fpeEfaDict[str(i)]),fpeEyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeEaDict[str(i)]+fpeEfaDict[str(i)])]
                                                #corners of sbend

                                                if fpeSaDict[str(i)] > fpeEaDict[str(i)]:
                                                        outerRadius = np.sqrt((fpeSxDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])**2 + (fpeSyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1])**2)
                                                        innerRadius = np.sqrt((fpeSxDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])**2 + (fpeSyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1])**2)
                                                else:
                                                        outerRadius = -np.sqrt((fpeSxDict[str(i)]-fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])**2 + (fpeSyDict[str(i)]+fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1])**2)
                                                        innerRadius = -np.sqrt((fpeSxDict[str(i)]+fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[0])**2 + (fpeSyDict[str(i)]-fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)]-fpeSfaDict[str(i)])-intersection[1])**2)
                                                #radii of sbend arc edges



                                                middleAngle = (fpeSaDict[str(i)]+fpeEaDict[str(i)])/2

                                                top = [intersection[0]-outerRadius*np.sin(middleAngle),intersection[1]+outerRadius*np.cos(middleAngle)]
                                                bottom = [intersection[0]-innerRadius*np.sin(middleAngle),intersection[1]+innerRadius*np.cos(middleAngle)]
                                                #midpoints of top and bottom arcs in an sbend

                                                topCP = [2*(top[0])-.5*(c1[0])-.5*(c2[0]),2*(top[1])-.5*(c1[1])-.5*(c2[1])]
                                                bottomCP = [2*(bottom[0])-.5*(c3[0])-.5*(c4[0]),2*(bottom[1])-.5*(c3[1])-.5*(c4[1])]
                                                #corresponding control points for a quadratic Bezier curve that passes through the corners and arc midpoint

                                                verts = [c1,topCP,c2,c4,bottomCP,c3,c1]
                                                codes = [Path.MOVETO,Path.CURVE3,Path.CURVE3,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CLOSEPOLY]
                                                pathDict[str(i)] = Path(verts,codes)

                                                '''patch = patches.PathPatch(Path(verts,codes),facecolor='green',alpha = .5)
                                                GraphDict['FloorPlan'].add_patch(patch)'''
                                                #visualize clickable regions
                                        #path approximating sbend region for clickable region on graph using lines and quadratic Bezier curves

                                        else: #for non sbend click detection
                                                corner1[str(i)] = [fpeSxDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeSyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)])]
                                                corner2[str(i)] = [fpeExDict[str(i)] - fpeY1Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeEyDict[str(i)] + fpeY1Dict[str(i)]*np.cos(fpeSaDict[str(i)])]
                                                corner3[str(i)] = [fpeSxDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeSyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])]
                                                corner4[str(i)] = [fpeExDict[str(i)] + fpeY2Dict[str(i)]*np.sin(fpeSaDict[str(i)]),fpeEyDict[str(i)] - fpeY2Dict[str(i)]*np.cos(fpeSaDict[str(i)])]
                                        #coordinates of corners of a floor plan element for clickable region


                                except KeyError:
                                        pass

                                if gInfoDict['draw_axes'].value == False:
                                        plt.axis('off')
                                #hides axes if draw_axes is turned off


                        '''Floor Plan Building Wall'''

                        try:

                                fbwInfo=pipe.cmd_in('python floor_building_wall '+gType,no_warn = True).splitlines()
                                #list of plotting parameter strings from tao command python floor_building_wall

                                fbwCurveList = []
                                for i in range(len(fbwInfo)):
                                        fbwCurveList.append(int(fbwInfo[i].split(';')[0]))

                                fbwCurveList = list(set(fbwCurveList)) #list of unique curve indices


                                def circle_intersection(x1,y1,x2,y2,r):
                                        '''takes centers and radius of circles, returns the 2 intersection points of overlapping circles with equal radii'''
                                        dx=x2-x1
                                        dy=y2-y1
                                        d=np.sqrt(dx**2 + dy**2)
                                        a = d/2
                                        h = np.sqrt(r**2 - a**2)
                                        xm = x1 + dx/2
                                        ym = y1 + dy/2
                                        xs1 = xm + h*dy/d
                                        xs2 = xm - h*dy/d
                                        ys1 = ym - h*dx/d
                                        ys2 = ym + h*dx/d
                                        return (xs1,ys1),(xs2,ys2)


                                bwn=pipe.cmd_in('python floor_building_wall '+gType+' name',no_warn = True).splitlines()
                                bwnTypeDict = {}
                                for i in range(len(bwn)):
                                        bwnTypeDict[bwn[i].split(';')[0]] = bwn[i].split(';')[1]
                                #dictionary where keys are wall indices and values are the corresponding building wall types

                                fps=pipe.cmd_in('python plot_shapes floor_plan',no_warn = True).splitlines()
                                fpsTypeDict = {} #building wall element types
                                fpsColorDict = {} #building wall segment colors
                                for i in range(len(fps)):
                                        fpsTypeDict[fps[i].split(';')[1].split(':')[0].lower()] = fps[i].split(';')[2].lower()
                                        if fps[i].split(';')[1].split(':')[0].lower() == 'building_wall':
                                                fpsColorDict[fps[i].split(';')[1].split(':')[2].lower()] = color(fps[i].split(';')[3].lower())
                                #dictionaries matching of the type building wall components to their color


                                for i in fbwCurveList:
                                        fbwIndexList = [] #index of point in curve
                                        fbwXList = [] #list of point x coordinates
                                        fbwYList = [] #list of point y coordinates
                                        fbwRadiusList = [] #straight line if element has 0 or missing radius
                                        for j in range(len(fbwInfo)):
                                                if i == int(fbwInfo[j].split(';')[0]):
                                                        fbwIndexList.append(int(fbwInfo[j].split(';')[1]))
                                                        fbwXList.append(float(fbwInfo[j].split(';')[2]))
                                                        fbwYList.append(float(fbwInfo[j].split(';')[3]))
                                                        try:
                                                                fbwRadiusList.append(float(fbwInfo[j].split(';')[4]))
                                                        except:
                                                                fbwRadiusList.append(0.0)
                                        k = max(fbwIndexList) #max line index

                                        while k > 1:
                                                kIndex = fbwIndexList.index(k)
                                                mIndex = fbwIndexList.index(k-1) #adjacent point to connect to

                                                if fbwRadiusList[kIndex] == 0: #draw building wall line
                                                        GraphDict['FloorPlan'].plot([fbwXList[kIndex],fbwXList[mIndex]],[fbwYList[kIndex],fbwYList[mIndex]],color=fpsColorDict[bwnTypeDict[str(i)]])

                                                else: #draw building wall arc
                                                        centerList = circle_intersection(fbwXList[mIndex],fbwYList[mIndex],fbwXList[kIndex],fbwYList[kIndex],abs(fbwRadiusList[kIndex]))
                                                        #radius and endpoints specify 2 possible circle centers for arcs
                                                        mpx = (fbwXList[mIndex] + fbwXList[kIndex])/2
                                                        mpy = (fbwYList[mIndex] + fbwYList[kIndex])/2
                                                        if np.arctan2((fbwYList[mIndex]-mpy),(fbwXList[mIndex]-mpx)) < np.arctan2(centerList[0][1],centerList[0][0]) < np.arctan2((fbwYList[mIndex]-mpy),(fbwXList[mIndex]-mpx)) and fbwRadiusList[kIndex]>0:
                                                                center = (centerList[1][0],centerList[1][1])
                                                        elif np.arctan2((fbwYList[mIndex]-mpy),(fbwXList[mIndex]-mpx)) < np.arctan2(centerList[0][1],centerList[0][0]) < np.arctan2((fbwYList[mIndex]-mpy),(fbwXList[mIndex]-mpx)) and fbwRadiusList[kIndex]<0:
                                                                center = (centerList[0][0],centerList[0][1])
                                                        elif fbwRadiusList[kIndex]>0:
                                                                center = (centerList[0][0],centerList[0][1])
                                                        else:
                                                                center = (centerList[1][0],centerList[1][1])
                                                        #find correct center

                                                        mAngle = 360 + conv*np.arctan2((fbwYList[mIndex]-center[1]),(fbwXList[mIndex]-center[0]))
                                                        kAngle = 360 + conv*np.arctan2((fbwYList[kIndex]-center[1]),(fbwXList[kIndex]-center[0]))

                                                        if abs(kAngle-mAngle) <= 180:
                                                                if kAngle > mAngle:
                                                                        t1=mAngle
                                                                        t2=kAngle
                                                                else:
                                                                        t1=kAngle
                                                                        t2=mAngle
                                                        else:
                                                                if kAngle > mAngle:
                                                                        t1=kAngle
                                                                        t2=mAngle
                                                                else:
                                                                        t1=mAngle
                                                                        t2=kAngle
                                                        #pick correct start and end angle for arc

                                                        GraphDict['FloorPlan'].add_patch(patches.Arc(center,fbwRadiusList[kIndex]*2,fbwRadiusList[kIndex]*2,theta1=t1,theta2=t2,color=fpsColorDict[bwnTypeDict[str(i)]]))
                                                        #draw building wall arc

                                                k = k - 1
                        except ValueError:
                                pass
                        #plot floor plan building walls




                        '''Floor Plan Orbit'''

                        if float(floInfoDict['floor_plan_orbit_scale'].value) != 0:
                                fpoInfo=pipe.cmd_in('python floor_orbit '+gType,no_warn = True).splitlines()

                                fpoIndexList = []
                                fpoXList = []
                                fpoYList = []
                                for i in range(len(fpoInfo)):
                                        fpoIndexList.append(int(fpoInfo[i].split(';')[1]))

                                        if fpoInfo[i].split(';')[2].lower() == 'x':
                                                for j in range(3,len(fpoInfo[i].split(';'))):
                                                        fpoXList.append(float(fpoInfo[i].split(';')[j]))
                                        if fpoInfo[i].split(';')[2].lower() == 'y':
                                                for j in range(3,len(fpoInfo[i].split(';'))):
                                                        fpoYList.append(float(fpoInfo[i].split(';')[j]))

                                GraphDict['FloorPlan'].plot(fpoXList,fpoYList,color=floInfoDict['floor_plan_orbit_color'].value.lower())
                        #Lists of floor plan orbit point indices, x coordinates, and y coordinates
                        #plot floor plan orbit




                        '''floor plan labels and axes'''

                        plt.xlabel(pgp_to_mpl(gInfoDict['x'].get_component('label')))
                        plt.ylabel(pgp_to_mpl(gInfoDict['y'].get_component('label')))
                        #plot floor plan axis labels

                        GraphDict['FloorPlan'].grid(gInfoDict['draw_grid'].value,which='major',axis='both')
                        plt.xlim(gInfoDict['x'].get_component('min'),gInfoDict['x'].get_component('max'))
                        plt.ylim(gInfoDict['y'].get_component('min'),gInfoDict['y'].get_component('max'))
                        GraphDict['FloorPlan'].set_axisbelow(True)
                        #plot floor plan grid




                '''''''''Output Data'''''''''
                #creates list of data to be returned with the figure, needed to make elements clickable on graphs

                if gInfoDict['graph^type'].value == 'lat_layout' or gInfoDict['graph^type'].value == 'floor_plan':
                        if gInfoDict['ix_universe'].value != -1:
                                gUniverse = gInfoDict[ix_universe].value

                        else:
                                gUniverse = 1

                        gBranch = gInfoDict['-1^ix_branch'].value
                        gComponent = 'model'

                else:
                        if gInfoDict['ix_universe'].value != -1:
                                gUniverse = gInfoDict[ix_universe].value

                        else:
                                gUniverse = 1

                        gBranch = gInfoDict['-1^ix_branch'].value
                        if gInfoDict['component'].value == 'model' or gInfoDict['component'].value == 'base' or gInfoDict['component'].value == 'design':
                                gComponent = gInfoDict['component'].value
                        else:
                                gComponent = 'model'
                #get universe, branch, and component

                if gInfoDict['graph^type'].value != 'floor_plan':
                        corner1=[]
                        corner2=[]
                        corner3=[]
                        corner4=[]
                        fpeIndexList = []
                        fpeShapeDict = []
                        fpeCenterDict = []
                        fpeRadiusDict = []
                if LatLayout != True:
                        eleIndexList = []
                        eleStartDict = []
                        eleEndDict = []
                        eleShapeDict = []
                        eleY1Dict = []
                #fills output list with blank lists if information does not apply to the selected graph type

                returnList = [gInfoDict['graph^type'].value, gUniverse, gBranch, gComponent, eleIndexList, eleStartDict, eleEndDict, fpeIndexList,fpeShapeDict,fpeCenterDict, fpeRadiusDict, corner1, corner2, corner3, corner4, pathDict, eleShapeDict, eleY1Dict, graph_list]
                #data to be returned with the figure to make elements clickable

                fig.tight_layout(pad=.5) #prevents graphs from overlapping



                return fig, returnList



