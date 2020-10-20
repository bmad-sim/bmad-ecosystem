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
    def __init__(self,pipe,PlotRegion):
        '''initializer, takes a tao interface and a graph region'''
        self.pipe = pipe #tao_interface object
        self.PlotRegion = PlotRegion #string describing region in tao of the desired plot

    ## @profile
    def plot(self, width):
        '''Returns a figure containing graphs using the data in the region PlotRegion of the tao instance in pipe, and plots a lat_layout below if applicable.
           Also returns information about the indices and locations of elements, width is used to modify the width of floor plan elements and the size of the lat layout.'''

        # Creates plotting figure
        fig = plt.figure()

        pipe = self.pipe # tao interface instance to plot in
        PlotRegion = self.PlotRegion # Region plot is in.

        # Determines if shapes are drawn on graph, changed to True later if needed
        LatLayout = False # draw a layout below data graphs?
        FloorPlan = False

        # Records information about element locations to be returned with the figure
        eleIndexList = []
        eleStartDict = {}
        eleEndDict = {}
        fpeIndexList = []
        fpeCenterDict = {}
        fpeRadiusDict = {}
        pathDict = {}

        def mpl_color(x):
            '''takes string containing pgplot color and returns corresponding matplotlib color'''
            x = x.lower()
            if x == 'yellow_green':
                return 'greenyellow'
            elif x == 'light_green':
                return 'limegreen'
            elif x == 'navy_blue':
                return 'navy'
            elif x == 'reddish_purple':
                return 'mediumvioletred'
            elif x == 'dark_grey':
                return 'gray'
            elif x == 'light_grey':
                return 'lightgray'
            elif x == 'Transparent' or x == 'transparent':
                return 'none'
            else:
                return x
        # End: def color

        def mpl_string(x):
            '''Takes string with pgplot characters and returns string with characters replaced with matplotlib equivalent.
               Raises NotImplementedError if an unknown pgplot character is used.'''
            x=x.replace('\\','\\\\')
            if '\\\\' in x:
                lx = '$'+x+'$'
                while lx.find('\\\\d') != -1 and lx.find('\\\\u') != -1:
                    if lx.find('\\\\d') <       lx.find('\\\\u'):
                        if lx.find('\\\\d') != -1:
                            start = lx.find('\\\\d')
                            end = lx.find('\\\\u')
                            sx = lx[start:end+3]
                            lx = lx.replace(sx,'_'+sx[3:-3])
                            sx = ''
                            start = ''
                    else:
                        if lx.find('\\\\u') != -1:
                            start = lx.find('\\\\u')
                            end = lx.find('\\\\d')
                            sx = lx[start:end+3]
                            lx = lx.replace(sx,'^'+sx[3:-3])
                            sx = ''
                            start = ''
                            end = ''

                lx = lx.replace(' ','\\ ')
                lx = lx.replace('%','\\%')

                lx = lx.replace('\\\\(2265)','\\partial')
                lx = lx.replace('\\\\ga','\\alpha')
                lx = lx.replace('\\\\gb','\\beta')
                lx = lx.replace('\\\\gg','\\gamma')
                lx = lx.replace('\\\\gd','\\delta')
                lx = lx.replace('\\\\ge','\\epsilon')
                lx = lx.replace('\\\\gz','\\zeta')
                lx = lx.replace('\\\\gy','\\eta')
                lx = lx.replace('\\\\gh','\\theta')
                lx = lx.replace('\\\\gi','\\iota')
                lx = lx.replace('\\\\gk','\\kappa')
                lx = lx.replace('\\\\gl','\\lambda')
                lx = lx.replace('\\\\gm','\\mu')
                lx = lx.replace('\\\\gn','\\nu')
                lx = lx.replace('\\\\gc','\\xi')
                lx = lx.replace('\\\\go','\\omicron')
                lx = lx.replace('\\\\gp','\\pi')
                lx = lx.replace('\\\\gr','\\rho')
                lx = lx.replace('\\\\gs','\\sigma')
                lx = lx.replace('\\\\gt','\\tau')
                lx = lx.replace('\\\\gu','\\upsilon')
                lx = lx.replace('\\\\gf','\\phi')
                lx = lx.replace('\\\\gx','\\chi')
                lx = lx.replace('\\\\gq','\\psi')
                lx = lx.replace('\\\\gw','\\omega')

                lx = lx.replace('\\\\gA','A')
                lx = lx.replace('\\\\gB','B')
                lx = lx.replace('\\\\gG','\\Gamma')
                lx = lx.replace('\\\\gD','\\Delta')
                lx = lx.replace('\\\\gE','E')
                lx = lx.replace('\\\\gZ','Z')
                lx = lx.replace('\\\\gY','H')
                lx = lx.replace('\\\\gH','\\Theta')
                lx = lx.replace('\\\\gI','I')
                lx = lx.replace('\\\\gK','\\Kappa')
                lx = lx.replace('\\\\gL','\\Lambda')
                lx = lx.replace('\\\\gM','M')
                lx = lx.replace('\\\\gN','N')
                lx = lx.replace('\\\\gC','\\Xi')
                lx = lx.replace('\\\\gO','O')
                lx = lx.replace('\\\\gP','\\Pi')
                lx = lx.replace('\\\\gR','P')
                lx = lx.replace('\\\\gS','\\Sigma')
                lx = lx.replace('\\\\gT','T')
                lx = lx.replace('\\\\gU','\\Upsilon')
                lx = lx.replace('\\\\gF','\\Phi')
                lx = lx.replace('\\\\gX','X')
                lx = lx.replace('\\\\gQ','\\Psi')
                lx = lx.replace('\\\\gW','\\Omega')
                lx = lx.replace('\\\\fn', '')
                if '\\\\' in lx:
                    raise NotImplementedError('unknown character in string, character not yet added to mpl_string: ' + repr(lx))
                return lx
            else:
                return x
        # End: def mpl_string(x)


        def circle_intersection(x1,y1,x2,y2,r):
            '''takes centers and radius of circles, returns the 2 intersection points of overlapping circles with equal radii'''
            dx = x2-x1
            dy = y2-y1
            d = np.sqrt(dx**2 + dy**2)
            a = d/2
            h = np.sqrt(r**2 - a**2)
            xm = x1 + dx/2
            ym = y1 + dy/2
            xs1 = xm + h*dy/d
            xs2 = xm - h*dy/d
            ys1 = ym - h*dx/d
            ys2 = ym + h*dx/d
            return (xs1,ys1),(xs2,ys2)


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

        #Dictionary with pgplot symbol strings as keys and corresponding matplotlib symbol strings as values
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

        # Plot Data...

        # List of plotting parameter strings from tao command python plot1
        pInfo = pipe.cmd_in('python plot1 '+PlotRegion,no_warn = True).splitlines()

        # Dictionary of Tao plot parameter name to the corresponding value.
        pInfoDict = {}
        for i in range(len(pInfo)):
            pInfoDict[pInfo[i].split(';')[0]] = str_to_tao_param(pInfo[i])

        # List of graph names and heights
        gNameList = []
        gHeightsList = []
        for i in range(pInfoDict['num_graphs'].value):
            gNameList.append(pInfoDict[('graph['+str(i+1)+']')].value)
            gHeightsList.append(1) #makes graphs the same size
        number_graphs = len(gNameList) + 1

        # Relative height of lat layout to graph if there is one graph
        layout_height = 0.2 * width

        # Makes layout approximately a fixed size as number of graphs changes
        gHeightsList.append(layout_height*len(gNameList) + layout_height*(len(gNameList)-1))

        gs = fig.add_gridspec(nrows = number_graphs,ncols = 1, height_ratios = gHeightsList)

        gSubPlotList = [None]  # List of graphs with the graph's matplotlib subplot as it's value. Indexed from 1.
        gFullNameList = [] # List of graph names, eg: r1.g or r3.x

        for gNum in range(1, len(gNameList)+1):

            # Create plots in figure, second line also makes x axes scale together
            if gNum == 1:
                gSubPlotList.append(fig.add_subplot(gs[gNum-1,0]))
            elif gNum > 1:
                gSubPlotList.append(fig.add_subplot(gs[gNum-1,0],sharex = gSubPlotList[1]))

            # Graph Data...

            # Graph "region.graph" full name, EG: "r13.g" or "top.x"
            gFullName = PlotRegion + '.' + gNameList[gNum-1]
            gFullNameList.append(gFullName)

            # List of graph parameters from tao command python plot_graph
            gInfo = pipe.cmd_in('python plot_graph '+gFullName,no_warn = True).splitlines()

            # Dictionary of graph parameters.
            gInfoDict = {}
            for i in range(len(gInfo)):
                gInfoDict[gInfo[i].split(';')[0]] = str_to_tao_param(gInfo[i])

            #------------------
            # Curve Data...

            # List of curve names
            cList = []
            for i in range(gInfoDict['num_curves'].value):
                cList.append(gInfoDict[('curve['+str(i+1)+']')].value)

            # List of curve parameters.
            cInfo = []
            for i in cList:
                cInfo.append(pipe.cmd_in('python plot_curve '+gFullName+'.'+i,no_warn = True).splitlines())

            cInfoDictList = []
            for i in range(len(cList)):
                cInfoDict = {}
                for j in range(len(cInfo[i])):
                    cInfoDict[cInfo[i][j].split(';')[0]] = str_to_tao_param(cInfo[i][j])
                cInfoDictList.append(cInfoDict)
                cInfoDict = {}


            #------------------
            # Line Data...

            # List of points from tao command python plot_line for each curve
            lInfo = []
            try:
                for i in cList:
                    lInfo.append(pipe.cmd_in('python plot_line '+gFullName+'.'+i,no_warn = True).splitlines())

                # List of lists of points used to draw each curve
                if len(lInfo) != 0:
                    PointsSuperList = []
                    for i in range(len(cList)):
                        pList = []
                        for j in range(len(lInfo[i])):
                            pList.append(lInfo[i][j].split(';'))
                        LineCoords = []
                        for j in range(len(pList)):
                            LineCoords.append([float(pList[j][1]),float(pList[j][2])])
                        PointsSuperList.append(LineCoords)
                        LineCoords = []
            except: #handle graph with no lines
                lInfo = []
                PointsSuperList = []
                for i in cList:
                    PointsSuperList.append([[0,0]])


            #------------------
            # Symbol Data

            # List  of points from tao command python plot_symbol for each curve
            try:
                sInfo = []
                for i in cList:
                    sInfo.append(pipe.cmd_in('python plot_symbol '+gFullName+'.'+i,no_warn = True).splitlines())
                SymbolSuperList = []
                # List of points used to draw symbols on each curve
                for i in range(len(cList)):
                    sList = []
                    for j in range(len(sInfo[i])):
                        sList.append(sInfo[i][j].split(';'))
                    SymCoords = []
                    for j in range(len(sList)):
                        SymCoords.append([float(sList[j][2]),float(sList[j][3])])
                    SymbolSuperList.append(SymCoords)
                    SymCoords = []
            except: # Handle graph with no symbols
                sInfo = []
                SymbolSuperList = []
                for i in cList:
                    SymbolSuperList.append([[0,0]])

            #------------------
            # Histogram Data...

            try:
                hInfo = []
                for i in cList:
                    hInfo.append(pipe.cmd_in('python plot_histogram '+gFullName+'.'+i,no_warn = True).splitlines())
                hInfoDictList = []
                # List of lists of dictionaries of plot_histogram data for each curve
                for i in range(len(cList)):
                    hInfoDict = {}
                    for j in range(len(hInfo[i])):
                        hInfoDict[hInfo[i][j].split(';')[0]] = str_to_tao_param(hInfo[i][j])
                    hInfoDictList.append(hInfoDict)
                    hInfoDict = {}

            except:
                hInfo = []

            # Plot Data

            # List of data needed to plot line and symbol graphs
            # Includes points, and line and symbol information for each curve
            CurvesList = []
            for i in range(len(cList)):
                CurveData = []
                CurveData.append(PointsSuperList[i]) #points for each curve
                CurveData.append(SymbolSuperList[i]) #symbols for each curve
                CurveData.append(mpl_color(cInfoDictList[i]['line'].get_component('color'))) #line color
                CurveData.append(StylesDict[cInfoDictList[i]['line'].get_component('pattern').lower()]) #line style
                if (cInfoDictList[i]['draw_line'].value == True): #line width if drawn
                    CurveData.append(cInfoDictList[i]['line'].get_component('width'))
                else:
                    CurveData.append(0)
                CurveData.append(mpl_color(cInfoDictList[i]['symbol'].get_component('color'))) #symbol color

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

            #------------------
            # Plotting...

            # Plots line graphs, histograms, phase space plots, dynamic aperture graphs, and plots for wave analysis
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

                if gInfoDict['graph^type'].value == 'data' or gInfoDict['graph^type'].value == 'wave.0' or \
                   gInfoDict['graph^type'].value == 'wave.a' or gInfoDict['graph^type'].value == 'wave.b':
                    LineList.append(gSubPlotList[gNum].plot(xpList,ypList,color = i[2],linestyle = i[3],linewidth = i[4]/2))
                    gSubPlotList[gNum].plot(xsList,ysList,color = i[5],linewidth = 0,markerfacecolor = i[6],markersize = i[7]/2,marker = i[8],mew = i[9]/2)

                    # Add lat_layout if x-axis is "s" (longitudinal postion).
                    if pInfoDict['x_axis_type'].value == 's':
                        LatLayout = True

                    # Wave region boundaries
                    if gInfoDict['graph^type'].value != 'data': #wave analysis rectangles
                        wInfo = pipe.cmd_in('python wave params',no_warn = True).splitlines()
                        a1 = float(wInfo[1].split(';')[3])
                        a2 = float(wInfo[2].split(';')[3])
                        b1 = float(wInfo[3].split(';')[3])
                        b2 = float(wInfo[4].split(';')[3])

                    if i[5] == 'blue' or i[5] == 'navy' or i[5] == 'cyan' or i[5] == 'green' or i[5] == 'purple':
                        waveColor = 'orange'
                    else:
                        waveColor = 'blue'
                    #wave analysis rectangle color

                    if gInfoDict['graph^type'].value == 'wave.0' or gInfoDict['graph^type'].value == 'wave.a':
                        gSubPlotList[gNum].add_patch(patches.Rectangle((a1,yMin),a2-a1,yMax-yMin,fill = False,color = waveColor))
                    if gInfoDict['graph^type'].value == 'wave.0' or gInfoDict['graph^type'].value == 'wave.b':
                        gSubPlotList[gNum].add_patch(patches.Rectangle((b1,yMin),b2-b1,yMax-yMin,fill = False,color = waveColor))
                #line and symbol graphs

                elif gInfoDict['graph^type'].value == 'dynamic_aperture':
                    LineList.append(gSubPlotList[gNum].plot(xpList,ypList,color = i[2],linestyle = i[3],linewidth = i[4]/2))
                    gSubPlotList[gNum].plot(xsList,ysList,color = i[5],linewidth = 0,markerfacecolor = i[6],markersize = i[7]/2,marker = i[8],mew = i[9]/2)
                #dynamic aperture graphs

                elif gInfoDict['graph^type'].value == 'phase_space':
                    if lInfo != []:
                        LineList.append(gSubPlotList[gNum].plot(xpList,ypList,color = i[2],linestyle = i[3],linewidth = i[4]/2))
                        gSubPlotList[gNum].plot(xsList,ysList,color = i[5],linewidth = 0,markerfacecolor = i[6],markersize = i[7]/2,marker = i[8],mew = i[9]/2)
                    else:
                        LineList.append(gSubPlotList[gNum].plot(xsList,ysList,color = i[5],linewidth = 0,markerfacecolor = i[6],markersize = i[7]/2,marker = i[8],mew = i[9]/2))
                #phase space graphs

                elif gInfoDict['graph^type'].value == 'histogram':
                    LineList.append(gSubPlotList[gNum].hist(xpList,bins = int(hInfoDictList[CurvesList.index(i)]['number'].value),weights = ypList,histtype = 'step',color = i[5]))
                #histogram

            if gInfoDict['graph^type'].value == 'key_table':
                raise NotImplementedError('key table is not available in the GUI')

            if gInfoDict['graph^type'].value == 'lat_layout':
                LatLayout = True
                gNameList = []
                plt.axis('off')
            #sets up lat layout plot

            if gInfoDict['graph^type'].value == 'floor_plan':
                FloorPlan = True
                gNameList = []
                plt.axis('off')
            #sets up floor plan plot

            if gInfoDict['draw_axes'].value == False:
                plt.axis('off')
            #hides axes if draw_axes is turned off

            if gInfoDict['why_invalid'].value != '':
                raise ValueError(gInfoDict['why_invalid'].value + ', make sure graph is properly initialized')

            #plots line and symbol graphs and histograms, lat layouts and floor plans are drawn later
            #LineList gives names of curves




            plt.title(mpl_string(gInfoDict['title'].value)+' '+gInfoDict['title_suffix'].value)
            #plot title

            LegendList = [] #legends for each graph
            LabelList = [] #labels in each legend
            try:
                for i in range(len(CurvesList)):
                    LegendList.append(LineList[i][0])
                    if mpl_string(cInfoDictList[i]['legend_text'].value) != '':
                        LabelList.append(mpl_string(cInfoDictList[i]['legend_text'].value))
                    elif mpl_string(cInfoDictList[i]['data_type'].value) == 'physical_aperture':
                        LabelList.append(mpl_string(cInfoDictList[i]['data_type'].value))
                    else:
                        LabelList.append('')
                #list of curves to be added to a legend and list of labels for each curve in the legend
            except IndexError:
                raise NotImplementedError('unknown graph type')


            if (gInfoDict['draw_curve_legend'].value == True and LabelList != ['']) and gInfoDict['graph^type'].value != 'lat_layout' and gInfoDict['graph^type'].value != 'floor_plan':
                gSubPlotList[gNum].legend(LegendList,LabelList)
            #plot legend

            plt.xlabel(mpl_string(gInfoDict['x'].get_component('label')))
            plt.ylabel(mpl_string(gInfoDict['y'].get_component('label')))
            #plot axis labels

            gSubPlotList[gNum].grid(gInfoDict['draw_grid'].value,which = 'major',axis = 'both')
            #plot grid

            plt.xlim(gInfoDict['x'].get_component('min'),gInfoDict['x'].get_component('max'))
            plt.ylim(gInfoDict['y'].get_component('min'),gInfoDict['y'].get_component('max'))
            #set axis limits

            gSubPlotList[gNum].set_axisbelow(True)
            #place graphs over grid lines

        #-----------------------
        # Plots lat layouts

        if LatLayout == True:
            if gInfoDict['graph^type'].value != 'lat_layout': #add space for lat layout below graph
                latLayoutSubPlot = fig.add_subplot(gs[len(gNameList),0],sharex = gSubPlotList[1])

            else: # Standalone lat layout graph
                latLayoutSubPlot = fig.add_subplot(len(gNameList)+1,1,len(gNameList)+1,sharex = gSubPlotList[1])

            # List of parameter strings from tao command python plot_graph
            layInfo = pipe.cmd_in('python plot_graph layout.g',no_warn = True).splitlines()

            # Tao_parameter object names from python plot_graph
            # Dictionary of tao_parameter name string keys to the corresponding tao_parameter object
            layInfoDict = {}
            for i in range(len(layInfo)):
                layInfoDict[layInfo[i].split(';')[0]] = str_to_tao_param(layInfo[i])



            # Makes lat layout only have horizontal axis for panning and zooming
            twinAxes = latLayoutSubPlot.axes.twinx()
            plt.xlim(gInfoDict['x'].get_component('min'),gInfoDict['x'].get_component('max'))
            plt.ylim(layInfoDict['y'].get_component('min'),layInfoDict['y'].get_component('max'))
            twinAxes.set_navigate(True)
            latLayoutSubPlot.axis('off')
            twinAxes.axis('off')

            # Sets axis limits and creates second axis to allow x panning and zooming
            latLayoutSubPlot.axes.set_navigate(False)
            latLayoutSubPlot.axhline(y = 0,xmin = 1.1*layInfoDict['x'].get_component('min'),xmax = 1.1*layInfoDict['x'].get_component('max'),color = 'Black')


            # Lat layout branch and universe information
            if layInfoDict['ix_universe'].value != -1:
                universe = layInfoDict[ix_universe].value
            else:
                universe = 1
            branch = layInfoDict['-1^ix_branch'].value


            # List of strings containing information about each element
            eleInfo = pipe.cmd_in('python plot_lat_layout '+str(universe)+'@'+str(branch),no_warn = True).splitlines()

            # All dict keys and entries are strings which match a lattice layout element index (eg: '1') string to the corresponding information
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
                eleStartDict[eleInfo[i].split(';')[0]] = float(eleInfo[i].split(';')[1])
                eleEndDict[eleInfo[i].split(';')[0]] = float(eleInfo[i].split(';')[2])
                eleLwDict[eleInfo[i].split(';')[0]] = float(eleInfo[i].split(';')[3])
                eleShapeDict[eleInfo[i].split(';')[0]] = eleInfo[i].split(';')[4].lower()
                eleY1Dict[eleInfo[i].split(';')[0]] = float(eleInfo[i].split(';')[5])
                eleY2Dict[eleInfo[i].split(';')[0]] = float(eleInfo[i].split(';')[6])
                eleColorDict[eleInfo[i].split(';')[0]] = mpl_color(eleInfo[i].split(';')[7].lower())
                eleNameDict[eleInfo[i].split(';')[0]] = eleInfo[i].split(';')[8]

            # Plotting line segments one-by-one can be slow if there are thousands of lattice elements.
            # So keep a list of line segments and plot all at once at the end.

            y_max = max(max(eleY1Dict.values()), max(eleY2Dict.values()))
            y2_floor = -max(eleY2Dict.values())  # Note negative sign
            lines = []
            widths = []
            colors = []

            for i in eleIndexList:
                s1 = eleStartDict[str(i)]
                s2 = eleEndDict[str(i)]
                y1 = eleY1Dict[str(i)]
                y2 = -eleY2Dict[str(i)]  # Note negative sign.
                wid = eleLwDict[str(i)]
                color = eleColorDict[str(i)]
                shape = eleShapeDict[str(i)]
                name = eleNameDict[str(i)]

                try:
                    # Normal case where element is not wrapped around ends of lattice.
                    if s2-s1 > 0:

                        # Draw box element
                        if shape == 'box':
                            latLayoutSubPlot.add_patch(patches.Rectangle((s1,y1), s2-s1, y2-y1, lw = wid,color = color,fill = False))

                        # Draw xbox element
                        elif shape == 'xbox':
                            latLayoutSubPlot.add_patch(patches.Rectangle((s1,y1), s2-s1, y2-y1, lw = wid,color = color,fill = False))
                            lines.extend([ [(s1,y1), (s2,y2)], [(s1,y2), (s2,y1)] ])
                            colors.extend([color, color])
                            widths.extend([wid, wid])

                        # Draw x element
                        elif shape == 'x':
                            lines.extend([ [(s1,y1), (s2,y2)], [(s1,y2), (s2,y1)] ])
                            colors.extend([color, color])
                            widths.extend([wid, wid])

                        # Draw bow_tie element
                        elif shape == 'bow_tie':
                            lines.extend([ [(s1,y1), (s2,y2)], [(s1,y2), (s2,y1)], [(s1,y1), (s1,y2)], [(s2,y1), (s2,y2)] ])
                            colors.extend([color, color, color, color])
                            widths.extend([wid, wid, wid, wid])

                        # Draw rbow_tie element
                        elif shape == 'rbow_tie':
                            lines.extend([ [(s1,y1), (s2,y2)], [(s1,y2), (s2,y1)], [(s1,y1), (s2,y1)], [(s1,y2), (s2,y2)] ])
                            colors.extend([color, color, color, color])
                            widths.extend([wid, wid, wid, wid])

                        # Draw diamond element
                        elif shape == 'diamond':
                            s_mid = (s1 + s2) / 2
                            lines.extend([ [(s1,0), (s_mid, y1)], [(s1,0), (s_mid, y2)], [(s2,0), (s_mid, y1)], [(s2,0), (s_mid, y2)] ])
                            colors.extend([color, color, color, color])
                            widths.extend([wid, wid, wid, wid])

                        # Draw circle element
                        elif shape == 'circle':
                            s_mid = (s1 + s2) / 2
                            latLayoutSubPlot.add_patch(patches.Ellipse((s_mid,0), y1-y2, y1-y2, lw = wid,color = color,fill = False))

                        # Draw element name
                        latLayoutSubPlot.text((s1+s2)/2, 1.1*y2_floor, name,ha = 'center',va = 'top',clip_on = True,color = color)

                    # Case where element is wrapped round the lattice ends.
                    else:
                        s_min = layInfoDict['x'].get_component('min')
                        s_max = layInfoDict['x'].get_component('max')

                        # Draw wrapped box element
                        if shape == 'box':
                            latLayoutSubPlot.plot([s1,s_max],[y1,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s_max],[y2,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[y1,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[y2,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s1],[y1,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s2,s2],[y1,y2],lw = wid,color = color)

                        # Draw wrapped xbox element
                        elif shape == 'xbox':
                            latLayoutSubPlot.plot([s1,s_max],[y1,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s_max],[y2,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s_max],[y1,0],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s_max],[y2,0],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[y1,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[y2,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[0,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[0,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s1],[y1,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s2,s2],[y1,y2],lw = wid,color = color)

                        # Draw wrapped x element
                        elif shape == 'x':
                            latLayoutSubPlot.plot([s1,s_max],[y1,0],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s_max],[y2,0],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[0,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[0,y2],lw = wid,color = color)

                        # Draw wrapped bow tie element
                        elif shape == 'bow_tie':
                            latLayoutSubPlot.plot([s1,s_max],[y1,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s_max],[y2,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s_max],[y1,0],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s_max],[y2,0],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[y1,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[y2,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[0,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[0,y2],lw = wid,color = color)

                        # Draw wrapped diamond element
                        elif shape == 'diamond':
                            latLayoutSubPlot.plot([s1,s_max],[0,y1],lw = wid,color = color)
                            latLayoutSubPlot.plot([s1,s_max],[0,y2],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[y1,0],lw = wid,color = color)
                            latLayoutSubPlot.plot([s_min,s2],[y2,0],lw = wid,color = color)

                        # Draw wrapped element name
                        latLayoutSubPlot.text(s_max, 1.1*y2_floor, name,ha = 'right',va = 'top',clip_on = True,color = color)
                        latLayoutSubPlot.text(s_min, 1.1*y2_floor, name,ha = 'left',va = 'top',clip_on = True,color = color)

                except KeyError:
                    pass

            # Draw all line segments
            latLayoutSubPlot.add_collection(mp.collections.LineCollection(lines, colors = colors, linewidths = widths))

            # Invisible line to give the lat layout enough vertical space.
            # Without this, the tops and bottoms of shapes could be cut off
            latLayoutSubPlot.plot([0, 0],[-1.7*y_max,1.3*y_max],lw = wid,color = color,alpha = 0)

        else:
            latLayoutSubPlot = fig.add_subplot(gs[len(gNameList),0],sharex = gSubPlotList[1])
            latLayoutSubPlot.remove()


        # Plots floor plans

        if FloorPlan == True:
            gSubPlotForFloorPlan = fig.add_subplot(len(gNameList)+1,1,len(gNameList)+1,sharex = gSubPlotList[1])

            floInfo = pipe.cmd_in('python plot_graph '+gFullName,no_warn = True).splitlines()
            #list of plotting parameter strings from tao command python plot_graph


            floInfoDict = {}
            for i in range(len(floInfo)):
                floInfoDict[floInfo[i].split(';')[0]] = str_to_tao_param(floInfo[i])
            #tao_parameter object names from python plot_graph for a floor plan
            #dictionary of tao_parameter name string keys to the corresponding tao_parameter object

            if floInfoDict['ix_universe'].value != -1:
                universe = floInfoDict[ix_universe].value

            else:
                universe = 1

            fpeInfo = pipe.cmd_in('python floor_plan '+gFullName,no_warn = True).splitlines()
            #list of plotting parameter strings from tao command python floor_plan


            fpeIndexList = [] # Contains lists of branch index then element index for each point
            fpeTypeDict = {}
            fpeSxDict = {} # Start x coordinate
            fpeSyDict = {} # Start y coordinate
            fpeSaDict = {} # Start angle
            fpeExDict = {} # End x coordinate
            fpeEyDict = {} # End y coordinate
            fpeEaDict = {} # End angle
            fpeLwDict = {} # Line width
            fpeShapeDict = {}
            fpeY1Dict = {} # Distance above
            fpeY2Dict = {} # Distance below
            fpeColorDict = {}
            fpeNameDict = {}
            fpeAlDict = {}  # Arc length
            fpeBaDict = {}  # Bend angle
            fpeSfaDict = {} # Relative angle of starting face to incoming line
            fpeEfaDict = {} # Relative angle of ending face to incoming line
            corner1 = {}    # Corner coordinates of objects
            corner2 = {}
            corner3 = {}
            corner4 = {}
            for i in range(len(fpeInfo)):
                fpeIndexList.append([int(fpeInfo[i].split(';')[0]),int(fpeInfo[i].split(';')[1])])
                fpeTypeDict[str(fpeIndexList[i])] = fpeInfo[i].split(';')[2].lower()
                fpeSxDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[3])
                fpeSyDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[4])
                fpeSaDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[5])
                fpeExDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[6])
                fpeEyDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[7])
                fpeEaDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[8])
                fpeLwDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[9])
                fpeShapeDict[str(fpeIndexList[i])] = fpeInfo[i].split(';')[10].lower()
                fpeY1Dict[str(fpeIndexList[i])] = width*float(fpeInfo[i].split(';')[11])
                fpeY2Dict[str(fpeIndexList[i])] = width*float(fpeInfo[i].split(';')[12])
                fpeColorDict[str(fpeIndexList[i])] = mpl_color(fpeInfo[i].split(';')[13].lower())
                fpeNameDict[str(fpeIndexList[i])] = fpeInfo[i].split(';')[14]
                try:
                    fpeAlDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[15])
                    fpeBaDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[16])
                    fpeSfaDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[17])
                    fpeEfaDict[str(fpeIndexList[i])] = float(fpeInfo[i].split(';')[18])
                except IndexError:
                    fpeAlDict[str(fpeIndexList[i])] = 0
                    fpeBaDict[str(fpeIndexList[i])] = 0
                    fpeSfaDict[str(fpeIndexList[i])] = 0
                    fpeEfaDict[str(fpeIndexList[i])] = 0

            #dict keys and entries are strings which match a floor plan element index (eg: '1') to the corresponding information

            conv = (180)/(np.pi) #radian to degree conversion
            for i in fpeIndexList:
                fpeCenterDict[str(i)] = ([fpeSxDict[str(i)] + (fpeExDict[str(i)]-fpeSxDict[str(i)])/2,fpeSyDict[str(i)] + (fpeEyDict[str(i)]-fpeSyDict[str(i)])/2])
                fpeRadiusDict[str(i)] = fpeY1Dict[str(i)] #for click detection
                ele_key = fpeTypeDict[str(i)]
                x1 = fpeSxDict[str(i)]
                x2 = fpeExDict[str(i)]
                y1 = fpeSyDict[str(i)]
                y2 = fpeEyDict[str(i)]
                off1 = fpeY1Dict[str(i)]
                off2 = fpeY2Dict[str(i)]
                width = fpeLwDict[str(i)]
                color = fpeColorDict[str(i)]
                angStart = fpeSaDict[str(i)]
                angEnd = fpeEaDict[str(i)]
                shape = fpeShapeDict[str(i)]
                ele_name = fpeNameDict[str(i)]
                relAngStart = fpeSfaDict[str(i)]
                relAngEnd = fpeEfaDict[str(i)]

                try:
                    if ele_key == 'drift' or ele_key == 'kicker':
                        gSubPlotForFloorPlan.plot([x1,x2],[y1,y2],color = 'black')
                    #draw drift element

                    if off1 == 0 and off2 == 0 and ele_key != 'sbend' and color != '':
                        gSubPlotForFloorPlan.plot([x1,x2],[y1,y2],lw = width,color = color)
                    #draw line element

                    elif shape == 'box' and ele_key != 'sbend' and color != '':
                        gSubPlotForFloorPlan.add_patch(patches.Rectangle((x1 + off2*np.sin(angStart), y1 - off2*np.cos(angStart)),
                                                       np.sqrt((x2-x1)**2 + (y2-y1)**2), off1+off2,
                                                       lw = width,color = color,fill = False,angle = angStart*conv))
                    #draw box element

                    elif shape == 'xbox' and ele_key != 'sbend' and color != '':
                        gSubPlotForFloorPlan.add_patch(patches.Rectangle((x1 + off2*np.sin(angStart), y1 - off2*np.cos(angStart)),
                                                                         np.sqrt((x2-x1)**2 + (y2-y1)**2), off1+off2,
                                                                         lw = width,color = color,fill = False,angle = angStart*conv))
                        gSubPlotForFloorPlan.plot([x1 + off2*np.sin(angStart),x2 - off1*np.sin(angStart)],[y1 - off2*np.cos(angStart),y2 + off1*np.cos(angStart)],lw = width,color = color)
                        gSubPlotForFloorPlan.plot([x1 - off1*np.sin(angStart),x2 + off2*np.sin(angStart)],[y1 + off1*np.cos(angStart),y2 - off2*np.cos(angStart)],lw = width,color = color)
                    #draw xbox element

                    elif shape == 'x' and ele_key != 'sbend' and color != '':
                        gSubPlotForFloorPlan.plot([x1 + off2*np.sin(angStart),x2 - off1*np.sin(angStart)],[y1 - off2*np.cos(angStart),y2 + off1*np.cos(angStart)],lw = width,color = color)
                        gSubPlotForFloorPlan.plot([x1 - off1*np.sin(angStart),x2 + off2*np.sin(angStart)],[y1 + off1*np.cos(angStart),y2 - off2*np.cos(angStart)],lw = width,color = color)
                    #draw x element


                    elif shape == 'bow_tie' and ele_key != 'sbend' and color != '':
                        gSubPlotForFloorPlan.plot([x1 + off2*np.sin(angStart),x2 - off1*np.sin(angStart)],[y1 - off2*np.cos(angStart),y2 + off1*np.cos(angStart)],lw = width,color = color)
                        gSubPlotForFloorPlan.plot([x1 - off1*np.sin(angStart),x2 + off2*np.sin(angStart)],[y1 + off1*np.cos(angStart),y2 - off2*np.cos(angStart)],lw = width,color = color)
                        gSubPlotForFloorPlan.plot([x1 - off1*np.sin(angStart),x2 - off1*np.sin(angStart)],[y1 + off1*np.cos(angStart),y2 + off1*np.cos(angStart)],lw = width,color = color)
                        gSubPlotForFloorPlan.plot([x1 + off2*np.sin(angStart),x2 + off2*np.sin(angStart)],[y1 - off2*np.cos(angStart),y2 - off2*np.cos(angStart)],lw = width,color = color)
                    #draw bow_tie element


                    elif shape == 'diamond' and ele_key != 'sbend' and color != '':
                        gSubPlotForFloorPlan.plot([x1,x1 + (x2-x1)/2 - off1*np.sin(angStart)],[y1,y1 + (y2-y1)/2 + off1*np.cos(angStart)],lw = width,color = color)
                        gSubPlotForFloorPlan.plot([x1 + (x2-x1)/2 - off1*np.sin(angStart),x2],[y1 + (y2-y1)/2 + off1*np.cos(angStart),y2],lw = width,color = color)
                        gSubPlotForFloorPlan.plot([x1,x1 + (x2-x1)/2 + off2*np.sin(angStart)],[y1,y1 + (y2-y1)/2 - off2*np.cos(angStart)],lw = width,color = color)
                        gSubPlotForFloorPlan.plot([x1 + (x2-x1)/2 + off2*np.sin(angStart),x2],[y1 + (y2-y1)/2 - off2*np.cos(angStart),y2],lw = width,color = color)
                    #draw diamond element


                    elif shape == 'circle' and ele_key != 'sbend' and color != '':
                        gSubPlotForFloorPlan.add_patch(patches.Circle((x1 + (x2-x1)/2,y1 + (y2-y1)/2),off1,lw = width,color = color,fill = False))
                    #draw circle element


                    elif shape == 'box' and ele_key == 'sbend' and color != '':
                        gSubPlotForFloorPlan.plot([x1-off1*np.sin(angStart-relAngStart), x1+off2*np.sin(angStart-relAngStart)],
                                                  [y1+off1*np.cos(angStart-relAngStart), y1-off2*np.cos(angStart-relAngStart)], lw = width, color = color)
                        gSubPlotForFloorPlan.plot([x2-off1*np.sin(angEnd+relAngEnd), x2+off2*np.sin(angEnd+relAngEnd)],
                                                  [y2+off1*np.cos(angEnd+relAngEnd), y2-off2*np.cos(angEnd+relAngEnd)], lw = width, color = color)
                        #draws straight sbend edges

                        intersection = intersect(line([x1-off1*np.sin(angStart), y1+off1*np.cos(angStart)], [x1+off2*np.sin(angStart), y1-off2*np.cos(angStart)]),
                                                 line([x2-off1*np.sin(angEnd), y2+off1*np.cos(angEnd)], [x2+off2*np.sin(angEnd), y2-off2*np.cos(angEnd+relAngEnd)]))
                        #center of circle used to draw arc edges of sbends

                        if intersection == False:
                            gSubPlotForFloorPlan.plot([x1-off1*np.sin(angStart-relAngStart), x2-off1*np.sin(angEnd+relAngEnd)],
                                                      [y1+off1*np.cos(angStart-relAngStart), y2+off1*np.cos(angEnd+relAngEnd)], lw = width, color = color)
                            gSubPlotForFloorPlan.plot([x1+off2*np.sin(angStart-relAngStart), x2+off2*np.sin(angEnd+relAngEnd)],
                                                      [y1-off2*np.cos(angStart-relAngStart), y2-off2*np.cos(angEnd+relAngEnd)], lw = width, color = color)
                        #draw sbend edges if bend angle is 0


                        elif intersection != False:
                            angle1 = 360 + conv*np.arctan2(y1+off1*np.cos(angStart-relAngStart)-intersection[1],x1-off1*np.sin(angStart-relAngStart)-intersection[0])
                            angle2 = 360 + conv*np.arctan2(y2+off1*np.cos(angEnd+relAngEnd)-intersection[1],x2-off1*np.sin(angEnd+relAngEnd)-intersection[0])
                            #angles of further curve endpoints relative to center of circle
                            angle3 = 360 + conv*np.arctan2(y1-off2*np.cos(angStart-relAngStart)-intersection[1],x1+off2*np.sin(angStart-relAngStart)-intersection[0])
                            angle4 = 360 + conv*np.arctan2(y2-off2*np.cos(angEnd+relAngEnd)-intersection[1],x2+off2*np.sin(angEnd+relAngEnd)-intersection[0])
                            #angles of closer curve endpoints relative to center of circle

                            if abs(angle1-angle2)<180:
                                a1 = min(angle1, angle2)
                                a2 = max(angle1, angle2)
                            else:
                                a1 = max(angle1, angle2)
                                a2 = min(angle1, angle2)

                            if abs(angle3-angle4)<180:
                                a3 = min(angle3, angle4)
                                a4 = max(angle3, angle4)
                            else:
                                a3 = max(angle3, angle4)
                                a4 = min(angle3, angle4)
                            #determines correct start and end angles for arcs

                            gSubPlotForFloorPlan.add_patch(patches.Arc((intersection[0],intersection[1]),
                                                            np.sqrt((x1-off1*np.sin(angStart-relAngStart)-intersection[0])**2 + (y1+off1*np.cos(angStart-relAngStart)-intersection[1])**2)*2,
                                                            np.sqrt((x1-off1*np.sin(angStart-relAngStart)-intersection[0])**2 + (y1+off1*np.cos(angStart-relAngStart)-intersection[1])**2)*2,
                                                            theta1 = a1,theta2 = a2,lw = width,color = color))
                            gSubPlotForFloorPlan.add_patch(patches.Arc((intersection[0],intersection[1]),
                                                            np.sqrt((x1+off2*np.sin(angStart-relAngStart)-intersection[0])**2 + (y1-off2*np.cos(angStart-relAngStart)-intersection[1])**2)*2,
                                                            np.sqrt((x1+off2*np.sin(angStart-relAngStart)-intersection[0])**2 + (y1-off2*np.cos(angStart-relAngStart)-intersection[1])**2)*2,
                                                            theta1 = a3,theta2 = a4,lw = width,color = color))
                            #draw sbend edges if bend angle is nonzero
                    #draw sbend element


                    if ele_name != '' and color != '' and np.sin(((angEnd+angStart)/2)) > 0:
                        gSubPlotForFloorPlan.text(x1+(x2-x1)/2 - 1.3*off1*np.sin(angStart), y1+(y2-y1)/2 + 1.3*off1*np.cos(angStart),
                                             ele_name, ha = 'right',va = 'center',color = 'black', rotation = -90+((angEnd+angStart)/2)*conv,clip_on = True,rotation_mode = 'anchor')

                    elif ele_name != '' and color != '' and np.sin(((angEnd+angStart)/2)) <= 0:
                        gSubPlotForFloorPlan.text(x1+(x2-x1)/2 - 1.3*off1*np.sin(angStart),y1+(y2-y1)/2 + 1.3*off1*np.cos(angStart),
                                             ele_name,ha = 'left',va = 'center',color = 'black',rotation = 90+((angEnd+angStart)/2)*conv,clip_on = True,rotation_mode = 'anchor')
                    #draw element name




                    '''floor plan click detection'''

                    if ele_key == 'sbend' and intersection != False: #for sbend click detection
                        c1 = [x1-off1*np.sin(angStart-relAngStart),y1+off1*np.cos(angStart-relAngStart)]
                        c2 = [x2-off1*np.sin(angEnd+relAngEnd),y2+off1*np.cos(angEnd+relAngEnd)]
                        c3 = [x1+off2*np.sin(angStart-relAngStart),y1-off2*np.cos(angStart-relAngStart)]
                        c4 = [x2+off2*np.sin(angEnd+relAngEnd),y2-off2*np.cos(angEnd+relAngEnd)]
                        #corners of sbend

                        if angStart > angEnd:
                            outerRadius = np.sqrt((x1-off1*np.sin(angStart-relAngStart)-intersection[0])**2 + (y1+off1*np.cos(angStart-relAngStart)-intersection[1])**2)
                            innerRadius = np.sqrt((x1+off2*np.sin(angStart-relAngStart)-intersection[0])**2 + (y1-off2*np.cos(angStart-relAngStart)-intersection[1])**2)
                        else:
                            outerRadius = -np.sqrt((x1-off1*np.sin(angStart-relAngStart)-intersection[0])**2 + (y1+off1*np.cos(angStart-relAngStart)-intersection[1])**2)
                            innerRadius = -np.sqrt((x1+off2*np.sin(angStart-relAngStart)-intersection[0])**2 + (y1-off2*np.cos(angStart-relAngStart)-intersection[1])**2)
                        #radii of sbend arc edges



                        middleAngle = (angStart+angEnd)/2

                        top = [intersection[0]-outerRadius*np.sin(middleAngle),intersection[1]+outerRadius*np.cos(middleAngle)]
                        bottom = [intersection[0]-innerRadius*np.sin(middleAngle),intersection[1]+innerRadius*np.cos(middleAngle)]
                        #midpoints of top and bottom arcs in an sbend

                        topCP = [2*(top[0])-.5*(c1[0])-.5*(c2[0]),2*(top[1])-.5*(c1[1])-.5*(c2[1])]
                        bottomCP = [2*(bottom[0])-.5*(c3[0])-.5*(c4[0]),2*(bottom[1])-.5*(c3[1])-.5*(c4[1])]
                        #corresponding control points for a quadratic Bezier curve that passes through the corners and arc midpoint

                        verts = [c1,topCP,c2,c4,bottomCP,c3,c1]
                        codes = [Path.MOVETO,Path.CURVE3,Path.CURVE3,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CLOSEPOLY]
                        pathDict[str(i)] = Path(verts,codes)

                        '''patch = patches.PathPatch(Path(verts,codes),facecolor = 'green',alpha = .5)
                        gSubPlotForFloorPlan.add_patch(patch)'''
                        #visualize clickable regions
                    #path approximating sbend region for clickable region on graph using lines and quadratic Bezier curves

                    else: #for non sbend click detection
                        corner1[str(i)] = [x1 - off1*np.sin(angStart),y1 + off1*np.cos(angStart)]
                        corner2[str(i)] = [x2 - off1*np.sin(angStart),y2 + off1*np.cos(angStart)]
                        corner3[str(i)] = [x1 + off2*np.sin(angStart),y1 - off2*np.cos(angStart)]
                        corner4[str(i)] = [x2 + off2*np.sin(angStart),y2 - off2*np.cos(angStart)]
                    #coordinates of corners of a floor plan element for clickable region


                except KeyError:
                    pass

                if gInfoDict['draw_axes'].value == False:
                    plt.axis('off')
                #hides axes if draw_axes is turned off


            '''Floor Plan Building Wall'''

            try:

                fbwInfo = pipe.cmd_in('python building_wall_graph '+gFullName,no_warn = True).splitlines()
                #list of plotting parameter strings from tao command python floor_building_wall

                fbwCurveList = []
                for i in range(len(fbwInfo)):
                    fbwCurveList.append(int(fbwInfo[i].split(';')[0]))

                fbwCurveList = list(set(fbwCurveList)) #list of unique curve indices

                bwn = pipe.cmd_in('python building_wall_list',no_warn = True).splitlines()
                bwnTypeDict = {}
                for i in range(len(bwn)):
                    bwnTypeDict[bwn[i].split(';')[0]] = bwn[i].split(';')[1]
                #dictionary where keys are wall indices and values are the corresponding building wall types

                fps = pipe.cmd_in('python shape_list floor_plan',no_warn = True).splitlines()
                fpsTypeDict = {} #building wall element types
                fpsColorDict = {} #building wall segment colors
                for i in range(len(fps)):
                    fpsTypeDict[fps[i].split(';')[1].split(':')[0].lower()] = fps[i].split(';')[2].lower()
                    if fps[i].split(';')[1].split(':')[0].lower() == 'building_wall':
                        fpsColorDict[fps[i].split(';')[1].split(':')[2].lower()] = mpl_color(fps[i].split(';')[3].lower())
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
                        if bwnTypeDict[str(i)] not in fpsColorDict.keys():
                            # TODO: This is a temporary fix to deal with building wall segments
                            # that don't have an associated floor_plan shape
                            # Currently this will fail to match to wild cards
                            # in the shape name (e.g. building_wall::* should match
                            # to every building wall segment, but currently it
                            # matches to none).  A more sophisticated way of getting the
                            # floor_plan shape settings for building walls will be required
                            # in the future, either through a python command in tao or
                            # with a method on the python to match wild cards to wall segment names
                            print("No floor_plan shape defined for building_wall segment " + bwnTypeDict[str(i)])
                            k -= 1
                            continue

                        if fbwRadiusList[kIndex] == 0: #draw building wall line
                            gSubPlotForFloorPlan.plot([fbwXList[kIndex],fbwXList[mIndex]],[fbwYList[kIndex],fbwYList[mIndex]],color = fpsColorDict[bwnTypeDict[str(i)]])

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
                                    t1 = mAngle
                                    t2 = kAngle
                                else:
                                    t1 = kAngle
                                    t2 = mAngle
                            else:
                                if kAngle > mAngle:
                                    t1 = kAngle
                                    t2 = mAngle
                                else:
                                    t1 = mAngle
                                    t2 = kAngle
                            #pick correct start and end angle for arc

                            gSubPlotForFloorPlan.add_patch(patches.Arc(center,fbwRadiusList[kIndex]*2,fbwRadiusList[kIndex]*2,theta1 = t1,theta2 = t2,color = fpsColorDict[bwnTypeDict[str(i)]]))
                            #draw building wall arc

                        k = k - 1
            except ValueError:
                pass
            #plot floor plan building walls




            '''Floor Plan Orbit'''

            if float(floInfoDict['floor_plan_orbit_scale'].value) != 0:
                fpoInfo = pipe.cmd_in('python floor_orbit '+gFullName,no_warn = True).splitlines()

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

                gSubPlotForFloorPlan.plot(fpoXList,fpoYList,color = floInfoDict['floor_plan_orbit_color'].value.lower())
            #Lists of floor plan orbit point indices, x coordinates, and y coordinates
            #plot floor plan orbit




            '''floor plan labels and axes'''

            plt.xlabel(mpl_string(gInfoDict['x'].get_component('label')))
            plt.ylabel(mpl_string(gInfoDict['y'].get_component('label')))
            #plot floor plan axis labels

            gSubPlotForFloorPlan.grid(gInfoDict['draw_grid'].value,which = 'major',axis = 'both')
            plt.xlim(gInfoDict['x'].get_component('min'),gInfoDict['x'].get_component('max'))
            plt.ylim(gInfoDict['y'].get_component('min'),gInfoDict['y'].get_component('max'))
            gSubPlotForFloorPlan.set_axisbelow(True)
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
            corner1 = []
            corner2 = []
            corner3 = []
            corner4 = []
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

        returnList = [gInfoDict['graph^type'].value, gUniverse, gBranch, gComponent, eleIndexList, eleStartDict, eleEndDict, fpeIndexList,fpeShapeDict,fpeCenterDict, fpeRadiusDict, corner1, corner2, corner3, corner4, pathDict, eleShapeDict, eleY1Dict, gFullNameList]
        #data to be returned with the figure to make elements clickable

        fig.tight_layout(pad = .5) #prevents graphs from overlapping

        return fig, returnList
