#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import re
from PyQt4 import QtGui
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt


class NavigationToolbar2(NavigationToolbar):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('Save', 'Home')]


def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out


class Window(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.reflections = []
        self.initUI()
        self.showMaximized()
        self.readFcoFile()


    def initUI(self):
        self.fcoLineRegExp = re.compile('^(\s+\d+){3}')
        self.newLayout = QtGui.QVBoxLayout(self)
        self.main_widget = QtGui.QWidget(self)

        exitAction = QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)

        openAction = QtGui.QAction(QtGui.QIcon('open.png'), '&Open File', self)
        openAction.setShortcut('Ctrl+O')
        openAction.setStatusTip('Open .fco or .fcf File')
        openAction.triggered.connect(lambda: self.readFcoFile())

        drkAction = QtGui.QAction('&DRKplot', self)
        drkAction.setObjectName('drk')
        drkAction.setShortcut('Ctrl+D')
        drkAction.setStatusTip('Make DRKplot')
        drkAction.triggered.connect(lambda: self.plotThis(drkAction.objectName(),self.drkPlot(0.05)))

        drk2Action = QtGui.QAction('&DRKplot equal binning', self)
        drk2Action.setObjectName('drk2')
        drk2Action.setShortcut('Ctrl+E')
        drk2Action.setStatusTip('Make DRKplot equal binning')
        drk2Action.triggered.connect(lambda: self.plotThis(drk2Action.objectName(),self.drkPlot2()))

        lnFoOverFCAction = QtGui.QAction('&ln(Fo^2/Fc^2)+1', self)
        lnFoOverFCAction.setObjectName('lnover')
        lnFoOverFCAction.setShortcut('Ctrl+L')
        lnFoOverFCAction.setStatusTip('Make ln Plot')
        lnFoOverFCAction.triggered.connect(lambda: self.plotThis(lnFoOverFCAction.objectName(),self.lnFoOverFCPlot()))

        FoOverFCAction = QtGui.QAction('&Fo^2/Fc^2', self)
        FoOverFCAction.setObjectName('FoOverFc')
        FoOverFCAction.setShortcut('Ctrl+F')
        FoOverFCAction.setStatusTip('Make Fo/Fc Plot')
        FoOverFCAction.triggered.connect(lambda: self.plotThis(FoOverFCAction.objectName(),self.FoOverFCPlot()))

        normPlotAction = QtGui.QAction('&Normal Propability Plot', self)
        normPlotAction.setObjectName('normplot')
        normPlotAction.setShortcut('Ctrl+N')
        normPlotAction.setStatusTip('Make Normal Propability Plot')
        normPlotAction.triggered.connect(lambda: self.plotThis(normPlotAction.objectName(),self.normplot()))

        menubar = self.menuBar()
        self.addToolBarBreak()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openAction)
        fileMenu.addAction(exitAction)
        plotMenu = menubar.addMenu('&Plot')
        plotMenu.addAction(drkAction)
        plotMenu.addAction(drk2Action)
        plotMenu.addAction(lnFoOverFCAction)
        plotMenu.addAction(FoOverFCAction)
        plotMenu.addAction(normPlotAction)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.statusBar().showMessage('Ready')
        self.addToolBarBreak()
#        self.matplotlibToolbar = NavigationToolbar(self.canvas, self)
#        self.matplotlibToolbar = self.figure.canvas.toolbar
        self.setWindowTitle('DRKPlot')
        self.newLayout.addWidget(self.canvas)
        self.cursor = QtGui.QTextCursor()
        self.main_widget.setLayout(self.newLayout)
        self.setCentralWidget(self.main_widget)
        self.show()

    def readFcoFile(self):
        self.reflections = []
        file = QtGui.QFileDialog.getOpenFileName(self)
        print file
        self.fcoFileContent = open(file).readlines()
        self.maxRes = 0.
        for i in self.fcoFileContent:
            if re.match(self.fcoLineRegExp,i):
                if not (float(i.split()[3])==0.0 or float(i.split()[4])==0.0):
                    if float(i.split()[6]) > self.maxRes:
                        self.maxRes = float(i.split()[6])
                    self.reflections.append([int(i.split()[0]),int(i.split()[1]),int(i.split()[2]),float(i.split()[3]),float(i.split()[4]),float(i.split()[5]),float(i.split()[6])])
        self.reflections = sorted(self.reflections,key=lambda x: x[6])
        self.plotThis('drk',self.drkPlot(0.05))

    def pickReflection(self,event):
        ind = event.ind
        self.statusBar().showMessage('H: {} K: {} L: {} F(calc)^2: {} F(obs)^2: {} Sigma(F^2): {} sin(theta)/lambda: {} F(obs)^2/F(calc)^2: {}'.format(self.reflections[ind[0]][0],self.reflections[ind[0]][1],self.reflections[ind[0]][2],self.reflections[ind[0]][3],self.reflections[ind[0]][4],self.reflections[ind[0]][5],self.reflections[ind[0]][6],float(self.reflections[ind[0]][4])/float(self.reflections[ind[0]][3])))

    def pickBin(self,event):
        ind = event.ind
        self.statusBar().showMessage('sum Fo^2 / sum Fc^2: {} sin(theta)/lambda: {} Number of Reflections in Bin: {}'.format(self.ydata[ind],self.xdata[ind],len(self.toAverageFc[self.xdata[ind]])))

    def plotThis(self,whichPlot,(x,y),textFieldContent=0.0):
        try:
            self.canvas.mpl_disconnect(self.cid)
        except AttributeError:
            pass
        try:
            self.removeToolBar(self.toolbar)
        except AttributeError:
            pass
        try:
            self.removeToolBar(self.toolbar2)
        except AttributeError:
            pass
        try:
            self.matplotlibToolbar.deleteLater()
            self.matplotlibToolbar = None
        except AttributeError:
            pass

        self.ax = self.figure.add_subplot(111)
        self.ax.hold(False)
        if whichPlot == 'normplot':
            self.matplotlibToolbar = NavigationToolbar2(self.canvas, self)
            self.newLayout.addWidget(self.matplotlibToolbar)
            #print x
            osm = stats.probplot(x)  #,plot=self.ax
            #
            # output of normal probability plot for testing purposes
            #
            #testOutput = open('normplot_test.log','w')
            #for i in x:
            #    testOutput.write('{}\n'.format(i))
            #testOutput.close()
            self.ax.plot([xdiag for xdiag in range(-20, 20)], [ydiag for ydiag in range(-20, 20)], 'b-',osm[0][0], osm[0][1], 'ro', linewidth=2,markersize=6, markevery=0.01, zorder = 10)  #
            self.ax.axhline(color='k', lw=1, zorder=2)
            self.ax.axvline(color='k', lw=1, zorder=2)
            self.ax.set_ylim(-4, 4)
            self.ax.set_xlim(-4, 4)
            self.ax.set_xlabel(r'$\mathit{Expected DR}$')
            self.ax.set_ylabel(r'$\mathit{Experimental DR}$')
        elif whichPlot == 'drk2':
            self.matplotlibToolbar = NavigationToolbar(self.canvas, self)
            self.newLayout.addWidget(self.matplotlibToolbar)
            self.toolbar2 = self.addToolBar('DRKPLOT2')
            self.toolbar2.setMovable(False)
            self.validator = QtGui.QIntValidator()
            self.numberLabel = QtGui.QLabel('Number of Bins')
            self.binNumber = QtGui.QLineEdit()
            self.binNumber.setText('400')
            self.binNumber.setFixedWidth(40)
            self.binNumber.setValidator(self.validator)
            self.setNumber = QtGui.QPushButton('Set')
            self.setNumber.clicked.connect(lambda: self.plotThis('drk2',self.drkPlot2(int(self.binNumber.text()))))
            self.toolbar2.addWidget(self.numberLabel)
            self.toolbar2.addWidget(self.binNumber)
            self.toolbar2.addWidget(self.setNumber)
            self.ax.plot(self.xdata, self.ydata, 'ro', [x for x in range(-1, 3)], [1 for x in range(-1, 3)], 'k-', picker = 5, zorder = 10)
            self.cid = self.canvas.mpl_connect('pick_event', self.pickBin)
            self.ax.set_xlim([-0.05, self.maxRes + 0.05])
            self.ax.set_ylim([0.8, 1.2])
            self.ax.set_xlabel(r'$\mathit{sin(\theta)/\lambda}$')
            self.ax.set_ylabel(r'$\mathit{\sum F(obs)^2/\sum F(calc)^2}$')
        elif whichPlot == 'drk':
            self.matplotlibToolbar = NavigationToolbar(self.canvas, self)
            self.newLayout.addWidget(self.matplotlibToolbar)
            self.toolbar = self.addToolBar('DRKPLOT')
            self.toolbar.setMovable(False)
            self.validator = QtGui.QDoubleValidator()
            self.incrementLabel = QtGui.QLabel('Stepsize for Bins')
            self.binIncrement = QtGui.QLineEdit()
            self.binIncrement.setText('0.05')
            self.binIncrement.setFixedWidth(40)
            self.binIncrement.setValidator(self.validator)
            self.setBins = QtGui.QPushButton('Set')
            self.setBins.clicked.connect(lambda: self.plotThis(('drk'),self.drkPlot(float(self.binIncrement.text()))))
            self.toolbar.addWidget(self.incrementLabel)
            self.toolbar.addWidget(self.binIncrement)
            self.toolbar.addWidget(self.setBins)
            self.ax.plot(self.xdata,self.ydata, 'ro', picker = 5, zorder=10)
            self.cid = self.canvas.mpl_connect('pick_event',self.pickBin)
            self.ax.set_xlim([-0.05,self.maxRes+0.05])
            self.ax.set_ylim([0.8,1.2])
            self.ax.axhline(0.95,0.0,color='r', lw=1,alpha=0.25, zorder=2)
            self.ax.axhline(1.05,0.0,color='r', lw=1,alpha=0.25, zorder=2)
            self.ax.axhline(1.0,0.0,color='k', lw=1,alpha=0.5, zorder=2)
            self.ax.set_xlabel(r'$\mathit{sin(\theta)/\lambda}$')
            self.ax.set_ylabel(r'$\mathit{\sum F(obs)^2/\sum F(calc)^2}$')
        else:
            self.matplotlibToolbar = NavigationToolbar(self.canvas, self)
            self.newLayout.addWidget(self.matplotlibToolbar)
            self.ax.hold(False)
            self.ax.plot(x,y, 'ro',  picker = 5, zorder = 1)
            self.ax.axhline(1.0,0.0,color='k', lw=1,alpha=0.5, zorder=2)
            self.ax.set_xlim([0.0,self.maxRes+0.05])
            self.cid = self.canvas.mpl_connect('pick_event',self.pickReflection)
            self.ax.set_xlabel(r'$\mathit{sin(\theta)/\lambda}$')
            if whichPlot == 'lnover':
                self.ax.set_ylabel(r'$\mathit{ln(F(obs)^2/F(calc)^2)+1}$')
            else:
                self.ax.set_ylabel(r'$\mathit{F(obs)^2/F(calc)^2}$')
        self.canvas.draw()
        self.ax.grid()
        self.ax.set_axisbelow(True)

    def normplot(self):
        self.distribution = []
        for i in self.reflections:
            if i[0] == 1 and i[1] == 0 and i[2] == 0:
                print i
                print (i[3]-i[4])/i[5]
            self.distribution.append((i[3]-i[4])/i[5])
        return self.distribution,[]

    def drkPlot2(self,nBins=400):
        test = chunkIt(self.reflections,nBins)
        self.toAverageFo = {}
        self.toAverageFc = {}
        for i in test:
            mean = np.mean([x[6] for x in i])
            for j in i:
                try:
                    self.toAverageFo[mean].append(j[4])
                    self.toAverageFc[mean].append(j[3])
                except KeyError:
                    self.toAverageFo[mean] = [j[4]]
                    self.toAverageFc[mean] = [j[3]]
        self.xdata = []
        self.ydata = []
        for i in self.toAverageFc.keys():
            self.ydata.append(np.sum(self.toAverageFo[i])/np.sum(self.toAverageFc[i]))
            self.xdata.append(i)
        return self.xdata,self.ydata

    def drkPlot(self, Increment):
        self.toAverageFo = {}
        self.toAverageFc = {}
        for i in self.reflections:
            for j in np.arange(0.0, self.maxRes, Increment):
                if i[6] <= float(j) and i[6] > float(j-Increment):
                    try:
                        self.toAverageFo[float(j)].append(i[4])
                        self.toAverageFc[float(j)].append(i[3])
                    except KeyError:
                        self.toAverageFo[float(j)] = [i[4]]
                        self.toAverageFc[float(j)] = [i[3]]
        self.xdata = []
        self.ydata = []
        for i in self.toAverageFc.keys():
            self.ydata.append(np.sum(self.toAverageFo[i])/np.sum(self.toAverageFc[i]))
            self.xdata.append(float(i))
        return self.xdata,self.ydata

    def FoOverFCPlot(self):
        xdata = []
        ydata = []
        for i in self.reflections:
            ydata.append(i[4]/i[3])
            xdata.append(i[6])
        print len(xdata),len(ydata)
        return xdata,ydata

    def lnFoOverFCPlot(self):
        xdata = []
        ydata = []
        for i in self.reflections:
            ydata.append(np.log(i[4]/i[3])+1)
            xdata.append(i[6])
        return xdata,ydata



if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    main = Window()
    #main.showMaximized()

    sys.exit(app.exec_())
