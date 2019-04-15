#!/usr/bin/env python
''' Interactive Image Viewer/Analyzer of streamed images.
 
Features:
+ Input source: 
    * file, 
    * HTTP link to image, 
    * RHIC ADO parameter (requires cns.py), 
    * EPICS PV (requires epics.py).
+ Wide range of image formats.
+ 16-bit/channel images supported.
+ Image orientation and rotation (program options: -o and -R).
+ Interactive zooming, panning, rotation.
+ Contrast control: Displays histogram of image data with movabled region defining the dark/light levels.
+ ROI and embedded plot for measuring image values.
+ Isocurves. The isocurve level defines the threshold for spot finding.
+ Fast multi-spot finder, reports and logs centroid position and integral of most intense spots in the ROI.
+ Gaussian fit of the brightest spot (can be disabled using -T option)
+ Export as PNG,TIFF, JPG..., SVG?, Matplotlib, CSV, HDF5.
+ Interactive python console with access to image data, graphics objects and shell commands (program option: -c).
+ Configuration and reporting in the parameter dock.
+ Ref. Images: save/retrieve image to/from a reference slots.
+ Binary operation on current image and a reference: addition, subtraction.
+ Background subtraction using a reference image.
+ User add-ons (-u option)
+ Fast browsing/cleanup of the image directories
'''
__version__ = 'v174 2019-04-13'# forked from imageViewer.py
import sys
import traceback
import time
import datetime
import struct 
import subprocess
import os
import threading
from collections import OrderedDict
from json import dumps
    
#from PyQt4 import QtGui, QtCore
from pyqtgraph.Qt import QtGui, QtCore
Gray_color_table = [QtGui.qRgb(i, i, i) for i in range(256)]
import pyqtgraph as pg
import pyqtgraph.dockarea
import pyqtgraph.console
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem
from pyqtgraph.parametertree.parameterTypes import (
    WidgetParameterItem, registerParameterType)

# Interpret image data as row-major instead of col-major
pg.setConfigOptions(imageAxisOrder='row-major')

import numpy as np
#import skimage.transform as st
from scipy import ndimage
import math

# if graphics is done in callback, then we need this:
QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)
# otherwise errors: [xcb] Unknown request in queue while dequeuing

#````````````````````````````Necessary explicit globals```````````````````````
pargs = None
imager = None
Prec = 3 # precision for logging and display
# define signal on data arrival
EventPocessingFinished = threading.Event()
#v93#cprintLock = threading.Lock()
X,Y = 0,1
#Addon = None
#````````````````````````````Helper Functions`````````````````````````````````        
#````````````````````````````Stuff for profiling``````````````````````````````
from timeit import default_timer as timer
profilingState = OrderedDict() # keeps processing times for diagnostics

def profile(state):
    # store the state
    #global profilingState
    profilingState[state] = timer()

def profStates(first,last):
    # returns text lines with time differences between intermediate states
    txt = ''
    l = timer()
    t = 0
    for key,value in profilingState.items():
        if key == first: 
            t = value
        elif key == last:
            break
        if t:
            d = value - t    
            txt += 'time of '+key+' :%0.3g'%(d)+'\n'
            t = value
    return txt

def profDif(first,last):
    return '%0.3g'%(profilingState[last] - profilingState[first])
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
def printi(msg): print('info: '+msg)
    
def printw(msg): print('WARNING: '+msg)
    
def printe(msg): print('ERROR: '+msg)

def printd(msg): 
    if pargs.dbg: print('dbg: '+msg)

gWidgetConsole = None
def cprint(msg):
    """Print info on console dock"""
    if gWidgetConsole:
        gWidgetConsole.write('#'+msg+'\n') # use it to inform the user

def cprinte(msg): # print to Console
    """Print error message on console dock"""
    if gWidgetConsole:
        gWidgetConsole.write('#ERROR: '+msg+'\n') # use it to inform the user
    printe(msg)

def cprintw(msg): # print to Console
    """Print warning message on console dock"""
    if gWidgetConsole:
        gWidgetConsole.write('#WARNING: '+msg+'\n') # use it to inform the user
    printw(msg)

def font():
    font = pg.QtGui.QFont()
    font.setPointSize(10)
    font.setBold(True)
    return font

def mainSpotText():
    o = pg.TextItem('',color=(0,150,255),anchor=(1,0),fill='w')
    o.setFont(font())
    return o

def pointerText():
    o = pg.TextItem('',color=(0,150,255),anchor=(0,0),fill='w')
    o.setFont(font())
    return o

def checkPath(path):
    """Check if path exists, if not, create directory"""
    try:
        if not os.path.exists(path):
            print('checkPath created new path:',path)
            os.makedirs(path)
    except Exception as e:
        cprinte('in checkPath '+path+' error: '+str(e))

def file_idx(fileList,positionPercent):
    l = len(fileList)
    ipos = int(l*positionPercent/100.)
    try:
        fn = fileList[ipos].rsplit('/',1)[1]
    except:
        fn = fileList[ipos]
    cprint('file at %.0f%%'%(100.*ipos/l)+': [%d]'%ipos+' of %d: '%l+fn)
    return ipos
    
def qMessage(text):
    """Message dialog in separate window"""
    ans = QtGui.QMessageBox.question(None, 'Confirm', text, 
      QtGui.QMessageBox.Yes, defaultButton = QtGui.QMessageBox.No)
    return 1 if ans == QtGui.QMessageBox.Yes else 0

def rgb2gray(data,graysum=True):
    """Convert RGB to Grayscale"""
    #st = timer()
    if len(data.shape) < 3:
        return data
    #'''
    r = data[:,:,0].astype('uint16')
    g = data[:,:,1].astype('uint16')
    b = data[:,:,2].astype('uint16')
    if graysum: # uniform sum
        #rc = (r/3.+ g/3. + b/3.).astype('uint16') #12ms
        rc = ((r + g + b)/3).astype('uint8') #3ms,  4 times faster
        #print('rc:',type(r),rc[100][100])
    else: # using perception-based weighted sum 
        rc = (0.2989 * r + 0.5870 * g + 0.1140 * b).astype('uint8')
    #'''
    #rc = np.mean(data,axis=2).astype('uint8') #this is slightly slower, check on large images
    #print('rgb2gray',timer() - st)
    return rc

def crop_width(array):
    '''crop 2d array horizontally to have the width divisible by 4
    this is important for PNG encoding'''
    wi = array.shape[1]
    wo = wi & -4
    if wi != wo:
        array = array[:,:wo,...]
    return array

def rotate(data,degree):
    if pargs.flip: degree = -degree
    degree += pargs.rotate
    fracp,n = math.modf(degree/90)
    #s = timer()
    if fracp == 0:
        data = np.rot90(data,n) # fast rotate by integral of 90 degree
    else: 
        #if degree: data = st.rotate(data, degree, resize = True, preserve_range = True)
        if degree: data = ndimage.rotate(data, degree, reshape = True, order=1)
    #printi('rotate time:'+str(timer()-s))
    if pargs.flip:
        if   pargs.flip == 'V': return data[::-1,...]
        elif pargs.flip == 'H': return data[:,::-1,...]
    return data

def blur(a,width=2):
    if len(a.shape) == 2:
        return ndimage.gaussian_filter(a,(width,width)) # 10 times faster than pg.gaussianFilter
    else:
        cprintw('blurring of color images is not implemented yet')
        return a       

def ellipse_prob05(sigmaX,sigmaY,pcc):
    """Parameters of the probability 0.5 ellipse
    (Contour ellipse with sum = 0.5 spot sum)"""
    ds = sigmaX**2 - sigmaY**2
    theta = math.atan(2.*pcc*sigmaX*sigmaY/ds)/2. if abs(ds) > 0.001 else math.pi/2.
    cos2 = math.cos(theta)**2
    sin2 = 1. - cos2
    cs = 2.*cos2 - 1. # cos2 - sin2
    sigmaX2 = sigmaX**2
    sigmaY2 = sigmaY**2
    sigmaU = math.sqrt((cos2*sigmaX2 - sin2*sigmaY2)/cs)
    sigmaV = math.sqrt((cos2*sigmaY2 - sin2*sigmaX2)/cs)
    return theta,sigmaU,sigmaV

#````````````````````````````Spot processing stuff````````````````````````````
def centroid(data):
    """Returns first (mean) an second (sigma) central moments and 
    Person Correlation Coefficient of 2D array data, 
    Caution on PCC calculation: depending on the numbers involved, it can 
    sometimes be numerically unstable."""
    s = np.zeros(2) # sum of samples along axis
    ss = np.zeros(2) # sum of squares of samples along axis
    sxy = 0. #sum of X[i]*Y[i]
    iax = [0.]*2 # index vectors along axis
    n = 0 # number of samples
    for axis in (X,Y):
        idx = range(data.shape[axis])
        iax[axis] = np.array(idx,dtype=float)
        oppositeAxis = int(not axis)
        projection = data.sum(axis = oppositeAxis).astype(float)
        s[axis] = np.dot(projection,iax[axis])
        ss[axis] = np.dot(projection,iax[axis]**2)
        if axis == 1:
            n += sum(projection)
            
            # get sum(xy) for correlation coeff.
            # note, that section is most time-consuming
            for i in idx:
                sxdot = i*np.dot(data[:,i],iax[oppositeAxis])
                sxy += sxdot
    
    sigman = np.sqrt(n*ss - s**2)
    means = s/n
    sigmas = sigman/n
    if min(sigman) != 0.:
        pcc = (n*sxy - s[0]*s[1]) / (sigman[0]*sigman[1])
    else: pcc = 0.
    return means, sigmas, pcc, n

FinalFit = True
RankBkg = 1 # Rank = 3 is troublesome for wide peaks, it favors the quadratic background 
if FinalFit:
    from scipy.optimize import curve_fit
    def peak_shape(xx, halfWidth):
        """Function, representing the peak shape,
        It should be in a performance-optimized form.
        The halfWidth = sqrt(2)*sigma.
        """
        try: r = np.exp(-0.5*(xx/halfWidth)**2) 
        except: r = np.zeros(len(x))
        return r
    def func_sum_of_peaks(xx, *par):
        # Fitting function: base and sum of peaks.
        #print('fsop:',par)
        s = np.zeros(len(xx)) + par[0] # if RankBkg = 1
        #s = par[0] + par[1]*xx + par[2]*xx**2 # if RankBkg = 3
        for i in range(RankBkg,len(par),3):
            s += par[i+2]*peak_shape(xx-par[i],par[i+1])
        return s    

def find_spots(region, threshold, maxSpots, fitBrightest = False,
               fitBaseBounds = (-np.inf,+np.inf)):
    """find up to maxSpots in the ndarray region and return its centroids and sums
    """
    rh,rw =  region.shape
    #print('>find_spots:','hw:',(rh,rw),(region[0,0],region[0,rw-1]),'sum:%.1f'%region.sum()+', mean:%.1f'%region.mean()+', std:%.1f'%region.std())
    profile('startFind')
    
    above_threshold = np.copy(region) 
    #above_threshold[above_threshold < threshold] = 0
    
    # Subtract the threshold, otherwise, if the pedestal is high,
    # we will get the sigma = ~width/sqrt(12)
    above_threshold = region - threshold
    above_threshold[above_threshold < 0.] = 0
    
    profile('thresholding') # 5% of processing time spent here
    
    # now find the objects
    labeled, number_of_objects = ndimage.label(above_threshold)
    if number_of_objects == 0:
        #print('no objects:',region.shape,np.max(region),np.sum(region))
        return []
    profile('labeling') # 12% of processing time spent here
    
    # sort the labels according to their sums
    sums = ndimage.sum(above_threshold,labeled,index=range(1,number_of_objects+1))
    sumsSorted = sorted(enumerate(sums),key=lambda idx: idx[1],reverse=True)
    labelIndexesSortedBySum = [i[0] for i in sumsSorted]
    profile('sums') #TODO: 30% of processing time spent here
    peak_slices = ndimage.find_objects(labeled)
    profile('find') #TODO: 15% of processing time spent here
    
    # calculate centroids
    centroids = []
    for spotIdxx,lbl in enumerate(labelIndexesSortedBySum[:maxSpots]):
        islice = peak_slices[lbl]
        ofsXY = islice[1].start,islice[0].start
        only_labeled = np.copy(above_threshold[islice])
        only_labeled[labeled[islice] != lbl+1] = 0 # zero all not belonging to label lbl+1
        p,w,pcc,pixSum = centroid(only_labeled.T)
        #pos = [ofsXY[0]+p[0],ofsXY[1]+p[1]] # the following is nicer:
        pos = p + ofsXY

        fittedPars = []
        fitRange = []
        #print('spotIdxx,w,fb',spotIdxx,w,fitBrightest)
        #too_small = w[0]*w[1] < 0.1
        too_small = min(w) < 1.
        #print('w',w,('','too small')[too_small])
        if too_small:
            continue
        if spotIdxx==0 and not too_small and fitBrightest:
            # correct position and width by fitting each projection with gaussian and a baseline
            for axis in (X,Y):
                projection = region.sum(axis).astype(float)
                lArr = len(projection)
                nP = region.shape[axis] # number of points summed i each column
                xyLetter = 'XY'[axis]
                #print('axis'+xyLetter+' l:%d'%lArr+',w:%.1f'%w[axis])
                posMax = pos

                seqx = range(lArr)
                iPosMax = int(round(posMax[axis]))
                
                # we have good estimation on position and width
                # correct the width using gaussian fit with linear base in 6*width range
                # IMPORTANT! using less that 6*width may cause underestimated baseline
                halfRange = int(round(w[axis]))*6
                il,ir =  max(0,iPosMax - halfRange), min(lArr,iPosMax+halfRange)-1
                #print('il,ir',il,ir,iPosMax,halfRange)
                if ir - il < 6:
                    printw('Not enough area to calculate background for '+xyLetter\
                      +'need 6, got %d'%(ir - il))
                    fittedPars.append([])
                    break
                fitRange.append((il,ir))
                
                # estimate the amplitude
                amp = max(projection)
                ilm,irm = iPosMax-2, iPosMax+2
                ## quadratic fit around top
                #r = np.polyfit(seqx[ilm:irm],projection[ilm:irm],2)
                #print('pftime',timer() - ts) # it may take longer than the gaussian iPosMax
                #print('pfit amp:%.1f'%amp+', [%d]'%iPosMax+'=%.1f'%projection[iPosMax],r)
                #amp = r[0]*iPosMax**2 + r[1]*iPosMax + r[2]
                # simple mean is good enough
                amp = projection[ilm:irm].mean()
                
                #base = threshold*region.shape[axis]
                #base = 5.*region.shape[axis]
                base = min(projection[il],projection[ir])               
                guess = [base,posMax[axis],w[axis],amp-base] # if  RankBkg = 1

                # bounds significantly slows the performance, from 7ms to 20ms
                gl = len(guess) - 1
                bounds = ([fitBaseBounds[0]*nP]+[-np.inf]*gl,
                  [fitBaseBounds[1]*nP]+[+np.inf]*gl)
                #bounds = ([-np.inf]+[-np.inf]*gl, [+np.inf]+[+np.inf]*gl)
                #print('bounds',bounds)
                
                xr = seqx[il:ir]
                yr = projection[il:ir]
                #print(il,projection[il:ir])
                #ts = timer()profile('find')
                try:
                    fp,pcov = curve_fit(func_sum_of_peaks,xr,yr,guess,bounds=bounds)#,method='dogbox')#, factor = 1.)#,bounds=bounds)
                    stdevPars = np.sqrt(np.diag(pcov))/abs(fp)
                    #print('stdevPars',stdevPars)
                    sigmaStdev = stdevPars[RankBkg+1]
                    if sigmaStdev > 0.2 or math.isnan(sigmaStdev):
                        printw('fitted sigma'+xyLetter+' is bad, stdev:%.4g'%stdevPars[RankBkg+1])
                        fittedPars.append([])
                        continue
                    #v152#fp[RankBkg+1] /= 1.41421356 # sigma is width/sqrt(2)
                    #print('curfit  time',timer()-ts) # ~2ms for 100 points
                    #print('fitted',fp)
                    pos[axis] = fp[RankBkg]#+ofsXY[axis]
                    w[axis] = fp[RankBkg+1]
                except Exception as e:
                    msg = 'Fit'+xyLetter+' failed: '+str(e)
                    printe(msg)
                    cprint('WARNING '+msg)
                    #fittedPars.append([])
                    break
                fittedPars.append(fp)
        centroids.append((pos, w, pcc, pixSum, fittedPars, fitRange))
    profile('centroids') # 15% of processing time spent here, 09/10:30%
    return centroids 
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#````````````````````````````Graphs outside of the main widget```````````````
class GraphRoiProj():
    def __init__(self,axis,pxPmm=10.,maxPix=640):
        self.oppositeAxis = {X:Y,Y:X}[axis]
        self.axis = axis
        self.widget = None
        self.pxPmm = pxPmm
        self.maxPix = maxPix
        #self.bottomAxis.setScale(abs(1./self.pxPmm))
        
    def update(self,xyArray,ofs=(0,0),peakPars=[]):
        l = xyArray.shape[:2][self.oppositeAxis]
        a = xyArray.mean(axis=self.axis)
        n = a.sum()
        offs = ofs[self.axis]
        x = np.arange(l+1,dtype=float)
        
        xmin = imager.pix2mm(offs,offs)[self.axis]
        xmax = imager.pix2mm(offs+l,offs+l)[self.axis]
        xPoints = np.linspace(xmin,xmax,l+1)        
        if self.widget is not None:
            self.widget.plot(xPoints,a,stepMode=True,pen='k',clear=True)
            self.win = self.widget.win
        else:
            self.widget = pg.plot(xPoints,a,stepMode=True,pen='k',clear=True)
            self.widget.showGrid(True,True)
            self.hv = {X:'Horizontal',Y:'Vertical'}
            self.widget.win.setWindowTitle(self.hv[self.axis]+' ROI Means')
            self.widget.win.resize(480,320)
            #self.bottomAxis = self.widget.getAxis('bottom')           
            #self.bottomAxis.setScale(abs(1./self.pxPmm))
            self.widget.setLabel('bottom','Position (mm)')
            self.widget.setLabel('left','Average intensity')
        s = np.dot(a,x[:-1])
        ss = np.dot(a,x[:-1]**2)
        mean = s/n + offs
        std = np.sqrt(n*ss - s**2)/n
        #self.widget.setLabel('top','Mean,std: %.f, %.f pix, '%(mean,std)\
        #  +'%.2f, %.3f mm'%(xf*(mean-self.maxPix/2)/abs(self.pxPmm),
        #  std/abs(self.pxPmm)))
        #self.widget.setLabel(label)
        if len(peakPars) == 0:
            self.widget.setLabel('top','Final fit is disabled')
            return
        pp = peakPars[self.axis]
        #print('pp',self.axis,peakPars)
        if len(pp) == 0:
            self.widget.setLabel('top','No '+self.hv[self.axis]+' fit data')
            return
        else:
            posPix = pp[RankBkg+0] + offs
            pos = imager.pix2mm(posPix,posPix)[self.axis]
            sigPix = abs(pp[RankBkg+1])
            sig = imager.mm_per_pixel(sigPix,sigPix)[self.axis]
            self.widget.setLabel('top','Fitted pos:%.3f mm'%pos
              +', sigma:%.2f mm'%sig)
            try:
                fitY = func_sum_of_peaks(x,*(pp))\
                  /xyArray.shape[self.axis]
                self.widget.plot(xPoints,fitY,
                  pen=pg.mkPen(color=(0, 0, 255), style=QtCore.Qt.DashLine))
            except: pass
          
    #def __del__(self):
    #    self.widget.win.close()
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#````````````````````````````Interactive console``````````````````````````````
Console = True
if Console: 
    #````````````````````````````Bug fix in pyqtgraph 0.10.0``````````````````
    import pickle
    class CustomConsoleWidget(pyqtgraph.console.ConsoleWidget):
        """ Fixing bugs in pyqtgraph 0.10.0:
        Need to rewrite faulty saveHistory()
        and handle exception in loadHistory() if history file is empty."""
        def loadHistory(self):
            """Return the list of previously-invoked command strings (or None)."""
            if self.historyFile is not None:
                try:
                    pickle.load(open(self.historyFile, 'rb'))
                except Exception as e:
                    printw('History file '+' not open: '+str(e))

        def saveHistory(self, history):
            """Store the list of previously-invoked command strings."""
            #TODO: no sense to provide history argument, use self.input.history instead
            #printd('>saveHistory')
            if self.historyFile is not None:
                #bug#pickle.dump(open(self.historyFile, 'wb'), history)
                pickle.dump(history,open(self.historyFile, 'wb'))
#`````````````````````````````````````````````````````````````````````````````
class CustomViewBox(pg.ViewBox):
    """ defines actions, activated on the right mouse click in the dock
    """
    def __init__(self, **kwds):
        self.dockName = kwds['name'] # cannot use name due to an issue in demo
        del kwds['name'] # the name in ViewBox.init fails in demo
        #printd('CustomViewBox: '+str(self.dockName)+', '+str(kwds))

        # call the init method of the parent class
        super(CustomViewBox, self).__init__()
        # the above is equivalent to:#pg.ViewBox.__init__(self, **kwds)

        # IMPORTANT: menu creation is deferred because it is expensive 
        # and often the user will never see the menu anyway.
        self.menu = None
           
    #v32#def mouseClickEvent(self, ev) removed, due to blank exports

    def raiseContextMenu(self, ev):
        # Let the scene add on to the end of our context menu
        menuIn = self.getContextMenus()        
        menu = self.scene().addParentContextMenus(self, menuIn, ev)
        menu.popup(ev.screenPos().toPoint())
        return True

    def getContextMenus(self, event=None):
        ''' This method will be called when this item's children want to raise
        a context menu that includes their parents' menus.
        '''
        if self.menu:
            printd('menu exist')
            return self.menu
        #print('getContextMenus for '+str(self.dockName))
        self.menu = ViewBoxMenu.ViewBoxMenu(self)
        self.menu.setTitle(str(self.dockName)+ " options..")
                   
        # showControls
        showControls = pg.QtGui.QWidgetAction(self.menu)
        showControlsGui = pg.QtGui.QCheckBox('ShowControls')
        showControlsGui.setChecked(not pargs.miniPanes)
        showControlsGui.stateChanged.connect(
          lambda x: imager.hideDocks(not x))
        showControls.setDefaultWidget(showControlsGui)
        self.menu.addAction(showControls)
        
        # AxesMM
        axesMM = pg.QtGui.QWidgetAction(self.menu)
        axesMMGui = pg.QtGui.QCheckBox('Millimeters')
        #if imager.axesInMm:
        axesMMGui.setChecked(imager.axesInMm)
        axesMMGui.stateChanged.connect(lambda x: imager.set_axes_scale(x))
        axesMM.setDefaultWidget(axesMMGui)
        self.menu.addAction(axesMM)

        # Isocurves
        isoCurves = pg.QtGui.QWidgetAction(self.menu)
        isoCurvesGui = pg.QtGui.QCheckBox('IsoCurves')
        isoCurvesGui.setChecked(1)
        isoCurvesGui.stateChanged.connect(
          lambda x: imager.show_isocurves(x))
        isoCurves.setDefaultWidget(isoCurvesGui)
        self.menu.addAction(isoCurves)
        
        return self.menu
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
def MySlot(a):
    """Global redirector of the SignalSourceDataReady"""
    #printd('>MySlot received event:'+str(a))
    if imager:
        imager.process_image()
    else:
        printe('Imager not defined yet')
#````````````````````````````Custom parameter widgets``````````````````````````
from pyqtgraph import ViewBoxMenu
class SliderParameterItem(WidgetParameterItem):
    """Slider widget"""
    def makeWidget(self):
        w = QtGui.QSlider(QtCore.Qt.Horizontal, self.parent())
        w.sigChanged = w.valueChanged
        w.sigChanged.connect(self._set_tooltip)
        self.widget = w
        return w

    def _set_tooltip(self):
        self.widget.setToolTip(str(self.widget.value()))

class SliderParameter(Parameter):
    itemClass = SliderParameterItem
    ''' no response on that
    sigActivated = QtCore.Signal(object)
    
    def activate(self):
        self.sigActivated.emit(self)
        self.emitStateChanged('activated', None)
    '''
        
registerParameterType('slider', SliderParameter, override=True)
    
class ButtonParameterItem(ParameterItem):
    """ Button with reduced spacing than 'action'"""
    def __init__(self, param, depth):
        ParameterItem.__init__(self, param, depth)
        self.layoutWidget = QtGui.QWidget()
        self.layout = QtGui.QHBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layoutWidget.setLayout(self.layout)
        self.button = QtGui.QPushButton(param.name())
        #self.layout.addSpacing(100)
        self.layout.addWidget(self.button)
        self.layout.addStretch()
        self.button.clicked.connect(self.buttonClicked)
        param.sigNameChanged.connect(self.paramRenamed)
        self.setText(0, '')
        
    def treeWidgetChanged(self):
        ParameterItem.treeWidgetChanged(self)
        tree = self.treeWidget()
        if tree is None:
            return
        
        tree.setFirstItemColumnSpanned(self, True)
        tree.setItemWidget(self, 0, self.layoutWidget)
        
    def paramRenamed(self, param, name):
        self.button.setText(name)
        
    def buttonClicked(self):
        self.param.activate()

class ButtonParameter(Parameter):
    """Used for displaying a button within the tree."""
    itemClass = ButtonParameterItem
    sigActivated = QtCore.Signal(object)
    
    def activate(self):
        self.sigActivated.emit(self)
        self.emitStateChanged('activated', None)

registerParameterType('button', ButtonParameter, override=True)

class ComplexParameter(pTypes.GroupParameter):
    def __init__(self, **opts):
        opts['type'] = 'bool'
        opts['value'] = True
        pTypes.GroupParameter.__init__(self, **opts)
        
        self.addChild({'name': 'X(pixels)', 'type': 'float', 
          'value': imager.ref_X, 'step':10.})
        self.addChild({'name': 'Y(pixels)', 'type': 'float', 
          'value': imager.ref_Y, 'step':10.})
        self.addChild({'name': 'D(pixels)', 'type': 'float', 
          'value': imager.ref_diameter, 'step':10.})
        self.addChild({'name': 'Target(mm)', 'type': 'float',
          'value': round(imager.ref_diameter/imager.pixelPmm[1],3)})
        self.addChild({'name': 'XPix/mm', 'type': 'float',
          'value': round(imager.pixelPmm[0],3),'readonly':True})
        self.addChild({'name': 'YPix/mm', 'type': 'float',
          'value': round(imager.pixelPmm[1],3),'readonly':True})
        if pargs.expert:
            self.addChild({'name':'Save', 'type': 'button'})
            self.p_save = self.param('Save')
            self.p_save.sigActivated.connect(self.p_save_changed)
        self.p_refd = self.param('D(pixels)')
        self.p_target = self.param('Target(mm)')
        self.p_XpixMm = self.param('XPix/mm')
        self.p_YpixMm = self.param('YPix/mm')
        self.p_X = self.param('X(pixels)')
        self.p_Y = self.param('Y(pixels)')
        self.p_refd.sigValueChanged.connect(self.p_refd_changed)
        self.p_target.sigValueChanged.connect(self.p_refd_changed)
        #self.p_XpixMm.sigValueChanged.connect(self.p_XpixMm_changed)
        self.p_X.sigValueChanged.connect(self.update_refRing)
        self.p_Y.sigValueChanged.connect(self.update_refRing)
        
    def p_save_changed(self):
        # update DB
        if True:#try:
            import Sybase
            mdb = Sybase.Connection('OPSYB1', 'ops', 'opsops', 'appConfig')
            mcrs = mdb.cursor()
            def set_camera_DB(db,cursor,row,col,value):
                updateCmd = 'UPDATE cameras SET '+col+' = '+str(value)\
                  +' WHERE camName = "'+row+'"'
                #print('executing DB',updateCmd)
                cursor.execute(updateCmd)
                db.commit()
            ppmY = self.p_YpixMm.value()
            cols = {'ref_X':self.p_X.value(),
                    'ref_Y':-self.p_Y.value(),# DB (0,0) is top left corner
                    'pixelPmm':ppmY,
                    'target':self.p_refd.value()/ppmY}
            if not qMessage('Are you sure to change DataBase %s?'%\
              ('appConfig/'+imager.cameraName[0])+'\nFields:\n'+str(cols)
              #'\nYou will need to restart imageMan for this camera using StartUp.'
              ):
                return
            for c,v in cols.items():
                set_camera_DB(mdb,mcrs,imager.cameraName[0],c,v)
            v1 = round(ppmY,3)
            v0 = round(ppmY*imager.pixelXScale,3)
            adoName = 'img.'+imager.cameraName[0]
            r = pvMonitor.httpDev.set(adoName,'pixelPmmM',[v0,v1])
            #print('httpDev',r)
            cprint('DataBase and Imager updated successfully')
        else:#except Exception as e:
            cprinte('in update DB :'+str(e))
        return
 
    def update_refRing(self):
        imager.ref_diameter = self.p_refd.value()
        imager.ref_X = self.p_X.value()
        imager.ref_Y = self.p_Y.value()
        imager.remove_ref_ring()
        imager.create_refRing()
    
    def p_refd_changed(self):
        ypmm = self.p_refd.value()/self.p_target.value()
        self.p_YpixMm.setValue(round(ypmm,3))#, blockSignal=self.p_XpixMm_changed)
        self.p_XpixMm.setValue(round(ypmm*imager.pixelXScale,3))
        imager.pixelPmm[0] = self.p_XpixMm.value()
        imager.pixelPmm[1] = self.p_YpixMm.value()
        self.update_refRing()
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#````````````````````````````Image decoders/coders````````````````````````````
# the PIL is best codec, as fast as QT, it takes any pixel depth.
# It can be substituted with even faster Pillow-SIMD 2.5 which is clamed to be
# 2.5 times faster than Pillow and 10 times faster than ImageMagick.
Image = None
class CodecPIL():
    def __init__(self):
        global Image, fromimage
        from PIL import Image

    def load(self,fname):
        #try:
        self.img = Image.open(fname)
        return np.asarray(self.img)
        
    def loadFromData(self,blob): #
        import io
        self.img = Image.open(io.BytesIO(blob))
        return np.asarray(self.img)

    def save(self,fileName,npArray=None):
        #if self.img:
        if npArray is None:
            self.img.save(fileName,'PNG')
        else:
            vertFlipped = npArray[::-1,...]
            Image.fromarray(vertFlipped).save(fileName,'PNG')
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       
#````````````````````````````PV Monitor Objects for different access systems``
''' The Monitor object is instantiated as:
pvm = PVMonitorXXX(sourceName, reader = readerName)
Derived class must override at two functions: getTimeStamp(), get_data_and_timestamp() and getSensImageShape()
'''
#````````````````````````````Base PVMonitor class`````````````````````````````
#class PVMonitor(object): # need to be derived from object for super to work
class PVMonitor(QtCore.QThread): # inheritance from QtCore.QThread is needed for qt signals
    SignalSourceDataReady = QtCore.pyqtSignal(object)
    def __init__(self):
        super(PVMonitor,self).__init__()
        self.eventNumber = 0
        self.data = np.array([])
        self.paused = False
        self.exit = False
        self.dbg = None # debugging flag
            
    #def getSensImageShape(self): return self.data.shape[:2]
    
    def set_refresh(self,r):
        rr = 0.99/r
        printi('Refresh rate changed to %0.4g'%rr)
        self.refreshTime = rr

    def monitor(self):
        """starts the monitoring"""
        printi('pvmonitor.monitor() is not instrumented') 
        
    def clear(self):
        """clears a monitor, stops related threads"""
        printi('pvmonitor.clear() is not instrumented') 

    def getTimeStamp(self):
        """returns timestamp, used for polling data delivery"""
        printi('pvmonitor.getTimeStamp() is not instrumented')
        
    def get_data_and_timestamp(self):
        """returns the image ndarray and timestamp used for polling data delivery
        if stream has been finished then the ndarray is empty, timestamp = 0
        """
        printi('pvmonitor.get_data_and_timestamp() is not instrumented')
        return [],0
        
    def get_threshold_and_ROI(self,cameraName):
        """returns threshold level, ROI and pedestal subtraction"""
        return 19.9,[0.5,0.5,0.9,0.9],0
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#````````````````````````````Image access objects`````````````````````````````            
#````````````````````````````PVMonitor for data from a file```````````````````
class PVMonitorFile(PVMonitor):
    def __init__(self,pvname,**kwargs):
        super(PVMonitorFile,self).__init__()
        #print('pvname:'+str(pvname))
        if len(pvname)==1:
            import glob
            if os.path.isdir(pvname[0]):
                globname = pvname[0]+'/*.png'
            else:
                globname = pvname[0]
            self.fileList = sorted(glob.glob(globname))
            #printd('glob size %d'%len(self.fileList))
        else:
            self.fileList = pvname
        nFiles = len(self.fileList)
        print('Files in the stream: %d'%nFiles)
        if nFiles == 0:
            # select files of today from the camera
            globname = '/operations/app_store/cameras/run_fy19/'\
              +pargs.controlSystem+'/'+pvname[0]+'/'+pargs.prefix+pvname[0]\
              +time.strftime('_%Y%m%d*')
            #printd('globname: '+str(globname))
            self.fileList = sorted(glob.glob(globname))
            #printd('today glob size %d'%len(self.fileList))
            if len(self.fileList) == 0:
                printe('No such files: '+str(globname))
                sys.exit()
                       
        try: self.cameraName = kwargs['camName']
        except: self.cameraName = 'No cameraName'
        try: r = kwargs['refreshRate']
        except: r = 1.
        self.set_refresh(r)
        self.pvsystem = 'File' # for Accelerator Device Objects, ADO
        self.qimg = QtGui.QImage() # important to have it persistent
        self.profTime = time.time()
        self.profN = len(self.fileList)
        self.curFileIdx = 0
        self.lastFileIdx = len(self.fileList)
                
        thread = threading.Thread(target=self.thread_process)
        thread.start()
        
        # thread safe data delivery
        self.SignalSourceDataReady.connect(MySlot)
                   
    def trash(self,fromIdx):
        trashDir = '/operations/app_store/cameras/Trash/'+self.cameraName+'/'
        if not os.path.exists(trashDir):
            os.makedirs(trashDir)
        print('created ',trashDir)
        toIdx = self.lastFileIdx
        nImages = toIdx - fromIdx
        try:
            firstFile = self.fileList[fromIdx]
            lastFile = self.fileList[toIdx-1]
            try:
                firstFile = firstFile.rsplit('/',1)[1]
                lastFile = lastFile.rsplit('/',1)[1]
            except: pass
            ok = qMessage('Moving %d files '%nImages+' out of %d, from '%len(self.fileList)\
              +firstFile+' through '+lastFile+' to '+trashDir)
            if not ok:
                return
            #print('deleting',fromIdx,toIdx)
            for fullname in self.fileList[fromIdx:toIdx]:
                #source,fn = self.fileIter.next().rsplit('/',1)
                try:
                    source,fn = fullname.rsplit('/',1)
                    source += '/' + fn
                except:
                    source, fn = fullname, fullname
                os.rename(source,trashDir+fn)
                #print('trashing '+str(source+'/'+fn))
            cprint('Moved %d files '%nImages+' out of %d, from '%len(self.fileList)\
              +firstFile+' through '+lastFile+' to '+trashDir)
            self.fileList[fromIdx:toIdx] = []
            self.curFileIdx = fromIdx
            #print('after:',len(self.fileList),self.curFileIdx)
        except Exception as e:
            cprinte('in trash(): '+str(e))
    
    def thread_process(self):
        #TODO: using cprint is not thread-safe as it calls GUI.
        while True:
            EventPocessingFinished.wait(2) #
            if not EventPocessingFinished.is_set():
                #print('Processing timeout')
                if self.exit:
                    cprint('FIXME.File streaming finished. You have to restart the program')
                    return
            else:
                EventPocessingFinished.clear()
                #print('EventPocessingFinished',self.eventNumber)
                time.sleep(self.refreshTime)
                #self.getTimeStamp()
                self.SignalSourceDataReady.emit('sourceDataReady')
            
    def clear(self):
        self.exit = True

    def get_data_and_timestamp(self):
        # get the next file
        try:
            fname = self.fileList[self.curFileIdx]
            if self.curFileIdx > self.lastFileIdx:
                raise IndexError
        except IndexError:
            #printd('No more files')
            # no more files, do not update timestamp
            #if not self.noMore:
            cprint('Processed %d images'%self.profN+' in %0.4g s'%(time.time()-self.profTime))
            #self.noMore = True
            if imager.addon: imager.addon.stop()
            try:
                imager.gl
                imager.gl.setBackground('y')
            except:
                printe('Processing first event')
            time.sleep(1)
            return [],0

        self.curFileIdx +=1
        self.pvname = fname
        #printd('image:'+self.pvname)
        # extract timestamp from the filename
        try:
            ds,dt = fname.split('_')[-2:]
            dsdt = (ds+dt)
            dsdt.strip('.png')
            seconds = dsdt[:14]
            #print(seconds)
            timestamp = time.mktime(datetime.datetime.strptime(seconds,"%Y%m%d%H%M%S").timetuple())
            timestamp += float(dsdt[14:20])*1.e-6 # add microsecond part
        except:
            timestamp = time.time()
        #print('timestamp:',timestamp)
        
        # load image
        self.data = codec.load(fname)
        self.eventNumber +=1
        return self.data,timestamp

#````````````````````````````PVMonitor for a HTTP image```````````````````````
class PVMonitorHTTP(PVMonitor):
    def __init__(self,pvname,**kwargs):
        super(PVMonitorHTTP,self).__init__()
        self.pvsystem = 'HTTP'
        self.qimg = QtGui.QImage() # important to have it persistent
        self.req = Backend.get(pvname)
        
    def get_data_and_timestamp(self):
        self.eventNumber +=1
        if pargs.dbg:
            print('>http.get_data')
            print('encoding:',self.req.encoding)
            print('status_code:',self.req.status_code)
            print('elapsed:',self.req.elapsed)
            #prin(':',self.req.url)
            print('history:',self.req.history)
            print('Content-Type:',self.req.headers['Content-Type'])
            print('cont:',type(self.req.content),len(self.req.content),self.req.content[:20])
        printw('FixIt: HTTP:get_data_and_timestamp')
        self.data = None
        #self.qimg.loadFromData(self.req.content)
        #self.data = convert_qimg_to_ndArray(self.qimg)
        return self.data,time.time()

    def clear(self):
        return

#````````````````````````````PVMonitor of a Process Variable from ADO system``
class PVMonitorAdo(PVMonitor):    
    def __init__(self,adoName,**kwargs):
        global pargs
        super(PVMonitorAdo, self).__init__() # for signal/slot paradigm we need to call the parent init
            
        from cad import pyado
        self.cadDev = pyado.useDirect()
        self.httpDev = pyado.useHttp()# http supports set history
        self.adoName = adoName
        #print('kwargs:',kwargs)
        try: r = kwargs['refreshRate']
        except: r = 1.
        self.set_refresh(r)
        self.pvsystem = 'ADO' # for Accelerator Device Objects, ADO
        self.qimg = QtGui.QImage() # important to have it persistent
        self.par = 'imageM' if pargs.fullsize else 'gmImageM'
        self.handle = Backend.CreateAdo(adoName)
        if self.handle is None:
            printe('cannot create '+adoName)
            sys.exit(1)
        #print('ado:',self.handle.systemName,self.handle.genericName)
        self.adoName = self.handle.genericName

        # check if parameter exists
        metaData = Backend.adoMetaData(self.handle)
        if not isinstance(metaData, dict):
            printe("while trying to get metadata"+str(metaData))
            sys.exit(2)
        if (self.par,'value') not in metaData:
            printw('ado '+adoName+' does not have '+self.par)
            self.par = 'minimageM'
            if (self.par,'value') not in metaData:
                printe('ado '+adoName+' does not have '+self.par)
                sys.exit(3)
        self.pvname = adoName+':'+self.par
        
        # get sensor size for calibration
        w = Backend.adoGet(self.handle,'imageWidthM','value')[0][0][0]
        h = Backend.adoGet(self.handle,'imageHeightM','value')[0][0][0]
        #self.sensImageShape = (h,w)
        
        if self.par is not 'gmImageM':
        # find width/height from other parameters
            pargs.width = self.get_image_shape()
            pargs.extract = 'raw'
            print('--width modified: '+pargs.width)
        
        # store the requests for future use
        self.dataRequest = [(self.handle, self.par, 'value'),
          (self.handle,self.par,'timestampSeconds'),
          (self.handle,self.par,'timestampNanoSeconds')]

        # check if parameter has timestamp property
        ts = int(str(Backend.adoGet(self.handle,self.par,'timestampSeconds')[0][0][0]))
        p = 'nImagesM' if ts == 0 else self.par 
        self.tsRequest = [(self.handle,p,'timestampSeconds'),          
          (self.handle,p,'timestampNanoSeconds'),
          (self.handle,p,'value')]
          
        self.asyncPar = 'nImagesM'
        self.tsKey = self.adoName+':'+self.asyncPar
        self.dataKey = self.adoName+':'+self.par
        self.ts = self.cadDev.get(self.adoName,self.asyncPar)[self.tsKey]['timestamp']
        r = self.cadDev.getAsync(self.callback,self.adoName,self.asyncPar)
        self.SignalSourceDataReady.connect(MySlot)
        
    def clear(self):
        self.cadDev.cancelAsync()
    
    def callback(self,*args):
        profile('start')
        #print('cb:',args)
        props = args[0][self.tsKey]
        #print(props)
        ts = props['timestampSeconds']+props['timestampNanoSeconds']*1.e-9
        self.sourceImageN = props['value']
        #print('source image %d, at %.4f, last:%.4f'%(self.sourceImageN,ts,self.ts))
        if ts >= self.ts + self.refreshTime:
            self.ts = ts
            if EventPocessingFinished.is_set() or self.eventNumber == 0:
                EventPocessingFinished.clear()
                #print('sourceDataReady')
                self.SignalSourceDataReady.emit('sourceDataReady')            
        
    def blobToNdArray(self,blob):
        if pargs.fullsize:
            data = np.array(blob,'u1')
        else:
            data = codec.loadFromData(blob)
        return data
        
    def get_data_and_timestamp(self):
        printd('>get_data')
        self.eventNumber +=1
        data = None
        #r = self.cadDev.get(self.adoName,self.par)[self.dataKey]
        
        # get data with minimal overhead
        items, status = Backend.adoGet(list = self.dataRequest)
        profile('adoGet')
        #print('dr:',self.dataRequest)
        if items is None or len(items) == 0:
            cprintw('no data from '+self.pvname+' is manager all right? '+str(status))
            return [],0
        else:
            blob = items[0]
            if len(blob) <= 0:
                return [],0
            # print('got from '+self.pvname+': blob['+str(len(blob))+']: '+repr(blob[:100]))
            data = self.blobToNdArray(blob)
            timestamp = items[1][0] + items[2][0]*1.e-09
        #print('timestamp:%.4f'%timestamp)
        #print('gdts:',data.shape)
        self.blob = blob
        return data,timestamp
        
    def get_threshold_and_ROI(self,cameraName):
        adoName = 'img.'+cameraName
        d = {'thresholdS':19.8, 'roiS':[0.06,0.06,0.9,0.9], 'subtractPedS':0}
        for par in d:
            adoKey = adoName+':'+par
            try:
                v = self.cadDev.get(adoName,par)[adoKey]['value']
                d[par]=v
            except Exception as e:
                cprintw('Cannot get '+adoKey+': '+str(e))
        return d['thresholdS'],d['roiS'],d['subtractPedS']

    def get_image_shape(self):
        # find width/height from other parameters
        bitsPerPixel = ''
        if self.par is 'imageM':
            #bpp = Backend.adoGet(self.handle,'imageBppM','value')[0][0][0]
            bpp = Backend.adoGet(self.handle,'bppM','value')[0][0][0]
            bitsPerPixel = ','+str(bpp)
        w = str(Backend.adoGet(self.handle,self.par[:-1]+'WidthM','value')[0][0][0])
        h = str(Backend.adoGet(self.handle,self.par[:-1]+'HeightM','value')[0][0][0])
        # substitute -w and -e options
        return w+','+h+bitsPerPixel

#````````````````````````````PVMonitor of a an EPICS Process Variable`````````
class PVMonitorEpics(PVMonitor):
    def __init__(self,pvname,**kwargs):
        super(PVMonitorEpics,self).__init__()
        self.pvsystem = 'Epics' # for Accelerator Device Objects, ADO
        #epics.ca.replace_printf_handler(self.handle_messages)
        self.pv  = Backend.PV(pvname)
        self.ts = 0
        print('Epics:',self.pv)
        
    def getTimeStamp(self):
        r = self.pv.get_timevars()
        if r is None:
            printe('No PV:'+str(self.pv))
            sys.exit(4)
        return r['timestamp']
        
    def get_data_and_timestamp(self):
        if pargs.dbg:
            print('>data, dt:%0.4g'%(self.pv.timestamp - self.ts))
            #print(pargs)
        self.ts = self.pv.timestamp
        data = self.pv.get()
        self.blob = blob
        return data,self.ts
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#```````````````````````````Image viewing/processing object```````````````````
class Imager(QtCore.QThread): # for signal/slot paradigm the inheritance from QtCore.QThread is necessary
    SignalCPrint = QtCore.pyqtSignal(object) # for thread-safe cprint
    def __init__(self,pvname,camName):
        super(Imager, self).__init__() # for signal/slot paradigm we need to call the parent init
        self.viewBoxBackgroundColor = 0.9
        self.pvname = pvname
        self.cameraName = camName
        self.qimg = QtGui.QImage()
        self.hwpb = [0]*4 # height width, number of planes, bits/channel
        self.dataLen = 0
        self.sizeFactor = 0.
        self.plot = None # ROI projection plot
        self.roi = None # Region Of Interest
        self.roiPlotEnabled = True
        self.contrast = None # Contrast histogram
        self.roiRect = []
        try: 
            pargs.roiRect = [float(i) for i in pargs.roiRect.split(',')]
            self.roiRect = pargs.roiRect
        except: pass
        if len(self.roiRect) != 4:
            cprinte('setting ROI '+pargs.roiRect)
            self.roiRect = [0.05,0.05,0.9,0.9] # set default
        self.threshold = pargs.threshold
        
        self.iso = None # Isocurve object
        self.isoInRoi = False
        self.data = np.array([])
        self.grayData = None
        self.mainWidget = None
        self.imageItem = None
        self.docks = {}
        self.needToRedrawAxes = False
        self.dockParRotate = 0 # rotation angle
        self.spotLog = None # enable logging of calculated spot parameters
        self.events = 0 # counter of processed events
        #see below in start()#self.threshold = pargs.threshold # threshold for image thresholding
        self.cleanImage = False # show image only
        self.zoomROI = False # zoom image to ROI
        self.maxSpots = pargs.maxSpots # number of spots to find
        self.spots = [] # calculated spot parameters
        self.spotShapes = [] # spot finder graphics objects
        self.marksColor = (0,170,0) # color of marks and spot contours
        self.mainSpotTxt = ''
        self.roiArray = [] # array of ROI-selected data
        self.isocurve_enabled = False
        self.saving = pargs.saving # enable the continuous saving of images
        self.nImagesToTrash = 0 
        self.refresh = pargs.refreshRate # refresh rate [Hz]
        self.sleep = 0 # debugging sleep
        self.timestamp = -1 # timestamp from the source
        self.paused = pargs.pause # pause processing    
        if not self.paused:
            self.timestamp = 0 # invalidate timestamp to get one event
        self.rawData = None # rawData from reader
        self.stopProcThread = False # to stop procThread
        self.refRing = None
        self.ref_diameter = None
        self.pixelXScale = 1. # used in calibrations
        self.axesInMm = True # convert spot parameters to milimeters
        self.pedestals = None
        self.averageWindow = pargs.average
        self.imageFile = 'noPath/noName'
        self.backend = pargs.backend
        self.borderPen = pg.mkPen('g')
        self.outsideGraphs = {}
        self.roiIntensityPlot = None
        #self.roiVerticalPlot = None
        self.imageMan = None # handle to imageMan ADO if 'UpdateMan' is on
        #self.update_imageMan_locked = False
        #self.is_host_priviledged = os.uname()[1][:6]=='acnmcr'
        self.is_host_priviledged = True
        self.blurWidth = 0.
        self.fitBaseBounds = [0., +np.inf]

        self.emit_slit_width = 0.
        self.emit_slit_space = 0.
        self.emit_drift_length = 0.
        
        # connect signal to slot
        #v76#self.signalDataArrived.connect(self.process_image) # use mySlot because cannot connect external slot
        #for thread safe use of cprint in threads
        self.SignalCPrint.connect(cprint)

        #destroyIssue#self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)
        #destroyIssue#self.destroyed.connect(self.closeEvent)

        profile('start')

    #destroyIssue#@staticmethod
    #destroyIssue#def closeEvent(self):
    #destroyIssue#    print("Closing")

    def start_imager(self):
        if self.cameraName is None:
            if self.backend in ['ado','epics']:
                self.cameraName = [self.pvname,'?']
        if self.cameraName is None:
            self.savePath = '/tmp/imageViewer/'
            self.cameraName = ['noName','?']
        else:
            self.savePath = pargs.logdir+pargs.controlSystem+'/'+self.cameraName[0]+'/'
        if self.backend == 'file':
            self.imageFile = self.pvname
        checkPath(self.savePath)
        self.refPath = self.savePath+'refs/'
        checkPath(self.refPath)
        self.logPath = self.savePath+'logs/'
        checkPath(self.logPath)

        if pargs.userAddon:
            addon = 'ivAddon_'+pargs.userAddon
            Addon = importlib.import_module(addon)
            self.addon = Addon.Addon(self)
        else:
            self.addon = None        
        
        # provide status of 'Subtract' checkbox
        if pargs.backend == 'ado':
            self.threshold,self.roiRect,subtractPeds =\
              pvMonitor.get_threshold_and_ROI(self.cameraName[0])
        else: subtractPeds = 0
        if pargs.subtract: subtractPeds = 1
        #print('subtractPeds',subtractPeds)
        self.subtractPeds = ('None','ref0')[subtractPeds]

        # get first image
        ts = timer()
        while len(self.data) == 0 and timer() - ts < 10:
            self.process_image()
        if len(self.data) == 0:
            printe('Could not get first event in 10 s')
            sys.exit(1)

        # first event available, subtraction can be activated.
        self.pedestals = self.enable_subtraction(self.subtractPeds)

        #print('threshold,roi,pix/mm:',(self.threshold,self.roiRect,self.pixelPmm))
        #if self.axesInMm:
        # Fixing the loss of the axis scalings when mainWidget was changed.
        # There should be a better way to keep axis scalings persistent
        self.mainWidget.sigRangeChanged.connect(self.mainWidget_changed)
        #print('<start')
        
        if pargs.refRing:
            self.create_refRing()

    def update_imageMan(self,newState=2):
        if self.backend != 'ado':
            return
        from cad import cns
        #if self.update_imageMan_locked: # should not happen
        #    return 0
        if newState == False:
            if self.imageMan is not None:
                print('Updating of the imageMan is disabled ')#+self.imageMan.genericName)
            self.imageMan = None
            return 0 # normal return
            
        elif newState == True:
            if self.imageMan is not None:
                print('Strange, '+self.imageMan.genericName+' is already open')
            adoName = 'img.'+self.cameraName[0]
            #self.imageMan = cns.CreateAdo(adoName)
            self.imageMan = adoName
            if self.imageMan is None:
                qMessage('ERROR: imageMan for '+adoName+' not running')
                return 2
         
        if self.imageMan is None: # no updating required
            return 0
            
        # update imageMan
        if self.addon:
            try:
                self.addon.update_imageMan(adoName)
            except Exception as e:
                cprinte(' in addon.update_imageMan: '+str(e))           
        d = {}
        for parval in (('thresholdS',self.threshold),
                       ('fitS',pargs.finalFit),
                       ('subtractPedS',self.pedestals is not None),
                       # roiS can be updated only in 'expert' mode
                       ('roiS',self.roiRect) if pargs.expert else ('',None),
                       ):
            par,val = parval
            if len(par) == 0:
                continue 
            try:
                if False:
                    r = cns.adoGet(self.imageMan, par, 'value')
                    old = r[0][0]
                    if len(old) == 1: 
                        old = old[0]
                    else:
                        old = list(old)
                    imageManName = self.imageMan.genericName
                imageManName = self.imageMan
                key = self.imageMan+':'+par
                r = pvMonitor.cadDev.get(self.imageMan,par)[key]['value']
                old = r
            except Exception as e:
                cprinte('getting '+str(parval)+',exception: '+str(e))
                continue
            d[par] = {'old':old}
            if parval[0] == 'fitS':
                val = 'FinalFit' if val else 'FitLess'
            if old == val:
                #printi('old value of '+par+' is the same, no need to change')
                continue
            txt = 'Are you going to modify '+par+' on imageMan for '+\
              imageManName+' from '+str(old)+' to '+str(val)
            ok = qMessage(txt)
            if not ok:
                continue
            status = pvMonitor.httpDev.set(self.imageMan,par,val) # http supports set history
            #status = pvMonitor.cadDev.set(self.imageMan,par,val)
            #print('status',status)
            d[par]['new'] = val
            #print('set succeeded')
        #printi('update_imageMan:'+str(d))

    def open_spotLog(self):
        pname = self.pvname
        if self.backend == 'file':
            # strip the path, leave only the filename
            pname = pname.split('/')[-1:][0].rsplit('.',1)[0]
        checkPath(self.logPath)
        fn = self.logPath+'sl_'+pname\
          +time.strftime('_%y%m%d%H%M%S.json')
        self.spotLog = open(fn,'w',1)
        cprint('file: '+self.spotLog.name+' opened')
        
        # dump static parameters
        logdata = OrderedDict([\
          ('version', __version__),
          ('refreshRate', self.refresh),
          ('roi', self.roiRect),
          ('threshold', self.threshold),
          ('imageHWPB', self.hwpb),
          ('sizeFactor', self.sizeFactor),
          ('pixelPmm', [round(r,4) for r in self.pixelPmm]),
          ('imageCount', self.events),         
          ('spotInfo', 'mean(x,y),sigma(x,y),pcc,sum'),
        ])
        self.spotLog.write(dumps(logdata))#,indent=2))
        self.spotLog.write('\n')
        #printd('open,written:'+str(logdata))      

    def iv_show(self):
        """ Display the widget, called once, when the first image is available"""
        #printd('>iv_show')
        # Update sizeFactorM only once at first image
        w = float(max(self.data.shape))# - 1 #-1 is a kludge
        #sensorShape = pvMonitor.getSensImageShape()
        #self.sizeFactor = round(max(sensorShape)/w)
        sensorShape = self.hwpb[:2]
        try:
            ss = pargs.sensorShape.split(',')
            sensorShape[1],sensorShape[0] = [int(z) for z in ss[:2]]
        except: printe('sensorShape wrong: '+str(pargs.sensorShape))
        self.sizeFactor = float(min(sensorShape))/min(self.hwpb[:2])
        # ratio of heights is correct even when the rows were padded (by a PNG formatter)
        print('sizeFactor:%.2f'%self.sizeFactor+', sensorShape,imageShape: '\
          +str((sensorShape,self.hwpb)))

        self.win = QtGui.QMainWindow()
        #printd('mainWindow created')
        #does not work#self.win.setAttribute(QtCore.Qt.WA_DeleteOnClose, True)
        #destroyIssue#self.win.destroyed.connect(self.closeEvent)

        self.area = pg.dockarea.DockArea()
        self.win.setCentralWidget(self.area)
        #destroyIssue#self.connect(self,QtCore.SIGNAL('triggered()'), self.closeEvent)
        
        # image info
        self.winTitle = 'image:'+self.pvname+' hwpb:'+str(self.hwpb)
        self.win.setWindowTitle(self.winTitle)

        #````````````````````````````parameters data``````````````````````````
        self.numRefs = 5
        refs = ['None']+['ref'+str(i) for i in range(self.numRefs)]
        refreshList = ['1Hz','0.1Hz','10Hz','Instant']
        #if self.backend == 'file': 
        #    refreshList.append('Instant')
        if self.addon:
            self.addonEntry = self.addon.entryName()
            paramAddon = {'name': self.addonEntry, 'type': 'group', 'children': self.addon.controlPane()}
            #paramAddon = {'name': 'Addon', 'type': 'group', 'children': self.addon.controlPane()}
        else:
            self.addonEntry = '??????'
            paramAddon = {'name': 'Addon', 'type': 'group','visible':False}
        params = [
            {'name': 'Control', 'type': 'group', 'children': [
                {'name':'Image#','type':'str','value':'0',
                  'tip':'Accepted events','readonly':True},
                {'name':'Pause', 'type': 'bool', 'value': self.paused,
                  'tip': 'Enable/disable receiving of images, '},
                {'name':'Next', 'type': 'button',
                  'tip':'Process one image from the stream'},
                {'name':'Fast', 'type': 'bool', 'value': 0,
                  'visible':True if self.backend == 'file' else False,
                  'tip': 'Fast browsing, many features disabled'},
                {'name':'Images/s','type':'float','value':0.,
                  'tip':'Performance','readonly':True},
                {'name':'FirstFile%', 'type':'slider', 'value':0,
                  'visible':True if self.backend == 'file' else False,
                  'tip':'Starting file relative position in directory.'},
                {'name':'LastFile%', 'type':'slider', 'value':100,
                  'visible':True if self.backend == 'file' else False,
                  'tip':'Ending file relative position in directory.'},
                {'name':'Saving', 'type': 'bool', 'value': self.saving,
                  'visible':False if self.backend == 'file' else True,
                  'tip': 'Enable/disable saving of images'},
                {'name':'Trash', 'type': 'button',
                  'visible':True if self.backend == 'file' else False,
                  'tip':'Move selected images to trash'},                
                {'name':'View saved', 'type': 'button',
                  'visible':False if self.backend == 'file' else True,
                  'tip':'View saved images from this camera using separate application'},
                {'name':'Threshold', 'type': 'float', 'value':self.threshold,
                  'tip': 'Threshold level for spot finding, changed with isoCurve level'},
                {'name':'Update Mgr', 'type': 'bool','value':False,
                  'visible':False if self.backend == 'file' else True,
                  #'readonly':True if self.is_host_priviledged else False,
                  'tip':'Update control parameters of the imageMan with local setting for continuous logging of the analysis'},                
                {'name':'LogView', 'type': 'button',
                  'visible':True if self.backend in ('ado','file') else True,
                  'tip':'Open LogView for related camera'},                
                {'name':'Gpm', 'type': 'button',
                  'visible':True if self.backend in ('ado','file') else True,
                  'tip':'Open Gpm for related camera'},
                #{'name':'Camera', 'type': 'button',
                #  'visible':True if self.backend in ('ado','file') else True,
                #  'tip':'Open Camera Management Pet page'},
                {'name':'Zoom to ROI', 'type': 'bool', 'value':self.zoomROI,
                  'tip': 'Zoom to ROI'},
                {'name':'Normalize', 'type': 'bool', 'value':pargs.normalize,
                  'tip': 'Normalize intensity'},
            ]},
            {'name': 'Configuration', 'type': 'group','expanded':False, 'children': [
                {'name':'RefRing', 'type': 'bool', 'value': False,
                  'tip':'Draw reference ring'},
                ComplexParameter(name='RefRing Adj',expanded=False),
                {'name':'Refresh Rate', 'type':'list',
                  'values':refreshList,
                  'value':str(self.refresh)+'Hz','tip':'Refresh rate'},
                {'name':'Reset ROI', 'type': 'button',
                  'tip': 'Reset ROI to original'},
                {'name':'Clean Image', 'type': 'bool', 'value':self.cleanImage,
                  'tip': 'Show image only'},
                {'name':'Rotate', 'type': 'float', 'value': 0,
                  'tip':'Rotate image view by degree clockwise'},
                #{'name':'Axes in mm', 'type': 'bool', 'value':self.axesInMm,
                #  'tip':'Convert coordinates to milimeters '},
            ]},
            {'name': 'Analysis', 'type': 'group','expanded':False, 'children': [
                {'name':'FinalFit', 'type': 'bool', 'value': pargs.finalFit,
                  'tip':'Do not fit the brightest spot for faster processing'},
                {'name':'View results', 'type': 'bool', 'value': False,
                  'tip':'Log the spots parameters to a file'},
                {'name':'Subtract', 'type':'list','values':refs,
                  'value':self.subtractPeds,
                  'tip':'Streaming subtraction of a reference image inside the ROI'},
                {'name':'Average', 'type': 'int', 'value': self.averageWindow,
                  'tip':'Moving average of several images'},
                {'name':'ROI Intensity', 'type': 'bool','value':False,
                 'tip':'Intesity histogram of all pixels in the ROI'},
                {'name':'ROI Vertical', 'type': 'bool','value':False,
                 'tip':'Vertical projection of the ROI'},
                {'name':'ROI Horizontal', 'type': 'bool','value':False,
                 'tip':'Horizontal projection of the ROI'},
            ]},
            {'name':'SpotFinder', 'type':'group','expanded':False, 'children': [
                {'name':'MaxSpots', 'type': 'int', 'value':self.maxSpots,
                  'limits':(0,pargs.maxSpots),
                  'tip': 'Max number of spots to find in the ROI'},
                {'name':'Found:', 'type': 'int', 'value':0,'readonroiArrayly':True,
                  'tip': 'Number of spots found in the ROI'},
                #{'name':'Spots', 'type':'str','value':'(0,0)',
                #  'readonly': True,'tip':'X,Y and integral of found spots'},
            ]},
            paramAddon,
            {'name':'Ref Images', 'type': 'group','expanded':False,'children': [
                {'name':'View', 'type':'list','values': refs,
                  'tip':'View reference image, use space/backspace for next/previous image'},
                {'name':'Store', 'type':'list','values': refs},
                {'name':'Retrieve', 'type':'list','values': refs},
                {'name':'Add', 'type':'list','values': refs},
                {'name':'Subtract', 'type':'list','values': refs},
            ]},
            {'name':'For Experts', 'type':'group','expanded':False,'children': [
                {'name':'Color', 'type':'list','values':['Native','Gray','Red','Green','Blue'],
                  'tip':'Convert image to grayscale or use only one color channel'},
                {'name':'Fast', 'type': 'bool', 'value': 0,
                   'tip': 'Fast and limited image processing'},
                {'name':'ContrastCtrl', 'type': 'bool', 'value': True,
                  'tip':'Show contrast control'},
                # ROI Plot is better to control using dock height
                {'name':'ROI Plot', 'type': 'bool', 'value': self.roiPlotEnabled,
                  'tip':'Update ROI plot, disable it for faster processing'},
                {'name':'SpotText', 'type': 'bool', 'value': True,
                  'tip':'Show coordinates and sigmas of the main spot'},
                {'name':'PointerText', 'type': 'bool', 'value': True,
                  'tip':'Show pointer coordinates and intensity '},
                #{'name':'Blur', 'type': 'button',
                # 'tip':'Convert the current image to gray and blur it using gaussian filter with sigma 2'},
                {'name':'Blur', 'type':'float','value':self.blurWidth,
                  'tip':'Streaming blurring of the ROI'},
                {'name':'FitBaseLo', 'type':'float','value':0.,
                  'tip':'Lover bound for fitted baseline'},
                {'name':'FitBaseHi', 'type':'float','value':1.e9,
                  'tip':'Upper bound for fitted baseline'},
                {'name':'Debug', 'type': 'bool', 'value': False},
                #{'name':'Sleep', 'type': 'float', 'value': 0},
                #{'name':'Test', 'type': 'str', 'value': 'abcd'},
                {'name':'Debug Action', 'type': 'button'},
            ]},
            {'name':'Help', 'type': 'button','tip':'User Instructiona'},
            {'name':'Exit', 'type': 'button','tip':'Exit imageViewer'},
        ]
        #```````````````````````````Create parameter tree`````````````````````````````
        ## Create tree of Parameter objects
        self.pgPar = Parameter.create(name='params', type='group', children=params)
        #printd('par tree created')
        
        # Handle any changes in the parameter tree
        def handle_change(param, changes):
            #printd('tree changes:')
            for param, change, itemData in changes:
                path = self.pgPar.childPath(param)
                if path is not None:
                    childName = '.'.join(path)
                else:
                    childName = param.name()
                #print('  parameter: %s'% childName)
                #print('  change:    %s'% change)
                if change == 'options': continue # do not print lengthy text
                #printd('  itemData:      %s'% str(itemData))
                #printd('  ----------')
            
                parGroupName,parItem = childName,''
                try: 
                    parGroupName,parItem = childName.split('.',1)
                except: None
                if parGroupName == 'Control':
                    if parItem == 'Pause':
                        #printd('Pause')
                        self.paused = itemData
                        if not self.paused:
                            self.show_isocurves(False)
                            if self.backend == 'file':
                                pvMonitor.curFileIdx = self.events
                            self.profTime = 0                            
                            EventPocessingFinished.set()
                    elif parItem == 'Next':
                        self.show_isocurves(False)
                        EventPocessingFinished.set()
                                                
                    # items for backed != 'file'
                    elif parItem == 'Saving':
                        self.saving = itemData
                        if self.saving:
                            self.save_image()
                            if self.addon: self.addon.start_saving()
                            cprint('Saving images to '+self.savePath)
                        else:
                            if self.addon: self.addon.stop_saving()
                            cprint('Stopped saving to '+self.savePath)
                    elif parItem == 'View saved':
                        if not os.path.exists(self.savePath):
                            cprinte('opening path: '+self.savePath)
                            return
                        try:
                            #cmd = ["gm",'display',self.savePath+'IV_'+self.pvname+'_*']
                            cmd = 'xterm -e imageViewer.py -p -T -m1 -bfile '+self.savePath
                            cprint('executing:'+str(cmd))
                            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True) #stderr=subprocess.PIPE)
                        except Exception as e: cprinte('in View saved'+str(e))
                    elif parItem == 'Update Mgr':
                        if self.is_host_priviledged:
                            self.update_imageMan(newState=itemData)
                        elif itemData:
                            qMessage('ERROR\nimageMan can be updated only from a MCR host')
                            # self.set_dockPar('Control','Update Mgr',False)#DNW
                    elif parItem == 'LogView':
                        cmd = None
                        try:
                            if pargs.controlSystem == 'LEREC':
                                cmd = ['LogView','-file','/operations/app_store/Gpm/LEReC/Cameras/img.'+self.cameraName[0]+'.logreq']
                            elif pargs.controlSystem == 'CEC':
                                cmd = ['LogView','-file','/operations/app_store/Gpm/RHIC/Systems/CeC/Cameras/img.'+self.cameraName[0]+'.logreq']
                            else:
                                cprinte('Control System '+str(pargs.controlSystem)+' not supported')
                            if cmd is not None:
                                p = subprocess.Popen(cmd)
                        except Exception as e: cprinte('in LogView: '+str(e))
                    elif parItem == 'Gpm':
                        cmd = None
                        try:
                            if pargs.controlSystem == 'LEREC':
                                cmd = ['Gpm','-file','/operations/app_store/Gpm/LEReC/Cameras/img.'+self.cameraName[0]+'.logreq']
                            elif pargs.controlSystem == 'CEC':
                                cmd = ['Gpm','-file','/operations/app_store/Gpm/RHIC/Systems/CeC/Cameras/img.'+self.cameraName[0]+'.logreq']
                            else:
                                cprinte('Control System '+str(pargs.controlSystem)+' not supported')
                            if cmd is not None:
                                p = subprocess.Popen(cmd)
                        except Exception as e: cprinte('in Gpm: '+str(e))
                        '''
                    elif parItem == 'Camera':
                        #cmd = ['pet','-f',
                        #  '/operations/acop/FECs/Instrumentation/EBICs/'\
                        #  +self.cameraName[1]+'/device_list.ado']
                        cmd = ['find',
                          '/operations/acop/FECs/Instrumentation/EBICs/',
                          '-name','"'+self.cameraName[1]+'"',
                          '-exec','pet -f','"{}"/device_list.ado','\;']
                        print('cmd:',cmd)
                        try:
                            p = subprocess.Popen(cmd)
                        except Exception as e: cprinte('in Gpm: '+str(e))
                        '''
                    elif parItem == 'Threshold':
                        self.threshold = itemData
                        if self.imageMan is not None: 
                            self.update_imageMan()
                        if self.isoInRoi:
                            #TODO: need to relocate the isocurves to ROI origin
                            self.iso.setData(blur(self.roiArray))
                        else:
                            self.iso.setData(blur(self.grayData))
                        self.show_isocurves(True)
                        self.isoLine.setValue(self.threshold)
                        #self.isoLine.setZValue(1000) # bring iso line above contrast controls
                        self.iso.setLevel(self.threshold)
                        if self.roi:
                            self.update_ROI()
                    elif parItem == 'Zoom to ROI':
                        self.zoomROI = itemData
                        if self.zoomROI:
                            self.setup_zoomROI()
                        else:
                            self.setup_main_widget()
                    elif parItem == 'Normalize':
                        pargs.normalize = itemData
                        self.update_imageItem_and_ROI()

                    # items for backend == 'file'
                    elif parItem == 'TrashImages':
                        self.nImagesToTrash = itemData
                    elif  parItem == 'Trash':
                        pvMonitor.trash(self.events)
                    elif parItem == 'FirstFile%':
                        #self.events = pvMonitor.jump_to(itemData)
                        self.events = file_idx(pvMonitor.fileList,itemData)
                        pvMonitor.curFileIdx = self.events
                        EventPocessingFinished.set()
                    elif parItem == 'LastFile%':
                        #self.events = pvMonitor.jump_to(itemData)
                        #EventPocessingFinished.set()
                        pvMonitor.lastFileIdx = file_idx(pvMonitor.fileList,
                          itemData)
                        pvMonitor.curFileIdx = pvMonitor.lastFileIdx
                        EventPocessingFinished.set()
                        
                if parGroupName == 'Configuration':               
                    if parItem == 'Reset ROI':
                        h,w = self.hwpb[:2]
                        self.roi.setPos(pargs.roiRect[0]*w,pargs.roiRect[1]*h)
                        self.roi.setSize(pargs.roiRect[2]*w,pargs.roiRect[3]*h)
                    elif parItem == 'Clean Image':
                        self.cleanImage = itemData
                        if self.cleanImage:
                            self.remove_spotShapes()
                            self.show_isocurves(False)
                            self.mainWidget.removeItem(self.roi)
                            self.mainWidget.removeItem(self.pointerText)
                            self.pointerText = None
                            self.mainWidget.removeItem(self.mainSpotText)
                            self.mainSpotText = None
                        else:
                            self.mainWidget.addItem(self.roi)
                            self.mainSpotText = mainSpotText()
                            self.mainSpotText.setPos(*self.top_right_corner())
                            #self.mainSpotText.setText(self.mainSpotTxt)
                            self.mainWidget.addItem(self.mainSpotText)
                            self.pointerText = pointerText()
                            self.mainWidget.addItem(self.pointerText)
                            
                    elif parItem == 'RefRing':
                        if itemData:
                            self.create_refRing()
                        else:
                            self.remove_ref_ring()
                    elif parItem == 'Rotate':
                        self.dockParRotate = float(itemData)
                        self.data = rotate(self.receivedData,self.dockParRotate)
                        self.grayData = rgb2gray(self.data)
                        self.update_imageItem_and_ROI()
                        self.update_isocurve()
                    elif parItem == 'Refresh Rate':
                        self.refresh = {'1Hz':1,'0.1Hz':0.1,'10Hz':10,'Instant':1000}[itemData]
                        pvMonitor.set_refresh(self.refresh)
                        cprint('Refresh period changed to '+str(self.refresh))
                    #elif  parItem == 'Axes in mm':
                    #    self.set_axes_scale(itemData)

                if parGroupName == 'Analysis':               
                    if  parItem == 'FinalFit':
                        pargs.finalFit = itemData
                        #print('finalFit',pargs.finalFit)
                        if self.imageMan is not None: 
                            self.update_imageMan()
                    elif parItem == 'View results':
                        if itemData:
                            self.open_spotLog()
                            if self.spotLog:
                                cmd = ['xterm','-e','tail -f '+str(self.spotLog.name)]
                                #print 'tail:',cmd
                                time.sleep(.5)
                                p = subprocess.Popen(cmd)
                        else:
                            try: 
                                self.spotLog.close()
                                cprint('file:'+str(self.spotLog.name)+' closed')
                                p = subprocess.Popen(['pkill','-9','-f',self.spotLog.name])
                            except Exception as e: cprinte('in spotLog '+str(e))
                            self.spotLog = None

                    elif parItem == 'Average':
                        self.averageWindow = itemData
                        if self.averageWindow:
                            self.start_averaging()

                    elif parItem == 'Subtract':
                        self.pedestals = self.enable_subtraction(itemData)

                    elif parItem in ('ROI Vertical','ROI Horizontal'):
                        if itemData:
                            axis = {'ROI Vertical':Y,'ROI Horizontal':X}
                            i = axis[parItem]
                            pxPmm = self.pixelPmm[i]/self.sizeFactor
                            maxPix = self.data.shape[1],self.data.shape[0]
                            self.outsideGraphs[parItem]\
                              = GraphRoiProj(i,pxPmm,maxPix[i])
                            self.update_ROI()
                        else:
                            try:    del self.outsideGraphs[parItem]
                            except: pass
                                          
                    elif parItem == 'ROI Intensity':
                        if itemData:
                            self.create_roiIntensity()
                            self.outsideGraphs[parItem] = \
                              self.roiIntensityPlot
                        else:
                            if self.roiIntensityPlot:
                                self.roiIntensityPlot.win.close()
                                self.roiIntensityPlot = None
                       
                if parGroupName == 'SpotFinder':               
                    if parItem == 'MaxSpots':
                        self.maxSpots = itemData
                    
                if parGroupName == 'Ref Images':
                    if itemData == 'None':
                        return
                    prefix = self.refPath+'_'
                    child = self.pgPar.child(parGroupName).child(parItem)
                    fileName = prefix+itemData+'.png'
                    if parItem == 'View':
                        #cmd = ["gm",'display',prefix+'*']
                        if not os.path.exists(fileName):
                           cprint('file does not exist: '+fileName)
                           return
                        cmd = 'imageViewer.py -M -r -p -o0 -C '+self.cameraName[0]+\
                          ' -bfile '+ fileName
                        #cprint('viewing from '+prefix+'*'+', use Backspace/Space for browsing')
                        cprint('executing: '+cmd)
                        p = subprocess.Popen(cmd.split(), 
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    else:
                        if parItem == 'Store':
                            if itemData == 'ref0':
                                #msg = 'cannot store to '+parItem+', it is reserved for fixed background'
                                #cprinte(msg)
                                if not qMessage('ref0 is reserved for streamed pedestal subtraction. Are you sure to overwrite it?'):
                                    return
                            if os.path.exists(fileName):
                                if not qMessage('Are you sure you want to overwrite '+itemData+'?'):
                                    return
                            codec.save(fileName,self.data)
                            #self.save_image(fileName)
                            cprint('Current image saved to '+fileName)
                        else:
                            self.reference_operation(fileName,parItem)

                if parGroupName == self.addonEntry:
                    self.addon.addon_clicked(parItem,itemData)
 
                if parGroupName in ('For Experts','Control'):
                    if parItem == 'Fast':
                        if itemData:
                            self.set_dockPar('For Experts','ROI Plot',0)
                            self.set_dockPar('For Experts','ContrastCtrl',0)
                            pvMonitor.set_refresh(1000)
                        else:
                            self.set_dockPar('For Experts','ROI Plot',1)
                            self.set_dockPar('For Experts','ContrastCtrl',1)
                            pvMonitor.set_refresh(self.refresh)

                if parGroupName == 'For Experts':
                    if parItem == 'Color':
                        if itemData == 'Gray':
                            pargs.gray = True
                            self.savedData = self.data
                            self.data = rgb2gray(self.data)
                            self.update_imageItem_and_ROI()
                                
                        elif itemData == 'Native':
                            pargs.gray = False
                            self.data = self.savedData
                            self.grayData = rgb2gray(self.data)
                            self.update_imageItem_and_ROI()
                        else:
                            cprintw('Color = '+itemData+' reserved for future updates')
                    elif parItem == 'ContrastCtrl':
                        if itemData:
                            self.gl.addItem(self.contrast)
                        else:
                            self.gl.removeItem(self.contrast)
                    elif parItem == 'ROI Plot':
                        self.roiPlotEnabled = itemData
                        self.plot.clear()
                    elif parItem == 'SpotText':
                        if itemData:
                            self.mainSpotText = mainSpotText()
                            self.mainSpotText.setPos(*self.top_right_corner())
                            self.mainSpotText.setText(self.mainSpotTxt)
                            self.mainWidget.addItem(self.mainSpotText)
                        else:
                            self.mainWidget.removeItem(self.mainSpotText)
                            self.mainSpotText = None
                    elif parItem == 'PointerText':
                        if itemData:
                            self.pointerText = pointerText()
                            self.mainWidget.addItem(self.pointerText)
                        else:
                            if self.pointerText:
                                self.mainWidget.removeItem(self.pointerText)
                            self.pointerText = None                     
                    elif parItem == 'Blur':
                        #TODO: blur color images
                        #self.data = blur(self.data)
                        self.blurWidth = itemData
                        self.update_imageItem_and_ROI()
                    elif parItem == 'FitBaseLo':
                        self.fitBaseBounds[0] = itemData
                        if self.fitBaseBounds[0] <= -1.e9:
                            self.fitBaseBounds[0] = -np.inf
                    elif parItem == 'FitBaseHi':
                        self.fitBaseBounds[1]= itemData
                        if self.fitBaseBounds[1] >= 1.e9:
                            self.fitBaseBounds[1]= np.inf
                    elif parItem == 'Debug':
                        pargs.dbg = itemData
                        printi('Debugging is '+('en' if pargs.dbg else 'dis')+'abled')
                    elif parItem == 'Sleep':
                        self.sleep = itemData
                    elif parItem == 'Debug Action':
                        printi('Debug Action pressed')
                if parGroupName == 'Help':
                    import webbrowser
                    webbrowser.open(pargs.userGuide)
                if parGroupName == 'Exit':
                    print('Exit')
                    if self.addon: 
                        self.addon.exit()
                    for item in self.outsideGraphs:
                        #self.outsideGraphs[item].__del__()
                        self.outsideGraphs[item].win.close()
                    self.win.close()
                    self.exit()
                         
        #QtGui.QSpinBox.setKeyboardTracking(False) # does not work
        self.pgPar.sigTreeStateChanged.connect(handle_change)
           
        # Too lazy for recursion:
        '''
        def valueChanging(param, value):
            printi('Value changing (not finalized):'+str((param, value)))

        for child in self.pgPar.children():
            child.sigValueChanging.connect(valueChanging)
            for ch2 in child.children():
                ch2.sigValueChanging.connect(valueChanging)
        '''        

        def valueChanged(param, value):
            printi('Value changed:'+str((param, value)))

        #printd('>connect(valueChanged)')
        for child in self.pgPar.children():
            child.sigValueChanged.connect(valueChanged)
            #child.setKeyboardTracking(False) # does not work
            for ch2 in child.children():
                if not ch2.readonly:
                    ch2.sigValueChanged.connect(valueChanged)
                    #ch2.setKeyboardTracking(False) # does not work 
                    
        #````````````````````Create docks, place them into the window`````````
        # Note that size arguments are only a suggestion; docks will still 
        # have to fill the entire dock area and obey the limits of their 
        # internal widgets.

        # Create ParameterTree widgets, both accessing the same data
        pgParTree = ParameterTree()
        pgParTree.setParameters(self.pgPar, showTop=False)
        pgParTree.setWindowTitle('Parameter Tree')
        self.registerDock('dockPar',pgParTree,(220,0),'left')

        #set dockImage: a plot area (ViewBox + axes) for displaying the image
        self.gl = pg.GraphicsLayoutWidget()
        self.viewBox = CustomViewBox(name='viewbox')        
        self.mainWidget = self.gl.addPlot(viewBox = self.viewBox)
        self.mainWidget.addItem(self.imageItem)
        self.setup_main_widget()        
        self.mainSpotText = mainSpotText()
        self.mainSpotText.setPos(*self.top_right_corner())
        self.mainWidget.addItem(self.mainSpotText)
        # mouse interaction stuff
        #self.mainWidget.getViewBox().setMouseMode(1)# RectMode
        self.viewBox.setMouseMode(1)# RectMode
        self.pointerText = pointerText()
        self.mainWidget.addItem(self.pointerText)
        self.proxy = pg.SignalProxy(self.mainWidget.scene().sigMouseMoved, 
        rateLimit=60, slot=self.mouseMoved)

        h,w = self.data.shape[:2]        
        imageHSize = float(w)/float(h)*pargs.vertSize
        winHSize = float(w+100)/float(h)*pargs.vertSize # correct for the with of the contrast hist
        self.registerDock('dockImage',self.gl,(imageHSize,pargs.vertSize),'left')

        if pargs.hist:
            # Contrast/color control

            self.contrast = pg.HistogramLUTItem()#fillHistogram=False) #,rgbHistogram=True)
            # fillHistogram=False is less CPU-hungry
            # rgbHistogram is not implemented
            # levelMode='rgba' not implemented

            self.contrast.setImageItem(self.imageItem)
            #print('image max intensity:',self.grayData.max())
            #maxIntensity = max(self.grayData.max(),255)
            # TODO: maxIntensity could be not optimal for high dynamic range images,  
            #maxIntensity = (1<<self.hwpb[3])-1
            #self.contrast.setLevels(0, maxIntensity) #-1 to get rid of possible saturation
            
            self.gl.addItem(self.contrast)
            #printd('<contrast control')

        if pargs.iso != 'Off':
        # Isocurve drawing
            if pargs.iso == 'ROI':
                self.isoInRoi = True
                printw('iso == ROI is not fully functional yet')
            self.iso = pg.IsocurveItem(level=0.8)
            self.iso.setParentItem(self.imageItem)
            self.iso.setZValue(5)
            # Draggable line for setting isocurve level
            self.isoLine = pg.InfiniteLine(angle=0, movable=True, pen='g')
            print('isoLine created')
            self.contrast.vb.addItem(self.isoLine)
            self.contrast.vb.setMouseEnabled(y=False) # makes user interaction a little easier
            self.isoLine.setValue(self.threshold)
            self.isoLine.setZValue(1000) # bring iso line above contrast controls
            # Connect callback to signal
            #self.isoLine.sigDragged.connect(self.update_isocurve)
            self.isoLine.sigPositionChangeFinished.connect(self.update_isocurve)
            #self.update_iso()

        if pargs.roi:
        # Custom ROI for selecting an image region
            self.plot = pg.PlotWidget()
            self.plot.showGrid(True,True)
            self.registerDock('dockPlot',self.plot,(1,100),'bottom')

            h,w = self.data.shape[:2]
            rect = [i*g for i,g in zip((w,h,w,h),self.roiRect)]
            self.roi = pg.RectROI(rect[:2], rect[2:], sideScalers=True, movable=False)
            self.roi.addScaleHandle([0, 0],[1, 1])
            self.mainWidget.addItem(self.roi)
            self.roi.setZValue(10)  # make sure pargs.roi is drawn above image
            
            # create max number of spot labels
            #self.spotLabels = [pg.TextItem('*',color='r',anchor=(0.5,0.5)) 
            self.spotLabels = [pg.TextItem('*',color=self.marksColor,
              anchor=(0.5,0.5)) for i in range(pargs.maxSpots)]
            for sl in self.spotLabels:
                self.mainWidget.addItem(sl)

            # Connect callback to signal
            self.roi.sigRegionChangeFinished.connect(self.cb_update_ROI)
            #printd('<roi')

        if pargs.console:
        # interactive python console
            global gWidgetConsole
    
            def sh(s): # console-available metod to execute shell commands
                print(subprocess.Popen(s,shell=True, stdout = None if s[-1:]=="&" else subprocess.PIPE).stdout.read())

            gWidgetConsole = CustomConsoleWidget(
                namespace={'pg': pg, 'np': np, 'plot': self.plot, 
                'roi':self.roi, 'data':self.data, 'image': self.qimg, 
                'imageItem':self.imageItem, 'pargs':pargs, 'sh':sh},
                historyFile='/tmp/pygpm_console.pcl',text="")
            self.widgetConsole = gWidgetConsole # could be needed for addons
            self.registerDock('dockConsole',gWidgetConsole,(0,10),'bottom')
            self.docks['dockConsole'][0].setStretch(0,0)
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

        self.win.resize(winHSize, pargs.vertSize)
        self.win.show()
        
        # set color map to grayClip: red color for saturated pixels
        if len(self.data.shape) == 2:
            cmap_grayClip = pg.ColorMap(pos=[0.,0.995,1.],color=
              [[0,0,0,255],[255,255,255,255],[255,0,0,255]])
            lut_grayClip = cmap_grayClip.getLookupTable()
            self.imageItem.setLookupTable(lut_grayClip)
        else:
            msg = 'Avoiding colormapping for color images'
            printw(msg)
            #self.stop()
            #raise Exception('ERROR '+msg)
        self.hideDocks(pargs.miniPanes)
        #print('<iv_show')
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

    def save_refring(self):
        try:
            import Sybase
            mdb = Sybase.Connection('OPSYB1', 'ops', 'opsops', 'appConfig')
            mcrs = mdb.cursor()
            def set_camera_DB(db,cursor,row,col,value):
                updateCmd = 'UPDATE cameras SET '+col+' = '+str(value)\
                  +' WHERE camName = "'+row+'"'
                #print('executing DB',updateCmd)
                cursor.execute(updateCmd)
                db.commit()
            cols = {'ref_X':self.ref_X,
                    'ref_Y':-self.ref_Y,# DB (0,0) is top left corner
                    #TODO:add pix/mm for X and Y
                    'target':self.ref_diameter/self.pixelPmm[1]}
            if not qMessage('Are you sure to change DataBase %s?'%\
              ('appConfig/'+self.cameraName[0])):
                                    return
            for c,v in cols.items():
                set_camera_DB(mdb,mcrs,self.cameraName[0],c,v)
            cprint('DataBase updated successfully')
        except Exception as e:
            cprint('in update DB :'+str(e))
        return
            
    def hideDocks(self,okToHide):
        for name,dock_size in self.docks.items():
            stretch = (0,0) if okToHide else dock_size[1]
            dock_size[0].setStretch(*stretch)

    def registerDock(self,name,widget,size,side):
        dock = pg.dockarea.Dock(name, size=size)
        self.docks[name] = [dock,size]
        dock.addWidget(widget)
        self.area.addDock(dock,side)

    def set_axes_scale(self,mm_pix):
        self.axesInMm = mm_pix
        self.redraw_axes()

    def roiIntensity(self):
        #y,x = np.histogram(self.roiArray.flatten(),bins=self.roiArray.max()-1)
        a = self.roiArray
        m = a.max()
        #mPos = np.unravel_index(np.argmax(a), a.shape)
        #print('max at ',mPos,m)
        y,x = np.histogram(a.flatten(),bins=np.linspace(0,m+1,m+2))
        return x,y
        
    def plot_roiIntensity(self):
        self.roiIntensityPlot.plot(*self.roiIntensity(),stepMode=True,pen='b',clear=True)
        mean = self.roiArray.mean()
        std = self.roiArray.std()
        self.roiIntensityPlot.setLabel('top','Mean:%.2f'%mean+', std:%.2f'%std)

    def create_roiIntensity(self):
        self.roiIntensityPlot = pg.plot(*self.roiIntensity(),stepMode=True,pen=None)
        self.roiIntensityPlot.win.resize(480,320)
        self.roiIntensityPlot.win.setWindowTitle('Pixel Amplitudes in ROI')
        self.roiIntensityPlot.showGrid(True,True)
        self.roiIntensityPlot.setLogMode(y=True)
        self.plot_roiIntensity()
        
    def set_dockPar(self,child,grandchild,value):
        self.pgPar.child(child).child(grandchild).setValue(value)

    def closeEvent(self, event):
        print("Closing")

    def setup_main_widget(self):
        #self.mainWidget.autoRange(padding=None) # remove default padding
        ##self.mainWidget.setrange(0.,w)
        ##self.mainWidget.setRange(range=(0,w),yRange=(0,h),padding=None)
        #self.mainWidget.setAutoVisible(x=True,y=True) # better use of graph area
        #self.mainWidget.setLimits(xMin=0.,yMin=0.) # limit autopan 
        self.mainWidget.setAspectLocked()
        h,w = self.data.shape[:2]
        self.mainWidget.setRange(rect=pg.QtCore.QRect(0,0,w,h),padding=None)
        
    def redraw_axes(self):
        self.needToRedrawAxes = False
        axx = self.mainWidget.getAxis('bottom')
        axy = self.mainWidget.getAxis('left')
        rx,ry = axx.range, axy.range
        #print('redraw_axes',axx.range, axy.range)
        if self.axesInMm:
            rx[0] = self.pix2mm(rx[0],0.)[0]
            rx[1] = self.pix2mm(rx[1],0.)[0]
            ry[0] = self.pix2mm(0.,ry[0])[1]
            ry[1] = self.pix2mm(0.,ry[1])[1]
            axx.setRange(*rx)
            axy.setRange(*ry)
            #print('axx',axx.range)
        else:
            # Not a best method to revert axes to pixels
            # by slightly resizing it...
            size = self.win.size()
            self.win.resize(size.width()+2,size.height())

    #v156#def set_calibs(self, pixPerMM=(10.,10.), ringDiameter=[100,100],
    def set_calibs(self, pixPerMM=10., ringDiameter=100.,
                    ringCenter=(0,0), xScale=1.):
        self.ref_diameter = ringDiameter
        self.pixelPmm = [round(pixPerMM*xScale,3), round(pixPerMM,3)]
        self.ref_X,self.ref_Y = ringCenter
        self.pixelXScale = xScale
        #if self.pixelXScale != 1.:
        #    printi('pixelXScale != 1: %.3f'%self.pixelXScale)

    def create_refRing(self):
        #color = 'b'
        if self.refRing is not None:
            #printi('RefRing exists')
            return
        color = (199,234,70) #lime
        pen = pg.mkPen(color=color, width=1, style=QtCore.Qt.DashLine)
        h,w = self.data.shape[:2]
        x = self.ref_X/self.sizeFactor + w/2.
        y = self.ref_Y/self.sizeFactor + h/2.
        dx = self.ref_diameter/self.sizeFactor
        dy = dx
        refring = pg.QtGui.QGraphicsEllipseItem(x-dx/2.,y-dy/2.,dx,dy)
        cross0 = pg.QtGui.QGraphicsLineItem(x-dx/2.,y,x+dx/2.,y)
        cross1 = pg.QtGui.QGraphicsLineItem(x,y-dy/2.,x,y+dy/2.)
        self.refRing = (refring,cross0,cross1)
        for i in self.refRing:
            i.setPen(pen)
            self.mainWidget.addItem(i)

    def remove_ref_ring(self):
        if self.refRing is not None:
            for i in self.refRing:
                self.mainWidget.removeItem(i)
            self.refRing = None
        
    def enable_saving(self,trueFalse):
        """Enable/disable saving."""
        self.saving = trueFalse
        
    def save_image(self,fileName=None):
        if fileName is None:
            fmt = '_%Y%m%d_%H%M%S%f.png'
            if self.timestamp:
                strtime = datetime.datetime.fromtimestamp(self.timestamp).strftime(fmt)
            else:
                strtime = time.strftime(fmt) # note, time.strftime does not handle the %f format
            fileName = self.savePath+'IV_'+self.cameraName[0]+strtime
            self.imageFile = fileName
        try:
            ts = timer()
            if pargs.fullsize:
                codec.save(fileName)
            else: # PNG image format
                printd('writing original data')
                with open(fileName,'wb') as f:
                    f.write(pvMonitor.blob)
            try:
                r = pvMonitor.cadDev.set('img.'+self.cameraName[0],
                  'lastSavedM',fileName)
            except Exception as e:
                printe('exception:'+str(e))
                r = 1
            if (r):
                printe('updating lastSavedM')
            printd('save_image time %.3f'%(timer()-ts)+'s for %.3fkB'%(len(self.data)/1000.)) 
        except Exception as e:
            printe('save_image exception: '+str(e)+'\n'+traceback.format_exc())

    def shape_data(self,data):
        if data is not None:
            if len(self.data.shape) > 2:
                data = data[:,:,0] # saved files are always color PNG
            data = data[::-1]# flip vertically
        return data
        
    def adjust_refData(self,data):
        try:
            if data.shape != self.data.shape:
                printw('reference image shape '+str(data.shape)+' != '\
                  +str(self.data.shape))
                h,w = self.data.shape
                if data.shape[1] > w:
                    data = data[:,:w]
                print('reference image reshaped '+str(data.shape))
            return data
        except:
            raise NameError('')

    def reference_operation(self,fileName,operation):
        """ binary operation with current and restored image """
        try: #if True:
            data = self.shape_data(codec.load(fileName))
            data = self.adjust_refData(data)
            if data is None:
                raise NameError(operation)
            print('refop',data.shape,self.data.shape)
            if operation == 'Retrieve':
                self.data = data
                cprint('retrieved image '+fileName)
            elif operation == 'Add':
                self.data = (self.data.astype(int) + data)/2
                cprint('added image '+fileName+' to current image')
            elif operation == 'Subtract':
                self.data = self.data.astype(int) - data
                cprint('subtracted image '+fileName+' from current image')
            else: pass
            self.grayData = rgb2gray(self.data)
            self.update_imageItem_and_ROI()
        except Exception as e:
            msg = 'in ref_op '+operation+': '+str(e)
            cprinte(msg)

    def enable_subtraction(self,ref):
        if ref == 'None':
            return None
        fn = self.refPath+'_'+ref+'.png'
        try:
            d = codec.load(fn)
            if d is None:
                raise NameError('')
            sd = self.shape_data(d)
            #print('subdata',self.data.shape)
            r = self.adjust_refData(sd)
            cprint('Subtraction of '+ref+' enabled')
            return r
        except Exception as e:
            cprinte('cannot load '+ref+':'+str(e))
            self.set_dockPar('Analysis','Subtract','None')# this has no effect
            return None
    
    def stop(self):
        self.stopProcThread = True
        pvMonitor.clear()
        printi('imager stopped')
        try: self.addon.stop()
        except: pass
        try: self.spotLog.close()
        except: pass

    def mm_per_pixel(self,xPx,yPx):
        # scaling of pixels to mm (usually used for converting widths).
        g = [self.sizeFactor/r for r in self.pixelPmm]
        #x = round(xPx*g, Prec)
        #y = round(yPx*g, Prec)
        x = xPx*g[0]
        y = yPx*g[1]
        return x,y
        
    def pix2mm(self,xPx,yPx):
        # pixel to mm conversion, zero mm is at the center of image
        h,w = self.data.shape[:2]
        return self.mm_per_pixel(xPx - w/2.,yPx - h/2.)
        
    def standard_analysis(self):
            #self.imageItem.setBorder(self.borderPen)
            ts = timer()
            self.gl.setBackground(self.viewBoxBackgroundColor)
            ox,oy = self.roiOrigin
            self.spots = find_spots(self.roiArray,self.threshold,
              self.maxSpots, pargs.finalFit, fitBaseBounds=self.fitBaseBounds)
            #print('spots',self.spots)
            profile('spotsFound')
            if pargs.profile:
                #print(profStates('initRoi','spotsFound'))
                print('FindSpot time: '+profDif('initRoi','spotsFound'))
            if self.sizeFactor == 0.:
                # should not be here
                printw('TODO logic error sizeFactor=0')
                return
            spen = pg.mkPen(self.marksColor)
            if self.spotLog: 
                logdata = OrderedDict([('time:',time.strftime('%y-%m-%d %H:%M:%S'))])
            self.stdPars = []
            for spotIdx,spot in enumerate(self.spots):
                p,wPx,pcc,psum,fittedPars = spot[:5]
                #print('fp',spotIdx,len(fittedPars))
                posPx = p + (ox,oy)
                
                theta,sigmaU,sigmaV = ellipse_prob05(wPx[0],wPx[1],pcc)
                # store the spot parameters, they could be useful for user analysis
                if wPx[0]*wPx[1] > 4.: # do not store too narrow spots 
                    tilt = theta if sigmaU > sigmaV else math.pi/2. + theta
                    self.stdPars.append([posPx,wPx,pcc,psum, # primary parameters
                      tilt,(sigmaU,sigmaV)]) # derivative parameters
                    #print('u,v,th,tilt',sigmaU,sigmaV,theta,tilt,wPx)
                # record the results
                if self.spotLog: 
                    if self.axesInMm:
                        posLog = self.pix2mm(*posPx)
                        wLog = self.mm_per_pixel(*wPx)
                    else:
                        posLog = posPx
                        wLog = wPx
                    #round and log the parameters
                    posLog = map(lambda z: round(z,Prec), posLog)
                    wLog = map(lambda z: round(z,Prec), wLog)
                    logdata['spot_'+str(spotIdx)] = (list(posLog),list(wLog),
                      round(pcc,Prec),round(psum,1))
                
                # plot ellipse and text labels with position and widths
                if not self.cleanImage:

                    # position the spot mark
                    self.spotLabels[spotIdx].setPos(*posPx)

                    too_small = max(wPx) < 0.1 # 0.1 is big enough for ellipse
                    if not too_small:                    
                        # draw ellipse
                        spotShape = pg.QtGui.QGraphicsEllipseItem(posPx[0]-sigmaU,
                          posPx[1]-sigmaV,sigmaU*2.,sigmaV*2.)
                        spotShape.setTransformOriginPoint(*posPx)
                        if theta:
                            spotShape.setRotation(np.rad2deg(theta))
                        spotShape.setPen(spen)
                        self.mainWidget.addItem(spotShape)
                        self.spotShapes.append(spotShape)

                    # draw mainSpotText
                    if self.mainSpotText is not None and spotIdx == 0:
                        if self.axesInMm:
                            posTxt = self.pix2mm(*posPx)
                            w = self.mm_per_pixel(*wPx)
                        else:
                            posTxt = tuple(posPx)
                            w = wPx
                        self.mainSpotTxt = 'posXY:(%.1f,%.1f)'%posTxt
                        if not too_small:
                            st = ' sigma' if len(fittedPars)>0\
                              else ' stDev'
                            self.mainSpotTxt +=\
                              st+'XY:(%.2f,%.2f)'%(abs(w[0]),abs(w[1]))
                        self.mainSpotText.setText(self.mainSpotTxt)
                        
                        # if text width too large, plot empty text
                        sbr = self.mainSpotText.boundingRect()
                        ibr = self.imageItem.boundingRect()
                        #scale them into viewbox data dimensions
                        sx, sy = self.viewBox.viewPixelSize() #scaling factors
                        if sx*sbr.width() > ibr.width():
                            self.mainSpotText.setText('')  
                                                                        
            # reset outstanding spotLabels
            for j in range(len(self.spots),len(self.spotLabels)):
                self.spotLabels[j].setPos(0,0)
            self.set_dockPar('SpotFinder','Found:',len(self.spots))
            #self.set_dockPar('SpotFinder','Spots',msg)
            if self.spotLog:
                #print('ld:',logdata)
                if self.saving:
                    #print('if:'+self.imageFile)
                    logdata['file'] = '_'.join(self.imageFile.split('_')[-2:])
                self.spotLog.write(dumps(logdata))#,indent=2))
                self.spotLog.write('\n')
                #printd('written:'+str(logdata))
                
    def remove_spotShapes(self):
        # remove previously found spots
        for s in self.spotShapes:
            self.mainWidget.removeItem(s)
        self.spotShapes = []
        
        # move labels to 0
        for i in range(len(self.spotLabels)):
            self.spotLabels[i].setPos(0,0)
                
    def cb_update_ROI(self):
        """callback for handling changed ROI"""
        self.update_ROI(report=True)
        if self.imageMan is not None:
            self.update_imageMan()
        if self.zoomROI:
            self.setup_zoomROI()

    def update_ROI(self,report=False):
        """handle changed ROI"""
        profile('initRoi')
        # the following is much faster than getArrayRegion
        slices = self.roi.getArraySlice(self.grayData,self.imageItem)[0][:2]
        self.roiOrigin = slices[1].start, slices[0].start
        self.roiArray = self.grayData[slices]
        self.colorRoiArray = self.data[slices]
        
        if self.blurWidth > 0.:
            #print('roi max pre:',self.roiArray.max())
            self.roiArray = blur(self.roiArray,self.blurWidth)
            self.grayData[slices] = self.roiArray
            #self.data[slices] = self.roiArray
            #print('roi max post:',self.roiArray.max())
        
        if report:
            relativeRect = [float(v1)/float(v2) 
              for v1,v2 in zip(self.roiOrigin,self.data.shape[::-1])]
            relativeRect += [float(v1)/float(v2) 
              for v1,v2 in zip(self.roiArray.shape,self.data.shape)][::-1]
            cprint('ROI changed: '+str((self.roiOrigin,self.roiArray.shape,
              ['%0.2f'%i for i in relativeRect])))
            self.roiRect = map(lambda z: round(z,3), relativeRect)

        self.remove_spotShapes()
        profile('roiArray')
        
        # find spots using isoLevel as a threshold
        if self.threshold>1 and self.maxSpots>0:
            self.standard_analysis()
            profile('standard processing')
            if self.addon:
                self.addon.process()
                profile('addon processing')
        
        if not self.roiPlotEnabled:
            return
        # plot the ROI histograms
        meansV = self.roiArray.mean(axis=X) # vertical means
        x = np.arange(len(meansV)); s = False
        pen = 'k' if not pargs.black else 'w'
        if len(self.data.shape) == 2: # gray image
            self.plot.plot(x,meansV,clear=True,pen=pen,stepMode=s)
        else: # color image
            # plot color intensities
            cMeansV = self.colorRoiArray.mean(axis=X)
            self.plot.plot(x,cMeansV[:,0],pen='r', clear=True,stepMode=s)#plot red
            self.plot.plot(x,cMeansV[:,1],pen='g',stepMode=s) # plot green
            self.plot.plot(x,cMeansV[:,2],pen='b',stepMode=s) # plot blue
            self.plot.plot(x,meansV,pen=pen,stepMode=s)
        fitPars = []
        if len(self.spots) > 0:
            if len(self.spots[0][4]):
                try: 
                    fitPars = self.spots[0][4][X], self.spots[0][4][Y]
                    #print('fp',fitPars)
                    #print(self.roiArray.shape[0])
                    il,ir = self.spots[0][5][X]
                    fitY = func_sum_of_peaks(x[il:ir],*fitPars[X])/self.roiArray.shape[0]
                    #print(fitY)
                    #print(x[il:ir])
                    self.plot.plot(x[il:ir],fitY[0:ir-il],pen=pg.mkPen(color=(0, 0, 255), style=QtCore.Qt.DashLine),stepMode=s)
                except Exception as e:
                    printe('no fit data ')
                    pass
        if self.roiIntensityPlot:       
            self.plot_roiIntensity()
        #if self.outsideGraphs['roiVerticalPlot']:        self.plot_roiVertical(fitParsY)
        for name,graph in self.outsideGraphs.items():
            try:
                graph.update(self.roiArray,self.roiOrigin,fitPars)
            except Exception as e:
                pass
                #TODO, divide by zero encountered in log10 for intensity plot
                #printw('exception in graph.update() for '+str(name))
        #if self.blurWidth > 0.:
        #    self.imageItem.setImage(self.data)        
        profile('roiPlot')
        
    def setup_zoomROI(self):
        x0,y0 = self.roiOrigin[0],self.roiOrigin[1]
        wy,wx = self.roiArray.shape
        self.mainWidget.setRange(rect=pg.QtCore.QRect(x0,y0,wx,wy))

    def update_isocurve(self):
    # callback for handling ISO
        print('update_isocurve')
        self.show_isocurves(True)
        #profile('init iso')
        v = self.isoLine.value()
        # inform imager on changed threshold
        self.set_dockPar('Control','Threshold',v)
         
    def show_isocurves(self,show=True): 
        #print('show_isocurves',show)
        self.isocurve_enabled = show
        if show:
            self.iso.setPen(pg.mkPen('g'))
        else:
            try:
                self.iso.setPen(pg.mkPen(None))
            except Exception as e:
                printe('in show_isocurve:'+str(e))
                pass
        
    def mainWidget_changed(self):
        # redrawing axes here have no effect as it will be overridden
        # we need to postpone the activation somehow.
        # simple delayed thread is not working.
        if self.axesInMm:
            self.needToRedrawAxes = True
        #print('need to redraw axes: ',self.needToRedrawAxes)
    
    def update_imageItem_and_ROI(self):
        #DNW#self.imageItem.setImage(self.data,autoLevels=pargs.normalize)
        if pargs.normalize:
            if self.contrast:
                mx = pargs.contrastLevel
                if mx == 0:
                    #mx = self.data.max() # too noisy
                    # very fast, 5ms/Mbyte:
                    mx = np.convolve(self.data.flatten(), np.ones((4,))/4, mode='same').max()
                self.contrast.setLevels(0,mx)
                self.contrast.setHistogramRange(0, mx)
        if self.roi:
            self.update_ROI()
        self.imageItem.setImage(self.data)
        if self.contrast is not None: self.contrast.regionChanged() # update contrast histogram
        profile('setImage')

    def process_image(self):
        global profilingState
        if not pvMonitor:
            printe('pvMonitor did not start, processing stopped')
            sys.exit(10)        
        profile('start process')
        
        data,self.timestamp = pvMonitor.get_data_and_timestamp()
        #print('d,t:'+str((len(data),self.timestamp)))
        profile('got image array')
        dataLen = len(data)
        if dataLen == 0:
            if self.timestamp == 0:
                print('No more data')
                if self.imageItem is None:# first event was not processed
                    printe('No valid images were found')
                    sys.exit(1)
                if self.backend == 'file':
                    self.set_pause(True)
            if not self.paused:
                EventPocessingFinished.set()
            return
        #print('d,ts:',data.shape,self.timestamp)
        #print('data from monitor:'+str(data.dtype)+':\n'+str(data[:100]))
        if pargs.width:
            #``````````The source is vector parameter with user-supplied shape
            #print('dl',dataLen,self.dataLen)
            if self.dataLen != dataLen: # data size changed
                print('size changed',self.dataLen,dataLen)
                pargs.width = pvMonitor.get_image_shape()
                tokens = pargs.width.split(',')
                w = int(tokens[0])
                h = int(tokens[1])
                bytesPerPixel = dataLen/w/h
                try:    bitsPerPixel = int(tokens[2])
                except: bitsPerPixel = 8*bytesPerPixel                
                # we cannot decide exactly about nPlanes and bytesPerChannel based on bitsPerPixel
                # here is assumption:
                nPlanes = 3 if bitsPerPixel > 16 else 1 #
                self.hwpb = [h, w, nPlanes,bitsPerPixel]
                #print('hwpb',self.hwpb)
                if self.roi:
                    #print('roi',self.roi)
                    self.roi.setPos(self.roiRect[0]*w,self.roiRect[1]*h)
                    self.roi.setSize(self.roiRect[2]*w,self.roiRect[3]*h)

            bytesPerChannel = ((self.hwpb[3]-1)/8)+1
            #print('bytesPerChannel:',bytesPerChannel,self.hwpb)
            if bytesPerChannel > 1: # we need to merge pairs of bytes to integers
                #data = np.array(struct.unpack('<'+str(dataLen/bytesPerChannel )+'H', data.data),'u2')
                #data = struct.unpack('<'+str(dataLen/bytesPerChannel )+'H', bytearray(data)) 
                data = struct.unpack('<'+str(dataLen/bytesPerChannel )+'H', data.data) 
                profile('merge')
            shape = (self.hwpb[:3] if self.hwpb[2]>1 else self.hwpb[:2])
            data = np.reshape(data,shape).astype('u2')
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        #````````````````````Got numpy array from data````````````````````````
        self.dataLen = dataLen # save data length for monitoring        
        data = data[::-1,...] # flip vertical axis
        data = crop_width(data) # correct the data width to be divisible by 4
        self.receivedData = data # store received data
        #printd('array: '+str((data.shape,data.dtype,data)))
        profile('array')
                                        
        self.data = rotate(self.receivedData,self.dockParRotate)
        #profile('rotate'); print('rotation time:'+profDif('array','rotate'))
        #print('after rotation:',self.data.shape)

        if pargs.gray:
            self.data = rgb2gray(self.data)
            self.grayData = self.data
        else:
            self.grayData = rgb2gray(self.data)
        profile('gray conversion')

        #````````````````````Data array is ready for analisys`````````````````
        h,w = self.data.shape[:2]
        if self.imageItem is None: 
            #````````````````First event, do the show() only once`````````````
            if self.hwpb[0] == 0: # get p,b: number of planes and bits/channel
                try: p = self.data.shape[2]
                except: p = 1
                b = self.data.dtype.itemsize*8
                self.hwpb = [h,w,p,b]
            #printd('hwpb:'+str(self.hwpb))
            #printd('self.array: '+str((self.data.shape,self.data)))
            if self.averageWindow:
                self.start_averaging()
            self.imageItem = pg.ImageItem(self.data)
            self.iv_show()
        else:  #`````````````udate data```````````````````````````````````````
            pass
        if True:
            #TODO: react on shape change correctly, cannot rely on self.hwpb because of possible rotation
            if self.saving: self.save_image() #TODO shouldn't it be after update
            if self.backend == 'file':
                try:    self.winTitle = pvMonitor.pvname.rsplit('/',1)[1]
                except: self.winTitle = pvMonitor.pvname
            else:
                self.winTitle = self.pvname+time.strftime('_%H%M%S')
            self.win.setWindowTitle(self.winTitle)

            if self.pedestals is not None:
                try:
                    self.data = self.data.astype(int) - self.pedestals
                    self.grayData = rgb2gray(self.data)
                except:
                    cprinte('Reference is not compatible, subtraction disabled')
                    self.pedestals = None
                    self.set_dockPar('Analysis','Subtract','None')

            if self.averageWindow:
                l = len(self.averageQueue)
                #print('len(a)',l,self.average.shape)
                dtype = self.data.dtype
                self.average += self.data
                if l >= self.averageWindow:
                    pl = self.averageQueue.popleft()
                    self.average -= pl
                else: l += 1
                self.averageQueue.append(self.data.astype(dtype))
                self.data = self.average/l
                self.grayData = rgb2gray(self.data)
            
        self.update_imageItem_and_ROI()
        if self.events % 100 == 0:
            try:
                dt = timer() - profilingState['100 events']
                #print('Performance: %.1f f/s'%(100./dt))
                self.set_dockPar('Control','Images/s',100./dt)
            except: pass
            profile('100 events')
        self.events += 1
        self.set_dockPar('Control','Image#',str(self.events)) # concern: time=0.5ms
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,            
        if pargs.profile:
            print('### Total time: '+'%0.3g' % (timer()-profilingState['start']))
            print(profStates('start','finish'))
        if self.needToRedrawAxes:
            self.redraw_axes()
        if not self.paused:
            EventPocessingFinished.set()
            
    def set_pause(self,pause=True):
        """Pause imager, used mainly in addons"""
        cprint('Setting pause '+str(pause))
        self.paused = pause
        self.set_dockPar('Control','Pause',pause)
        if not pause: # wake up processing
            EventPocessingFinished.set()
            
    def mouseMoved(self,evt):
        if self.pointerText is None:
            return
        absPos = evt[0]  ## using signal proxy turns original arguments into a tuple
        #print('pos',absPos)
        #print('vrect',vr,top,left)
        #print('sbr',self.mainWidget.sceneBoundingRect())
        #print('geom',self.mainWidget.viewGeometry())
        vrange = self.mainWidget.viewRange()
        #print('vrange',vrange)
        left = max(0,vrange[0][0])
        top = min(self.data.shape[0],vrange[1][1])
        txt = ''
        if self.mainWidget.sceneBoundingRect().contains(absPos):
            mousePoint = self.mainWidget.vb.mapSceneToView(absPos)
            pos = mousePoint.x(), mousePoint.y()
            ipos = int(pos[1]),int(pos[0])
            pos = self.pix2mm(*pos) if self.axesInMm else ipos
            try:
                intensity = np.mean(self.data[ipos])
                txt = 'x:%.1f y:%.1f'%pos+' i:%d'%intensity
            except:
                pass
        self.pointerText.setPos(left,top)
        self.pointerText.setText(txt)
            
    def top_right_corner(self):
        vrange = self.mainWidget.viewRange()
        right = min(self.data.shape[1],vrange[0][1])
        top = min(self.data.shape[0],vrange[1][1])
        return right,top

    def start_averaging(self):
        from collections import deque
        h,w = self.data.shape
        self.average = np.zeros(h*w).astype(self.data.dtype).reshape(h,w)
        self.averageQueue = deque()
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#````````````````````````````Addon template```````````````````````````````````
class AddonBase():
    def __init__(self,imager):
        """Constructor, inherits the Imager object imager"""
        #self.imager = imager # that should be executed in the derived class
            
    def cprint(self,msg):
        """Print in imageViewer console dock"""
        #print('cprint',msg)
        '''The following is not thread safe
        with cprintLock:
            try:
                self.imager.widgetConsole.write('#'+time.strftime('%H:%M:%S:')+msg+'\n')
            except:
                IV.printe('_cprint not functional')
        '''
        # Thread-safe printing
        self.imager.SignalCPrint.emit('_'+time.strftime('%H:%M:%S: ')+msg)

    def entryName(self):
        """Group name in the parameter tree"""
        return 'Addon'

    def controlPane(self):
        """Define user items in the control pane. For example:
        return [\
            {'name':'Dist', 'type': 'int', 'value':2.0,
                  'tip': 'Distance'},
            {'name':'Where My Results?', 'type': 'button', 'value':False},
            ]
       """
    def addon_clicked(self,parItem,itemData):
        """Called when user item in control pane is clicked. For example:
        if DBG: print('Addon clicked',parItem,itemData)
        if parItem == 'Where My Results?':
            w = IV.QtGui.QWidget()
            msg = 'The results are in /tmp/'
            IV.QtGui.QMessageBox.information(w,'Message',msg)
        """
    def process(self):
        """Called when data region is updated"""
        print('Empty addon.process()')
              
    def start_saving(self):
        """Called when the saving is enabled in the imageViewer"""
        print('Obsolete addon.start_saving()')
              
    def stop_saving(self):
        """Called when the saving is disabled in the imageViewer"""
        
    def stop(self):
        """Called when imager is stopped"""
        print('Empty addon.stop()')
        
    def exit(self):
        print('Empty addon.exit()')       
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#````````````````````````````Main Program`````````````````````````````````````
pvMonitor = None
Backend = None
#png = None
import importlib
def main():
    global pargs, imager, pvMonitor, Backend, codec#, png
    #`````````````````````````````````````````````````````````````````````````
    # parse program arguments
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(
      description = 'Interactive analysis of images from ADO or file',
      formatter_class=RawTextHelpFormatter)
    parser.add_argument('-a','--refreshRate',type=float,
      help='Refresh rate [Hz]')
    parser.add_argument('-A','--average',type=int,default=0,
      help='Averaging, number of images to average')
    parser.add_argument('-B','--black', action='store_true', help=
      'Black background, white foreground for all graphics')
    parser.add_argument('-b','--backend', default = 'ado', help=
      'Data access backend: file/ado/epics/http')
    parser.add_argument('-c','--console', action='store_false', help=
      'Disable interactive python console')
    parser.add_argument('-C','--cameraName', default=None, help=
      'Camera name, used mainly for backend = file when camera name is not recognized from the filename')
    parser.add_argument('-d','--dbg', action='store_true', help=
      'Turn on debugging')
    parser.add_argument('-D','--noDatabase', action='store_true', help=
      'Dont use database')
    #parser.add_argument('-e','--extract',default='png',help=
    #  'Image extractor: qt for QT (default), png for pyPng (for 16-bit+ images)') #cv for OpenCV, 
    parser.add_argument('-F','--flip', help=
      "Flip image, 'V':vertically, 'H':horizontally")
    parser.add_argument('-f','--fullsize', action='store_true', help=
      'Use full-size full-speed imageM parameter')
    parser.add_argument('-g','--gray', action='store_false', help=
      'Set it for non-gray images, not functional in latest versions')
    #parser.add_argument('-G','--graysum', action='store_false', help=
    #  'Use perceptional color-to-gray conversion, rather than uniform')
    parser.add_argument('-H','--hist', action='store_false', help=
      'Disable histogram with contrast and isocurve contol')
    #parser.add_argument('-i','--iso', action='store_false', help=
    #  'Disable Isocurve drawing')
    parser.add_argument('-i','--iso',default='Image',help=
      '''Isocurve drawing options: ROI - only in ROI (default), 
      Image - in full image, 
      Off - no isocurves''')
    parser.add_argument('-j','--contrastLevel', type=int, default=0,help=
      'Contrast control level, 0 for auto.')    
    #parser.add_argument('-l','--logdir', default = '/operations/app_store/imageViewer/',
    logdir = '/operations/app_store/cameras/run_fy19/'
    parser.add_argument('-l','--logdir', default = logdir,
      help='When DB is not available, this defines a directory for saving images, logging and references, i.e: '+logdir)
    parser.add_argument('-m','--maxSpots',type=int,default=4,
      help='Maximum number of spots to find')
    parser.add_argument('-M','--miniPanes', action='store_true', help=
      'Start with minimized panes (i.e. control, histogram and console panes')
    defaultROI = '0.05,0.05,0.9,0.9'
    parser.add_argument('-n','--normalize', action='store_false', help=
      'Disable intensity normalization')
    parser.add_argument('-O','--roiRect',default=defaultROI,
      help='ROI rectangle: posX,pozY,sizeX,sizeY, i.e -O0.05,0.05,0.9,0.9')
    parser.add_argument('-o','--orientation', type=int, default=None, 
      help='''Specifies the order in which pixels are drawn,
to adjust for camera mounting orientation, mirrors, etc.
0 = As captured
1 = Flip Vertical
2 = Flip Horizontal
3 = Flip Vertical & Horizontal
4 = Rotate 90 deg clockwise
5 = Rotate & Flip Vertical
6 = Rotate & Flip Horizontal
7 = Rotate, Flip Horizontal & Vertical''')
    parser.add_argument('-P','--profile', action='store_true', help=
      'Enable code profiling')
    parser.add_argument('-p','--pause', action='store_true', help=
      'Start in paused state')
    parser.add_argument('--prefix',default='IV_')
    parser.add_argument('-R','--rotate', type=float, default=0, help=
      'Rotate image by ROTATE degree')
    parser.add_argument('-r','--roi', action='store_false', help=
      'Disable Region Of Interest analysis')
    parser.add_argument('--refRing', action='store_true', help=
      'Show reference ring')
    parser.add_argument('-s','--saving', action='store_true', help=
      'Enable image saving with evry new image.')
    parser.add_argument('-S','--controlSystem',
      help='Control system, when DB is not available this defines the path for logging directory  i.e TEST or CEC, or LEReC')
    parser.add_argument('--sensorShape',default=None,
      help='Sensor shape, width,height')
    parser.add_argument('-t','--threshold',type=float,default=0.,
      help='Threshold for spot finding')
    parser.add_argument('-T','--finalFit',action='store_true',
      help='Do not fit the brightest spot for faster performance. The calculated spot witdth will be less accurate and more sensitive to the threshold setting')
    parser.add_argument('-v','--vertSize',type=float,default=800,
      help='Vertical size of the display window')
    parser.add_argument('-u','--userAddon', help=
      '''User addon, the addon script name should be prefixed with ivAddon_,
      i.e -uNewAddon will try to import ivAddon_NewAddon.py''')
    parser.add_argument('-U','--userGuide', help=
      'Location of the user instructions', default=
      'http://www.cadops.bnl.gov/Controls/ControlsWiki/index.php/ImageViewer')
    parser.add_argument('-w','--width', help=
      '''For blob data: width,height,bits/pixel i.e 1620,1220,12. 
      The bits/pixel may be omitted for standard images''')
    #parser.add_argument('-W','--white', action='store_true', help=
    #  'White background, black foreground for all graphics')
    #parser.add_argument('-x','--pixPerMM',type=float,
    #  help='Force pixel/mm conversion with applied value.')
    parser.add_argument('-x','--expert',action='store_true', help=
      'Expert/Admin mode')
    parser.add_argument('-z','--subtract',action='store_true', help=
      'Subtract reference image ref0')
    defaultPname = 'ebic.avt24'
    parser.add_argument('pname', nargs='*', 
      default=[defaultPname],
      help='''Image stream source. i.e: -b ado ebic.avt29
or -b file /operations/app_store/cameras/run_fy19/TEST/g2-lerec.laser-relay-cam/*.png -t100 -a100,
or -b http https://cdn.spacetelescope.org/archives/images/news/heic1523b.jpg.
for epics:
    ssh acnlinhc
    source /home/cfsd/laster/setup-epics
    imageViewer.py -b epics -w1000,800,8 13PS1:image1:ArrayData
''')
    pargs = parser.parse_args()
    #print('pargs.refRing',pargs.refRing)
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    #`````````````````````````````````````````````````````````````````````````
    if not pargs.black:
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
    if pargs.refreshRate is None:
        pargs.refreshRate = 1000. if pargs.backend == 'file' else 1.

    pname = pargs.pname

    extractor = {'raw':'raw', 'pil':'Pillow'}
    if pargs.width: # if width is provided, the pargs.extract should be set to 'raw'
        pargs.extract = 'raw'
    else:
        pargs.extract = 'pil'  
        codec = CodecPIL()

    print(parser.prog+' items:%i, first:'%len(pname)+pname[0]+\
      ' using '+extractor[pargs.extract]+', version '+__version__)
    if not pargs.hist: pargs.iso = 'Off'
    
    pixScale = 1.
    ref_diameter = 0
    ref_X, ref_Y = 0, 0
    pixelPmm = 1. #1.001 is for Addon to distinguish pixels from mm
    emit_slit_width, emit_slit_space, emit_drift_length = 0.15, 1.8, 2.0
    
    if pargs.cameraName: 
        camNm = pargs.cameraName
    else:
        if pargs.backend == 'file':
            camNm = pargs.pname[0]
            if camNm[-1] == '/':
                camNm = camNm[:-1]
            print('<'+camNm+'>')
            try: # get camNm from folder name: /xxx/.../camNm/xxx_camNm_xxx
                #camNm = pargs.pname[0].split('/')[-2]
                camNm = camNm.split('/')[-1]
                camNm = camNm.split('_')[1]
            except: 
                #if '_' in camNm:
                #    try: # assuming we are at the camera folder already
                #       camNm = camNm.split('_')[1]
                #    except: 
                pass
            print('camNm: '+camNm)
        else:
            camNm = pargs.pname[0]
    camName = 'NoName'
    xScale = 1.
    if not pargs.noDatabase:
        #````````````````````````Database stuff```````````````````````````````````
        try:
            import Sybase
            def get_tables( dbname, tablename=None, command=None ):
                '''Returns dictionary of fields and list of records from the table 
                tablename of the database dbname '''
                rc = None,None
                try:
                    dbR = Sybase.connect('OPSYB1', 'harmless', 'harmless', dbname)
                except Sybase.DatabaseError as msg:
                    printe('connecting Database: '+str(msg))
                    return rc
                cR = dbR.cursor()
                if command:
                    sql = command
                    if sql == 'getall':
                        sql = 'select * from ' + tablename
                else:
                    sql = "sp_help"
                    if tablename is not None: sql += ' ' + tablename
                try:
                    cR.execute( sql )
                except Sybase.DatabaseError as msg:
                    printe('database Error:'+str(msg))
                    return rc
                
                fields = OrderedDict()
                for i,f in enumerate(cR.description): fields[f[0]]=i
                return fields,cR.fetchall()

            # get configuration from database
            dbF,dbT = get_tables( 'appConfig', 'cameras', 'getall' )
            #print('dbF,dbT:'+str((dbF,dbT)))
            
            # substitute orientation from database
            try: name = camNm.split('.')[1]
            except: name = ''
            #print('db:'+str((dbF['camNm'],camNm,name)))
            for row in dbT:
                if row[dbF['camName']] == camNm or row[dbF['name']] == name:
                    #print('row:',row)
                    camName = [row[dbF['camName']],row[dbF['name']]]
                    if pargs.controlSystem is None:
                        pargs.controlSystem = row[dbF['system']]
                    if pargs.orientation is None:
                        # The saved images are already properly oriented
                        #if pargs.backend != 'file':
                        pargs.orientation = row[dbF['orientNum']]
                            
                    if pargs.roiRect == defaultROI:
                        pargs.roiRect = str(row[dbF['ROI_Origin_X']])
                        pargs.roiRect += ','+str(row[dbF['ROI_Origin_Y']])
                        pargs.roiRect += ','+str(row[dbF['ROI_Width']])
                        pargs.roiRect += ','+str(row[dbF['ROI_Height']])

                    try:
                        xScale = row[dbF['pixelXPmm']]
                        pixelPmm = row[dbF['pixelPmm']]
                    except:
                        printe('not in db: pixelXPmm,pixelPmm')
                        sys.exit(1)
                    if None in (xScale,pixelPmm):
                        printe('not defined in db pixelXPmm,pixelPmm')
                        sys.exit(1)

                    if pargs.sensorShape is None:
                        try :
                            pargs.sensorShape = str(row[dbF['sensorWidth']])\
                              +','+str(row[dbF['sensorHeight']])\
                              +','+str(row[dbF['calibScale']])
                        except:
                            printe('sensorShape is missing in DB for '+camname[0])
                    try: pixScale = float(pargs.sensorShape.split(',')[2])
                    except: pixScale = 1.
                    print('db sensorShape:'+str((pargs.sensorShape,pixScale)))
                    ref_diameter = row[dbF['target']]*pixelPmm
                    ref_X = row[dbF['ref_X']]
                    ref_Y = -row[dbF['ref_Y']] # the Y direction in DB is down
                    # correct for missing size factor in db
                    #corrDict = {'cs2-gun.cath-cam':0.5,
                    #            'lecs1-inj.cath-cam1':0.5,
                    #            }
                    #if camName in corrDict:
                    #    pixScale = corrDict[camName]
                    # correct the parameters 
                    ref_diameter,ref_X,ref_Y,pixelPmm =\
                      [x/pixScale for x in \
                      (ref_diameter,ref_X,ref_Y,pixelPmm)]
                    printi('cam:'+camName[0]+', orient:'+str(pargs.orientation)\
                   +',\nrefRing: d=%.2f, x=%.2f, y=%.2f'%(ref_diameter,ref_X,ref_Y))
                    emit_slit_width = row[dbF['emit_slit_width']]
                    emit_slit_space = row[dbF['emit_slit_space']]
                    emit_drift_length = row[dbF['emit_drift_length']]
        except Exception as e:
            cprintw('No database support. '+str(e))
            if pargs.backend != 'ado':
                cprintw('You may try to specify camera name explicitly using -C')
            camName = [camNm,'?']
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    if camName == 'NoName':
        printw('Camera not recognized: '+str(camNm))
        print('po',pargs.orientation)
        if pargs.orientation != None:
            sys.exit(1)
    if pargs.controlSystem is None: pargs.controlSystem = 'TEST'
    printi('Control System:'+pargs.controlSystem)
    if pargs.orientation is None: pargs.orientation = 0
    #printi('orientation: '+str(pargs.orientation))
    # convert -orientation option to -R and -F combination 
    odict = {0:(pargs.rotate,pargs.flip), 1:(0,'V'), 2:(0,'H'),
      3:(180,None), 4:(90,None), 5:(90,'V'), 6:(90,'H'), 7:(270,None)}
    pargs.rotate, pargs.flip = odict[pargs.orientation]

    # instantiate the imager
    imager = Imager(pname[0],camName)
    #if pargs.pixPerMM:
    #    pixelPmm = pargs.pixPerMM
    #print('pix/mm x:%.2f, y:%.2f'%(pixelPmm[0],pixelPmm[1]))
    imager.set_calibs(pixPerMM=pixelPmm,
      ringDiameter=ref_diameter, ringCenter=(ref_X,ref_Y), xScale=xScale)

    imager.emit_slit_width = emit_slit_width
    imager.emit_slit_space = emit_slit_space
    imager.emit_drift_length = emit_drift_length
    
    # instantiate the data monitor
    # note, only backend file accepts list in pname, all others should use pname[0]
    # TODO: arguments for PVMonitorFile could be omitted, use pargs instead
    if pargs.backend == 'file':
        pvMonitor = PVMonitorFile(pname,camName=camName[0],refreshRate=pargs.refreshRate)
    elif pargs.backend == 'http':
        import requests as Backend
        if pname[0] == defaultPname:
            #pname = 'https://upload.wikimedia.org/wikipedia/commons/thumb/5/5a/Hubble_deep_field.jpg/584px-Hubble_deep_field.jpg'
            #pname = 'https://www.ifa.hawaii.edu/~kaiser/pictures/ntt/a1689.gif'
            #pname = 'https://www.hep.shef.ac.uk/research/dm/images/hubbleDeepField.jpg'
            pname = ['http://www.dlr.de/dlr/en/Portaldata/1/Resources/bilder/portal/portal_2012_3/scaled/GalaxyClusterAbell1689_sn_l.jpg']
        pvMonitor = PVMonitorHTTP(pname[0])
    elif pargs.backend == 'ado':
        from cad import cns as Backend
        #print('cns version '+Backend.__version__)
        pvMonitor = PVMonitorAdo(pname[0],refreshRate=pargs.refreshRate)
    elif pargs.backend == 'epics':
        import epics as Backend
        pvMonitor = PVMonitorEpics(pname[0])
    else:
        print('Unknown backend:',pargs.backend)
        exit(8)
        
    #imager.start_imager()
    #time.sleep(2)
    pvMonitor.dbg = pargs.dbg
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
if __name__ == '__main__':

    # enable Ctrl-C to kill application
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)    
    
    app = QtGui.QApplication([])
    main()
    imager.start_imager()
    try:
        # Start Qt event loop unless running in interactive mode or using pyside
        #print('starting QtGUI'+('.' if pargs.file else ' Waiting for data...'))
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()
    except KeyboardInterrupt:
        print('keyboard interrupt: exiting')
    print('Application exit')
    imager.stop()
    #print('exit')
    #sys.exit(app.exec_())



