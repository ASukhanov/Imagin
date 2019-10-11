#!/usr/bin/env python
from __future__ import print_function
''' Interactive Image Viewer/Analyzer of streamed images.
 
Features:
+ Input source: 
    * file, 
    * HTTP link to image, 
    * RHIC ADO parameter (requires cns.py), 
    * EPICS PV (requires epics.py).
    * USB
+ Wide range of image formats.
+ 16+ bit/channel images supported.
+ Image orientation and rotation (program options: -o and -R).
+ Interactive zooming, panning, rotation.
+ Contrast control: Displays histogram of image data with movabled region defining the dark/light levels.
+ ROI and embedded plot for measuring image values.
+ Isocurve. The isocurve level defines the threshold for spot finding.
+ Fast multi-spot finder, reports and logs centroid position and integral of most intense spots in the ROI.
+ Gaussian (1D and 2D) fit of the brightest spot.
+ Export as PNG,TIFF, JPG..., SVG?, Matplotlib, CSV, HDF5.
+ Interactive python console with access to image data, graphics objects and shell.
+ Configuration and reporting in the parameter dock.
+ Reference Images: save/retrieve image to/from a reference slots.
+ Binary operation on current image and a reference: addition, subtraction.
+ Background subtraction using a reference image.
+ User add-ons (-u option)
+ Fast browsing/cleanup of image directories
'''
#
__version__ = 'v211 2019-10-11'# Improved profiling.
#TODO: Thread entry time is 40 ms at 90%CPU. Consider using qt timer insted of event-slot mechanism

import sys, os, subprocess, time, datetime, struct
import traceback
import threading
from collections import OrderedDict
from json import dumps
    
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

# The MKL library take lots of CPU when no analysis is carried on.
# It seen using pstack on the busy thread, which shows that the time is spent
# in libiomp5.
# Limit the number of threads to 1 will disable the MKL.
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np

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
X,Y = 0,1
#````````````````````````````Helper Functions`````````````````````````````````        
#````````````````````````````Stuff for profiling``````````````````````````````
from timeit import default_timer as timer
profilingState = {}
def profile(state):
    profilingState[state] = timer()

def profStates(first,last):
    """Return text lines with time differences between intermediate states
    Do not print entry, prefixed with '>'"""
    txt = ''
    l = timer()
    t = 0
    sortedt = sorted(profilingState.items(), key=lambda kv: kv[1])
    #print('sortedt',sortedt)
    for key,value in sortedt:
        if key == first: 
            t = value
        elif key == last:
            break
        if t:
            d = value - t
            if key[0] != '>':
                txt += 'time of '+key+' :%.4f'%(d)+'\n'
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
    """Check if path exists, if not, then create it"""
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
    if len(data.shape) < 3:
        return data
    r = data[:,:,0].astype('uint16')
    g = data[:,:,1].astype('uint16')
    b = data[:,:,2].astype('uint16')
    if graysum: # uniform sum
        #rc = (r/3.+ g/3. + b/3.).astype('uint16') #12ms
        rc = ((r + g + b)/3).astype('uint8') #3ms,  4 times faster
    else: # using perception-based weighted sum 
        rc = (0.2989 * r + 0.5870 * g + 0.1140 * b).astype('uint8')
    #rc = np.mean(data,axis=2).astype('uint8') #this is slightly slower, check on large images
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
        datao = np.rot90(data,int(n)) # fast rotate by integral of 90 degree
    else: 
        #if degree: data = st.rotate(data, degree, resize = True, preserve_range = True)
        if degree: 
            datao = ndimage.rotate(data, degree, reshape = True, order=1)
    #printi('rotate time:'+str(timer()-s))
    if pargs.flip:
        if   pargs.flip == 'V': return datao[::-1,...]
        elif pargs.flip == 'H': return datao[:,::-1,...]
    return datao

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
def stdev(x,weights): #using np.average
    '''Returns sum, mean and standard deviation of array of weights
    using np.average. Slightly slower than the dot-based but 
    mathematically stable'''
    mean = np.average(x,weights=weights)
    stDev = np.sqrt(np.average((x - mean)**2, weights=weights))
    return weights.sum(), mean, stDev 
def stdev_dot_based(x,weights):
    '''Returns mean and standard deviation of array of weights
    using np.dot. Slightly faster than the average-based.
    Math. stability is not concern on 64-bit machine'''
    s = np.dot(weights,x)
    ss = np.dot(weights,x**2)
    ws = weights.sum()
    return ws, s/ws, np.sqrt(ws*ss - s**2)/ws
  
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
        
#````````````````````````````Fitting helpers``````````````````````````````````
def calculate_fitRange(region,iPosMax,width):
    # calculate range of the fitting area
    # IMPORTANT! using less that 6*width may cause underestimated baseline
    lArr = region.shape[1],region.shape[0]
    fitRange = []
    for axis in (X,Y):   
        halfRange = int(round(width[axis]))*6
        il = max(0,iPosMax[axis] - halfRange)
        ir = min(lArr[axis],iPosMax[axis]+halfRange)-1
        #print('il,ir',il,ir,iPosMax,halfRange)
        if ir - il < 6:
            printw('Not enough area to calculate background for '+'XY'[axis]\
              +', need 6, got %d'%(ir - il))
            fitRange.append((0, lArr[axis]-1))
        else:
            fitRange.append((il,ir))
    return fitRange
        
def estimate_spot_amplitude(region,iPosMax):
    subRegion = region[iPosMax[1]-1:iPosMax[1]+2,iPosMax[0]-1:iPosMax[0]+2]
    #print('esa',subRegion.shape,subRegion.max(),subRegion.mean())
    return subRegion.mean()

def estimate_base(a):
    # calculate the base as a mean along the boundary
    h,w = a.shape
    return np.mean(list(a[0]) + list(a[:,w-1]) + list(a[h-1]) + list(a[:,0]))

from scipy.optimize import curve_fit
RankBkg = 1 # Rank = 3 is troublesome for wide peaks, it favors the quadratic background 
PBase = 0           # Function parameter defining the baseline
PPos = RankBkg      # Function parameter defining the peak position
PWidth = RankBkg+1  # Function parameter defining the peak width
PAmp = RankBkg+2    # Function parameter defining the peak amplitude
#````````````````````````````Fit X and Y projections separately```````````````
# Strange, it takes 230% CPU while 2dgauss only 50%
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
def fit_gaus1D(xr, yr, guess, bounds=(-np.inf, np.inf),axis='?'):
    try:
        fp,pcov = curve_fit(func_sum_of_peaks,xr,yr,
          guess)
          #
        stdevPars = np.sqrt(np.diag(pcov))/abs(fp)
        sigmaStdev = stdevPars[PWidth]
        if sigmaStdev > 0.2 or math.isnan(sigmaStdev):
            msg = 'fitted sigma'+axis+' is bad, stdev:%.4g'%stdevPars[PWidth]
            cprint('WARNING: '+msg)
            return []
    except Exception as e:
        msg = 'Fit'+axis+' failed: '+str(e)
        cprint('WARNING: '+msg)
        return []
    return fp
def fit_region_1D(fRegion,posGuess,width,maxLevel,base,fitBaseBounds):
    fittedPars = []
    pos = []
    sigmas = []
    #
    for axis in (X,Y):
        projection = np.sum(fRegion,axis=axis,dtype=np.float)
        nP = fRegion.shape[axis] # number of points summed in each column
        xyLetter = 'XY'[axis]
        #
        
        amp = (maxLevel - base)*fRegion.shape[axis]/2.5066
        axBase = base*fRegion.shape[axis]
        guess = [axBase,posGuess[axis],width[axis],amp] # if  RankBkg = 1
        #

        # bounds may significantly slow the performance, from 7ms to 20ms
        gl = len(guess) - 1
        bounds = ([fitBaseBounds[0]*nP]+[-np.inf]*gl,
          [fitBaseBounds[1]*nP]+[+np.inf]*gl)        
        #
        
        lArr = len(projection)
        xr = range(lArr)
        yr = projection
        
        ts = timer()
        fp = fit_gaus1D(xr,yr,guess,bounds,axis=xyLetter)
        #print( 'fitted: '+', '.join(['%.2f'%i for i in fp])+'. time:%.3f'%(timer()-ts))
        fittedPars.append(fp)
    try:
        pos = np.array((fittedPars[0][PPos],fittedPars[1][PPos]))
        sigmas = np.array((fittedPars[0][PWidth],fittedPars[1][PWidth]))
    except Exception as e:
        #printw('exception in fit_region_1D:'+str(e))
        pass
    return pos,sigmas,fittedPars

#````````````````````````````Gauss2D``````````````````````````````````````````
# it is actually less CPU consuming than 1D (50% vs 250%) for small area
def twoD_Gaussian_tilted(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x,y = xy
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def oneD_Gaussian(x, amplitude, xo, sigma):
    return amplitude*np.exp(-(0.5*((x-xo)/sigma)**2))

def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, offset):
    x,y = xy   
    xo = float(xo)
    yo = float(yo)    
    a = 1./(2*sigma_x**2)
    c = 1./(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + c*((y-yo)**2)))
    return g.ravel()

def fit_region_2D(fRegion,pos,width,amp,base,fitBaseBounds):
    #
    fittedPars = []
    h,w = fRegion.shape
    x = np.arange(w)
    y = np.arange(h)
    x, y = np.meshgrid(x, y)

    g = 1.5 # moments width is usually underestimated
    guess = [amp-base,pos[0],pos[1],width[0]*g,width[1]*g,base]
    #
    
    if guess[0] < base/2.:
        printw('too dim, guess amp:%.2f'%guess[0]+' < base/2:%.2f'%(base/2))
        return pos,width,[]

    try:
        ts = timer()
        fp,pcov = curve_fit(twoD_Gaussian, (x, y), fRegion.ravel(), guess)
        #
        stdevPars = np.sqrt(np.diag(pcov))/abs(fp)
        if stdevPars.max() > 0.2:
            i = np.argmax(stdevPars)
            pnames = ['amp','posX','posY','sX','sY','theta','base']
            msg = 'fitted parameter '+pnames[i]+' bad stdev:%.4g'\
              %stdevPars[i]
            printw(msg)
            return pos,width,[]
    except Exception as e:
        msg = '2DFit failed: '+str(e)
        printe(msg)
        #
        return pos,width,[]
    sX,sY = abs(fp[3]),abs(fp[4])
    rfp = [(fp[5],fp[1],sX,fp[0]),(fp[5],fp[2],sY,fp[0])]
    return np.array((fp[1],fp[2])), np.array((sX,sY)), rfp
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
def find_spots(region, threshold, maxSpots, fitBrightest = 'None',
               fitBaseBounds = (-np.inf,+np.inf), ofs=(0.,0.)):
    """Find up to maxSpots in the ndarray region and return its centroids 
    and sums"""
    #print( 'fitBaseBounds',fitBaseBounds)
    profile('startFind')
    
    above_threshold = np.copy(region) 
    #above_threshold[above_threshold < threshold] = 0
    
    # Subtract the threshold, otherwise, if the pedestal is high,
    # we will get the sigma = ~width/sqrt(12)
    above_threshold = region - threshold
    # zero all the pixels with negative values
    above_threshold[above_threshold < 0.] = 0
    
    profile('thresholding') # 5% of processing time spent here
    
    # now find the objects
    labeled, number_of_objects = ndimage.label(above_threshold)
    if number_of_objects == 0:
        #
        return []
    profile('sums') #TODO: 30% of processing time spent here
    peak_slices = ndimage.find_objects(labeled)
    profile('find') #TODO: 15% of processing time spent here
    
    # sort the labels according to their sums
    sums = ndimage.sum(above_threshold,labeled,
      index=range(1,number_of_objects+1))
    sumsSorted = sorted(enumerate(sums),key=lambda idx: idx[1],reverse=True)
    labelIndexesSortedBySum = [i[0] for i in sumsSorted]
    profile('sort') #TODO: 15% of processing time spent here
    
    # calculate centroids
    centroids = []
    for spotIdx,lbl in enumerate(labelIndexesSortedBySum[:maxSpots]):
        islice = peak_slices[lbl]
        ofsXY = islice[1].start, islice[0].start
        only_labeled = np.copy(above_threshold[islice])
        only_labeled[labeled[islice] != lbl+1] = 0 # zero all not belonging to label lbl+1
        p,w,pcc,pixSum = centroid(only_labeled.T)
        pos = p + ofsXY

        too_small = min(w) < 1.
        #
        if too_small:
            continue

        fittedPars = []
        fitRange = []
        if spotIdx==0 and fitBrightest != 'None':
            iPosMax = [int(round(i)) for i in pos]
            fitRange = calculate_fitRange(region,iPosMax,w)
            amp = estimate_spot_amplitude(region,iPosMax)
            fRegion = region[fitRange[1][0]:fitRange[1][1],
              fitRange[0][0]:fitRange[0][1]]
            posGuess = [pos[0] - fitRange[0][0], pos[1] - fitRange[1][0]]
            #
            base = estimate_base(fRegion)
            brightness = pixSum/w[0]/w[1]
            #
            minBrightness = 1.*base
            if brightness < minBrightness:
                #
                break

            # fit the region
            if fitBrightest == '2D':
                pos,w,fittedPars = \
                  fit_region_2D(fRegion,posGuess,w,amp,base,fitBaseBounds)
            else:
                pos,w,fittedPars = \
                  fit_region_1D(fRegion,posGuess,w,amp,base,fitBaseBounds)
            if len(pos) == 0:
                break
            pos[X] += ofs[X] + fitRange[X][0]
            pos[Y] += ofs[Y] + fitRange[Y][0]
        else:
            pos[X] += ofs[X]
            pos[Y] += ofs[Y]            

        centroids.append((pos, w, pcc, pixSum, fittedPars, fitRange))
    profile('centroids') # 15% of processing time spent here, 09/10:30%
    return centroids 
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
class MainWindow(pg.QtGui.QMainWindow):
    # closeEvent() overridden to exit application properly
    def closeEvent(self, event):
        # recursion may happen if it is called from myExit, but that is fine
        print('>MainWindow.closeEvent')
        imager.myExit()
        
#````````````````````````````Graphs outside of the main widget````````````````
class GraphRoiProj():
    def __init__(self,axis,pxPmm=10.,maxPix=640):
        self.oppositeAxis = {X:Y,Y:X}[axis]
        self.axis = axis
        self.widget = None
        self.pxPmm = pxPmm
        self.maxPix = maxPix
        self.gausIntegral = 2.5067 #(sqrt(2.*3.1415))
        
    def update(self,xyArray,ofs=(0,0),peakPars=(),fitRange=()):
        #
        l = xyArray.shape[:2][self.oppositeAxis]
        a = xyArray.mean(axis=self.axis)
        offs = ofs[self.axis]
        x = np.arange(l+1,dtype=float)
        
        xmin = imager.pix2mm(offs,offs)[self.axis]
        xmax = imager.pix2mm(offs+l,offs+l)[self.axis]
        xPoints = np.linspace(xmin,xmax,l+1)        
        if self.widget is not None:
            self.widget.plot(xPoints,a,stepMode=True,pen='k',clear=True)
        else:
            self.widget = pg.plot(xPoints,a,stepMode=True,pen='k',clear=True)
            self.widget.showGrid(True,True)
            self.hv = {X:'Horizontal',Y:'Vertical'}
            self.widget.win.setWindowTitle(self.hv[self.axis]+' ROI Means')
            self.widget.win.resize(480,320)
            self.widget.setLabel('bottom','Position (mm)')
            self.widget.setLabel('left','Average intensity')
        wsum,mean,std = stdev(x[:-1],a)
        mean += offs
        std = abs(imager.mm_per_pixel(std,std)[self.axis])
        pos = imager.pix2mm(mean,mean)[self.axis]
        topLine = 'Mean:%.2f, StDev:%.2f, Sum:%i'%(pos,std,wsum)
        try:    pp = peakPars[self.axis]
        except: pp = []
        #
        if len(pp) == 0:
            self.widget.setLabel('top',topLine)
            return
        posPix = pp[RankBkg+0] + offs
        pos = imager.pix2mm(posPix,posPix)[self.axis]
        sigPix = abs(pp[PWidth])
        sig = imager.mm_per_pixel(sigPix,sigPix)[self.axis]
        txt2d = '2D' if pargs.finalFit=='2D' else '1D'
        self.widget.setLabel('top',txt2d+'-fitted pos:%.3f mm'\
        %pos+', sigma:%.2f mm.<br>'%sig+topLine)
        
        # calculate scaling factor to adjust to projection plot
        scale = 1./xyArray.shape[self.axis]
        il,ir = fitRange[self.axis]
        try:
            if pargs.finalFit=='2D':
                sigOppAxis = abs(peakPars[self.oppositeAxis][PWidth])
                g = sigOppAxis*self.gausIntegral*scale
                a,xo,sig = pp[PAmp], pp[PPos], pp[PWidth]
                #fitY = g*oneD_Gaussian(x,a,xo,sig) + pp[PPos]
                fitY = g*(oneD_Gaussian(range(ir-il),a,xo,sig)) + pp[PBase]
            else:
                fitY = func_sum_of_peaks(range(ir-il),*(pp))*scale
            self.widget.plot(xPoints[il:ir],fitY[0:ir-il],
              pen=pg.mkPen(color=(0, 0, 255), style=QtCore.Qt.DashLine))
        except Exception as e:
            printw('exception in plotting:'+str(e))
          
    def __del__(self):
        printi('deleting GraphRoiProj:'+str(GraphRoiProj))
        self.widget.win.close()

class GraphRoiIntensity():
    def __init__(self):
        self.widget = None

    def roiIntensity(self,a):
        m = a.max()
        y,x = np.histogram(a.flatten(),bins=np.linspace(0,m+1,m+2))
        return x,y

    def update(self,array,ofs=(0,0),peakPars=(),fitRange=()):
        if self.widget is not None:
            self.widget.plot(*self.roiIntensity(array),stepMode=True,pen='b'\
            ,clear=True)
            mean = array.mean()
            std = array.std()
            self.widget.setLabel('top','Mean:%.2f'%mean+', std:%.2f'%std)
        else:
            self.widget = pg.plot(*self.roiIntensity(array),stepMode=True\
            ,pen='b')
            self.widget.win.resize(480,320)
            self.widget.win.setWindowTitle('Pixel Amplitudes in ROI')
            self.widget.showGrid(True,True)
            self.widget.setLabel('bottom','Pixel amplitude')
            self.widget.setLabel('left','Count')
            self.widget.setLogMode(y=True)

    def __del__(self):
        printi('deleting GraphRoiIntensity')
        self.widget.win.close()
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

        # call the init method of the parent class
        super(CustomViewBox, self).__init__()
        # the above is equivalent to:#pg.ViewBox.__init__(self, **kwds)

        # IMPORTANT: menu creation is deferred because it is expensive 
        # and often the user will never see the menu anyway.
        self.menu = None
           
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
            #printd('menu exist')
            return self.menu
        #
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
        axesMMGui.setChecked(imager.axesInMm)
        axesMMGui.stateChanged.connect(lambda x: imager.set_axes_scale(x))
        axesMM.setDefaultWidget(axesMMGui)
        self.menu.addAction(axesMM)

        # isocurve
        isocurve = pg.QtGui.QWidgetAction(self.menu)
        isocurveGui = pg.QtGui.QCheckBox('Isocurve')
        isocurveGui.setChecked(1)
        isocurveGui.stateChanged.connect(
          lambda x: imager.show_isocurve(x))
        isocurve.setDefaultWidget(isocurveGui)
        self.menu.addAction(isocurve)
        
        return self.menu
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
def MySlot(a):
    """Global redirector of the SignalSourceDataReady"""
    #printd('>MySlot received event:'+str(a))
    if imager is not None:
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
        #
        self.p_X.sigValueChanged.connect(self.update_refRing)
        self.p_Y.sigValueChanged.connect(self.update_refRing)
        
    def p_save_changed(self):
        # update DB 
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
        global Image
        from PIL import Image
        self.img = None

    def load(self,fname):
        # takes 40ms for 2.4Mpixel image
        self.img = Image.open(fname)
        return np.asarray(self.img)
        
    def loadFromData(self,blob): # not used
        import io
        if isinstance(blob,tuple):
            blob = bytearray(blob)
        self.img = Image.open(io.BytesIO(blob))
        return np.asarray(self.img)

    def save(self,fileName,npArray):
        if True:#try:
            #ts = timer()
            if self.img is None:
                vertFlipped = npArray[::-1,...].astype('int32')
                img = Image.fromarray(vertFlipped)
            else:
                img = self.img
            img.save(fileName,'PNG')
            #print('saving time',timer()-ts)
        else:#except Exception as e:
            printe('in CodecPIL.save():'+str(e))
codec = CodecPIL()
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       
#````````````````````````````Back-ends: image delivery interface``````````````
''' The Monitor object is instantiated as:
pvm = PVMonitorXXX(sourceName, reader = readerName)
Derived class must override at two functions: getTimeStamp(), get_data_and_timestamp()
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
        self.refreshTime = 1
            
    def set_refresh(self,r):
        rr = 0.99/r
        printi('Refresh rate changed to %0.4g'%rr)
        self.refreshTime = rr

    def monitor(self):
        """starts the monitoring"""
        printi('pvmonitor.monitor() is not instrumented') 
        
    def start_thread_process(self):
        # start the event thread 
        thread = threading.Thread(target=self.thread_process)
        thread.start()        
        # thread safe data delivery
        self.SignalSourceDataReady.connect(MySlot)

    def thread_process(self):
        # dummy thread to operate EventPocessingFinished and self.SignalSourceDataReady
        #TODO: using cprint is not thread-safe as it calls GUI.
        while True:
            EventPocessingFinished.wait(2) #
            if EventPocessingFinished.is_set():
                EventPocessingFinished.clear()
                ct = timer()
                time.sleep(self.refreshTime)
                if pargs.profile:
                    profile('>Flag cleared')
                    print('Thread entry time: %.4f'%\
                      #(profilingState['Flag cleared']-profilingState['Finish']))
                      (ct-profilingState['Finish']))
                self.SignalSourceDataReady.emit('sourceDataReady')
            else:
                #print('Processing timeout')
                if self.exit:
                    cprint('FIXME.Streaming finished. You have to restart the program')
                    return

    def clear(self):
        """clears a monitor, stops related threads"""
        printi('pvmonitor.clear() is not instrumented') 

    def get_data_and_timestamp(self):
        """returns the image ndarray and timestamp used for polling data delivery
        if stream has been finished then the ndarray is empty, timestamp = 0
        """
        printi('pvmonitor.get_data_and_timestamp() is not instrumented')
        return [],0
        
    def get_initial_setting(self,cameraName=None):
        """Usually calls once, during initialization, returns threshold level,
        ROI, fitting mode and pedestal reference name"""
        return OrderedDict((('thresholdS',None),('roiS',None),
          ('subtractPedS','None'),('fitS','None'),('deSpeckleS',0)))
#````````````````````````````PVMonitor for data from a file```````````````````
class PVMonitorFile(PVMonitor):
    def __init__(self,pvname,**kwargs):
        super(PVMonitorFile,self).__init__()
        if len(pvname)==1:
            import glob
            if os.path.isdir(pvname[0]):
                globname = pvname[0]+'/'+'*.png'
            else:
                globname = pvname[0]
            self.fileList = sorted(glob.glob(globname))
        else:
            self.fileList = pvname
        nFiles = len(self.fileList)
        print('Files in the stream: %d'%nFiles)
        if nFiles == 0:
            # select files of today from the camera
            globname = '/operations/app_store/cameras/'+pargs.year+'/'\
              +pargs.controlSystem+'/'+pvname[0]+'/'+pargs.prefix+pvname[0]\
              +time.strftime('_%Y%m%d*')
            self.fileList = sorted(glob.glob(globname))
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
        
        self.start_thread_process()

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
                try:
                    source,fn = fullname.rsplit('/',1)
                    source += '/' + fn
                except:
                    source, fn = fullname, fullname
                os.rename(source,trashDir+fn)
            cprint('Moved %d files '%nImages+' out of %d, from '%len(self.fileList)\
              +firstFile+' through '+lastFile+' to '+trashDir)
            self.fileList[fromIdx:toIdx] = []
            self.curFileIdx = fromIdx
            #print('after:',len(self.fileList),self.curFileIdx)
        except Exception as e:
            cprinte('in trash(): '+str(e))
            
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
            cprint('Processed %d images'%self.profN+' in %0.4g s'%(time.time()-self.profTime))
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
        # extract timestamp from the filename
        try:
            ds,dt = fname.split('_')[-2:]
            dsdt = (ds+dt)
            dsdt.strip('.png')
            seconds = dsdt[:14]
            timestamp = time.mktime(datetime.datetime.strptime(seconds,"%Y%m%d%H%M%S").timetuple())
            timestamp += float(dsdt[14:20])*1.e-6 # add microsecond part
        except:
            timestamp = time.time()
        #print('timestamp:',timestamp)
        
        # load image
        profile('>il')
        self.data = codec.load(fname)
        profile('image loading')
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
            print('history:',self.req.history)
            print('Content-Type:',self.req.headers['Content-Type'])
            print('cont:',type(self.req.content),len(self.req.content),self.req.content[:20])
        self.data = codec.loadFromData(self.req.content)
        return self.data,time.time()

    def clear(self):
        return

#
#````````````````````````````PVMonitor of a an EPICS Process Variable`````````
class PVMonitorEpics(PVMonitor):
    def __init__(self,pvname,**kwargs):
        super(PVMonitorEpics,self).__init__()
        self.pvsystem = 'Epics'
        #epics.ca.replace_printf_handler(self.handle_messages)
        print('epics version: '+str(Backend.__version__))
        self.pv  = Backend.PV(pvname)
        time.sleep(.1) # give it some time to initiate
        self.ts = 0
        print('Epics:',self.pv)        
        try: r = kwargs['refreshRate']
        except: r = 1.
        self.set_refresh(r)
        self.start_thread_process()
        
    def getTimeStamp(self):
        r = self.pv.get_timevars()
        if r is None:
            printe('No PV:'+str(self.pv))
            sys.exit(4)
        return r['timestamp']
        
    def get_data_and_timestamp(self):
        if pargs.dbg:
            print('>data, dt:%0.4g'%(self.pv.timestamp - self.ts))
        self.ts = self.pv.timestamp
        data = self.pv.get()
        #print( 'pv.get:',len(data),data[:100])
        self.blob = data
        return data,self.ts
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#````````````````````````````PVMonitor of an USB camera``````````````````````````
class PVMonitorUSB(PVMonitor):
    def __init__(self,pvname,**kwargs):
        super(PVMonitorUSB,self).__init__()
        self.pvsystem = 'USB'
        # if usb port is not numeric, search ports from 9 to 0
        cams = [int(pvname)] if unicode(pvname).isnumeric() else range(9,-1,0)
        print('Trying to open one of the usb cameras: [%s]'%pvname)
        #````````````````````camera initialization
        # capture from the LAST camera in the system
        # presumably, if the system has a built-in webcam it will be the first
        for i in cams:
            #printd("Testing for presence of camera %i..."%i)
            cv2_cap = Backend.VideoCapture(i)
            if cv2_cap.isOpened():
                print('USB camera %i opened'%i)
                break
        
        if not cv2_cap.isOpened():
            print("Camera not found!")
            exit(1)
        self.videoCapture = cv2_cap
        try: r = kwargs['refreshRate']
        except: r = 1.
        self.set_refresh(r)
        self.start_thread_process()
            
    def get_data_and_timestamp(self):
        ret, img = self.videoCapture.read()
        ts = time.time()
        if not ret:
            printw("Error reading image")
            return [],ts
        #
        return img, ts

    def clear(self):
        self.videoCapture.release()
        print('USB camera released')
        self.exit = True
             
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
        self.cleanImage = False # show image only
        self.zoomROI = False # zoom image to ROI
        self.maxSpots = pargs.maxSpots # number of spots to find
        self.spots = [] # calculated spot parameters
        self.spotShapes = [] # spot finder graphics objects
        self.marksColor = (0,170,0) # color of marks and spot contours
        self.mainSpotTxt = ''
        self.roiArray = np.array([]) # array of ROI-selected data
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
        #
        self.outsideGraphs = {}
        #self.roiIntensityPlot = None
        self.imageMan = None # handle to imageMan ADO if 'UpdateMan' is on
        self.is_host_priviledged = True
        self.blurWidth = 0.
        self.fitBaseBounds = [-np.inf, +np.inf]
        self.normalize = pargs.pixLimit == 0.
        self.stdPars = []
        self.addon = None
        
        self.emit_slit_width = 0.
        self.emit_slit_space = 0.
        self.emit_drift_length = 0.
        
        # connect signal to slot
        #
        #for thread safe use of cprint in threads
        self.SignalCPrint.connect(cprint)

        #
        profile('>start')

    def start_imager(self):
        print('>start imager')
        if self.cameraName is None:
            if self.backend in ['ado','epics']:
                self.cameraName = [self.pvname,'?']
        if self.cameraName is None:
            self.imagePath = '/tmp/imageViewer/'
            self.cameraName = ['noName','?']
        else:
            self.imagePath = pargs.logdir+pargs.controlSystem+'/'+self.cameraName[0]+'/'
        if self.backend in ('file','http'):
            self.imageFile = self.pvname
            self.refPath = self.imagePath+'refs/'
            self.logPath = self.imagePath+'logs/'
        else:
            checkPath(self.imagePath)
            self.refPath = self.imagePath+'refs/'
            checkPath(self.refPath)
            self.logPath = self.imagePath+'logs/'
            checkPath(self.logPath)
        print('imagePath:'+str(self.imagePath))
        
        # try to obtain threshold,roiRect,subtractPeds from the manager
        threshold,roiRect,subtractPeds,fitS, despeckleKernel=\
          pvMonitor.get_initial_setting(self.cameraName[0]).values()
        #print( 'threshold,roiRect,subtractPeds,fitS, ixLimitS',threshold,roiRect,subtractPeds,fitS,deSpeckleS)
        if subtractPeds == 0: subtractPeds = 'None' # to deal with old imageMans
                
        # evaluate ROI Rectangle
        self.roiRect = None
        try:
            rr = pargs.roiRect
            if rr is None: 
                rr = pargs.roiRectDb
                self.roiRect = roiRect
            if self.roiRect is None:
                self.roiRect = [float(i) for i in rr.split(',')]
        except Exception as e:
            #printe('wrong -O format: '+str(pargs.roiRect)+', '+str(e))
            printw('ROI is not recognized, use default')
            self.roiRect = 0.05,0.05,0.9,0.9# default
        self.initialRoiRect = self.roiRect
        
        # evaluate threshold
        self.threshold = threshold if pargs.threshold is None\
          else pargs.threshold
        if self.threshold is None: self.threshold = 0

        # evaluate pedestal reference
        self.subtractPeds = subtractPeds if pargs.subtract == 'None'\
          else pargs.subtract
        #print('subtractPeds',self.subtractPeds)

        # evaluate fitting mode
        if pargs.finalFit == 'None':
            # take care of old managers
            try:    fitS = {'FitLess':'None', 'FinalFit':'1D'}[fitS]
            except: pass
            pargs.finalFit = fitS
        printi('Fitting mode:'+str(pargs.finalFit))
        
        ## evaluate pixLimit
        self.pixLimit = pargs.pixLimit
        #if self.pixLimit == 0:
        #    self.pixLimit = pixLimitS
            
        # evaluate despeckling
        self.despeckleKernel = pargs.despeckle
        if self.despeckleKernel == 0:
            self.despeckleKernel = despeckleKernel

        # initate addon
        if pargs.userAddon:
            addon = 'ivAddon_'+pargs.userAddon
            Addon = importlib.import_module(addon)
            try:
                self.addon = Addon.Addon(self)
            except Exception as e:
                printw('could not start the Addon:'+str(e))
                self.addon = None
        else:
            self.addon = None        

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
        self.cb_update_ROI()
        self.started = True

    def myExit(self):
        print('>Exit')
        if self.addon: 
            self.addon.exit()
        for item in self.outsideGraphs:
            #self.outsideGraphs[item].__del__()
            self.outsideGraphs[item].widget.win.close()
        self.win.close()
        #self.exit() # has no effect

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
                       ('subtractPedS',self.subtractPeds),
                       # roiS can be updated only in 'expert' mode
                       ('roiS',self.roiRect) if pargs.expert else ('',None),
                       ('deSpeckleS',self.despeckleKernel),
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
            #if parval[0] == 'fitS':
            #    val = 'FinalFit' if val else 'FitLess'
            if old == val:
                #printi('old value of '+par+' is the same, no need to change')
                continue
            txt = 'Are you going to modify '+par+' on imageMan for '+\
              imageManName+' from '+str(old)+' to '+str(val)
            ok = qMessage(txt)
            if not ok:
                continue
            status = pvMonitor.httpDev.set(self.imageMan,par,val) # http supports set history
            #print( 'setting '+self.imageMan+':'+par+' to '+str(val)+': '+str(status))
            d[par]['new'] = val
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
        # Update sizeFactorM only once at first image
        w = float(max(self.data.shape))# - 1 #-1 is a kludge
        sensorShape = self.hwpb[:2]
        if pargs.sensorShape is not None:
            try:
                ss = pargs.sensorShape.split(',')
                sensorShape[1],sensorShape[0] = [int(z) for z in ss[:2]]
            except: printe('sensorShape not avaialble: '+str(pargs.sensorShape))
        self.sizeFactor = float(min(sensorShape))/min(self.hwpb[:2])
        # ratio of heights is correct even when the rows were padded (by a PNG formatter)
        print('sizeFactor:%.2f'%self.sizeFactor+', sensorShape,imageShape: '\
          +str((sensorShape,self.hwpb)))

        #self.win = QtGui.QMainWindow()
        # MainWindow is QtGui.QMainWindow() but with proper handling of exit.
        self.win = MainWindow()

        self.area = pg.dockarea.DockArea()
        self.win.setCentralWidget(self.area)
        
        # image info
        self.winTitle = 'image:'+self.pvname+' hwpb:'+str(self.hwpb)
        self.win.setWindowTitle(self.winTitle)

        #````````````````````````````parameters data``````````````````````````
        self.partitioning = None if pargs.part == 'None' else pargs.part
        self.numRefs = 5
        refs = ['None']+['_ref'+str(i) for i in range(self.numRefs)]
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
                {'name':'Images/s','type':'float','value':0.,
                  'tip':'Performance','readonly':True},
                {'name':'Pause', 'type': 'bool', 'value': self.paused,
                  'tip': 'Enable/disable receiving of images, '},
                {'name':'Fast', 'type': 'bool', 'value': 0,
                  'visible':True if self.backend == 'file' else False,
                  'tip': 'Fast browsing, many features disabled'},
                {'name':'Next', 'type': 'button',
                  'tip':'Process next image from the stream'},
                {'name':'Prev', 'type': 'button',
                  'visible':True if self.backend == 'file' else False,
                  'tip':'Process previous image from the stream'},
                {'name':'FirstFile%', 'type':'slider', 'value':0,
                  'visible':True if self.backend == 'file' else False,
                  'tip':'Starting file relative position in directory.'},
                {'name':'LastFile%', 'type':'slider', 'value':100,
                  'visible':True if self.backend == 'file' else False,
                  'tip':'Ending file relative position in directory.'},
                {'name':'Saving', 'type': 'bool', 'value': self.saving,
                  'visible':False if self.backend == 'file' else True,
                  'tip': 'Enable/disable saving of images'},
                {'name':'Delete first to last', 'type': 'button',
                  'visible':True if self.backend == 'file' else False,
                  'tip':'Move selected images to trash'},                
                {'name':'View saved', 'type': 'button',
                  'visible':False if self.backend == 'file' else True,
                  'tip':'View saved images from this camera using separate application'},
                {'name':'Threshold', 'type': 'float', 'value':self.threshold,
                  'tip': 'Threshold level for spot finding, changed with isoCurve level'},
                #
                #
                {'name':'Zoom to ROI', 'type': 'bool', 'value':self.zoomROI,
                  'tip': 'Zoom to ROI'},
                {'name':'Normalize', 'type': 'bool', 'value':self.normalize,
                  'tip': 'Normalize intensity'},
                {'name':'PixLimit', 'type': 'int', 'value':self.pixLimit,
                  'tip': 'Pixel amplitude limit, if not zero, then ROI pixels will be limited to that value'},
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
                {'name':'FinalFit', 'type': 'list', 'value': pargs.finalFit,
                  'values':['None','1D','2D'],
                  'tip':'Fitting of found spots: 1D- one-dimensional on projections, 2D- two-dimensional'},
                {'name':'View results', 'type': 'bool', 'value': False,
                  'tip':'Log the spots parameters to a file'},
                {'name':'Despeckle', 'type':'int','value':self.despeckleKernel,
                  'tip':'Kernel size of a median filter for removing speckles in the ROI'},
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
                {'name':'Partitioning', 'type': 'list',
                    'value':pargs.part,     
                    'values':['None','Vertical'],
                  'tip': 'ROI prtitioning'},
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
                {'name':'FitBaseLo', 'type':'float','value':-1.e9,
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
                            self.show_isocurve(False)
                            if self.backend == 'file':
                                pvMonitor.curFileIdx = self.events
                            self.profTime = 0                            
                            EventPocessingFinished.set()
                    elif parItem == 'Next':
                        self.show_isocurve(False)
                        EventPocessingFinished.set()
                                                
                    # items for backed != 'file'
                    elif parItem == 'Saving':
                        self.saving = itemData
                        if self.saving:
                            self.save_image()
                            if self.addon: self.addon.start_saving()
                            cprint('Saving images to '+self.imagePath)
                        else:
                            if self.addon: self.addon.stop_saving()
                            cprint('Stopped saving to '+self.imagePath)
                    elif parItem == 'View saved':
                        if not os.path.exists(self.imagePath):
                            cprinte('opening path: '+self.imagePath)
                            return
                        try:
                            #cmd = 'xterm -e imageViewer.py -p -m1 -a100 -bfile '+self.imagePath
                            cmd = 'imageViewer.py -p -m1 -a100 -bfile '+self.imagePath
                            cprint('executing:'+str(cmd))
                            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True) #stderr=subprocess.PIPE)
                        except Exception as e: cprinte('in View saved'+str(e))
                    elif parItem == 'Update Mgr':
                        if self.is_host_priviledged:
                            self.update_imageMan(newState=itemData)
                        elif itemData:
                            qMessage('ERROR\nimageMan can be updated only from a MCR host')
                            # self.set_dockPar('Control','Update Mgr',False)#DNW
                    elif parItem in ('LogView','Gpm'):
                        cmd = None
                        try:
                            if pargs.controlSystem == 'LEREC':
                                cmd = [parItem,'-file','/operations/app_store/Gpm/LEReC/Cameras/img.'+self.cameraName[0]+'.logreq']
                            elif pargs.controlSystem == 'CEC':
                                cmd = [parItem,'-file','/operations/app_store/Gpm/RHIC/Systems/CeC/Cameras/img.'+self.cameraName[0]+'.logreq']
                            elif pargs.controlSystem == 'RHIC':
                                cmd = [parItem,'-file','/operations/app_store/Gpm/RHIC/Instrumentation/Cameras/img.'+self.cameraName[0]+'.logreq']
                            else:
                                cprinte('Control System '+str(pargs.controlSystem)+' not supported')
                            if cmd is not None:
                                cprint('Executing:'+str(cmd))
                                p = subprocess.Popen(cmd)
                        except Exception as e: cprinte('in LogView: '+str(e))
                    elif parItem == 'Threshold':
                        self.threshold = itemData
                        if self.imageMan is not None: 
                            self.update_imageMan()
                        if self.isoInRoi:
                            #TODO: need to relocate the isocurve to ROI origin
                            self.iso.setData(blur(self.roiArray))
                        else:
                            self.iso.setData(blur(self.grayData))
                        self.show_isocurve(True)
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
                        self.normalize = itemData
                        if itemData:
                            self.pixLimit = 0
                            self.set_dockPar('Control','PixLimit',self.pixLimit,
                              deferred = True)
                        print('update in normalize')
                        self.update_imageItem_and_ROI()
                    elif parItem == 'PixLimit':
                        self.set_pixLimit(itemData)
                    # items for backend == 'file'
                    elif  parItem == 'Prev':
                        self.events -= 2
                        pvMonitor.curFileIdx = self.events
                        EventPocessingFinished.set()                        
                    elif  parItem == 'Delete first to last':
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
                        rr = self.initialRoiRect
                        self.roi.setPos(rr[0]*w,rr[1]*h)
                        self.roi.setSize(rr[2]*w,rr[3]*h)
                    elif parItem == 'Clean Image':
                        self.cleanImage = itemData
                        if self.cleanImage:
                            self.remove_spotShapes()
                            self.show_isocurve(False)
                            self.mainWidget.removeItem(self.roi)
                            self.mainWidget.removeItem(self.pointerText)
                            self.pointerText = None
                            self.mainWidget.removeItem(self.mainSpotText)
                            self.mainSpotText = None
                        else:
                            self.mainWidget.addItem(self.roi)
                            self.mainSpotText = mainSpotText()
                            self.update_mainSpotText()
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
                        print('update in rotate()')
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
                        self.update_ROI()
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

                    elif parItem == 'Despeckle':
                        self.despeckleKernel = itemData
                        self.data = self.dataOriginal
                        self.update_imageItem_and_ROI()

                    elif parItem == 'Subtract':
                        self.pedestals = self.enable_subtraction(itemData)
                        if self.pedestals is not None:
                            self.subtract_pedestals()
                            print('update in subtract()')
                            self.update_imageItem_and_ROI()

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
                            self.outsideGraphs[parItem] = GraphRoiIntensity()
                        else:
                            try:    del self.outsideGraphs[parItem]
                            except: pass
                       
                if parGroupName == 'SpotFinder':               
                    if parItem == 'MaxSpots':
                        self.maxSpots = itemData
                    elif parItem == 'Partitioning':
                        self.partitioning = None if itemData == 'None'\
                          else itemData
                        if self.partitioning != 'Vertical':
                            try:
                                for pl in self.partLines:
                                    self.mainWidget.removeItem(pl)
                            except:
                                pass
                            return
                    self.cb_update_ROI()
                    
                if parGroupName == 'Ref Images':
                    if itemData == 'None':
                        return
                    prefix = self.refPath
                    child = self.pgPar.child(parGroupName).child(parItem)
                    fileName = prefix+itemData+'.png'
                    if parItem == 'View':
                        #cmd = ["gm",'display',prefix+'*']
                        if not os.path.exists(fileName):
                           cprint('file does not exist: '+fileName)
                           return
                        cmd = 'imageViewer.py -p -o0 -C '+self.cameraName[0]+\
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
                            self.update_mainSpotText()
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
                        self.data = self.dataOriginal
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
                        self.imageItem.setLookupTable(self.lut_greyClip)
                if parGroupName == 'Help':
                    import webbrowser
                    webbrowser.open(pargs.userGuide)
                if parGroupName == 'Exit':
                    self.myExit()
                         
        #QtGui.QSpinBox.setKeyboardTracking(False) # does not work
        self.pgPar.sigTreeStateChanged.connect(handle_change)
           
        #        

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
        self.update_mainSpotText()
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
            #

            self.contrast.setImageItem(self.imageItem)
            #
            self.gl.addItem(self.contrast)
            
            # Connect callback to signal
            #self.contrast.sigLevelsChanged.connect(self.cb_contrast_levels_changed)

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
            #print('isoLine created')
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
        cmap_greyClip = pg.ColorMap(pos=[0.,0.995,1.],color=
          [[0,0,0,255],[255,255,255,255],[255,0,0,255]])
        self.lut_greyClip = cmap_greyClip.getLookupTable()
        
        self.hideDocks(pargs.miniPanes)
        #print('<iv_show')
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

    def set_pixLimit(self,itemData):
        # drawback: due to Normalize the analysis will be done twice
        self.set_dockPar('Control','Normalize',itemData==0.)
        if itemData == 0.:
            return
        self.pixLimit = itemData
        self.contrast.setLevels(0,itemData)
        self.contrast.setHistogramRange(0,itemData)
        print('update in pixLimit()')
        # we have to update because self.pixLimit have changed
        self.update_imageItem_and_ROI()
        if self.imageMan is not None: 
            self.update_imageMan()
    
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
        
    def set_dockPar(self,child,grandchild,value,deferred = False):
        """Change dock parameter, if it is called inside the parameter 
        processing then deferred should be set True, in that case the execution
        will be done in a separate thred"""
        if deferred:
            thread = threading.Thread(target=self.set_dockPar_thread,
              args = (child,grandchild,value))
            thread.start()
        else:  
            self.pgPar.child(child).child(grandchild).setValue(value)
        
    def set_dockPar_thread(self,child,grandchild,value):
        #print('>set_dockPar_thread',child,grandchild,value)
        time.sleep(.1)
        self.pgPar.child(child).child(grandchild).setValue(value)

    def closeEvent(self, event):
        print("Closing")

    def setup_main_widget(self):
        #
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

    def set_calibs(self, pixPerMM=10., ringDiameter=100.,
                    ringCenter=(0,0), xScale=1.):
        self.ref_diameter = ringDiameter
        self.pixelPmm = [round(pixPerMM*xScale,3), round(pixPerMM,3)]
        self.ref_X,self.ref_Y = ringCenter
        self.pixelXScale = xScale

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
                strtime = datetime.datetime.fromtimestamp(self.timestamp)\
                .strftime(fmt)
            else:
                strtime = time.strftime(fmt)
            folderName = strtime.split('_')[1]
            path = self.imagePath+'/'+folderName
            checkPath(path)
            fileName = path+'/IV_'+self.cameraName[0]+strtime
            self.imageFile = fileName
        try:
            ts = timer()
            if pargs.fullsize:
                codec.save(fileName,self.data)
            else: # PNG image format
                #printd('writing original data')
                with open(fileName,'wb') as f:
                    f.write(pvMonitor.blob)
            if self.backend == 'ado':
                try:
                    r = pvMonitor.cadDev.set('img.'+self.cameraName[0],
                      'lastSavedM',fileName)
                except Exception as e:
                    printe('exception:'+str(e))
                    r = 1
                if (r):
                    printe('updating lastSavedM')
            #printd('save_image time %.3f'%(timer()-ts)+'s for %.3fkB'%(len(self.data)/1000.)) 
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
            #print('refop',data.shape,self.data.shape)
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
            print('update in refs')
            self.update_imageItem_and_ROI()
        except Exception as e:
            msg = 'in ref_op '+operation+': '+str(e)
            cprinte(msg)

    def enable_subtraction(self,ref):
        if ref == 'None' or ref is None:
            return None
        fn = self.refPath+str(ref)+'.png'
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
            cprinte('cannot load '+str(ref)+':'+str(e))
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
        x = xPx*g[0]
        y = yPx*g[1]
        return x,y
        
    def pix2mm(self,xPx,yPx):
        # pixel to mm conversion, zero mm is at the center of image
        h,w = self.data.shape[:2]
        return self.mm_per_pixel(xPx - w/2.,yPx - h/2.)

    def find_partitioned_spots(self):
        """Find one spot per partition"""
        h,w = self.roiArray.shape
        nS = self.maxSpots
        if self.partitioning != 'Vertical':
            return
        stripHeight = int(h/nS)
        croppedHeight = nS*stripHeight
        #
        strips = self.roiArray[:croppedHeight].reshape((nS,stripHeight,w))
        centroids = []
        for stripIdx,strip in enumerate(strips):
            foundSpots = find_spots(strip,self.threshold,1,
              fitBrightest = pargs.finalFit, fitBaseBounds=self.fitBaseBounds,
              ofs=(0,(stripIdx)*stripHeight))
            #print('foundSpots:',foundSpots)
            if len(foundSpots) == 0:
                centroids.append([])
            else:
                centroids.append(foundSpots[0])
        return centroids
    def standard_analysis(self):
            if self.events == 0:
                # do not analyze first image as it will be processed during update
                return
            #
            ts = timer()
            self.gl.setBackground(self.viewBoxBackgroundColor)
            ox,oy = self.roiOrigin
            self.stdPars = [] #fix:v199
            if self.partitioning is None:
                self.spots = find_spots(self.roiArray,self.threshold,
                  self.maxSpots, fitBrightest = pargs.finalFit,
                  fitBaseBounds=self.fitBaseBounds)
            else:
                self.spots = self.find_partitioned_spots()
            self.mainSpotTxt = ''
            if len(self.spots) == 0:
                self.update_mainSpotText()
                return
            #print( 'spots',self.spots)
            profile('spotsFound')
            if self.sizeFactor == 0.:
                # should not be here
                printw('TODO logic error sizeFactor=0')
                return
            spen = pg.mkPen(self.marksColor)
            if self.spotLog: 
                logdata = OrderedDict([('time:',time.strftime('%y-%m-%d %H:%M:%S'))])
            for spotIdx,spot in enumerate(self.spots):
                if len(spot) == 0:
                    continue
                p,wPx,pcc,psum,fittedPars = spot[:5]
                #print( 'fp',spotIdx,len(fittedPars))
                posPx = p + (ox,oy)
                #print( 'posPx',posPx)
                try:
                    theta,sigmaU,sigmaV = ellipse_prob05(wPx[0],wPx[1],pcc)
                except:
                    printw('in ellipse_prob05 for spot:'+str(spot))
                    self.update_mainSpotText()
                    continue
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
                            txt2d = '2D' if pargs.finalFit=='2D' else ''
                            self.mainSpotTxt +=\
                              st+txt2d+'XY:(%.2f,%.2f)'%(abs(w[0]),abs(w[1]))
                        self.update_mainSpotText()
                        
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
        if self.imageMan is not None:
            self.update_imageMan()
        if self.zoomROI:
            self.setup_zoomROI()
        self.data = self.dataOriginal
        self.update_imageItem_and_ROI(report=True)
        if self.partitioning != 'Vertical':
            return
        orig = self.roi.pos()
        size = self.roi.size()
        ds = size[1] / self.maxSpots
        try:
            for pl in self.partLines:
                self.mainWidget.removeItem(pl)
        except:
            pass
        self.partLines = []
        for i in range(self.maxSpots-1):
            y = orig[1]+ds*(i+1)
            pl = self.mainWidget.addLine(y=y)
            self.partLines.append(pl)
            
    def cb_contrast_levels_changed(self):
        #print('Levels changed',self.contrast.getLevels())
        #self.contrastLevels = self.contrast.getLevels())
        return

    def update_ROI(self,report=False, analysis=True):
        """handle changed ROI"""
        if self.roi is None:
            return

        if report:
            relativeRect = [float(v1)/float(v2) 
              for v1,v2 in zip(self.roiOrigin,self.data.shape[::-1])]
            relativeRect += [float(v1)/float(v2) 
              for v1,v2 in zip(self.roiArray.shape,self.data.shape)][::-1]
            cprint('ROI changed: '+str((self.roiOrigin,self.roiArray.shape,
              ['%0.2f'%i for i in relativeRect])))
            self.roiRect = map(lambda z: round(z,3), relativeRect)

        # the following is much faster than getArrayRegion
        slices = self.roi.getArraySlice(self.grayData,self.imageItem)[0][:2]
        self.roiOrigin = slices[1].start, slices[0].start
        self.roiArray = self.grayData[slices]
        self.colorRoiArray = self.data[slices]
        
        try:
            if self.pixLimit:
                tmp = self.roiArray.copy()
                tmp[tmp > self.pixLimit] = self.pixLimit
                self.roiArray = tmp
        except Exception as e:
            self.pixLimit = 0
        if self.blurWidth > 0. or self.despeckleKernel > 0.:
            #ts = timer()
            #print('roi max pre:',self.roiArray.max())
            # blurring will modify data, they cannot be read_only
            self.data = self.data.copy()
            self.grayData = rgb2gray(self.data)
            if self.blurWidth > 0.:
                self.roiArray = blur(self.roiArray,self.blurWidth).astype('i4')
            if self.despeckleKernel > 0.:
                self.roiArray = ndimage.median_filter(self.roiArray\
                  , size=self.despeckleKernel)
            try:
                self.grayData[slices] = self.roiArray
            except:
                printw('Exception line 2819, grayData is read_only')
            #print('blurring/despeckling time: %.4fs'%(timer()-ts))                
            #self.data[slices] = self.roiArray
            #print('roi max post:',self.roiArray.max())

        self.remove_spotShapes()
        profile('roiArray')
        
        # find spots using isoLevel as a threshold
        if self.threshold>1 and self.maxSpots>0:
            self.standard_analysis()
            profile('standard processing')
            if self.addon:
                # call addon even no spots are found, it may need this information
                self.addon.process()
                profile('addon processing')
        
        if not self.roiPlotEnabled:
            return
        # plot the ROI histograms
        yV = self.roiArray.mean(axis=X) # vertical means
        x = np.arange(len(yV)); s = False
        pen = 'k' if not pargs.black else 'w'
        if len(self.data.shape) == 2: # gray image
            self.plot.plot(x,yV,clear=True,pen=pen,stepMode=s)
        else: # color image
            # plot color intensities
            cyV = self.colorRoiArray.mean(axis=X)
            self.plot.plot(x,cyV[:,0],pen='r', clear=True,stepMode=s)#plot red
            self.plot.plot(x,cyV[:,1],pen='g',stepMode=s) # plot green
            self.plot.plot(x,cyV[:,2],pen='b',stepMode=s) # plot blue
            self.plot.plot(x,yV,pen=pen,stepMode=s)
        fitPars = []
        fitRange = []
        if len(self.spots) > 0:
          if len(self.spots[0]) > 0:
            if len(self.spots[0][4]):
                scale = 1./self.roiArray.shape[0]
                try: 
                    fitPars = self.spots[0][4][X], self.spots[0][4][Y]
                    fitRange = self.spots[0][5][X], self.spots[0][5][Y]
                    il,ir = fitRange[0]
                    #print( 'update ROI fpars',fitPars)
                    pp = fitPars[X]
                    a,x0,w = pp[PAmp],pp[PPos],pp[PWidth]
                    if pargs.finalFit=='2D':
                        sigOppAxis = abs(fitPars[Y][PWidth])
                        g = sigOppAxis*2.5066*scale
                        fitY = g*oneD_Gaussian(range(ir-il),a,x0,w)+pp[PBase]
                    else:
                        fitY = func_sum_of_peaks(range(ir-il),*pp)*scale
                    #print('il,ir',il,ir,len(x[il:ir]),len(fitY[0:ir-il]))
                    self.plot.plot(x[il:ir],fitY[0:ir-il],
                      pen=pg.mkPen(color=(0, 0, 255), 
                      style=QtCore.Qt.DashLine),stepMode=s)
                except Exception as e:
                    printe('in ROI plot: '+str(e))
                    pass
        for name,graph in self.outsideGraphs.items():
            try:
                graph.update(self.roiArray,self.roiOrigin,fitPars,fitRange)
            except Exception as e:
                print('exception in graph.update:'+str(e))
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
        #print('update_isocurve')
        self.show_isocurve(True)
        #profile('init iso')
        v = self.isoLine.value()
        # inform imager on changed threshold
        self.set_dockPar('Control','Threshold',v)
         
    def show_isocurve(self,show=True): 
        #print('show_isocurve',show)
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
        if len(self.data.shape) == 2:
            #self.imageItem.setLookupTable(self.lut_greyClip)
            pass
        else:
            msg = 'Avoiding colormapping for color images'
            printw(msg)

        # redrawing axes here have no effect as it will be overridden
        # we need to postpone the activation somehow.
        # simple delayed thread is not working.
        if self.axesInMm:
            self.needToRedrawAxes = True
        #print('need to redraw axes: ',self.needToRedrawAxes)
        if self.mainSpotText: self.update_mainSpotText()
    
    def update_imageItem_and_ROI(self,report=False):
        #DNW#self.imageItem.setImage(self.data,autoLevels=self.normalize)
        profile('>norm')
        if self.normalize:
            if self.contrast:
                # fast, 5ms/Mbyte:
                if len(self.roiArray) > 10: 
                    profile('>convolve')
                    mx = np.convolve(self.roiArray.flatten(), np.ones((4,))/4, mode='same').max()
                    profile('convolve')
                else:
                    try:
                        mx = self.roiArray.max()
                    except:
                        mx = self.data.max()
                self.contrast.setLevels(0,mx)
                self.contrast.setHistogramRange(0, mx)
        elif self.pixLimit and self.contrast is not None:
                self.contrast.setLevels(0,self.pixLimit)
                self.contrast.setHistogramRange(0, self.pixLimit)
        profile('>update_ROI')
        if self.roi:
            self.update_ROI(report)
        self.imageItem.setImage(self.data)
        if self.contrast is not None: self.contrast.regionChanged() # update contrast histogram
        profile('setImage')

    def subtract_pedestals(self):
        try:
            self.data = self.data.astype('i4') - self.pedestals
            self.grayData = rgb2gray(self.data)
        except:
            cprinte('Reference is not compatible, subtraction disabled')
            self.pedestals = None
            self.set_dockPar('Analysis','Subtract','None')
    def process_image(self):
        global profilingState
        if not pvMonitor:
            printe('pvMonitor did not start, processing stopped')
            sys.exit(10)        
        profile('start process')
        
        data,self.timestamp = pvMonitor.get_data_and_timestamp()
        #
        profile('get_data')
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
        #
        #print('data from monitor:'+str(data.dtype)+':\n'+str(data[:100]))
        if pargs.width:
            #``````````The source is vector parameter with user-supplied shape
            if self.dataLen != dataLen: # data size changed
                print('size changed',self.dataLen,dataLen)
                try:
                    pargs.width = pvMonitor.get_image_shape()
                except:
                    printw('get_image_shape() not awailable')
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
                    print('roi',self.roi)
                    self.roi.setPos(self.roiRect[0]*w,self.roiRect[1]*h)
                    self.roi.setSize(self.roiRect[2]*w,self.roiRect[3]*h)

            bytesPerChannel = ((self.hwpb[3]-1)/8)+1
            #
            if bytesPerChannel > 1: # we need to merge pairs of bytes to integers
                # 
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
        #
        profile('>gray conv')
        dataRotated = rotate(self.receivedData,self.dockParRotate)
        #

        if pargs.gray:
            self.data = rgb2gray(dataRotated)
            self.grayData = self.data
        else:
            self.data = dataRotated
            self.grayData = rgb2gray(self.data)
        self.dataOriginal = self.data.copy()
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
            if self.averageWindow:
                self.start_averaging()
            self.imageItem = pg.ImageItem(self.data)
            self.iv_show()
        #````````````````````udate data```````````````````````````````````````
        if True:
            #TODO: react on shape change correctly, cannot rely on self.hwpb because of possible rotation
            if self.saving:
                self.save_image() #TODO should not it be after update
            if self.backend == 'file':
                try:    self.winTitle = pvMonitor.pvname.rsplit("/",1)[1]
                except: self.winTitle = pvMonitor.pvname
            else:
                self.winTitle = self.pvname+time.strftime('_%H%M%S')
            self.win.setWindowTitle(self.winTitle)

            if self.pedestals is not None:
                self.subtract_pedestals()

            if self.averageWindow:
                l = len(self.averageQueue)
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
        if self.isocurve_enabled: self.show_isocurve(False)
        if self.events % 100 == 0:
            try:
                dt = timer() - profilingState['>100 events']
                self.set_dockPar('Control','Images/s',100./dt)
            except: pass
            profile('>100 events')
        self.events += 1
        self.set_dockPar('Control','Image#',str(self.events)) # concern: time=0.5ms
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,            
        if self.needToRedrawAxes:
            self.redraw_axes()
        if pargs.profile:
            ct = timer()
            print('### Total time: %.4f @ %.4f'%\
                (ct-profilingState['start process'],ct)\
                +' t/ev=%.4f'%((ct-profilingState['>start'])/self.events))
            profile('Finish')
            print(profStates('>start','Finish'))
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
        absPos = evt[0]  # using signal proxy turns original arguments into a tuple
        #
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
        #print('vrange',vrange,self.data.shape)
        right = min(self.data.shape[1],vrange[0][1])
        top = min(self.data.shape[0],vrange[1][1])
        return right,top
        
    def update_mainSpotText(self):
        self.mainSpotText.setText(self.mainSpotTxt)
        self.mainSpotText.setPos(*self.top_right_corner())

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
        print('addon.controlPanel()')
       
    def addon_clicked(self,parItem,itemData):
        """Called when user item in control pane is clicked. For example:
        if parItem == 'Where My Results?':
            w = IV.QtGui.QWidget()
            msg = 'The results are in /tmp/'
            IV.QtGui.QMessageBox.information(w,'Message',msg)
        """
        print('addon.addon_clicked',parItem,itemData)
        
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
    year = time.strftime('run_fy%y')
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
    parser.add_argument('-e','--despeckle', type=int, default=0, help=
      'median filtering of the ROI to remove speckles, -e3 will do the best cleaning')
    parser.add_argument('-F','--flip', help=
      "Flip image, 'V':vertically, 'H':horizontally")
    parser.add_argument('-f','--fullsize', action='store_true', help=
      'Use full-size full-speed imageM parameter')
    parser.add_argument('-g','--gray', action='store_false', help=
      'Set it for non-gray images')
    parser.add_argument('-H','--hist', action='store_false', help=
      'Disable histogram with contrast and isocurve contol')
    parser.add_argument('-i','--iso',default='Image',help=
      '''Isocurve drawing options: ROI - only in ROI (default), 
      Image - in full image, 
      Off - no isocurve''')
    parser.add_argument('-j','--pixLimit', type=int, default=0,help=
      'Contrast control level, 0 for auto.')    
    #parser.add_argument('-l','--logdir', default = '/operations/app_store/imageViewer/',
    logdir = '/operations/app_store/cameras/'+year+'/'
    parser.add_argument('-l','--logdir', default = logdir,
      help='When DB is not available, this defines a directory for saving images, logging and references, i.e: '+logdir)
    parser.add_argument('-m','--maxSpots',type=int,default=4,
      help='Maximum number of spots to find')
    parser.add_argument('-M','--miniPanes', action='store_true', help=
      'Start with minimized panes (i.e. control, histogram and console panes')
    #parser.add_argument('-n','--normalize', action='store_false', help=
    #  'Disable intensity normalization')
    parser.add_argument('-O','--roiRect',
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
    parser.add_argument('-t','--threshold',type=float,
      help='Threshold for spot finding')
    parser.add_argument('-T','--finalFit',default='None',
      help='Fitting of the brightest spots in regions, "None", "1D" or "2D"')
    parser.add_argument('-v','--vertSize',type=float,default=800,
      help='Vertical size of the display window')
    parser.add_argument('-V','--part',default='None', help=
      'ROI Partitioning, Vertical/Horizontal') 
    parser.add_argument('-u','--userAddon', help=
      '''User addon, the addon script name should be prefixed with ivAddon_,
      i.e -uNewAddon will try to import ivAddon_NewAddon.py''')
    ug = 'https://github.com/ASukhanov/Imagin'#
    
    parser.add_argument('-U','--userGuide',default=ug, help=
      'Location of the user instructions')
    parser.add_argument('-w','--width', help=
      '''For blob data: width,height,bits/pixel i.e 1620,1220,12. 
      The bits/pixel may be omitted for standard images''')
    parser.add_argument('-X','--pixPerMM',type=float,
      help='Force pixel/mm conversion with applied value.')
    parser.add_argument('-x','--expert',action='store_true', help=
      'Expert/Admin mode')
    parser.add_argument('-y','--year', default=year, help=
      'Year, applied for backend=file')
    parser.add_argument('-z','--subtract', nargs='?', const='_ref0', 
      default='None', help=
      'Name of the pedestal file to subtract, eg: -z_ref0')
    defaultPname = 'ebic.avt24'
    parser.add_argument('pname', nargs='*', 
      default=[defaultPname],
      help='''Image stream source. i.e: -b ado ebic.avt29
or -b file /operations/app_store/cameras/run_fy20/TEST/g2-lerec.laser-relay-cam -t100 -a100,
or -b http https://cdn.spacetelescope.org/archives/images/news/heic1523b.jpg.
for epics:
    ssh acnlinhc
    source /home/cfsd/laster/setup-epics
    imageViewer.py -b epics -w1000,800,8 13PS1:image1:ArrayData
''')
    pargs = parser.parse_args()
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    #`````````````````````````````````````````````````````````````````````````
    if not pargs.black:
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
    
    pargs.backend = pargs.backend.lower()
    if pargs.refreshRate is None:
        pargs.refreshRate = 1000. if pargs.backend == 'file' else 1.

    pname = pargs.pname

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
            try: # get camNm from folder name: /xxx/.../camNm/xxx_camNm_xxx
                #camNm = pargs.pname[0].split('/')[-2]
                pre,camNm = camNm.rsplit('/',1)
                if camNm[:2] == pargs.year[:2]: # the folder is date, use the previus one
                    camNm = pre.split('/')[-1]
                try:    camNm = camNm.split('_')[1]
                except: pass
            except: 
                printw('Could not get camera name from '+str(pargs.pname[0]))
                pass
            print('camNm: '+camNm)
        else:
            camNm = pargs.pname[0]
    camName = [camNm,'?']
    xScale = 1.
    #
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    if camName == 'NoName':
        printw('Camera not recognized: '+str(camNm))
        if pargs.orientation != None:
            sys.exit(1)
    if pargs.controlSystem is None: pargs.controlSystem = 'TEST'
    printi('Control System:'+pargs.controlSystem)
    if pargs.orientation is None: pargs.orientation = 0
    # convert -orientation option to -R and -F combination 
    odict = {0:(pargs.rotate,pargs.flip), 1:(0,'V'), 2:(0,'H'),
      3:(180,None), 4:(90,None), 5:(90,'V'), 6:(90,'H'), 7:(270,None)}
    pargs.rotate, pargs.flip = odict[pargs.orientation]

    # instantiate the imager
    imager = Imager(pname[0],camName)
    if pargs.pixPerMM:
        pixelPmm = pargs.pixPerMM
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
        pvMonitor = PVMonitorAdo(pname[0],refreshRate=pargs.refreshRate)
    elif pargs.backend == 'epics':
        import epics as Backend
        pvMonitor = PVMonitorEpics(pname[0],refreshRate=pargs.refreshRate)
        pargs.fullsize = True
    elif pargs.backend == 'usb':
        try:
            import cv2 as Backend
            print('cv2 imported')
        except ImportError:
            print("ERROR python-opencv must be installed")
            exit(1)        
        pvMonitor = PVMonitorUSB(pname[0],refreshRate=pargs.refreshRate)
        pargs.fullsize = True
    else:
        print('Unknown backend:',pargs.backend)
        exit(8)
    pvMonitor.dbg = pargs.dbg
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
if __name__ == '__main__':

    # enable Ctrl-C to kill application
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)    
    
    app = QtGui.QApplication([])
    main()
    imager.start_imager()
    #
    try:
        # Start Qt event loop unless running in interactive mode or using pyside
        #print('starting QtGUI'+('.' if pargs.file else ' Waiting for data...'))
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            app.exec_()
    except KeyboardInterrupt:
        print('keyboard interrupt: exiting')
    print('Application exit')
    imager.stop()

