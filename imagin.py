#!/usr/bin/python
''' Interactive Image Viewer/Analyzer from streamed images or files 
 
Features:
+ Input source: 
    * file, 
    * HTTP link to image, 
    * EPICS PV (requires epics.py).
+ Accepts two-dimensional data, user have to define the array shape (program options: -w --width,height,bits/channel).
+ Wide range of image formats.
+ 16-bit/channel images supported (requires PyPNG or OpenCV).
+ Streaming rotation and flipping.
+ Interactive zooming, panning, rotation.
+ Contrast control: Displays histogram of image data with movable region defining the dark/light levels.
+ ROI and embedded plot for measuring image values.
+ Isocurves. The isocurve level defines the threshold for spot finding.
+ Fast multi-spot finder, reports and logs centroid position and integral of most intense spots in the ROI.
+ Fitted ellipses plotted: the contour ellipses which covers the area of 50% of spot brightness.  
+ Export as PNG,TIFF, JPG..., SVG?, Matplotlib, CSV, HDF5.
+ Interactive python console with access to image data, graphics objects and shell commands (program option: -c).
+ Configuration and reporting in the parameter dock.
+ Image references: save/retrieve image to/from a reference slots.
+ Binary operation on current image and a reference: addition, subtraction.
+ The background subtraction can be achieved by subtraction of a blurred image.

The pyqtgraph is fast in most cases but it only supports 8-bits/channel. 

The OpenCV (option -e cv) is as fast as pyqtgraph, supports 16-bits/channel, 
it is a large package, not widely available. 

The PyPNG (option -e png) supports 16-bit and more per channel for PNG images, 
it is pure python can be downloaded from: https://github.com/drj11/pypng.
The PyPNG is slow on color images.
'''
#__version__ = 'v01 2018-04-05' # created
__version__ = 'v02 2018-04-30' #faster find_spots(), contour ellipses
 
import io
import sys
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
import pyqtgraph as pg
import pyqtgraph.dockarea
import pyqtgraph.console
# Interpret image data as row-major instead of col-major
pg.setConfigOptions(imageAxisOrder='row-major')

import numpy as np
#import skimage.transform as st
from scipy import ndimage
import math

#````````````````````````````Stuff for profiling``````````````````````````````
from timeit import default_timer as timer
import collections
profilingState = collections.OrderedDict() # keeps processing times for diagnostics

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
    
# if graphics is done in callback, then we need this:
QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)

app = QtGui.QApplication([])

#necessary explicit globals
pargs = None
#imager = None
MaxSpotLabels = 16
Prec = 2 # precision for logging and display
DocksShrinked = False
#````````````````````````````Helper Functions`````````````````````````````````        
def printi(msg): print('info: '+msg)
    
def printw(msg): print('WARNING: '+msg)
    
def printe(msg): print('ERROR: '+msg)

def printd(msg): 
    if pargs.dbg: print('dbg: '+msg)

gWidgetConsole = None
def cprint(msg): # print to Console
    if gWidgetConsole:
        gWidgetConsole.write('#'+msg+'\n') # use it to inform the user
    #print(msg)

def cprinte(msg): # print to Console
    if gWidgetConsole:
        gWidgetConsole.write('#ERROR: '+msg+'\n') # use it to inform the user
    printe(msg)

def cprintw(msg): # print to Console
    if gWidgetConsole:
        gWidgetConsole.write('#WARNING: '+msg+'\n') # use it to inform the user
    printw(msg)

def checkPath(path):
# check if path exists, if not, create directory
    try:
        if not os.path.exists(path):
            print('checkPath created new path:',path)
            os.makedirs(path)
    except Exception as e:
        cprinte('in checkPath '+path+' error: '+str(e))

def rgb2gray(data):
    # convert RGB to Grayscale
    if len(data.shape) < 3:
        return data
    else:
        r,g,b = data[:,:,0], data[:,:,1], data[:,:,2]
        if pargs.graysum: # uniform sum
            return r/3 + g/3 + b/3
        else: # using perception-based weighted sum 
            return 0.2989 * r + 0.5870 * g + 0.1140 * b

def imageToArray(img, copy=False, transpose=True):
    """ Corrected pyqtgraph function, supporting Indexed formats.
    Convert a QImage into numpy array. The image must have format RGB32, ARGB32, or ARGB32_Premultiplied.
    By default, the image is not copied; changes made to the array will appear in the QImage as well (beware: if 
    the QImage is collected before the array, there may be trouble).
    The array will have shape (width, height, (b,g,r,a)).
    &RA: fix for Indexed8, take care of possible padding
    """
    nplanes = int(img.byteCount()/img.height()/img.width())
    fmt = img.format()
    ptr = img.bits()
    bpl = img.bytesPerLine() # the bpl is width + len(padding). The padding area is not used for storing anything,
    dtype = np.ubyte
    USE_PYSIDE = False
    if USE_PYSIDE:
        arr = np.frombuffer(ptr, dtype=dtype)
    else:
        ptr.setsize(img.byteCount())
        #arr = np.asarray(ptr)
        arr = np.frombuffer(ptr, dtype=dtype) # this is 30% faster than asarray
    if fmt in (img.Format_Indexed8, 24):
        arr = arr.reshape(img.height(), bpl)
    else:
        arr = arr.reshape(img.height(), img.width(),nplanes)
    
    if copy:
        arr = arr.copy()
        
    if transpose:
        return arr.transpose((1,0,2))
    else:
        return arr

def convert_qimg_to_ndArray(qimg):
    w,h = qimg.width(),qimg.height()
    if w == 0:
        printe('Width unknown, use -w to specify it')
        exit(5)
    planes = qimg.byteCount()/w/h
    t = False
    if planes == 4: 
        return imageToArray(qimg,transpose=t)[...,[2,1,0,3]] # convert BGRA to RGBA
    else:
        return imageToArray(qimg,transpose=t)

def rotate(data,degree):
    import math    
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

def blur(a):
    if len(a.shape) == 2:
        return ndimage.gaussian_filter(a,(2,2)) # 10 times faster than pg.gaussianFilter
    else:
        cprintw('blurring of color images is not implemented yet')
        return a       
    
def sh(s): # console-available metod to execute shell commands
    print(subprocess.Popen(s,shell=True, stdout = None if s[-1:]=="&" else subprocess.PIPE).stdout.read())

#````````````````````````````Spot processing stuff````````````````````````````
def moments(image):
    '''Calculates image moments M00, M10, M01, M20, M02 and M11,
    returns n,s,ss,sxy where n=M00, s=(M10,M01), ss=(M20,M02), sxy=M11'''
    s = np.zeros(2) # sum of samples along axis
    ss = np.zeros(2) # sum of squares of samples along axis
    sxy = 0. #sum of X[i]*Y[i]
    iax = [0.]*2 # index vectors along axis
    n = 0 # number of samples
    for axis in (0,1):
        idx = range(image.shape[axis])
        iax[axis] = np.array(idx,dtype=float)
        oppositeAxis = int(not axis)
        projection = image.sum(axis = oppositeAxis).astype(float)
        psum = projection.sum()
        s[axis] = np.dot(projection,iax[axis])
        ss[axis] = np.dot(projection,iax[axis]**2)
        if axis == 1:
            n += sum(projection)
            for i in idx:
                sxdot = i*np.dot(image[:,i],iax[oppositeAxis])
                sxy += sxdot
    return n,s,ss,sxy

def centroid(data):
    '''Returns first (mean) an second (sigma) moments of 2D array data, 
    if PCC = True, then also returns the Person Correlation Coefficient (PCC).
    Caution on PCC calculation: depending on the numbers involved, it can 
    sometimes be numerically unstable'''
    
    n,s,ss,sxy = moments(data)
    #n = np.sum(data) # this is equally fast
    sigman = np.sqrt(n*ss - s**2)
    means = s/n
    sigmas = sigman/n
    pcc = (n*sxy - s[0]*s[1]) / (sigman[0]*sigman[1])
    return means, sigmas, pcc, n

def find_spots(region, threshold, maxSpots):
    '''find up to maxSpots in the ndarray region and return its centroids and sums
    '''
    rh,rw =  region.shape
    #print('>find_spots:','hw:',(rh,rw),(region[0,0],region[0,rw-1]))

    profile('startFind')
    above_threshold = np.copy(region)
    above_threshold[above_threshold < threshold] = 0
    profile('thresholding')
    
    # now find the objects
    labeled, number_of_objects = ndimage.label(above_threshold)
    profile('labeling')
    
    # sort the labels according to their sums
    sums = ndimage.sum(above_threshold,labeled,index=range(1,number_of_objects+1))
    sumsSorted = sorted(enumerate(sums),key=lambda idx: idx[1],reverse=True)
    labelIndexesSortedBySum = [i[0] for i in sumsSorted]
    profile('sums')
    peak_slices = ndimage.find_objects(labeled)
    profile('find')
    
    # calculate centroids
    centroids = []
    for i in labelIndexesSortedBySum[:maxSpots]:
        islice = peak_slices[i]
        y,x = islice[0].start, islice[1].start
        only_labeled = np.copy(above_threshold[islice])
        only_labeled[labeled[islice] != i+1] = 0 # zero all not belonging to label i+1
        p,w,pcc,s = centroid(only_labeled.T)
        centroids.append(((x+p[0], y+p[1]), w, pcc, s))
    profile('centroids')
    return centroids
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,    
Console = True # needed only when the interactive console is used
if Console: 
    #````````````````````````````Bug fix in pyqtgraph 0.10.0`````````````````````
    import pickle
    class CustomConsoleWidget(pyqtgraph.console.ConsoleWidget):
        ''' Fixing bugs in pyqtgraph 0.10.0:
        Need to rewrite faulty saveHistory()
        and handle exception in loadHistory() if history file is empty.'''
        def loadHistory(self):
            '''Return the list of previously-invoked command strings (or None).'''
            if self.historyFile is not None:
                try:
                    pickle.load(open(self.historyFile, 'rb'))
                except Exception as e:
                    printw('History file '+' not open: '+str(e))

        def saveHistory(self, history):
            '''Store the list of previously-invoked command strings.'''
            #TODO: no sense to provide history argument, use self.input.history instead
            printd('>saveHistory')
            if self.historyFile is not None:
                #bug#pickle.dump(open(self.historyFile, 'wb'), history)
                pickle.dump(history,open(self.historyFile, 'wb'))
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#````````````````````````````PV Monitor Objects for different access systems``
''' The Monitor object is instantiated as:
pvm = PVMonitorXXX(sourceName, reader = readerName)
Derived class must override at two functions: getTimeStamp(), getData() and getSensImageShape()
'''
#````````````````````````````Base PVMonitor class`````````````````````````````
#class PVMonitor(object): # need to be derived from object for super to work
class PVMonitor(QtCore.QThread): # inheritance from QtCore.QThread is needed for qt signals
    def __init__(self):
        super(PVMonitor,self).__init__()
        self.firstEvent = True
        self.data = None
        self.paused = False
        self.timestamp = 0
        self.dbg = None # debugging flag
    
    def convert_png_to_ndArray(self,pngReader):
        width,height,pix,meta = pngReader.read_flat() # slow for RGB images
        printd('meta: '+str(meta))    
        if self.firstEvent:
            self.firstEvent = False
            #self.bitsPerChannel = meta['bitdepth']
            p = meta['planes']
            self.hwp = [height,width,p]
        printd('hwp:'+str(self.hwp))
        return np.reshape(pix, self.hwp if self.hwp[2]>1 else self.hwp[:2])
        
    def getSensImageShape(self): return self.data.shape[:2]
    
    def nextImage(self):
        return None

    def monitor(self):
        '''starts a monitor on the named PV by pvmonitor().'''
        printi('pvmonitor.monitor() is not instrumented') 
        
    def clear(self):
        '''clears a monitor set on the named PV by pvmonitor().'''
        printi('pvmonitor.clear() is not instrumented') 

    def getTimeStamp(self):
        '''returns timestamp, used for polling data delivery'''
        printi('pvmonitor.getTimeStamp() is not instrumented')
        
    def getData(self):
        '''returns the image ndarray, used for polling data delivery'''
        printi('pvmonitor.getData() is not instrumented')
        return []
            
#````````````````````````````PVMonitor for data from a file```````````````````
class PVMonitorFile(PVMonitor):
    def __init__(self,pvname,**kwargs):
        #import glob
        super(PVMonitorFile,self).__init__()
        self.pvsystem = 'File'
        self.qimg = QtGui.QImage() # important to have it persistent
        self.reader = pargs.extract
        #print('pvname:',pvname)
        #fileList = glob.glob(pvname)
        fileList = pvname
        #print('flist:'+str(fileList))
        self.fileIter = iter(fileList)
        self.ts = time.time()
        self.profTime = self.ts
        self.profN = len(fileList)
        self.noMore = False

    def getTimeStamp(self):
        try:
            fname = next(self.fileIter)
        except Exception as e:
            # no more files, do not update timestamp
            if not self.noMore:
                printi('Processed %d images'%self.profN+' in %0.4g s'%(time.time()-self.profTime))
            self.noMore = True
            return self.ts
        printd('image:'+fname)
        self.pvname = fname
        if self.reader == 'qt':
            if not self.qimg.load(fname):
                printe('Loading image '+fname)
                self.ts = time.time()
                self.data = []
                return self.ts
            self.data = convert_qimg_to_ndArray(self.qimg)
        elif self.reader == 'png':
            pngReader = png.Reader(fname)
            self.data = self.convert_png_to_ndArray(pngReader)
        else:
            printw('Reader '+self.reader+' not yet implemented for --access file')
        self.ts = time.time()
        return self.ts
        
    def getData(self):
        return self.data

#````````````````````````````PVMonitor for a HTTP image```````````````````````
class PVMonitorHTTP(PVMonitor):
    def __init__(self,pvname,**kwargs):
        super(PVMonitorHTTP,self).__init__()
        self.pvsystem = 'HTTP'
        self.qimg = QtGui.QImage() # important to have it persistent
        self.req = Backend.get(pvname)
        
    def getTimeStamp(self):
        if not self.timestamp: self.timestamp = time.time()
        return self.timestamp
        
    def getData(self):
        if pargs.dbg:
            print('>http.getdata')
            print('encoding:',self.req.encoding)
            print('status_code:',self.req.status_code)
            print('elapsed:',self.req.elapsed)
            #print(':',self.req.url)
            print('history:',self.req.history)
            print('Content-Type:',self.req.headers['Content-Type'])
            print('cont:',type(self.req.content),len(self.req.content),self.req.content[:20])
        self.qimg.loadFromData(self.req.content)
        data = convert_qimg_to_ndArray(self.qimg)
        return data
        
#````````````````````````````PVMonitor of a an EPICS Process Variable`````````
class PVMonitorEpics(PVMonitor):
    def __init__(self,pvname,**kwargs):
        super(PVMonitorEpics,self).__init__()
        self.pvsystem = 'Epics'
        #epics.ca.replace_printf_handler(self.handle_messages)
        self.pv  = Backend.PV(pvname)
        print('Epics:',self.pv)
        
    def getTimeStamp(self):
        r = self.pv.get_timevars()
        if r == None:
            printe('No PV:'+str(self.pv))
            sys.exit(1)
        return r['timestamp']
        
    def getData(self):
        if pargs.dbg:
            print('>data, dt:%0.4g'%(self.pv.timestamp - self.timestamp))
            #print(pargs)
        self.timestamp = self.pv.timestamp
        data = self.pv.get()
        return data
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#```````````````````````````Imager````````````````````````````````````````````
class Imager(QtCore.QThread): # for signal/slot paradigm the inheritance from QtCore.QThread is necessary
    # define signal on data arrival
    signalDataArrived = QtCore.pyqtSignal(object)
    def __init__(self,pvname):
        super(Imager, self).__init__() # for signal/slot paradigm we need to call the parent init
        self.pvname = pvname
        self.qimg = QtGui.QImage()
        self.hwpb = [0]*4 # height width, number of planes, bits/channel
        self.plot = None # ROI projection plot
        self.roi = None # Region Of Interest
        self.contrast = None # Contrast histogram
        try: self.roiRect = [float(i) for i in pargs.ROI.split(',')]
        except: self.roiRect = None
        self.iso = None # Isocurve object
        self.isoInRoi = False
        self.data = None
        self.grayData = None
        self.mainWidget = None
        self.imageItem = None
        self.dockParRotate = 0 # rotation angle
        self.spotLog = None # enable logging of calculated spot parameters
        self.events = 0 # counter of processed events
        self.threshold = pargs.threshold # threshold for image thresholding
        self.cleanImage = False # show image only
        self.maxSpots = pargs.maxSpots # number of spots to find
        self.spots = [] # calculated spot parameters
        self.spotShapes = [] # spot finder graphics objects
        self.roiArray = [] # array of ROI-selected data
        self.save = False # enable the continuous saving of imahes 
        self.controlSystem = 'TEST' # part of the savePat
        self.refresh = 1./pargs.refreshRate # refresh period in seconds
        self.sleep = 0 # debugging sleep
        self.timestamp = -1 # timestamp from the source
        self.paused = pargs.pause # pause processing    
        if not self.paused:
            self.timestamp = 0 # invalidate timestamp to get one event
        self.nextImage = False
        self.blocked = False # to synchronize event loops in GUI and procThread 
        self.rawData = None # rawData from reader
        self.stopProcThread = False # to stop procThread
        self.ref_diameter = None
        self.ref_X = None
        self.ref_Y = None
        self.pixelPmm = None
        self.spotParToMM = False # convert spot parameters to milimeters
        self.subtraction = None
        self.sizeFactor = 1.
        # connect signal to slot
        self.signalDataArrived.connect(self.process_image) # use mySlot because cannot connect external slot
        profile('start')

    def start(self):
        path = '/tmp/'
        self.savePath = path.replace('TEST',self.controlSystem)+self.pvname+'/'

        checkPath(self.savePath)
        self.refPath = self.savePath+'refs/'
        checkPath(self.refPath)
        self.logPath = self.savePath+'logs/'
        checkPath(self.logPath)
        
        print('Processing thread for '+self.pvname+' started')
        thread = threading.Thread(target=self.procThread)
        thread.start()

    def qMessage(self,text):
        ans = QtGui.QMessageBox.question(self.gl, 'Confirm', text)
        #,                QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        return 1 if ans == QtGui.QMessageBox.Ok else 0

    def open_spotLog(self):
        pname = self.pvname
        if pargs.backend == 'file':
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
          ('pixelPmm', round(self.pixelPmm,4)),
          ('imageCount', self.events),         
          ('spotInfo', 'mean(x,y),sigma(x,y),pcc,sum'),
        ])
        self.spotLog.write(dumps(logdata))#,indent=2))
        self.spotLog.write('\n')
        printd('open,written:'+str(logdata))      

    def show(self):
        ''' Display the widget, called once, when the first image is available'''
        self.win = QtGui.QMainWindow()
        area = pg.dockarea.DockArea()
        self.win.setCentralWidget(area)
        
        # image info
        self.winTitle = 'image:'+self.pvname+' hwpb:'+str(self.hwpb)
        self.win.setWindowTitle(self.winTitle)

        #````````````````````````````parameters data``````````````````````````
        import pyqtgraph.parametertree.parameterTypes as pTypes
        from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
        self.numRefs = 5
        refs = ['ref'+str(i) for i in range(self.numRefs)]
        refreshList = ['1Hz','0.1Hz','10Hz']
        refresh = str(self.refresh)+'Hz'
        if pargs.backend == 'file': refreshList.append('Instant')
        if pargs.userAddon:
            paramAddon = {'name': 'Addon', 'type': 'group', 'children': Addon.controlPane}
        else:
            paramAddon = {'name': 'Addon', 'type': 'group'}
        params = [
            {'name': 'Control', 'type': 'group', 'children': [
                {'name':'Event','type':'int','value':0,
                  'tip':'Accepted events','readonly':True},
                {'name':'Pause', 'type': 'bool', 'value': self.paused,
                  'tip': 'Enable/disable receiving of images, '},
                {'name':'Next', 'type': 'action',
                  'tip':'Process one image from stream'},                
                {'name':'Saving', 'type': 'bool', 'value': False,
                  'tip': 'Enable/disable saving of images'},
                {'name':'View saved images', 'type': 'action',
                  'tip':'View saved images using GraphicMagic, use space/backspace for next/previous image'},
                {'name':'Threshold', 'type': 'float', 'value':self.threshold,
                  'tip': 'Threshold level for spot finding, changed with isoCurve level'},
                {'name':'Clean Image', 'type': 'bool', 'value':self.cleanImage,
                  'tip': 'Show image only'},
            ]},
            {'name': 'Configuration', 'type': 'group', 'children': [
                {'name':'Color', 'type':'list','values':['Native','Gray','Red','Green','Blue'],
                  'tip':'Convert image to grayscale or use only one color channel'},
                {'name':'Refresh Rate', 'type':'list','values':refreshList,
                  'value':refresh,'tip':'Refresh rate'},
                {'name':'Rotate', 'type': 'float', 'value': 0,
                  'tip':'Rotate image view by degree clockwise'},
                {'name':'Subtract', 'type':'list','values':['None',] + refs,
                  'tip':'Streaming subtraction of a reference image inside the ROI'},
                {'name':'RefRing', 'type': 'bool', 'value': False,
                  'tip':'Draw reference ring'},
                 {'name':'Axes in mm', 'type': 'bool', 'value': False,
                  'tip':'convert coordinates to milimeters'},
            ]},
            {'name':'SpotFinder', 'type':'group', 'children': [
                {'name':'MaxSpots', 'type': 'int', 'value':self.maxSpots,
                  'limits':(0,pargs.maxSpots),
                  'tip': 'Max number of spots to find in the ROI'},
                {'name':'Found:', 'type': 'int', 'value':0,'readonroiArrayly':True,
                  'tip': 'Number of spots found in the ROI'},
                {'name':'SpotLog', 'type': 'bool', 'value': False,
                  'tip':'Log the spots parameters to a file'},
            ]},
            paramAddon,
            {'name':'Reference images', 'type': 'group','children': [
                {'name':'Slot','type':'list','values': refs,
                  'tip':'Slot to store/retrieve/ reference image to/from local file "slot#.png"'},
                {'name':'View', 'type': 'action',
                  'tip':'View reference image, use space/backspace for next/previous image'},
                {'name':'Store', 'type': 'action'},
                {'name':'Retrieve', 'type': 'action'},
                {'name':'Add', 'type': 'action'},
                {'name':'Subtract', 'type': 'action'},
            ]},
            {'name':'For Experts', 'type':'group', 'children': [
                {'name':'Blur', 'type': 'action',
                  'tip':'Convert the current image to gray and blur it using gaussian filter with sigma 2'},
                #{'name':'Debug', 'type': 'bool', 'value': False},
                #{'name':'Sleep', 'type': 'float', 'value': 0},
                #{'name':'Test', 'type': 'str', 'value': 'abcd'},
                #{'name':'Debug Action', 'type': 'action'},
            ]},
        ]
        #```````````````````````````Create parameter tree`````````````````````````````
        ## Create tree of Parameter objects
        self.pgPar = Parameter.create(name='params', type='group', children=params)
        
        # Handle any changes in the parameter tree
        def handle_change(param, changes):
            global args
            printd('tree changes:')
            for param, change, itemData in changes:
                path = self.pgPar.childPath(param)
                if path is not None:
                    childName = '.'.join(path)
                else:
                    childName = param.name()
                printd('  parameter: %s'% childName)
                printd('  change:    %s'% change)
                if change == 'options': continue # do not print lengthy text
                printd('  itemData:      %s'% str(itemData))
                printd('  ----------')
            
                parGroupName,parItem = childName,''
                try: 
                    parGroupName,parItem = childName.split('.')
                except: None
                if parGroupName == 'Control':
                    if parItem == 'Pause':
                        printd('Pause')
                        self.paused = itemData
                        self.updateTitle()
                        if not self.paused:
                            self.timestamp = 0 # invalidate timestamp to get one event
                    elif parItem == 'Next':
                        self.nextImage = True
                    elif parItem == 'Saving':
                        self.save = itemData
                        if self.save:
                            self.saveImage()
                            cprint('Saving images to '+self.savePath)
                        else:
                            cprint('Stopped saving to '+self.savePath)
                    elif parItem == 'View saved images':
                        if not os.path.exists(self.savePath):
                            cprinte('opening path: '+self.savePath)
                            return
                        #cmd = ["gm",'display',self.savePath+'IV_'+self.pvname+'_*']
                        cmd = 'xterm -e python imageFrames.py '+self.savePath+'IV_'+self.pvname+'_*'
                        cprint('executing:'+str(cmd))
                        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True) #stderr=subprocess.PIPE)
                    elif parItem == 'Threshold':
                        self.threshold = itemData
                    elif parItem == 'Clean Image':
                        self.cleanImage = itemData
                        if self.cleanImage:
                            self.remove_spotShapes()
                            self.hide_isocurve()
                            self.mainWidget.removeItem(self.roi)
                        else:
                            self.mainWidget.addItem(self.roi)
                        
                elif parGroupName == 'Configuration':               
                    if parItem == 'Color':
                        if itemData == 'Gray':
                            pargs.gray = True
                            self.savedData = self.data
                            self.data = rgb2gray(self.data)
                            self.updateImageItemAndRoi()
                        elif itemData == 'Native':
                            pargs.gray = False
                            self.data = self.savedData
                            self.grayData = rgb2gray(self.data)
                            self.updateImageItemAndRoi()
                        else:
                            cprintw('Color = '+itemData+' reserved for future updates')
                    elif parItem == 'Rotate':
                        self.dockParRotate = float(itemData)
                        self.data = rotate(self.receivedData,self.dockParRotate)
                        self.grayData = rgb2gray(self.data)
                        self.updateImageItemAndRoi()
                        self.updateIsocurve()
                    elif parItem == 'Refresh Rate':
                        self.refresh = {'1Hz':1,'0.1Hz':0.1,'10Hz':10,'Instant':1000}[itemData]
                    elif parItem == 'Subtract':
                        if itemData == 'None':
                            self.subraction = None
                        else:
                            fn = self.refPath+'_'+itemData+'.png'
                            self.subraction = self.shapeData(self.load(fn))
                    elif parItem == 'RefRing':
                        if itemData:
                            self.refRing = self.createRefRing()
                            self.mainWidget.addItem(self.refRing)
                        else:
                            self.mainWidget.removeItem(self.refRing)
                    elif  parItem == 'Axes in mm':
                        self.spotParToMM = itemData
                        self.redraw_axes()                        
                        
                if parGroupName == 'SpotFinder':               
                    if parItem == 'MaxSpots':
                        self.maxSpots = itemData
                    elif parItem == 'mm':
                        self.spotParToMM = itemData
                    elif parItem == 'SpotLog':
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
                            except Exception as e: printe('in spotLog '+str(e))
                            self.spotLog = None
                    
                if parGroupName == 'Reference images':
                    prefix = pargs.logdir+self.pvname+'_'
                    if parItem == 'Slot': pass
                    elif parItem == 'View':
                        cmd = ["gm",'display',prefix+'*']
                        cprint('viewing from '+prefix+'*'+', use Backspace/Space for browsing')
                        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    else:
                        child = self.pgPar.child(parGroupName).child(parItem)
                        slot = self.pgPar.child(parGroupName).child('Slot').value()
                        fn = prefix+slot+'.png'
                        if parItem == 'Store':
                            if slot == 'ref0':
                                self.qMessage('cannot store to '+slot+', it is reserved for fixed background')
                            else:
                                if os.path.exists(fn):
                                    if not self.qMessage('Are you sure you want to overwrite '+slot+'?'):
                                        return
                                img = self.imageItem.qimage
                                if img.mirrored().save(fn,"PNG"): 
                                    cprint('Current image stored to '+fn)
                                else:    cprinte('saving '+fn)
                        else:
                            self.referenceOperation(fn,parItem)
                           
                elif parGroupName == 'Addon':
                    Addon.addonClicked(parItem,itemData)

                elif parGroupName == 'For Experts':
                    if parItem == 'Blur':
                        #TODO: blur color images
                        self.data = blur(self.data)
                        self.updateImageItemAndRoi()
                    elif parItem == 'Debug':
                        pargs.dbg = itemData
                        printi('Debugging is '+('en' if pargs.dbg else 'dis')+'abled')
                    elif parItem == 'Sleep':
                        self.sleep = itemData
                    elif parItem == 'Debug Action':
                        printi('Debug Action pressed')
                        ## print all graphics objects:
                        #print('mw:',self.mainWidget)
                        ##print(dir(self.mainWidget))
                        #items = self.mainWidget.items
                        #for i in items: 
                        #    print(i)

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
        
        for child in self.pgPar.children():
            child.sigValueChanged.connect(valueChanged)
            #child.setKeyboardTracking(False) # does not work
            for ch2 in child.children():
                if not ch2.readonly:
                    ch2.sigValueChanged.connect(valueChanged)
                    #ch2.setKeyboardTracking(False) # does not work 

        ## Create ParameterTree widgets, both accessing the same data
        pgParTree = ParameterTree()
        pgParTree.setParameters(self.pgPar, showTop=False)
        pgParTree.setWindowTitle('Parameter Tree')
        s = (120,0) if not DocksShrinked else (0,0)
        dockPar = pg.dockarea.Dock('dockPar', size=s)
        #dockPar = pg.dockarea.Dock('dockPar', size=(0,0))
        dockPar.addWidget(pgParTree)
        area.addDock(dockPar, 'left')
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

        ## Create docks, place them into the window one at a time.
        ## Note that size arguments are only a suggestion; docks will still have to
        ## fill the entire dock area and obey the limits of their internal widgets.
        h,w = self.data.shape[:2]
        
        imageHSize = float(w)/float(h)*pargs.vertSize
        winHSize = float(w+100)/float(h)*pargs.vertSize # correct for the with of the contrast hist
        dockImage = pg.dockarea.Dock('dockImage - Image', size=(imageHSize,pargs.vertSize))
        area.addDock(dockImage, 'left')
        dockImage.hideTitleBar()
                
        #````````````````````Add widgets into each dock```````````````````````
        # dockImage: a plot area (ViewBox + axes) for displaying the image
        self.gl = pg.GraphicsLayoutWidget()
        self.mainWidget = self.gl.addPlot()
        dockImage.addWidget(self.gl)
        # Item for displaying image data
        #print 'adding imageItem:',self.imageItem.width(),self.imageItem.height()
        self.mainWidget.addItem(self.imageItem)
        self.mainWidget.autoRange(padding=0) # remove default padding
        self.mainWidget.setAspectLocked()
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

        if pargs.hist:
            # Contrast/color control
            self.contrast = pg.HistogramLUTItem()
            self.contrast.setImageItem(self.imageItem)
            #grayData = rgb2gray(self.data)
            
            self.contrast.setLevels(self.grayData.min(), self.grayData.max())
            
            #TODO dockHist#
            #hw = pg.GraphicsView()
            #hw.addItem(self.contrast)
            #dockHist.addWidget(hw)
            
            self.gl.addItem(self.contrast)

        if pargs.iso != 'Off':
        # Isocurve drawing
            if pargs.iso == 'ROI':
                self.isoInRoi = True
                printw('iso == ROI is not fully functional yet')
            self.iso = pg.IsocurveItem(level=0.8, pen='g')
            self.iso.setParentItem(self.imageItem)
            self.iso.setZValue(5)
            # Draggable line for setting isocurve level
            self.isoLine = pg.InfiniteLine(angle=0, movable=True, pen='g')
            self.contrast.vb.addItem(self.isoLine)
            self.contrast.vb.setMouseEnabled(y=False) # makes user interaction a little easier
            self.isoLine.setValue(0.8)
            self.isoLine.setZValue(1000) # bring iso line above contrast controls
            # Connect callback to signal
            #self.isoLine.sigDragged.connect(self.updateIsocurve)
            self.isoLine.sigPositionChangeFinished.connect(self.updateIsocurve)
            #self.updateIso()

        if pargs.roi:
        # Custom ROI for selecting an image region
            s = (1,100) if not DocksShrinked else (0,0)
            dockPlot = pg.dockarea.Dock('dockPlot', size=s)
            area.addDock(dockPlot, 'bottom')
            dockPlot.hideTitleBar()
            self.plot = pg.PlotWidget()
            dockPlot.addWidget(self.plot)

            h,w,p,b = self.hwpb
            rect = [i*g for i,g in zip((w,h,w,h),self.roiRect)]
            self.roi = pg.RectROI(rect[:2], rect[2:], sideScalers=True)
            self.mainWidget.addItem(self.roi)
            self.roi.setZValue(10)  # make sure pargs.roi is drawn above image
            
            # create max number of spot labels
            self.spotLabels = [pg.TextItem('*',color='r',anchor=(0.5,0.5)) 
              for i in range(pargs.maxSpots)]
            for sl in self.spotLabels:
                self.mainWidget.addItem(sl)

            # Connect callback to signal
            self.roi.sigRegionChangeFinished.connect(self.updateRoi)
            self.updateRoi()

        if pargs.console:
        # interactive python console
            s = (0,10) if not DocksShrinked else (0,0)
            dockConsole = pg.dockarea.Dock('dockConsole - Console', size=s, closable=True)
            area.addDock(dockConsole, 'bottom')
            ## Add the console widget
            global gWidgetConsole
            gWidgetConsole = CustomConsoleWidget(
                namespace={'pg': pg, 'np': np, 'plot': self.plot, 'roi':self.roi, #'roiData':meansV,
          'data':self.data, 'image': self.qimg, 'imageItem':self.imageItem, 'pargs':pargs, 'sh':sh},
                historyFile='/tmp/pygpm_console.pcl',text="")
            dockConsole.addWidget(gWidgetConsole)

        self.win.resize(winHSize, pargs.vertSize) 
        self.win.show()
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        
    def redraw_axes(self):
        axx = self.mainWidget.getAxis('bottom')
        axy = self.mainWidget.getAxis('left')
        rx,ry = axx.range, axy.range
        if self.spotParToMM:
            rx[0] = self.pix2mm(rx[0],0.)[0]
            rx[1] = self.pix2mm(rx[1],0.)[0]
            ry[0] = self.pix2mm(0.,ry[0])[1]
            ry[1] = self.pix2mm(0.,ry[1])[1]
            axx.setRange(*rx)
            axy.setRange(*ry)
        else:
            # Not a best method to revert axes to pixels
            # by slightly resizing it...
            size = self.win.size()
            
    def setCalibs(self,pixPerMM=10.,ringDiameter=(100,100),ringCenter=(0,0)):
        self.ref_diameter = ringDiameter
        self.pixelPmm = pixPerMM
        self.ref_X,self.ref_Y = ringCenter

    def createRefRing(self):
        #color = 'b'
        color = (199,234,70) #lime
        pen = pg.mkPen(color=color, width=1, style=QtCore.Qt.DashLine)
        x = self.ref_X/self.sizeFactor + self.hwpb[1]/2.
        y = self.ref_Y/self.sizeFactor + self.hwpb[0]/2.
        dx = self.ref_diameter[0]/self.sizeFactor
        dy = self.ref_diameter[1]/self.sizeFactor
        refring = pg.QtGui.QGraphicsEllipseItem(x-dx/2.,y-dy/2.,dx,dy)
        refring.setPen(pen)
        return refring
        
    def load(self,fn):
    # load image from file
        data = None
        img = QtGui.QImage()
        if img.load(fn):
            #printd('image: '+str((img.width(),img.height())))
            data = convert_qimg_to_ndArray(img.mirrored())
            #printd('loaded:'+str(data))
            #cprint('Loaded image from '+fn)
        else: cprinte('loading '+fn)        
        #print 'load:',data.dtype,self.data.dtype
        return data
        
    def checkPath(self):
    # check if path exists, if not, create directory
        try:
            if not os.path.exists(self.savePath):
                os.makedirs(self.savePath)
        except Exception as e:
            cprinte('creating '+self.savePath+' error: '+str(e))
            self.savePath = '/dev/null/'
    
    def saveImage(self):
        fmt = '_%Y%m%d_%H%M%S.png'
        if self.timestamp:
            strtime = datetime.datetime.fromtimestamp(self.timestamp).strftime(fmt)
        else:
            strtime = time.strftime(fmt)
        self.checkPath()
        fn = self.savePath+'IV_'+self.pvname+strtime
        img = self.imageItem.qimage
        if img == None:
            return
        if not img.mirrored().save(fn,"PNG"): 
            cprinte('saving '+fn)
            # does not work: self.set_dockPar('Control','Saving',False)
    
    def shapeData(self,data):
        if data is not None:
            if len(self.data.shape) == 2:
                data = data[:,:,0] # saved files are always color PNG
        return data

    def referenceOperation(self,fn,operation):
        ''' binary operation with current and restored image '''
        try: 
            data = self.shapeData(self.load(fn))
            if data is None: return
            if operation == 'Retrieve':
                self.data = data
                cprint('retrieved image '+fn)
            elif operation == 'Add':
                self.data = (self.data.astype(int) + data)/2
                cprint('added image '+fn+' to current image')
            elif operation == 'Subtract':
                self.data = self.data.astype(int) - data
                cprint('subtracted image '+fn+' from current image')
            else: pass
            self.grayData = rgb2gray(self.data)
            self.updateImageItemAndRoi()
        except Exception as e: printe(str(e))
        
    def set_dockPar(self,child,grandchild,value):
        self.pgPar.child(child).child(grandchild).setValue(value)

    def stop(self):
        self.stopProcThread = True
        printi('imager stopped')
        try: self.spotLog.close()
        except: pass

    def updateTitle(self):
        self.win.setWindowTitle(('Waiting','Paused')[self.paused]+' '+self.winTitle)
                
    def standard_process(self,offset):
        # standard analysis
        ox,oy = offset
        self.spots = find_spots(self.roiArray,self.threshold,self.maxSpots)
        profile('spotsFound')
        if pargs.profile:
            print(profStates('initRoi','spotsFound'))
            print('FindSpot time: '+profDif('initRoi','spotsFound'))
        spen = pg.mkPen('r')
        logdata = OrderedDict([('time:',time.strftime('%y-%m-%d %H:%M:%S'))])
        #Prec = 2
        for i,spot in enumerate(self.spots):
            #print 'l%i:(%0.4g,%0.4g)'%(i,spot[0],spot[1])
            p,wPx,pcc,s = spot
            xPx,yPx = p[0]+ox, p[1]+oy
            #print('widths:',w,'pcc:',pcc,pcc**2)
            
            # convert to mm if needed and record results
            #x,y = [round((pos-c[0])*c[1],Prec) for pos,c in zip((xPx,yPx),conv)]
            #w = [round(wPx[0]*conv[0][1],Prec), round(wPx[1]*conv[1][1],Prec)]  
            if self.spotParToMM:
                x,y = self.pix2mm(xPx,yPx)
                w = self.pix2mmScale(*wPx)
            else:
                x,y = xPx,yPx
                w = wPx
            x,y = round(x,Prec), round(y,Prec)
            w = [round(w[0],Prec), round(w[1],Prec)]
            # record the results
            logdata['spot_'+str(i)] = ((x,y),w,round(pcc,Prec),s)
            
            if not self.cleanImage:
                # show spotShapes
                self.spotLabels[i].setPos(xPx,yPx)
                # rotate ellipse only when the pcc is significant
                # otherwise the angle is not reliable
                ew,eh = wPx

                def ellipse_prob05(sigmaX,sigmaY,pcc):
                    # Parameters of the probability 0.5 ellipse
                    # Contour ellipse with sum = 0.5 spot sum
                    tan = 2.*pcc*sigmaX*sigmaY/(sigmaX**2 - sigmaY**2)
                    theta = math.atan(tan)/2.
                    cos2 = math.cos(theta)**2
                    sin2 = 1. - cos2
                    cs = 2*cos2 - 1. # cos2 - sin2
                    sigmaX2 = sigmaX**2
                    sigmaY2 = sigmaY**2
                    sigmaU = math.sqrt((cos2*sigmaX2 - sin2*sigmaY2)/cs)
                    sigmaV = math.sqrt((cos2*sigmaY2 - sin2*sigmaX2)/cs)
                    return theta,sigmaU,sigmaV

                if min(wPx) < 0.5:
                    continue
            
                theta,sigmaU,sigmaV = ellipse_prob05(ew,eh,pcc)
                    
                spotShape = pg.QtGui.QGraphicsEllipseItem(xPx-sigmaU,yPx-sigmaV,sigmaU*2.,sigmaV*2.)
                spotShape.setTransformOriginPoint(xPx,yPx)
                if theta:
                    spotShape.setRotation(np.rad2deg(theta))
                
                spotShape.setPen(spen)
                self.mainWidget.addItem(spotShape)
                self.spotShapes.append(spotShape)
                            
        # reset outstanding spotLabels
        for j in range(len(self.spots),len(self.spotLabels)):
            self.spotLabels[j].setPos(0,0)
        self.set_dockPar('SpotFinder','Found:',len(self.spots))
        #self.set_dockPar('SpotFinder','Spots',msg)
        if self.spotLog:
            #print('ld:',logdata)
            if self.save:
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

    def updateRoi(self):
    # callback for handling ROI
        profile('initRoi')
        # the following is much faster than getArrayRegion
        slices = self.roi.getArraySlice(self.grayData,self.imageItem)[0][:2]
        self.roiArray = self.grayData[slices]
        oy,ox = slices[0].start, slices[1].start
        profile('roiArray')
        self.remove_spotShapes()
        
        # find spots using isoLevel as a threshold
        if self.threshold>1 and self.maxSpots>0:
            self.standard_process((ox,oy))
            if pargs.userAddon:
                Addon.process(self.roiArray,self.threshold,(ox,oy),
                  self.mainWidget,logPath=self.logPath,
                  cameraName=self.cameraName,imageFile=self.imageFile)
        
        # plot the ROI histograms
        meansV = self.data[slices].mean(axis=0) # vertical means
        #x = range(len(meansV)+1); s = True
        x = range(len(meansV)); s = False
        #if self.hwpb[2] == 1: # gray image
        if len(self.data.shape) == 2: # gray image
            self.plot.plot(x,meansV,clear=True,stepMode=s)
        else: # color image
            # plot color intensities
            self.plot.plot(x,meansV[:,0],pen='r', clear=True,stepMode=s) # plot red
            self.plot.plot(x,meansV[:,1],pen='g',stepMode=s) # plot green
            self.plot.plot(x,meansV[:,2],pen='b',stepMode=s) # plot blue
            meansVG = self.grayData[slices].mean(axis=0)
            pen = 'k' if pargs.white else 'w'
            self.plot.plot(x,meansVG,pen=pen,stepMode=s) # plot white
        #profile('roiPlot')
        
    def updateIsocurve(self):
    # callback for handling ISO
        printd('>uIsoCurve')
        #profile('init iso')
        #if len(self.roiArray):
        if self.isoInRoi:
            #TODO: need to relocate the isocurves to ROI origin
            self.iso.setData(blur(self.roiArray))
        else:
            self.iso.setData(blur(self.grayData))
        self.threshold = self.isoLine.value()
        #printi('isolevel:'+str(self.threshold))
        self.iso.setLevel(self.threshold)
        #profile('iso')
        self.set_dockPar('Control','Threshold',self.threshold)
        if self.roi: 
            self.updateRoi()
         
    def updateImageItemAndRoi(self):
        self.imageItem.setImage(self.data)
        if self.contrast is not None: self.contrast.regionChanged() # update contrast histogram
        if self.roi: 
            self.updateRoi()

    def process_image(self, **kargs):
        global profilingState
        profile('image')
        printd('uimage:'+str(kargs))
        data = self.rawData
        printd('data:'+str(data[:100]))

        if pargs.width:
            #``````````The source is vector parameter with user-supplied shape
            l = len(data)
            if self.imageItem == None: # do it once
                tokens = pargs.width.split(',')
                w = int(tokens[0])
                h = int(tokens[1])
                bytesPerPixel = l/w/h
                try:    bitsPerPixel = int(tokens[2])
                except: bitsPerPixel = 8*bytesPerPixel                
                # we cannot decide exactly about nPlanes and bytesPerChannel based on bitsPerPixel
                # here is assumption:
                nPlanes = 3 if bitsPerPixel > 16 else 1 #
                bytesPerChannel = bytesPerPixel / nPlanes
                self.hwpb = [h, w, nPlanes,bytesPerChannel]
            
            if self.hwpb[3] > 1: # we need to merge pairs of bytes to integers
                #data = np.array(struct.unpack('<'+str(l/self.hwpb[3] )+'H', data.data),'u2')
                data = struct.unpack('<'+str(l/self.hwpb[3] )+'H', data.data) 
                #data = struct.unpack('<'+str(l/self.hwpb[3] )+'H', bytearray(data)) 
                profile('merge')
            shape = (self.hwpb[:3] if self.hwpb[2]>1 else self.hwpb[:2])
            data = np.reshape(data,shape)
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        #````````````````````Got numpy array from data````````````````````````
        data = data[::-1,...] # flip vertical axis            
        self.receivedData = data # store received data
        printd('array: '+str((data.shape,data.dtype,data)))
        profile('array')
                                        
        self.data = rotate(self.receivedData,self.dockParRotate)
        #profile('rotate'); print('rotation time:'+profDif('array','rotate'))

        if pargs.gray: 
            self.data = rgb2gray(self.data)
            self.grayData = data
        else:
            self.grayData = rgb2gray(self.data)

        #````````````````````Data array is ready for analisys`````````````````
        h,w = self.data.shape[:2]
        if self.imageItem == None: # first event, do the show() only once
            if self.hwpb[0] == 0: # get p,b: number of planes and bytes/channel
                try: p = self.data.shape[2]
                except: p = 1
                b = self.data.dtype.itemsize
                self.hwpb = [h,w,p,b]
            printd('hwpb:'+str(self.hwpb))
            printd('self.array: '+str((self.data.shape,self.data)))
            self.imageItem = pg.ImageItem(self.data)
            self.show()
        else:  # udate data
            #TODO: react on shape change correctly, cannot rely on self.hwpb because of possible rotationg2-cec.laser-relay-cam_ref1.png
            if self.subtraction is not None:
                self.data = self.data.astype(int) - self.subtraction
                self.grayData = rgb2gray(self.data)
            if self.save: self.saveImage() #TODO shouldn't it be after update
            if pargs.backend == 'file':
                self.win.setWindowTitle(pvMonitor.pvname)
            self.updateImageItemAndRoi()
        self.events += 1
        self.set_dockPar('Control','Event',self.events) # concern: time=0.5ms
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,            
        if pargs.profile:
            print('### Total time: '+'%0.3g' % (timer()-profilingState['start']))
            print(profStates('start','finish'))
        self.blocked = False
            
    def procThread(self):
    # data processing thread
        if not pvMonitor:
            printe('pvMonitor did not start, processing stopped')
            sys.exit(10)
        pvMon = pvMonitor
        printi('Processing started, refresh: '+str(self.refresh)+'s, pause='+str(self.paused))
        paused = False
        while not self.stopProcThread:
            profile('start')
            profile('tmp')
            sleepIdle = self.refresh/2.
            if paused and not self.nextImage:
                paused = self.paused
                time.sleep(sleepIdle)
                continue
            paused = self.paused
            self.nextImage = False
            t = pvMon.getTimeStamp() if pvMon else 0.
            profile('getTS')
            dt = t - self.timestamp
            printd('paused:'+str((paused,self.paused,dt,self.timestamp)))
            if dt > self.refresh:
                self.timestamp = t
                self.rawData = pvMon.getData()
                if len(self.rawData) > 0:
                    profile('getData')
                    self.blocked = True
                    self.signalDataArrived.emit('newData')
                    # more correctly, we need to wait for an event here
                    while self.blocked:
                        time.sleep(0.001)
                    #print('blocked for '+str(timer()-s))
                    time.sleep(sleepIdle)
                #time.sleep(self.sleep) # adjustable sleep for debugging
                #print('Data processed')
            else:
                time.sleep(0.5) # time to sleep when no events
        printi('Processing finished')
        
#````````````````````````````Main Program`````````````````````````````````````
pvMonitor = None
Backend = None
#png = None
def main():
    global pargs, imager, pvMonitor, png, Backend
    # parse program arguments
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(
      description = 'Interactive analysis of images from streaming source or files',
      formatter_class=RawTextHelpFormatter)
    parser.add_argument('-d','--dbg', action='store_true', help=
      'Turn on debugging')
    parser.add_argument('-R','--rotate', type=float, default=0, help=
      'Rotate image by ROTATE degree')
    parser.add_argument('-F','--flip', help=
      "Flip image, 'V':vertically, 'H':horizontally")
    #parser.add_argument('-i','--iso', action='store_false', help=
    #  'Disable Isocurve drawing')
    parser.add_argument('-i','--iso',default='Image',help=
      '''Isocurve drawing options: ROI - only in ROI (default), 
      Image - in full image, 
      Off - no isocurves''')
    parser.add_argument('-r','--roi', action='store_false', help=
      'Disable Region Of Interest analysis')
    parser.add_argument('-a','--refreshRate',type=float,default=1.,
      help='Refresh rate [Hz]')    
    parser.add_argument('-f','--fullsize', action='store_true', help=
      'Use full-size full-speed imageM parameter')
    parser.add_argument('-c','--console', action='store_false', help=
      'Disable interactive python console')
    parser.add_argument('-H','--hist', action='store_false', help=
      'Disable histogram with contrast and isocurve contol')
    parser.add_argument('-w','--width', help=
      'For blob data: width,height,bits/pixel i.e 1620,1220,12. The bits/pixel may be omitted for standard images')
    parser.add_argument('-W','--white', action='store_true', help=
      'White background, black foreground for all graphics')
    parser.add_argument('-P','--profile', action='store_true', help=
      'Enable code profiling')
    parser.add_argument('-p','--pause', action='store_true', help=
      'Start in paused state')    
    parser.add_argument('-e','--extract',default='qt',help=
      'Image extractor: qt for QT (default), png for pyPng (for 16-bit+ images)') #cv for OpenCV, 
    parser.add_argument('-g','--gray', action='store_true', help=
      'Show gray image')
    #parser.add_argument('-s','--spot', action='store_true', help='Enable Spot Finder to estimate spot and background parameters inside the ROI. It could be slow on some images')
    parser.add_argument('-G','--graysum', action='store_false', help=
      'Use perceptional color-to-gray conversion, rather than uniform')
    parser.add_argument('-b','--backend', default = 'http', help=
      'Data access backend: file/epics/http') 
    parser.add_argument('-l','--logdir', default = '/tmp/',
      help='Directory for logging and references')
    parser.add_argument('-m','--maxSpots',type=int,default=16,
      help='Maximum number of spots to find')
    parser.add_argument('-t','--threshold',type=float,default=50,
      help='Threshold for spot finding')
    parser.add_argument('-O','--ROI',default='0.05,0.05,0.9,0.9',
      help='ROI rectangle: posX,pozY,sizeX,sizeY, i.e -O0.05,0.05,0.9,0.9')
    parser.add_argument('-v','--vertSize',type=float,default=800,
      help='Vertical size of the display window')
    parser.add_argument('-u','--userAddon', help=
      '''User addon, the addon script name should be prefixed with ivAddon_,
      i.e -uNewAddon will try to import ivAddon_NewAddon.py''')          
    defaultPname = 'heic1523b.jpg'
    parser.add_argument('pname', nargs='*', 
      default=[defaultPname],
      help='''Image stream source. i.e:
or -bepics 13PS1:image1:ArrayData,
or -bfile /tmp/*.png -t100 -a100,
or -bhttp https://cdn.spacetelescope.org/archives/images/news/heic1523b.jpg.''')

    pargs = parser.parse_args()
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    #`````````````````````````````````````````````````````````````````````````
    if pargs.white:
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')

    pname = pargs.pname

    extractor = {'qt':'QT','cv':'OpenCV','png':'PyPng','raw':'raw'}
    if pargs.width: # if width is provided, the pargs.extract should be set to 'raw'
        pargs.extract = 'raw'
                
    print(parser.prog+' items:%i, first:'%len(pname)+pname[0]+\
      ' using '+extractor[pargs.extract]+', version '+__version__)
    if not pargs.hist: pargs.iso = 'Off'
                
    if pargs.extract == 'png':
        import png

    elif pargs.extract == 'cv':
        import cv
               
    # instantiate the imager
    imager = Imager(pname[0])

    # instantiate the data monitor
    # note, only backend file accepts list in pname, all others should use pname[0]
    # TODO: arguments for PVMonitorFile could be omitted, use pargs instead
    if pargs.backend == 'file':
        pvMonitor = PVMonitorFile(pname)
    elif pargs.backend == 'http':
        import requests as Backend
        if pname[0] == defaultPname:
            #pname = 'https://upload.wikimedia.org/wikipedia/commons/thumb/5/5a/Hubble_deep_field.jpg/584px-Hubble_deep_field.jpg'
            #pname = 'https://www.ifa.hawaii.edu/~kaiser/pictures/ntt/a1689.gif'
            #pname = 'https://www.hep.shef.ac.uk/research/dm/images/hubbleDeepField.jpg'
            #pname = ['http://www.dlr.de/dlr/en/Portaldata/1/Resources/bilder/portal/portal_2012_3/scaled/GalaxyClusterAbell1689_sn_l.jpg']
            pname = ['https://cdn.spacetelescope.org/archives/images/news/heic1523b.jpg']
        pvMonitor = PVMonitorHTTP(pname[0])
    elif pargs.backend == 'epics':
        import epics as Backend
        pvMonitor = PVMonitorEpics(pname[0])
    else:
        print('Unknown backend:',pargs.backend)
        exit(8)        
        
    imager.start()
    imager.setCalibs()
    #time.sleep(2)
    pvMonitor.dbg = pargs.dbg
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
if __name__ == '__main__':

    # enable Ctrl-C to kill application
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)    
    
    main()

    try:
        # Start Qt event loop unless running in interactive mode or using pyside
        #print('starting QtGUI'+('.' if pargs.file else ' Waiting for data...'))
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()
    except KeyboardInterrupt:
        print('keyboard interrupt: exiting')
    print('Done')
    imager.stop()




