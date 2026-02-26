"""Image analysis methods, most of the are deployed in the imgMan"""
#TODO: design it as a class, cprint should be class feature
#__version__ = 'v0.1.5 2020-11-13'# tested with iv and imgMan
__version__ = 'v0.1.6a 2021-02-02'# parent_cprint
print(f'imageas {__version__}')

import math
import numpy as np
from scipy import ndimage

X,Y = 0,1

# Dummies so far
def profile(state):
    return
def profileReport():
    return

def printe(msg): print('ERROR in imageas: '+msg)

#reporting warnings through externally supplied functions, if needed 
def _cprint(msg): print('cprint in imageas: '+msg)
CPrint = _cprint
def set_parent_cprint(cprint_func):
    global CPrint
    CPrint = cprint_func
def cprintw(msg):
    CPrint('WARNING  in imageas: '+msg)

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
    """crop 2d array horizontally to have the width divisible by 4
    this is important for PNG encoding"""
    wi = array.shape[1]
    wo = wi & -4
    if wi != wo:
        array = array[:,:wo,...]
    return array

def blur(a,width=2):
    """Blur array a using gaussian filter with width=width"""
    if len(a.shape) == 2:
        return ndimage.gaussian_filter(a,(width,width)) # 10 times faster than pg.gaussianFilter
    else:
        cprintw('blurring of color images is not implemented yet')
        return a       

def filter_prominence(input, size=None, blurwidth=None, footprint=None):
    """Calculate a two-dimensional prominence filter:
output[x,y] = input[x,y] - min(image(footprint(x,y)))
Where  footprint is a boolean array that specifies (implicitly) a shape, 
but also which of the elements within this shape will get passed to the 
filter function.
"""
    blurred = blur(input, blurwidth) if blurwidth else input
    background = ndimage.minimum_filter(blurred, size, footprint)
    return input - background

def ellipse_prob05(sigmaX,sigmaY,pcc):
    """Parameters of the probability 0.5 ellipse
    (Contour ellipse with sum = 0.5 spot sum)"""
    ds = sigmaX**2 - sigmaY**2
    theta = math.atan(2.*pcc*sigmaX*sigmaY/ds)/2. if abs(ds) > 0.001\
      else math.pi/2.
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
    """Returns sum, mean and standard deviation of array of weights
    using np.average. Slightly slower than the dot-based but 
    mathematically stable"""
    mean = np.average(x,weights=weights)
    stDev = np.sqrt(np.average((x - mean)**2, weights=weights))
    return weights.sum(), mean, stDev 
def stdev_dot_based(x,weights):
    """Returns mean and standard deviation of array of weights
    using np.dot. Slightly faster than the average-based.
    Math. stability is not concern on 64-bit machine"""
    s = np.dot(weights,x)
    ss = np.dot(weights,x**2)
    ws = weights.sum()
    return ws, s/ws, np.sqrt(ws*ss - s**2)/ws

def movingAverage (values, window): 
    weights = np.repeat(1.0, window)/window 
    return np.convolve(values, weights, 'valid') 

from scipy.interpolate import splrep, sproot, splev
class FWHM(Exception): pass
def fwhms(y):
    """Determine full-with-half-maximum of a peaked set of points, x and y.
    The function uses a spline interpolation. Reliable after 6-point smoothing
    Performance:0.3ms/800pts"""
    #ts = timer()
    l = len(y)
    #wl = l//20
    wl = 6
    if wl > 1:
        y = movingAverage(y,wl)
    halfMax = (max(y)+min(y))/2.0
    x = np.arange(len(y))
    s = splrep(x, y - halfMax, k=3)
    roots = sproot(s)
    nr = len(roots)
    if nr == 2: # success
        r = abs(roots[1] - roots[0])
        #print('fwhms timing:%.5f'%(timer()-ts))
        return r
    # failure
    #msg = 'fwhms, number of peaks !=2 %i'%nr
    #printw(msg)
    #raise FWHM(msg)
    return 0.

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
        idx = list(range(data.shape[axis]))
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
    """Calculate range of the fitting area
    IMPORTANT! using less than 6*width may cause underestimated baseline."""
    lArr = region.shape[1],region.shape[0]
    fitRange = []
    for axis in (X,Y):   
        halfRange = int(round(width[axis]))*6
        il = max(0,iPosMax[axis] - halfRange)
        ir = min(lArr[axis],iPosMax[axis]+halfRange)-1
        #print('il,ir',il,ir,iPosMax,halfRange)
        if il >= ir:
            return []
        if ir - il < 6:
            cprintw('Not enough area to calculate background for '+'XY'[axis]\
              +', need 6, got %d'%(ir - il))
            fitRange.append((0, lArr[axis]-1))
        else:
            fitRange.append((il,ir))
    return fitRange
        
def estimate_spot_amplitude(region,iPosMax):
    subRegion = region[iPosMax[1]-1:iPosMax[1]+2,iPosMax[0]-1:iPosMax[0]+2]
    return subRegion.mean()

def estimate_base(a):
    """Calculate the base as a mean along the boundary"""
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
    #s = par[0] + par[1]*xx + par[2]*xx**2 # if RankBkg = 3
    s = par[0] # if RankBkg = 1
    for i in range(RankBkg,len(par),3):
        s += par[i+2]*peak_shape(xx-par[i],par[i+1])
    return s
def fit_gaus1D(xr, yr, guess, bounds=(-np.inf, np.inf),axis='?'):
    try:
        fp,pcov = curve_fit(func_sum_of_peaks,xr,yr,
          guess)
          #/*,bounds=bounds)
          #,method='dogbox')#, factor = 1.)#,bounds=bounds)*/
        stdevPars = np.sqrt(np.diag(pcov))/abs(fp)
        sigmaStdev = stdevPars[PWidth]
        if sigmaStdev > 0.2 or math.isnan(sigmaStdev):
            msg = 'fitted sigma'+axis+' is bad, stdev:%.4g'%stdevPars[PWidth]
            cprintw(msg)
            return []
    except Exception as e:
        msg = 'Fit'+axis+' failed: '+str(e)
        cprintw(msg)
        return []
    return fp
def fit_region_1D(fitRegion,posGuess,width,maxLevel,base,fitBaseBounds):
    """One-dimensional gaussian fit of the fitRegion"""
    fittedPars = []
    #v213#pos = []
    #sigmas = []
    pos = np.array(posGuess)
    sigmas = np.array(width)
    #/*print( 'fr',fitRegion.shape)
    for axis in (X,Y):
        projection = np.sum(fitRegion,axis=axis,dtype=np.float)
        nP = fitRegion.shape[axis] # number of points summed in each column
        xyLetter = 'XY'[axis]
        #/*print( 'lp',axis,len(projection),np.argmax(projection))
        
        amp = (maxLevel - base)*fitRegion.shape[axis]/2.5066
        axBase = base*fitRegion.shape[axis]
        guess = [axBase,posGuess[axis],width[axis],amp] # if  RankBkg = 1
        #/*print( 'axBase,posGuess,width,amp',axBase,posGuess,width,amp)
        #print( 'guess'+xyLetter+': '+', '.join(['%.2f'%i for i in guess]))*/

        # bounds may significantly slow the performance, from 7ms to 20ms
        gl = len(guess) - 1
        bounds = ([fitBaseBounds[0]*nP]+[-np.inf]*gl,
          [fitBaseBounds[1]*nP]+[+np.inf]*gl)        
        #/*bounds = ([-np.inf]+[-np.inf]*gl, [+np.inf]+[+np.inf]*gl)
        #print('bounds',bounds)*/
        
        lArr = len(projection)
        xr = list(range(lArr))
        yr = projection
        
        #ts = timer()
        fp = fit_gaus1D(xr,yr,guess,bounds,axis=xyLetter)
        #print( 'fitted: '+', '.join(['%.2f'%i for i in fp])+'. time:%.3f'%(timer()-ts))
        fittedPars.append(fp)
    try:
        pos = np.array((fittedPars[0][PPos],fittedPars[1][PPos]))
        sigmas = np.array((fittedPars[0][PWidth],fittedPars[1][PWidth]))
    except Exception as e:
        #printw('exception in fit_region_1D:'+str(e))
        pass
    return pos,np.abs(sigmas),fittedPars

#````````````````````````````Gauss2D``````````````````````````````````````````
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

def fit_region_2D(fitRegion,pos,width,amp,base,fitBaseBounds):
    """Two-dimensional gaussian fit of the fitRegion"""
    #/*print( '>f2G',fitRegion.shape,pos,width,fitBaseBounds)*/
    fittedPars = []
    pos = np.array(pos)
    width = np.array(width)
    h,w = fitRegion.shape
    x = np.arange(w)
    y = np.arange(h)
    x, y = np.meshgrid(x, y)

    g = 1.5 # moments width is usually underestimated
    guess = [amp-base,pos[0],pos[1],width[0]*g,width[1]*g,base]
    #/*print( 'guess :',['%.2f,'%i for i in guess])*/
    
    if guess[0] < base/2.:
        cprintw('too dim, guess amp:%.2f'%guess[0]+' < base/2:%.2f'%(base/2))
        return pos, width, []

    try:
        ts = timer()
        fp,pcov = curve_fit(twoD_Gaussian, (x, y), fitRegion.ravel(), guess)
        #/*print( 'fitted:',['%.2f,'%i for i in fp],'time:%.3f'%(timer()-ts))*/
        stdevPars = np.sqrt(np.diag(pcov))/abs(fp)
        if stdevPars.max() > 0.2:
            i = np.argmax(stdevPars)
            pnames = ['amp','posX','posY','sX','sY','theta','base']
            msg = 'fitted parameter '+pnames[i]+' bad stdev:%.4g'\
              %stdevPars[i]
            cprintw(msg)
            return pos,width,[]
    except Exception as e:
        msg = '2DFit failed: '+str(e)
        cprintw(msg)
        return pos,width,[]
    sX,sY = abs(fp[3]),abs(fp[4])
    fittedPars = [[fp[5],fp[1],sX,fp[0]],[fp[5],fp[2],sY,fp[0]]]
    return np.array((fp[1],fp[2])), np.array((sX,sY)), fittedPars

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
        #/*print('no objects:',region.shape,np.max(region),np.sum(region))*/
        return [[],[]]
    profile('sums') #TODO: 30% of processing time spent here
    peak_slices = ndimage.find_objects(labeled)
    profile('find') #TODO: 15% of processing time spent here
    
    # sort the labels according to their sums
    sums = ndimage.sum(above_threshold,labeled,
      index=list(range(1,number_of_objects+1)))
    sumsSorted = sorted(enumerate(sums),key=lambda idx: idx[1],reverse=True)
    labelIndexesSortedBySum = [i[0] for i in sumsSorted]
    profile('sort') #TODO: 15% of processing time spent here
    
    # calculate centroids
    centroids = []
    fitRegion = region
    for spotIdx,lbl in enumerate(labelIndexesSortedBySum[:maxSpots]):
        islice = peak_slices[lbl]
        ofsXY = islice[1].start, islice[0].start
        only_labeled = np.copy(above_threshold[islice])
        # zero all not belonging to label lbl+1
        only_labeled[labeled[islice] != lbl+1] = 0
        p,w,pcc,pixSum = centroid(only_labeled.T)
        pos = p + ofsXY

        too_small = min(w) < 0.01# this is, actually check for 0.
        if too_small:
            cprintw( 'Spot too small '+str(w)+', '+str(('','too small')[too_small]))
            continue

        fittedPars = []
        fitRange = []
        if spotIdx==0 and fitBrightest != 'None':
            iPosMax = [int(round(i)) for i in pos]
            fitRange = calculate_fitRange(region,iPosMax,w)
            if len(fitRange) == 0:
                break
            amp = estimate_spot_amplitude(region,iPosMax)
            fitRegion = region[fitRange[1][0]:fitRange[1][1],
              fitRange[0][0]:fitRange[0][1]]
            posGuess = [pos[0] - fitRange[0][0], pos[1] - fitRange[1][0]]
            #/*print( 'fitRange',fitRange,fitRegion.shape,region.shape)
            #print('pos,posg',pos,posGuess)*/
            base = estimate_base(fitRegion)
            brightness = pixSum/w[0]/w[1]
            #/*print( 'pixSum',pixSum,w,brightness)*/
            minBrightness = 1.*base
            if brightness < minBrightness:
                cprintw('spot[%i] brightness %.1f < %.1f'\
                  %(spotIdx,brightness,minBrightness))
                break

            # fit the region
            if fitBrightest == '2D':
                posFit,wFit,fittedPars = \
                  fit_region_2D(fitRegion,posGuess,w,amp,base,fitBaseBounds)
            else:
                posFit,wFit,fittedPars = \
                  fit_region_1D(fitRegion,posGuess,w,amp,base,fitBaseBounds)
            #print( 'fitting time %.4f'%(timer() - ts))
            if len(posFit) == 0:
                break
        profile('fitting')        
        pos[X] += ofs[X]
        pos[Y] += ofs[Y]

        centroids.append((pos, w, pcc, pixSum, fittedPars, fitRange))
    profile('centroids') # 15% of processing time spent here, 09/10:30%
    return centroids, fitRegion

class CodecPIL():
    Image = None
    """The PIL is best codec, as fast as QT, it takes any pixel depth.
    It can be substituted with even faster Pillow-SIMD 2.5 which is clamed to be
    2.5 times faster than Pillow and 10 times faster than ImageMagick.
    """
    def __init__(self):
        from PIL import Image
        self.Image = Image
        self.img = None

    def load(self,fname):
        # takes 40ms for 2.4Mpixel image
        self.img = self.Image.open(fname)
        if self.img.mode == 'P':# palette mode if GIF Image
            self.img = self.img.convert('RGB')
        return np.asarray(self.img)
        
    def loadFromData(self,blob): # not used
        import io
        if isinstance(blob,tuple):
            blob = bytearray(blob)
        self.img = self.Image.open(io.BytesIO(blob))
        if self.img.mode == 'P':# palette mode if GIF Image
            self.img = self.img.convert('RGB')
        return np.asarray(self.img)

    def save(self,fileName,npArray=None):
        #print(f'codec dmax {npArray.max(), npArray.shape, npArray.dtype}')
        if True:#try:
            #ts = timer()
            if npArray is not None:
                vertFlipped = npArray[::-1,...].astype('int32')
                img = self.Image.fromarray(vertFlipped)
            else:
                img = self.img
            img.save(fileName,'PNG')
            #print('saving time',timer()-ts)
        #else:#except Exception as e:
        #    printe('in CodecPIL.save():'+str(e))

Codec = CodecPIL
'''to show an image in ipython:
%gui qt
import pyqtgraph as pg
import imageas
fname = '/cfs/e/LoggerData/run_fy20/fullRun/CEC/Cameras/cs2-inj.yag2-cam/scans/20200710/IM_cs2-inj.yag2-cam_20200710_090501797054.png'
img = imageas.codec.load(fname)
pg.show(img)
'''
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       

