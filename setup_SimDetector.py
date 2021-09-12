# Set up the simulated camera
import argparse
from imagin import epicsAccess_caproto as epicsInterface

parser = argparse.ArgumentParser(description=__doc__
,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('camera', nargs ='?', default='sim:det', help=\
'EPICS camera name, e.g. 13SIM1')
pargs = parser.parse_args()

dev = pargs.camera
epicsInterface.set((dev,'cam1:GainX',1))
epicsInterface.set((dev,'cam1:GainY',1))
epicsInterface.set((dev,'cam1:Gain',124))
epicsInterface.set((dev,'cam1:GainRed',1))
epicsInterface.set((dev,'cam1:GainGreen',1))
epicsInterface.set((dev,'cam1:GainBlue',1))
epicsInterface.set((dev,'cam1:SimMode','Peaks'))
epicsInterface.set((dev,'cam1:PeakStartX',300))
epicsInterface.set((dev,'cam1:PeakStartY',300))
epicsInterface.set((dev,'cam1:PeakWidthX',30))
epicsInterface.set((dev,'cam1:PeakWidthY',20))
epicsInterface.set((dev,'cam1:PeakNumX',4))
epicsInterface.set((dev,'cam1:PeakNumY',4))
epicsInterface.set((dev,'cam1:PeakStepX',128))
epicsInterface.set((dev,'cam1:PeakStepY',127))
epicsInterface.set((dev,'cam1:PeakVariation',100))
epicsInterface.set((dev,'cam1:Noise',50))
epicsInterface.set((dev,'cam1:AcquirePeriod',1))
#Does not work#epicsInterface.set((dev,'cam1:Acquire',1))

# Enable data taking
epicsInterface.set((dev,'image1:EnableCallbacks',1))

