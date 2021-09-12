"""Basic EPICS access API using caproto threading.
The API is similar to cad_io."""
#__version__ = 'v01 2020-06-03'# get() added
#__version__ = 'v02 2020-06-04'# cleanup
#__version__ = 'v03 2020-08-18'#
#__version__ = 'v04b 2020-08-18'# major functions address parameters as (devName,parName)
#__version__ = 'v05a 2020-10-01'# #set() args corrected
#__version__ = 'v06 2020-10-14'# set() is corrected for legalValues
#__version__ = 'v07 2020-11-02'# feature 'W' translated to 'WE'
#TODO: unsubscribe all subscriptions.
#__version__ = 'v3.4.7 2020-11-11'# unsubscribe() corrected, cleanup
__version__ = 'v3.5.1 2021-01-11'# get() accepts *args, *kwargs for copmpatibility with ADO
print(f'{__file__}, {__version__}')

#from timeit import default_timer as timer

from caproto.threading.client import Context
Ctx = Context()

PVCache = {}# cache of PVs
Subscriptions = []

PVC_PV, PVC_CB, PVC_Props = 0, 1, 2
dbg = False
CA_data_type_STRING = 14

def printd(msg):
    if dbg:
        print(f'CPAccess:{msg}')

def init():
    return

def _get_pv(pvName):
    r = PVCache.get(pvName)
    if r is not None:
        pv,*_ = r
    else:
        #printd( f'register pv {pvName}')
        pv,*_ = Ctx.get_pvs(pvName, timeout=2)
        PVCache[pvName] = [pv, None, None]
        _fill_PVCacheProps(pvName)
    return pv

def _fill_PVCacheProps(pvName):
    pv = PVCache[pvName][PVC_PV]
    if pv is None:
        return None
    pvData = pv.read(data_type='time')
    #printd(f'pvData:{pvData}')
    val = pvData.data
    #printd(f'val:{val}')
    if len(val) == 1:
        try:
            # treat it as numpy
            val = pvData.data[0].item()
        except:
            # data is not numpy
            val = pvData.data[0]
    
    # get properties
    #ISSUE: the caproto reports timestamp as float, the precision for float64 
    #presentation is ~300ns
    featureBit2Letter = {1:'R', 2:'WE'}
    featureCode = pv.access_rights
    features = ''
    for bit, letter in featureBit2Letter.items():
        if bit & featureCode:
            features += letter
    pvControl = pv.read(data_type='control')
    #printd(f'pvcontrol {pvName}: {pvControl}')
    datatype = pvControl.data_type
    #printd(f'data_type:{datatype}')
    if datatype == CA_data_type_STRING:# convert text bytearray to str
        val = val.decode()
    props = {'value':val}
    props['timestamp'] = pvData.metadata.timestamp
    props['count'] = len(pvData.data)
    props['features'] = features
    try:    
        props['units'] = pvControl.metadata.units.decode()
        if props['units'] == '':   props['units'] = None
    except: pass

    try:    props['engLow'] = pvControl.metadata.lower_ctrl_limit
    except: pass

    try:    
        props['engHigh'] = pvControl.metadata.upper_ctrl_limit
        if props['engHigh'] == 0.0 and props['engHigh'] == 0.0:
            props['engHigh'], props['engLow'] = None, None
    except: pass

    try:
        props['alarm'] = pvControl.metadata.severity 
        #printd(f'status {pvControl.metadata.severity}')
        if props['alarm'] == 17:# UDF
            props['alarm'] = None
    except: pass

    try:    # legalValues
        enum_strings = pvControl.metadata.enum_strings
        props['legalValues'] = [i.decode() for i in enum_strings]
        props['value'] = props['legalValues'][val]
    except:
        #props['legalValues'] = None
        if 'legalValues' in props:
            del props['legalValues']
    #printd(f'_props {props}')
    PVCache[pvName][PVC_Props] = props

def info(devParName):
    """Abridged PV info"""
    pvName = ':'.join(devParName)
    pv = _get_pv(pvName)
    return {pvName:PVCache[pvName][PVC_Props]}

def get(devParName, *args, **kwargs):
    pvName = ':'.join(devParName)
    pv = _get_pv(pvName)
    pvData = pv.read(data_type='time')
    rDict = _unpack_ReadNotifyResponse(pvName, pvData)
    return rDict

def set(devParValue):
    dev, par, value = devParValue
    #print(f'epicsAccess.set({dev,par,value})')
    pvName = ':'.join((dev,par))
    pv = _get_pv(pvName)
    try: # if PV has legalValues then the value should be index of legalValues
        value = PVCache[pvName][PVC_Props]['legalValues'].index(value)
        #lv = PVCache[pvName][PVC_Props]['legalValues']
        #print(f'lv:{lv}')
    except Exception as e:
        #print(f'in epicsAccess.set. Value not in legalValues: {e}')
        pass
    pv.write(value)
    return 1

def _unpack_ReadNotifyResponse(pvName, pvData):
    #printd('>uRNR')
    val = pvData.data
    if len(val) == 1:
        try:    #it as numpy
            val = pvData.data[0].item()
        except: # it is not numpy
            val = pvData.data[0]
    #printd(f'pvData:{pvData}')
    #printd(f'val:{val}, {pvData.data_type}')
    if pvData.data_type == CA_data_type_STRING:
        val = val.decode()
    #rDict = {'pvname':pvName, 'value':val}
    
    #rDict['timestamp'] = pvData.metadata.timestamp
    ##printd(f'uRNR1:{rDict}')

    legalValues = PVCache[pvName][PVC_Props].get('legalValues')
    #printd(f'uRNR2:{legalValues}')
    if legalValues is not None:
        #rDict['value'] = legalValues[int(val)]
        val = legalValues[int(val)]
    #printd('uRNR3')
    alarm = pvData.metadata.severity
    #rDict['alarm'] = alarm
    #printd(f'pvName:{pvName}')
    key = tuple(pvName.rsplit(':',1))
    #printd(f'key:{key}')
    rDict = {key: {'value':val\
    , 'timestamp':pvData.metadata.timestamp, 'alarm': alarm}}
    #printd(f'<uRNR:{rDict}')
    return rDict

def _callback(subscription, pvData):
    #print(f'>epicsAccess._callback: {pvData})')
    #tMark = [timer(), 0., 0.]
    pvName = subscription.pv.name
    rDict = _unpack_ReadNotifyResponse(pvName, pvData)
    #tMark[1] = timer()
    #printd(f'rDict for {pvName}:{rDict}')
    ##printd(f'PVCache:{PVCache}')
    cache = PVCache.get(pvName)
    #printd(f'cache[{len(cache)}]:{cache}')
    ##printd(f'PVC_CB:{PVC_CB}')
    cb = cache[PVC_CB]
    ##printd(f'c0:{cache[0]}')
    ##printd(f'c1:{cache[1]}')
    #printd(f'cb:{str(cb)}')
    if cb:
        #print(f'epicsAccess call {cb}({rDict})')
        cb(rDict)
    #tMark[2] = timer() - tMark[1]
    #tMark[1] -= tMark[0]
    ##printd(f'caproto cb times {tMark}')# 20-30 uS
    #printd('<callback')
    
def subscribe(callback, devParName):
    pvName = ':'.join(devParName)
    #print(f'>epicsAccess subs({pvName})')
    if not isinstance(pvName, str):
        msg = f'ERROR: Second argument of subscribe() should be a string, not {type(pvName)}'
        raise SystemError(msg)
    pv = _get_pv(pvName)
    subscription = pv.subscribe(data_type='time')
    PVCache[pvName][PVC_CB] = callback
    #printd('>add_callback')
    subscription.add_callback(_callback)
    Subscriptions.append(subscription)
    #print('<subs')

def unsubscribe():
    global Subscriptions
    for subscription in Subscriptions:
        #print(f'>epicsAccess clear subs: {subscription}')
        subscription.clear()
    Subscriptions = []
    
