# $Id$
# $Date$
# $Author$
# $Revision$

import synthlisa
import lisawp

import sys
import os.path
import time
import string
import re
import math

import Numeric

import pyRXP
# require xmlutils from pyRXP examples
import xmlutils

# begin definitions encoding synthLISA syntax

argumentList = {}
outputList = {}

# give parameter name, default unit, default value (or None if parameter is required)
# this list is used to generate xml from the "initsave" arglist of the SWIGged Wave classes
# and also to generate a Wave class from an xml PlaneWave structure 


# LISA types

argumentList['OriginalLISA'] = ( ('Armlength1','Second','16.6782'),
                                 ('Armlength2','Second','16.6782'),
                                 ('Armlength3','Second','16.6782') )

outputList['OriginalLISA'] = argumentList['OriginalLISA']

                            
argumentList['ModifiedLISA'] = argumentList['OriginalLISA']

outputList['ModifiedLISA'] = argumentList['ModifiedLISA']                            
                            

argumentList['CircularRotating'] = ( ('InitialEta','Radian','0'),
                                     ('InitialXi','Radian','0'),
                                     ('ArmSwitch','1','1'),
                                     ('TimeOffset','Second','0') )

outputList['CircularRotating'] =  ( ('TimeOffset','Second',None),
                                    ('InitialPosition','Radian',None),
                                    ('InitialRotation','Radian',None),
                                    ('Armlength','Second','16.6782'),
                                    ('OrbitRadius','Second','499.004'),
                                    ('OrbitApproximation','String','CircularRigid') )


argumentList['EccentricInclined'] = ( ('InitialEta','Radian','0'),
                                      ('InitialXi','Radian','0'),
                                      ('ArmSwitch','1','1'),
                                      ('TimeOffset','Second','0') )

outputList['EccentricInclined'] = ( ('TimeOffset','Second',None),
                                    ('InitialPosition','Radian',None),
                                    ('InitialRotation','Radian',None),
                                    ('Armlength','Second','16.6782') )


# PyLISA, AllPyLISA, CacheLISA (special treatment of initargs), 
# SimpleLISA, CacheLengthLISA?


# noises

argumentList['PowerLawNoise'] = ( ('Cadence','Second',None),
                                  ('TimeOffset','Second',None),
                                  ('PowerSpectralDensity','(f/Hz)^n/Hz',None),
                                  ('SpectralType','1',None),
                                  ('InterpolatorLength','1',None),
                                  ('PseudoRandomSeed','1',None) )

# need to convert SpectralType from Exponent to String,
# InterpolatorType to Interpolator and InterpolatorWindow

outputList['PowerLawNoise'] = ( ('SpectralType','String',None),
                                ('Cadence','Second',None),
                                ('TimeOffset','Second',None),
                                ('PowerSpectralDensity','(f/Hz)^n/Hz',None),
                                ('PseudoRandomGenerator','String','taus2-gsl1.4'),
                                ('PseudoRandomSeed','1',None),
                                ('Interpolator','String',None),
                                ('InterpolatorWindow','1','None') )

# Wave objects

argumentList['SimpleBinary'] = ( ('Frequency','Hertz',None),
                                 ('InitialPhase','1',None),
                                 ('Inclination','Radian',None),
                                 ('Amplitude','1',None),
                                 ('EclipticLatitude','Radian',None),
                                 ('EclipticLongitude','Radian',None),
                                 ('Polarization','Radian',None) )

outputList['SimpleBinary'] = ( ('EclipticLatitude','Radian',None),
                               ('EclipticLongitude','Radian',None),
                               ('Polarization','Radian',None),
                               ('Frequency','Hertz',None),
                               ('InitialPhase','1',None),
                               ('Inclination','Radian',None),
                               ('Amplitude','1',None) )


argumentList['PNBinary'] = ( ('Mass1','SolarMass',None),
                             ('Mass2','SolarMass',None),
                             ('Frequency','Hertz',None),
                             ('InitialPhase','1',None),
                             ('Distance','Parsec',None),
                             ('Inclination','Radian',None),
                             ('EclipticLatitude','Radian',None),
                             ('EclipticLongitude','Radian',None),
                             ('Polarization','Radian',None),
                             ('IntegrationStep','Second','10'),
                             ('TimeOffset','Second','650') )

outputList['PNBinary'] = ( ('EclipticLatitude','Radian',None),
                           ('EclipticLongitude','Radian',None),
                           ('Polarization','Radian',None),
                           ('TimeOffset','Second',None),
                           ('Mass1','SolarMass',None),
                           ('Mass2','SolarMass',None),
                           ('Frequency','Hertz',None),
                           ('InitialPhase','1',None),
                           ('Distance','Parsec',None),
                           ('Inclination','Radian',None),
                           ('IntegrationStep','Second',None) )


# give translations between synthlisa and XML, and backwards

ObjectToXML = {
    'OriginalLISA': 'OriginalLISA',
    'CircularRotating': 'PseudoLISA',
    'EccentricInclined': 'PseudoLISA',

    'PowerLawNoise': 'PseudoRandomNoise',

    'SimpleBinary': 'GalacticBinary',
    'PNBinary': 'BlackHoleBinary'
}

XMLToObject = {
    # Synthetic LISA objects

    'OriginalLISA': ('OriginalLISA',synthlisa.OriginalLISA),
    'CircularRotating': ('CircularRotating',synthlisa.CircularRotating),
    'EccentricInclined': ('EccentricInclined',synthlisa.EccentricInclined),
    
    'SimpleBinary': ('SimpleBinary',synthlisa.SimpleBinary),
    'PNBinary': ('PNBinary',lisawp.PNBinary),

    # standard lisaXML objects

    'PseudoLISA': ('EccentricInclined',synthlisa.EccentricInclined),    

    'PseudoRandomNoise': ('PowerLawNoise',synthlisa.PowerLawNoise),

    'GalacticBinary': ('SimpleBinary',synthlisa.SimpleBinary),
    'BlackHoleBinary': ('PNBinary',lisawp.PNBinary)
}

# begin definitions encoding XML syntax

minimumParameterSet = {}
optionalParameterSet = {}

# for the moment, no minimum or optional parameter sets for synthLISA objects
# that are not default lisaXML objects
# we'll let the code discover if any parameters should be given that are not

minimumParameterSet['OriginalLISA'] = []
optionalParameterSet['OriginalLISA'] = []

minimumParameterSet['CircularRotating'] = []
optionalParameterSet['CircularRotating'] = []

minimumParameterSet['EccentricInclined'] = []
optionalParameterSet['EccentricInclined'] = []

# LISA objects (aren't these duplications of what we can get from the above?)

minimumParameterSet['PseudoLISA'] = [
    lambda s: 'InitialPosition' in s or 'Polarization',
    lambda s: 'InitialRotation' in s or 'Amplitude',
    lambda s: 'TimeOffset' in s or 'TimeOffset'
]

optionalParameterSet['PseudoLISA'] = [
    lambda s: 'Armlength' in s or ('Armlength',(16.6782,'Second')),
    lambda s: 'ArmSwitch' in s or ('ArmSwitch',(-1.0,'1'))
]

minimumParameterSet['PseudoRandomNoise'] = []
optionalParameterSet['PseudoRandomNoise'] = []

# this construction would replicate the outputList as lambda tests
#  map(lambda p: lambda s: p[0] in s or p[0], outputList['PowerLawNoise'])

# for PlaneWave sources only...

standardSourceParameterSet = [
    lambda s: 'SourceType' in s or 'SourceType',
    lambda s: ( (('EclipticLatitude' in s) and ('EclipticLongitude' in s)) or (('RightAscension' in s) and ('Declination' in s))
                or 'EclipticLatitude/EclipticLongitude or RightAscension/Declination' )
]

# PlaneWave objects

minimumParameterSet['GalacticBinary'] = [
    lambda s: 'Polarization' in s or 'Polarization',
    lambda s: 'Amplitude' in s or 'Amplitude',
    lambda s: 'Inclination' in s or 'Inclination',
    lambda s: 'InitialPhase' in s or 'InitialPhase',
    lambda s: 'Frequency' in s or 'Frequency'
]

optionalParameterSet['GalacticBinary'] = [
    lambda s: 'TimeOffset' in s or ('TimeOffset',('0.0','Second')),
    lambda s: 'FrequencyDot' in s or ('FrequencyDot',('0.0','Hertz/Second')),
    lambda s: 'FrequencyDotDot' in s or ('FrequencyDotDot',('0.0','Hertz/Second^2')),
    lambda s: 'Eccentricity' in s or ('Eccentricity',('0.0','1'))
]

minimumParameterSet['BlackHoleBinary'] = [
    lambda s: 'Mass1' in s or 'Mass1',
    lambda s: 'Mass2' in s or 'Mass2',
    lambda s: 'Frequency' in s or 'Frequency',
    lambda s: 'InitialPhase' in s or 'InitialPhase',
    lambda s: 'Distance' in s or 'Distance',
    lambda s: 'Inclination' in s or 'Inclination'
]

optionalParameterSet['GalacticBinary'] = [
    lambda s: 'IntegrationStep' in s or ('IntegrationStep',('10','Second')),
    lambda s: 'TimeOffset' in s or ('TimeOffset',('650','Second'))
]

# default units

defaultUnits = {
    'EclipticLatitude': 'Radian',
    'EclipticLongitude': 'Radian',
    'Polarization': 'Radian',
    'Amplitude': '1',
    'Inclination': 'Radian',
    'InitialPhase': 'Radian',
    'Frequency': 'Hertz',
    'TimeOffset': 'Second',
    'FrequencyDot': 'Hertz/Second',
    'FrequencyDotDot': 'Hertz/Second^2',
    'Eccentricity': '1',
    'InitialPosition': 'Radian',
    'InitialRotation': 'Radian',
    'Armlength': 'Second'
}

# so far we can convert
#
# - Degree to Radian
# - Degree/Minute/Second (DMS), space-separated, to Radian
# - Hour/Minute/Second (HMS), space-separated, to Radian

def convertUnit(param,unitin,unitout,paramname=''):
    if unitout == unitin:
        return (param,unitin)

    if unitout == "Radian":
        if unitin == "Degree":
            return ((math.pi/180.0) * float(param),'Radian')
        elif unitin == "DMS":
            # space separated
            dms = map(float, param.split())

            return ((math.pi/180.0) * (dms[0] + dms[1]/60.0 + dms[2]/3600.0),'Radian')
        elif unitin == "HMS":
            # space separated
            hms = map(float, param.split())
            
            return ((15.0*math.pi/180.0) * (hms[0] + hms[1]/60.0 + hms[2]/3600.0),'Radian')

    if unitout == 'String':
        if paramname == 'SpectralType':
            if unitin == '1':
                if float(param) == 0.0:
                    return ('WhiteFrequency','String')
                elif float(param) == 2.0:
                    return ('WhitePhase','String')
                elif float(param) == -2.0:
                    return ('WhiteAcceleration','String')

    if unitout == '1':
        if paramname == 'SpectralType':
            if unitin == 'String':
                if param == 'WhiteFrequency':
                    return ('0.0','1')
                elif param == 'WhitePhase':
                    return ('2.0','1')
                elif param == 'WhiteAcceleration':
                    return ('-2.0','1')                

    raise NotImplementedError, "convertUnit(): cannot convert %s from %s to %s" % (param,unitin,unitout)

# so far we can convert:
#
# - RightAscension/Declination to EclipticLatitude/EclipticLongitude 
# - InitialPosition/InitialRotation to InitialEta/InitialXi
# - InitialEta/InitialXi to InitialPosition/InitialRotation

# we may later move this to a more abstract coding structure (unify
# the parsing and unit conversion of parameters, call externally stored
# conversion routines...)

def convertParameters(param,sourceparams):
    if param[0] in ('EclipticLatitude','EclipticLongitude'):    
        try:
            alpha = float(convertUnit(sourceparams['RightAscension'][0],
                                      sourceparams['RightAscension'][1],
                                      'Radian')[0])
                                      
            delta = float(convertUnit(sourceparams['Declination'][0],
                                      sourceparams['Declination'][1],
                                      'Radian')[0])

            # J2000 epoch
            epsilon = (math.pi/180.0) * 23.439291

            beta = math.asin( math.cos(epsilon) * math.sin(delta) -
                              math.sin(epsilon) * math.cos(delta) * math.sin(alpha) )
                              
            coslambd = math.cos(alpha) * math.cos(delta) / math.cos(beta)
            sinlambd = ( (math.sin(delta) - math.cos(epsilon) * math.sin(beta)) /
                         (math.sin(epsilon) * math.cos(beta)) )
                         
            lambd = math.atan2(sinlambd,coslambd)
            if lambd < 0:
                lambd = lambd + 2.0*math.pi

            if param[0] == 'EclipticLatitude':
                return (str(beta),'Radian')
            elif param[0] == 'EclipticLongitude':
                return (str(lambd),'Radian')
        except KeyError:
            raise AttributeError, "convertParameters(): need RightAscension and Declination (in the right units) to return EclipticLatitude and EclipticLongitude"
    elif param[0] in ('InitialPosition','InitialRotation'):
        try:
            eta = float(convertUnit(sourceparams['InitialEta'][0],
                                    sourceparams['InitialEta'][1],
                                    'Radian')[0])

            xi  = float(convertUnit(sourceparams['InitialXi'][0],
                                    sourceparams['InitialXi'][1],
                                    'Radian')[0])

            sw  = float(sourceparams['ArmSwitch'][0])

            if sw > 0:
                raise NotImplementedError, "convertParameters(): LISA eta/xi configuration with sw > 0 not compatible with PseudoLISA"

            initpos = eta
            initrot = xi + initpos - 1.5*math.pi
    
            if param[0] == 'InitialPosition':
                return (str(initpos),'Radian')
            elif param[0] == 'InitialRotation':
                return (str(initrot),'Radian')
        except KeyError:
            raise AttributeError, "convertParameters(): need LISA eta/xi/sw (in the right units) to return InitialPosition and InitialRotation"
    elif param[0] in ('InitialEta','InitialXi','ArmSwitch'):
        try:
            kappa = float(convertUnit(sourceparams['InitialPosition'][0],
                                      sourceparams['InitialPosition'][1],
                                      'Radian')[0])

            lambd = float(convertUnit(sourceparams['InitialRotation'][0],
                                      sourceparams['InitialRotation'][1],
                                      'Radian')[0])
            eta = kappa
            xi = lambd - kappa + 1.5*math.pi

            if param[0] == 'InitialEta':
                return (str(eta),'Radian')
            elif param[0] == 'InitialXi':
                return (str(xi),'Radian')
            elif param[0] == 'ArmSwitch':
                return ('-1.0','1')
        except KeyError:
            raise AttributeError, "convertParameters(): need PseudoLISA InitialPosition and InitialRotation (in the right units) to return eta/xi/sw"
    elif param[0] == 'InterpolatorLength':
        try:
            interptype = convertUnit(sourceparams['Interpolator'][0],
                                     sourceparams['Interpolator'][1],
                                     'String')[0]
    
            if interptype == 'NearestNeighbor':
                return ('0','1')
            elif interptype == 'LinearExtrapolator':
                return ('-1','1')
            elif interptype == 'Linear':
                return ('1','1')
            elif interptype == 'Lagrange': 
                interplen = convertUnit(sourceparams['InterpolatorWindow'][0],
                                        sourceparams['InterpolatorWindow'][1],
                                        '1')[0]
                return (interplen,'1')
            else:
                raise AttributeError
        except KeyError:
            raise AttributeError, "convertParameters(): need Interpolator/InterpolatorWindow (if Interpolator == 'Lagrange') to return InterpolatorLength"
    elif param[0] in ('Interpolator','InterpolatorWindow'):
        try:
            interplen = float(convertUnit(sourceparams['InterpolatorLength'][0],
                                          sourceparams['InterpolatorLength'][1],
                                          '1')[0])

            if param[0] == 'Interpolator':
                 if interplen == 0.0:
                    return ('NearestNeighbor','String')
                 elif interplen == -1.0:
                    return ('LinearExtrapolator','String')
                 elif interplen == 1.0:
                    return ('Linear','String')
                 elif interplen > 1.0:
                    return ('Lagrange','String')            
            elif param[0] == 'InterpolatorWindow':
                if interplen > 1.0:
                    return (str(int(interplen)),'1')
                else:
                    return ('None','1')
        except KeyError:
            raise AttributeError, "convertParameters(): need InterpolatorLength to return Interpolator/InterpolatorWindow"
    else:
        raise AttributeError, "convertParameters(): cannot obtain %s from provided parameters %s" % (param[0],sourceparams)


class writeXML:
    def __init__(self,filename):
        if filename:
            self.f = open(filename,'w')
            self.opened = 1
        else:
            self.f = sys.stdout
            self.opened = -1

    def __del__(self):
        if self.opened == 1:
            self.close()

    def close(self):
        if self.opened == 1:
            self.f.close()
            self.opened = 0

    # handling of XML indentation
    
    stdindent = 4
    indent = ""

    def incind(self):
        self.indent += " " * self.stdindent
    
    def decind(self):
        self.indent = self.indent[0:len(self.indent)-self.stdindent]

    def iprint(self,s):
        if self.opened != 0:
            print >> self.f, self.indent + s

    def doattrs(self,attrs):
        string = ''
    
        if len(attrs) > 0:
            # always put 'Name' first
            
            if 'Name' in attrs.keys():
                string += ' Name="' + attrs['Name'] +'"'
        
            for attr in attrs.keys():
                if attr != 'Name':
                    string += ' ' + attr + '="' + str(attrs[attr]) + '"'
                    
        return string

    def opentag(self,tag,attrs):
        """Open an XML element; take dictionary for attributes""" 
        
        string = '<' + tag + self.doattrs(attrs) + '>'

        self.iprint(string)
        self.incind()
    
    def singletag(self,tag,attrs):
        """Do an XML singleton; take dictionary for attributes"""

        string = '<' + tag + self.doattrs(attrs) + '/>'

        self.iprint(string)

    def coupletag(self,tag,attrs,thevalue):
        """Do inline XML open/close tag"""
    
        string = '<' + tag + self.doattrs(attrs) + '>'     
        string += str(thevalue)
        string += '</' + tag + '>'        
        
        self.iprint(string)

    def closetag(self,tag):
        """Close an XML element"""
    
        string = '</' + tag + '>'

        self.decind()
        self.iprint(string)    

    def content(self,thevalue):
        """Output XML characters"""
        
        # try to keep indentation for multiline content
        for line in str(thevalue).split('\n'):
            # look here for content serialization concerns
            self.iprint(line)

    def outputrxp(self,rxpexp):
        """Output RXP tuple-based expression"""
        
        if rxpexp[2]:
            # have children!

            self.opentag(rxpexp[0],rxpexp[1])

            for elem in rxpexp[2]:
                if type(elem) in (tuple,list):
                    self.outputrxp(elem)
                else:
                    self.content(elem)

            self.closetag(rxpexp[0])
        else:
            # I am a singleton
            self.singletag(rxpexp[0],rxpexp[1])


class typeTimeSeries:
    pass

class typeFrequencySeries:
    pass


def dumpXML(object):
    stdoutXML = lisaXML('')
    
    stdoutXML.outputrxp(stdoutXML.ProcessObjectData(object))
    

class lisaXML(writeXML):
    """Create lisaXML file with metadata (author,date,comments);
       date should be given as ISO-8601, or will be set to today"""

    
    def __init__(self,filename,author='',date='',comments=''):
        """Create lisaXML file with metadata [author,date,comments];
           date should be given as ISO-8601."""
    
        self.filename = filename
    
        # add ".xml" to file if necessary
        if not re.match('.*\.xml$',self.filename) and self.filename != '':
            self.filename += '.xml'
        
        writeXML.__init__(self,self.filename)
        
        self.author = author
        
        # ??? the date should be validated before using it!
        if date:
            self.date = date
        else:
            self.date = time.strftime('%Y-%m-%dT%H:%M:%S%Z',time.localtime())

        self.comments = comments
        
        self.binaryfiles = 0
        self.theLISAData = []
        self.theNoiseData = []
        self.theSourceData = []
        self.theTDIData = []

    
    def doComment(self,comments):
        return ('Comment', {}, [comments])


    def ProcessObjectData(self,object,name='',comments=''):
        """Add an object (LISA, Wave,Noise) to an XML output file."""
    
        # get init arguments and type

        if hasattr(object,'xmlargs'):
            objectarglist = object.xmlargs
        else:
            objectarglist = object.initargs

        if hasattr(object,'xmltype'):
            objecttype = object.xmltype
        else:
            objecttype = object.__class__.__name__

        try:
            defaultarglist = argumentList[objecttype]
        except KeyError:
            raise KeyError, 'lisaXML.ObjectData: unknown object type %s' % objecttype

        # the synthlisa constructor parameters are contained (as a list)
        # in object.xmlargs (preferred) or object.initargs
        # in the order and with the units given by argumentList[objecttype]
        
        # assign the standard parameters to my structure

        objectpars = {}        

        for i in range(len(defaultarglist)):
            param = defaultarglist[i]

            try:
                if param[1]:
                    # if argumentList[objecttype] specifies a standard Unit, use it

                    objectpars[param[0]] = (objectarglist[i],param[1])
                else:
                    # otherwise get it from the element in xmlargs/initargs,
                    # which is a tuple

                    objectpars[param[0]] = objectarglist[i]
            except IndexError:
                if param[2] != None:
                    # if a parameter is missing in xmlargs/initargs, see if
                    # we have a default value and use it
                
                    objectpars[param[0]] = (param[2],param[1])
                else:
                    raise AttributeError, 'lisaXML.ProcessObjectData(): missing internal parameter %s in object %s' % (param[0],object)

        # translate the object to an XML type

        xmltype = ObjectToXML[objecttype]

        # the Params to be output in the XSIL element are given
        # in outputList[objecttype]; this allows for parameter reordering,
        # for converting to the right units, and for setting some default values

        params = []

        for param in outputList[objecttype]:
            # first, see if we have the parameter...

            if not param[0] in objectpars:
                # try obtaining the parameter from a transformation of the other ones

                try:
                    thisparam = convertParameters(param,objectpars)
                except AttributeError:
                    # try using a default output value if there is one
                    
                    if param[2] != None:
                        thisparam = (param[2],param[1])
                    else:
                        raise AttributeError, 'readXML.ProcessObjectData(): missing external parameter(s) %s for object %s' % (param[0],object)
            else:
                thisparam = objectpars[param[0]]

            # convert to the correct units (if we know how to)

            if param[1]:
                thisparam = convertUnit(thisparam[0],thisparam[1],param[1],param[0])

            # add in pyRXP format

            params.append(('Param', {'Name': param[0], 'Unit': thisparam[1]}, [thisparam[0]]))
            
        if comments:
            params.append(self.doComment(comments))

        # create the right kind of XSIL element

        if isinstance(object,synthlisa.Wave):
            xsildict = {'Type': 'PlaneWave'}
            params.append(('Param',{'Name': 'SourceType'},[xmltype]))
        else:
            xsildict = {'Type': xmltype}
           
        if name:
            xsildict['Name'] = name
            
        return ('XSIL', xsildict, params)

        
    # the following calls piggyback on ProcessObjectData
    
    def SourceData(self,source,name='',comments=''):
        """Add a SourceData entry describing a synthlisa source to
        a lisaXML file object."""

        self.theSourceData.append(self.ProcessObjectData(source,name,comments))


    def LISAData(self,lisa,comments=''):
        """Add a LISAData entry describing a synthlisa LISA object to
        a lisaXML file object."""
        
        # allow only one LISA object
        self.theLISAData = [self.ProcessObjectData(lisa,'LISA',comments)]


    def NoiseData(self,object,name='',comments=''):
        """Add a NoiseData entry describing a synthlisa Noise object to
        a lisaXML file object; if called for a TDInoise object, will
        loop over all component noises and dump an XML description for
        each of them."""
        
        if object.__class__.__name__ == 'TDInoise':
            # if we get a TDInoise, do the component noises one by one

            if name:
                comments = comments + '\n(part of %s) ' % name 

            if not self.theLISAData:
                self.LISAData(object.initargs[0],comments)

            # proof-mass noise
            
            map(lambda noise, name: self.NoiseData(noise,name,comments),
                object.initargs[1], 
                ['pm1','pm1s','pm2','pm2s','pm3','pm3s'])
            
            # photodetector noise; note mapping per abstract-format.txt
            
            map(lambda noise, name: self.NoiseData(noise,name,comments),
                object.initargs[2], 
                ['pdm3','pd3','pdm1','pd1','pdm2','pd2'])
                
            # laser noise
                
            map(lambda noise, name: self.NoiseData(noise,name,comments),
                object.initargs[3],
                ['C1','C1s','C2','C2s','C3','C3s'])
        else:
            self.theNoiseData.append(self.ProcessObjectData(object,name,comments))


    def TDIData(self,data,length,cadence,description,offset=0,encoding='Binary',comments=''):
        """Add a TimeSeries object to a lisaXML file object. Here
        'data' is the Numeric array containing the time series
        (simultaneous entries on the same row); 'length' is the desired
        length to be written in the array; 'cadence' is its nominal
        cadence in seconds; 'description' should be a comma-separated
        string listing the TDI observables represented in the time
        series; 'offset' is the nominal initial time for the data
        (currently in seconds); 'encoding' can be 'Binary' for storage
        in a separate binary file, or 'Text' for inline storage in the
        XML file; 'comments' is added to the TimeSeries entry. To
        skip some records at the beginning of the array, use a slicing
        syntax such as data[1:].
        
        Each lisaXML file can contain several TimeSeries objects,
        all contained in the TDIData block; if binary storage is
        requested, each TimeSeries (or FrequencySeries) object
        corresponds to a separate binary file, named file-0.bin,
        file-1.bin, etc., if the main XML file is file.xml."""
           
        # mimick the getobsc call, and get description as comma-separated string
        # ??? provide way to read variable names directly from list
        # Should use different names for phase and frequency TDI variables

        TimeSeries = typeTimeSeries()
        
        try:
            TimeSeries.data = data
            TimeSeries.dim = len(Numeric.shape(data))
            TimeSeries.alength = Numeric.shape(data)[0]
            
            if data.typecode() != 'd':
                raise TypeError
        except:
            print "lisaXML::TDIData: data must be a proper Numeric array of doubles"
            raise TypeError
        
        if TimeSeries.dim == 1:
            TimeSeries.records = 1
        elif TimeSeries.dim == 2:        
            TimeSeries.records = Numeric.shape(data)[1]
        else:
            print "lisaXML::TDIData: behavior undetermined for arrays of this dimension"
            raise NotImplementedError

        TimeSeries.length = length

        if length > TimeSeries.alength:
            print "lisaXML::TDIData: data shorter than requested output length"
            raise IndexError

        TimeSeries.cadence = cadence

        TimeSeries.duration = cadence * length

        TimeSeries.description = description
        TimeSeries.start = offset
        
        TimeSeries.encoding = encoding
   
        TimeSeries.comments = comments
   
        self.theTDIData.append(TimeSeries)

    def writearray(self,data,length,records,description,encoding):
        self.opentag('Array',{'Name': description,'Type': 'double'})
        
        self.coupletag('Dim',{'Name': 'Length'},length)        
        self.coupletag('Dim',{'Name': 'Records'},records)

        # ??? support only remote binary, local ascii
    
        if 'Binary' in encoding:
            # determine the system binary encoding
            
            if sys.byteorder == 'big':
                encoding = 'Binary,BigEndian'
            elif sys.byteorder == 'little':
                encoding = 'Binary,LittleEndian'
            else:
                print 'lisaXML::writedata: cannot determine binary encoding'
                raise NotImplementedError

            # defaulting to remote storage
            # determine binary filename (base filename + ordinal + '.bin')

            binaryfilename = (re.sub('\.xml$','',self.filename) +
                              '-' + str(self.binaryfiles) + '.bin')
            self.binaryfiles += 1

            self.coupletag('Stream',{'Type': 'Remote','Encoding': encoding},
                                    os.path.basename(binaryfilename))

            bfile = open(binaryfilename, 'w')

            if len(data) != length:
                bfile.write(data[0:length].tostring())
            else:
                bfile.write(data.tostring())
                
            bfile.close()
        elif 'Text' in encoding:
            # defaulting to inline storage
        
            self.opentag('Stream',{'Type': 'Local', 'Encoding': 'Text', 'Delimiter': ' '})

            linepattern = '%le ' * records

            for line in range(0,length):
                self.content(linepattern % tuple(data[line]))

            self.closetag('Stream')
        else:
            print 'lisaXML::writedata: encoding not implemented'
            raise NotImplementedError
            
        self.closetag('Array')

    def writeTimeSeries(self,TimeSeries):
        # write out the TimeSeries defined in TimeSeries
        
        self.opentag('XSIL',{'Type': 'TimeSeries',
                             'Name': TimeSeries.description})

        if TimeSeries.comments:
            self.opentag('Comment',{})
            self.content(TimeSeries.comments)
            self.closetag('Comment')
        
        # ??? fix the time types to "s" (XSIL extension) for the moment
        
        self.coupletag('Time',{'Name': 'StartTime','Type': 's'},
                              str(TimeSeries.start))

        self.coupletag('Param',{'Name': 'Cadence','Unit': 's'},
                               str(TimeSeries.cadence))        

        self.coupletag('Param',{'Name': 'Duration','Unit': 's'},
                                  str(TimeSeries.duration))
       
        # ??? use <Column> to define columns (not in XSIL, but in BFD)?

        self.writearray(TimeSeries.data,
                        TimeSeries.length,
                        TimeSeries.records,
                        TimeSeries.description,
                        TimeSeries.encoding)
       
        self.closetag('XSIL')

    def TDISpectraSelfDescribed(self,data,description,encoding='Binary',comments=''):
        """Add a FrequencySeries object to a lisaXML file object.
        Here 'data' is the Numeric array containing the time
        series (simultaneous entries on the same row); 'description'
        is a comma-separated string listing the TDI observables
        represented in the time series; 'encoding' can be 'Binary'
        for storage in a separate binary file, or 'Text' for inline
        storage in the XML file; 'comments' is added to the
        FrequencySeries entry.
        
        The other parameters to TDISpectra are obtained by examining the first
        column of 'data', which is assumed to contain frequencies in Hz."""

        return self.TDISpectra(data,
                               len(data),
                               data[1,0]-data[0,0],
                               description,
                               data[0,0],
                               encoding,
                               comments)

    def TDISpectra(self,data,length,deltaf,description,offset=0,encoding='Binary',comments=''):
        """Add a FrequencySeries object to a lisaXML file object.
        Here 'data' is the Numeric array containing the time series
        (simultaneous entries on the same row); 'length' is the desired
        length of the array to be written; 'deltaf' is the difference
        between successive records [in Hz]; comma-separated string
        listing the TDI observables represented in the time series;
        'offset' is the initial frequency for the spectra (in Hz);
        'encoding' can be 'Binary' for storage in a separate binary
        file, or 'Text' for inline storage in the XML file; 'comments'
        is added to the FrequencySeries entry. To skip some records
        at the beginning of the array, use a slicing syntax such as
        data[1:].
        
        Each lisaXML file can contain several FrequencySeries
        objects, all contained in the TDIData block; if binary storage
        is requested, each FrequencySeries (or TimeSeries)
        object corresponds to a separate binary file, named file-0.bin,
        file-1.bin, etc., if the main XML file is file.xml."""

        # mimick the getobsc call, and get description as comma-separated string
        # ??? provide way to read variable names directly from list
        # Should use different names for phase and frequency TDI variables

        FrequencySeries = typeFrequencySeries()
        
        try:
            FrequencySeries.data = data
            FrequencySeries.dim = len(Numeric.shape(data))
            FrequencySeries.alength = Numeric.shape(data)[0]
            
            if data.typecode() != 'd':
                raise TypeError
        except:
            print "lisaXML::TDISpectra: data must be a proper Numeric array of doubles"
            raise TypeError
        
        if FrequencySeries.dim == 1:
            FrequencySeries.records = 1
        elif FrequencySeries.dim == 2:        
            FrequencySeries.records = Numeric.shape(data)[1]
        else:
            print "lisaXML::TDISpectra: behavior undetermined for arrays of this dimension"
            raise NotImplementedError

        FrequencySeries.length = length

        if length > FrequencySeries.alength:
            print "lisaXML::TDISpectra: data shorter than requested output length"
            raise IndexError

        FrequencySeries.minf = offset * deltaf
        FrequencySeries.maxf = (offset + length - 1) * deltaf
        FrequencySeries.deltaf = deltaf

        FrequencySeries.description = description
        
        FrequencySeries.encoding = encoding
   
        FrequencySeries.comments = comments
   
        self.theTDIData.append(FrequencySeries)

    def writeFrequencySeries(self,FrequencySeries):
        # write out the FrequencySeries defined in FrequencySeries
        
        self.opentag('XSIL',{'Type': 'FrequencySeries',
                             'Name': FrequencySeries.description})

        if FrequencySeries.comments:
            self.opentag('Comment',{})
            self.content(FrequencySeries.comments)
            self.closetag('Comment')
        
        # ??? fix the frequency types to "Hz" (XSIL extension) for the moment
        # ??? provide facility for automatic step specification
        
        self.coupletag('Param',{'Name': 'MinFreq','Type': 'Hz'},
                                      str(FrequencySeries.minf))
        self.coupletag('Param',{'Name': 'MaxFreq','Type': 'Hz'},
                                      str(FrequencySeries.maxf))
        self.coupletag('Param',{'Name': 'DeltaFreq','Type': 'Hz'},
                                      str(FrequencySeries.deltaf))
       
        # ??? use <Column> to define columns (not in XSIL, but in BFD)?
       
        self.writearray(FrequencySeries.data,
                        FrequencySeries.length,
                        FrequencySeries.records,
                        FrequencySeries.description,
                        FrequencySeries.encoding)
       
        self.closetag('XSIL')
    
    
    def close(self):
        """Write the XML file to disk. This happens also on destruction of
        the lisaXML object."""
        
        if self.opened == 0:
            print "lisaXML::close(): File is closed already"
            raise IOError
        
        self.content('<?xml version="1.0"?>')
        self.content('<!DOCTYPE XSIL SYSTEM "bfd.dtd">')

        self.content('<?xml-stylesheet type="text/xsl" href="http://www.vallis.org/lisa-xml.xsl"?>')

        self.opentag('XSIL',{})

        self.coupletag('Param',{'Name': 'Author'},self.author)
        self.coupletag('Param',{'Name': 'GenerationDate', 'Type': 'ISO-8601'},
                               self.date)
        
        if self.comments:
            self.comments += '\n\n'

        self.comments += 'This file produced by Synthetic LISA v. %s\n' % synthlisa.version_short
        self.comments += '(c) 2006 Michele Vallisneri, California Institute of Technology\n'
        self.comments += '---------------------------------------------------------------\n'
        self.comments += synthlisa.version_full

        self.outputrxp(self.doComment(self.comments))

        if self.theLISAData:
            self.opentag('XSIL',{'Type': 'LISAData'})
    
            for object in self.theLISAData:
                self.outputrxp(object)

            self.closetag('XSIL')

        if self.theNoiseData:
            self.opentag('XSIL',{'Type': 'NoiseData'})
            
            for object in self.theNoiseData:
                self.outputrxp(object)
            
            self.closetag('XSIL')

        if self.theSourceData:
            self.opentag('XSIL',{'Type': 'SourceData'})
    
            for object in self.theSourceData:
                self.outputrxp(object)

            self.closetag('XSIL')

        # do the TDIdata objects (first supported)
        
        if len(self.theTDIData) > 0:
            self.opentag('XSIL',{'Type': 'TDIData'})
    
            for object in self.theTDIData:
                if isinstance(object,typeTimeSeries):
                    self.opentag('XSIL',{'Type': 'TDIObservable',
                                 'Name': object.description})
                    self.coupletag('Param',{'Name': 'DataType'},'FractionalFrequency')
                    self.writeTimeSeries(object)
                    self.closetag('XSIL')
                elif isinstance(object,typeFrequencySeries):
                    self.opentag('XSIL',{'Type': 'TDIObservable',
                                 'Name': object.description})
                    self.coupletag('Param',{'Name': 'DataType'},'FractionalFrequency')
                    self.writeFrequencySeries(object)
                    self.closetag('XSIL')                

            self.closetag('XSIL')
                
        self.closetag('XSIL')

        # do the actual writing
        
        writeXML.close(self)


class readXML:
    def __init__(self,filename):
        p = pyRXP.Parser()        

        f = open(filename)
        tree = p(f.read())
        f.close()

        if tree[0] != 'XSIL':
            print 'Not a LISA XSIL file!'
            raise TypeError

        self.directory = os.path.dirname(filename)

        self.tw = xmlutils.TagWrapper(tree)

    def getTime(self,node):
        try:
            # keep Time as string, get Type if provided
            return (str(node),node.Type)
        except AttributeError:
            return (str(node),)
        
    def getParam(self,node):
        try:
            # convert Param to float, get Unit if provided
            return [str(node),node.Unit]
        except AttributeError:
            return [str(node),None]

    def getDim(self,node):
        return int(str(node))

    def processSeries(self,node):
        timeseries = {}
    
        timeseries['Type'] = node.Type
    
        # I suppose 'Name' must be provided!
        timeseries['Name'] = node.Name
        timeseries['Vars'] = node.Name.split(',')
   
        for node2 in node:
            if node2.tagName == 'Time':
                timeseries[node2.Name] = self.getTime(node2)
            elif node2.tagName == 'Param':
                timeseries[node2.Name] = self.getParam(node2)
            elif node2.tagName == 'Array':
                for node3 in node2:
                    if node3.tagName == 'Dim':
                        timeseries[node3.Name] = self.getDim(node3)
                    elif node3.tagName == 'Stream':
                        timeseries['Encoding'] = node3.Encoding
                        
                        if node3.Type == 'Remote':
                            timeseries['Filename'] = str(node3)
                            
                            if 'Binary' in timeseries['Encoding']:
                                # assume length of doubles is 8 (generic?)
                                readlength = 8*timeseries['Length']*timeseries['Records']
    
                                # need to catch reading errors here
                                binaryfile = open(self.directory + '/' + timeseries['Filename'],'r')
                                readbuffer = Numeric.fromstring(binaryfile.read(readlength),'double')
                                binaryfile.close()
                
                                if (('BigEndian' in timeseries['Encoding'] and sys.byteorder == 'little') or
                                    ('LittleEndian' in timeseries['Encoding'] and sys.byteorder == 'big')):
                                    readbuffer = readbuffer.byteswapped()
     
                                if timeseries['Records'] == 1:
                                    timeseries['Data'] = readbuffer
                                else:
                                    timeseries['Data'] = Numeric.reshape(readbuffer,
                                                                         [timeseries['Length'],timeseries['Records']])
                            else:
                                # remote data, not binary
                                raise NotImplementedError
                        elif node3.Type == 'Local':
                            if 'Text' in timeseries['Encoding']:
                                timeseries['Delimiter'] = node3.Delimiter
                                
                                datastring = str(node3)
                                
                                for delchar in timeseries['Delimiter']:
                                    datastring = string.join(datastring.split(delchar),' ')
    
                                # there may be a more efficient way to initialize an array
                                datavalues = map(float,datastring.split())
    
                                if timeseries['Records'] == 1:
                                    timeseries['Data'] = Numeric.array(datavalues,'d')
                                else:
                                    timeseries['Data'] = Numeric.reshape(Numeric.array(datavalues,'d'),
                                                                         [timeseries['Length'],timeseries['Records']])
    
                                # should try different delimiters
                            else:
                                # local data, not textual
                                raise NotImplementedError
    
        return timeseries

    def getTDITimeSeries(self):
        result = []
        
        for node in self.tw:
            # outermost XSIL level container
            
            if node.tagName == 'XSIL':
                if node.Type == 'TDIData':
                    # inside TDIData
        
                    for node2 in node:
                        if node2.tagName == 'XSIL':
                            if node2.Type == 'TDIObservable':
                                for node3 in node2:
                                    if node3.tagName == 'XSIL':
                                        if node3.Type == 'TimeSeries':
                                            # got a TimeSeries!
        
                                            result.append(self.processSeries(node3))
    
        return result

    def getTDIFrequencySeries(self):
        result = []
        
        for node in self.tw:
            # outermost XSIL level container
            
            if node.tagName == 'XSIL':
                if node.Type == 'TDIData':
                    # inside TDIData
                    
                    for node2 in node:
                        if node2.tagName == 'XSIL':
                            if node2.Type == 'TDIObservable':
                                for node3 in node2:
                                    if node3.tagName == 'XSIL':
                                        if node3.Type == 'FrequencySeries':
                                            # got a FrequencySeries!
                                
                                            result.append(self.processSeries(node3))
    
        return result

    def getLISASampledNoise(self):
        result = []
        
        for node in self.tw:
            # outermost XSIL level container
            
            if node.tagName == 'XSIL':
                if node.Name == 'NoiseData':
                    # inside NoiseData
                    
                    for node2 in node:
                        if node2.tagName == 'XSIL':
                            if node2.Type == 'TimeSeries':
                                # got a TimeSeries!
                            
                                result.append(self.processSeries(node2))

        return result

    def getLISASources(self,returnFactory=False):
        result = []

        for node in self.tw:
            if node.tagName == 'XSIL':
                if node.Type == 'SourceData':        
                    # inside SourceData
                    
                    for node2 in node:
                        if node2.tagName == 'XSIL':
                            if node2.Type == 'PlaneWave':
                                r = self.processObject(node2)

                                if returnFactory:
                                    result.append(r)
                                else:
                                    result.append( (r[0])(*(r[1])) )
                            elif node2.Type == 'SampledPlaneWave':
                                raise NotImplementedError, 'readXML.getLISASources(): cannot currently handle SampledPlaneWave objects'
                                
                                # but the idea is to replace r[0] with a class instance that references r[0], and that will
                                # process all the data loading implicit in the TimeSeries parameters, etc., and will return
                                # r[0] applied correctly...

        if returnFactory:
            return LISASourceFactory(result)
        else:
            return result


    def getLISAGeometry(self):
        result = None

        for node in self.tw:
            if node.tagName == 'XSIL':
                if node.Type == 'LISAData':        
                    for node2 in node:
                        if node2.tagName == 'XSIL':
                            if node2.Type == 'PseudoLISA':
                                r = self.processObject(node2)

                                result = (r[0])(*(r[1]))

        return result


    def getLISANoise(self):
        result = []
        
        for node in self.tw:
            if node.tagName == 'XSIL':
                if node.Type == 'NoiseData':        
                    for node2 in node:
                        if node2.tagName == 'XSIL':
                            if node2.Type == 'PseudoRandomNoise':
                                r = self.processObject(node2)

                                ri = (r[0])(*(r[1]))
                                ri.name = r[2]                   

                                result.append(ri)

        return result


    def getTDInoise(self):
        try:
            lisa = self.getLISAGeometry()
            
            if not lisa:
                raise AttributeError
        except:
            raise AttributeError, 'readXML.getTDInoise(): problems reading LISA geometry'

        try:
            noises = self.getLISANoise()
        except:
            raise AttributeError, 'readXML.getTDInoise(): problems reading LISA noises'

        noisedict = {}

        for noise in noises:
            noisedict[noise.name] = noise

        try:
            pmnoises = [noisedict[id] for id in ('pm1','pm1s','pm2','pm2s','pm3','pm3s')]
            pdnoises = [noisedict[id] for id in ('pdm3','pd3','pdm1','pd1','pdm2','pd2')]
            lsnoises = [noisedict[id] for id in ('C1','C1s','C2','C2s','C3','C3s')]
        except KeyError, k:
            raise KeyError, 'readXML.getTDInoise(): undefined noise %s needed for TDInoise' % k

        return synthlisa.TDInoise(lisa,pmnoises,pdnoises,lsnoises)
        

    def processObject(self,node):
        objectparams = {}

        try:
            objectname = node.Name
        except AttributeError:
            objectname = ''

        objectparams = {}

        for node2 in node:
            if node2.tagName == 'Param':
                objectparams[node2.Name] = self.getParam(node2)

        if node.Type == 'PlaneWave':
            for test in standardSourceParameterSet:
                if test(objectparams) != True:
                    raise KeyError, 'readXML.processObject(): need standard parameter(s) %s in source %s' % (test(objectparams),objectname)

            objecttype = objectparams['SourceType'][0]
        else:
            objecttype = node.Type
    
        if not objecttype in minimumParameterSet:
            raise NotImplementedError, 'readXML.processObject(): unknown object type %s for object %s' % (objecttype,objectname)

        # check that minimum parameter set is included

        for test in minimumParameterSet[objecttype]:
            if test(objectparams) != True:
                raise KeyError, 'readXML.processObject(): need parameter(s) %s for object %s of type %s' % (test(objectparams),objectname,objecttype)

        # add the default units if not included and if defined
        
        for param in objectparams:
            if (not param[1]) and (param[0] in defaultUnits):
                param[1] = defaultUnits[param[0]]

        # add default value of optional parameters if not defined

        for test in optionalParameterSet[objecttype]:
            t = test(objectparams)

            if t != True:
                objectparams[t[0]] = t[1]                

        # now convert to a synthlisa object; see if we have it defined
        
        if not objecttype in XMLToObject:
            raise NotImplementedError, 'readXML.processObject(): unknown object %s of type %s' % (objectname,objecttype)

        synthlisatype = XMLToObject[objecttype][0]
            
        # assemble the argument list

        arglist = []

        for param in argumentList[synthlisatype]:
            # first, see if we have the parameter...

            if not param[0] in objectparams:
                # try obtaining the parameter from the other ones

                try:
                    thisparam = convertParameters(param,objectparams)            
                except AttributeError:         
                    raise AttributeError, 'readXML.processObject(): need parameter(s) %s in object %s of type %s' % (param[0],objectname,objecttype)
            else:
                thisparam = objectparams[param[0]]

            # convert to the correct units (if we know how to)

            thisparam = convertUnit(thisparam[0],thisparam[1],param[1],param[0])
                        
            try:
                # we get a string; first try converting to an int...
                
                evalparam = int(thisparam[0])
            except ValueError:
                # if it doesn't work, try a float...

                try:
                    evalparam = float(thisparam[0])
                except ValueError:
                    # if the float does not work, default to str...
                
                    evalparam = str(thisparam[0])

            # this will accept ints and floats, but not strings...

            arglist.append(evalparam)

        return (XMLToObject[objecttype][1],arglist,objectname)


class LISASourceFactory:
    def __init__(self,sourcelist):
        self.sourcelist = sourcelist
        self.sourcenumber = len(sourcelist)
        
    def __getitem__(self,index):
        i = self.sourcelist[index]
        return (i[0])(*(i[1]))

    def __getslice__(self,index1,index2):
        i = self.sourcelist[index]
        return map((i[0])(*(i[1])),self.sourcelist[index1:index2])
    
    def __len__(self):
        return self.sourcenumber

    def __iter__(self):
        return LISASourceIterator(self.sourcelist)

        
class LISASourceIterator:
    def __init__(self,sourcelist):
        self.sourcelist = sourcelist
        self.sourcenumber = len(sourcelist)
        self.last = 0

    def __iter__(self):
        return self

    def next(self):
        if self.last == self.sourcenumber:
            raise StopIteration
        else:
            i = self.sourcelist[self.last]
            self.last = self.last + 1

            return (i[0])(*(i[1]))
