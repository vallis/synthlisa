# $Id$
# $Date$
# $Author$
# $Revision$

import synthlisa

import sys
import os.path
import time
import string
import re
import math

import numpy

import pyRXPU as pyRXP
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
                                 ('InitialPhase','Radian',None),
                                 ('ThetaInclination','Radian',None),
                                 ('Amplitude','1',None),
                                 ('EclipticLatitude','Radian',None),
                                 ('EclipticLongitude','Radian',None),
                                 ('SLPolarization','Radian',None) )

outputList['SimpleBinary'] = ( ('EclipticLatitude','Radian',None),
                               ('EclipticLongitude','Radian',None),
                               ('Polarization','Radian',None),
                               ('Frequency','Hertz',None),
                               ('InitialPhase','Radian',None),
                               ('Inclination','Radian',None),
                               ('Amplitude','1',None) )
    
# let's not support normalization, right now...

argumentList['SampledWave'] = ( ('hp','numpy',None),
                                ('hc','numpy',None),
                                ('Length','1',None),
                                ('Cadence','Second',None),
                                ('Prebuffer','Second',None),
                                ('Normalization','1','1.0'),
                                ('Filtering','Filter','None'),
                                ('InterpolatorLength','1','1'),
                                ('EclipticLatitude','Radian',None),
                                ('EclipticLongitude','Radian',None),
                                ('SLPolarization','Radian',None) )

outputList['SampledWave'] = ( ('EclipticLatitude','Radian',None),
                              ('EclipticLongitude','Radian',None),
                              ('Polarization','Radian',None),
                              ('Interpolator','String',None),
                              ('InterpolatorWindow','1','None') )

# this is special...

outputList['TimeSeries'] = ( ('TimeOffset','Second',None),
                             ('Cadence','Second',None),
                             ('Length','1',None),
                             ('hc','numpy',None),
                             ('hp','numpy',None) )

# give translations between synthlisa and XML, and backwards

ObjectToXML = {
    'OriginalLISA': 'OriginalLISA',
    'CircularRotating': 'PseudoLISA',
    'EccentricInclined': 'PseudoLISA',

    'PowerLawNoise': 'PseudoRandomNoise',

    'SimpleBinary': 'GalacticBinary',
    
    'SampledWave': 'SampledPlaneWave'
}

XMLToObject = {
    # Synthetic LISA objects

    'OriginalLISA': ('OriginalLISA',synthlisa.OriginalLISA),
    'CircularRotating': ('CircularRotating',synthlisa.CircularRotating),
    'EccentricInclined': ('EccentricInclined',synthlisa.EccentricInclined),
    
    'SimpleBinary': ('SimpleBinary',synthlisa.SimpleBinary),

    # standard lisaXML objects

    'PseudoLISA': ('EccentricInclined',synthlisa.EccentricInclined),    

    'PseudoRandomNoise': ('PowerLawNoise',synthlisa.PowerLawNoise),

    'GalacticBinary': ('SimpleBinary',synthlisa.SimpleBinary),
    
    'SampledPlaneWave': ('SampledWave',synthlisa.SampledWave)
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

def makeminimum(parlist):
    return map(lambda p: lambda s: p in s or p,parlist) 

def makeoptional(parlist):
    return map(lambda p: lambda s: p[0] in s or p,parlist) 


minimumParameterSet['PseudoLISA'] = makeminimum(['InitialPosition',
                                                 'InitialRotation',
                                                 'TimeOffset'])

optionalParameterSet['PseudoLISA'] = makeoptional([('Armlength',(16.6782,'Second')),
                                                   ('ArmSwitch',(-1.0,'1'))])
minimumParameterSet['PseudoRandomNoise'] = []
optionalParameterSet['PseudoRandomNoise'] = []

# this construction would replicate the outputList as lambda tests
#  map(lambda p: lambda s: p[0] in s or p[0], outputList['PowerLawNoise'])

# for PlaneWave sources only...

standardSourceParameterSet = [
    lambda s: ( (('EclipticLatitude' in s) and ('EclipticLongitude' in s)) or (('RightAscension' in s) and ('Declination' in s))
                or 'EclipticLatitude/EclipticLongitude or RightAscension/Declination' ),
    lambda s: 'Polarization' in s or 'Polarization'
]

# PlaneWave objects

minimumParameterSet['GalacticBinary'] = makeminimum(['Polarization',
                                                     'Amplitude',
                                                     'Inclination',
                                                     'InitialPhase',
                                                     'Frequency'])

optionalParameterSet['GalacticBinary'] = makeoptional([('TimeOffset',('0.0','Second')),
                                                       ('FrequencyDot',('0.0','Hertz/Second')),
                                                       ('FrequencyDotDot',('0.0','Hertz/Second^2')),
                                                       ('Eccentricity',('0.0','1'))])

minimumParameterSet['SampledPlaneWave'] = makeminimum(['TimeOffset',
                                                       'Cadence',
                                                       'Duration'])

optionalParameterSet['SampledPlaneWave'] = []

# default units

defaultUnits = {
    'EclipticLatitude': 'Radian',
    'EclipticLongitude': 'Radian',
    'Polarization': 'Radian',
    'SLPolarization': 'Radian',
    'Amplitude': '1',
    'Inclination': 'Radian',
    'InitialPhase': 'Radian',
    'Frequency': 'Hertz',
    'TimeOffset': 'Second',
    'FrequencyDot': 'Hertz/Second',
    'FrequencyDotDot': 'Hertz/Second^2',
    'Eccentricity': '1',
    'InitialAngularOrbitalPhase': 'Radian',
    'CoalescenceTime': 'Second',
    'InitialPosition': 'Radian',
    'InitialRotation': 'Radian',
    'Armlength': 'Second'
}

from convertunit import convertParameters, convertUnit

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

    # it would be better to redefine this function in lisaXML
    # to implement the binaryData routine

    def content(self,thevalue):
        """Output XML characters"""
        
        if isinstance(thevalue,binaryData):
            filename = re.sub('\.xml$','',self.filename) + '-' + str(self.binaryfiles) + '.bin'
            self.iprint(os.path.basename(filename))

            thevalue.dumpdata(filename)          
            self.binaryfiles += 1
        else:
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


class binaryData:
    def __init__(self,thedata,length):
        self.thedata = thedata

        self.records = len(thedata)
        self.length = length

    def dumpdata(self,filename):
        buffer = numpy.zeros([self.length,self.records],'d')

        for i in range(self.records):
            buffer[:,i] = self.thedata[i][0:self.length]

        bfile = open(filename,'w')
        bfile.write(buffer.tostring())              
        bfile.close()


class fileData:
    def __init__(self,filename,basename,length,records,encoding,index):
        self.filename = filename
        self.basename = basename
        
        self.length = length
        self.records = records
        self.encoding = encoding
        
        self.index = index



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

        if isinstance(object,synthlisa.Wave) and not xmltype == 'SampledPlaneWave':
            xsildict = {'Type': 'PlaneWave'}
            params.append(('Param',{'Name': 'SourceType'},[xmltype]))
        else:
            xsildict = {'Type': xmltype}
           
        if name:
            xsildict['Name'] = name
        
        # add a TimeSeries if needed
        
        if objecttype == 'SampledWave':
            params.append(self.doSampledWaveTimeSeries(objectpars))
        
        return ('XSIL', xsildict, params)


    def doSampledWaveTimeSeries(self,myobjectpars):
        objectpars = {}

        for param in outputList['TimeSeries']:
            if not param[0] in myobjectpars:
                try:
                    thisparam = convertParameters(param,myobjectpars)
                except AttributeError:
                    if param[2] != None:
                        thisparam = (param[2],param[1])
                    else:
                        raise AttributeError, 'readXML.doSampledWaveTimeSeries(): missing external parameter(s) %s for object %s' % (param[0],object)
            else:
                thisparam = myobjectpars[param[0]]

            if param[1]:
                thisparam = convertUnit(thisparam[0],thisparam[1],param[1],param[0])

            objectpars[param[0]] = thisparam

        children = []

        for name in ('TimeOffset','Cadence'):
            children.append( ('Param',
                              {'Name': name, 'Unit': objectpars[name][1]},
                              [objectpars[name][0]]) )

        children.append( ('Param',
                         {'Name': 'Duration', 'Unit': objectpars['Cadence'][1]},
                         [str(float(objectpars['Cadence'][0]) * float(objectpars['Length'][0]))]) )


        arraycontent = []
        
        arraycontent.append( ( 'Dim', {'Name': 'Length'}, [objectpars['Length'][0]] ) )
        arraycontent.append( ( 'Dim', {'Name': 'Records'}, [str(2)] ) )        

        # check bigendian, littleendian

        if sys.byteorder == 'big':
            encoding = 'Binary,BigEndian'
        elif sys.byteorder == 'little':
            encoding = 'Binary,LittleEndian'

        arraycontent.append( ('Stream',
                             {'Type': 'Remote', 'Encoding': encoding},
                             [ binaryData( (objectpars['hp'][0], objectpars['hc'][0]),
                                           int(objectpars['Length'][0]) ) ]
                           ) )
        
        children.append( ('Array',
                         {'Name': 'hp,hc', 'Type': 'double', 'Unit': '1'},
                         arraycontent) ) 

        return ('XSIL', {'Name': 'hp,hc', 'Type': 'TimeSeries'}, children)
        
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


    def NoiseData(self,object,name='',comments='',hideseed=0):
        """Add a NoiseData entry describing a synthlisa Noise object to
        a lisaXML file object; if called for a TDInoise object, will
        loop over all component noises and dump an XML description for
        each of them."""
        
        if object.__class__.__name__ == 'TDInoise':
            # if we get a TDInoise, do the component noises one by one

            if name:
                if comments:
                    comments = comments + '\n(part of %s) ' % name 
                else:
                    comments = '(part of %s) ' % name

            if not self.theLISAData:
                self.LISAData(object.initargs[0],comments)

            # proof-mass noise
            
            map(lambda noise, name: self.NoiseData(noise,name,comments,hideseed),
                object.initargs[1], 
                ['pm1','pm1s','pm2','pm2s','pm3','pm3s'])
            
            # photodetector noise; note mapping per abstract-format.txt
            
            map(lambda noise, name: self.NoiseData(noise,name,comments,hideseed),
                object.initargs[2], 
                ['pdm3','pd3','pdm1','pd1','pdm2','pd2'])
                
            # laser noise
                
            map(lambda noise, name: self.NoiseData(noise,name,comments,hideseed),
                object.initargs[3],
                ['C1','C1s','C2','C2s','C3','C3s'])
        elif object.__class__.__name__ == 'SumSignal':
            self.NoiseData(object.initargs[0],name,comments,hideseed)
            self.NoiseData(object.initargs[1],name,comments,hideseed)
        elif object.__class__.__name__ == 'NoSignal':
            return
        else:
            if hasattr(object,'xmltype') and object.xmltype == 'PowerLawNoise' and hideseed == 1:
                # hide pseudorandom seed if so requested  
                object.xmlargs[5] = 0
            
            self.theNoiseData.append(self.ProcessObjectData(object,name,comments))


    def TDIData(self,data,length,cadence,description,offset=0,encoding='Binary',comments=''):
        """Add a TimeSeries object to a lisaXML file object. Here
        'data' is the numpy array containing the time series
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
            TimeSeries.dim = len(numpy.shape(data))
            TimeSeries.alength = numpy.shape(data)[0]
            
            if data.dtype.char != 'd':
                raise TypeError
        except:
            print "lisaXML::TDIData: data must be a proper numpy array of doubles"
            raise TypeError
        
        if TimeSeries.dim == 1:
            TimeSeries.records = 1
        elif TimeSeries.dim == 2:        
            TimeSeries.records = numpy.shape(data)[1]
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
                
        self.coupletag('Param',{'Name': 'TimeOffset','Type': 'Second'},
                              str(TimeSeries.start))

        self.coupletag('Param',{'Name': 'Cadence','Unit': 'Second'},
                               str(TimeSeries.cadence))        

        self.coupletag('Param',{'Name': 'Duration','Unit': 'Second'},
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
        Here 'data' is the numpy array containing the time
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
        Here 'data' is the numpy array containing the time series
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
            FrequencySeries.dim = len(numpy.shape(data))
            FrequencySeries.alength = numpy.shape(data)[0]
            
            if data.dtype.char != 'd':
                raise TypeError
        except:
            print "lisaXML::TDISpectra: data must be a proper numpy array of doubles"
            raise TypeError
        
        if FrequencySeries.dim == 1:
            FrequencySeries.records = 1
        elif FrequencySeries.dim == 2:        
            FrequencySeries.records = numpy.shape(data)[1]
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
        self.content('<!DOCTYPE XSIL SYSTEM "http://www.vallis.org/lisa-xml.dtd">')

        self.content('<?xml-stylesheet type="text/xsl" href="lisa-xml.xsl"?>')
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
        lines = f.read()
        f.close()
        
        try:
            tree = p(lines)
        except pyRXP.error:
            print "XML validation error! (Or perhaps I couldn't access the DTD)."
            print "I'll try to use the file anyway by removing the DTD..."

            lines = re.sub('<!DOCTYPE XSIL SYSTEM ".*">','',lines)
            tree = p(lines)

        if tree[0] != 'XSIL':
            print 'Not a LISA XSIL file!'
            raise TypeError

        self.directory = os.path.dirname(filename)

        self.tw = xmlutils.TagWrapper(tree)

    def close(self):
        pass

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
                                if self.directory:
                                    binaryfile = open(self.directory + '/' + timeseries['Filename'],'r')
                                else: 
                                    binaryfile = open(timeseries['Filename'],'r')
                                                            
                                readbuffer = numpy.fromstring(binaryfile.read(readlength),'double')
                                binaryfile.close()
                
                                if (('BigEndian' in timeseries['Encoding'] and sys.byteorder == 'little') or
                                    ('LittleEndian' in timeseries['Encoding'] and sys.byteorder == 'big')):
                                    readbuffer = readbuffer.byteswap()
     
                                if timeseries['Records'] == 1:
                                    timeseries['Data'] = readbuffer
                                else:
                                    timeseries['Data'] = numpy.reshape(readbuffer,
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
                                    timeseries['Data'] = numpy.array(datavalues,'d')
                                else:
                                    timeseries['Data'] = numpy.reshape(numpy.array(datavalues,'d'),
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
                            if node2.Type in ('PlaneWave','SampledPlaneWave'):
                                r = self.processObject(node2)

                                if returnFactory:
                                    result.append(r)
                                else:
                                    result.append(MakeObject(r))

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

                                result = MakeObject(r)

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

                                result.append(MakeObject(r))

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
            if noise.name in noisedict:
                # handle composite noise
                
                noisedict[noise.name] = synthlisa.SumSignal(noisedict[noise.name],noise)
            else:
                noisedict[noise.name] = noise

        # undefined noises are replaced with NoSignal()

        checknoise = lambda id: id in noisedict and noisedict[id] or synthlisa.NoSignal()

        pmnoises = [checknoise(id) for id in ('pm1','pm1s','pm2','pm2s','pm3','pm3s')]
        pdnoises = [checknoise(id) for id in ('pdm3','pd3','pdm1','pd1','pdm2','pd2')]
        lsnoises = [checknoise(id) for id in ('C1','C1s','C2','C2s','C3','C3s')]

        return synthlisa.TDInoise(lisa,pmnoises,pdnoises,lsnoises)
        

    def processArray(self,node,objectparams):
        for node2 in node:
            if node2.tagName == 'Dim':
                objectparams[node2.Name] = (self.getDim(node2),'1')

        for node2 in node:
            if node2.tagName == 'Stream':
                try:
                    encoding = node2.Encoding
                except:
                    raise AttributeError, 'readXML.processArray(): encoding/type not specified in stream'

                if node2.Type == 'Remote':
                    vars = map(lambda s: string.lstrip(string.rstrip(s)),string.split(node.Name,','))

                    for v in range(len(vars)):
                        objectparams[vars[v]] = (fileData(self.directory + '/' + str(node2),str(node2),
                                                          int(objectparams['Length'][0]),int(objectparams['Records'][0]),encoding,v),
                                                 'numpy')

                # should handle local data here...

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

        if node.Type in ('PlaneWave','SampledPlaneWave'):
            for test in standardSourceParameterSet:
                if test(objectparams) != True:
                    raise KeyError, 'readXML.processObject(): need standard parameter(s) %s in source %s' % (test(objectparams),objectname)

        if node.Type == 'PlaneWave':
            try:
                objecttype = objectparams['SourceType'][0]
            except:
                raise KeyError, 'readXML.processObject(): need SourceType for PlaneWave source %s' % objectname
        else:
            objecttype = node.Type

        # look inside included TimeSeries to get more parameters there
        # can only handle the first TimeSeries within each SampledPlaneWave

        # look inside included Array to get more parameters there
        # can only handle the first Array

        if node.Type == 'SampledPlaneWave':
            for node2 in node:
                if node2.tagName == 'XSIL' and node2.Type == 'TimeSeries':
                    for node3 in node2:
                        if node3.tagName == 'Param':
                            objectparams[node3.Name] = self.getParam(node3)                   
                        elif node3.tagName == 'Array':
                            self.processArray(node3,objectparams)
                    break             
    
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
                    if param[2]:
                        thisparam = (param[2],param[1])
                    else:
                        raise AttributeError, 'readXML.processObject(): need parameter(s) %s in object %s of type %s' % (param[0],objectname,objecttype)
            else:
                thisparam = objectparams[param[0]]

            # convert to the correct units (if we know how to)

            thisparam = convertUnit(thisparam[0],thisparam[1],param[1],param[0])

            if param[1] == 'String' or param[1] == 'numpy':
                evalparam = thisparam[0]
            else:
                try:
                    # first try converting to an int...
                    
                    evalparam = int(thisparam[0])
                except ValueError:
                    # if it doesn't work, try a float...
    
                    try:
                        evalparam = float(thisparam[0])
                    except ValueError:
                        # if the float does not work, try calling Python...

                        evalparam = eval(thisparam[0])

            arglist.append(evalparam)

        return (XMLToObject[objecttype][1],arglist,objectname,objectparams)


# for the moment handle only remote files...

def handleFiles(args):
    files = {}
    
    for i in range(0,len(args)):
        if isinstance(args[i],fileData):
            if not args[i].filename in files:
                readlength = 8 * args[i].length * args[i].records

                try:
                    binaryfile = open(args[i].filename,'r')
                    readbuffer = numpy.fromstring(binaryfile.read(readlength),'double')
                    binaryfile.close()
                except:
                    try:
                        binaryfile = open(args[i].basename,'r')
                        readbuffer = numpy.fromstring(binaryfile.read(readlength),'double')
                        binaryfile.close()
                    except:
                        raise IOError, 'handleFiles(): problems reading file %s' % args[i].filename

                if (('BigEndian' in args[i].encoding and sys.byteorder == 'little') or
                    ('LittleEndian' in args[i].encoding and sys.byteorder == 'big')):
                    readbuffer = readbuffer.byteswap()

                if args[i].records == 1:
                    files[args[i].filename] = [readbuffer]
                else:
                    readbuffer = numpy.reshape(readbuffer,[args[i].length,args[i].records])
                
                    files[args[i].filename] = [readbuffer[:,j].copy() for j in range(0,args[i].records)]

            args[i] = files[args[i].filename][args[i].index].copy()

def MakeObject(s):
    # if any data needs to be loaded, do so!

    for arg in s[1]:
        if isinstance(arg,fileData):
            handleFiles(s[1])
            break
    
    ret = (s[0])(*(s[1]))
    ret.name = s[2]
    
    for param in s[3]:
        if s[3][param][1] != 'numpy':
            ret.__setattr__(param,(s[3][param]))
    
    return ret


class LISASourceFactory:
    def __init__(self,sourcelist):
        self.sourcelist = sourcelist
        self.sourcenumber = len(sourcelist)
        
    def __getitem__(self,index):
        return MakeObject(self.sourcelist[index])

    def __getslice__(self,index1,index2):
        return map(MakeObject,self.sourcelist[index1:index2])
    
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

            return MakeObject(i)
