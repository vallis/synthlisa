# $Id$
# $Date$
# $Author$
# $Revision$

import sys
import os.path
import time
import string
import re

import Numeric

import pyRXP
# require xmlutils from pyRXP examples
import xmlutils

class writeXML:
    def __init__(self,filename):
        self.f = open(filename,'w')
        self.opened = 1

    def __del__(self):
        if self.opened == 1:
            self.close()

    def close(self):
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
        if self.opened:
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

class typeLISATimeSeries:
    pass

class typeLISAFrequencySeries:
    pass

class lisaXML(writeXML):
    """Create lisaXML file with metadata (author,date,comments);
       date should be given as ISO-8601."""
    # actual writing of file
    
    def __init__(self,filename,author='',date='',comments=''):
        """Create lisaXML file with metadata [author,date,comments];
           date should be given as ISO-8601."""
    
        self.filename = filename
    
        # add ".xml" to file if necessary
        if not re.match('.*\.xml$',self.filename):
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
        self.theTDIData = []

    def TDIData(self,data,length,cadence,description,offset=0,encoding='Binary',comments=''):
        """Add a LISATimeSeries object to a lisaXML file object. Here
        'data' is the Numeric array containing the time series
        (simultaneous entries on the same row); 'length' is the desired
        length to be written in the array; 'cadence' is its nominal
        cadence in seconds; 'description' should be a comma-separated
        string listing the TDI observables represented in the time
        series; 'offset' is the nominal initial time for the data
        (currently in seconds); 'encoding' can be 'Binary' for storage
        in a separate binary file, or 'Text' for inline storage in the
        XML file; 'comments' is added to the LISATimeSeries entry. To
        skip some records at the beginning of the array, use a slicing
        syntax such as data[1:].
        
        Each lisaXML file can contain several LISATimeSeries objects,
        all contained in the TDIData block; if binary storage is
        requested, each LISATimeSeries (or LISAFrequencySeries) object
        corresponds to a separate binary file, named file-0.bin,
        file-1.bin, etc., if the main XML file is file.xml."""
           
        # mimick the getobsc call, and get description as comma-separated string
        # ??? provide way to read variable names directly from list
        # Should use different names for phase and frequency TDI variables

        LISATimeSeries = typeLISATimeSeries()
        
        try:
            LISATimeSeries.data = data
            LISATimeSeries.dim = len(Numeric.shape(data))
            LISATimeSeries.alength = Numeric.shape(data)[0]
            
            if data.typecode() != 'd':
                raise TypeError
        except:
            print "lisaXML::TDIData: data must be a proper Numeric array of doubles"
            raise TypeError
        
        if LISATimeSeries.dim == 1:
            LISATimeSeries.records = 1
        elif LISATimeSeries.dim == 2:        
            LISATimeSeries.records = Numeric.shape(data)[1]
        else:
            print "lisaXML::TDIData: behavior undetermined for arrays of this dimension"
            raise NotImplementedError

        LISATimeSeries.length = length

        if length > LISATimeSeries.alength:
            print "lisaXML::TDIData: data shorter than requested output length"
            raise IndexError

        LISATimeSeries.cadence = cadence

        LISATimeSeries.duration = cadence * length

        LISATimeSeries.description = description
        LISATimeSeries.start = offset
        
        LISATimeSeries.encoding = encoding
   
        LISATimeSeries.comments = comments
   
        self.theTDIData.append(LISATimeSeries)

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

    def writeTimeSeries(self,LISATimeSeries):
        # write out the LISATimeSeries defined in LISATimeSeries
        
        self.opentag('XSIL',{'Type': 'LISATimeSeries',
                             'Name': LISATimeSeries.description})

        if LISATimeSeries.comments:
            self.opentag('Comment',{})
            self.content(LISATimeSeries.comments)
            self.closetag('Comment')
        
        # ??? fix the time types to "s" (XSIL extension) for the moment
        
        self.coupletag('Time',{'Name': 'StartTime','Type': 's'},
                              str(LISATimeSeries.start))

        self.coupletag('Param',{'Name': 'Cadence','Unit': 's'},
                               str(LISATimeSeries.cadence))        

        self.coupletag('Param',{'Name': 'Duration','Unit': 's'},
                                  str(LISATimeSeries.duration))
       
        # ??? use <Column> to define columns (not in XSIL, but in BFD)?

        self.writearray(LISATimeSeries.data,
                        LISATimeSeries.length,
                        LISATimeSeries.records,
                        LISATimeSeries.description,
                        LISATimeSeries.encoding)
       
        self.closetag('XSIL')

    def TDISpectraSelfDescribed(self,data,description,encoding='Binary',comments=''):
        """Add a LISAFrequencySeries object to a lisaXML file object.
        Here 'data' is the Numeric array containing the time
        series (simultaneous entries on the same row); 'description'
        is a comma-separated string listing the TDI observables
        represented in the time series; 'encoding' can be 'Binary'
        for storage in a separate binary file, or 'Text' for inline
        storage in the XML file; 'comments' is added to the
        LISAFrequencySeries entry.
        
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
        """Add a LISAFrequencySeries object to a lisaXML file object.
        Here 'data' is the Numeric array containing the time series
        (simultaneous entries on the same row); 'length' is the desired
        length of the array to be written; 'deltaf' is the difference
        between successive records [in Hz]; comma-separated string
        listing the TDI observables represented in the time series;
        'offset' is the initial frequency for the spectra (in Hz);
        'encoding' can be 'Binary' for storage in a separate binary
        file, or 'Text' for inline storage in the XML file; 'comments'
        is added to the LISAFrequencySeries entry. To skip some records
        at the beginning of the array, use a slicing syntax such as
        data[1:].
        
        Each lisaXML file can contain several LISAFrequencySeries
        objects, all contained in the TDIData block; if binary storage
        is requested, each LISAFrequencySeries (or LISATimeSeries)
        object corresponds to a separate binary file, named file-0.bin,
        file-1.bin, etc., if the main XML file is file.xml."""

        # mimick the getobsc call, and get description as comma-separated string
        # ??? provide way to read variable names directly from list
        # Should use different names for phase and frequency TDI variables

        LISAFrequencySeries = typeLISAFrequencySeries()
        
        try:
            LISAFrequencySeries.data = data
            LISAFrequencySeries.dim = len(Numeric.shape(data))
            LISAFrequencySeries.alength = Numeric.shape(data)[0]
            
            if data.typecode() != 'd':
                raise TypeError
        except:
            print "lisaXML::TDISpectra: data must be a proper Numeric array of doubles"
            raise TypeError
        
        if LISAFrequencySeries.dim == 1:
            LISAFrequencySeries.records = 1
        elif LISAFrequencySeries.dim == 2:        
            LISAFrequencySeries.records = Numeric.shape(data)[1]
        else:
            print "lisaXML::TDISpectra: behavior undetermined for arrays of this dimension"
            raise NotImplementedError

        LISAFrequencySeries.length = length

        if length > LISAFrequencySeries.alength:
            print "lisaXML::TDISpectra: data shorter than requested output length"
            raise IndexError

        LISAFrequencySeries.minf = offset * deltaf
        LISAFrequencySeries.maxf = (offset + length - 1) * deltaf
        LISAFrequencySeries.deltaf = deltaf

        LISAFrequencySeries.description = description
        
        LISAFrequencySeries.encoding = encoding
   
        LISAFrequencySeries.comments = comments
   
        self.theTDIData.append(LISAFrequencySeries)

    def writeFrequencySeries(self,LISAFrequencySeries):
        # write out the LISAFrequencySeries defined in LISAFrequencySeries
        
        self.opentag('XSIL',{'Type': 'LISAFrequencySeries',
                             'Name': LISAFrequencySeries.description})

        if LISAFrequencySeries.comments:
            self.opentag('Comment',{})
            self.content(LISAFrequencySeries.comments)
            self.closetag('Comment')
        
        # ??? fix the frequency types to "Hz" (XSIL extension) for the moment
        # ??? provide facility for automatic step specification
        
        self.coupletag('Param',{'Name': 'MinFreq','Type': 'Hz'},
                                      str(LISAFrequencySeries.minf))
        self.coupletag('Param',{'Name': 'MaxFreq','Type': 'Hz'},
                                      str(LISAFrequencySeries.maxf))
        self.coupletag('Param',{'Name': 'DeltaFreq','Type': 'Hz'},
                                      str(LISAFrequencySeries.deltaf))
       
        # ??? use <Column> to define columns (not in XSIL, but in BFD)?
       
        self.writearray(LISAFrequencySeries.data,
                        LISAFrequencySeries.length,
                        LISAFrequencySeries.records,
                        LISAFrequencySeries.description,
                        LISAFrequencySeries.encoding)
       
        self.closetag('XSIL')
    
    def close(self):
        """Write the XML file to disk. This happens also on destruction of
        the lisaXML object."""
        
        if self.opened == 0:
            print "lisaXML::close(): File is closed already"
            raise IOError
        
        self.content('<?xml version="1.0"?>')
        self.content('<!DOCTYPE XSIL SYSTEM "bfd.dtd">')

        self.opentag('XSIL',{})

        self.coupletag('Param',{'Name': 'Author'},self.author)
        self.coupletag('Param',{'Name': 'GenerationDate', 'Type': 'ISO-8601'},
                               self.date)
        
        if self.comments:
            self.opentag('Comment',{})
            self.content(self.comments)
            self.closetag('Comment')

        # ??? do the LISAModel (LISAGeometry and LISANoise)

        # ??? do the SourceModel objects

        # do the TDIdata objects (first supported)
        
        if len(self.theTDIData) > 0:
            self.opentag('XSIL',{'Type': 'TDIData','Name': 'TDIData'})
    
            for object in self.theTDIData:
                if isinstance(object,typeLISATimeSeries):
                    self.writeTimeSeries(object)
                elif isinstance(object,typeLISAFrequencySeries):
                    self.writeFrequencySeries(object)

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
            return (float(str(node)),node.Unit)
        except AttributeError:
            return (float(str(node)),)     

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
                                raise NotImplemented
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
                                raise NotImplemented
    
        return timeseries

    def getLISATimeSeries(self):
        result = []
        
        for node in self.tw:
            # outermost XSIL level container
            
            if node.tagName == 'XSIL':
                if node.Name == 'TDIData':
                    # inside TDIData
        
                    for node2 in node:
                        if node2.tagName == 'XSIL':
                            if node2.Type == 'LISATimeSeries':
                                # got a TimeSeries!
        
                                result.append(self.processSeries(node2))
    
        return result

    def getLISAFrequencySeries(self):
        result = []
        
        for node in self.tw:
            # outermost XSIL level container
            
            if node.tagName == 'XSIL':
                if node.Name == 'TDIData':
                    # inside TDIData
        
                    for node2 in node:
                        if node2.tagName == 'XSIL':
                            if node2.Type == 'LISAFrequencySeries':
                                # got a FrequencySeries
                                
                                result.append(self.processSeries(node2))
    
        return result
