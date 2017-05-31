grep -Rn "print" . |grep .py -->
--------------------------------

./examples/paper-examples/plotxml.py:45:            print "multiplot: need more output files!"
./examples/paper-examples/test-binary.py:88:print "S/N = ", math.sqrt(sn(signalA,noiseA,stime,patches)**2 +
./examples/paper-examples/test-binary.py:119:    print "Spectra plotted in eps/binary.pdf"
./examples/paper-examples/test-binary.py:121:    print "Could not output PyX graphics"
./examples/paper-examples/test-tdiequal.py:65:# getobsc (as opposed to getobs) prints status information during evaluation
./examples/paper-examples/test-tdiequal.py:133:print "You can plot the results of this script by running"
./examples/paper-examples/test-tdiequal.py:134:print "  ./plotxml.py data/tdiequal.xml eps/tdiequal-spectra.eps"
./examples/paper-examples/test-tdiequal.py:135:print "or"
./examples/paper-examples/test-tdiequal.py:136:print "  ./plotxml.py data/tdiequal.xml eps/tdiequal-spectra.pdf"

x./lisasim/lisapar.py:94:                print "LISApar.getobsp(...): third parameter must be a 4-tuple containing a",
x./lisasim/lisapar.py:95:                print "LISA instance, a Wave factory, an array of parameters for the factory,",
x./lisasim/lisapar.py:96:                print "and a set of TDI observables given as class methods (such as synthlisa.TDI.X)."
x./lisasim/lisapar.py:101:                print "LISApar.getobsp(...): needs a list of parameters to feed to the factory!"
x./lisasim/lisapar.py:106:                print "LISApar.getobsp(...): must be run with more than one cpu!"
x./lisasim/lisapar.py:111:                print "LISApar.getobsp(...): needs to run with more sources than cpus!"
x./lisasim/lisapar.py:123:            print "Standard block: ", blocksize,
x./lisasim/lisapar.py:124:            print "; root block: ", len(parameters) - blocksize * (size-1)
x./lisasim/lisapar.py:128:                print "Preparing for parallel execution..."
x./lisasim/lisapar.py:147:            print "CPU ", myrank, " received ", len(mypars), " source parameters ", mypars
x./lisasim/lisapar.py:159:                print "LISApar.getobsp(...): srcfunc must return a synthlisa.Wave when applied",
x./lisasim/lisapar.py:160:                print "to each element of the parameter list"
x./lisasim/lisapar.py:164:            print "CPU ", myrank, " created sources ", sources
x./lisasim/lisapar.py:170:                print "LISApar.getobsp(...): lisa must be an instance of synthlisa.LISA."
x./lisasim/lisapar.py:199:            print "Completed in %d s [%d (multi)samples/s]." % (int(currenttime),int(vel))

x./lisasim/lisautils.py:284:        print "\r...%d/%d (%d%%) done [%d (multi)samples/s], est. %dm%ds left..." % (i,snum,percdone,vel,minleft,secleft),
x./lisasim/lisautils.py:298:    print "Processing...",
x./lisasim/lisautils.py:317:        print "lisautils::getobsc: I have trouble accessing time ", zerotime+i*stime,
x./lisasim/lisautils.py:318:        print "; you may try to reset your objects and repeat..."
x./lisasim/lisautils.py:333:        print "\r...completed in %d s [%d (multi)samples/s].                           " % (int(currenttime),int(vel))
x./lisasim/lisautils.py:335:        print "\r...completed.                                                         "

x./lisasim/lisaxml.py:728:            print "lisaXML::TDIData: data must be a proper Numeric array of doubles"
x./lisasim/lisaxml.py:736:            print "lisaXML::TDIData: behavior undetermined for arrays of this dimension"
x./lisasim/lisaxml.py:742:            print "lisaXML::TDIData: data shorter than requested output length"
x./lisasim/lisaxml.py:774:                print 'lisaXML::writedata: cannot determine binary encoding'
x./lisasim/lisaxml.py:807:            print 'lisaXML::writedata: encoding not implemented'
x./lisasim/lisaxml.py:897:            print "lisaXML::TDISpectra: data must be a proper Numeric array of doubles"
x./lisasim/lisaxml.py:905:            print "lisaXML::TDISpectra: behavior undetermined for arrays of this dimension"
x./lisasim/lisaxml.py:911:            print "lisaXML::TDISpectra: data shorter than requested output length"
x./lisasim/lisaxml.py:963:            print "lisaXML::close(): File is closed already"
x./lisasim/lisaxml.py:1051:            print "XML validation error! (Or perhaps I couldn't access the DTD)."
x./lisasim/lisaxml.py:1052:            print "I'll try to use the file anyway by removing the DTD..."
x./lisasim/lisaxml.py:1058:            print 'Not a LISA XSIL file!'
x all the raise xyz...


x./setup.py:91:print >> version_py, "version_full = \"\"\"%s\"\"\"\n" % idcatalog
x./setup.py:92:print >> version_py, "version_short = \"%s\"\n" % versiontag
x./setup.py:115:            print 'Sorry, I am unable to swig the modified ' + lisasim_isource
x./setup.py:157:print >> setdir_sh, """if [ -z "${PYTHONPATH}" ]
x./setup.py:165:print >> setdir_csh, """if !($?PYTHONPATH) then
x./setup.py:172:print >> recompile_sh, """#!/bin/sh
x./setup.py:271:                print "No GSL, skipping " + contrib_swigfile
"

grep -Rn "raise" . |grep .py -->
--------------------------------

x./lisasim/convertunit.py:56:    raise NotImplementedError("convertUnit(): cannot convert {} from {} to {} (parameter {})".format(param,unitin,unitout,paramname))
x./lisasim/convertunit.py:117:        raise NotImplementedError("convertParameters(): LISA eta/xi configuration with sw > 0 not compatible with PseudoLISA")
x./lisasim/convertunit.py:161:        raise NotImplementedError, "convertParameters(): unknown interpolator type %s" % interptype
x./lisasim/convertunit.py:179:        raise NotImplementedError, "convertParameters(): unknown interpolator length %s" % interplen

x./lisasim/lisaxml.py:482:            raise KeyError, 'lisaXML.ObjectData: unknown object type %s' % objecttype
x./lisasim/lisaxml.py:512:                    raise AttributeError, 'lisaXML.ProcessObjectData(): missing internal parameter %s in object %s' % (param[0],object)
x./lisasim/lisaxml.py:538:                        raise AttributeError, 'readXML.ProcessObjectData(): missing external parameter(s) %s for object %s' % (param[0],object)
x./lisasim/lisaxml.py:584:                        raise AttributeError, 'readXML.doSampledWaveTimeSeries(): missing external parameter(s) %s for object %s' % (param[0],object)
x./lisasim/lisaxml.py:1285:            raise AttributeError, 'readXML.getTDInoise(): problems reading LISA geometry'
x./lisasim/lisaxml.py:1290:            raise AttributeError, 'readXML.getTDInoise(): problems reading LISA noises'
x./lisasim/lisaxml.py:1323:                    raise AttributeError, 'readXML.processArray(): encoding/type not specified in stream'
x./lisasim/lisaxml.py:1352:                    raise KeyError, 'readXML.processObject(): need standard parameter(s) %s in source %s' % (test(objectparams),objectname)
x./lisasim/lisaxml.py:1358:                raise KeyError, 'readXML.processObject(): need SourceType for PlaneWave source %s' % objectname
x./lisasim/lisaxml.py:1379:            raise NotImplementedError, 'readXML.processObject(): unknown object type %s for object %s' % (objecttype,objectname)
x./lisasim/lisaxml.py:1385:                raise KeyError, 'readXML.processObject(): need parameter(s) %s for object %s of type %s' % (test(objectparams),objectname,objecttype)
x./lisasim/lisaxml.py:1404:            raise NotImplementedError, 'readXML.processObject(): unknown object %s of type %s' % (objectname,objecttype)
x./lisasim/lisaxml.py:1424:                        raise AttributeError, 'readXML.processObject(): need parameter(s) %s in object %s of type %s' % (param[0],objectname,objecttype)
x./lisasim/lisaxml.py:1474:                        raise IOError, 'handleFiles(): problems reading file %s' % args[i].filename
x./lisasim/lisaxml.py:1536:            raise StopIteration

x./lisasim/xmlutils.py:99:            raise AttributeError, msg
x./lisasim/xmlutils.py:116:            raise IndexError, '%s no index %s' % (self.__repr__(), `idx`)

x./lisasim/lisasim-swig.i:565:        raise NotImplementedError, "getInterpolator: undefined interpolator length %s (lisasim-swig.i)." % interplen
x./lisasim/lisasim-swig.i:573:        raise NotImplementedError, "getDerivativeInterpolator: undefined interpolator length %s (lisasim-swig.i)." % interplen
x./lisasim/lisasim-swig.i:593:        raise NotImplementedError, "PowerLawNoise: undefined PowerLaw exponent %s (lisasim-swig.i)." % exponent
x./lisasim/lisasim-swig.i:632:        raise NotImplementedError, "SampledSignal: need Numeric array or filename as first argument (lisasim-swig.i)."
