grep -Rn "print" . |grep .py
-> only round 55 places to edit !
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
./lisasim/lisapar.py:94:                print "LISApar.getobsp(...): third parameter must be a 4-tuple containing a",
./lisasim/lisapar.py:95:                print "LISA instance, a Wave factory, an array of parameters for the factory,",
./lisasim/lisapar.py:96:                print "and a set of TDI observables given as class methods (such as synthlisa.TDI.X)."
./lisasim/lisapar.py:101:                print "LISApar.getobsp(...): needs a list of parameters to feed to the factory!"
./lisasim/lisapar.py:106:                print "LISApar.getobsp(...): must be run with more than one cpu!"
./lisasim/lisapar.py:111:                print "LISApar.getobsp(...): needs to run with more sources than cpus!"
./lisasim/lisapar.py:123:            print "Standard block: ", blocksize,
./lisasim/lisapar.py:124:            print "; root block: ", len(parameters) - blocksize * (size-1)
./lisasim/lisapar.py:128:                print "Preparing for parallel execution..."
./lisasim/lisapar.py:147:            print "CPU ", myrank, " received ", len(mypars), " source parameters ", mypars
./lisasim/lisapar.py:159:                print "LISApar.getobsp(...): srcfunc must return a synthlisa.Wave when applied",
./lisasim/lisapar.py:160:                print "to each element of the parameter list"
./lisasim/lisapar.py:164:            print "CPU ", myrank, " created sources ", sources
./lisasim/lisapar.py:170:                print "LISApar.getobsp(...): lisa must be an instance of synthlisa.LISA."
./lisasim/lisapar.py:199:            print "Completed in %d s [%d (multi)samples/s]." % (int(currenttime),int(vel))
./lisasim/lisautils.py:284:        print "\r...%d/%d (%d%%) done [%d (multi)samples/s], est. %dm%ds left..." % (i,snum,percdone,vel,minleft,secleft),
./lisasim/lisautils.py:298:    print "Processing...",
./lisasim/lisautils.py:317:        print "lisautils::getobsc: I have trouble accessing time ", zerotime+i*stime,
./lisasim/lisautils.py:318:        print "; you may try to reset your objects and repeat..."
./lisasim/lisautils.py:333:        print "\r...completed in %d s [%d (multi)samples/s].                           " % (int(currenttime),int(vel))
./lisasim/lisautils.py:335:        print "\r...completed.                                                         "
./lisasim/lisaxml.py:293:    def iprint(self,s):
./lisasim/lisaxml.py:295:            print >> self.f, self.indent + s
./lisasim/lisaxml.py:317:        self.iprint(string)
./lisasim/lisaxml.py:325:        self.iprint(string)
./lisasim/lisaxml.py:334:        self.iprint(string)
./lisasim/lisaxml.py:342:        self.iprint(string)
./lisasim/lisaxml.py:352:            self.iprint(os.path.basename(filename))
./lisasim/lisaxml.py:360:                self.iprint(line)
./lisasim/lisaxml.py:728:            print "lisaXML::TDIData: data must be a proper Numeric array of doubles"
./lisasim/lisaxml.py:736:            print "lisaXML::TDIData: behavior undetermined for arrays of this dimension"
./lisasim/lisaxml.py:742:            print "lisaXML::TDIData: data shorter than requested output length"
./lisasim/lisaxml.py:774:                print 'lisaXML::writedata: cannot determine binary encoding'
./lisasim/lisaxml.py:807:            print 'lisaXML::writedata: encoding not implemented'
./lisasim/lisaxml.py:897:            print "lisaXML::TDISpectra: data must be a proper Numeric array of doubles"
./lisasim/lisaxml.py:905:            print "lisaXML::TDISpectra: behavior undetermined for arrays of this dimension"
./lisasim/lisaxml.py:911:            print "lisaXML::TDISpectra: data shorter than requested output length"
./lisasim/lisaxml.py:963:            print "lisaXML::close(): File is closed already"
./lisasim/lisaxml.py:1051:            print "XML validation error! (Or perhaps I couldn't access the DTD)."
./lisasim/lisaxml.py:1052:            print "I'll try to use the file anyway by removing the DTD..."
./lisasim/lisaxml.py:1058:            print 'Not a LISA XSIL file!'


x./setup.py:91:print >> version_py, "version_full = \"\"\"%s\"\"\"\n" % idcatalog
x./setup.py:92:print >> version_py, "version_short = \"%s\"\n" % versiontag
x./setup.py:115:            print 'Sorry, I am unable to swig the modified ' + lisasim_isource
x./setup.py:157:print >> setdir_sh, """if [ -z "${PYTHONPATH}" ]
x./setup.py:165:print >> setdir_csh, """if !($?PYTHONPATH) then
x./setup.py:172:print >> recompile_sh, """#!/bin/sh
x./setup.py:271:                print "No GSL, skipping " + contrib_swigfile
