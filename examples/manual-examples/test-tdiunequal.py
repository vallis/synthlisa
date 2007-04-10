#!/usr/bin/env python

# $Id$
# $Date$
# $Author$
# $Revision$

# test of first-generation TDI noises for a stationary unequal-arm interferometer

# this script demonstrates the generation of a TDI time series, and
# in particular:
# - creating a LISA geometry object
# - creating a TDInoise object based on standard pseudorandom noises
# - defining a short hand for a basic Doppler observable y_slr,d1...dn expression
# - calling getobs to get a 2D array of several TDI observables at equispaced times
# - getting the spectra of the time series
# - writing the spectra to an XML file

# import all the libraries that are needed

from synthlisa import *
from numpy.oldnumeric import transpose

# we create the stationary LISA (with unequal arms) and the standard
# TDI noise objects

originallisa = OriginalLISA(15.0,16.0,17.0)

stime = 4.0

originalTDI = TDInoise(originallisa,
                       stime, 2.5e-48, # proof-mass noise parameters
                       stime, 1.8e-37, # optical-path noise parameters
                       stime, 1.1e-26) # laser frequency noise parameters

# we get all the observables (see test-tdiequal.py)

samples = 2**18  # 2**18 takes 60 s on a 1.25GHz G4; 2**21 used for plot

patches = 256

originalTDI.y231 = lambda t : originalTDI.y(2,3,1,0,0,0,0,0,0,0,t)

[noiseX,
 noisea,
 noisez,
 noiseu,
 noisey] = transpose(getobsc(samples,stime,[originalTDI.X,
                                            originalTDI.alpha,
                                            originalTDI.zeta,
                                            originalTDI.U,
                                            originalTDI.y231]))

# compute spectra

myspecX = spect(noiseX,stime,patches)
myspeca = spect(noisea,stime,patches)
myspecz = spect(noisez,stime,patches)
myspecu = spect(noiseu,stime,patches)
myspecy = spect(noisey,stime,patches)

# create the output XML file and write spectra (see test-tdiequal.py)

outputXML = lisaXML('data/tdiunequal',
                    author='Michele Vallisneri',
                    comments='OriginalLISA (stationary unequal-arm) TDI noises')

allspec = transpose([myspecX[:,0],myspecX[:,1],
                                  myspeca[:,1],
                                  myspecz[:,1],
                                  myspecu[:,1],
                                  myspecy[:,1]])

outputXML.TDISpectraSelfDescribed(allspec[1:],'f,Xf,alphaf,zetaf,Uf,y231f',encoding='Text')

outputXML.close()

print "You can plot the results of this script by running"
print "  ./plotxml.py data/tdiunequal.xml eps/tdiunequal-spectra.eps"
print "or"
print "  ./plotxml.py data/tdiunequal.xml eps/tdiunequal-spectra.pdf"
