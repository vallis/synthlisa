#!/usr/bin/env python

# $Id$
# $Date$
# $Author$
# $Revision$

# test of first-generation TDI noises for an equal-arm, stationary interferometer

# this script demonstrates the generation of a TDI time series, and
# in particular:
# - creating a LISA geometry object
# - creating a TDInoise object based on standard pseudorandom noises
# - defining a short hand for a basic Doppler observable y_slr,d1...dn expression
# - calling getobsc to get a 2D array of several TDI observables at equispaced times
# - writing the time series to an XML file
# - getting the spectra of the time series
# - writing the spectra to an XML file

# import all the libraries that are needed

from synthlisa import *
from numpy import transpose

# we create a LISA geometry object corresponding to a stationary LISA
# with equal armlengths

originallisa = OriginalLISA(16.6782,16.6782,16.6782)

# we create a TDInoise object based on the LISA geometry object, and
# on 18 (6+6+6) pseudorandom noise objects with standard parameters
# (the exponents are implicitly -2.0, 2.0, and 0.0; the interpolation
# windows are 1)

# need to use the same sampling time (or lower) for noise generation
# w.r.t. the later TDI sampling

stime = 4.0

originalTDI = TDInoise(originallisa,
                       stime, 2.5e-48, # proof-mass noise parameters
                       stime, 1.8e-37, # optical-path noise parameters
                       stime, 1.1e-26) # laser frequency noise parameters

# get "samples" values of the noises, at times separated by "stime"

samples = 2**18   # 2**18 takes 60 s on a 1.25GHz G4; 2**21 used for plot

# number of averaging patches for the spectra

patches = 64

# getting the value of one of the basic Doppler observables requires a
# little trick: here send = 2, link = 3, recv = 1, and all the seven
# possible delays are set to 0 with a similar syntax (but eight delays)
# you could also get the z's

originalTDI.y231 = lambda t : originalTDI.y(2,3,1,0,0,0,0,0,0,0,t)

# notice the call to getobs with a list "[...]" of observables
# (originalTDI.time, or originalTDI.t is just the time); the result is
# a 2D array with six columns, corresponding to time and the five
# observables

# getobsc (as opposed to getobs) prints status information during evaluation

observables = getobsc(samples,stime,[originalTDI.time,
                                     originalTDI.Xm,
                                     originalTDI.Ym,
                                     originalTDI.Zm])
 
# observables = getobsc(samples,stime,[originalTDI.time,
#                                     originalTDI.Xm,
#                                     originalTDI.alpha,
#                                     originalTDI.zeta,
#                                     originalTDI.U,
#                                     originalTDI.y231])

# create the output XML file

outputXML = lisaXML('data/tdiequal',
                    author='Michele Vallisneri',
                    comments='OriginalLISA (stationary equal-arm) TDI noises')

# save time series to the XML file

# outputXML.TDIData(observables,samples,stime,'t,Xmf,Ymf,Zmf')
# outputXML.TDIData(observables,samples,stime,'t,Xf,alphaf,zetaf,Uf,y231f')

# collect the columns of 'observables' in single-column arrays
# for each observable

[time,noiseX,noiseY,noiseZ] = transpose(observables)

# [time,noiseX,noisea,noisez,noiseu,noisey] = transpose(observables)

# compute the spectra and save them to the XML file

myspecX = spect(noiseX,stime,patches)
myspecY = spect(noiseY,stime,patches)
myspecZ = spect(noiseZ,stime,patches)

# myspeca = spect(noisea,stime,patches)
# myspecz = spect(noisez,stime,patches)
# myspecu = spect(noiseu,stime,patches)
# myspecy = spect(noisey,stime,patches)

allspec = transpose([myspecX[:,0],myspecX[:,1],
                                  myspecY[:,1],
                                  myspecZ[:,1]])

#                                  myspeca[:,1],
#                                  myspecz[:,1],
#                                  myspecu[:,1],
#                                  myspecy[:,1]])

# the first column of the spectra contains the frequency abscissa, so
# information about minimum and maximum frequency and about frequency
# cadence may be gotten there; we omit the DC component, which may be
# problematic to plot log-log

outputXML.TDISpectraSelfDescribed(allspec[1:],'f,Xmf,Ymf,Zmf',encoding='Text')
outputXML.TDISpectraSelfDescribed(allspec[1:],'f,Xmf2,Ymf2,Zmf2',encoding='Text')

# outputXML.TDISpectraSelfDescribed(allspec[1:],'f,Xf,alphaf,zetaf,Uf,y231f',encoding='Text')

# all data is written only on closing the outputXML; thus the arrays
# referenced in the TDIData and TDISpectra calls should not be altered
# in the meantime

outputXML.close()

print "You can plot the results of this script by running"
print "  ./plotxml.py data/tdiequal.xml eps/tdiequal-spectra.eps"
print "or"
print "  ./plotxml.py data/tdiequal.xml eps/tdiequal-spectra.pdf"
