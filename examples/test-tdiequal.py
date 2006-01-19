#!/usr/bin/python

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
# - calling getobs to get a 2D array of several TDI observables at equispaced times
# - getting the spectra of the time series
# - writing the spectra to disk

# import all the libraries that are needed

from synthlisa import *
from Numeric import *

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

samples = 2**16    # 2**18 takes 60 s on a 1.25GHz G4; 2**21 used for plot

# number of averaging patches for the spectra

patches = 256

# getting the value of one of the basic Doppler observables requires a little trick:
# here send = 2, link = 3, recv = 1, and all the seven possible delays are set to 0
# with a similar syntax (but eight delays) you could also get the z's

originalTDI.y231 = lambda t : originalTDI.y(2,3,1,0,0,0,0,0,0,0,t)

# notice the call to getobs with a list "[...]" of observables
# the result is a 2D array with four columns, corresponding to the four observables
 
observables = getobsc(samples,stime,[originalTDI.X,
                                     originalTDI.alpha,
                                     originalTDI.zeta,
                                     originalTDI.U,
                                     originalTDI.y231])

outputXML = lisaXML('data/tdiequal',
                    author='Michele Vallisneri',
                    comments='OriginalLISA (stationary equal-arm) TDI noises')

# save time series to disk

outputXML.TDIData(observables,samples,stime,'Xf,alphaf,zetaf,Uf,y231f')

# collect the columns of 'observables' in single-column arrays
# for each observable

[noiseX,noisea,noisez,noiseu,noisey] = transpose(observables)

# compute the spectra and write them to disk

myspecX = spect(noiseX,stime,patches)
myspeca = spect(noisea,stime,patches)
myspecz = spect(noisez,stime,patches)
myspecu = spect(noiseu,stime,patches)
myspecy = spect(noisey,stime,patches)

outputXML.TDISpectraSelfDescribed(myspecX[1:],'f,Xf')
outputXML.TDISpectraSelfDescribed(myspeca[1:],'f,alphaf')
outputXML.TDISpectraSelfDescribed(myspecz[1:],'f,zetaf')
outputXML.TDISpectraSelfDescribed(myspecu[1:],'f,uf')
outputXML.TDISpectraSelfDescribed(myspecy[1:],'f,y231f')

# all data is written only on closing the outputXML; thus the arrays referenced
# in the TDIData and TDISpectra calls should not be altered in the meantime

outputXML.close()

# if PyX is available, plot the spectra

try:
    import pyx

    # ??? how to choose the plotting range?

    g1 = pyx.graph.graphxy(width=8,
                           x=pyx.graph.axis.log(),
                           y=pyx.graph.axis.log())

    g1.plot(pyx.graph.data.list(map(tuple,myspecX[1:]),x=1,y=2),
            [pyx.graph.style.line([pyx.color.rgb.black])])
    g1.plot(pyx.graph.data.list(map(tuple,myspeca[1:]),x=1,y=2),
            [pyx.graph.style.line([pyx.color.rgb.red])])
    g1.plot(pyx.graph.data.list(map(tuple,myspecz[1:]),x=1,y=2),
            [pyx.graph.style.line([pyx.color.rgb.green])])
    g1.plot(pyx.graph.data.list(map(tuple,myspecu[1:]),x=1,y=2),
            [pyx.graph.style.line([pyx.color.rgb.blue])])

    g1.writeEPSfile("eps/tdiequal-spectra")

    print "Spectra plotted in eps/tdiequal-spectra.eps"

    # ??? encapsulate in lisautils.py?
    # ??? legend?

except ImportError:
    print "Cannot find the plotting module PyX; skipping plotting"
