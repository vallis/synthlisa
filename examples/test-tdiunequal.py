#!/usr/bin/python

# test of first-generation TDI noises for a stationary unequal-arm interferometer

# this script demonstrates the generation of a TDI time series, and
# in particular:
# - creating a LISA geometry object
# - creating a TDInoise object based on standard pseudorandom noises
# - defining a short hand for a basic Doppler observable y_slr,d1...dn expression
# - calling getobs to get a 2D array of several TDI observables at equispaced times
# - getting the spectra of the time series
# - writing the spectra to disk

# import all the libraries that are needed

from lisaswig import *
from Numeric import *
from lisautils import *

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

tdiy231 = lambda t : originalTDI.y(2,3,1,0,0,0,0,0,0,0,t)

[noiseX,
 noisea,
 noisez,
 noiseu,
 noisey] = transpose(getobs(samples,stime,[originalTDI.X,
                                           originalTDI.alpha,
                                           originalTDI.zeta,
                                           originalTDI.U,
                                           tdiy231]))

# compute spectra, and write to disk

myspecX = spect(noiseX,stime,patches)
writearray('data/tdiunequal-X-freq.txt',    myspecX[1:])

myspeca = spect(noisea,stime,patches)
writearray('data/tdiunequal-alpha-freq.txt',myspeca[1:])

myspecz = spect(noisez,stime,patches)
writearray('data/tdiunequal-zeta-freq.txt', myspecz[1:])

myspecu = spect(noiseu,stime,patches)
writearray('data/tdiunequal-U-freq.txt', myspecu[1:])

myspecy = spect(noisey,stime,patches)
writearray('data/tdiunequal-y231-freq.txt', myspecy[1:])


