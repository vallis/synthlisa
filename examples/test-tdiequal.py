#!/usr/bin/python

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

from lisaswig import *
from Numeric import *
from lisautils import *

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

samples = 2**18    # 2**18 takes 60 s on a 1.25GHz G4; 2**21 used for plot

patches = 256

# getting the value of one of the basic Doppler observables requires a little trick:
# here send = 2, link = 3, recv = 1, and all the seven possible delays are set to 0
# with a similar syntax (but eight delays) you could also get the z's

tdiy231 = lambda t : originalTDI.y(2,3,1,0,0,0,0,0,0,0,t)

# notice the call to getobs with a list "[...]" of observables
# the result is a 2D array with four columns, corresponding to the four observables
# we use transpose to get the four columns into four separate 1D arrays

[noiseX,
 noisea,
 noisez,
 noiseu,
 noisey] = transpose(getobs(samples,stime,[originalTDI.X,
                                           originalTDI.alpha,
                                           originalTDI.zeta,
                                           originalTDI.U,
                                           tdiy231]))

# you could also use the observable "time" to get the time at which the others
# are evaluated

# compute the spectra and write them to disk

myspecX = spect(noiseX,stime,patches)
writearray('data/tdiequal-X-freq.txt',myspecX[1:])

myspeca = spect(noisea,stime,patches)
writearray('data/tdiequal-alpha-freq.txt',myspeca[1:])

myspecz = spect(noisez,stime,patches)
writearray('data/tdiequal-zeta-freq.txt',myspecz[1:])

myspecu = spect(noiseu,stime,patches)
writearray('data/tdiequal-U-freq.txt',myspecu[1:])

myspecy = spect(noisey,stime,patches)
writearray('data/tdiequal-y231-freq.txt',myspecy[1:])
