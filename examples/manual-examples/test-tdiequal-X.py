#!/usr/bin/env python

# test of first-generation TDI X noise against the theoretical expression

# this script demonstrates the generation of a TDI X time series, and
# in particular:
# - creating a LISA geometry object
# - creating a TDInoise object based on standard pseudorandom noises
# - calling getobs to get an array of X values at equispaced times
# - getting the spectrum of the time series
# - writing the spectrum to disk
# - using some Python/Numeric magic to create a theoretical X spectrum

# import all the libraries that are needed

from synthlisa import *

# we create a LISA geometry object corresponding to a stationary LISA
# with equal armlengths

L = 16.6782
originallisa = OriginalLISA(L)

# we create a TDInoise object based on the LISA geometry object, and
# on 18 (6+6+6) pseudorandom noise objects with standard parameters
# (the exponents are implicitly -2.0, 2.0, and 0.0; the interpolation
# windows are 1; the sampling time of the noise needs to be the same,
# or lower than the sampling time of the final TDI time series, to
# prevent aliasing)

originalTDI = TDInoise(originallisa,
                       1.0, 2.5e-48, # proof-mass noise parameters
                       1.0, 1.8e-37, # optical-path noise parameters
                       1.0, 1.1e-26) # laser frequency noise parameters

# get "samples" values of the noise, at times separated by "stime"
# note that X is a "method" of the TDInoise object

samples = 2**18 # 2**18 takes 11 s on a 1.25GHz
stime = 1

patches = 256

noiseX = getobsc(samples,stime,originalTDI.Xm)

# the result is a 1D Numeric array, that we feed to "spect" to get a
# triangle-windowed, averaged spectrum, using "patches" averaging
# periods

myspecX = spect(noiseX,stime,patches)

# ...and write the spectrum to disk

writearray('data/tdiequal-X-sampled.txt',myspecX[1:])

# now let's do the theoretical expression. we need pi...

from math import pi

# we get the minimum and maximum frequency from the first column of
# the spectrum computed above to create the range of frequencies where
# we're computing our theoretical spectrum (as an array index, -1 is
# shorthand for the last element)

fmin = myspecX[1,0]
fmax = myspecX[-1,0]

f = numpy.arange(fmin,fmax,(fmax-fmin)/999,'d')

# theoretical expression from Estabrook, Tinto, and Armstrong,
# Phys. Rev. D 62, 042002 (2000); just math follows

Syproof = 2.5e-48 * f**-2
Syopt = 1.8e-37 * f**2

theoryX = ( 8*numpy.sin(4*pi*f*L)**2 + 32*numpy.sin(2*pi*f*L)**2 ) * Syproof + 16*numpy.sin(2*pi*f*L)**2 * Syopt

# we want to write a 2-column file

writearray('data/tdiequal-X-theory.txt',numpy.transpose([f, theoryX]))
