#!/usr/bin/env python

# test of first-generation TDI X noise with the baseline configuration,
# and with a bad proof-mass noise

# this script demonstrates the generation of a TDI X time series with
# custom noise objects, and in particular:
# - creating a LISA geometry object
# - creating a TDInoise object based on standard pseudorandom noises
# - creating a TDInoise object based on custom spectral density values
#   for the noises
# - calling getobs to get an array of X values at equispaced times
# - getting the spectrum of the time series and writing it to disk

# import all the libraries that are needed

from synthlisa import *
import numpy

# we create a LISA geometry object corresponding to a stationary LISA
# with equal armlengths

originallisa = OriginalLISA(16.6782,16.6782,16.6782)

# create a TDInoise object with standard pseudorandom noises

stime = 4.0

originalnoise = TDInoise(originallisa,
                         stime, 2.5e-48, # all six proof masses get this
                         stime, 1.8e-37, # all six optical-path noises get this
                         stime, 1.1e-26) # all six laser noises get this

# create a TDInoise object with custom parameters for the eighteen
# here pm_1 (unstarred) is a hundred times worse, in power

badnoise = TDInoise(originallisa,
                    [PowerLawNoise(stime,256.0,2.5e-46,-2.0,1)] + [PowerLawNoise(stime,256.0,2.5e-48,-2.0,1) for i in range(5)],
                    [PowerLawNoise(stime,256.0,1.8e-37,2.0,1) for i in range(6)],
                    [PowerLawNoise(stime,256.0,1.1e-26,0.0,1) for i in range(6)])

# get the time series; since they're from different TDInoise objects,
# there's no point in doing them together

samples = 2**19  # 2**18 takes 22 s on a 1.25GHz G4

patches = 1024

noisegood = getobsc(samples,stime,originalnoise.Xm)

[noisebad,noisebad2,noisebad3] = numpy.transpose(getobsc(samples,stime,[badnoise.Xm,badnoise.Ym,badnoise.Zm]))

# compute spectra, and write to disk

myspecgood = spect(noisegood,stime,patches)
writearray('data/tdibadmass-good.txt',myspecgood[1:])

myspecbad  = spect(noisebad,stime,patches)
writearray('data/tdibadmass-bad.txt', myspecbad[1:])

myspecbad2 = spect(noisebad2,stime,patches)
writearray('data/tdibadmass-bad2.txt',myspecbad2[1:])

myspecbad3 = spect(noisebad3,stime,patches)
writearray('data/tdibadmass-bad3.txt',myspecbad3[1:])
