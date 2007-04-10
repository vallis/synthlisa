#!/usr/bin/env python

# test of the standard LISA proof-mass noise

# this script demonstrates the generation of pseudorandom noise with
# an InterpolateNoise object, and in particular:
# - creating an InterpolationNoise object
# - calling getobs to get an array of noise values at equispaced times
#   (with a comment on how a standard "for" would get the same result)
# - calling spect to get the spectrum of the time series
# - writing the spectrum to disk

# import all the libraries that are needed

from synthlisa import *

# here we create a proof-mass noise object with
# - noise-generation sampling time 1.0 s
# - zero offset 256 s (quite arbitrary here)
# - noise spectral density S_n = 2.5e-48 (f/Hz)^2 1/Hz
# - noise spectral density exponent -2.0
# - linear interpolation (1)

proofnoise = PowerLawNoise(1.0,256.0,2.5e-48,-2.0,interplen=1)

# get "samples" values of the noise, at times separated by "sstep"

samples = 2**20  # 2**20 takes 13s on a 1.25GHz G4; 2*22 used for plot
sstep = 0.1

patches = 64

noises = getobsc(samples,sstep,proofnoise)

# "getobs" is really a short hand for creating an empty Numeric array,
# and then running through a standard Python cycle to call
# proofnoise.noise(time), as follows:
#
# noises = zeros(samples,typecode='d') # we want doubles
# for i in range(0,samples):
#     noises[i] = proofnoise.noise(i*sstep)

# the result is a 1D Numeric array, that we feed to "spect" to get a
# triangle-windowed, averaged spectrum, using "patches" averaging
# periods; use detrending to fix DC component

myspec = spect(noises,sstep,patches,detrend=1)

# the result is a 2D Numeric array (for instance, myspec[0,0] = 0 Hz,
# and myspec[0,1] is the spectral density at DC; myspec[1,0] is the
# frequency of the first bin, and myspec[1,1] is the spectral density
# at that frequency

# we write it to disk as an ASCII file with (frequency,density) pairs
# on each line; we skip the DC component which is confusing to gnuplot

writearray('data/proofnoise-freq.txt',myspec[1:])
