#!/usr/bin/python

# test of the standard LISA laser frequency noise

# this script demonstrates the generation of pseudorandom noise with
# an InterpolateNoise object, and in particular:
# - creating an InterpolationNoise object
# - calling getobs to get an array of noise values at equispaced times
# - calling spect to get the spectrum of the time series
# - writing the spectrum to disk

# import all the libraries that are needed

from lisaswig import *
from Numeric import *
from lisautils import *

# here we create a proof-mass noise object with
# - noise-generation sampling time 1.0 s
# - zero offset 256 s (quite arbitrary here)
# - noise spectral density S_n = 1.1e-26 1/Hz
# - noise spectral density exponent 0.0
# - linear interpolation (1)

lasernoise = InterpolateNoise(1.0,256.0,1.1e-26,0.0,1)

# get "samples" values of the noise, at times separated by "sstep"

samples = 2**20 # 2**20 takes 13s on a 1.25GHz G4; 2*22 used for plot
sstep = 0.1

patches = 64

noises = getobs(samples,sstep,lasernoise.noise)

# the result is a 1D Numeric array, that we feed to "spect" to get a
# triangle-windowed, averaged spectrum, using "patches" averaging
# periods

myspec = spect(noises,sstep,patches)

# the result is a 2D Numeric array (for instance, myspec[0,0] = 0 Hz,
# and myspec[0,1] is the spectral density at DC; myspec[1,0] is the
# frequency of the first bin, and myspec[1,1] is the spectral density
# at that frequency

# we write it to disk as an ASCII file with (frequency,density) pairs
# on each line; we skip the DC component which is confusing to gnuplot

writearray('data/lasernoise-freq.txt',myspec[1:])
