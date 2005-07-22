#!/usr/bin/python

# test of first-generation TDI signals for a monochromatic-binary GW source

# this script demonstrates the generation of a TDI time series for a
# GW signal, and in particular:
#
# - creating a LISA geometry object
# - creating a Wave GW object
# - creating a TDIsignal object based on the Wave object
# - creating a standard TDInoise object
# - calling getobs to get an array of the TDI X observable at equispaced times
# - calling getobs to get an array of the TDI X noise at equispaced times
# - computing the optimal S/N for the source

# import all the libraries that are needed

from synthlisa import *
from Numeric import *


###################################################
# creating the GW source from physical parameters #
###################################################

# the following parameters are taken from the standard test source
# for the LISA simulator, v. 2.0 (Newtonian.c)

import math

# masses of the binary components (kg)

Msun = 1.9889e30

M1 = 0.5*Msun
M2 = 0.033*Msun

# luminosity distance (m)

kpc = 3.0856675807e19
r = 0.1*kpc

# orbital frequency (Hz)

forb = 1.0/1028.76
montanafgw = 2.0*forb

# find overall amplitude (G is mks, c is m/s)

G = 6.67259e-11
c = 2.99792458e8

R = (G*(M1+M2)*(1.0/(2.0*math.pi*forb))**2)**(1.0/3.0)

montanaamplitude = (2.0*M1*M2/(r*R) * (G/(c*c))**2)

# initial phase, inclination, position in the sky, and polarization of the source

montanaphi0 = 0.0
montanainc = 1.53921

montanacolatitude = 1.57
montanalongitude = 0.0

montanapsi = 0.0

# translation between Synthetic LISA and LISA Simulator GW source parameters:
# - the frequency needs to be adjusted for the different definition of year
#   in the two codes *1.000702365550482
# - difference in the definition of inclination
# - Synthetic LISA uses SSB ecliptic latitude rather than colatitude
# - difference in the definition of polarization

# by the way, the following is a Python dictionary: a bit like an
# array, but you can retrieve values by "key" (in this case, parameter
# name)

mypars = {'frequency': montanafgw,
          'phase':     montanaphi0,
          'inc':       math.pi - montanainc,
          'amplitude': montanaamplitude,
          'beta':      0.5*math.pi - montanacolatitude,
          'lambda':    montanalongitude,
          'psi':       -montanapsi}

#####################################
# creating the LISA and TDI objects #
#####################################

# create a SimpleBinary Wave object

mybinary = SimpleBinary(mypars['frequency'],
                        mypars['phase'],
                        mypars['inc'],
                        mypars['amplitude'],
                        mypars['beta'],
                        mypars['lambda'],
                        mypars['psi'])

# these parameters yield the same initial LISA orientation and position
# as in the LISA Simulator, v. 2.0

circularlisa  = CircularRotating(0.0,3.0*math.pi/2.0,-1.0)
eccentriclisa = EccentricInclined(0.0,3.0*math.pi/2.0,-1.0)

# TDIsignal and TDInoise structures

circularsignal = TDIsignal(circularlisa,mybinary)
eccentricsignal = TDIsignal(circularlisa,mybinary)

# need to know at what frequency we're sampling TDI not to alias too badly

stime = 16

# let's forget laser noise here, and pretend X can do a good job at canceling it

circularnoise = TDInoise(circularlisa,stime,2.5e-48,stime,1.8e-37,stime,0*1.1e-26)

################################
# computing signals and noises #
################################

# one year of data is approx 2**25 samples; we use getobs to produce the signals
# and the noises

samples = 2**17  # 2**17 takes 42s on a 1.25GHz G4
                 # the graph however was done with stime = 16 and samples = 2**21
                 # and test-montana.py also needs 2**21
patches = 16

# "getobsc" (as opposed to "getobs") will display completion and speed date

signalXcirc = getobsc(samples,stime,circularsignal.X)
noiseXcirc  = getobsc(samples,stime,circularnoise.Xm)

# write the binary time series

writebinary('data/tdibinary-X-circ.bin',signalXcirc)
writebinary('data/tdibinary-X-noise.bin',noiseXcirc)

# compute signal spectrum without windowing or averaging

sspec = spect(signalXcirc,stime,0)

# on the other hand, average the noise profile

nspec = spect(noiseXcirc ,stime,patches)

# interpolate the noise to the frequencies of the signal spectrum
# (arrayfns is part of Numeric)

ispec = zeros(shape(sspec),typecode='d')
ispec[:,0] = sspec[:,0]

import arrayfns
ispec[:,1] = arrayfns.interp(nspec[:,1],nspec[:,0],ispec[:,0])

# the (S/N)^2 is be given by 2T times the integrated ratio
# of the spectral densities (the factor of 2 because the spectral
# density is one-sided); notice however that the df is 1/T,
# so we need only to sum up the array containing the ratio,
# and multiply by two

sratio = zeros(shape(sspec)[0],typecode='d')
sratio[1:] = sspec[1:,1] / ispec[1:,1]
sn2 = 2.0 * sum(sratio[1:])

print "S/N = %f over %f seconds" % (sqrt(sn2),1.0*stime*samples)

# it is hard to plot the signal meaningfully over noise: plotting
# twice the signal spectral density on top of the noise spectral
# density would allow us to read (as a ratio between the curves) the
# (S/N)^2 IN EACH BIN; however, given the LISA amplitude modulation,
# S/N is spread on a very short frequency scale (a few bins) around
# the central frequency of the monochromatic signal

writearray('data/tdibinary-X-noise.txt',nspec[1:])
writearray('data/tdibinary-X-circ.txt', 2.0*sspec[1:])
