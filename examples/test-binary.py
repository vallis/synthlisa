#!/usr/bin/python

# test of first-generation TDI signals for a monochromatic-binary GW source

# this script demonstrates the generation of a TDI time series for a
# GW signal, and in particular:
#
# - creating a LISA geometry object
# - creating a Wave GW object
# - creating a TDIsignal object based on the Wave object
# - calling getobs to get an array of the TDI X observable at equispaced times
# - creating a standard TDInoise object
# - calling getobs to get an array of the TDI X noise at equispaced times
# - summing the noise and signal time series
# - getting spectra and writing them to disk

# import all the libraries that are needed

from lisaswig import *
from Numeric import *
from lisautils import *

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
#   in the two codes
# - difference in the definition of inclination
# - Synthetic LISA uses SSB ecliptic latitude rather than colatitude
# - difference in the definition of polarization

# by the way, the following is a Python dictionary: a bit like an
# array, but you can retrieve values by "key" (in this case, parameter
# name)

mypars = {'frequency': montanafgw*1.000702365550482,
          'phase':     montanaphi0,
          'inc':       math.pi - montanainc,
          'amplitude': montanaamplitude,
          'lambda':    0.5*math.pi - montanacolatitude,
          'beta':      montanalongitude,
          'psi':       -montanapsi}

#####################################
# creating the LISA and TDI objects #
#####################################

# create a SimpleBinary Wave object

mybinary = SimpleBinary(mypars['frequency'],
                        mypars['phase'],
                        mypars['inc'],
                        mypars['amplitude'],
                        mypars['lambda'],
                        mypars['beta'],
                        mypars['psi'])

# these parameters yield the same initial LISA orientation and position
# as in the LISA Simulator, v. 2.0

circularlisa  = CircularRotating(0.0,3.0*math.pi/2.0,-1.0)
eccentriclisa = EccentricInclined(0.0,3.0*math.pi/2.0,-1.0)

# TDIsignal and TDInoise structures

circularsignal = TDIsignal(circularlisa,mybinary)
eccentricsignal = TDIsignal(circularlisa,mybinary)

# need to know at what frequency we're sampling TDI not to alias too badly

stime = 32

# let's forget laser noise here, and pretend X can do a good job at canceling it

circularnoise = TDInoise(circularlisa,stime,2.5e-48,stime,1.8e-37,stime,0*1.1e-26)

################################
# computing signals and noises #
################################

# one year of data is approx 2**25 samples; we use getobs to produce the signals
# and the noises

samples = 2**17  # we want to do a full year, but 2**17 takes 42s on a 1.25GHz G4
                 # the graph however was done with stime = 16 and samples = 2**21

patches = 16

signalXcirc = getobs(samples,stime,circularsignal.X)
signalXecc  = getobs(samples,stime,eccentricsignal.X)

noiseXcirc  = getobs(samples,stime,circularnoise.Xm)

# write the binary time series

writebinary('data/tdibinary-X-circ.bin',signalXcirc)
writebinary('data/tdibinary-X-noise.bin',signalXcirc)

# compute the spectra and write them to disk

myspecXcirc = spect(signalXcirc,stime,patches)
writearray('data/tdibinary-X-circ.txt',myspecXcirc[1:])

myspecXecc = spect(signalXecc,stime,patches)
writearray('data/tdibinary-X-ecc.txt',myspecXecc[1:])

myspecXnoise = spect(noiseXcirc,stime,patches)
writearray('data/tdibinary-X-noise.txt',myspecXnoise)

# for closeups, we do not average the signal spectra

myspecXcirc = spect(signalXcirc,stime)
writearray('data/tdibinary-X-circ-closeup.txt',myspecXcirc[1:])

myspecXecc = spect(signalXecc,stime)
writearray('data/tdibinary-X-ecc-closeup.txt',myspecXecc[1:])

# compute the S/N of the source

mysn = sn(signalXcirc,stime,myspecXnoise)

print 'The S/N for this source is', mysn
print 'Extrapolated to a year: ', sqrt(31536000.0 / (samples*stime)) * mysn

mysn2 = sn2(stime*samples,myspecXcirc,myspecXnoise)

print 'Alternative computation', mysn2

# if you want to see the S/N graphically, you have to multiply the fourier
# transform by frequency, and plot if against the noise S_n

# it seems that the LISA Simulator source has a very high S/N of about 36 per year
