#!/usr/bin/python

# test of laser-noise cancellation for first- and second-generation TDI observables

# note: the plotting scripts test-tdi2nd.plt, test-tdi2nd-2.plt,
# test-tdi2nd-3.plt, test-tdi2nd-4.plt, all use the output of this Python script
# to show different aspects of the result

# this script demonstrates the generation of 1st- and 2nd-generation
# TDI time series, and in particular:
# - creating the appropriate LISA geometry objects
# - creating TDInoise objects based on standard pseudorandom noises
# - calling getobs to get the time series
# - getting the spectra of the time series
# - writing the spectra to disk

# import all the libraries that are needed

from lisaswig import *
from Numeric import *
from lisautils import *

# we use the most realistic LISA geometry model, with eccentric and inclined
# orbits that give time-changing and direction-dependent armlengths

# add a time offset at the end to see a different period

oneyear = 31536000.0
eccentriclisa = EccentricInclined(0.0,0.0,1.0,0.0*oneyear)

# we use the same sampling time for noises and for TDI time series

stime = 8

# create TDInoise object

eccentricTDI = TDInoise(eccentriclisa,
                       stime, 2.5e-48, # proof-mass noise parameters
                       stime, 1.8e-37, # optical-path noise parameters
                       stime, 1.1e-26) # laser frequency noise parameters

# for comparison, we'll take out the laser noises

nolaserTDI = TDInoise(eccentriclisa,
                      stime, 2.5e-48, # proof-mass noise parameters
                      stime, 1.8e-37, # optical-path noise parameters
                      stime, 0.0)     # laser frequency noise parameters

# how many samples do we need? to do one year at 1 s we need 31,536,000
# samples, or about 2**25. It's a bit much, both memory wise and
# CPU wise (on a Mac G4 1.25 GHz, I get about 13,000 samples/s)

samples = 2**25 / stime
#samples = 2**18 / stime

patches = 256

# ok, let's compute everything! modified X and X1

[noiseX, noiseX1] = transpose(getobs(samples,stime,[eccentricTDI.Xm,eccentricTDI.X1]))

myspecX =  spect(noiseX, stime,patches)
writearray('data/tdi2nd-X.txt', myspecX[1:])

myspecX1 = spect(noiseX1,stime,patches)
writearray('data/tdi2nd-X1.txt',myspecX1[1:])

# same, with no laser noise; we reuse the arrays to save some space

[noiseX, noiseX1] = transpose(getobs(samples,stime,[nolaserTDI.Xm,nolaserTDI.X1]))

myspecXc =  spect(noiseX, stime,patches)
writearray('data/tdi2nd-X-nolaser.txt', myspecXc[1:])

myspecX1c = spect(noiseX1,stime,patches)
writearray('data/tdi2nd-X1-nolaser.txt',myspecX1c[1:])

# now, I'd like to see noise cancellation for X with reduced laser noise

somelaserTDI = TDInoise(eccentriclisa,
                        stime, 2.5e-48, # proof-mass noise parameters
                        stime, 1.8e-37, # optical-path noise parameters
                        stime, 0.1*1.1e-26) # laser frequency noise parameters

noiseX = getobs(samples,stime,somelaserTDI.Xm)

myspecXb =  spect(noiseX, stime,patches)
writearray('data/tdi2nd-X-somelaser.txt', myspecXb[1:])

# change it again?

lesslaserTDI = TDInoise(eccentriclisa,
                        stime, 2.5e-48, # proof-mass noise parameters
                        stime, 1.8e-37, # optical-path noise parameters
                        stime, 0.01*1.1e-26) # laser frequency noise parameters

noiseX = getobs(samples,stime,lesslaserTDI.Xm)

myspecXl =  spect(noiseX, stime,patches)
writearray('data/tdi2nd-X-lesslaser.txt', myspecXl[1:])

# and now, we'd like to see the ratio between the "dirty" and "clean" noise curves

freqs = myspecX[1:,0]
ratio1 = sqrt(myspecX[1:,1]  / myspecXc[1:,1])
ratio2 = sqrt(myspecXb[1:,1] / myspecXc[1:,1])
ratio3 = sqrt(myspecXl[1:,1] / myspecXc[1:,1])

writearray('data/tdi2nd-comparison.txt',transpose([freqs, ratio1, ratio2, ratio3]))
