#!/usr/bin/python

from lisaswig import *
from Numeric import *
from lisautils import *

# create the stationary LISA and standard TDI noise objects

oneyear = 31536000.0
mylisa = EccentricInclined(0.0,0.0,1.0,0.0*oneyear)

stime = 8

cleanTDI = TDInoise(mylisa,
                    stime, 2.5e-48, # proof-mass noise parameters
                    stime, 1.8e-37, # optical-path noise parameters
                    stime, 1.1e-26) # laser frequency noise parameters

# set the measurement noise (white) from its variance, and bandlimit frequency
#
#       dx^2 = (10 m)^2 = (3.33e-8 s)^2 = S_h f_b
#
# hence S_h = 1.11e-15 (dx/10 m)^2 (f_b/1 Hz)^-1 s^2/Hz
#
# using dx = 50 m, f_b = 1/(2*16) s, S_h = 8.88e-13 s^2/Hz

noisylisa = NoisyLISA(mylisa,16.0,8.88e-13)

noisyTDI = TDInoise(noisylisa,
                    stime, 2.5e-48, # proof-mass noise parameters
                    stime, 1.8e-37, # optical-path noise parameters
                    stime, 1.1e-26) # laser frequency noise parameters

cleanerlisa = NoisyLISA(mylisa,16.0,0.1*8.88e-13)

cleanerTDI = TDInoise(cleanerlisa,
                      stime, 2.5e-48, # proof-mass noise parameters
                      stime, 1.8e-37, # optical-path noise parameters
                      stime, 1.1e-26) # laser frequency noise parameters

cleanerTDI.lock(1)

# ok, let's compute everything!

samples = 2**25/stime # on a 1.25GHz G4, 2**18 samples take 36 s
windows = 512

[noiseclean,
 noisenoisy,
 noisecleaner] = transpose(getobs(samples,stime,[cleanTDI.X1,noisyTDI.X1,cleanerTDI.X1]))

myspecclean = spect(noiseclean,stime,windows)
writearray('data/tdinoisy-2nd-clean.txt',myspecclean[1:])

myspecnoisy = spect(noisenoisy,stime,windows)
writearray('data/tdinoisy-2nd-noisy.txt',myspecnoisy[1:])

myspeccleaner = spect(noisecleaner,stime,windows)
writearray('data/tdinoisy-2nd-cleaner.txt',myspeccleaner[1:])
