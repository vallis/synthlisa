#!/usr/bin/python

# test of laser-noise cancellation with imperfect knowledge of the
# armlengths: second-generation TDI, additional error budgets

from lisaswig import *
from Numeric import *
from lisautils import *

# create the stationary LISA and standard TDI noise objects

oneyear = 31536000.0
mylisa = EccentricInclined(0.0,0.0,1.0,0.0*oneyear)

stime = 8

prebuffer = stime * 32

proofnoises = [InterpolateNoise(stime,prebuffer,2.5e-48,-2.0) for i in range(0,6)]
shotnoises  = [InterpolateNoise(stime,prebuffer,1.8e-37, 2.0) for i in range(0,6)]
lasernoises = [InterpolateNoise(stime,prebuffer,1.1e-26, 0.0) for i in range(0,6)]

# set the measurement noise (white) from its variance, and bandlimit frequency
#
#       dx^2 = (10 m)^2 = (3.33e-8 s)^2 = S_h f_b
#
# hence S_h = 1.11e-15 (dx/10 m)^2 (f_b/1 Hz)^-1 s^2/Hz
#
# using dx = 50 m, f_b = 1/(2*16) s, S_h = 8.88e-13 s^2/Hz

intl = 16.0
band = 1/(2*intl)
ersh = lambda errm : (errm / 3.0e8)**2 / band

noisylisa = NoisyLISA(mylisa,intl,ersh(2.0/3.0))

noisyTDI = TDIaccurate(noisylisa,proofnoises,shotnoises,lasernoises)

cleanerlisa = NoisyLISA(mylisa,intl,ersh(1.0/3.0))

cleanerTDI = TDIaccurate(cleanerlisa,proofnoises,shotnoises,lasernoises)

# ok, let's compute everything!

samples = 2**25/stime # on a 1.25GHz G4, 2**18 samples take 36 s
windows = 512

[noisenoisy,
 noisecleaner] = transpose(getobsc(samples,stime,[noisyTDI.X1,cleanerTDI.X1]))

writebinary('data/tdinoisy-2nd-cleanerer.bin',noisenoisy)
writebinary('data/tdinoisy-2nd-cleanest.bin',noisecleaner)

myspecnoisy = spect(noisenoisy,stime,windows)
writearray('data/tdinoisy-2nd-cleanerer.txt',myspecnoisy[1:])

myspeccleaner = spect(noisecleaner,stime,windows)
writearray('data/tdinoisy-2nd-cleanest.txt',myspeccleaner[1:])
