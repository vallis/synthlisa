#!/usr/bin/python

# test of first-generation TDI X noise against theory

from lisaswig import *
from Numeric import *
from lisautils import *

stime = 1
samples = 2**18 / stime

patches = 64

# create the stationary LISA and standard TDI noise objects

originallisa = OriginalLISA(16.6782,16.6782,16.6782)

originalTDI = TDInoise(originallisa,
                       stime, 2.5e-48, # proof-mass noise parameters
                       stime, 1.8e-37, # optical-path noise parameters
                       stime, 1.1e-26) # laser frequency noise parameters

lockedTDI = TDInoise(originallisa,
                     stime, 2.5e-48, # proof-mass noise parameters
                     stime, 1.8e-37, # optical-path noise parameters
                     stime, 1.1e-26) # laser frequency noise parameters

lockedTDI.lock(3)

# ok, let's compute everything!

#tdiy231 = lambda t : lockedTDI.y(2,-1,3,0,0,0,0,0,0,0,t)
#noisey231 = getobs(samples,stime,tdiy231)

[noiseX, noiseXlock, noiseXsimp] = transpose(getobs(samples,stime,[originalTDI.Xm,lockedTDI.Xm,lockedTDI.Xmlock3]))

myspecX =  spect(noiseX, stime,patches)
writearray('tdi-X-nonlocked.txt', myspecX[1:])

myspecXlock = spect(noiseXlock,stime,patches)
writearray('tdi-X-locked.txt',myspecXlock[1:])

myspecXsimp = spect(noiseXsimp,stime,patches)
writearray('tdi-X-locksimp.txt',myspecXsimp[1:])

#myspecy231 = spect(noisey231,stime,patches)
#writearray('tdi-y-zero.txt',myspecy231[1:])
