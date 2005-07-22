#!/usr/bin/python

import pypar

from synthlisa import *
from Numeric import *

from lisapar import *

import math
import random

numsources = 6

stime = 16
samples = 2**12

# define sources (root only)

lofreq = 1e-3 # 1 mHz
hifreq = 2e-3 # 2 mHz

def randomph():
    return 2.0 * math.pi * random.random()

def randomth():
    cosz = 2.0 * random.random() - 1.0
    return math.acos(cosz) - 0.5 * math.pi

parameters = []

if pypar.rank() == 0:
    for i in range(0,numsources):
        parameters.append( [ lofreq + (hifreq - lofreq) * random.random(),
                                          # frequency
                             randomph(),  # phase
                             randomth(),  # inclination
                             1.0,         # amplitude
                             randomth(),  # ecliptic latitude
                             randomph(),  # ecliptic longitude
                             randomph()   # polarization
                             ] )

circularlisa  = CircularRotating(0.0,0.0,1.0)

# I have a problem with the persistence of objects: I cannot just return TDIsignal(...,SimpleBinary(...))
# because SimpleBinary will be deallocated and is needed by TDIsignal; fix in SWIG?

def sourcefn(pars):
    lisa = circularlisa
    source = SimpleBinary(*pars)
    tdi = TDIsignal(lisa,source)

    return (lisa,source,tdi)

[signalX,signalY] = transpose(getobsp(samples,stime,[sourcefn,parameters,(TDIsignal.X,TDIsignal.Y)]))

if pypar.rank() == 0:
    writearray('data/tdibinary-parallel-X.txt',signalX)
    writearray('data/tdibinary-parallel-Y.txt',signalY)

    # signalsumX = zeros(shape(signalX),typecode='d')
    # signalsumY = zeros(shape(signalY),typecode='d')

    # for par in parameters:
    #    (lisa,source,tdi) = sourcefn(par)

    #    [signalX,signalY] = transpose(getobsc(samples,stime,(tdi.X,tdi.Y)))
    #    signalsumX += signalX
    #    signalsumY += signalY

    # writearray('data/tdibinary-parallel-X2.txt',signalsumX)
    # writearray('data/tdibinary-parallel-Y2.txt',signalsumY)

pypar.finalize()
