from Numeric import *
from FFT import *

# estimate spectrum
# patches = 0 gives the unwindowed spectrum
# patches = 1 gives the triangle-windowed spectrum
# patches > 1 gives the triangle-windowed, averaged over "patches" periods

# detrend controls the subtraction of DC components

def spect(series,sampling,patches=1,detrend=0):
    nyquistf = 0.5 / sampling

    if patches==0:
        period = pdg(series)
    elif patches==1:
        period = wpdg(series,detrend)
    else:
        period = opwpdg(series,patches,detrend)

    pdlen = shape(period)[0]-1

    freqs = arange(0,pdlen+1,typecode='d') * (nyquistf / pdlen)

    deltaf = nyquistf / pdlen
   
    period[0] = 2 * period[0] / deltaf
    period[1:pdlen] = period[1:pdlen] / deltaf
    period[pdlen] = 2 * period[pdlen] / deltaf

    spectrum = zeros((pdlen+1,2),typecode='d')

    spectrum[:,0] = freqs[:]
    spectrum[:,1] = period[:]

    return spectrum

# simple periodogram

def pdg(series):
    samples = shape(series)[0]

    pdlen = samples/2

    fourier = abs(fft(series))

    pdgram = zeros(pdlen+1,typecode='d')

    pdgram[0] = fourier[0]**2
    pdgram[1:pdlen] = fourier[1:pdlen]**2 + fourier[-1:pdlen:-1]**2
    pdgram[pdlen] = fourier[pdlen]**2

    pdgram /= (1.0*samples**2)

    return pdgram

# triangle-windowed periodogram
# - since series is a reference to successive overlapping slices,
#   we should not modify its value

def wpdg(series,detrend=0):
    samples = shape(series)[0]
    pdlen = samples/2

    window = 1.0 - abs(arange(0,samples,typecode='d') - pdlen) / (pdlen)
    weight = samples * sum(window ** 2)

    # detrending
    if detrend==0:
        mean = 0.0
    else:
        mean = sum(series) / (1.0*samples)

    wseries = window * (series - mean)

    wpdgram = pdg(wseries) * (1.0 * samples**2 / weight)

    return wpdgram

# overlapping-averaged, triangle-windowed periodogram

def opwpdg(series,patches,detrend=0):
    samples = shape(series)[0]
    serlen = samples - (samples % (4*patches))

    patlen = serlen/patches
    pdlen = patlen/2

    opwpdgram = zeros(pdlen+1,typecode='d')

    for cnt in range(0,2*patches-1):
        opwpdgram[:] += wpdg(series[cnt*pdlen:(cnt+2)*pdlen],detrend)

    opwpdgram[:] /= (2.0*patches - 1.0)

    return opwpdgram
    
# make it so that "time" is an observable

def time(t):
    return t

# return a 1 or 2 dimensional Numeric array, according to whether
# "observables" is a function or a list of functions

def getobs(snum,stime,observables):
    if len(shape(observables)) == 0:
        array = zeros(snum,typecode='d')
        for i in arange(0,snum):
            array[i] = observables(i*stime)
    else:
        obslen = shape(observables)[0]
        array = zeros((snum,obslen),typecode='d')
        for i in arange(0,snum):
            for j in range(0,obslen):
                array[i,j] = observables[j](i*stime)
    return array

# used by getobsc (hoping time.time() will work on all platforms...)

import sys
from time import time

def dotime(i,snum,inittime,lasttime):
    currenttime = int(time()) - inittime

    if currenttime - lasttime > 2:
        percdone = int((100.0*i)/snum)
        timeleft = int(((1.0 * currenttime) / i) * (snum-i))

        vel = ((1.0*i)/currenttime)

        minleft = timeleft / 60
        secleft = timeleft - minleft*60

        print "\r...%d/%d (%d%%) done [%d (multi)samples/s], est. %dm%ds left..." % (i,snum,percdone,vel,minleft,secleft),
        sys.stdout.flush()

        return currenttime

    return lasttime

# the next version, getobsc, will display a countdown to completion

def getobsc(snum,stime,observables):
    inittime = int(time())
    lasttime = 0
    
    print "Processing...",
    sys.stdout.flush()

    if len(shape(observables)) == 0:
        array = zeros(snum,typecode='d')
        for i in arange(0,snum):
            array[i] = observables(i*stime)
            if i % 1024 == 0:
                lasttime = dotime(i,snum,inittime,lasttime)
    else:
        obslen = shape(observables)[0]
        array = zeros((snum,obslen),typecode='d')
        for i in arange(0,snum):
            for j in range(0,obslen):
                array[i,j] = observables[j](i*stime)
            if i % 1024 == 0:
                lasttime = dotime(i,snum,inittime,lasttime)

    currenttime = int(time()) - inittime
    print "\r...completed in %d s.                                                 " % currenttime
    
    return array

# this version is less efficient, probably because of the conversion
# between the result of map (a sequence) and the numeric array

def getobs2(snum,stime,observables):
    if len(shape(observables)) == 0:
        return asarray(map(observables,arange(0,snum*stime,stime)),typecode='d')
    else:
        def mapfun(x):
            return map(lambda y: y(x),observables)
        return asarray(map(mapfun,arange(0,snum*stime,stime)),typecode='d')

def writeobs(filename,snum,stime,observables):
    writearray(filename,getobs(snum,stime,observables))

def writearray(filename,a):
    file = open(filename, 'w')
    if len(a.shape) == 1:
	a = a[:, Numeric.NewAxis]
    for line in a:
	for element in line:
	    file.write('%le ' % element)
	file.write('\n')
    file.close()

# still not satisfactory; should get the length from the file,
# not from the user, and should allocate its own buffer

def readarray(buffer,length,filename):
    file = open(filename, 'r')
    for index in range(0,length):
        buffer[index] = float(file.readline())
    file.close()

# how is the record order going to be?

def writebinary(filename,a):
    file = open(filename, 'w')
    file.write(a.tostring())
    file.close()

# in a future version, this should also get the length from the file

def readbinary(filename,length):
    file = open(filename,'r')
    buffer = Numeric.fromstring(file.read(length*8),'double')
    file.close()
    # then reshape the buffer if needed
    return buffer
