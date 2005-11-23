# these might be replaced with simple imports; otherwise all of
# Numeric is available within the lisaswig context

import Numeric
import FFT
import math

# estimate spectrum
# patches = 0 gives the unwindowed spectrum
# patches = 1 gives the triangle-windowed spectrum
# patches > 1 gives the triangle-windowed, averaged over "patches" periods

# overlap controls the use of overlapping or nonoverlapping averaging periods

# detrend controls the subtraction of DC components

def spect(series,sampling,patches=1,detrend=0,overlap=1):
    nyquistf = 0.5 / sampling

    if patches==0:
        period = pdg(series)
    elif patches==1:
        period = wpdg(series,detrend)
    else:
        if overlap==0:
            period = nopwpdg(series,patches,detrend)
        else:
            period = opwpdg(series,patches,detrend)

    pdlen = Numeric.shape(period)[0]-1

    freqs = Numeric.arange(0,pdlen+1,typecode='d') * (nyquistf / pdlen)

    deltaf = nyquistf / pdlen
   
    period[0] = 2 * period[0] / deltaf
    period[1:pdlen] = period[1:pdlen] / deltaf
    period[pdlen] = 2 * period[pdlen] / deltaf

    spectrum = Numeric.zeros((pdlen+1,2),typecode='d')

    spectrum[:,0] = freqs[:]
    spectrum[:,1] = period[:]

    return spectrum

# simple periodogram

def pdg(series):
    samples = Numeric.shape(series)[0]

    pdlen = samples/2

    fourier = abs(FFT.fft(series))

    pdgram = Numeric.zeros(pdlen+1,typecode='d')

    pdgram[0] = fourier[0]**2
    pdgram[1:pdlen] = fourier[1:pdlen]**2 + fourier[-1:pdlen:-1]**2
    pdgram[pdlen] = fourier[pdlen]**2

    pdgram /= (1.0*samples**2)

    return pdgram

# triangle-windowed periodogram
# - since series is a reference to successive overlapping slices,
#   we should not modify its value

def wpdg(series,detrend=0):
    samples = Numeric.shape(series)[0]

    #pdlen = (samples-1)/2.0
    #window = 1.0 - abs(arange(0,samples,typecode='d') - pdlen) / (pdlen)

    # try Blackman
    wrange = Numeric.arange(0,samples,typecode='d') / (samples - 1.0);
    window = 0.42 - 0.5 * Numeric.cos(2*math.pi*wrange) + 0.08 * Numeric.cos(4*math.pi*wrange)

    weight = samples * Numeric.sum(window ** 2)

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
    samples = Numeric.shape(series)[0]
    serlen = samples - (samples % (4*patches))

    patlen = serlen/patches
    pdlen = patlen/2

    opwpdgram = Numeric.zeros(pdlen+1,typecode='d')

    for cnt in range(0,2*patches-1):
        opwpdgram[:] += wpdg(series[cnt*pdlen:(cnt+2)*pdlen],detrend)

    opwpdgram[:] /= (2.0*patches - 1.0)

    return opwpdgram

# nonoverlapping-averaged, triangle-windowed periodogram

def nopwpdg(series,patches,detrend=0):
    samples = Numeric.shape(series)[0]
    serlen = samples - (samples % (4*patches))

    patlen = serlen/patches
    pdlen = patlen/2

    opwpdgram = Numeric.zeros(pdlen+1,typecode='d')

    for cnt in range(0,patches):
        opwpdgram[:] += wpdg(series[cnt*patlen:(cnt+1)*patlen],detrend)

    opwpdgram[:] /= 1.0*patches

    return opwpdgram

def whiten(series,patches=1):
    samples = Numeric.shape(series)[0]

    patlen = samples / patches
    
    for j in range(0,patches-1):
        integ = 0.0

        for i in range(0,patlen):
            next = series[j*patlen+i]
            series[j*patlen+i] = integ
            integ = integ+next

def darken(spectrum,stime):
    samples = Numeric.shape(spectrum)[0]

    for cnt in range(1,samples):
        spectrum[cnt,1] = spectrum[cnt,1] * (2*Numeric.sin(math.pi*spectrum[cnt,0]*stime))**2.0

# make it so that "time" is an observable

def time(t):
    return t

# return a 1 or 2 dimensional Numeric array, according to whether
# "observables" is a function or a list of functions

def getobs(snum,stime,observables,zerotime=0.0):
    if len(Numeric.shape(observables)) == 0:
        array = Numeric.zeros(snum,typecode='d')
        for i in Numeric.arange(0,snum):
            array[i] = observables(zerotime+i*stime)
    else:
        obslen = Numeric.shape(observables)[0]
        array = Numeric.zeros((snum,obslen),typecode='d')
        for i in Numeric.arange(0,snum):
            for j in range(0,obslen):
                array[i,j] = observables[j](zerotime+i*stime)
    return array

def testobs(snum,stime,observables,zerotime=0.0):
    if len(Numeric.shape(observables)) == 0:
        array = Numeric.zeros(snum,typecode='d')
        for i in Numeric.arange(0,snum):
            array[i] = observables(zerotime+i*stime)
    else:
        obslen = Numeric.shape(observables)[0]
        array = Numeric.zeros((snum,obslen),typecode='d')
        for i in Numeric.arange(0,snum):
            for j in range(0,obslen):
                array[i,j] = observables[j](zerotime+i*stime)
    return None

# used by getobsc (hoping time.time() will work on all platforms...)

import sys
from time import time

def dotime(i,snum,inittime,lasttime):
    currenttime = int(time()) - inittime

    if currenttime - lasttime > 2 and i > 0:
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

def getobsc(snum,stime,observables,zerotime=0.0):
    inittime = int(time())
    lasttime = 0
    
    print "Processing...",
    sys.stdout.flush()

    try:
        if len(Numeric.shape(observables)) == 0:
            array = Numeric.zeros(snum,typecode='d')
            for i in Numeric.arange(0,snum):
                array[i] = observables(zerotime+i*stime)
                if i % 1024 == 0:
                    lasttime = dotime(i,snum,inittime,lasttime)
        else:
            obslen = Numeric.shape(observables)[0]
            array = Numeric.zeros((snum,obslen),typecode='d')
            for i in Numeric.arange(0,snum):
                for j in range(0,obslen):
                    array[i,j] = observables[j](zerotime+i*stime)
                if i % 1024 == 0:
                    lasttime = dotime(i,snum,inittime,lasttime)
    except IndexError:
        print "lisautils::getobsc: I have trouble accessing time ", zerotime+i*stime,
        print "; you may try to reset your objects and repeat..."

        raise

    # there was good stuff here, but it's better to let the user deal with this
    # in particular, if observables is a noise-like variable, then it should
    # have the 'reset' method ['reset' in dir(observables)]; if observables is
    # a method of a TDI class, which should have the 'reset' method, the test
    # is ['reset' in dir(observables.im_self)], where "im_self" returns the
    # class instance for a given instance method

    currenttime = int(time()) - inittime

    if currenttime > 0:
        vel = ((1.0*snum)/currenttime)
        print "\r...completed in %d s [%d (multi)samples/s].                           " % (currenttime,vel)  
    else:
        print "\r...completed.                                                         "

    return array

# this version is less efficient, probably because of the conversion
# between the result of map (a sequence) and the numeric array

def getobs2(snum,stime,observables):
    if len(Numeric.shape(observables)) == 0:
        return Numeric.asarray(map(observables,Numeric.arange(0,snum*stime,stime)),typecode='d')
    else:
        def mapfun(x):
            return map(lambda y: y(x),observables)
        return Numeric.asarray(map(mapfun,Numeric.arange(0,snum*stime,stime)),typecode='d')

def writeobs(filename,snum,stime,observables):
    writearray(filename,getobs(snum,stime,observables))

# should be able to convert from a sequence to an array,
# if an array is not given

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
# can probably do it with append()

def readarray(filename):
    file = open(filename, 'r')
    lines = [line.split() for line in file.readlines() if line[0] != '#']    
    file.close()
    
    bshape = (len(lines),len(lines[0]))

    if bshape[1] == 1:
        buffer = Numeric.zeros(bshape[0],'d')

        for index in range(0,bshape[0]):
            buffer[index] = float(lines[index][0])
    else:
        buffer = Numeric.zeros(bshape,'d')

        for index in range(0,buffer.shape[0]):
            buffer[index,:] = map(float,lines[index])

    return buffer

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


# healpix (need to sort out the license)

import healpix

def hnpix(nside):
    return healpix.healpix.nside2npix(nside)

def hn2ec(nside,ipix):
    (th,ph) = healpix.healpix.pix2ang_nest(nside,ipix)
    return (math.pi-th,ph)
    
def hr2ec(nside,ipix):
    (th,ph) = healpix.healpix.pix2ang_ring(nside,ipix)
    return (math.pi-th,ph)


# lisa positions from Ted Sweetser's file

import os
import os.path

def stdLISApositions():
    """Returns four Numeric arrays corresponding to times [in seconds] along
ten years and to the corresponding positions [the three SSB coordinates
in Earth Mean Equator and J2000 Equinox, also given in seconds] of the
three LISA spacecraft according to simulations run by JPL's Ted Sweetser
on 2005-07-23 (from Excel spreadsheet states_baseline2.xls). The times
are spaced by 1 day (86400 seconds), and they begin from zero instead of
Ted's Julian date 2457023.5"""

    pos = readarray(os.path.join(os.environ['SYNTHLISABASE'],'share/synthlisa','positions.txt'))

    t = pos[:,0].copy()

    p1 = pos[:,1:4].copy()
    p2 = pos[:,4:7].copy()
    p3 = pos[:,7:10].copy()

    secondsperday = 86400

    t -= t[0]
    t *= secondsperday

    speedoflight = 299792.458

    p1 /= speedoflight
    p2 /= speedoflight
    p3 /= speedoflight

    return t,p1,p2,p3


import lisaswig

def stdSampledLISA(interp=1):
"""Returns an interpolated SampledLISA object based on the position arrays
returned by stdLISApositions(); the argument interp sets the semilength
of the interpolation window, and prebuffering is set to interp times 1
day (the spacing of the stdLISApositions() data."""

    [t,p1,p2,p3] = stdlisapositions()
    
    return lisaswig.SampledLISA(p1,p2,p3,86400,86400*interp,interp)
