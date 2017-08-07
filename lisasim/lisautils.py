# $Id$
# $Date$
# $Author$
# $Revision$

import lisaswig
import numpy
import numpy
import numpy.fft as FFT
import math

# estimate spectrum
# patches = 0 gives the unwindowed spectrum
# patches = 1 gives the triangle-windowed spectrum
# patches > 1 gives the triangle-windowed, averaged over "patches" periods

# overlap controls the use of overlapping or nonoverlapping averaging periods

# detrend controls the subtraction of DC and linear components

def spect(series,sampling,patches=1,detrend=0,overlap=1,win='triangle'):
    nyquistf = 0.5 / sampling

    if patches==0:
        period = pdg(series)
    elif patches==1:
        period = wpdg(series,detrend,win)
    else:
        if overlap==0:
            period = nopwpdg(series,patches,detrend,win)
        else:
            period = opwpdg(series,patches,detrend,win)

    pdlen = numpy.shape(period)[0]-1

    freqs = numpy.arange(0,pdlen+1,dtype='d') * (nyquistf / pdlen)

    deltaf = nyquistf / pdlen
   
    period[0] = 2 * period[0] / deltaf
    period[1:pdlen] = period[1:pdlen] / deltaf
    period[pdlen] = 2 * period[pdlen] / deltaf

    spectrum = numpy.zeros((pdlen+1,2),dtype='d')

    spectrum[:,0] = freqs[:]
    spectrum[:,1] = period[:]

    return spectrum

# simple periodogram

def pdg(series):
    samples = numpy.shape(series)[0]

    pdlen = samples/2

    fourier = abs(FFT.fft(series))

    pdgram = numpy.zeros(pdlen+1,dtype='d')

    pdgram[0] = fourier[0]**2
    pdgram[1:pdlen] = fourier[1:pdlen]**2 + fourier[-1:pdlen:-1]**2
    pdgram[pdlen] = fourier[pdlen]**2

    pdgram /= (1.0*samples**2)

    return pdgram

# quick-and-dirty least squares fit
# returns the coefficients of the model a + b x
# with x = [0,...,N-1]
# if detrend = 1, also detrends the data

def leastsquares(array,detrend=0):
    S = len(array)
    x = numpy.arange(0,S,dtype='d')
    
    Sx = numpy.sum(x)
    Sy = numpy.sum(array)

    Sxx = numpy.sum(x**2)
    Sxy = numpy.sum(x * array)

    delta = S * Sxx - Sx**2
    
    a = (Sxx * Sy - Sx * Sxy) / delta
    b = (S * Sxy - Sx * Sy) / delta

    if detrend:
        array -= a + b * x
    
    return (a,b)

# triangle-windowed periodogram
# - since series is a reference to successive overlapping slices,
#   we should not modify its value

def wpdg(series,detrend=0,win='triangle'):
    samples = numpy.shape(series)[0]

    wrange = numpy.arange(0,samples,dtype='d') / (samples - 1.0);

    if win == 'blackman':
        window = 0.42 - 0.5 * numpy.cos(2*math.pi*wrange) + 0.08 * numpy.cos(4*math.pi*wrange)
    elif win == 'sin4':
        window = numpy.sin(math.pi*wrange)**4.0
    else:
        # if we don't recognize a window, default to triangle
        pdlen = (samples - 1) / 2.0
        window = 1.0 - abs(numpy.arange(0,samples,dtype='d') - pdlen) / (pdlen)

    wseries = series.copy()

    if detrend == 1:
        leastsquares(wseries,detrend=1)
    
    wseries *= window

    weight = samples * numpy.sum(window ** 2)
    wpdgram = pdg(wseries) * (1.0 * samples**2 / weight)

    return wpdgram

# overlapping-averaged, triangle-windowed periodogram

def opwpdg(series,patches,detrend=0,win='triangle'):
    samples = numpy.shape(series)[0]
    serlen = samples - (samples % (4*patches))

    patlen = serlen/patches
    pdlen = patlen/2

    opwpdgram = numpy.zeros(pdlen+1,dtype='d')

    for cnt in xrange(0,2*patches-1):
        opwpdgram[:] += wpdg(series[cnt*pdlen:(cnt+2)*pdlen],detrend,win)

    opwpdgram[:] /= (2.0*patches - 1.0)

    return opwpdgram

# nonoverlapping-averaged, triangle-windowed periodogram

def nopwpdg(series,patches,detrend=0,win='triangle'):
    samples = numpy.shape(series)[0]
    serlen = samples - (samples % (4*patches))

    patlen = serlen/patches
    pdlen = patlen/2

    opwpdgram = numpy.zeros(pdlen+1,dtype='d')

    for cnt in xrange(0,patches):
        opwpdgram[:] += wpdg(series[cnt*patlen:(cnt+1)*patlen],detrend,win)

    opwpdgram[:] /= 1.0*patches

    return opwpdgram

def whitentime(series,patches=1):
    samples = numpy.shape(series)[0]

    patlen = samples / patches
    
    for j in xrange(0,patches-1):
        integ = 0.0

        for i in xrange(0,patlen):
            next = series[j*patlen+i]
            series[j*patlen+i] = integ
            integ = integ+next

def darkenspectrum(spectrum,stime):
    samples = numpy.shape(spectrum)[0]

    for cnt in xrange(1,samples):
        spectrum[cnt,1] = spectrum[cnt,1] * (2*numpy.sin(math.pi*spectrum[cnt,0]*stime))**2.0

def darkentime(series,patches=1):
    samples = numpy.shape(series)[0]

    patlen = samples / patches
    
    for j in xrange(0,patches-1):
        for i in xrange(0,patlen-1):
            series[j*patlen+i] = series[j*patlen+i+1] - series[j*patlen+i]

def whitenspectrum(spectrum,stime):
    samples = numpy.shape(spectrum)[0]

    for cnt in xrange(1,samples):
        spectrum[cnt,1] = spectrum[cnt,1] / (2*numpy.sin(math.pi*spectrum[cnt,0]*stime))**2.0

# make it so that "time" is an observable

def time(t):
    return t

# return a 1 or 2 dimensional numpy array, according to whether
# "observables" is a function or a list of functions

import operator
import types

def checkobs(observables):
    retobs = []

    for obs in observables:
        if isinstance(obs,lisaswig.Signal):
            retobs.append(obs)
        elif isinstance(obs,types.MethodType) and isinstance(obs.im_self,lisaswig.TDI):
            try:
                retobs.append(obs())
            except:
                return None
        else:
            return None
            
    return retobs

def getobsc(snum,stime,observables,zerotime=0.0,forcepython=0):
    return getobs(snum,stime,observables,zerotime,display=1,forcepython=forcepython)

def getobs(snum,stime,observables,zerotime=0.0,display=0,forcepython=0):
    if len(numpy.shape(observables)) == 0:
        obsobj = checkobs([observables])
                
        if obsobj and (not forcepython):
            array = numpy.zeros(snum,dtype='d')

            if display:
                lisaswig.fastgetobsc(array,snum,stime,obsobj,zerotime)
            else:
                lisaswig.fastgetobs(array,snum,stime,obsobj,zerotime)
        else:
            if display:
                return getobscount(snum,stime,observables,zerotime)
            else:
                array = numpy.zeros(snum,dtype='d')
            
                for i in numpy.arange(0,snum):
                    array[i] = observables(zerotime+i*stime)
    else:
        obsobj = checkobs(observables)

        if obsobj and (not forcepython):
            obslen = numpy.shape(observables)[0]
            array = numpy.zeros((snum,obslen),dtype='d')

            if display:
                lisaswig.fastgetobsc(array,snum,stime,obsobj,zerotime)
            else:
                lisaswig.fastgetobs(array,snum,stime,obsobj,zerotime)
        else:
            if display:
                return getobscount(snum,stime,observables,zerotime)
            else:
                obslen = numpy.shape(observables)[0]
                array = numpy.zeros((snum,obslen),dtype='d')
            
                for i in numpy.arange(0,snum):
                    for j in xrange(0,obslen):
                        array[i,j] = observables[j](zerotime+i*stime)
    return array

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

def getobscount(snum,stime,observables,zerotime=0.0):
    fullinittime = time()
    inittime = int(fullinittime)
    lasttime = 0
    
    print "Processing...",
    sys.stdout.flush()

    try:
        if len(numpy.shape(observables)) == 0:
            array = numpy.zeros(snum,dtype='d')
            for i in numpy.arange(0,snum):
                array[i] = observables(zerotime+i*stime)
                if i % 1024 == 0:
                    lasttime = dotime(i,snum,inittime,lasttime)
        else:
            obslen = numpy.shape(observables)[0]
            array = numpy.zeros((snum,obslen),dtype='d')
            for i in numpy.arange(0,snum):
                for j in xrange(0,obslen):
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

    currenttime = time() - fullinittime

    if currenttime > 0:
        vel = snum/currenttime
        print "\r...completed in %d s [%d (multi)samples/s].                           " % (int(currenttime),int(vel))  
    else:
        print "\r...completed.                                                         "

    return array

# this version is less efficient, probably because of the conversion
# between the result of map (a sequence) and the numpy array

def getobs2(snum,stime,observables):
    if len(numpy.shape(observables)) == 0:
        return numpy.asarray(map(observables,numpy.arange(0,snum*stime,stime)),dtype='d')
    else:
        def mapfun(x):
            return map(lambda y: y(x),observables)
        return numpy.asarray(map(mapfun,numpy.arange(0,snum*stime,stime)),dtype='d')

def writeobs(filename,snum,stime,observables):
    writearray(filename,getobs(snum,stime,observables))

# should be able to convert from a sequence to an array,
# if an array is not given

def writearray(filename,a):
    file = open(filename, 'w')
    if len(a.shape) == 1:
        a = a[:, numpy.newaxis]
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
        buffer = numpy.zeros(bshape[0],'d')

        for index in xrange(0,bshape[0]):
            buffer[index] = float(lines[index][0])
    else:
        buffer = numpy.zeros(bshape,'d')

        for index in xrange(0,buffer.shape[0]):
            buffer[index,:] = map(float,lines[index])

    return buffer

# how is the record order going to be?

def writebinary(filename,a,append=0):
    if append:
        file = open(filename, 'a')
    else:
        file = open(filename, 'w')

    file.write(a.tostring())

    file.close()

# in a future version, this should also get the length from the file

def readbinary(filename,length):
    file = open(filename,'r')
    buffer = numpy.fromstring(file.read(length*8),'double')
    file.close()
    # then reshape the buffer if needed
    return buffer

# S/N utility functions

# this is a slow replacement for arrayfns.interp, which is not present
# in numpy

def linearinterpolate(y,x,xi):
    yi = xi.copy()

    j = 0

    for i in range(0,len(xi)):
        while j + 2 < len(x) and xi[i] > x[j+1]:
            j = j + 1

        yi[i] = y[j] + (y[j+1] - y[j]) * (xi[i] - x[j]) / (x[j+1] - x[j])

    return yi

# import arrayfns

def sn(signal,noise,stime,npatches):
    """Compute the optimal S/N for signal, sampled at intervals of
    stime, and for the total duration represented in the array, against noise
    represented by the time series noise; npatches overlapping periods are used
    to estimate the PSD of the noise."""

    # compute signal spectrum without windowing or averaging
    sspec = spect(signal,stime,0)

    # compute the noise spectrum, using segment averaging
    nspec = spect(noise,stime,npatches)

    # interpolate the noise to be defined on the same frequencies
    # of the signal's spectrum

    ispec = numpy.zeros(numpy.shape(sspec),dtype='d')
    ispec[:,0] = sspec[:,0]

    # ispec[:,1] = arrayfns.interp(nspec[:,1],nspec[:,0],ispec[:,0])
    ispec[:,1] = linearinterpolate(nspec[:,1],nspec[:,0],ispec[:,0])

    # the (S/N)^2 is given by 2T times the integrated ratio
    # of the spectral densities (the factor of 2 because the spectral
    # density is one-sided); notice however that the df is 1/T,
    # so we need only to sum up the array containing the ratio,
    # and multiply by two

    sratio = numpy.zeros(numpy.shape(sspec)[0],dtype='d')
    sratio[1:] = sspec[1:,1] / ispec[1:,1]
    sn2 = 2.0 * sum(sratio[1:])

    return math.sqrt(sn2)

import numpy.fft as FFT

def real(number):
	try:
		return number.real
	except AttributeError:
		return abs(number)

def noiseproduct(signal1,signal2,noise,stime,npatches):
    """Compute the noise inner product for signal1 and signal2, where
    the two signals are sampled at intervals stime, and the product is computed
    for the total duration represented in the array, against noise
    represented by the time series noise; npatches overlapping periods are used
    to estimate the PSD of the noise."""
    
    # compute signal FFT without windowing or averaging
    # this definition of the FFT satisfied Parseval's theorem, with
    # sum(signal**2) * stime == sum(abs(sfour)**2) / (stime*length(signal))
    # [since deltaf = 1 / (totaltime)]

    sfour1 = stime * FFT.fft(signal1)
    sfour2 = stime * FFT.fft(signal2)
    
    # compute the noise spectrum, using segment averaging,
    # and interpolate the noise to be defined on the same frequencies
    # of the signal's spectrum

    nspec = spect(noise,stime,npatches)

    siglen = len(signal1)
    deltaf = 1.0 / (stime * siglen)

    fourlen = siglen/2 + 1
    ispec = numpy.zeros([fourlen,2],dtype='d')

    ispec[:,0] = deltaf * numpy.arange(0,fourlen)
    # ispec[:,1] = arrayfns.interp(nspec[:,1],nspec[:,0],ispec[:,0])
    ispec[:,1] = linearinterpolate(nspec[:,1],nspec[:,0],ispec[:,0])
   
    return 4.0 * real(sum(sfour1[1:fourlen] * conjugate(sfour2[1:fourlen]) / ispec[1:,1])) * deltaf
	
# lisa positions from Ted Sweetser's file

import os
import os.path

datadir = os.path.join(os.path.dirname(__file__),'data')

def getLISApositions(filename):
    """Get LISA positions in seconds from an ASCII file, and returns four
    numpy arrays, corresponding to times [in seconds] and SSB coordinates
    of the LISA spacecraft [also in seconds]. 
    
    The file should be formatted as follows: column 1 gives the Julian date
    in seconds; columns 2-4, 5-7, 8-10 give the three SSB coordinates for each
    spacecraft, given as kms in an Earth Mean Equator, J2000 Equinox frame.
    It's assumed that timestamps are equally spaced."""
    
    try:
        pos = numpy.genfromtxt(filename)
    except IOError:
        pos = numpy.genfromtxt(os.path.join(datadir,filename))
    
    t = pos[:,0].copy()
    p1, p2, p3 = pos[:,1:4].copy(), pos[:,4:7].copy(), pos[:,7:10].copy()
    
    t -= t[0]
    secondsperday = 86400
    t *= secondsperday
    
    speedoflight = 299792.458
    p1 /= speedoflight; p2 /= speedoflight; p3 /= speedoflight
    
    return t,p1,p2,p3


def stdLISApositions():
    """Returns four numpy arrays corresponding to times [in seconds] along
    ten years and to the corresponding positions [the three SSB coordinates
    in Earth Mean Equator and J2000 Equinox, also given in seconds] of the
    three LISA spacecraft according to simulations run by JPL's Ted Sweetser
    on 2005-07-23 (from Excel spreadsheet states_baseline2.xls). The times
    are spaced by 1 day (86400 seconds), and they begin from zero instead of
    Ted's Julian date 2457023.5."""
    
    return getLISApositions(os.path.join(datadir,'positions.txt'))
    

def makeSampledLISA(filename,interp=2):
    """Returns an interpolated and cached SampledLISA object based on the position
    arrays returned by getLISApositions(); the argument interp sets the semilength
    of the interpolation window, and prebuffering is set to interp times the spacing
    of the stdLISApositions() data."""
    
    [t,p1,p2,p3] = getLISApositions(filename)
    
    dt = t[1] - t[0]
    slisa = lisaswig.SampledLISA(p1,p2,p3,dt,dt*interp,interp)
    
    return lisaswig.CacheLengthLISA(slisa,len(t),dt,interp)


def stdSampledLISA(interp=1):
    """Calls makeSampledLISA with the standard file given by stdLISApositions()."""
    
    return makeSampledLISA(os.path.join(datadir,'positions.txt'),interp)

