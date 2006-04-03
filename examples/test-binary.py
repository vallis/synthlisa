#!/usr/bin/env python

from synthlisa import *
from Numeric import *

import math

# this function returns a PhysicalBinary object of given frequency, 

def PhysicalBinary(freq,phi0,inc,amp,elat,elon,psi):
    """Returns a SimpleBinary Wave object of frequency freq [Hz],
    initial phase phi0 [0,2pi], inclination inc [0,2pi], ecliptic latitude
    and longitude elat [-pi/2,pi/2] and elon [0,2pi], and polarization
    angle psi [0,pi]. The amplitude is compute from the two masses [given
    in Msun] and the distance [given in Kpc]; these three parameters are
    passed as the tuple amp = (m1,m2,dist)."""

    G = 6.67259e-11
    c = 2.99792458e8

    kpc = 3.0856675807e19
    Msun = 1.9889e30

    M1 = amp[0]*Msun
    M2 = amp[1]*Msun

    r = amp[2]*kpc

    R = (G*(M1+M2)*(1.0/(2.0*math.pi*freq/2))**2)**(1.0/3.0)
    ampl = (2.0*M1*M2/(r*R) * (G/(c*c))**2)

    return SimpleBinary(freq, phi0, inc, ampl, elat, elon, psi)


# creates a standard LISA geometry object

lisa = EccentricInclined()

# creates the SimpleBinary object

mysystem = PhysicalBinary(freq = 1.0e-3,
                          phi0 = 0.0,
                          inc = math.pi/4.2,
                          amp = (0.5,0.033,0.05),
                          elat = math.pi/3.0,
                          elon = 0.0, 
                          psi = math.pi/10.0)
               
# creates a TDIsignal object with lisa and mysystem
               
signalTDI = TDIsignal(lisa,mysystem)

stime = 32.0

# creates a TDInoise object with standard proof-mass and shot noise
# no laser noise since it would not be canceled with first-generation TDI

noiseTDI =  TDInoise(lisa,stime,2.5e-48,
                          stime,1.8e-37,
                          stime,0.0)

# get one year's worth of samples

samples = int(2**25/stime)
patches = 256
 
[alphas,betas,gammas] = transpose(getobsc(samples,stime,[signalTDI.alpha,
                                                         signalTDI.beta,
                                                         signalTDI.gamma]))

[alphan,betan,gamman] = transpose(getobsc(samples,stime,[noiseTDI.alpha,
                                                         noiseTDI.beta,
                                                         noiseTDI.gamma]))

# form the optimal TDI observables

signalA = (gammas - alphas) / math.sqrt(2.0)
signalE = (alphas - 2.0*betas + gammas) / math.sqrt(6.0)
signalT = (alphas + betas + gammas) / math.sqrt(3.0)

noiseA = (gamman - alphan) / math.sqrt(2.0)
noiseE = (alphan - 2.0*betan + gamman) / math.sqrt(6.0)
noiseT = (alphan + betan + gamman) / math.sqrt(3.0)

# compute the optimal signal to noise (for short time series, this is somewhat
# sensitive to the fluctuations in the noise, which is used to compute S_n)

print "S/N = ", math.sqrt(sn(signalA,noiseA,stime,patches)**2 + 
                          sn(signalE,noiseE,stime,patches)**2 + 
                          sn(signalT,noiseT,stime,patches)**2), "over T = ", stime * samples, " s."

# output signal+noise time series and spectra to disk

outputXML = lisaXML('data/binary',
                    author='Michele Vallisneri',
                    comments='EccentricInclined 1st-gen. optimal-TDI observables for 1 mHz binary')

outputXML.TDIData(transpose([signalA + noiseA,
                             signalE + noiseE,
                             signalT + noiseT]),samples,stime,'AOf,EOf,TOf')

specA = spect(signalA + noiseA,stime,patches)
specE = spect(signalE + noiseE,stime,patches)
specT = spect(signalT + noiseT,stime,patches)

outputXML.TDISpectraSelfDescribed(transpose([specA[1:,0],specA[1:,1],
                                                         specE[1:,1],
                                                         specT[1:,1]]),'f,AOf,EOf,TOf',encoding='Text')

outputXML.close()

# plot the spectra

try:
    import pyx    
    import os
    
    os.system('./plotxml.py ' + 'data/binary.xml' + ' ' + 'eps/binary.pdf')
    print "Spectra plotted in eps/binary.pdf"
except:
    print "Could not output PyX graphics"
