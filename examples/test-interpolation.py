#!/usr/bin/env python

# $Id$
# $Date$
# $Author$
# $Revision$

# test of noise distortion with interpolation schemes of increasing order

# this script demonstrates the generation of pseudorandom noise time
# series with interpolation schemes of increasing order

# import all the libraries that are needed

from synthlisa import *
from Numeric import *

# generate the same noise object with:

# no interpolation (noise defaults to nearest sample)
lasernoise0 = PowerLawNoise(1.0,256.0,1.0,0.0,0)

# linear interpolation
lasernoise1 = PowerLawNoise(1.0,256.0,1.0,0.0,1)

# Lagrange interpolation with order = 4 (semiwindow of 2)
lasernoise2 = PowerLawNoise(1.0,256.0,1.0,0.0,2)

# Lagrange interpolation with order = 8 (semiwindow of 4)
lasernoise4 = PowerLawNoise(1.0,256.0,1.0,0.0,4)

# Lagrange interpolation with order = 32 (semiwindow of 16)
lasernoise8 = PowerLawNoise(1.0,256.0,1.0,0.0,16)

# get "samples" values of the noises, at times separated by "stime"

samples = 2**18 # 2**20 takes 23s on a 1.25GHz G4
sstep = 0.1

patches = 512

# produce the noises

noise0 = getobs(samples,sstep,lasernoise0)
noise1 = getobs(samples,sstep,lasernoise1)
noise2 = getobs(samples,sstep,lasernoise2)
noise4 = getobs(samples,sstep,lasernoise4)
noise8 = getobs(samples,sstep,lasernoise8)

# compute the spectra and write to disk

spec0 = spect(noise0,sstep,patches)
spec1 = spect(noise1,sstep,patches)
spec2 = spect(noise2,sstep,patches)
spec4 = spect(noise4,sstep,patches)
spec8 = spect(noise8,sstep,patches)

allspec = transpose([spec0[:,0],spec0[:,1],
                                spec1[:,1],
                                spec2[:,1],
                                spec4[:,1],
                                spec8[:,1]])

outputXML = lisaXML('data/interpolation',
                    author='Michele Vallisneri',
                    comments='Test of interpolated white noise')
                    
outputXML.TDISpectraSelfDescribed(allspec[1:],'f,i0,i1,i2,i4,i8',encoding='Text')

outputXML.close()

print "You can plot the results of this script by running"
print "  ./plotxml.py data/interpolation.xml eps/interpolation-spectra.eps"
print "or"
print "  ./plotxml.py data/interpolation.xml eps/interpolation-spectra.pdf"
