#!/usr/bin/env python

# test of noise distortion with interpolation schemes of increasing order

# this script demonstrates the generation of pseudorandom noise time
# series with interpolation schemes of increasing order

# import all the libraries that are needed

from synthlisa import *

# generate the same noise object with:

# no interpolation (noise defaults to nearest sample)
lasernoise0 = PowerLawNoise(1.0,256.0,1.0,0.0,interplen=0)

# linear interpolation
lasernoise1 = PowerLawNoise(1.0,256.0,1.0,0.0,interplen=1)

# Lagrange interpolation with order = 4 (semiwindow of 2)
lasernoise2 = PowerLawNoise(1.0,256.0,1.0,0.0,interplen=2)

# Lagrange interpolation with order = 8 (semiwindow of 4)
lasernoise4 = PowerLawNoise(1.0,256.0,1.0,0.0,interplen=4)

# Lagrange interpolation with order = 32 (semiwindow of 16)
lasernoise8 = PowerLawNoise(1.0,256.0,1.0,0.0,interplen=16)

# get "samples" values of the noises, at times separated by "stime"

samples = 2**18 # 2**20 takes 23s on a 1.25GHz G4
sstep = 0.1

patches = 512

# produce the noises

noise0 = getobsc(samples,sstep,lasernoise0)
noise1 = getobsc(samples,sstep,lasernoise1)
noise2 = getobsc(samples,sstep,lasernoise2)
noise4 = getobsc(samples,sstep,lasernoise4)
noise8 = getobsc(samples,sstep,lasernoise8)

# compute the spectra and write to disk

writearray('data/interpnoise0-freq.txt',spect(noise0,sstep,patches)[1:])
writearray('data/interpnoise1-freq.txt',spect(noise1,sstep,patches)[1:])
writearray('data/interpnoise2-freq.txt',spect(noise2,sstep,patches)[1:])
writearray('data/interpnoise4-freq.txt',spect(noise4,sstep,patches)[1:])
writearray('data/interpnoise8-freq.txt',spect(noise8,sstep,patches)[1:])
