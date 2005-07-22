#!/usr/bin/python

# test of noise distortion with interpolation schemes of increasing order

# this script demonstrates the generation of pseudorandom noise time
# series with interpolation schemes of increasing order

# import all the libraries that are needed

from synthlisa import *
from Numeric import *


# generate the same noise object with:

# no interpolation (noise defaults to nearest sample)
lasernoise0 = InterpolateNoise(1.0,256.0,1.0,0.0,0)

# linear interpolation
lasernoise1 = InterpolateNoise(1.0,256.0,1.0,0.0,1)

# Lagrange interpolation with order = 4 (semiwindow of 2)
lasernoise2 = InterpolateNoise(1.0,256.0,1.0,0.0,2)

# Lagrange interpolation with order = 8 (semiwindow of 4)
lasernoise4 = InterpolateNoise(1.0,256.0,1.0,0.0,4)

# Lagrange interpolation with order = 32 (semiwindow of 16)
lasernoise8 = InterpolateNoise(1.0,256.0,1.0,0.0,16)

# get "samples" values of the noises, at times separated by "stime"

samples = 2**18 # 2**20 takes 23s on a 1.25GHz G4
sstep = 0.1

patches = 512

# produce the noises

noise0 = getobs(samples,sstep,lasernoise0.noise)
noise1 = getobs(samples,sstep,lasernoise1.noise)
noise2 = getobs(samples,sstep,lasernoise2.noise)
noise4 = getobs(samples,sstep,lasernoise4.noise)
noise8 = getobs(samples,sstep,lasernoise8.noise)

# compute the spectra and write to disk

writearray('data/interpnoise0-freq.txt',spect(noise0,sstep,patches)[1:])
writearray('data/interpnoise1-freq.txt',spect(noise1,sstep,patches)[1:])
writearray('data/interpnoise2-freq.txt',spect(noise2,sstep,patches)[1:])
writearray('data/interpnoise4-freq.txt',spect(noise4,sstep,patches)[1:])
writearray('data/interpnoise8-freq.txt',spect(noise8,sstep,patches)[1:])
