#!/usr/bin/python

# test of the Synthetic LISA TDI output for a monochromatic binary
# against the output of the LISA Simulator (www.physics.montana.edu/lisa)
# note: this test needs the data file "tdibinary-X-circ.bin",
# to be downloaded separately

# it also needs "test-binary.py" to have been run with samples = 2**21
# and stime = 16

# import all the libraries that are needed

from lisaswig import *
from Numeric import *
from lisautils import *

# read the product of test-binary.py

signalXcirc = readbinary('data/tdibinary-X-circ.bin',2**21)

# we have stored the LISA Simulator signal in the examples directory

signalXmont = readbinary('tdibinary-X-montana.bin',2**20)

# set the sampling rates

ssrate = 16
msrate = 3.15581498e7/1048576

# use only a year of synthetic LISA data

yearlen = 365 * 3600 * 24 / ssrate
signalXcirc = signalXcirc[0:yearlen]

# take a clean spectrum of my data (no windowing or averaging)

sspec = spect(signalXcirc,ssrate,0)

# convert the montana spectrum to my units...

# first take a clean spectrum (no windowing or averaging)

mspec = spect(signalXmont,msrate,0)

# multiply by (10^10 m)^2

mspec[:,1] = mspec[:,1] * (1.0e10)**2

# multiply by (2 pi f)^2

mspec[:,1] = mspec[:,1] * (2.0 * math.pi * mspec[:,0])**2

# divide by c^2

mspec[:,1] = mspec[:,1] / (2.99792458e8)**2

# write spectra to disk
# there are roughly twice as many data points for Synthetic LISA
# as for the LISA Simulator

writearray('data/tdibinary-X-circ-closeup.txt',sspec[61208:61409])
writearray('data/tdibinary-X-montana-closeup.txt',mspec[61301:61402])
