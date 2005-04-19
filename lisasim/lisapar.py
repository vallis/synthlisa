import pypar

from lisaswig import *
from Numeric import *

from time import time

# triad must be given as [sourcefn, parameters, observables]

def getobsp(snum,stime,triad,zerotime=0.0,debug=0):
    inittime = int(time())

    (srcfunc, parameters, obs) = triad

    if type(obs) == list or type(obs) == tuple:
        multobs = True
    else:
        multobs = False

    # I seem to have trouble using broadcast...
    # well, at least I could avoid sending everything to everybody

    if pypar.rank() == 0:
        if debug != 0:
            print "Preparing for parallel execution...",

        for i in range(1,pypar.size()):
            pypar.send(parameters,i)
    else:
        parameters = pypar.receive(0)

    if multobs:
        obslen = len(obs)
        array = zeros((snum,obslen),typecode='d')
    else:
        array = zeros(snum,typecode='d')

    sourcesdone = 0

    for mysrc in range(pypar.rank(),len(parameters),pypar.size()):
        currenttime = int(time())
        
        if debug != 0:
            print "   ...running source ",mysrc," (",parameters[mysrc][0],") on CPU ",pypar.rank()

        if multobs:
            [lisa, source, tdi] = srcfunc(parameters[mysrc])

            for i in arange(0,snum):
                for j in range(0,obslen):
                    array[i,j] += obs[j](tdi,zerotime+i*stime)
        else:
            [lisa, source, tdi] = srcfunc(parameters[mysrc])

            for i in arange(0,snum):
                array[i] += obs(tdi,zerotime+i*stime)

        sourcesdone += 1

    currenttime = int(time()) - inittime
    print "CPU %d: %d sources, %d s [%d (multi)samples/s]" % (pypar.rank(),sourcesdone,
                                                              currenttime,int(sourcesdone*snum/currenttime))
    
    if pypar.rank() == 0:
        if debug != 0:
            print "   ...root ready to collect results"

        sumarray = array.copy()

        for i in range(1,pypar.size()):
            array,status = pypar.receive(pypar.any_source,return_status=True)

            if debug != 0:
                print "   ...root received results from ",status.source

            sumarray += array

        currenttime = int(time()) - inittime
        print "Total elapsed time for %d sources, %d s [%d (multi)samples/s]" % (len(parameters),
                                                                                 currenttime,int(sourcesdone*snum/currenttime))
        
        return sumarray
    else:
        if debug != 0:
            print "   ...CPU ",pypar.rank()," sending results to root"
        
        pypar.send(array,0)

        if multobs:
            return zeros((1,obslen),typecode='d')
        else:
            return zeros(1,typecode='d')
