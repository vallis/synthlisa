#!/usr/bin/python

import os
import os.path
import glob
import operator

def findpackage(packagename,dirname):
    packagetargz = glob.glob(packagename + '-*.tar.gz')
    packagetgz   = glob.glob(packagename + '-*.tgz')
    packagetar   = glob.glob(packagename + '-*.tar')    
    packagezip   = glob.glob(packagename + '-*.zip')

    if packagetargz:
        print "Unpacking " + packagetargz[-1]
        os.system('tar zxvf ' + packagetargz[-1] + ' > /dev/null')
    elif packagetgz:
        print "Unpacking " + packagetgz[-1]
        os.system('tar zxvf ' + packagetgz[-1] + ' > /dev/null')
    elif packagetar:
        print "Unpacking " + packagetar[-1]
        os.system('tar xvf ' + packagetar[-1] + ' > /dev/null')
    elif packagezip:
        print "Unpacking " + packagezip[-1]
        os.system('unzip ' + packagezip[-1] + ' > /dev/null')

    dir = filter(os.path.isdir,glob.glob(dirname + '-*') + glob.glob(packagename + '-*'))

    print "Using unpacking dir " + dir[-1]
    return dir[-1]

thisdir = os.getcwd()

os.chdir('packages')

pckgdir = os.getcwd()

numericdir = findpackage('Numeric','Numeric')
swigdir    = findpackage('swig','SWIG')
rxpdir     = findpackage('pyRXP','pyRXP')
pyxdir     = findpackage('PyX','PyX')

if numericdir:
    try:
        os.chdir(numericdir)    
        assert os.system('python setup.py install --prefix=' + thisdir) == 0
        os.chdir(pckgdir)
    except AssertionError:
        print "Error while compiling and installing Numeric"
        raise

if swigdir:
    try:
        os.chdir(swigdir)
        assert os.system('./configure --prefix=' + thisdir) == 0
        assert os.system('make') == 0
        assert os.system('make install') == 0
        os.chdir(pckgdir)
    except AssertionError:
        print "Error while compiling and installing SWIG"
        raise

if pyxdir:
    try:
        os.chdir(pyxdir)
        assert os.system('python setup.py install --prefix=. --root=' + thisdir) == 0
        os.chdir(pckgdir)
    except AssertionError:
        print "Error while compiling and installing PyX"
        raise
        
if rxpdir:
    try:
        os.chdir(rxpdir + '/pyRXP')
        assert os.system('python setup.py install --prefix=' + thisdir) == 0
        os.chdir(pckgdir)
    except AssertionError:
        print "Error while compiling and installing pyRXP"
        raise

os.chdir(thisdir)
        
try:
    assert os.system(('python setup.py install --with-numeric=' + thisdir +
                      ' --with-swig=' + thisdir + '/bin/swig --prefix=' + thisdir)) == 0
except AssertionError:
    print "Error while compiling and installing synthLISA"
    raise
