#!/usr/bin/env python

import sys
import os
import os.path
import glob
import operator
import re

def escapespace(dirname):
	return re.sub(' ','\ ',dirname)

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

# if we're given a command-line argument, try to parse it to decide what to
# install; otherwise install everything that we find

if len(sys.argv) > 1:
	donumeric,doswig,dorxp,dopyx,dosynthlisa,dommpi = 0,0,0,0,0,0

	for arg in sys.argv:
		if arg == 'Numeric' or arg == 'numeric':
			donumeric = 1
		elif arg == 'swig' or arg == 'SWIG':
			doswig = 1
		elif arg == 'rxp' or arg == 'RXP' or arg == 'pyRXP':
			dorxp = 1
		elif arg == 'pyx' or arg == 'PyX':
			dopyx = 1
		elif arg == 'synthlisa' or arg == 'synthLISA':
			dosynthlisa = 1
		elif arg == 'mmpi' or arg == 'mpi' or arg == 'MMPI' or arg == 'MPI':
			dommpi = 1
else:
	# don't do MMPI by default
	donumeric,doswig,dorxp,dopyx,dosynthlisa,dommpi = 1,1,1,1,1,0

thisdir = os.getcwd()

os.chdir('packages')

pckgdir = os.getcwd()

if donumeric:
	numericdir = findpackage('Numeric','Numeric')
if doswig:
	swigdir    = findpackage('swig','SWIG')
if dorxp:
	rxpdir     = findpackage('pyRXP','pyRXP')
if dopyx:
	pyxdir     = findpackage('PyX','PyX')
if dommpi:
	mmpidir    = findpackage('mikempi','mikempi')

if donumeric and numericdir:
    try:
        os.chdir(numericdir)    
        assert os.system('python setup.py install --prefix=' + escapespace(thisdir)) == 0
        os.chdir(pckgdir)
    except AssertionError:
        print "Error while compiling and installing Numeric"
        raise

if doswig and swigdir:
    try:
        os.chdir(swigdir)
        assert os.system('./configure --prefix=' + escapespace(thisdir)) == 0
        assert os.system('make') == 0
        assert os.system('make install') == 0
        os.chdir(pckgdir)
    except AssertionError:
        print "Error while compiling and installing SWIG"
        raise

if dopyx and pyxdir:
    try:
        os.chdir(pyxdir)
        assert os.system('python setup.py install --prefix=. --root=' + escapespace(thisdir)) == 0
        os.chdir(pckgdir)
    except AssertionError:
        print "Error while compiling and installing PyX"
        raise
        
if dorxp and rxpdir:
    try:
        os.chdir(rxpdir + '/pyRXP')
        assert os.system('python setup.py install --prefix=' + escapespace(thisdir)) == 0
        os.chdir(pckgdir)
    except AssertionError:
        print "Error while compiling and installing pyRXP"
        raise

if dommpi and mmpidir:
	try:
		assert sys.version_info[0] + 0.1*sys.version_info[1] >= 2.4, "Need at least Python 2.4 to compile MMPI!"
		assert os.system('mpicc --version > /dev/null') == 0, "Need to have mpicc in your path!"
		
		os.chdir(mmpidir)
		assert os.system('python setup.py install --prefix=' + escapespace(thisdir) + ' --with-numeric=' + escapespace(thisdir)) == 0 
		os.chdir(pckgdir)
	except AssertionError:
		print "Error while compiling and installing MMPI"
		raise

os.chdir(thisdir)
        
try:
    if dosynthlisa:
        assert os.system(('python setup.py install --with-numeric=' + escapespace(thisdir) +
                          ' --with-swig=' + escapespace(thisdir) + '/bin/swig --prefix=' + escapespace(thisdir))) == 0

        print """Now (and before every synthLISA session) you should run the command

source bin/synthlisa-setdir.sh  [if you use bash or]
source bin/synthlisa-setdir.csh [if you use t/csh].

Even better, add these commands to your .profile (for bash) or
.cshrc/.tcshrc (for t/csh)."""
                 
except AssertionError:
    print "Error while compiling and installing synthLISA"
    raise
