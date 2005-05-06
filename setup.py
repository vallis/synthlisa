#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc, get_python_lib
from distutils.dep_util import newer_group
from distutils.spawn import spawn

import sys
import os
import glob

synthlisa_prefix = ''
numeric_prefix = ''
swig_bin = 'swig'

# At the moment, this setup script does not deal with --home.
# I should also modify the --help text to discuss these options

argv_replace = []

for arg in sys.argv:
    if arg.startswith('--prefix='):
        synthlisa_prefix = arg.split('=', 1)[1]
        argv_replace.append(arg)
    elif arg.startswith('--with-numeric='):
        numeric_prefix = arg.split('=', 1)[1]
    elif arg.startswith('--with-swig='):
        swig_bin = arg.split('=', 1)[1]
    else:
        argv_replace.append(arg)

sys.argv = argv_replace

gsl_cfiles = glob.glob('lisasim/GSL/*.c')
gsl_hfiles = glob.glob('lisasim/GSL/*.h')

lisasim_cppfiles = glob.glob('lisasim/*.cpp')
lisasim_hppfiles = glob.glob('lisasim/*.h')

source_files = lisasim_cppfiles + gsl_cfiles
header_files = lisasim_hppfiles + gsl_hfiles

# Swig away!

lisasim_isource = 'lisasim/lisasim-swig.i'
lisasim_iheader = 'lisasim/lisasim-typemaps.i'

lisasim_cppfile = 'lisasim/lisasim-swig_wrap.cpp'
lisasim_pyfile = 'lisasim/lisaswig.py'

# Should I clean the swig-generated files when setup.py clean is issued?
# Better not, since it could engender problems if swig is not available.

if not os.path.isfile(lisasim_cppfile) or not os.path.isfile(lisasim_pyfile) \
       or newer_group(header_files + [lisasim_isource,lisasim_iheader],
                      lisasim_cppfile) \
       or newer_group(header_files + [lisasim_isource,lisasim_iheader],
                      lisasim_pyfile):
    try:
        spawn([swig_bin,
               '-w402','-c++','-python',
               '-o',lisasim_cppfile,
               lisasim_isource])
    except:
        print 'Sorry, I am unable to swig the modified ' + lisasim_isource
        sys.exit(cmd)

if lisasim_cppfile not in source_files:
    source_files.append(lisasim_cppfile)

# Create the setdir.sh and setdir.csh scripts; they will be recreated
# each time setup.py is run, but it does not matter.

setdir_sh = open('lisasim/synthlisa-setdir.sh','w')
setdir_csh = open('lisasim/synthlisa-setdir.csh','w')

pythonpath = ''

if synthlisa_prefix:
    pythonpath = get_python_lib(prefix=synthlisa_prefix)

if numeric_prefix:
    if pythonpath:
        pythonpath = pythonpath + ':'

    pythonpath = pythonpath + get_python_lib(prefix=numeric_prefix) + '/Numeric'

print >> setdir_sh, """
if [ -z "${PYTHONPATH}" ]
then
    PYTHONPATH="%s"; export PYTHONPATH
else
    PYTHONPATH="%s:$PYTHONPATH"; export PYTHONPATH
fi
""" % (pythonpath, pythonpath)

print >> setdir_csh, """
if !($?PYTHONPATH) then
    setenv PYTHONPATH %s
else
    setenv PYTHONPATH %s:$PYTHONPATH
endif
""" % (pythonpath, pythonpath)

setdir_csh.close()
setdir_sh.close()

# Ready to setup, build, install!

if numeric_prefix:
    numeric_hfiles = get_python_inc(prefix=numeric_prefix)
else:
    numeric_hfiles = get_python_inc()

if pythonpath:
    setdir_scripts = ['lisasim/synthlisa-setdir.sh','lisasim/synthlisa-setdir.csh']
else:
    setdir_scripts = []

setup(name = 'synthLISA',
      version = '1.2 (2005/05/06)',
      description = 'Synthetic LISA Simulator',

      author = 'Michele Vallisneri',
      author_email = 'vallis@vallis.org',
      url = 'http://www.vallis.org/syntheticlisa',
      
      package_dir = {'' : 'lisasim'},
      py_modules = ['lisautils','lisaswig'],

      scripts = setdir_scripts,

      ext_modules = [Extension('_lisaswig',
                               source_files,
                               include_dirs = [numeric_hfiles],
                               depends = header_files
                               )]
     )
