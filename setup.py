#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc, get_python_lib
from distutils.dep_util import newer_group
from distutils.dir_util import mkpath
from distutils.util import get_platform
from distutils.spawn import spawn
from distutils.command import install_lib
from distutils.file_util import copy_file

import sys
import os
import glob
import re

def my_get_python_lib(prefixdir=None):
    if prefixdir and sys.platform == 'darwin':
        # this hack needed because of a bug in OS X Leopard's stock Python 2.5.1
        return get_python_lib(standard_lib=True,prefix=prefixdir) + '/site-packages'
    else:
        return get_python_lib(prefix=prefixdir)


versiontag = '1.3.7'

synthlisa_prefix = ''
numpy_prefix = ''
swig_bin = 'swig'
gsl_prefix = ''
make_clib = False

# At the moment, this setup script does not deal with --home.
# I should also modify the --help text to discuss these options

argv_replace = []

for arg in sys.argv:
    if arg.startswith('--prefix='):
        synthlisa_prefix = arg.split('=', 1)[1]
        argv_replace.append(arg)
    elif arg.startswith('--with-numpy='):
        numpy_prefix = arg.split('=', 1)[1]
    elif arg.startswith('--with-gsl='):
        gsl_prefix = arg.split('=', 1)[1]
    elif arg.startswith('--with-swig='):
        swig_bin = arg.split('=', 1)[1]
    elif arg.startswith('--make-clib'):
        make_clib = True
    else:
        argv_replace.append(arg)

sys.argv = argv_replace

gsl_cfiles = glob.glob('lisasim/GSL/*.c')
gsl_hfiles = glob.glob('lisasim/GSL/*.h')

lisasim_cppfiles = glob.glob('lisasim/*.cpp')
lisasim_hppfiles = glob.glob('lisasim/*.h')
lisasim_pyfiles  = glob.glob('lisasim/*.py')

# remove lisasim/lisasim-swig_wrap.h from headers

lisasim_hppfiles = filter(lambda s: s != 'lisasim/lisasim-swig_wrap.h',
                          lisasim_hppfiles)

source_files = lisasim_cppfiles + gsl_cfiles
header_files = lisasim_hppfiles + gsl_hfiles

# compile Id catalog for lisasim

idcatalog = ""

signature = re.compile('.*\$(Id.*)\$.*')

for file in source_files + header_files + lisasim_pyfiles:
    handle = open(file)

    for i in range(0,5):
        line = handle.readline()

        if signature.match(line):
            if idcatalog:
                idcatalog += '\n' + signature.match(line).group(1)
            else:
                idcatalog = signature.match(line).group(1)

    handle.close()

version_py = open('lisasim/version.py','w')
print >> version_py, "version_full = \"\"\"%s\"\"\"\n" % idcatalog
print >> version_py, "version_short = \"%s\"\n" % versiontag
version_py.close()

# Swig away!

lisasim_isource = 'lisasim/lisasim-swig.i'
lisasim_iheader = 'lisasim/lisasim-typemaps.i'

lisasim_cppfile = 'lisasim/lisasim-swig_wrap.cpp'
lisasim_pyfile = 'lisasim/lisaswig.py'

# Should I clean the swig-generated files when setup.py clean is issued?
# Better not, since it could engender problems if swig is not available.

def runswig(source,cppfile,pyfile,deps,cpp=1):
    if not os.path.isfile(cppfile) or not os.path.isfile(pyfile) \
           or newer_group(deps,cppfile) or newer_group(deps,pyfile):
        try:
            if cpp:
                spawn([swig_bin,'-w402','-c++','-python','-o',cppfile,source])
            else:
                spawn([swig_bin,'-w402','-python','-o',cppfile,source])
        except:
            print 'Sorry, I am unable to swig the modified ' + lisasim_isource
            sys.exit(1)

runswig(lisasim_isource,lisasim_cppfile,lisasim_pyfile,
        header_files + [lisasim_isource,lisasim_iheader])

if lisasim_cppfile not in source_files:
    source_files.append(lisasim_cppfile)

# Create the setdir.sh and setdir.csh scripts; they will be recreated
# each time setup.py is run, but it does not matter.

scripttempdir = 'build/temp.' + get_platform() + '-%s.%s' % sys.version_info[0:2]
mkpath(scripttempdir)

setdir_sh     = open(scripttempdir + '/synthlisa-setdir.sh','w')
setdir_csh    = open(scripttempdir + '/synthlisa-setdir.csh','w')
recompile_sh  = open(scripttempdir + '/synthlisa-recompile.sh','w')

pythonpath = ''
installpath = sys.exec_prefix

if synthlisa_prefix:
    pythonpath = get_python_lib(prefix=synthlisa_prefix)
    installpath = synthlisa_prefix

# not needed with numpy, which can find the module from site-packages

# if numpy_prefix:
#     if pythonpath:
#         pythonpath = pythonpath + ':'
#
#     pythonpath = pythonpath + get_python_lib(prefix=numpy_prefix) + '/numpy'

mpi_prefix = numpy_prefix

if mpi_prefix:
    if pythonpath:
        pythonpath = pythonpath + ':'
        
    pythonpath = pythonpath + get_python_lib(prefix=mpi_prefix) + '/mpi'

print >> setdir_sh, """if [ -z "${PYTHONPATH}" ]
then
    PYTHONPATH="%s"; export PYTHONPATH
else
    PYTHONPATH="%s:$PYTHONPATH"; export PYTHONPATH
fi
""" % (pythonpath, pythonpath)

print >> setdir_csh, """if !($?PYTHONPATH) then
    setenv PYTHONPATH %s
else
    setenv PYTHONPATH %s:$PYTHONPATH
endif
""" % (pythonpath, pythonpath)

print >> recompile_sh, """#!/bin/sh
pushd %s
python setup.py install --prefix=%s $*
popd""" % (os.getcwd(),installpath)

recompile_sh.close()
setdir_csh.close()
setdir_sh.close()

# Ready to setup, build, install!

# this used to get the location of the Numeric C headers

# if numpy_prefix:
#     numpy_hfiles = get_python_inc(prefix=numpy_prefix)
# else:
#     numpy_hfiles = get_python_inc()

# get the location of the numpy C headers

if numpy_prefix:
	numpy_hfiles = get_python_lib(prefix=numpy_prefix) + '/numpy/core/include'
else:
    try:
        import numpy
        numpy_hfiles = numpy.__path__[0] + '/core/include'
    except ImportError:
	    numpy_hfiles = get_python_lib() + '/numpy/core/include'

if pythonpath:
    setdir_scripts = [scripttempdir + '/synthlisa-setdir.sh',
                      scripttempdir + '/synthlisa-setdir.csh',
                      scripttempdir + '/synthlisa-recompile.sh']
else:
    setdir_scripts = []

synthlisapackages = ['synthlisa']

synthlisapackage_dir = {'synthlisa' : 'lisasim'}

# do contribs (CPP only for the moment...)

contribs = []

# Duncan's recipe to install static c libraries

class qm_install_lib(install_lib.install_lib):
    def run(self):
        install_lib.install_lib.run(self)
        
        if self.distribution.has_c_libraries():
            build_clib = self.get_finalized_command('build_clib')
            libs = build_clib.get_library_names()

            clib_dir = build_clib.build_clib

            for lib in libs:
                clib = 'lib' + lib + '.a'
            
                src_file = os.path.join(clib_dir, clib)
                dest_file = os.path.join(self.install_dir, clib)
                    
                copy_file(src_file, dest_file)
                
                if sys.platform[:6] == "darwin":
                    spawn(['ranlib'] + [dest_file])

for entry in glob.glob('contrib/*'):
    if os.path.isdir(entry):
        contrib_packagename = os.path.basename(entry)

        contrib_source_files = glob.glob(entry + '/*.cpp') + glob.glob(entry + '/*.cc')
        contrib_header_files = glob.glob(entry + '/*.h') + glob.glob(entry + '/*.hh')

        contrib_swigfiles = glob.glob(entry + '/*.i')

        extensions = 0

        # each SWIG file creates a separate extension
        for contrib_swigfile in contrib_swigfiles:

            # let's check if this is a C or C++ extension
            iscpp = 1

            for line in open(contrib_swigfile).readlines():
                if 'is C extension' in line:
                    # do C extension
                    iscpp = 0
                    contrib_source_files = glob.glob(entry + '/*.c')

            # let's check if GSL is needed
            gsl_required = 0

            for line in open(contrib_swigfile).readlines():
                if 'requires GSL' in line:
                    # the magic words are "requires GSL"
                    gsl_required = 1

            if gsl_required == 1 and gsl_prefix == '':
                print "No GSL, skipping " + contrib_swigfile
                break

            contrib_basefile = re.sub('\.i$','',contrib_swigfile)
            contrib_basename = os.path.basename(contrib_basefile)
            
            # assume SWIG file has the same name of the module
            if iscpp == 1:
                contrib_wrapfile = contrib_basefile + '_wrap.cpp'
            else:
                contrib_wrapfile = contrib_basefile + '_wrap.c'
                
            contrib_pyfile = contrib_basefile + '.py'
            
            runswig(contrib_swigfile,contrib_wrapfile,contrib_pyfile,
                    contrib_header_files + [contrib_swigfile],iscpp)
    
            contrib_extname = contrib_packagename + '/_' + contrib_basename

            if not contrib_wrapfile in contrib_source_files: 
                contrib_source_files.append(contrib_wrapfile)

            if gsl_required == 1:
                contrib_extension = Extension(contrib_extname,
                                              contrib_source_files,
                                              include_dirs = [entry,numpy_hfiles,gsl_prefix + '/include'],
                                              library_dirs = [gsl_prefix + '/lib'],
                                              runtime_library_dirs = [gsl_prefix + '/lib'],
                                              libraries=['gsl', 'gslcblas'],
                                              depends = contrib_header_files)            
            else:
                contrib_extension = Extension(contrib_extname,
                                              contrib_source_files,
                                              include_dirs = [entry,numpy_hfiles],
                                              depends = contrib_header_files)
    
            contribs.append(contrib_extension)
            
            extensions = extensions + 1

        if extensions:
            synthlisapackages.append(contrib_packagename)
            synthlisapackage_dir[contrib_packagename] = entry

# if we're asked to make a static .a library, set that up, and remove the old already-built library
# which causes trouble to OS X's Universal Python...

if make_clib == True:
    clibrary = [('synthlisa',{'sources': filter(lambda s: '_wrap.cpp' not in s,source_files),
                              'depends': header_files,
                              'include_dirs': [numpy_hfiles, get_python_inc()]} )]

    # this hack needed on OS X for universal binaries since ar fails if the .a is already present...

    try:
        os.remove('build/temp.' + get_platform() + '-%s.%s' % sys.version_info[0:2] + '/libsynthlisa.a')
    except:
        pass
else:
    clibrary = []

# find all synthlisa include files for installation under lib/pythonx.x/site-packages/synthlisa

installincludes = [re.sub('lisasim/','',hfile) for hfile in header_files]

# do the actual setup

setup(name = 'synthLISA',
      version = versiontag,
      description = 'Synthetic LISA Simulator',
      long_description = idcatalog,

      author = 'Michele Vallisneri',
      author_email = 'vallis@vallis.org',
      url = 'http://www.vallis.org/syntheticlisa',

      packages = synthlisapackages,
      
      package_dir = synthlisapackage_dir,

      scripts = setdir_scripts,

      # TO DO these data files should be installed within the pythonx.x/site-packages/synthlisa directory
      package_data = {'synthlisa': ['data/positions.txt',
                                    'data/lisa-xml.dtd',
                                    'data/lisa-xml.xsl',
                                    'data/lisa-xml.css'] + installincludes},

      cmdclass = {'install_lib' : qm_install_lib},

      libraries = clibrary,

      ext_modules = [Extension('synthlisa/_lisaswig',
                               source_files,
                               include_dirs = [numpy_hfiles],
                               depends = header_files
                               )] + contribs
      )
