#!/usr/bin/env python

from __future__ import print_function

import sys, os, platform

from setuptools import setup
from setuptools import Extension
import distutils.sysconfig

from Cython.Build import cythonize

import numpy

lisasources = ['lisasim/lisasim-lisa.cpp',
               'lisasim/lisasim-signal.cpp',
               'lisasim/lisasim-tens.cpp',
               'lisasim/lisasim-wave.cpp']

setup(name = 'synthlisa',
      version = '2.1.0', # remember to change it in __init__.py.in
      description = 'Synthetic LISA Simulator',

      author = 'Michele Vallisneri',
      author_email = 'vallis@vallis.org',
      url = 'https://github.com/vallis/synthlisa',

      packages = ['synthlisa'],
      package_dir = {'synthlisa': 'lisasim'},

      py_modules = [],

      ext_modules = cythonize(Extension('synthlisa.synthlisa',
                                        ['lisasim/synthlisa.pyx'] + lisasources,
                                        language = "c++",
                                        include_dirs = [numpy.get_include()],
                                        libraries = ['gsl', 'gslcblas'],
                                        # library_dirs = ['/usr/local'],
                                        # extra_link_args = [],
                                        extra_compile_args = ["-Wno-unused-function"]
                                        ))
      )
