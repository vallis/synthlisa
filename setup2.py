#!/usr/bin/env python

from __future__ import print_function

import sys, os, platform

from setuptools import setup
from setuptools import Extension
import distutils.sysconfig

from Cython.Build import cythonize

import numpy

setup(name = 'synthlisa',
      version = '2.1.0', # remember to change it in __init__.py.in
      description = 'Synthetic LISA Simulator',

      author = 'Michele Vallisneri',
      author_email = 'vallis@vallis.org',
      url = 'https://github.com/vallis/synthlisa',

      packages = ['synthlisa'],
      package_dir = {'synthlisa': 'lisasim'},

      py_modules = [],

      ext_modules = cythonize(Extension('synthlisa.synthlisa',['lisasim/synthlisa.pyx'],
                                        language = "c++",
                                        include_dirs = [numpy.get_include()],
                                        # libraries = [],
                                        # library_dirs = [],
                                        # extra_link_args = linkArgs,
                                        extra_compile_args = ["-Wno-unused-function"]
                                        ))
      )
