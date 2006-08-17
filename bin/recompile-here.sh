#!/bin/sh

# $Id$

python setup.py install --with-numpy=`pwd` --with-swig=`pwd`/bin/swig --prefix=`pwd` $1
