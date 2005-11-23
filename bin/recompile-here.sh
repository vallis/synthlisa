#!/bin/sh

# $Id$

python setup.py install --with-numeric=`pwd` --with-swig=`pwd`/bin/swig --prefix=`pwd`
