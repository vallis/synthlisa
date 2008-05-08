#!/bin/sh

# $Id$

export COPY_EXTENDED_ATTRIBUTES_DISABLE=true
python setup.py sdist --with-swig=`pwd`/bin/swig --with-numpy=`pwd`
