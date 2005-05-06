#!/bin/sh

NUMERIC=`ls Numeric-*.tar.gz | tail -1`
SWIG=`ls swig-*.tar.gz | tail -1`

tar zxvf $NUMERIC
tar zxvf $SWIG

NUMERICDIR=`ls | grep Numeric | head -1`
SWIGDIR=`ls | grep SWIG | head -1`
THISDIR=`pwd`

cd $NUMERICDIR

python setup.py install --prefix=$THISDIR

cd ..

cd $SWIGDIR

./configure --prefix=$THISDIR
make
make install

cd ..

python setup.py install --with-numeric=$THISDIR --with-swig=$THISDIR/bin/swig --prefix=$THISDIR
