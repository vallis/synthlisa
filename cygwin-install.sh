#!/bin/sh

NUMERIC=`ls Numeric-*.tar.gz | tail -1`

tar zxvf $NUMERIC

NUMERICDIR=`ls | grep Numeric | head -1`
THISDIR=`pwd`

cd $NUMERICDIR

python setup.py install --prefix=$THISDIR

cd ..

python setup.py install --with-numeric=$THISDIR --prefix=$THISDIR
