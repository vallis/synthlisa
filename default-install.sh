#!/bin/sh

# $Id$

NUMERIC=`ls Numeric-*.tar.gz | tail -1`

if [ -z $NUMERIC ]
then
    NUMERIC=`ls Numeric-*.tar | tail -1`

    if [ -z $NUMERIC ]
    then
        echo "Can't find Numeric tar in this directory!"
        exit 1
    else
        tar xvf $NUMERIC
    fi
else
    tar zxvf $NUMERIC
fi

SWIG=`ls swig-*.tar.gz | tail -1`

if [ -z $SWIG ]
then
    SWIG=`ls swig-*.tar | tail -1`

    if [ -z $SWIG ]
    then
        echo "Can't find SWIG in this directory!"
        exit 1
    else
        tar xvf $SWIG
    fi
else
    tar zxvf $SWIG
fi

NUMERICDIR=`ls | grep Numeric | grep -v tar | head -1`

SWIGDIR=`ls | grep swig | grep -v tar | head -1`

if [ -z $SWIGDIR ]
then
    SWIGDIR=`ls | grep SWIG | head -1`
fi

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
