#!/bin/sh

# $Id$

NUMERIC=`ls Numeric-*.tar.gz | tail -1`

if [ -z $NUMERIC ]
then
    NUMERIC=`ls Numeric-*.tar | tail -1`

    if [ -z $NUMERIC ]
    then
        echo "Can't find Numeric tar in this directory!"
        return 1
    else
        tar xvf $NUMERIC
    fi
else
    tar zxvf $NUMERIC
fi

NUMERICDIR=`ls | grep Numeric | head -1`
THISDIR=`pwd`

cd $NUMERICDIR

python setup.py install --prefix=$THISDIR

cd ..

python setup.py install --with-numeric=$THISDIR --prefix=$THISDIR
