#!/bin/sh

PLATFORM=`uname -s`-`uname -r`

PLADIR=`pwd`/$PLATFORM

LISASWIGDIR=$PLADIR/lib
PYTHONDIR=$PLADIR/lib/`ls $PLADIR/lib | grep python | head -1`/site-packages
NUMERICDIR=$PLADIR/lib/`ls $PLADIR/lib | grep python | head -1`/site-packages/Numeric

export PYTHONPATH=$LISASWIGDIR:$NUMERICDIR:$PYTHONDIR
