#!/bin/sh

set PLATFORM=`uname -s`-`uname -r`

set PLADIR=`pwd`/${PLATFORM}

set LISASWIGDIR=${PLADIR}/lib
set PYTHONDIR=${PLADIR}/lib/`ls $PLADIR/lib | grep python | head -1`/site-packages
set NUMERICDIR=${PLADIR}/lib/`ls $PLADIR/lib | grep python | head -1`/site-packages/Numeric

setenv PYTHONPATH ${LISASWIGDIR}:${NUMERICDIR}:${PYTHONDIR}
