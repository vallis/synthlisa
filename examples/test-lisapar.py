#!/usr/bin/env python

# $Id$
# $Date$
# $Author$
# $Revision$

import synthlisa

def myfunc(x):
    return synthlisa.SimpleBinary(x*1e-3,0,0,1/x,0.5,0.5,0.5)

lisa = synthlisa.OriginalLISA()

samples = 2**14
stime = 0.25

lp = synthlisa.lisapar()

observables = lp.getobsp(samples,stime,(lisa,myfunc,[0.5,2.0,4.0],(synthlisa.TDI.time,synthlisa.TDI.Xm)),0,0)

if lp.rank == 0:
    outputXML = synthlisa.lisaXML('data/tdipar',
                                  author='Michele Vallisneri',
                                  comments='Parallel generation of SimpleBinary waves')
    outputXML.TDIData(observables,samples,stime,'t,Xmf',encoding='Binary')
    outputXML.close()
