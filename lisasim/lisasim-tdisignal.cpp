#include "lisasim-tdisignal.h"

TDIsignal::TDIsignal(LISA *mylisa, Wave *mywave) {
    lisa = mylisa;
    phlisa = mylisa;
    
    wave = mywave;
}

TDIsignal::TDIsignal(LISA *mylisa, LISA *physlisa, Wave *mywave) {
    lisa = mylisa;
    phlisa = physlisa;

    wave = mywave;
}

double TDIsignal::psi(Vector &lisan, double t, double twave) {
    Tensor cwave;
    wave->putwave(cwave,twave);
    
    Vector tmp;
    tmp.setproduct(cwave,lisan);

    return(0.5 * lisan.dotproduct(tmp));
}

double TDIsignal::retard(int craft, double t) {
    Vector lisap;
    lisa->putp(lisap,craft,t);
    
    return(-lisap.dotproduct(wave->k));
}

// the y as computed below should now be fully covariant
// the time of evaluation (reception) is retarded according to the nominal lisa
// then the time of emission is retarded according to physical lisa
// the wave retardation for the emitting link is taken at the correct retarded time

// note that the call to putn is already computing the armlength that we request later
// we could make putn return the corresponding armlength

double TDIsignal::y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t) {
    double retardedtime = t;

    if(ret7 != 0) retardedtime -= lisa->armlength(ret7,retardedtime);
    if(ret6 != 0) retardedtime -= lisa->armlength(ret6,retardedtime);
    if(ret5 != 0) retardedtime -= lisa->armlength(ret5,retardedtime);    
    if(ret4 != 0) retardedtime -= lisa->armlength(ret4,retardedtime);

    return y(send,link,recv,ret1,ret2,ret3,retardedtime);
}

double TDIsignal::y(int send, int slink, int recv, int ret1, int ret2, int ret3, double t) {
    int link = abs(slink);

    double retardedtime = t;

    if(ret3 != 0) retardedtime -= lisa->armlength(ret3,retardedtime);
    if(ret2 != 0) retardedtime -= lisa->armlength(ret2,retardedtime);    
    if(ret1 != 0) retardedtime -= lisa->armlength(ret1,retardedtime);

    Vector linkn;
    double denom, retardsignal;

    // for the moment, let's not trust the sign of slink, and recompute it
    // however, the linkn returned by the "modern" versions of putn is oriented,
    // therefore the sign in denom should be the same (-) for both cases
    // there's no problem in psi, because n is dotted twice into h
 
    if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2)) {
        lisa->putn(linkn,link,retardedtime);    
        denom = 1.0 - linkn.dotproduct(wave->k);

        retardsignal = retardedtime - phlisa->armlength(link,retardedtime);
    } else {
        lisa->putn(linkn,-link,retardedtime);    
        denom = 1.0 - linkn.dotproduct(wave->k);

        retardsignal = retardedtime - phlisa->armlength(-link,retardedtime);
    }

    return (( psi(linkn, retardsignal, retardsignal + retard(send, retardsignal)) -
              psi(linkn, retardedtime, retardedtime + retard(recv, retardedtime)) ) / denom);
}

// the next three names are not standard!

double TDIsignal::M(double t) {
    return( y(1, 2, 3, 2, 0, 0, t) -
            y(1, 3, 2, 3, 0, 0, t) + 
            y(3, 2, 1, 0, 0, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) );
}

double TDIsignal::N(double t) {
    return( y(2, 3, 1, 3, 0, 0, t) -
            y(2, 1, 3, 1, 0, 0, t) + 
            y(1, 3, 2, 0, 0, 0, t) -
            y(3, 1, 2, 0, 0, 0, t) );
}

double TDIsignal::O(double t) {
    return( y(3, 1, 2, 1, 0, 0, t) -
            y(3, 2, 1, 2, 0, 0, t) + 
            y(2, 1, 3, 0, 0, 0, t) -
            y(1, 2, 3, 0, 0, 0, t) );
}

