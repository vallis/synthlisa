#include "lisasim-tdisignal.h"

TDIsignal::TDIsignal(LISA *mylisa, Wave *mywave) {
    lisa = mylisa;
    wave = mywave;
}

double TDIsignal::psi(int arm, double t, double twave) {
    Vector lisan;
    lisa->putn(lisan,arm,t);
    
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

double TDIsignal::y(int send, int slink, int recv, int ret1, int ret2, int ret3, double t) {
    int link = abs(slink);

    double retardedtime = t;

    // this scheme needs to be updated in accordance with lisasim-tdinoise.cpp
    
    if(ret1) retardedtime -= lisa->armlength(ret1,t);
    if(ret2) retardedtime -= lisa->armlength(ret2,t);
    if(ret3) retardedtime -= lisa->armlength(ret3,t);

    Vector linkn;
    lisa->putn(linkn,link,t);
    
    double denom = linkn.dotproduct(wave->k);
    if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2))
        denom = 1.0 - denom;
    else
        denom = 1.0 + denom;

    return (( psi(link, retardedtime - lisa->armlength(link,t), retardedtime + retard(send, t) - lisa->armlength(link,t)) -
              psi(link, retardedtime, retardedtime + retard(recv, t)) ) / denom );
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

