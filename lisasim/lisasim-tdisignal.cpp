#include "lisasim-tdisignal.h"

TDIsignal::TDIsignal(LISA *mylisa, Wave *mywave) {
    lisa = mylisa->thislisa();
    phlisa = mylisa;
    
    wave = mywave;
}

void TDIsignal::setphlisa(LISA *mylisa) {
    phlisa = mylisa;
}

// the argument t is actually redundant here!

double TDIsignal::psi(Vector &lisan, double t) {
    Tensor cwave;
    wave->putwave(cwave,t);
    
    Vector tmp;
    tmp.setproduct(cwave,lisan);

    return(0.5 * lisan.dotproduct(tmp));
}

// the lisa used here should be the physical one, not the nominal one (02/13/04)

double TDIsignal::retard(int craft, double t) {
    Vector lisap;
    phlisa->putp(lisap,craft,t);
    
    return(-lisap.dotproduct(wave->k));
}

// the y as computed below should now be fully covariant
// the time of evaluation (reception) is retarded according to the nominal lisa
// then the time of emission is retarded according to physical lisa
// the wave retardation for the emitting link is taken at the correct retarded time

// note that the call to putn is already computing the armlength that we request later
// we could make putn return the corresponding armlength

double TDIsignal::y(int send, int slink, int recv, int ret1, int ret2, int ret3, double t) {
    return y(send,slink,recv,ret1,ret2,ret3,0,0,0,0,t);
}

double TDIsignal::y(int send, int slink, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t) {
    int link = abs(slink);

    retardtime myrt(lisa,t);

    myrt.retard(ret7); myrt.retard(ret6); myrt.retard(ret5);
    myrt.retard(ret4); myrt.retard(ret3); myrt.retard(ret2); myrt.retard(ret1);

    double retardedtime = myrt.retardedtime();

    // for the moment, let's not trust the sign of slink, and recompute it

    if( (link == 3 && recv == 2) || (link == 1 && recv == 3) || (link == 2 && recv == 1) )
	link = -link;

    // since the linkn returned by the "modern" versions of putn is oriented,
    // therefore the sign in denom is the same (-) for both positive and negative links
    // there's no problem in psi, because n is dotted twice into h

    Vector linkn;
    lisa->putn(linkn,link,retardedtime);    
    
    double denom = 1.0 - linkn.dotproduct(wave->k);
    double retardsignal = retardedtime - phlisa->armlength(link,retardedtime);

    return (( psi(linkn, retardsignal + retard(send, retardsignal)) -
              psi(linkn, retardedtime + retard(recv, retardedtime)) ) / denom);
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

