#include "lisasim-tdi.h"

// double TDI::psi(int arm, double t) {
double TDI::psi(int arm, double t, double twave) {
    Vector lisan;
    lisa->putn(lisan,arm,t);
    
    Tensor cwave;
    wave->putwave(cwave,twave);
    
    Vector tmp;
    tmp.setproduct(cwave,lisan);

    return(0.5 * lisan.dotproduct(tmp));
}

double TDI::retard(int craft, double t) {
    Vector lisap;
    lisa->putp(lisap,craft,t);
    
    return(-lisap.dotproduct(wave->k));
}

double TDI::y(int send, int link, int recv, int ret1, int ret2, int ret3, double t) {
    double retardedtime = t;
    
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

//    return (( psi(link, retardedtime + retard(send, t) - lisa->armlength(link,t)) -
//              psi(link, retardedtime + retard(recv, t)) ) / denom );
    
    return (( psi(link, retardedtime - lisa->armlength(link,t), retardedtime + retard(send, t) - lisa->armlength(link,t)) -
              psi(link, retardedtime, retardedtime + retard(recv, t)) ) / denom );
}

double TDI::M(double time) {
    double t = time * 3.17098E-8;
    
    return( y(1, 2, 3, 2, 0, 0, t) -
            y(1, 3, 2, 3, 0, 0, t) + 
            y(3, 2, 1, 0, 0, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) );
}

double TDI::N(double time) {
    double t = time * 3.17098E-8;
        
    return( y(2, 3, 1, 3, 0, 0, t) -
            y(2, 1, 3, 1, 0, 0, t) + 
            y(1, 3, 2, 0, 0, 0, t) -
            y(3, 1, 2, 0, 0, 0, t) );
}

double TDI::O(double time) {
    double t = time * 3.17098E-8;

    return( y(3, 1, 2, 1, 0, 0, t) -
            y(3, 2, 1, 2, 0, 0, t) + 
            y(2, 1, 3, 0, 0, 0, t) -
            y(1, 2, 3, 0, 0, 0, t) );
}

double TDI::X(double time) {
    double t = time * 3.17098E-8;
    
    return( y(1, 3, 2, 3, 2, 2, t) -
            y(1, 2, 3, 2, 3, 3, t) +
            y(2, 3, 1, 2, 2, 0, t) -
            y(3, 2, 1, 3, 3, 0, t) +
            y(1, 2, 3, 2, 0, 0, t) -
            y(1, 3, 2, 3, 0, 0, t) + 
            y(3, 2, 1, 0, 0, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) );
}

double TDI::Y(double time) {
    double t = time * 3.17098E-8;
        
    return( y(2, 1, 3, 1, 3, 3, t) -
            y(2, 3, 1, 3, 1, 1, t) +
            y(3, 1, 2, 3, 3, 0, t) -
            y(1, 3, 2, 1, 1, 0, t) +
            y(2, 3, 1, 3, 0, 0, t) -
            y(2, 1, 3, 1, 0, 0, t) + 
            y(1, 3, 2, 0, 0, 0, t) -
            y(3, 1, 2, 0, 0, 0, t) );
}

double TDI::Z(double time) {
    double t = time * 3.17098E-8;

    return( y(3, 2, 1, 2, 1, 1, t) -
            y(3, 1, 2, 1, 2, 2, t) +
            y(1, 2, 3, 1, 1, 0, t) -
            y(2, 1, 3, 2, 2, 0, t) +
            y(3, 1, 2, 1, 0, 0, t) -
            y(3, 2, 1, 2, 0, 0, t) + 
            y(2, 1, 3, 0, 0, 0, t) -
            y(1, 2, 3, 0, 0, 0, t) );
}

double TDI::alpha(double time) {
    double t = time * 3.17098E-8;
        
    return( y(3, 2, 1, 0, 0, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) +
            y(2, 1, 3, 2, 0, 0, t) -
            y(3, 1, 2, 3, 0, 0, t) +
            y(1, 3, 2, 1, 2, 0, t) -
            y(1, 2, 3, 1, 3, 0, t)  );
}

double TDI::beta(double time) {
    double t = time * 3.17098E-8;
        
    return( y(1, 3, 2, 0, 0, 0, t) -
            y(3, 1, 2, 0, 0, 0, t) +
            y(3, 2, 1, 3, 0, 0, t) -
            y(1, 2, 3, 1, 0, 0, t) +
            y(2, 1, 3, 2, 3, 0, t) -
            y(2, 3, 1, 2, 1, 0, t)  );
}

double TDI::gamma(double time) {
    double t = time * 3.17098E-8;
        
    return( y(2, 1, 3, 0, 0, 0, t) -
            y(1, 2, 3, 0, 0, 0, t) +
            y(1, 3, 2, 1, 0, 0, t) -
            y(2, 3, 1, 2, 0, 0, t) +
            y(3, 2, 1, 3, 1, 0, t) -
            y(3, 1, 2, 3, 2, 0, t)  );
}

double TDI::zeta(double time) {
    double t = time * 3.17098E-8;
        
    return( y(1, 3, 2, 2, 0, 0, t) -
            y(1, 2, 3, 3, 0, 0, t) +
            y(2, 1, 3, 3, 0, 0, t) -
            y(2, 3, 1, 1, 0, 0, t) +
            y(3, 2, 1, 1, 0, 0, t) -
            y(3, 1, 2, 2, 0, 0, t)  );
}

