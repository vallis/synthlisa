#include "lisasim.h"
#include "lisasim-noise.h"
#include "lisasim-tdinoise.h"

TDInoise::TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot) {
    lisa = mylisa;

    // to estimate size of noisebuffer, take an armlength at time zero,
    // and add 10%; need to convert from years to seconds (we use the factor from lisasim-tdi.cpp)

    double lighttime = 1.10 * (lisa->armlength(1,0.0) / 3.17098E-8);

    // create InterpolateNoise objects for proof-mass noises
    // we need quadruple retardations for the V's appearing in the z's

    double pbtproof = 4.0 * lighttime;

    for(int craft = 1; craft <= 3; craft++) {
        pm[craft] = new InterpolateNoise(stproof, pbtproof, sdproof, -2.0);
        pms[craft] = new InterpolateNoise(stproof, pbtproof, sdproof, -2.0);
    }
    
    // create InterpolateNoise objects for optical-path noises
    // we need only triple retardations for the shot's appearing in the y's
    
    double pbtshot = 3.0 * lighttime;
    
    for(int craft1 = 1; craft1 <= 3; craft1++) {
        for(int craft2 = 1; craft2 <= 3; craft2++) {
            if(craft1 != craft2)
                shot[craft1][craft2] = new InterpolateNoise(stshot, pbtshot, sdshot, 2.0);
        }
    }
}

TDInoise::~TDInoise() {
    // remove proof-mass-noise InterpolateNoise objects

    for(int craft = 1; craft <= 3; craft++) {
        delete pm[craft];
        delete pms[craft];
    }
 
    // remove optical-path-noise InterpolateNoise objects

    for(int craft1 = 1; craft1 <= 3; craft1++) {
        for(int craft2 = 1; craft2 <= 3; craft2++) {
            if(craft1 != craft2)
                delete shot[craft1][craft2];
        }
    }
}

double TDInoise::y(int send, int link, int recv, int ret1, int ret2, int ret3, double t) {
    double retardedtime = t;
    
    // need to convert retardations to seconds; conversion factor from lisasim-wave
    
    if(ret1) retardedtime -= 3.1536E7 * lisa->armlength(ret1,t);
    if(ret2) retardedtime -= 3.1536E7 * lisa->armlength(ret2,t);
    if(ret3) retardedtime -= 3.1536E7 * lisa->armlength(ret3,t);

    if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2)) {
        // cyclic combination

        return( -2.0 * (*pm[recv])[retardedtime] + (*shot[send][recv])[retardedtime] );
    } else {
        // anticyclic combination

        return( 2.0 * (*pms[recv])[retardedtime] + (*shot[send][recv])[retardedtime] );
    }
}

double TDInoise::z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t) {
    double retardedtime = t;
    
    // need to convert retardations to seconds; conversion factor from lisasim-wave
    
    if(ret1) retardedtime -= 3.1536E7 * lisa->armlength(ret1,t);
    if(ret2) retardedtime -= 3.1536E7 * lisa->armlength(ret2,t);
    if(ret3) retardedtime -= 3.1536E7 * lisa->armlength(ret3,t);
    if(ret4) retardedtime -= 3.1536E7 * lisa->armlength(ret4,t);

    if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2)) {
        // cyclic combination

        return( -2.0 * (*pms[recv])[retardedtime] );
    } else {
        // anticyclic combination

        return( 2.0 * (*pm[recv])[retardedtime] );
    }
}

// time in seconds is OK here

double TDInoise::X(double t) {
    return( y(1, 3, 2, 3, 2, 2, t) -
            y(1, 2, 3, 2, 3, 3, t) +
            y(2, 3, 1, 2, 2, 0, t) -
            y(3, 2, 1, 3, 3, 0, t) +
            y(1, 2, 3, 2, 0, 0, t) -
            y(1, 3, 2, 3, 0, 0, t) + 
            y(3, 2, 1, 0, 0, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) +
    0.5 * (-z(3, 2, 1, 2, 2, 3, 3, t) +
            z(3, 2, 1, 3, 3, 0, 0, t) +
            z(3, 2, 1, 2, 2, 0, 0, t) -
            z(3, 2, 1, 0, 0, 0, 0, t) ) +
    0.5 * ( z(2, 3, 1, 2, 2, 3, 3, t) -
            z(2, 3, 1, 3, 3, 0, 0, t) -
            z(2, 3, 1, 2, 2, 0, 0, t) +
            z(2, 3, 1, 0, 0, 0, 0, t) ) );
}

double TDInoise::Y(double t) {
    return( y(2, 1, 3, 1, 3, 3, t) -
            y(2, 3, 1, 3, 1, 1, t) +
            y(3, 1, 2, 3, 3, 0, t) -
            y(1, 3, 2, 1, 1, 0, t) +
            y(2, 3, 1, 3, 0, 0, t) -
            y(2, 1, 3, 1, 0, 0, t) + 
            y(1, 3, 2, 0, 0, 0, t) -
            y(3, 1, 2, 0, 0, 0, t) + 
    0.5 * (-z(1, 3, 2, 3, 3, 1, 1, t) +
            z(1, 3, 2, 1, 1, 0, 0, t) +
            z(1, 3, 2, 3, 3, 0, 0, t) -
            z(1, 3, 2, 0, 0, 0, 0, t) ) +
    0.5 * ( z(3, 1, 2, 3, 3, 1, 1, t) -
            z(3, 1, 2, 1, 1, 0, 0, t) -
            z(3, 1, 2, 3, 3, 0, 0, t) +
            z(3, 1, 2, 0, 0, 0, 0, t) ) );
}

double TDInoise::Z(double t) {
    return( y(3, 2, 1, 2, 1, 1, t) -
            y(3, 1, 2, 1, 2, 2, t) +
            y(1, 2, 3, 1, 1, 0, t) -
            y(2, 1, 3, 2, 2, 0, t) +
            y(3, 1, 2, 1, 0, 0, t) -
            y(3, 2, 1, 2, 0, 0, t) + 
            y(2, 1, 3, 0, 0, 0, t) -
            y(1, 2, 3, 0, 0, 0, t) +
    0.5 * (-z(2, 1, 3, 1, 1, 2, 2, t) +
            z(2, 1, 3, 2, 2, 0, 0, t) +
            z(2, 1, 3, 1, 1, 0, 0, t) -
            z(2, 1, 3, 0, 0, 0, 0, t) ) +
    0.5 * ( z(1, 2, 3, 1, 1, 2, 2, t) -
            z(1, 2, 3, 2, 2, 0, 0, t) -
            z(1, 2, 3, 1, 1, 0, 0, t) +
            z(1, 2, 3, 0, 0, 0, 0, t) ) );
}

double TDInoise::alpha(double t) {
    return( y(3, 2, 1, 0, 0, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) +
            y(2, 1, 3, 2, 0, 0, t) -
            y(3, 1, 2, 3, 0, 0, t) +
            y(1, 3, 2, 1, 2, 0, t) -
            y(1, 2, 3, 1, 3, 0, t) -
    0.5 * ( z(2, 1, 3, 2, 0, 0, 0, t) +
            z(2, 1, 3, 1, 3, 0, 0, t) +
            z(3, 2, 1, 0, 0, 0, 0, t) +
            z(3, 2, 1, 1, 2, 3, 0, t) +
            z(1, 3, 2, 3, 0, 0, 0, t) +
            z(1, 3, 2, 1, 2, 0, 0, t) ) +
    0.5 * ( z(1, 2, 3, 2, 0, 0, 0, t) +
            z(1, 2, 3, 1, 3, 0, 0, t) +
            z(2, 3, 1, 0, 0, 0, 0, t) +
            z(2, 3, 1, 1, 2, 3, 0, t) + 
            z(3, 1, 2, 3, 0, 0, 0, t) +
            z(3, 1, 2, 1, 2, 0, 0, t) ) );
}

double TDInoise::beta(double t) {
    return( y(1, 3, 2, 0, 0, 0, t) -
            y(3, 1, 2, 0, 0, 0, t) +
            y(3, 2, 1, 3, 0, 0, t) -
            y(1, 2, 3, 1, 0, 0, t) +
            y(2, 1, 3, 2, 3, 0, t) -
            y(2, 3, 1, 2, 1, 0, t) -
    0.5 * ( z(3, 2, 1, 3, 0, 0, 0, t) +
            z(3, 2, 1, 2, 1, 0, 0, t) +
            z(1, 3, 2, 0, 0, 0, 0, t) +
            z(1, 3, 2, 2, 3, 1, 0, t) +
            z(2, 1, 3, 1, 0, 0, 0, t) +
            z(2, 1, 3, 2, 3, 0, 0, t) ) +
    0.5 * ( z(2, 3, 1, 3, 0, 0, 0, t) +
            z(2, 3, 1, 2, 1, 0, 0, t) +
            z(3, 1, 2, 0, 0, 0, 0, t) +
            z(3, 1, 2, 2, 3, 1, 0, t) + 
            z(1, 2, 3, 1, 0, 0, 0, t) +
            z(1, 2, 3, 2, 3, 0, 0, t) ) );
}

double TDInoise::gamma(double t) {
    return( y(2, 1, 3, 0, 0, 0, t) -
            y(1, 2, 3, 0, 0, 0, t) +
            y(1, 3, 2, 1, 0, 0, t) -
            y(2, 3, 1, 2, 0, 0, t) +
            y(3, 2, 1, 3, 1, 0, t) -
            y(3, 1, 2, 3, 2, 0, t) -
    0.5 * ( z(1, 3, 2, 1, 0, 0, 0, t) +
            z(1, 3, 2, 3, 2, 0, 0, t) +
            z(2, 1, 3, 0, 0, 0, 0, t) +
            z(2, 1, 3, 3, 1, 2, 0, t) +
            z(3, 2, 1, 2, 0, 0, 0, t) +
            z(3, 2, 1, 3, 1, 0, 0, t) ) +
    0.5 * ( z(3, 1, 2, 1, 0, 0, 0, t) +
            z(3, 1, 2, 3, 2, 0, 0, t) +
            z(1, 2, 3, 0, 0, 0, 0, t) +
            z(1, 2, 3, 3, 1, 2, 0, t) + 
            z(2, 3, 1, 2, 0, 0, 0, t) +
            z(2, 3, 1, 3, 1, 0, 0, t) ) );
}

double TDInoise::zeta(double t) {
    return( y(1, 3, 2, 2, 0, 0, t) -
            y(1, 2, 3, 3, 0, 0, t) +
            y(2, 1, 3, 3, 0, 0, t) -
            y(2, 3, 1, 1, 0, 0, t) +
            y(3, 2, 1, 1, 0, 0, t) -
            y(3, 1, 2, 2, 0, 0, t) +
    0.5 * (-z(2, 1, 3, 2, 1, 0, 0, t) +
            z(1, 2, 3, 1, 2, 0, 0, t) -
            z(3, 2, 1, 2, 3, 0, 0, t) +
            z(2, 3, 1, 2, 3, 0, 0, t) -
            z(1, 3, 2, 1, 3, 0, 0, t) +
            z(3, 1, 2, 1, 3, 0, 0, t) ) +
    0.5 * (-z(1, 3, 2, 2, 0, 0, 0, t) +
            z(3, 1, 2, 2, 0, 0, 0, t) -
            z(2, 1, 3, 3, 0, 0, 0, t) +
            z(1, 2, 3, 3, 0, 0, 0, t) -
            z(3, 2, 1, 1, 0, 0, 0, t) +
            z(2, 3, 1, 1, 0, 0, 0, t) ) );
}

double TDInoise::P(double t) {
    return( y(1, 3, 2, 2, 0, 0, t) -
            y(1, 2, 3, 3, 0, 0, t) - 
            y(3, 1, 2, 2, 0, 0, t) +
            y(2, 1, 3, 3, 0, 0, t) +
            y(3, 1, 2, 1, 3, 0, t) -
            y(2, 1, 3, 1, 2, 0, t) + 
            y(1, 2, 3, 3, 1, 1, t) -
            y(1, 3, 2, 2, 1, 1, t) +
    0.5 * (-z(3, 2, 1, 2, 3, 0, 0, t) +
            z(3, 2, 1, 1, 1, 2, 3, t) + 
            z(2, 3, 1, 2, 3, 0, 0, t) -
            z(2, 3, 1, 1, 1, 2, 3, t) ) +
    0.5 * (-z(1, 3, 2, 2, 0, 0, 0, t) + 
            z(1, 3, 2, 1, 1, 2, 0, t) +
            z(3, 1, 2, 2, 0, 0, 0, t) -
            z(3, 1, 2, 1, 1, 2, 0, t) ) +
    0.5 * (-z(2, 1, 3, 3, 0, 0, 0, t) +
            z(2, 1, 3, 1, 1, 3, 0, t) +
            z(1, 2, 3, 3, 0, 0, 0, t) -
            z(1, 2, 3, 1, 1, 3, 0, t) ) );
}

double TDInoise::E(double t) {
    return( y(3, 1, 2, 2, 1, 0, t) -
            y(2, 1, 3, 3, 1, 0, t) -
            y(3, 1, 2, 3, 0, 0, t) +
            y(2, 1, 3, 2, 0, 0, t) +
            y(2, 3, 1, 1, 1, 0, t) -
            y(3, 2, 1, 1, 1, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) +
            y(3, 2, 1, 0, 0, 0, t) -
    0.5 * ( z(2, 1, 3, 2, 0, 0, 0, t) +
            z(3, 2, 1, 0, 0, 0, 0, t) +
            z(1, 3, 2, 3, 0, 0, 0, t) -
            z(2, 1, 3, 1, 1, 2, 0, t) +
            z(1, 2, 3, 1, 1, 2, 0, t) -
            z(1, 3, 2, 1, 1, 3, 0, t) ) +
    0.5 * ( z(1, 2, 3, 2, 0, 0, 0, t) +
            z(2, 3, 1, 0, 0, 0, 0, t) +
            z(3, 1, 2, 3, 0, 0, 0, t) -
            z(3, 1, 2, 1, 1, 3, 0, t) +
            z(3, 2, 1, 1, 1, 0, 0, t) -
            z(2, 3, 1, 1, 1, 0, 0, t) ) );
};
