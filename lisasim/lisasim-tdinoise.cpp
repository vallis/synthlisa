#include "lisasim-lisa.h"
#include "lisasim-noise.h"
#include "lisasim-tdinoise.h"
#include <time.h>

TDInoise::TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser) {
    lisa = mylisa;
    phlisa = mylisa;

    initialize(stproof, sdproof, stshot, sdshot, stlaser, sdlaser, claser);
}

TDInoise::TDInoise(LISA *mylisa, LISA *physlisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser) {
    lisa = mylisa;
    phlisa = physlisa;

    initialize(stproof, sdproof, stshot, sdshot, stlaser, sdlaser, claser);
}

void TDInoise::initialize(double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser) {
    // to estimate size of noisebuffer, take an armlength at time zero,
    // and add 10%; need to convert from years to seconds (we use the factor from lisasim-tdi.cpp)
    // ah, this might not work for OriginalLISA and strange geometries

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

    // create laser noise objects
    // quadruple retardations are needed for the C's

    double pbtlaser = 4.0 * lighttime;
    
     for(int craft = 1; craft <= 3; craft++) {
        c[craft] = new ExpGaussNoise(stlaser,pbtlaser,1/claser,sdlaser);
        cs[craft] = new ExpGaussNoise(stlaser,pbtlaser,1/claser,sdlaser);
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

    // remove laser-noise ExpGaussNoise objects

    for(int craft = 1; craft <= 3; craft++) {
        delete c[craft];
        delete cs[craft];
    }
}

void TDInoise::reset() {
    // initialize random-number-generator seed
    // here we're actually passing a long (on PPC); should work as long
    // as "long" is the same as "int"

    idum = -time(0);

    for(int craft = 1; craft <= 3; craft++) {
        pm[craft]->reset();
        pms[craft]->reset();
    }
 
    // remove optical-path-noise InterpolateNoise objects

    for(int craft1 = 1; craft1 <= 3; craft1++) {
        for(int craft2 = 1; craft2 <= 3; craft2++) {
            if(craft1 != craft2)
                shot[craft1][craft2]->reset();
        }
    }

    // remove laser-noise ExpGaussNoise objects

    for(int craft = 1; craft <= 3; craft++) {
        c[craft]->reset();
        cs[craft]->reset();
    }

    // reset also LISA, in case it includes noise of some kind
    
    lisa->reset();
    if(phlisa != lisa) phlisa->reset();
}

double TDInoise::y(int send, int slink, int recv, int ret1, int ret2, int ret3, double t) {
//    cout << "Doing y(" << send << "," << slink << "," << recv << "," << ret1 << "," << ret2 << "," << ret3 << ")" << endl;

    int link = abs(slink);

    // need to convert retardations to seconds; conversion factor from lisasim-wave
    // but armlength takes a time in years (when it does)

    // this recursive retardation procedure assumes smart TDI...

    double retardedtime = t;

    if(ret3 != 0) retardedtime -= 3.1536E7 * lisa->armlength(ret3, 3.17098E-8 * retardedtime);
    if(ret2 != 0) retardedtime -= 3.1536E7 * lisa->armlength(ret2, 3.17098E-8 * retardedtime);    
    if(ret1 != 0) retardedtime -= 3.1536E7 * lisa->armlength(ret1, 3.17098E-8 * retardedtime);

    if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2)) {
        // cyclic combination

        // if introducing error in the determination of the armlengths, it should not enter
        // the following (physical) retardation of the laser noise

        double retardlaser = retardedtime - 3.1536E7 * phlisa->armlength(link,3.17098E-8 * retardedtime);

        cout.precision(16);

        if(t==0.00 && recv==1) {
//            cout << "y(" << send << "," << slink << "," << recv << "," << ret1 << "," << ret2 << "," << ret3 << "); " << "-c(1," << retardedtime << "): " << -(*c[recv])[retardedtime] << endl; 
        }

        return( (*cs[send])[retardlaser] - 2.0 * (*pm[recv])[retardedtime]  - (*c[recv])[retardedtime]  + 
                (*shot[send][recv])[retardedtime] );
    } else {
        // anticyclic combination

        // ditto here

        double retardlaser = retardedtime - 3.1536E7 * phlisa->armlength(-link,3.17098E-8 * retardedtime);

        cout.precision(16);

        if(t==0.00 && send==1) {
//            cout << "y(" << send << "," << slink << "," << recv << "," << ret1 << "," << ret2 << "," << ret3 << "); " << "c(1," << retardlaser << "): " << (*c[send])[retardlaser] << endl; 
        }

        return( (*c[send])[retardlaser]  + 2.0 * (*pms[recv])[retardedtime] - (*cs[recv])[retardedtime] +
                (*shot[send][recv])[retardedtime] );
    }
}

double TDInoise::z(int send, int slink, int recv, int ret1, int ret2, int ret3, int ret4, double t) {
    int link = abs(slink);

//    cout << "Doing z(" << send << "," << slink << "," << recv << "," << ret1 << "," << ret2 << "," << ret3 << "," << ret4 << ")" << endl;

    // need to convert retardations to seconds; conversion factor from lisasim-wave

    // this recursive retardation procedure assumes smart TDI...
    // (and the correct order in the retardation expressions)

    double retardedtime = t;

    if(ret4 != 0) retardedtime -= 3.1536E7 * lisa->armlength(ret4,3.17098E-8 * retardedtime);
    if(ret3 != 0) retardedtime -= 3.1536E7 * lisa->armlength(ret3,3.17098E-8 * retardedtime);
    if(ret2 != 0) retardedtime -= 3.1536E7 * lisa->armlength(ret2,3.17098E-8 * retardedtime);    
    if(ret1 != 0) retardedtime -= 3.1536E7 * lisa->armlength(ret1,3.17098E-8 * retardedtime);

    cout.precision(16);
    
    if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2)) {
        // cyclic combination

        if(t==0.00 && recv==1) {
//            cout << "z(" << send << "," << slink << "," << recv << "," << ret1 << "," << ret2 << "," << ret3 << "," << ret4 << "); " << "-c(1," << retardedtime << "): " <<  -(*c[recv])[retardedtime] << endl; 
        }

        return( (*cs[recv])[retardedtime] - 2.0 * (*pms[recv])[retardedtime] - (*c[recv])[retardedtime] );
    } else {
        // anticyclic combination

        if(t==0.00 && recv==1) {
//            cout << "z(" << send << "," << slink << "," << recv << "," << ret1 << "," << ret2 << "," << ret3 << "," << ret4 << "); " << "c(1," << retardedtime << "):" << (*c[recv])[retardedtime] << endl; 
        }

        return( (*c[recv])[retardedtime]  + 2.0 * (*pm[recv])[retardedtime]  - (*cs[recv])[retardedtime] );
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

double TDInoise::Xm(double t) {
    return( y(1,-3, 2, 3, 2,-2, t) -
            y(1, 2, 3,-2,-3, 3, t) +
            y(2, 3, 1, 2,-2, 0, t) -
            y(3,-2, 1,-3, 3, 0, t) +
            y(1, 2, 3,-2, 0, 0, t) -
            y(1,-3, 2, 3, 0, 0, t) + 
            y(3,-2, 1, 0, 0, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) +
    0.5 * (-z(3,-2, 1, 2,-2,-3, 3, t) +
            z(3,-2, 1,-3, 3, 0, 0, t) +
            z(3,-2, 1, 2,-2, 0, 0, t) -
            z(3,-2, 1, 0, 0, 0, 0, t) ) +
    0.5 * ( z(2, 3, 1, 2,-2,-3, 3, t) -
            z(2, 3, 1,-3, 3, 0, 0, t) -
            z(2, 3, 1, 2,-2, 0, 0, t) +
            z(2, 3, 1, 0, 0, 0, 0, t) ) );
}
