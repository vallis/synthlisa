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
    // to estimate size of noisebuffer, take maximum armlength at time zero,
    // and add 10% for uplink-downlink uncertainty and flexing

    double arm1 = lisa->armlength(1,0.0);
    double arm2 = lisa->armlength(2,0.0);
    double arm3 = lisa->armlength(3,0.0);

    double maxarm = arm1 > arm2 ? arm1 : arm2;
    maxarm = maxarm > arm3 ? maxarm : arm3;

    double lighttime = 1.10 * maxarm;

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
        c[craft] = new ExpGaussNoise(stlaser,pbtlaser,claser,sdlaser);
        cs[craft] = new ExpGaussNoise(stlaser,pbtlaser,claser,sdlaser);
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
    int link = abs(slink);

    // this recursive retardation procedure assumes smart TDI...

    double retardedtime = t;

    if(ret3 != 0) retardedtime -= lisa->armlength(ret3,retardedtime);
    if(ret2 != 0) retardedtime -= lisa->armlength(ret2,retardedtime);    
    if(ret1 != 0) retardedtime -= lisa->armlength(ret1,retardedtime);

    if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2)) {
        // cyclic combination

        // if introducing error in the determination of the armlengths, it should not enter
        // the following (physical) retardation of the laser noise, so we use the phlisa object

        double retardlaser = retardedtime - phlisa->armlength(link,retardedtime);

        return( (*cs[send])[retardlaser] - 2.0 * (*pm[recv])[retardedtime]  - (*c[recv])[retardedtime]  + 
                (*shot[send][recv])[retardedtime] );
    } else {
        // anticyclic combination

        // ditto here

        double retardlaser = retardedtime - phlisa->armlength(-link,retardedtime);

        return( (*c[send])[retardlaser]  + 2.0 * (*pms[recv])[retardedtime] - (*cs[recv])[retardedtime] +
                (*shot[send][recv])[retardedtime] );
    }
}

double TDInoise::z(int send, int slink, int recv, int ret1, int ret2, int ret3, int ret4, double t) {
    int link = abs(slink);

    // this recursive retardation procedure assumes smart TDI...
    // (and the correct order in the retardation expressions)

    double retardedtime = t;

    if(ret4 != 0) retardedtime -= lisa->armlength(ret4,retardedtime);
    if(ret3 != 0) retardedtime -= lisa->armlength(ret3,retardedtime);
    if(ret2 != 0) retardedtime -= lisa->armlength(ret2,retardedtime);    
    if(ret1 != 0) retardedtime -= lisa->armlength(ret1,retardedtime);
    
    if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2)) {
        // cyclic combination

        return( (*cs[recv])[retardedtime] - 2.0 * (*pms[recv])[retardedtime] - (*c[recv])[retardedtime] );
    } else {
        // anticyclic combination

        return( (*c[recv])[retardedtime]  + 2.0 * (*pm[recv])[retardedtime]  - (*cs[recv])[retardedtime] );
    }
}

