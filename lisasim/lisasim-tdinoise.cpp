#include "lisasim-tdinoise.h"

#include <time.h>

// this version takes the parameters of the basic noises and lets us allocate objects as needed

TDInoise::TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser) {
    lisa = mylisa->thislisa();
    phlisa = mylisa;

    // allocate noise objects

    for(int craft = 1; craft <= 3; craft++) {
        pm[craft] = stdproofnoise(lisa,stproof,sdproof);
        pms[craft] = stdproofnoise(lisa,stproof,sdproof);
    }
        
    for(int craft1 = 1; craft1 <= 3; craft1++) {
        for(int craft2 = 1; craft2 <= 3; craft2++) {
            if(craft1 != craft2)
                shot[craft1][craft2] = stdopticalnoise(lisa,stshot,sdshot);
        }
    }

    for(int craft = 1; craft <= 3; craft++) {
        c[craft] = stdlasernoise(lisa,stlaser,sdlaser,claser);
        cs[craft] = stdlasernoise(lisa,stlaser,sdlaser,claser);
    }

    allocated = 1;
}

// this version takes arrays of basic-noise parameters, allowing for different noises on different objects,
// and lets us allocate objects as needed

TDInoise::TDInoise(LISA *mylisa, double *stproof, double *sdproof, double *stshot, double *sdshot, double *stlaser, double *sdlaser, double *claser) {
    lisa = mylisa->thislisa();
    phlisa = mylisa;

    // allocate noise objects
    // the convention is {1,1*,2,2*,3,3*}, and {12,21,23,32,31,13}

    for(int craft = 1; craft <= 3; craft++) {
        pm[craft]  = stdproofnoise(lisa,stproof[2*(craft-1)],  sdproof[2*(craft-1)]);
        pms[craft] = stdproofnoise(lisa,stproof[2*(craft-1)+1],sdproof[2*(craft-1)+1]);
    }
        
    for(int craft1 = 1; craft1 <= 3; craft1++) {
        for(int craft2 = 1; craft2 <= 3; craft2++) {
            if(craft1 != craft2) {
		if( (craft1 == 1 && craft2 == 2) || (craft1 == 2 && craft2 == 3) || (craft1 == 3 && craft2 == 1) )
		    shot[craft1][craft2] = stdopticalnoise(lisa,stshot[2*(craft1-1)],  sdshot[2*(craft1-1)]);
		else
		    shot[craft1][craft2] = stdopticalnoise(lisa,stshot[2*(craft2-1)+1],sdshot[2*(craft2-1)+1]);
	    }
        }
    }

    for(int craft = 1; craft <= 3; craft++) {
        c[craft]  = stdlasernoise(lisa,stlaser[2*(craft-1)],  sdlaser[2*(craft-1)],  claser[2*(craft-1)]);
        cs[craft] = stdlasernoise(lisa,stlaser[2*(craft-1)+1],sdlaser[2*(craft-1)+1],claser[2*(craft-1)+1]);
    }

    allocated = 1;
}

// this version takes pointers to noise objects, allowing for user-specified noises on different objects

TDInoise::TDInoise(LISA *mylisa, Noise *proofnoise[6],Noise *shotnoise[6],Noise *lasernoise[6]) {
    lisa = mylisa->thislisa();
    phlisa = mylisa;

    // set noise objects
    // the convention is {1,1*,2,2*,3,3*}, and {12,21,23,32,31,13}

    for(int craft = 1; craft <= 3; craft++) {
        pm[craft] = proofnoise[2*(craft-1)];
        pms[craft] = proofnoise[2*(craft-1)+1];
    }
        
    for(int craft1 = 1; craft1 <= 3; craft1++) {
        for(int craft2 = 1; craft2 <= 3; craft2++) {
            if(craft1 != craft2) {
		if( (craft1 == 1 && craft2 == 2) || (craft1 == 2 && craft2 == 3) || (craft1 == 3 && craft2 == 1) )
		    shot[craft1][craft2] = shotnoise[2*(craft1-1)];
		else
		    shot[craft1][craft2] = shotnoise[2*(craft2-1)+1];
	    }
        }
    }

    for(int craft = 1; craft <= 3; craft++) {
        c[craft] = lasernoise[2*(craft-1)];
        cs[craft] = lasernoise[2*(craft-1)];
    }

    allocated = 0;
}

void TDInoise::setphlisa(LISA *mylisa) {
    phlisa = mylisa;
}

TDInoise::~TDInoise() {
    if(allocated) {
	// allow for one noise object to be contained in multiple pointers
	// without calling delete twice

	// remove proof-mass-noise InterpolateNoise objects

	for(int craft = 1; craft <= 3; craft++) {
	    if(pm[craft])  {delete pm[craft]; pm[craft]=0;}
	    if(pms[craft]) {delete pms[craft]; pms[craft]=0;}
	}
 
	// remove optical-path-noise InterpolateNoise objects

	for(int craft1 = 1; craft1 <= 3; craft1++) {
	    for(int craft2 = 1; craft2 <= 3; craft2++) {
		if(craft1 != craft2)
		    if(shot[craft1][craft2]) {delete shot[craft1][craft2]; shot[craft1][craft2]=0;}
	    }
	}
	
	// remove laser-noise ExpGaussNoise objects

	for(int craft = 1; craft <= 3; craft++) {
	    if(c[craft])  {delete c[craft]; c[craft]=0;}
	    if(cs[craft]) {delete cs[craft]; cs[craft]=0;}
	}
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
 
    // reset optical-path-noise InterpolateNoise objects

    for(int craft1 = 1; craft1 <= 3; craft1++) {
        for(int craft2 = 1; craft2 <= 3; craft2++) {
            if(craft1 != craft2)
                shot[craft1][craft2]->reset();
        }
    }

    // reset laser-noise ExpGaussNoise objects

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

// standard noises for TDI, with utility function

double lighttime(LISA *lisa) {
    // to estimate size of noisebuffer, take maximum armlength at time zero,
    // and add 10% for uplink-downlink uncertainty and flexing

    double arm1 = lisa->armlength(1,0.0);
    double arm2 = lisa->armlength(2,0.0);
    double arm3 = lisa->armlength(3,0.0);

    double maxarm = arm1 > arm2 ? arm1 : arm2;
    maxarm = maxarm > arm3 ? maxarm : arm3;

    return(1.10 * maxarm);
}

Noise *stdproofnoise(LISA *lisa,double stproof, double sdproof) {
    // create InterpolateNoise objects for proof-mass noises
    // we need quadruple retardations for the V's appearing in the z's

    double pbtproof = 4.0 * lighttime(lisa);

    return new InterpolateNoise(stproof, pbtproof, sdproof, -2.0);
}

Noise *stdopticalnoise(LISA *lisa,double stshot, double sdshot) {
    // create InterpolateNoise objects for optical-path noises
    // we need only triple retardations for the shot's appearing in the y's

    double pbtshot = 3.0 * lighttime(lisa);

    return new InterpolateNoise(stshot, pbtshot, sdshot, 2.0);
}

Noise *stdlasernoise(LISA *lisa,double stlaser, double sdlaser, double claser) {
    // create laser noise objects
    // quadruple retardations are needed for the C's

    double pbtlaser = 4.0 * lighttime(lisa);

    return new ExpGaussNoise(stlaser,pbtlaser,claser,sdlaser);
}

