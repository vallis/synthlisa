#include "lisasim-tdinoise.h"

#include <time.h>

// this version takes the parameters of the basic noises and lets us allocate objects as needed

TDInoise::TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser) {
    phlisa = mylisa->thislisa();
    lisa = mylisa;

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
	c[craft] = stdlasernoise(lisa,stlaser,sdlaser);
        cs[craft] = stdlasernoise(lisa,stlaser,sdlaser);
    }

    allocated = 1;
}

// this version takes arrays of basic-noise parameters, allowing for different noises on different objects,
// and lets us allocate objects as needed

TDInoise::TDInoise(LISA *mylisa, double *stproof, double *sdproof, double *stshot, double *sdshot, double *stlaser, double *sdlaser) {
    phlisa = mylisa->thislisa();
    lisa = mylisa;

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
        c[craft]  = stdlasernoise(lisa,stlaser[2*(craft-1)],  sdlaser[2*(craft-1)]);
        cs[craft] = stdlasernoise(lisa,stlaser[2*(craft-1)+1],sdlaser[2*(craft-1)+1]);
    }

    allocated = 1;
}

// this version takes pointers to noise objects, allowing for user-specified noises on different objects

TDInoise::TDInoise(LISA *mylisa, Noise *proofnoise[6],Noise *shotnoise[6],Noise *lasernoise[6]) {
    phlisa = mylisa->thislisa();
    lisa = mylisa;

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

class zLockNoise : public Noise {
 private:
    int slave;

    Noise *masterpm, *slavepm, *masterc, *slavec;

 public:
    zLockNoise(int recv,Noise *mpm,Noise *spm,Noise *mc,Noise *sc)
	: slave(recv), masterpm(mpm), slavepm(spm), masterc(mc), slavec(sc) {};

    virtual ~zLockNoise() {
	delete slavec;
    };

    double operator[](double t) {
	if (slave > 0) {
	    // in this case master is starred, slave is not
	    return (*masterc)[t] - (*masterpm)[t] - (*slavepm)[t];
	} else {
	    // in this case master is not starred, slave is
	    return (*masterc)[t] + (*masterpm)[t] + (*slavepm)[t];
	}
    };

    double noise(double time) {
	return (*this)[time];
    };
};

class yLockNoise : public Noise {
 private:
    int slave, arm;

    LISA *lisa;

    Noise *slavepm, *shot, *masterc, *slavec;

 public:
    yLockNoise(int recv,int link,LISA *l,Noise *spm,Noise *sh,Noise *mc,Noise *sc)
	: slave(recv), arm(link), lisa(l), slavepm(spm), shot(sh), masterc(mc), slavec(sc) {};

    virtual ~yLockNoise() {
	delete slavec;
    };

    double operator[](double t) {
	if (slave > 0) {
	    // since the slave is not starred, the link is positive (cyclic)
	    return (*masterc)[t - lisa->armlength(arm,t)] - 2.0 * (*slavepm)[t] + (*shot)[t];
	} else {
	    // the slave is starred, the link is anticyclic
	    return (*masterc)[t - lisa->armlength(arm,t)] + 2.0 * (*slavepm)[t] + (*shot)[t];
	}
    };

    double noise(double time) {
	return (*this)[time];
    };
};

// locking procedure: use negative numbers for starred lasers

void TDInoise::lock(int master) {
    int mastera = abs(master);
    int slaveb = (mastera % 3) + 1;
    int slavec = (slaveb % 3) + 1;

    // first lock the laser on the same bench

    if(master > 0) {
	cs[mastera] = new zLockNoise(-mastera,pm[ mastera],pms[mastera],c[ mastera],cs[mastera]);
    } else {
	c[ mastera] = new zLockNoise( mastera,pms[mastera],pm[ mastera],cs[mastera],c[ mastera]);
    }

    // now lock across to the other benches

    cs[slaveb] = new yLockNoise(-slaveb,-slavec,phlisa,pms[slaveb],shot[mastera][slaveb],c[ mastera],cs[slaveb]);
    c[ slavec] = new yLockNoise( slavec, slaveb,phlisa,pm[ slavec],shot[mastera][slavec],cs[mastera],c[ slavec]);

    // finally, lock the lasers on the back of the other benches

    c[ slaveb] = new zLockNoise( slaveb,pms[slaveb],pm[slaveb],cs[slaveb],c[slaveb]);
    cs[slavec] = new zLockNoise(-slavec,pm[slavec],pms[slavec],c[slavec],cs[slavec]);
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
	
	// remove laser-noise InterpolateNoise objects

	for(int craft = 1; craft <= 3; craft++) {
	    if(c[craft])  {delete c[craft]; c[craft]=0;}
	    if(cs[craft]) {delete cs[craft]; cs[craft]=0;}
	}
    }
}

void TDInoise::reset() {
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

// this is a debugging function, which appears in lisasim-swig.i

double retardation(LISA *lisa,int ret1,int ret2,int ret3,int ret4,int ret5,int ret6,int ret7,int ret8,double t) {
    lisa->newretardtime(t);

    lisa->retard(ret8); lisa->retard(ret7); lisa->retard(ret6); lisa->retard(ret5);
    lisa->retard(ret4); lisa->retard(ret3); lisa->retard(ret2); lisa->retard(ret1);

    return lisa->retardedtime();
}

double TDInoise::y(int send, int slink, int recv, int ret1, int ret2, int ret3, double t) {
    return y(send,slink,recv,ret1,ret2,ret3,0,0,0,0,t);
}

double TDInoise::y(int send, int slink, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t) {
    int link = abs(slink);

    // this recursive retardation procedure assumes smart TDI...

    lisa->newretardtime(t);

    lisa->retard(ret7); lisa->retard(ret6); lisa->retard(ret5);
    lisa->retard(ret4); lisa->retard(ret3); lisa->retard(ret2); lisa->retard(ret1);

    double retardedtime = lisa->retardedtime();

    if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2)) {
        // cyclic combination
        // if introducing error in the determination of the armlengths, it should not enter
        // the following (physical) retardation of the laser noise, so we use the phlisa object

	lisa->retard(phlisa,link);
	double retardlaser = lisa->retardedtime();

        return( (*cs[send])[retardlaser] - 2.0 * (*pm[recv])[retardedtime]  - (*c[recv])[retardedtime]  + 
                (*shot[send][recv])[retardedtime] );
    } else {
        // anticyclic combination
        // ditto here

	lisa->retard(phlisa,-link);
	double retardlaser = lisa->retardedtime();

        return( (*c[send])[retardlaser]  + 2.0 * (*pms[recv])[retardedtime] - (*cs[recv])[retardedtime] +
                (*shot[send][recv])[retardedtime] );
    }
}

double TDInoise::z(int send, int slink, int recv, int ret1, int ret2, int ret3, int ret4, double t) {
    return z(send,slink,recv,ret1,ret2,ret3,ret4,0,0,0,0,t);
}

double TDInoise::z(int send, int slink, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t) {
    int link = abs(slink);

    // this recursive retardation procedure assumes smart TDI...
    // (and the correct order in the retardation expressions)

    lisa->newretardtime(t);

    lisa->retard(ret8); lisa->retard(ret7); lisa->retard(ret6); lisa->retard(ret5);
    lisa->retard(ret4); lisa->retard(ret3); lisa->retard(ret2); lisa->retard(ret1);

    double retardedtime = lisa->retardedtime();

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
    // we need quadruple retardations for the V's appearing in the z's (octuple for 2nd-gen TDI)

    double pbtproof = 8.0 * lighttime(lisa);

    return new InterpolateNoise(stproof, pbtproof, sdproof, -2.0);
}

Noise *stdopticalnoise(LISA *lisa,double stshot, double sdshot) {
    // create InterpolateNoise objects for optical-path noises
    // we need only triple retardations for the shot's appearing in the y's (septuple for 2nd-gen TDI)

    double pbtshot = 7.0 * lighttime(lisa);

    return new InterpolateNoise(stshot, pbtshot, sdshot, 2.0);
}

Noise *stdlasernoise(LISA *lisa,double stlaser, double sdlaser) {
    // create laser noise objects
    // quadruple retardations are needed for the C's (octuple for 2nd-gen TDI)

    double pbtlaser = 8.0 * lighttime(lisa);

    return new InterpolateNoise(stlaser, pbtlaser, sdlaser, 0.0);
}

TDInoise *stdnoise(LISA *mylisa) {
    return new TDInoise(mylisa,1.0,2.5e-48,1.0,1.8e-37,1.0,1.1e-26);
}
