#ifndef _LISASIM_RETARD_H_
#define _LISASIM_RETARD_H_

#include "lisasim-lisa.h"

/** Length of the cache buffer. It should be at least 2^(3 + maximum
    number of retardations - 1). */

static const unsigned long buflength = 1024;

/** Cached LISA class. Much like NoisyLISA, this class works by
    routing armlength and putp-putn calls to a "basiclisa" object that
    is passed upon construction. However, CacheLISA interposes a layer
    of its own making for "retard" calls: it maintains a cache of
    retardations for the most recently accessed time, and will return
    cached information if it has it. */

class CacheLISA : public LISA {
 private:
    /// Initial reference time for retardation.
    double it;  

    /// Retarded time so far.
    double rt;  

    /// Cumulative baseline retardation.
    double trb; 

    /// Cumulative additional accurate retardation.
    double tra; 

    /** Retardation key and hashed retardation key. The key is coded
	by expressing each retarding arm as a 3-bit signed integer,
	and shifting to the right as retarding arms are added. The
	hashed key is coded by doing the first retarding arm as for
	the key, but then coding each additional retarding arm as 1 if
	it is equal (in absolute value) to the preceding retarding
	arm, and 0 otherwise; therefore the hashed key is shifted to
	the right only by 1 before adding each new retarding arm. This
	hashed coding should be unique for causal retardations
	(corresponding to geometric paths), where it is essentially
	equivalent to a cw-ccw coding. */

    unsigned long rts, hash;

    /// Keeps track of the last retardation (needed to form hash).
    int lastarm;

    /// Basic LISA object.
    LISA *basiclisa;

    /// Cache keys.
    unsigned long keys[buflength];

    /// Caches for it, rt, trb, and tra.
    double its[buflength], rtis[buflength], trbs[buflength], tras[buflength];

 public:
    /// Default constructor. Sets caches and counters to zero.
    CacheLISA(LISA *l) : basiclisa(l) {
	for(int i=0;i<buflength;i++) {
	    keys[i] = 0;
	    its[i] = 0.0; rtis[i] = 0.0; trbs[i] = 0.0; tras[i] = 0.0;
	}

	newretardtime(0.0);
    }

    /// Reset function. Sets caches and counters to zero, and calls
    /// reset for basiclisa.
    void reset() {
	for(int i=0;i<buflength;i++) {
	    keys[i] = 0;
	    its[i] = 0.0; rtis[i] = 0.0; trbs[i] = 0.0; tras[i] = 0.0;
	}

	newretardtime(0.0);

	basiclisa->reset();
    }

    // The following is all standard for "encapsulating" LISA objects.

    /// Returns underlying PHYSICAL LISA object.
    LISA *thislisa() {
	return basiclisa->thislisa();
    }

    /// Calls underlying LISA for armlength
    double armlength(int arm, double t) {
	return basiclisa->armlength(arm,t);
    }

    /// Calls underlying LISA for armlengthbaseline
    double armlengthbaseline(int arm, double t) {
	return basiclisa->armlengthbaseline(arm,t);
    }

    /// Calls underlying LISA for armlengthaccurate
    double armlengthaccurate(int arm, double t) {
	return basiclisa->armlengthaccurate(arm,t);
    }

    /// Calls underlying LISA for putn
    void putn(Vector &n, int arm, double t) {
	basiclisa->putn(n,arm,t);
    }

    /// Calls underlying LISA for putp        
    void putp(Vector &p, int craft, double t) {
	basiclisa->putp(p,craft,t);
    }

    /// Calls underlying LISA for putn[]
    void putn(double n[], int arm, double t) {
	basiclisa->putn(n,arm,t);
    }

    /// Calls underlying LISA for putp[]
    void putp(double p[], int craft, double t) {
	basiclisa->putp(p,craft,t);
    }

    // The following is specific to CacheLISA.

    /// Sets up counters for a new retardation.
    void newretardtime(double t) {
	it = t; rt = t;
	trb = 0.0; tra = 0.0;
	rts = 0; hash = 0;
	lastarm = 0;
    }

    /// Returns the retarded time so far.
    double retardedtime() {
	return rt;
    }

    /// Computes a retardation. Contains the caching mechanism.
    void retard(int ret) {
	if(ret == 0)
	    return;

	// Update the retardation key and the hash.
	if(rts == 0) {
	    rts = (rts << 3) | (ret > 0 ? ret : (4 | -ret));
	    hash = rts;
	} else {
	    rts = (rts << 3) | (ret > 0 ? ret : (4 | -ret));
	    
	    hash = hash << 1;
	    if(abs(ret) == abs(lastarm))
		rts = rts | 1;
	}
	
	lastarm = ret;

	// We cannot handle such as a large hash. Should throw exception.
	if(hash >= buflength) {
	    cout << "CacheLISA: Hash out of range in retard()!" << endl;
	    abort();
	}

	if(its[hash] == it && keys[hash] == rts) {
	    // Got it!

	    rt = rtis[hash];
	    trb = trbs[hash]; tra = tras[hash];
	} else {
	    // Nah, will have to compute it.

	    trb += basiclisa->armlengthbaseline(ret,rt);
	    tra += basiclisa->armlengthaccurate(ret,rt);
		
	    rt = (it - trb) - tra;
	}

	if(its[hash] < it) {
	    // We can replace the old value in the cache with the
	    // newly computed value. But we don't do this if the time
	    // is the same, but the key is different (a cache collision).

	    its[hash] = it;
	    keys[hash] = rts;

	    rtis[hash] = rt;
	    trbs[hash] = trb; tras[hash] = tra;
	}
    }

    /** Computes a retardation with a different LISA geometry
	[generally it will be basiclisa->thislisa()]. Note that
	calling this form of retard will (probably) make the cache
	state undefined, so we should never do it unless it's the last
	retard call */

    void retard(LISA *anotherlisa,int ret) {
	if (ret != 0) {
	    trb += anotherlisa->armlengthbaseline(ret,rt);  
	    tra += anotherlisa->armlengthaccurate(ret,rt);

	    rt = (it - trb) - tra;
	}
    }
};

#endif /* _LISASIM_RETARD_H_ */
