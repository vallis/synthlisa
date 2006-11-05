/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#ifndef _LISASIM_RETARD_H_
#define _LISASIM_RETARD_H_

#include "lisasim-tens.h"
#include "lisasim-lisa.h"
#include "lisasim-except.h"

#include <iostream>

// ??? Might be good to move all these function definitions in a cpp file

/** Length of the cache buffer. It should be at least 2^(3 + maximum
    number of retardations). */

static const unsigned long buflength = 2048;

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

    /// p keys and Cache
    
    double pts[buflength];
    int pis[buflength];
    
    Vector ps[buflength];

 public:
    /// Default constructor. Sets caches and counters to zero.
    CacheLISA(LISA *l);

    /// Reset function. Sets caches and counters to zero, and calls reset for basiclisa.
    void reset();

    // The following is all standard for "encapsulating" LISA objects.

    LISA *physlisa() { return basiclisa->physlisa(); };

    double armlength(int arm, double t) { return basiclisa->armlength(arm,t); };

    double armlengthbaseline(int arm, double t) { return basiclisa->armlengthbaseline(arm,t); };
    double armlengthaccurate(int arm, double t) { return basiclisa->armlengthaccurate(arm,t); };

    void putn(Vector &n, int arm, double t) { basiclisa->putn(n,arm,t); };

    void putp(Vector &p, int craft, double t);
    void putp(LISA *anotherlisa,Vector &p, int craft, double t);

    // The following is specific to CacheLISA.

    /// Sets up counters for a new retardation.
    void newretardtime(double t);
    
    /// Returns the retarded time and retardation so far.
    double retardedtime();
    double retardation();

    /// Computes a retardation. Contains the caching mechanism.
    void retard(int ret);

    /** Computes a retardation with a different LISA geometry
	[generally it will be basiclisa->physlisa()]. Note that
	calling this form of retard will (probably) make the cache
	state undefined, so we should never do it unless it's the last
	retard call */

    void retard(LISA *anotherlisa,int ret);
};

#endif /* _LISASIM_RETARD_H_ */
