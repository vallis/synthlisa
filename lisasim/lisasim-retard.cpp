/* $Id: $
 * $Date: $
 * $Author: $
 * $Revision: $
 */

#include "lisasim-retard.h"

CacheLISA::CacheLISA(LISA *l) : basiclisa(l) {
   for(unsigned int i=0;i<buflength;i++) {
       keys[i] = 0;
       its[i] = 0.0; rtis[i] = 0.0; trbs[i] = 0.0; tras[i] = 0.0;

       pts[i] = 0.0; pis[i] = 0;
       ps[i][0] = 0.0; ps[i][1] = 0.0; ps[i][2] = 0.0;
   }

   newretardtime(0.0);
}

void CacheLISA::reset() {
   for(unsigned int i=0;i<buflength;i++) {
       keys[i] = 0;
       its[i] = 0.0; rtis[i] = 0.0; trbs[i] = 0.0; tras[i] = 0.0;
   }

   newretardtime(0.0);

   basiclisa->reset();
}

void CacheLISA::newretardtime(double t) {
   it = t; rt = t;
   trb = 0.0; tra = 0.0;
   rts = 0; hash = 0;
   lastarm = 0;
}

double CacheLISA::retardedtime() {
   return rt;
}

double CacheLISA::retardation() {
   return trb + tra;
}

void CacheLISA::retard(int ret) {
    if(ret == 0) return;

    // Update the retardation key and the hash.
    if(rts == 0) {
        rts = ret > 0 ? ret : (4 | -ret);

        // use a "guard bit" at the left side of the hash
        hash = 8 | rts;
    } else {
        rts = (rts << 3) | (ret > 0 ? ret : (4 | -ret));
    
        hash = hash << 1;
        if(ret > 0) hash = hash | 1;
    }

    // printf("ret %d, looking for hash %d... ",ret,hash);

    lastarm = ret;

    // We cannot handle such as a large hash. Should throw exception.
    if(hash >= buflength) {
        std::cerr << "CacheLISA::retard: Hash " << hash << " out of range"
                  << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;

		ExceptionOutOfBounds e;
       	throw e;
    }

    if(its[hash] == it && keys[hash] == rts) {
        // Got it!

        // printf("...found!\n");

        rt = rtis[hash];
        trb = trbs[hash]; tra = tras[hash];
    } else {
        // Nah, will have to compute it.

        // printf("...nah.\n");
        // printf("Getting arm %d at time %f\n",ret,rt);
        // printf("rts %d hash %d ",rts,hash);

        trb += basiclisa->armlengthbaseline(ret,rt);
        tra += basiclisa->armlengthaccurate(ret,rt);
	
        rt = (it - trb) - tra;

        // printf("setting hash %d\n",hash);

        // note: in case of cache collisions, we're overwriting the cache

        its[hash] = it;
        keys[hash] = rts;

        rtis[hash] = rt;
        trbs[hash] = trb; tras[hash] = tra;
    }
}

void CacheLISA::retard(LISA *anotherlisa,int ret) {
    if (anotherlisa == basiclisa) {
        retard(ret);
    } else if (ret != 0) {
        trb += anotherlisa->armlengthbaseline(ret,rt);  
        tra += anotherlisa->armlengthaccurate(ret,rt);

        rt = (it - trb) - tra;
    }
}

void CacheLISA::putp(Vector &p, int craft, double t) {
    // printf("Looking for hash %d (%f,%d)... ",hash,t,craft);

    if(pts[hash] == t && pis[hash] == craft) {
        // printf("found!\n");

        p[0] = ps[hash][0]; p[1] = ps[hash][1]; p[2] = ps[hash][2];
    } else {
        // printf("no, have %f,%d\n",pts[hash],pis[hash]);

        basiclisa->putp(p,craft,t);
        
        pts[hash] = t;
        pis[hash] = craft;
        ps[hash][0] = p[0]; ps[hash][1] = p[1]; ps[hash][2] = p[2];
    }
}

void CacheLISA::putp(LISA *anotherlisa,Vector &p, int craft, double t) {
    if (anotherlisa == basiclisa) {
        putp(p,craft,t);
    } else {
        anotherlisa->putp(p,craft,t);
    }
}
