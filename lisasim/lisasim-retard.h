#ifndef _LISASIM_RETARD_H_
#define _LISASIM_RETARD_H_

#include <map>

#include "lisasim-lisa.h"

/// define retardation key

class retardpos {
 public:
    double time;
    long arms;

    retardpos(double t,long a) : time(t), arms(a) {};
};

/// define weak comparison operator for retardation keys

class poscomp {
 public:
    bool operator()(const retardpos &a1,const retardpos &a2) {
	if (a1.time != a2.time)  // maybe should use fabs
	    return a1.time < a2.time;
	else
	    return a1.arms < a2.arms;
    }
};

/// define retardation value

class retardval {
 public:
    double rt, trb, tra;

    retardval(): rt(0.0), trb(0.0), tra(0.0) {};

    retardval(double r,double b,double a): rt(r), trb(b), tra(a) {};
};

typedef map<retardpos,retardval,poscomp> retardmap;

class CacheLISA : public LISA {
 private:
    double it;  // initial time
    double rt;  // retarded time so far

    double trb; // baseline retardation
    double tra; // additional accurate retardation

    long rts;   // all retardations so far

    LISA *basiclisa;

    double buffertime, oldesttime;

    retardmap *mymap;

 public:
    CacheLISA(LISA *l,double bt) : basiclisa(l), buffertime(bt), oldesttime(0.0) {
	mymap = new retardmap();
    }

    ~CacheLISA() {
	delete mymap;
    }

    void reset() {
	mymap->clear();

	oldesttime = 0.0;

	basiclisa->reset();
    }

    LISA *thislisa() {
	return basiclisa;
    }

    double armlength(int arm, double t) {
	return basiclisa->armlength(arm,t);
    }

    double armlengthbaseline(int arm, double t) {
	return basiclisa->armlengthbaseline(arm,t);
    }

    double armlengthaccurate(int arm, double t) {
	return basiclisa->armlengthaccurate(arm,t);
    }

    void putn(Vector &n, int arm, double t) {
	basiclisa->putn(n,arm,t);
    }
        
    void putp(Vector &p, int craft, double t) {
	basiclisa->putp(p,craft,t);
    }

    void putn(double n[], int arm, double t) {
	basiclisa->putn(n,arm,t);
    }
        
    void putp(double p[], int craft, double t) {
	basiclisa->putp(p,craft,t);
    }

    void newretardtime(double t) {
	it = t; rt = t;
	trb = 0.0; tra = 0.0;
	rts = 0;

	// clean buffer; alternative algorithm

	if(t != oldesttime) {
	    oldesttime = t;
	    mymap->clear();
	}

	return;

	// clean buffer; simplest algorithm

	if(t - buffertime > oldesttime) {
	    oldesttime = t - buffertime;

	    for (retardmap::iterator iter = mymap->begin(); iter != mymap->end(); /* empty */) {
		double thistime = (*iter).first.time;

		if( thistime < oldesttime )
		    mymap->erase(iter++); // must use postfix increment
		else
		    break; // was ++iter; but we are betting on ordered traversal
	    }
	}
    }

    double retardedtime() {
	return rt;
    }

    void retard(int ret) {
	if (ret != 0) {
	    rts = (rts << 3) | (ret > 0 ? ret : (4 | -ret));
	    retardpos rpos(it,rts);

	    retardmap::iterator iter = mymap->find(rpos);

	    if(iter == mymap->end()) { // nah, don't have it...
		trb += basiclisa->armlengthbaseline(ret,rt);
		tra += basiclisa->armlengthaccurate(ret,rt);
		
		rt = (it - trb) - tra;

		retardval rval(rt,trb,tra);

		(*mymap)[rpos] = rval;
	    } else {                // yup, got it!
		rt  = (*iter).second.rt;

		trb = (*iter).second.trb;
		tra = (*iter).second.tra;
	    }
	}
    }

    void retard(LISA *anotherlisa,int ret) {
	if (ret != 0) {
	    trb += anotherlisa->armlengthbaseline(ret,rt);  
	    tra += anotherlisa->armlengthaccurate(ret,rt);

	    rt = (it - trb) - tra;
	}
    }
};

#endif /* _LISASIM_RETARD_H_ */
