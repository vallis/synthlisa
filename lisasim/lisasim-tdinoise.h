#ifndef _LISASIM_TDINOISE_H_
#define _LISASIM_TDINOISE_H_

#include "lisasim-tdi.h"
#include "lisasim-lisa.h"
#include "lisasim-noise.h"

class TDInoise : public TDI {
 protected:
    LISA *lisa, *phlisa;

    Noise *pm[4], *pms[4];
        
    // I label shot noises by sending and receiving spacecraft, not by link and receiving
        
    Noise *shot[4][4];

    Noise *c[4], *cs[4];

    // set this to one if we are allocating noise objects

    int allocated;
    
 public:
    // standard noises for everybody, same levels

    TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser);

    // provide arrays of noise parameters
    // the convention is {1, 1*, 2, 2*, 3, 3*}, and {12,21,23,32,31,13} (sending and receiving)

    TDInoise(LISA *mylisa, double *stproof, double *sdproof, double *stshot, double *sdshot, double *stlaser, double *sdlaser);

    // provide arrays of pointers to noise objects

    TDInoise(LISA *mylisa, Noise *proofnoise[6],Noise *shotnoise[6],Noise *lasernoise[6]);

    // the destructor will delete all the noise objects only if they were created by the constructor

    virtual ~TDInoise();

    // change the physical LISA

    void setphlisa(LISA *mylisa);

    // lock all the laser noises to one of them; use negative "master" for starred lasers

    void lock(int master);

    // reset all noises

    void reset();

    // basic TDI observables

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
    double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

class TDIaccurate : public TDInoise {
 public:
    TDIaccurate(LISA *mylisa, Noise *proofnoise[6],Noise *shotnoise[6],Noise *lasernoise[6]) : TDInoise(mylisa,proofnoise,shotnoise,lasernoise) {};
    
    ~TDIaccurate() {};

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

// return approx lighttime, for estimation of noise buffer size

extern double lighttime(LISA *lisa);

// return standard elementary noises

extern Noise *stdproofnoise(LISA *lisa,double stproof, double sdproof);
extern Noise *stdopticalnoise(LISA *lisa,double stshot, double sdshot);

extern Noise *stdlasernoise(LISA *lisa,double stlaser, double sdlaser);
extern Noise *newstdlasernoise(LISA *lisa,double stlaser, double sdlaser, int window);

extern TDInoise *stdnoise(LISA *mylisa);

extern void retardone(LISA *lisa,int ret,double t,double *retardedtime,double *totalretardbaseline,double *totalretardaccurate);

#endif /* _LISASIM_TDINOISE_H_ */
