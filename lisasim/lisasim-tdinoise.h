#ifndef _LISASIM_TDINOISE_H_
#define _LISASIM_TDINOISE_H_

#include "lisasim-tdi.h"
#include "lisasim-lisa.h"
#include "lisasim-noise.h"

class TDInoise : public TDI {
    private:
        LISA *lisa, *phlisa;

        Noise *pm[4], *pms[4];
        
        // I label shot noises by sending and receiving spacecraft, not by link and receiving
        
        Noise *shot[4][4];

        Noise *c[4], *cs[4];

	// set this to one if we are allocating noise objects

	int allocated;

    public:
        // standard noises for everybody, same levels

        TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser);

	// provide arrays of noise parameters
	// the convention is {1, 1*, 2, 2*, 3, 3*}, and {12,21,23,32,31,13} (sending and receiving)

        TDInoise(LISA *mylisa, double *stproof, double *sdproof, double *stshot, double *sdshot, double *stlaser, double *sdlaser, double *claser);

	// provide arrays of pointers to noise objects

	TDInoise(LISA *mylisa, Noise **proofnoise,Noise **shotnoise,Noise **lasernoise);

	// the destructor will delete all the noise objects only if they were created by the constructor

        ~TDInoise();

	// change the physical LISA

	void setphlisa(LISA *mylisa);

	// reset all noises

        void reset();

	// basic TDI observables

        double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
        double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);
};

// return approx lighttime, for estimation of noise buffer size

extern double lighttime(LISA *lisa);

// return standard elementary noises

extern Noise *stdproofnoise(LISA *lisa,double stproof, double sdproof);
extern Noise *stdopticalnoise(LISA *lisa,double stshot, double sdshot);
extern Noise *stdlasernoise(LISA *lisa,double stlaser, double sdlaser, double claser);

#endif /* _LISASIM_TDINOISE_H_ */
