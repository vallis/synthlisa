#ifndef _LISASIM_TDINOISE_H_
#define _LISASIM_TDINOISE_H_

#include "lisasim-tdi.h"
#include "lisasim-lisa.h"
#include "lisasim-noise.h"

class TDInoise : public TDI {
    private:
        LISA *lisa, *phlisa;

        InterpolateNoise *pm[4], *pms[4];
        
        // I label shot noises by sending and receiving spacecraft, not by link and receiving
        
        InterpolateNoise *shot[4][4];

        ExpGaussNoise *c[4], *cs[4];

        void initialize(double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser);
        
    public:
        // the first constructor sets mylisa = physlisa; the real construction is done by initialize

        TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser);
        TDInoise(LISA *mylisa, LISA *physlisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser);

        ~TDInoise();

        void reset();

        double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
        double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);
};

#endif /* _LISASIM_TDINOISE_H_ */
