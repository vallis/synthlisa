#ifndef _LISASIM_TDINOISE_H_
#define _LISASIM_TDINOISE_H_

#include "lisasim-lisa.h"
#include "lisasim-noise.h"

class TDInoise {
    private:
        LISA *lisa;

        InterpolateNoise *pm[4], *pms[4];
        
        // I label shot noises by sending and receiving spacecraft, not by link and receiving
        
        InterpolateNoise *shot[4][4];

        ExpGaussNoise *c[4], *cs[4];
        
    public:

        // claser is a correlation e-folding time

        TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser);
        ~TDInoise();

        // leave these here so we can show the cancellation of laser noise

        double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
        double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);

        double X(double t);
        double Y(double t);
        double Z(double t);

        double alpha(double t);
        double beta(double t);
        double gamma(double t);

        double zeta(double t);

        double P(double t);
        
        double E(double t);

        double Xm(double t);
};

#endif /* _LISASIM_TDINOISE_H_ */
