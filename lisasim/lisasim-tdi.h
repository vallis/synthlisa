#ifndef _LISASIM_TDI_H_
#define _LISASIM_TDI_H_

#include "lisasim-lisa.h"
#include "lisasim-wave.h"

class TDI {
    public:
        LISA *lisa;
        Wave *wave;

        TDI(LISA *mylisa, Wave *mywave) {
            lisa = mylisa;
            wave = mywave;
        }

        double M(double t);
        double N(double t);
        double O(double t);
    
        double X(double t);
        double Y(double t);
        double Z(double t);
    
        double alpha(double t);
        double beta(double t);
        double gamma(double t);
    
        double zeta(double t);

	double P(double t);
	double E(double t);
	double U(double t);

    private:
	double psi(int arm, double t, double twave);
//        double psi(int arm, double t);
    
        double retard(int craft, double t);
    
        double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
};

#endif /* _LISASIM_TDI_H_ */
