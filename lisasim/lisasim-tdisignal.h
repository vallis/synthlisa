#ifndef _LISASIM_TDISIGNAL_H_
#define _LISASIM_TDISIGNAL_H_

#include "lisasim-tdi.h"
#include "lisasim-lisa.h"
#include "lisasim-wave.h"

class TDIsignal : public TDI {
    private:
        LISA *lisa;
        Wave *wave;
        
	double psi(int arm, double t, double twave);
        
        double retard(int craft, double t);
        
    public:
        TDIsignal(LISA *mylisa, Wave *mywave);

        double M(double t);
        double N(double t);
        double O(double t);

        double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
};

#endif /* _LISASIM_TDISIGNAL_H_ */
