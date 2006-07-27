/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#ifndef _LISASIM_TDI_H_
#define _LISASIM_TDI_H_

#include "lisasim-lisa.h"
#include "lisasim-wave.h"

#include <math.h>

class TDI {
 public:
    TDI() {};
    virtual ~TDI() {};

    virtual void reset() {};

    virtual double X(double t);
    virtual double Y(double t);
    virtual double Z(double t);
    
    virtual double alpha(double t);
    virtual double beta(double t);
    virtual double gamma(double t);

    virtual double alpham(double t);
    virtual double betam(double t);
    virtual double gammam(double t);

    virtual double alpha1(double t);
    virtual double alpha2(double t);
    virtual double alpha3(double t);

    virtual double zeta(double t);

    virtual double zeta1(double t);
    virtual double zeta2(double t);
    virtual double zeta3(double t);

    virtual double P(double t);
    virtual double E(double t);
    virtual double U(double t);
    
    virtual double Xm(double t);
    virtual double Ym(double t);
    virtual double Zm(double t);

    virtual double Xmlock1(double t);
    virtual double Xmlock2(double t);
    virtual double Xmlock3(double t);

    virtual double X1(double t);
    virtual double X2(double t);
    virtual double X3(double t);

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t) {return 0.0;};
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t) {return 0.0;};

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t) {return 0.0;};
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t) {return 0.0;};

    double y123(double t) {return y(1,2,3,0,0,0,0,0,0,0,t);};
    double y231(double t) {return y(2,3,1,0,0,0,0,0,0,0,t);};
    double y312(double t) {return y(3,1,2,0,0,0,0,0,0,0,t);};

    double y321(double t) {return y(3,-2,1,0,0,0,0,0,0,0,t);};
    double y132(double t) {return y(1,-3,2,0,0,0,0,0,0,0,t);};
    double y213(double t) {return y(2,-1,3,0,0,0,0,0,0,0,t);};
    
    double z123(double t) {return z(1,2,3,0,0,0,0,0,0,0,0,t);};
    double z231(double t) {return z(2,3,1,0,0,0,0,0,0,0,0,t);};
    double z312(double t) {return z(3,1,2,0,0,0,0,0,0,0,0,t);};
    
    double z321(double t) {return z(3,-2,1,0,0,0,0,0,0,0,0,t);};
    double z132(double t) {return z(1,-3,2,0,0,0,0,0,0,0,0,t);};
    double z213(double t) {return z(2,-1,3,0,0,0,0,0,0,0,0,t);};
};

class TDIquantize : public TDI {
 private:
    TDI *basetdi;

    double qlevel, qover;

    double quantize(double var) {
	   if(fabs(var) > qover) {
	       // hey, there is no "sign" function in C/C++ 
	       return (var > 0.0) ? qover : -qover;
	   } else {
	       return qlevel * (long long)(var/qlevel);
	   }
    };

 public:
    TDIquantize(TDI *bt,double qlev,int qbits,int qsat)
    	: basetdi(bt), qlevel(qlev/pow(2.0,qbits)), qover(qlev*pow(2.0,qsat)) {};
    
    virtual ~TDIquantize() {};

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t) {
    	return quantize(basetdi->y(send, link, recv, ret1, ret2, ret3, t));
    };

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t) {
    	return quantize(basetdi->y(send, link, recv, ret1, ret2, ret3, ret4, ret5, ret6, ret7, t));
    };

    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t) {
    	return quantize(basetdi->z(send, link, recv, ret1, ret2, ret3, ret4, t));
    };

    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t) {
    	return quantize(basetdi->z(send, link, recv, ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, t));
    };
};

class TDIalpha : public Signal {
 private:
    TDI *tdi;

 public:
    TDIalpha(TDI *t) : tdi(t) {};
    
    double value(double t) { return tdi->alpha(t); };
};

class TDIbeta : public Signal {
 private:
    TDI *tdi;

 public:
    TDIbeta(TDI *t) : tdi(t) {};
    
    double value(double t) { return tdi->beta(t); };
};

class TDIgamma : public Signal {
 private:
    TDI *tdi;

 public:
    TDIgamma(TDI *t) : tdi(t) {};
    
    double value(double t) { return tdi->gamma(t); };
};

class SampledTDI : public TDI {
 private:
    LISA *lisa;

    Signal *yobj[4][4], *zobj[4][4];

 public:
    SampledTDI(LISA *lisa,Noise *yijk[6],Noise *zijk[6]);
    ~SampledTDI() {};

    void reset(unsigned long seed = 0);

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
    double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

#endif /* _LISASIM_TDI_H_ */
