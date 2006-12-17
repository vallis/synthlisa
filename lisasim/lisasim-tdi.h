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

/* These objects could be returned by the main TDI class if the TDI functions are called without arguments! */
/* Could also try an implementation with method pointers... but should check whether it's slower */

class TDI;

class TDIobject : public Signal {
 protected:
   TDI *tdi;    
    
 public:
    TDIobject(TDI *t) : tdi(t) {};
    virtual ~TDIobject() {};
        
    virtual double value(double t) = 0;
};

class TDIobjectpnt : public TDIobject {
 private:
    double (TDI::*obs)(double t);
    
 public:
    TDIobjectpnt(TDI *t,double (TDI::*o)(double t)) : TDIobject(t), obs(o) {};
    ~TDIobjectpnt() {};
    
    double value(double t) { return (tdi->*obs)(t); };
};

class timeobject : public Signal {
 public:
    timeobject() {};
    
    double value(double t) { return t; };
};

class TDI {
 public:
    TDI() {};
    virtual ~TDI() {};

    virtual void reset() {};
    
    virtual double alpham(double t);
    TDIobject *alpham() { return new TDIobjectpnt(this,&TDI::alpham); };
    virtual double betam(double t);
    TDIobject *betam()  { return new TDIobjectpnt(this,&TDI::betam);  };
    virtual double gammam(double t);
    TDIobject *gammam() { return new TDIobjectpnt(this,&TDI::gammam); };

    virtual double zetam(double t);
    TDIobject *zetam()  { return new TDIobjectpnt(this,&TDI::zetam);  };

    virtual double alpha1(double t);
    TDIobject *alpha1() { return new TDIobjectpnt(this,&TDI::alpha1); };
    virtual double alpha2(double t);
    TDIobject *alpha2() { return new TDIobjectpnt(this,&TDI::alpha2); };
    virtual double alpha3(double t);
    TDIobject *alpha3() { return new TDIobjectpnt(this,&TDI::alpha3); };

    virtual double zeta1(double t);
    TDIobject *zeta1() { return new TDIobjectpnt(this,&TDI::zeta1); };
    virtual double zeta2(double t);
    TDIobject *zeta2() { return new TDIobjectpnt(this,&TDI::zeta2); };
    virtual double zeta3(double t);
    TDIobject *zeta3() { return new TDIobjectpnt(this,&TDI::zeta3); };

    // P, E, U still have non-signed delays
    
    virtual double P(double t);
    TDIobject *P() { return new TDIobjectpnt(this,&TDI::P); };
    virtual double E(double t);
    TDIobject *E() { return new TDIobjectpnt(this,&TDI::E); };
    virtual double U(double t);
    TDIobject *U() { return new TDIobjectpnt(this,&TDI::U); };
    
    virtual double Xm(double t);
    TDIobject *Xm() { return new TDIobjectpnt(this,&TDI::Xm); };
    virtual double Ym(double t);
    TDIobject *Ym() { return new TDIobjectpnt(this,&TDI::Ym); };
    virtual double Zm(double t);
    TDIobject *Zm() { return new TDIobjectpnt(this,&TDI::Zm); };

    virtual double Xmlock1(double t);
    TDIobject *Xmlock1() { return new TDIobjectpnt(this,&TDI::Xmlock1); };
    virtual double Xmlock2(double t);
    TDIobject *Xmlock2() { return new TDIobjectpnt(this,&TDI::Xmlock2); };
    virtual double Xmlock3(double t);
    TDIobject *Xmlock3() { return new TDIobjectpnt(this,&TDI::Xmlock3); };
    
    virtual double X1(double t);
    TDIobject *X1() { return new TDIobjectpnt(this,&TDI::X1); };
    virtual double X2(double t);
    TDIobject *X2() { return new TDIobjectpnt(this,&TDI::X2); };
    virtual double X3(double t);
    TDIobject *X3() { return new TDIobjectpnt(this,&TDI::X3); };

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t) { return 0.0; };
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t) { return 0.0; };

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t) { return 0.0; };
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t) { return 0.0; };

    double y123(double t) { return y(1,2,3,0,0,0,0,0,0,0,t); };
    TDIobject *y123() { return new TDIobjectpnt(this,&TDI::y123); };
    double y231(double t) { return y(2,3,1,0,0,0,0,0,0,0,t); };
    TDIobject *y231() { return new TDIobjectpnt(this,&TDI::y231); };
    double y312(double t) { return y(3,1,2,0,0,0,0,0,0,0,t); };
    TDIobject *y312() { return new TDIobjectpnt(this,&TDI::y312); };

    double y321(double t) { return y(3,-2,1,0,0,0,0,0,0,0,t); };
    TDIobject *y321() { return new TDIobjectpnt(this,&TDI::y321); };
    double y132(double t) { return y(1,-3,2,0,0,0,0,0,0,0,t); };
    TDIobject *y132() { return new TDIobjectpnt(this,&TDI::y132); };
    double y213(double t) { return y(2,-1,3,0,0,0,0,0,0,0,t); };
    TDIobject *y213() { return new TDIobjectpnt(this,&TDI::y213); };
    
    double z123(double t) { return z(1,2,3,0,0,0,0,0,0,0,0,t); };
    TDIobject *z123() { return new TDIobjectpnt(this,&TDI::z123); };
    double z231(double t) { return z(2,3,1,0,0,0,0,0,0,0,0,t); };
    TDIobject *z231() { return new TDIobjectpnt(this,&TDI::z231); };
    double z312(double t) { return z(3,1,2,0,0,0,0,0,0,0,0,t); };
    TDIobject *z312() { return new TDIobjectpnt(this,&TDI::z312); };
    
    double z321(double t) { return z(3,-2,1,0,0,0,0,0,0,0,0,t); };
    TDIobject *z321() { return new TDIobjectpnt(this,&TDI::z321); };
    double z132(double t) { return z(1,-3,2,0,0,0,0,0,0,0,0,t); };
    TDIobject *z132() { return new TDIobjectpnt(this,&TDI::z132); };
    double z213(double t) { return z(2,-1,3,0,0,0,0,0,0,0,0,t); };
    TDIobject *z213() { return new TDIobjectpnt(this,&TDI::z213); };
    
    double time(double t) { return t;};
    timeobject *time() { return new timeobject(); };
    double t(double t)    { return t;};
    timeobject *t()    { return new timeobject(); };
};

extern void fastgetobs(double *buffer,long length,long samples,double stime,Signal **thesignals,int signals,double inittime);
extern void fastgetobsc(double *buffer,long length,long samples,double stime,Signal **thesignals,int signals,double inittime);

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

class SampledTDI : public TDI {
 protected:
    LISA *lisa;

    Signal *yobj[4][4], *zobj[4][4];

 public:
    SampledTDI(LISA *lisa,Noise *yijk[6],Noise *zijk[6]);
    ~SampledTDI() {};

    void reset(unsigned long seed = 0);

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
    double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

class SampledTDIaccurate : public SampledTDI {
 public:
    SampledTDIaccurate(LISA *lisa,Noise *yijk[6],Noise *zijk[6]) : SampledTDI(lisa,yijk,zijk) {};
    ~SampledTDIaccurate() {};

    void reset(unsigned long seed = 0);

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

#endif /* _LISASIM_TDI_H_ */
