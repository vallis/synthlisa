#ifndef _LISASIM_TDI_H_
#define _LISASIM_TDI_H_

#include "lisasim-lisa.h"
#include "lisasim-wave.h"

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

};

class retardtime {
 private:
    LISA *lisa;

    double it, rt;
    double trb, tra;

 public:
    retardtime(LISA *l, double t) : lisa(l), it(t), rt(t), trb(0.0), tra(0.0) {};
    
    double retardedtime() {
	return rt;
    };

    void retard(int ret) {
	if (ret != 0) {
	    trb += lisa->armlengthbaseline(ret,rt);  
	    tra += lisa->armlengthaccurate(ret,rt);

	    rt = (it - trb) - tra;
	}
    };

    void retard(LISA *anotherlisa, int ret) {
	if (ret != 0) {
	    trb += anotherlisa->armlengthbaseline(ret,rt);  
	    tra += anotherlisa->armlengthaccurate(ret,rt);

	    rt = (it - trb) - tra;
	}
    };
};

#endif /* _LISASIM_TDI_H_ */
