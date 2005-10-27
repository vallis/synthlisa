#ifndef _LISASIM_LISA_EXTRA_H_
#define _LISASIM_LISA_EXTRA_H_

#include "lisasim-tens.h"
#include "lisasim-noise.h"
#include "lisasim-lisa.h"

#include <math.h>

/// Takes any LISA and adds noise to the TDI (nominal) armlengths.

class NoisyLISA : public LISA {
    private:
        LISA *cleanlisa;

        InterpolateNoise *uperror[4], *downerror[4];

    public:
        NoisyLISA(LISA *clean,double starm,double sdarm);
        ~NoisyLISA();
        
        void reset();

       	LISA *physlisa() {
	    return cleanlisa;
	}

        double armlength(int arm, double t);

	double armlengthbaseline(int arm, double t);
	double armlengthaccurate(int arm, double t);

        void putn(Vector &n, int arm, double t) {
            cleanlisa->putn(n,arm,t);
        }
        
        void putp(Vector &p, int craft, double t) {
            cleanlisa->putp(p,craft,t);
        }

        void putn(double n[], int arm, double t) {
            cleanlisa->putn(n,arm,t);
        }
        
        void putp(double p[], int craft, double t) {
            cleanlisa->putp(p,craft,t);
	}
};

/** Takes any LISA for putp, putn, and physical armlengths, but uses a
    parametrized model similar to EccentricInclined for the TDI
    (nominal) armlengths. */

class NominalLISA : public LISA {
    private:
        LISA *reallisa;

	double L, swi;

	double dL;
	double pdelmod, mdelmod, delmod3;
	double delmodph[4], delmodph2;
	double toffset;

    public:
	double Lnom, emod, cmod, toff;

	NominalLISA(double eta0,double xi0,double sw,double t0);
	~NominalLISA();

	void setparameters(double l,double cm,double em,double toff);
	void setparameters(double cm,double em,double toff);
	void setparameters3(double l,double cm,double em);

       	LISA *physlisa() {
	    return reallisa;
	}

        double armlength(int arm, double t);

	double armlengthbaseline(int arm, double t);
	double armlengthaccurate(int arm, double t);

        void putn(Vector &n, int arm, double t) {
            reallisa->putn(n,arm,t);
        }
        
        void putp(Vector &p, int craft, double t) {
            reallisa->putp(p,craft,t);
        }

        void putn(double n[], int arm, double t) {
            reallisa->putn(n,arm,t);
        }
        
        void putp(double p[], int craft, double t) {
            reallisa->putp(p,craft,t);
	}
};

// --- LinearLISA class ---------------------------------------------------

class LinearLISA : public LISA {
 private:
    LISA *reallisa;

    double L;
    double dL[6], dLdt[6];

    double toffset;

 public:
    LinearLISA(double eta0,double xi0,double sw,double t0);
    ~LinearLISA();

    void settimeoffset(double toff);
    void setparameters(double dl[6],double dldt[6]);
    
    double armlengtherror(int arm, double t);

    LISA *physlisa() {
	return reallisa;
    }

    double armlength(int arm, double t);
    
    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);

    void putn(Vector &n, int arm, double t) {
	reallisa->putn(n,arm,t);
    }
        
    void putp(Vector &p, int craft, double t) {
	reallisa->putp(p,craft,t);
    }

    void putn(double n[], int arm, double t) {
	reallisa->putn(n,arm,t);
    }
        
    void putp(double p[], int craft, double t) {
	reallisa->putp(p,craft,t);
    }
};

class LISAMakeNoise: public MakeNoise {
 private:
    LISA *mylisa;
    int myarm;

    double mynorm, mystime, myptime;

    MakeNoise *mynoise;

 public:
    LISAMakeNoise(LISA *lisa,int arm,double norm,double stime,double ptime,MakeNoise *noise)
	: mylisa(lisa), myarm(arm), mynorm(norm), mystime(stime), myptime(ptime), mynoise(noise) {};

    virtual ~LISAMakeNoise() {
	delete mynoise;
    };

    void reset() {
	mynoise->reset();
    };

    double operator[](long pos) {
	return mylisa->armlength(myarm,pos*mystime - myptime) + mynorm * (*mynoise)[pos];
    };
};

class LISANoise: public InterpolateNoise {
 public:
    LISANoise(LISA *cleanlisa,int link,double sampletime,double prebuffer,double density,int swindow)
	: InterpolateNoise(sampletime,prebuffer,density,0.0,swindow) {
	thenoise = new LISAMakeNoise(cleanlisa,link,density,sampletime,prebuffer,thenoise);
	normalize = 1.0;
    };
};

class MeasureLISA : public LISA {
 private:
    LISA *cleanlisa;

    LISANoise *uperror[4], *downerror[4];

 public:
    MeasureLISA(LISA *clean,double starm,double sdarm,int swindow = 1);
    ~MeasureLISA();
        
    void reset();

    LISA *physlisa() {
	return cleanlisa;
    }

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);

    void putn(Vector &n, int arm, double t) {
	cleanlisa->putn(n,arm,t);
    }
        
    void putp(Vector &p, int craft, double t) {
	cleanlisa->putp(p,craft,t);
    }

    void putn(double n[], int arm, double t) {
	cleanlisa->putn(n,arm,t);
    }
        
    void putp(double p[], int craft, double t) {
	cleanlisa->putp(p,craft,t);
    }
};

#endif /* _LISASIM_LISA_EXTRA_H_ */






