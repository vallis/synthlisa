#ifndef _LISASIM_LISA_H_
#define _LISASIM_LISA_H_

#include "lisasim-tens.h"
#include "lisasim-noise.h"

// ---- define file-wide numerical constants

// seconds per year; this is just 60x60x24x365 (it differs from the LISA Simulator)

const double secondsperyear = 31536000.0;
const double yearspersecond = 3.1709791983764586e-8;

// LISA's angular velocity

const double Omega = 2.0 * M_PI * yearspersecond;

// orbital radius of the guiding center (yrs)
// LISA simulator has Rgc = 1.49597870660e11 m = 499.005 s
 
static const double Rgc = 499.004;

// mean arm length of the LISA detector (yrs)
// LISA simulator has L = 5.0e9 m = 16.6782 s;

static const double Lstd = 16.6782;

// eccentricity of the LISA spacecraft orbits (for MontanaEccentric)
// LISA simulator has L/(2.0*sqrt(3.0)*Rgc)
// with L = 5.0e9, equal to 0.00964838

static const double ecc = 0.00964838;

// fit parameter for annual modulation of armlengths (1/sec, from armlength.nb)
// this might need fine tuning

const double delmodconst = 0.17647;

// generic LISA class

class LISA {
    protected:
        // guessL is needed by the generic version of armlength()
        // it needs to be initialized by the constructor of all the derived
        // LISA classes that do not define armlength
    
        double guessL[4];

    public:
        LISA() {};
        virtual ~LISA() {};

        virtual void reset() {};

	// will return pointer to itself; overridden by NoisyLISA, which returns the basic LISA

	virtual LISA *thislisa() {
	    return this;
	}

        // generic version of putn; uses (delayed) differences of putp, calling
        // armlength to get the right delay

        virtual void putn(Vector &n, int arm, double t);

        // putp should never be called for base LISA

        virtual void putp(Vector &p, int craft, double t) = 0;

        // generic version of armlength: will use putp iteratively to find
        // the delay corresponding to a photon trajectory backward from t
        // along "arm"

        virtual double armlength(int arm, double t);

        // these extra methods are needed to load arrays (not Vector objects)
        // with the spacecraft positions and links

        virtual void putn(double n[], int arm, double t) {
            Vector temp;
            
            putn(temp, arm, t);
            
            n[0] = temp[0];
            n[1] = temp[1];
            n[2] = temp[2];
	}
        
        virtual void putp(double p[], int craft, double t) {
            Vector temp;
            
            putp(temp, craft, t);
            
            p[0] = temp[0];
            p[1] = temp[1];
            p[2] = temp[2];
	}
};

class OriginalLISA : public LISA {
    protected:
        double L[4];
        
        Vector initn[4];
        Vector initp[4];
        
    public:
        // accept the armlength in seconds
    
        OriginalLISA(double arm1,double arm2,double arm3);

        // OriginalLISA defines optimized (look-up) versions of putn and armlength
    
        virtual void putn(Vector &n, int arm, double t);
        virtual void putp(Vector &p, int craft, double t);
    
        virtual double armlength(int arm, double t);
};

// for the moment, use the following class only for TDInoise

class ModifiedLISA : public OriginalLISA {
    private:
        double sagnac[4];
        double Lc[4], Lac[4];
        
    public:
        // accept the armlength in seconds
    
        ModifiedLISA(double arm1,double arm2,double arm3);

        // OriginalLISA defines a computed version of armlength
        // however, it uses the base putn
    
        void putp(Vector &p, int craft, double t);
    
        double armlength(int arm, double t);
};

// The following classes model the geometry of LISA and can be used also for "signal" TDI

class CircularRotating : public LISA {
    private:
        double scriptl;
        double R;
        double L;

        double eta0;
        double xi0;

        // Trick: we use 1-3 indexing for LISA positions and vectors, so we need to allocate 4 of everything

        double delmodamp;
        double delmodph[4];

        Vector initn[4];
        Vector initp[4];
        
        double rotationtime;
        Vector center;
        Tensor rotation;

        void settime(double t);

	void initialize(double e0, double x0, double sw);

    public:
        CircularRotating(double eta0 = 0.0,double xi0 = 0.0,double sw = 1.0);
	CircularRotating(double myL, double e0, double x0, double sw);

        // CircularRotating defines a computed (fitted) version of armlength
        // however, putn defaults to the base version
        // use genarmlength to get the base (exact) armlength,
        // use oldputn to get the nondelayed (rotated-only) n's

        void putp(Vector &p,int craft,double t);
        
        double armlength(int arm, double t);

        void oldputn(Vector &n,int arm,  double t);
        double genarmlength(int arm, double t);
};

// all the static constants below should have file scope and be used throughout

class MontanaEccentric : public LISA {
    private:
        // LISA armlength 

        double L;

        // initial azimuthal position of the guiding center
    
        double kappa;
        
        // initial orientation of the spacecraft
        
        double lambda;
        
        // caching the positions of spacecraft
        
        Vector cachep[4];
        double cachetime[4];
        
        void settime(int craft,double t);
        
    public:
        MontanaEccentric(double kappa0 = 0.0,double lambda0 = 0.0);

        // MontanaEccentric defaults to the base versions of putn and armlength

        void putp(Vector &p,int craft,double t);
};

class NoisyLISA : public LISA {
    private:
        LISA *cleanlisa;

        InterpolateNoise *uperror[4], *downerror[4];

    public:
        NoisyLISA(LISA *clean,double starm,double sdarm);
        ~NoisyLISA();
        
        void reset();

       	LISA *thislisa() {
	    return cleanlisa;
	}

        double armlength(int arm, double t);

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

extern LISA *stdlisa();

#endif /* _LISASIM_LISA_H_ */
