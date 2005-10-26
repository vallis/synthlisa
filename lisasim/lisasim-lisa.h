#ifndef _LISASIM_LISA_H_
#define _LISASIM_LISA_H_

#include "lisasim-tens.h"
#include "lisasim-noise.h"

#include <math.h>

/* Documentation rules: in the header, describe only objects that are
   actually defined, not just declared. A single line will become a
   Doxygen short description; otherwise, the first dot+space will
   separate the short from the long description. */

// ---- define file-wide numerical constants

/// Seconds per year. This is just 60x60x24x365 (it differs from the LISA Simulator).
const double secondsperyear = 31536000.0;
const double yearspersecond = 3.1709791983764586e-8;

/* Using 365.25 days, we'd have
   const double secondsperyear = 31557600.0;
   const double yearspersecond = 3.168808781402895e-8; */

/// Angular velocity of the LISA array orbit.
const double Omega = 2.0 * M_PI * yearspersecond;
const double Omega3 = 3.0 * Omega;

/** Orbital radius of the guiding center (yrs).
    The LISA simulator has Rgc = 1.49597870660e11 m = 499.005 s. */
static const double Rgc = 499.004;

/** Mean arm length of the LISA detector (yrs).
    The LISA simulator has L = 5.0e9 m = 16.6782 s. */
static const double Lstd = 16.6782;

/** Eccentricity of the LISA spacecraft orbits (used by
    EccentricInclined). LISA simulator has L/(2.0*sqrt(3.0)*Rgc);
    with L = 5.0e9, e = 0.00964838. */
static const double ecc = 0.00964838;

/// Base LISA geometry class.

class LISA {
 private:
    double it, rt;
    double trb, tra;

 protected:
    /** Initial armlength guess for the generic version of
	armlength(). It should be initialized by the constructor of
	all the derived LISA classes that do not define armlength. */    
    double guessL[4];
	
 public:
    LISA() {};

    virtual ~LISA() {};

    /// Resets LISA classes that have something to reset.
    virtual void reset() {};

    /** Returns a pointer to the TDI (nominal) LISA. Unless
	overridden, returns just "this". */
    virtual LISA *physlisa() { return this; }

    /*  Fills n with the photon direction vector along "arm" for
	reception at time t. */
    virtual void putn(Vector &n, int arm, double t);

    /** Fill p with the position of "craft" at time t. It is not
	defined for base LISA, which has no geometry, and thus becomes
	an abstract class. */
    virtual void putp(Vector &p, int craft, double t) = 0;

    /* Generic light propagation time along "arm" for reception at
       time t. */
    virtual double armlength(int arm, double t);

    /** Baseline value of the armlength (used to enhance precision in
	chain retardations). If we don't make a distinction between
	baseline and correction ("accurate"), return just the armlength. */
    virtual double armlengthbaseline(int arm, double t) {
	return armlength(arm,t);
    }

    /** Correction to the baseline armlength (used to enhance
	precision in chain retardations). If we don't make a
	distinction between baseline and correction ("accurate"),
	return 0. */
    virtual double armlengthaccurate(int arm, double t) {
	return 0.0;
    }

    /// FIX: I should remove all of these overloaded functions, not
    /// needed anymore now that Jeff's lisafast is not active anymore.

    /// Load an array rather than a vector with n (was used in LISAfast).
    virtual void putn(double n[], int arm, double t) {
	Vector temp;
            
	putn(temp, arm, t);
            
	n[0] = temp[0];
	n[1] = temp[1];
	n[2] = temp[2];
    }
        
    /// Load an array rather than a vector with p (was used in LISAfast).
    virtual void putp(double p[], int craft, double t) {
	Vector temp;
            
	putp(temp, craft, t);
            
	p[0] = temp[0];
	p[1] = temp[1];
	p[2] = temp[2];
    }

    virtual void newretardtime(double t);

    virtual double retardedtime();
    virtual double retardation();

    virtual void retard(int ret);
    virtual void retard(LISA *anotherlisa,int ret);
};


/// stationary LISA geometry

class OriginalLISA : public LISA {
 protected:
    /// Hold the three LISA armlengths in L[1]-L[3]
    
    double L[4];
        
    /// Hold the three LISA arms
    Vector initn[4];
    
    /// Hold the three LISA spacecraft positions
    Vector initp[4];
        
 public:
    // accept the armlength in seconds
    
    OriginalLISA(double arm1 = Lstd,double arm2 = Lstd,double arm3 = Lstd);

    // OriginalLISA defines optimized (look-up) versions of putn and armlength
    
    virtual void putn(Vector &n, int arm, double t);
    virtual void putp(Vector &p, int craft, double t);
	
    virtual double armlength(int arm, double t);
};

/// rotating LISA geometry

class ModifiedLISA : public OriginalLISA {
    private:
        double sagnac[4];
        double Lc[4], Lac[4];
        
    public:
        // accept the armlength in seconds
    
        ModifiedLISA(double arm1 = Lstd,double arm2 = Lstd,double arm3 = Lstd);

        // OriginalLISA defines a computed version of armlength
        // however, it uses the base putn
    
        void putp(Vector &p, int craft, double t);
    
        double armlength(int arm, double t);
        double genarmlength(int arm, double t);
};

/// Rigidly rotating, orbiting LISA.

class CircularRotating : public LISA {
    private:
        double scriptl;
        double R;
        double L;

	double toffset;

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
        CircularRotating(double eta0 = 0.0,double xi0 = 0.0,double sw = 1.0,double t0 = 0.0);
	CircularRotating(double myL,double e0,double x0,double sw,double t0);

        // CircularRotating defines a computed (fitted) version of armlength
        // use genarmlength to get the exact armlength,
        // use oldputn to get the nondelayed (rotated-only) n's

        // however, putn defaults to its base version

        void putp(Vector &p,int craft,double t);
        
        double armlength(int arm, double t);

	double armlengthbaseline(int arm, double t);
	double armlengthaccurate(int arm, double t);

        double genarmlength(int arm, double t);
};

/** Orbiting LISA with eccentric orbits. In the future it would be
    nice to have the eccentricity as a parameter */

class EccentricInclined : public LISA {
    private:
        // LISA armlength 

        double L;

        // initial azimuthal position of the guiding center
        // and initial orientation of the spacecraft
	// used internally, but computed from eta0 and xi0

	double toffset;

        double kappa; 
        double lambda;
	double swi;

	// constants needed in the approximation to the armlength

	double pdelmod, mdelmod, delmod3;

	double delmodph[4], delmodph2;
        
        // caching the positions of spacecraft
        
        Vector cachep[4];
        double cachetime[4];
        
        void settime(int craft,double t);
        
    public:
        EccentricInclined(double eta0 = 0.0,double xi0 = 0.0,double sw = 1.0,double t0 = 0.0);

        void putp(Vector &p,int craft,double t);

        // EccentricInclined defines a computed (leading order) version of armlength
        // use genarmlength to get the exact armlength

	double armlength(int arm, double t);

	double armlengthbaseline(int arm, double t);
	double armlengthaccurate(int arm, double t);

        double genarmlength(int arm, double t);
};

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

#include <Python.h>

class PyLISA : public LISA {
    public:
        LISA *baseLISA;

    // FIX: this ugliness removes a warning about armfunc being initialized before baseLISA
    //      but there should be a better way...

    private:
	PyObject *armfunc;

    public: 
        PyLISA(LISA *base,PyObject *func) : baseLISA(base), armfunc(func) {};
	virtual ~PyLISA() {};

        void reset() {
	    return baseLISA->reset();
	}

       	LISA *physlisa() {
	    return baseLISA;
	}

        double armlength(int arm, double t) {
	    PyObject *arglist, *result;

	    double dres = 0.0;
   
	    arglist = Py_BuildValue("(id)",arm,t);        // Build argument list
	    result = PyEval_CallObject(armfunc,arglist);  // Call Python
	    Py_DECREF(arglist);                           // Trash arglist
	    if (result) dres = PyFloat_AsDouble(result);  // If no errors, return double
	    Py_XDECREF(result);                           // Trash result
	    return dres;
	}

	double armlengthbaseline(int arm, double t) {
	    return armlength(arm,t);
	}

	double armlengthaccurate(int arm, double t) {
	    return 0.0;
	}
	
        void putn(Vector &n, int arm, double t) {
            baseLISA->putn(n,arm,t);
        }
        
        void putp(Vector &p, int craft, double t) {
            baseLISA->putp(p,craft,t);
        }

        void putn(double n[], int arm, double t) {
            baseLISA->putn(n,arm,t);
        }
        
        void putp(double p[], int craft, double t) {
            baseLISA->putp(p,craft,t);
	}
};

#endif /* _LISASIM_LISA_H_ */
