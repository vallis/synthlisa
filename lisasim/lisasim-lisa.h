/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#ifndef _LISASIM_LISA_H_
#define _LISASIM_LISA_H_

#include "lisasim-tens.h"
#include "lisasim-signal.h"
#include "lisasim-except.h"

#include <iostream>

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
static const double eccstd = 0.00964838;

/// Base LISA geometry class.


// --- spacecraft/link indexing ---

// check index validity

inline void assertArm(int arm) {
	if(arm < -3 || arm > 3 || arm == 0) {
		std::cerr << "synthlisa::assertArm: invalid arm index "
				  << arm << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;
	
		ExceptionUndefined e;
		throw e;
	}
}

inline void assertCraft(int sc) {
	if(sc < 1 || sc > 3) {
		std::cerr << "synthlisa::assertCraft invalid spacecraft index "
				  << sc << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;
	
		ExceptionUndefined e;
		throw e;
	}
}

// convert link to spacecraft
// assume we're given good arm indexes

inline int getSend(int arm) {
    if(arm > 0) {
        return (arm ==  1) ? 3 : (arm - 1);
    } else {
        return (arm == -3) ? 1 : (1 - arm);
    }
}

inline int getRecv(int arm) {
    if(arm > 0) {
        return (arm ==  3) ? 1 : ( arm + 1);
    } else {
        return (arm == -1) ? 3 : (-arm - 1);
    }
}


class LISA {
 private:
    double it, rt;
    double trb, tra;

 protected:
    /** Initial armlength guess for the generic version of
	armlength(). It should be initialized by the constructor of
	all the derived LISA classes that do not define armlength. */    
    double guessL[4];
	
	void setguessL(double time = 0.0);

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
    virtual void putp(LISA *anotherlisa,Vector &p, int craft, double t);

    /* Simple implementation of spacecraft velocity by first-order
    one-second finite-difference expression */
    virtual void putv(Vector &v, int craft, double t);

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

    virtual double dotarmlength(int arm, double t);

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

class ApproxLISA {
    public:
      virtual ~ApproxLISA() {};

      virtual double geteta0() = 0;
      virtual double getxi0() = 0;
      virtual double getsw() = 0;
};

class CircularRotating : public LISA, public ApproxLISA {
 private:
    double scriptl;
    double R;
    double L;

    double toffset;

    double eta0;
    double xi0;
    double sw;

    Vector initn[4];
    Vector initp[4];

    Vector center;
    Tensor rotation;
    double rotationtime;

    // Trick: we use 1-3 indexing for LISA positions and vectors, so we need to allocate 4 of everything

    double delmodamp;
    double delmodph[4];
    
    void initialize(double e0, double x0, double sw);
    void settime(double t);
    
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
    
    double geteta0() {return eta0;};
    double getxi0() {return xi0;};
    double getsw() {return sw;};
};


/** Lagrange-point halo orbits. */

const double OmegaO = Omega;
const double OmegaR = 2.0*Omega;

class HaloAnalytic: public LISA {
 private:
    double scriptl;
    double R;
    double L;

    double toffset;

    Vector initn[4];
    Vector initp[4];

    Vector center;
    Tensor rotation;
    double rotationtime;

    double delmodamp, delconstamp;
    double delmodph[4];

    void settime(double t);

 public:   
    HaloAnalytic(double myL,double t0 = 0.0);
    
    void putp(Vector &p,int craft,double t);
    
    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);

    double genarmlength(int arm, double t);
};


/** Orbiting LISA with eccentric orbits. */

class EccentricInclined : public LISA, public ApproxLISA {
 private:
    // LISA armlength 

    double L, ecc;

    // initial azimuthal position of the guiding center
    // and initial orientation of the spacecraft
    // used internally, but computed from eta0 and xi0

    double toffset;

    double eta0;
    double xi0;

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

    void initialize(double e0, double x0, double sw);
    
 public:
    EccentricInclined(double eta0 = 0.0,double xi0 = 0.0,double sw = 1.0,double t0 = 0.0);
    EccentricInclined(double myL,double eta0,double xi0,double sw,double t0);

    void putp(Vector &p,int craft,double t);

    // EccentricInclined defines a computed (leading order) version of armlength
    // use genarmlength to get the exact armlength

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);

    double genarmlength(int arm, double t);
    
    double geteta0() {return eta0;};
    double getxi0() {return xi0;};
    double getsw() {return swi;};
};


// --- SampledLISA ---

class SampledLISA : public LISA {
 private:
    double *buffer[4][3];
    
    SampledSignal *sampledp[4][3];

 public:
    SampledLISA(double *sc1,long length1,double *sc2,long length2,double *sc3,long length3,double deltat,double prebuffer,int interplen = 1);    
    ~SampledLISA();

    void putp(Vector &p,int craft,double t);
};


// --- CacheLengthLISA ---

/* Instead of passing LISASource to InterpolatedSignal, It would be
   more elegant to have a generic class BufferedFunctionLISASource
   that could take any pointer to function, but there's a question
   about managing back-conversion from time to long pos.
*/

class LISASource : public BufferedSignalSource {
 private:
	double deltat, prebuffer;

	LISA *basiclisa;
	int arm;

 public:
	LISASource(long len,double dt,double pbt,LISA *lisa,int l)
		: BufferedSignalSource(len), deltat(dt), prebuffer(pbt), basiclisa(lisa), arm(l) {};

	double getvalue(long pos);
};

class CacheLengthLISA : public LISA {
 private:
	LISA *basicLISA, *physLISA;

	// use {1,2,3,-1,-2,-3} = {1,2,3,4,5,6} indexing for armlengths

	LISASource *lisafuncs[7];

	Interpolator *interp, *dinterp;

	InterpolatedSignal *armlengths[7], *dotarmlengths[7];

 public:
    CacheLengthLISA(LISA *lisa,long length,double deltat,int interplen = 4);
	~CacheLengthLISA();

	LISA *physlisa();

	void reset(); 

	double armlength(int arm, double t);

	double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);

    double dotarmlength(int arm, double t);

    void putn(Vector &n, int arm, double t);
    void putp(Vector &p, int craft, double t);
};


class ZeroLISA : public OriginalLISA {
 public:
    ZeroLISA() {};
    virtual ~ZeroLISA() {};

    double armlength(int arm, double t) {return 0.0;};
    double dotarmlength(int arm, double t) {return 0.0;};
};


// --- PyLISA ---

#include <Python.h>

class PyLISA : public LISA {
 private:
    PyObject *armfunc;

 public: 
    LISA *baseLISA;

    PyLISA(LISA *base,PyObject *func) : armfunc(func), baseLISA(base) {};

	void reset();

	LISA *physlisa();

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);

    void putn(Vector &n, int arm, double t);
    void putp(Vector &p, int craft, double t);
};


class AllPyLISA : public LISA {
 private:
    PyObject *craftfunc, *armlengthfunc;

 public: 
    AllPyLISA(PyObject *cfunc,PyObject *afunc = 0);

	void reset();

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);

    void putp(Vector &p, int craft, double t);
};

#endif /* _LISASIM_LISA_H_ */









