/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#include "lisasim-lisa.h"
#include "lisasim-signal.h"
#include "lisasim-except.h"

#include <iostream>


// --- generic LISA class --------------------------------------------------------------

/** Fills the Vector n with "arm" for reception at time t. The base
    LISA version of putn uses delayed differences of putp; calls
    armlength to get the right delay */

void LISA::putn(Vector &n,int arms,double t) {
	assertArm(arms);

    // get spacecraft indices for b -> arms -> a
    
    int crafta = getRecv(arms);
    int craftb = getSend(arms);

    // this is consistent with pa(t) = pb(t-L(t;b->a)) + L(t;b->a) n
    // for instance, p_1 = p_2 + n_3

    Vector pa, pb;

    putp(pa,crafta,t);
    putp(pb,craftb,t-armlength(arms,t));

    n.setdifference(pa,pb);
    n.setnormalized();    
}

void LISA::putv(Vector &v, int craft, double t) {
    Vector pp, pm;
    
    putp(pp,craft,t + 0.5);
    putp(pm,craft,t - 0.5);
    
    v.setdifference(pp,pm);
}

/** Generic version of armlength. Will use putp iteratively to
    find the delay corresponding to a photon trajectory backward
    from t along "arm" */

double LISA::armlength(int arms, double t) {
	assertArm(arms);

    // get spacecraft indices for b -> arms -> a

    int crafta = getRecv(arms);
    int craftb = getSend(arms);

    // this is consistent with pa(t) = pb(t-L(t;b->a)) + L(t;b->a) n
    // for instance, p_1 = p_2 + n_3

    Vector pa, pb, n;

    putp(pa,crafta,t);

    // implement a simple bisection search for the correct armlength
    // use a 10% initial bracket

    const double tol = 1e-14;

    double newguess = guessL[abs(arms)];
    double hi = 1.10 * newguess, lo = 0.90 * newguess;

    double norm, guess;

    do {
        guess = newguess;
    
        putp(pb,craftb,t-guess);

        n.setdifference(pa,pb);

        norm = n.dotproduct() - guess*guess;

        if(norm > 0) {
            lo = guess;   // distance is spacelike; increase delay
        } else if(norm < 0) {
            hi = guess;   // distance is timelike; reduce delay
        } else {
            return guess; // should almost never happen
        }
        
        newguess = 0.5 * (hi + lo);
    } while( fabs(newguess - guess) > tol );

    return newguess;
}

double LISA::dotarmlength(int arm, double t) {
    return (armlength(arm,t + 0.5) - armlength(arm,t - 0.5));
}

void LISA::newretardtime(double t) {
    it = t;
    rt = t;
	
    trb = 0.0;
    tra = 0.0;
}
    
double LISA::retardedtime() {
    return rt;
};

double LISA::retardation() {
    return trb + tra;
};

void LISA::retard(int ret) {
    if (ret != 0) {
		trb += armlengthbaseline(ret,rt);  
		tra += armlengthaccurate(ret,rt);

		rt = (it - trb) - tra;
    }
};

void LISA::retard(LISA *anotherlisa, int ret) {
    if (ret != 0) {
		trb += anotherlisa->armlengthbaseline(ret,rt);  
		tra += anotherlisa->armlengthaccurate(ret,rt);

		rt = (it - trb) - tra;
    }
}

void LISA::setguessL(double t) {
	Vector pa,pb,n;

	for(int arm = 1;arm < 4;arm++) {
		putp(pa,getSend(arm),t);
		putp(pb,getRecv(arm),t);

		n.setdifference(pa,pb);
		guessL[arm] = sqrt(n.dotproduct());
	}
}


// --- OriginalLISA LISA class ---------------------------------------------------------

OriginalLISA::OriginalLISA(double L1,double L2,double L3) {
    // take armlengths in seconds (do not convert to years)

    L[1] = L1; L[2] = L2; L[3] = L3;
        
    // construct the p vectors with the required lengths;
                
    double dL = sqrt(2.0 * (L[2]*L[2] + L[3]*L[3]) - L[1]*L[1]);
    double dth = -acos((L[2]*L[2] + L[3]*L[3] - L[1]*L[1])/(2.0 * L[2] * L[3]));
    double th3 = acos((3.0*L[2]*L[2] + L[3]*L[3] - L[1]*L[1])/(2.0 * L[2] * dL));
    
    Vector tp1, tp2, tp3;
    
    tp1[0] = 0.0;                   tp1[1] = 0.0;                 tp1[2] = 0.0;
    tp2[0] = L[3] * cos(dth + th3); tp2[1] = L[3] * sin(dth+th3); tp2[2] = 0.0;
    tp3[0] = L[2] * cos(th3);       tp3[1] = L[2] * sin(th3);     tp3[2] = 0.0;

    // align them so that p1 + p2 + p3 = 0
    // need to do only the x axis, the others are guaranteed
    
    double xoffset = (tp1[0] + tp2[0] + tp3[0]) / 3.0;

    // now assign them to initp; reflect axes and exchange 2<>3 so that for equal arms we get the
    // initp values used in CircularRotating
    
    initp[1][0] = -(tp1[0] - xoffset); initp[1][1] = -tp1[1]; initp[1][2] = -tp1[2];
    initp[2][0] = -(tp3[0] - xoffset); initp[2][1] = -tp3[1]; initp[2][2] = -tp3[2];
    initp[3][0] = -(tp2[0] - xoffset); initp[3][1] = -tp2[1]; initp[3][2] = -tp2[2];

    // now compute the corresponding n's as differences of p's, and normalize
    
    initn[1].setdifference(initp[2],initp[3]);
    initn[1].setnormalized();

    initn[2].setdifference(initp[3],initp[1]);
    initn[2].setnormalized();

    initn[3].setdifference(initp[1],initp[2]);                
    initn[3].setnormalized();
}

// OriginalLISA does not move; hence the length (but not the
// direction!) of positive and negative arms is the same

void OriginalLISA::putn(Vector &n,int arm,double t) {
	assertArm(arm);

	double sign = arm > 0 ? 1.0 : -1.0;

    n.setproduct(sign,initn[abs(arm)]);
}

void OriginalLISA::putp(Vector &p,int craft,double t) {
	assertCraft(craft);

    p = initp[craft];
}

double OriginalLISA::armlength(int arms, double t) {
	assertArm(arms);

	return L[abs(arms)];
}

// --- ModifiedLISA LISA class ---------------------------------------------------------

ModifiedLISA::ModifiedLISA(double arm1,double arm2,double arm3) : OriginalLISA(arm1,arm2,arm3) {
    // compute everything in seconds

    for(int i=1;i<4;i++) {
        // get spacecraft indices for b -> arms -> a

        int crafta = getRecv(i);
        int craftb = getSend(i);

        double La = sqrt(initp[crafta].dotproduct());
        double Lb = sqrt(initp[craftb].dotproduct());

        sagnac[i] = La * Lb * sqrt(1.0 - (La*La + Lb*Lb - L[i]*L[i])*(La*La + Lb*Lb - L[i]*L[i]) / (4.0 * La*La * Lb*Lb)) * Omega;

		sagnac[i] = Omega * (initp[craftb][0]*initp[crafta][1] - initp[crafta][0]*initp[craftb][1]);

        Lc[i]  = L[i] + sagnac[i];
        Lac[i] = L[i] - sagnac[i];

    	guessL[i] = L[i];
    }
}

// implement simple LISA rotation in xy plane with period of a year
// take input time in seconds

void ModifiedLISA::putp(Vector &p,int craft,double t) {
	assertCraft(craft);

	p[0] = cos(Omega*t) * initp[craft][0] - sin(Omega*t) * initp[craft][1];
	p[1] = sin(Omega*t) * initp[craft][0] + cos(Omega*t) * initp[craft][1];
	p[2] = initp[craft][2];
}

// positive arms are corotating (have longer arms), negative arms are counterrotating (shorter arms)

double ModifiedLISA::armlength(int arm, double t) {
	assertArm(arm);

    if(arm > 0) {
        return(Lc[arm]);
	} else {
        return(Lac[-arm]);
	}
}

double ModifiedLISA::genarmlength(int arm, double t) {
    return LISA::armlength(arm,t);
}

// ??? modernized up to here

// --- CircularRotating LISA class -----------------------------------------------------

// full constructor; initialize positions and arm vectors
// measure everything in seconds
// if the last argument is negative, switch 2 and 3; needed for coherence with the Montana simulator

CircularRotating::CircularRotating(double myL, double e0, double x0, double s, double t0)
    : L(myL), toffset(t0) {

	initialize(e0,x0,s);
}

CircularRotating::CircularRotating(double e0, double x0, double s, double t0)
    : L(Lstd), toffset(t0) {

	initialize(e0,x0,s);
}

void CircularRotating::initialize(double e0, double x0, double s) {
    // set distance from the Sun (s)

    R = Rgc;

    // set distance from guiding center (s)

    scriptl = L / sqrt(3.0);

    // define guessL, in case we wish to call the base-LISA armlength()
    // for comparison with our fitted expressions

    guessL[1] = L; guessL[2] = L; guessL[3] = L;

    // set initial phases

    eta0 = e0;
    xi0 = x0;
    sw = s;

    // construct the arms

    initn[1][0] =  0.0; initn[1][1] = -1.0; initn[1][2] = 0.0;
    initn[2][0] = -cos(M_PI/6.0); initn[2][1] = sin(M_PI/6.0); initn[2][2] = 0.0;
    initn[3][0] =  cos(M_PI/6.0); initn[3][1] = sin(M_PI/6.0); initn[3][2] = 0.0;
    
    initp[1][0] = scriptl; initp[1][1] = 0.0; initp[1][2] = 0.0;
    initp[2][0] = scriptl * cos(2.*M_PI/3.0); initp[2][1] = -scriptl * sin(2.*M_PI/3.0); initp[2][2] = 0.0;
    initp[3][0] = scriptl * cos(4.*M_PI/3.0); initp[3][1] = -scriptl * sin(4.*M_PI/3.0); initp[3][2] = 0.0;

    // switch the arms 2 and 3 if required; needed for comparison with Montana simulator

    if (sw < 0.0) {
        double tmp;
        
        for(int i=0;i<3;i++) {
            tmp = initp[2][i];
            initp[2][i] = initp[3][i];
            initp[3][i] = tmp;

            // the arm switch is untested
        
            initn[1][i] = -initn[1][i];
            tmp = initn[2][i];
            initn[2][i] = -initn[3][i];
            initn[3][i] = -tmp;
        }
    }
    
    // set amplitude and phase of delay modulations
    
    delmodamp = (sw > 0.0 ? 1.0 : -1.0) * R * L * Omega;
    
    delmodph[1] = xi0;
    delmodph[2] = sw > 0.0 ? xi0 + 4.*M_PI/3.0 : xi0 + 2.*M_PI/3.0;
    delmodph[3] = sw > 0.0 ? xi0 + 2.*M_PI/3.0 : xi0 + 4.*M_PI/3.0;

    settime(0);
}

// Set the correct components of the rotation matrix for the current time
// set also the position of the guiding center

void CircularRotating::settime(double t) {
    const double zeta = -M_PI/6.0;

    // {eta -> Omega time + eta0, xi -> -Omega time + xi0}
    // {elat(beta) -> zeta, elon(lambda) -> eta, psi -> xi}
    // time is measured in seconds
    
    rotation.seteuler(zeta,Omega*(t+toffset)+eta0,-Omega*(t+toffset)+xi0);
    
    // R{Cos[Omega t + eta0], Sin[Omega t + eta0], 0}; leave center[2] = 0.0
    
    center[0] = R * cos(Omega*(t+toffset)+eta0);
    center[1] = R * sin(Omega*(t+toffset)+eta0);
}

// Return the position of "craft" at time t in vector p
// does so by multiplying the initial position of the craft (with respect to the guiding center)
// by a rotation matrix computed at the current time (if it is not already cached)

void CircularRotating::putp(Vector &p,int craft,double t) {
    if (t != rotationtime) settime(t);

    p.setproduct(rotation, initp[craft]);
    
    // forget the z axis for the LISA center
    
    p[0] += center[0];
    p[1] += center[1];
}

// fit to armlength modulation from armlength.nb

double CircularRotating::armlength(int arm, double t) {
	assertArm(arm);

    if(arm > 0) {
		return L + delmodamp * sin(Omega*(t+toffset) - delmodph[arm]);
    } else {
		return L - delmodamp * sin(Omega*(t+toffset) - delmodph[-arm]);
    }
}

double CircularRotating::armlengthbaseline(int arm, double t) {
    return L;
}

double CircularRotating::armlengthaccurate(int arm, double t) {
	assertArm(arm);
	
    if(arm > 0) {
		return delmodamp * sin(Omega*(t+toffset) - delmodph[arm]);
    } else {
		return -delmodamp * sin(Omega*(t+toffset) - delmodph[-arm]);
    }
}

// call the generic version

double CircularRotating::genarmlength(int arm, double t) {
    return LISA::armlength(arm,t);
}

// --- EccentricInclined LISA class -----------------------------------------------------

EccentricInclined::EccentricInclined(double e0,double x0,double sw,double t0) {
    L = Lstd;
    toffset = t0;

    // initialize guessL to the initial guess for the length of arms

    guessL[1] = L; guessL[2] = L; guessL[3] = L;

    eta0 = e0;
    xi0 = x0;

    kappa = e0;
    lambda = x0 + e0 - 3.*M_PI/2.0;
    swi = sw;

    // initialize constants for approximation to armlength

    pdelmod = (swi > 0.0 ? 1.0 : -1.0)*(Omega * Rgc * L) - (15.0/32.0) * (ecc*L);
    mdelmod = (swi > 0.0 ? -1.0 : 1.0)*(Omega * Rgc * L) - (15.0/32.0) * (ecc*L);
    delmod3 = (1.0/32.0) * (ecc*L);

    delmodph[1] = xi0;
    delmodph[2] = (swi > 0.0) ? 4.*M_PI/3.0 + xi0 : 2.*M_PI/3.0 + xi0;
    delmodph[3] = (swi > 0.0) ? 2.*M_PI/3.0 + xi0 : 4.*M_PI/3.0 + xi0;
    
    delmodph2 = 3.0*xi0;

    // Initialize the position cache

    settime(1,0.0); settime(2,0.0); settime(3,0.0);
}

// positions of spacecraft according to the LISA simulator

void EccentricInclined::settime(int craft, double t) {
    const double sqecc = ecc*ecc;
    const double sqrt3 = sqrt(3.0);

    double alpha = Omega*(t + toffset) + kappa;

    double beta;

	assertCraft(craft);

    switch(craft) {
	case 1:
		beta = lambda;
		break;
	case 2:
		beta = (swi > 0.0) ? 4.0*M_PI/3.0 + lambda : 2.0*M_PI/3.0 + lambda;
		break;
	case 3:
		beta = (swi > 0.0) ? 2.0*M_PI/3.0 + lambda : 4.0*M_PI/3.0 + lambda;
		break;
	default:
		beta = 0.0; // should never get here, having the assertCraft above
		break;
    }
	    
    cachep[craft][0] =   0.5 * Rgc * ecc * ( cos(2.0*alpha-beta) - 3.0*cos(beta) )
                       + 0.125 * Rgc * sqecc * ( 3.0*cos(3.0*alpha-2.0*beta) - 5.0*( 2.0*cos(alpha)+cos(alpha-2.0*beta) ) )
                       + Rgc * cos(alpha);
                          
    cachep[craft][1] =   0.5 * Rgc * ecc * ( sin(2.0*alpha-beta) - 3.0*sin(beta) )
                       + 0.125 * Rgc * sqecc * ( 3.0*sin(3.0*alpha-2.0*beta) - 5.0*( 2.0*sin(alpha)-sin(alpha-2.0*beta) ) )
                       + Rgc * sin(alpha);
           
    cachep[craft][2] = - sqrt3 * Rgc * ecc * cos(alpha-beta) 
                       + sqrt3 * Rgc * sqecc * ( cos(alpha-beta)*cos(alpha-beta) + 2.0*sin(alpha-beta)*sin(alpha-beta) );                                               

    cachetime[craft] = t;
}

void EccentricInclined::putp(Vector &p, int craft, double t) {
    if (t != cachetime[craft]) settime(craft,t);

    for(int i=0;i<3;i++)
        p[i] = cachep[craft][i];
}

double EccentricInclined::armlength(int arm, double t) {
	assertArm(arm);

    if(arm > 0) {
		return L + pdelmod * sin(Omega*(t+toffset) - delmodph[arm]) 
			     + delmod3 * sin(Omega3*(t+toffset) - delmodph2);
    } else {
		return L + mdelmod * sin(Omega*(t+toffset) - delmodph[-arm])
			     + delmod3 * sin(Omega3*(t+toffset) - delmodph2);
    }
}

double EccentricInclined::armlengthbaseline(int arm, double t) {
    return L;
}

double EccentricInclined::armlengthaccurate(int arm, double t) {
	assertArm(arm);

    if(arm > 0) {
		return pdelmod * sin(Omega*(t+toffset) - delmodph[arm]) 
			 + delmod3 * sin(Omega3*(t+toffset) - delmodph2);
    } else {
		return mdelmod * sin(Omega*(t+toffset) - delmodph[-arm])
			 + delmod3 * sin(Omega3*(t+toffset) - delmodph2);
    }
}

double EccentricInclined::genarmlength(int arm, double t) {
    return LISA::armlength(arm,t);
}


// --- SampledLISA ---

SampledLISA::SampledLISA(double *sc1,long length1,double *sc2,long length2,double *sc3,long length3,
                         double deltat,double prebuffer,int interplen) {
    for(int i=0;i<3;i++) {
        buffer[1][i] = new double[length1/3];
        buffer[2][i] = new double[length2/3];
        buffer[3][i] = new double[length3/3];
            
        for(int k=0;k<length1/3;k++) buffer[1][i][k] = sc1[k*3 + i];
        for(int k=0;k<length2/3;k++) buffer[2][i][k] = sc2[k*3 + i];
        for(int k=0;k<length3/3;k++) buffer[3][i][k] = sc3[k*3 + i];
        
        sampledp[1][i] = new SampledSignal(buffer[1][i],length1/3,deltat,prebuffer,1.0,0,interplen);
        sampledp[2][i] = new SampledSignal(buffer[2][i],length2/3,deltat,prebuffer,1.0,0,interplen);
        sampledp[3][i] = new SampledSignal(buffer[3][i],length3/3,deltat,prebuffer,1.0,0,interplen);
    }
    
    for(int c=1;c<4;c++) {
        guessL[c] = Lstd;
    }
}

SampledLISA::~SampledLISA() {
    for(int i=0;i<3;i++) {
        delete sampledp[1][i];
        delete sampledp[2][i];
        delete sampledp[3][i];
        
        delete [] buffer[1][i];
        delete [] buffer[2][i];
        delete [] buffer[3][i];
    }
}

void SampledLISA::putp(Vector &p,int craft,double t) {
	assertCraft(craft);

	for(int i=0;i<3;i++) {
		p[i] = sampledp[craft][i]->value(t);
	}
}


// --- PyLISA ---

void PyLISA::reset() {
	return baseLISA->reset();
}

LISA* PyLISA::physlisa() {
	return baseLISA;
}

double PyLISA::armlength(int arm, double t) {
	PyObject *arglist, *result;
  
	double dres = 0.0;
	
	arglist = Py_BuildValue("(id)",arm,t);        // Build argument list
	result = PyEval_CallObject(armfunc,arglist);  // Call Python
	Py_DECREF(arglist);                           // Trash arglist
	if (result) dres = PyFloat_AsDouble(result);  // If no errors, return double
	Py_XDECREF(result);                           // Trash result
	return dres;
}

double PyLISA::armlengthbaseline(int arm, double t) {
	return armlength(arm,t);
}

double PyLISA::armlengthaccurate(int arm, double t) {
	return 0.0;
}

void PyLISA::putn(Vector &n, int arm, double t) {
	baseLISA->putn(n,arm,t);
}
    
void PyLISA::putp(Vector &p, int craft, double t) {
	baseLISA->putp(p,craft,t);
}

// --- AllPyLISA ---

// problem here: setLguesses needs the virtual putp, which may not be ready

AllPyLISA::AllPyLISA(PyObject *cfunc,PyObject *afunc) : craftfunc(cfunc), armlengthfunc(afunc) {
	setguessL();
}

void AllPyLISA::reset() {
	setguessL();
}

double AllPyLISA::armlength(int arm, double t) {
	if (armlengthfunc != 0) {
		PyObject *arglist, *result;
    
		double dres = 0.0;
    
		arglist = Py_BuildValue("(id)",arm,t);              // Build argument list
		result = PyEval_CallObject(armlengthfunc,arglist);  // Call Python
		Py_DECREF(arglist);                                 // Trash arglist
		if (result) dres = PyFloat_AsDouble(result);        // If no errors, return double
		// no type checking!
		Py_XDECREF(result);                                 // Trash result
		return dres;
	} else {
		return LISA::armlength(arm,t);
	}
}

double AllPyLISA::armlengthbaseline(int arm, double t) {
	return armlength(arm,t);
}

double AllPyLISA::armlengthaccurate(int arm, double t) {
	return 0.0;
}
    
void AllPyLISA::putp(Vector &p, int craft, double t) {
	PyObject *arglist, *result;
        
	arglist = Py_BuildValue("(id)",craft,t);        // Build argument list
	result = PyEval_CallObject(craftfunc,arglist);  // Call Python
	Py_DECREF(arglist);                             // Trash arglist
	if (result) {                                   // If no errors, get results
		p[0] = PyFloat_AsDouble(PyTuple_GetItem(result,0));
		p[1] = PyFloat_AsDouble(PyTuple_GetItem(result,1));
		p[2] = PyFloat_AsDouble(PyTuple_GetItem(result,2));
	}
	Py_XDECREF(result);                             // Trash result
}


// --- CacheLengthLISA (including LISASource) ---

// always use the computed armlengths, since they are guaranteed to exist
// for all LISAs

double LISASource::getvalue(long pos) {
	return basiclisa->LISA::armlength(arm,pos*deltat - prebuffer);
};

CacheLengthLISA::CacheLengthLISA(LISA *l,long length,double deltat,int interplen)
    : basicLISA(l) {
	try {
		interp = getInterpolator(interplen);
        dinterp = getDerivativeInterpolator(interplen);
	} catch (ExceptionUndefined &e) {
		std::cerr << "CacheLengthLISA::CacheLengthLISA(...): undefined interpolator length "
				  << interplen << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;
			
		throw e;		
	}

	double prebuffer = interplen * deltat;

	for(int i=1;i<4;i++) {
		lisafuncs[i] = new LISASource(length,deltat,prebuffer,basicLISA,i);
		lisafuncs[i+3] = new LISASource(length,deltat,prebuffer,basicLISA,-i);
	}

	for(int i=1;i<7;i++) {
		armlengths[i] = new InterpolatedSignal(lisafuncs[i],interp,deltat,prebuffer);
	}

	for(int i=1;i<7;i++) {
		dotarmlengths[i] = new InterpolatedSignal(lisafuncs[i],dinterp,deltat,prebuffer,1.0/deltat);
	}
		
	if(l->physlisa() == l) {
	   physLISA = this;
	} else {
	   physLISA = new CacheLengthLISA(l->physlisa(),length,deltat,interplen);
	}
}
    
CacheLengthLISA::~CacheLengthLISA() {
    if(physLISA != this) delete physLISA;

	for(int i=1;i<7;i++) delete dotarmlengths[i];
	for(int i=1;i<7;i++) delete armlengths[i];

	for(int i=1;i<7;i++) delete lisafuncs[i];

    delete dinterp;
	delete interp;
}

void CacheLengthLISA::reset() {
    if(physLISA != this) physLISA->reset();

	for(int i=1;i<7;i++) dotarmlengths[i]->reset();
    for(int i=1;i<7;i++) armlengths[i]->reset();
}

LISA* CacheLengthLISA::physlisa() {
    return physLISA;
}

double CacheLengthLISA::armlength(int arm, double t) {
	assertArm(arm);

	if(arm > 0) {
		return armlengths[arm]->value(t);
	} else {
		return armlengths[3-arm]->value(t);
	}
}

double CacheLengthLISA::armlengthbaseline(int arm, double t) {
	return armlength(arm,t);
}

double CacheLengthLISA::armlengthaccurate(int arm, double t) {
	return 0.0;
}

double CacheLengthLISA::dotarmlength(int arm, double t) {
	assertArm(arm);

	if(arm > 0) {
		return dotarmlengths[arm]->value(t);
	} else {
		return dotarmlengths[3-arm]->value(t);
	}
}

// need this because basicLISA's putn will call basicLISA's armlength

void CacheLengthLISA::putn(Vector &n,int arms,double t) {
	assertArm(arms);

    int crafta = getRecv(arms);
    int craftb = getSend(arms);

    Vector pa, pb;

    basicLISA->putp(pa,crafta,t);

	if(arms > 0) {
        basicLISA->putp(pb,craftb,t - armlengths[arms]->value(t));
	} else {
        basicLISA->putp(pb,craftb,t - armlengths[3-arms]->value(t));
	}

    n.setdifference(pa,pb);
    n.setnormalized();
}

void CacheLengthLISA::putp(Vector &p, int craft, double t) {
	basicLISA->putp(p,craft,t);
}






