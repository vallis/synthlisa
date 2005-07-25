#include "lisasim-lisa.h"
#include "lisasim-noise.h"

using namespace std;
#include <iostream>

// --- generic LISA class --------------------------------------------------------------

/** Fills the Vector n with "arm" for reception at time t. The base
    LISA version of putn uses delayed differences of putp; calls
    armlength to get the right delay */

void LISA::putn(Vector &n,int arms,double t) {
    int arm = abs(arms);

    int crafta = (arm + 1 < 4) ? (arm + 1) : (arm - 2);
    int craftb = (arm + 2 < 4) ? (arm + 2) : (arm - 1);

    if(arms < 0) {
        int swap = crafta;
        crafta = craftb;
        craftb = swap;
    }

    // this is consistent with pa(t) = pb(t-L(t;b->a)) + L(t;b->a) n
    // for instance, p_1 = p_2 + n_3

    Vector pa, pb;

    putp(pa,crafta,t);
    putp(pb,craftb,t-armlength(arms,t));

    for(int i=0;i<3;i++)
        n[i] = pa[i] - pb[i];
    
    // normalize to a unit vector

    double norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);    
                
    for(int i=0;i<3;i++)
        n[i] /= norm;
}

/** Generic version of armlength. Will use putp iteratively to
    find the delay corresponding to a photon trajectory backward
    from t along "arm" */

double LISA::armlength(int arms, double t) {
    int arm = abs(arms);

    int crafta = (arm + 1 < 4) ? (arm + 1) : (arm - 2);
    int craftb = (arm + 2 < 4) ? (arm + 2) : (arm - 1);

    if(arms < 0) {
        int swap = crafta;
        crafta = craftb;
        craftb = swap;
    }

    // this is consistent with pa(t) = pb(t-L(t;b->a)) + L(t;b->a) n
    // for instance, p_1 = p_2 + n_3

    Vector pa, pb, n;

    putp(pa,crafta,t);

    // implement a simple bisection search for the correct armlength
    // use a 10% initial bracket

    double newguess = guessL[arm];
    double hi = 1.10 * newguess, lo = 0.90 * newguess;

    double guess;
    
    const double tol = 1e-14;

    do {
        guess = newguess;
    
        putp(pb,craftb,t-guess);

        for(int i=0;i<3;i++)
            n[i] = pa[i] - pb[i];

        // compute the invariant relativistic distance between the events
        // of emission and reception

        double norm = sqrt(-guess*guess + n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

        if(norm >= 0) {
            // distance is spacelike; increase delay

            lo = guess;
        } else {
            // distance is timelike; reduce delay
        
            hi = guess;
        }
        
        newguess = 0.5 * (hi + lo);
        
    } while( fabs(newguess - guess) > tol );

    return newguess;
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
    
    for(int i=0;i<3;i++) {
        initn[1][i] = initp[2][i] - initp[3][i];
        initn[2][i] = initp[3][i] - initp[1][i];
        initn[3][i] = initp[1][i] - initp[2][i];                
    }
        
    for(int j=1;j<4;j++) {
        double norm = sqrt(initn[j][0]*initn[j][0] + initn[j][1]*initn[j][1] + initn[j][2]*initn[j][2]);
        
        for(int i=0;i<3;i++)
            initn[j][i] /= norm;
    }
}

// OriginalLISA does not move; hence the length (but not the
// direction!) of positive and negative arms is the same

void OriginalLISA::putn(Vector &n,int arm,double t) {
    double sign = arm > 0 ? 1.0 : -1.0;

    for(int i=0;i<3;i++)
        n[i] = sign * initn[abs(arm)][i];
}

void OriginalLISA::putp(Vector &p,int craft,double t) {
    for(int i=0;i<3;i++)
        p[i] = initp[craft][i];
}

double OriginalLISA::armlength(int arm, double t) {
    return( L[abs(arm)] );
}

// --- ModifiedLISA LISA class ---------------------------------------------------------

ModifiedLISA::ModifiedLISA(double arm1,double arm2,double arm3) : OriginalLISA(arm1,arm2,arm3) {
    // compute everything in seconds

    for(int i=1;i<4;i++) {
        int crafta = (i + 1 < 4) ? i + 1 : i - 2;
        int craftb = (i + 2 < 4) ? i + 2 : i - 1;

        double La = sqrt(initp[crafta].dotproduct(initp[crafta]));
        double Lb = sqrt(initp[craftb].dotproduct(initp[craftb]));

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
    p[0] = cos(Omega*t) * initp[craft][0] - sin(Omega*t) * initp[craft][1];
    p[1] = sin(Omega*t) * initp[craft][0] + cos(Omega*t) * initp[craft][1];
    p[2] = initp[craft][2];
}

// positive arms are corotating (have longer arms), negative arms are counterrotating (shorter arms)

double ModifiedLISA::armlength(int arm, double t) {
    if(arm > 0)
        return(Lc[arm]);
    else
        return(Lac[-arm]);
}

double ModifiedLISA::genarmlength(int arm, double t) {
    return LISA::armlength(arm,t);
}

// --- CircularRotating LISA class -----------------------------------------------------

// full constructor; initialize positions and arm vectors
// measure everything in seconds
// if the last argument is negative, switch 2 and 3; needed for coherence with the Montana simulator

CircularRotating::CircularRotating(double myL, double e0, double x0, double sw, double t0)
    : L(myL), toffset(t0) {
	initialize(e0,x0,sw);
    }

CircularRotating::CircularRotating(double e0, double x0, double sw, double t0)
    : L(Lstd), toffset(t0) {
	initialize(e0,x0,sw);
    }

void CircularRotating::initialize(double e0, double x0, double sw) {
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
    // {psi -> xi, asc -> eta, dec -> zeta}
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

EccentricInclined::EccentricInclined(double eta0,double xi0,double sw,double t0) {
    L = Lstd;
    toffset = t0;

    // initialize guessL to the initial guess for the length of arms

    guessL[1] = L; guessL[2] = L; guessL[3] = L;

    kappa = eta0;
    lambda = xi0 + eta0 - 3.*M_PI/2.0;
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
	cout << "EccentricInclined::settime: invalid spacecraft index" << endl;
	abort(); 
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


// --- NoisyLISA class -----------------------------------------------------

// the error on the arms is defined in seconds

NoisyLISA::NoisyLISA(LISA *clean,double starm,double sdarm) {
    cleanlisa = clean;
    
    // to estimate size of noisebuffer, take the largest armlength at time zero, and add 10%

    double arm1 = cleanlisa->armlength(1,0.0);
    double arm2 = cleanlisa->armlength(2,0.0);
    double arm3 = cleanlisa->armlength(3,0.0);

    double maxarm = arm1 > arm2 ? arm1 : arm2;
    maxarm = maxarm > arm3 ? maxarm : arm3;

    double lighttime = 1.10 * maxarm;

    // create InterpolateNoise objects for the arm determination noises
    // we need triple retardations (septuple for 2nd-generation TDI)

    double pbtarm = 8.0 * lighttime + 2.0 * starm;

    // we use plain uncorrelated white noise, for the moment

    for(int link = 1; link <= 3; link++) {
        uperror[link] = new InterpolateNoise(starm, pbtarm, sdarm, 0.0);
        downerror[link] = new InterpolateNoise(starm, pbtarm, sdarm, 0.0);
    }
}

NoisyLISA::~NoisyLISA() {
    for(int link = 1; link <= 3; link++) {
        delete uperror[link];
        delete downerror[link];
    }
}

void NoisyLISA::reset() {
    for(int link = 1; link <= 3; link++) {
        uperror[link]->reset();
        downerror[link]->reset();
    }    
}

double NoisyLISA::armlength(int arm, double t) {
    if(arm > 0)
        return cleanlisa->armlength(arm,t) + (*uperror[arm])[t];
    else
        return cleanlisa->armlength(arm,t) + (*downerror[-arm])[t];
}

double NoisyLISA::armlengthbaseline(int arm, double t) {
    return cleanlisa->armlengthbaseline(arm,t);
}

double NoisyLISA::armlengthaccurate(int arm, double t)  {
    if(arm > 0)
        return cleanlisa->armlengthaccurate(arm,t) + (*uperror[arm])[t];
    else
        return cleanlisa->armlengthaccurate(arm,t) + (*downerror[-arm])[t];
}

// --- NominalLISA class ---------------------------------------------------

NominalLISA::NominalLISA(double eta0,double xi0,double sw,double t0) {
    reallisa = new EccentricInclined(eta0,xi0,sw,t0);

    L = Lstd;
    swi = sw;

    // these are the correct values of the user-set parameters
    // they won't be changed by setparameters, which works directly
    // on dL, toffset, pdelmod, mdelmod, delmod3

    Lnom = L;
    cmod = Omega * Rgc * L;
    emod = ecc * L;
    toff = t0;

    setparameters(Lnom,cmod,emod,toff);

    // we make these exact for the moment

    delmodph[1] = xi0;
    delmodph[2] = (swi > 0.0) ? 4.*M_PI/3.0 + xi0 : 2.*M_PI/3.0 + xi0;
    delmodph[3] = (swi > 0.0) ? 2.*M_PI/3.0 + xi0 : 4.*M_PI/3.0 + xi0;

    delmodph2 = 3.0*xi0;
}

void NominalLISA::setparameters(double l,double cm,double em,double to) {
    dL = l - L;

    toffset = to;
    
    pdelmod = (swi > 0.0 ? 1.0 : -1.0) * cm - (15.0/32.0) * em;
    mdelmod = (swi > 0.0 ? -1.0 : 1.0) * cm - (15.0/32.0) * em;

    delmod3 = (1.0/32.0) * em;
}

void NominalLISA::setparameters(double cm,double em,double to) {
    // setting cmod, emod, toff

    setparameters(Lnom,cm,em,to);
}

void NominalLISA::setparameters3(double p1,double p2,double p3) {
    // setting l, em, to
    setparameters(p1,cmod,p2,p3);

    // setting l, cm, to
    // setparameters(p1,p2,emod,p3);

    // setting l, cm, em
    // setparameters(p1,p2,p3,toff);
}

NominalLISA::~NominalLISA() {
    delete reallisa;
}

double NominalLISA::armlength(int arm, double t) {
    if(arm > 0) {
	return L + dL + pdelmod * sin(Omega*(t+toffset) - delmodph[arm]) 
	    + delmod3 * sin(Omega3*(t+toffset) - delmodph2);
    } else {
	return L + dL + mdelmod * sin(Omega*(t+toffset) - delmodph[-arm])
	    + delmod3 * sin(Omega3*(t+toffset) - delmodph2);
    }
}

double NominalLISA::armlengthbaseline(int arm, double t) {
    return L;
}

double NominalLISA::armlengthaccurate(int arm, double t) {
    if(arm > 0) {
	return dL + pdelmod * sin(Omega*(t+toffset) - delmodph[arm]) 
	    + delmod3 * sin(Omega3*(t+toffset) - delmodph2);
    } else {
	return dL + mdelmod * sin(Omega*(t+toffset) - delmodph[-arm])
	    + delmod3 * sin(Omega3*(t+toffset) - delmodph2);
    }
}

// --- LinearLISA class ---------------------------------------------------

LinearLISA::LinearLISA(double eta0,double xi0,double sw,double t0) {
    reallisa = new EccentricInclined(eta0,xi0,sw,t0);

    L = Lstd;

    for(int i=0;i<6;i++) {
	dL[i] = 0.0;
	dLdt[i] = 0.0;
    }
}

LinearLISA::~LinearLISA() {
    delete reallisa;
}

void LinearLISA::settimeoffset(double toff) {
    toffset = toff;
}

void LinearLISA::setparameters(double dl[6],double dldt[6]) {
    for(int i=0;i<6;i++) {
	dL[i] = dl[i];
	dLdt[i] = dldt[i];
    }
}

double LinearLISA::armlengtherror(int arm, double t) {
    return armlength(arm,t) - reallisa->armlength(arm,t);
}

// arm = 1 -> 0, 2 -> 1, 3 -> 2
// arm =-1 -> 3,-2 -> 4,-3 -> 5

double LinearLISA::armlength(int arm, double t) {
    if(arm > 0) {
	return L + dL[arm-1] + dLdt[arm-1]*(t-toffset);
    } else {
	return L + dL[2-arm] + dLdt[2-arm]*(t-toffset);
    }
}

double LinearLISA::armlengthbaseline(int arm, double t) {
    return L;
}

double LinearLISA::armlengthaccurate(int arm, double t) {
    if(arm > 0) {
	return dL[arm-1] + dLdt[arm-1]*(t-toffset);
    } else {
	return dL[2-arm] + dLdt[2-arm]*(t-toffset);
    }
}

// --- MeasureLISA class -----------------------------------------------------

// we're doing something funny here: we're defining Noise-like objects
// which will actually provide the LISA armlength measures, based on a regular
// (but noisy) sampling of the actual armlengths

MeasureLISA::MeasureLISA(LISA *clean,double starm,double sdarm,int swin) {
    cleanlisa = clean;
    
    // to estimate size of noisebuffer, take the largest armlength at time zero, and add 10%

    double arm1 = cleanlisa->armlength(1,0.0);
    double arm2 = cleanlisa->armlength(2,0.0);
    double arm3 = cleanlisa->armlength(3,0.0);

    double maxarm = arm1 > arm2 ? arm1 : arm2;
    maxarm = maxarm > arm3 ? maxarm : arm3;

    double lighttime = 1.10 * maxarm;

    // create InterpolateNoise objects for the arm determination noises
    // we need triple retardations (septuple for 2nd-generation TDI)

    if (swin == -1) swin = 2;

    double pbtarm = 8.0 * lighttime + (2.0 * swin + 1.0) * starm;

    // we use plain uncorrelated white noise, for the moment

    for(int link = 1; link <= 3; link++) {
        uperror[link] = new LISANoise(cleanlisa, link, starm, pbtarm, sdarm, swin);
        downerror[link] = new LISANoise(cleanlisa, -link, starm, pbtarm, sdarm, swin);
    }
}

MeasureLISA::~MeasureLISA() {
    for(int link = 1; link <= 3; link++) {
        delete uperror[link];
        delete downerror[link];
    }
}

void MeasureLISA::reset() {
    for(int link = 1; link <= 3; link++) {
        uperror[link]->reset();
        downerror[link]->reset();
    }    
}

// Unfortunately it is hard to make the distinction between the
// baseline and the accurate armlength within the LISANoise class. But
// this should be good enough.

double MeasureLISA::armlength(int arm, double t) {
    if(arm > 0)
        return (*uperror[arm])[t];
    else
        return (*downerror[-arm])[t];
}

double MeasureLISA::armlengthbaseline(int arm, double t) {
    if(arm > 0)
        return (*uperror[arm])[t];
    else
        return (*downerror[-arm])[t];
}

double MeasureLISA::armlengthaccurate(int arm, double t)  {
    return 0;
}
