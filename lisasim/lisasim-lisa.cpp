#include "lisasim-lisa.h"
#include "lisasim-noise.h"

// --- generic LISA class --------------------------------------------------------------

// This is the generic version of "armlength", which takes differences between positions
// generalized to take negative arguments, but it will still give the same answer

double LISA::armlength(int arm, double t) {
    Vector arm1, arm2;
    
    switch(abs(arm)) {
        case 1:
            putp(arm1,2,t); putp(arm2,3,t);
            break;
        case 2:
            putp(arm1,3,t); putp(arm2,1,t);
            break;
        case 3:
            putp(arm1,1,t); putp(arm2,2,t);
            break;            
    }
    
    Vector diff;
    for(int i=0;i<3;i++)
        diff[i] = arm1[i] - arm2[i];
        
    return(sqrt(diff.dotproduct(diff)));
}

// --- OriginalLISA LISA class ---------------------------------------------------------

OriginalLISA::OriginalLISA(double L1,double L2,double L3) {
    // convert from seconds to years

    L[1] = 3.17098E-8 * L1;
    L[2] = 3.17098E-8 * L2;
    L[3] = 3.17098E-8 * L3;
    
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

void OriginalLISA::putn(Vector &n,int arm,double t) {
    for(int i=0;i<3;i++)
        n[i] = initn[arm][i];
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
    // computing everything in years...

    double Omega = 2.0 * M_PI;

    for(int i=1;i<4;i++) {
        int crafta = (i + 1 < 4) ? i + 1 : i - 2;
        int craftb = (i + 2 < 4) ? i + 2 : i - 1;

        double La = sqrt(initp[crafta].dotproduct(initp[crafta]));
        double Lb = sqrt(initp[craftb].dotproduct(initp[craftb]));

        // see Mathematica file armlength.nb for this expression

        sagnac[i] = La * Lb * sqrt(1.0 - (La*La + Lb*Lb - L[i]*L[i])*(La*La + Lb*Lb - L[i]*L[i]) / (4.0 * La*La * Lb*Lb)) * Omega;

        Lc[i]  = L[i] + sagnac[i];
        Lac[i] = L[i] - sagnac[i];

//      cout.precision(8);
//      cout << "Light times on arm " << i << ": " << Lc[i] / 3.17098E-8 << " " << Lac[i] / 3.17098E-8 << endl;
    }
}

// this does not have to be so general, because the arm vectors are static up to a rotation.
// Still, ModifiedLISA should not need the n's often, and this putn can probably be adapted
// for a generic LISA

void ModifiedLISA::putn(Vector &n,int arms,double t) {
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

// implement simple LISA rotation in xy plane with period of a year

void ModifiedLISA::putp(Vector &p,int craft,double t) {
    const double twopi = 2.0*M_PI;
    
    p[0] = cos(twopi*t) * initp[craft][0] - sin(twopi*t) * initp[craft][1];
    p[1] = sin(twopi*t) * initp[craft][0] + cos(twopi*t) * initp[craft][1];
    p[2] = initp[craft][2];
}

double ModifiedLISA::armlength(int arm, double t) {
    if(arm > 0)
        return(Lc[arm]);
    else
        return(Lac[-arm]);
}

// --- CircularRotating LISA class -----------------------------------------------------

// full constructor; initialize positions and arm vectors; everything is measured in years
// if the last argument is negative, switch 2 and 3; needed for coherence with the Montana simulator

CircularRotating::CircularRotating(double e0, double x0, double sw) {
    scriptl = 3.0533885108232977E-7;
    R = 0.0000158233;
    L = 5.288624035993024E-7;

    eta0 = e0;
    xi0 = x0;

    initn[1][0] =  0.0; initn[1][1] = -1.0; initn[1][2] = 0.0;
    initn[2][0] = -cos(M_PI/6.0); initn[2][1] = sin(M_PI/6.0); initn[2][2] = 0.0;
    initn[3][0] =  cos(M_PI/6.0); initn[3][1] = sin(M_PI/6.0); initn[3][2] = 0.0;
    
    initp[1][0] = scriptl; initp[1][1] = 0.0; initp[1][2] = 0.0;
    initp[2][0] = scriptl * cos(2.*M_PI/3.0); initp[2][1] = -scriptl * sin(2.*M_PI/3.0); initp[2][2] = 0.0;
    initp[3][0] = scriptl * cos(4.*M_PI/3.0); initp[3][1] = -scriptl * sin(4.*M_PI/3.0); initp[3][2] = 0.0;

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
    
    settime(0);
}

// Set the correct components of the rotation matrix for the current time
// set also the position of the guiding center

void CircularRotating::settime(double t) {
    const double zeta = -M_PI/6.0;
    const double twopi = 2.0*M_PI;

    // {eta -> Omega time + eta0, xi -> -Omega time + xi0}
    // {psi -> xi, asc -> eta, dec -> zeta}
    // assume time measured in years 
    
    rotation.seteuler(zeta,twopi*t+eta0,-twopi*t+xi0);
    
    // R{Cos[Omega t + eta0], Sin[Omega t + eta0], 0}; leave center[2] = 0.0
    
    center[0] = R * cos(twopi*t+eta0);
    center[1] = R * sin(twopi*t+eta0);
}

// Return the unit vector along "arm" at time t in vector n
// does so by multiplying the initial vectors
// by a rotation matrix computed at the current time (if it is not already cached)

void CircularRotating::putn(Vector &n,int arm,double t) {
    if (t != rotationtime) settime(t);

    n.setproduct(rotation, initn[arm]);
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

// The length of arms is fixed to L for this model
// It will work also with negative arm arguments

double CircularRotating::armlength(int arm, double t) {
    return L;
}

// --- MontanaEccentric LISA class -----------------------------------------------------

// full constructor
// define the LISA Simulator kappa and lambda constants

MontanaEccentric::MontanaEccentric(double k0, double l0) {
    kappa = k0;
    lambda = l0;

    // Initialize the position cache

    settime(1,0.0); settime(2,0.0); settime(3,0.0);
}

// arm vectors as differences of position vectors; this sign convention
// is consistent with the one used for CircularRotating

void MontanaEccentric::putn(Vector &n, int arm, double t) {
    Vector arm1, arm2;
    
    switch(arm) {
        case 1:
            putp(arm1,2,t); putp(arm2,3,t);
            break;
        case 2:
            putp(arm1,3,t); putp(arm2,1,t);
            break;
        case 3:
            putp(arm1,1,t); putp(arm2,2,t);
            break;            
    }
    
    Vector diff;
    for(int i=0;i<3;i++)
        diff[i] = arm1[i] - arm2[i];

    // normalize

    double norm = 1.0/sqrt(diff.dotproduct(diff));
    
    for(int i=0;i<3;i++)
        n[i] = norm * diff[i];
}

// positions of spacecraft according to the LISA simulator

void MontanaEccentric::settime(int craft, double t) {
    const double twopi = 2.0*M_PI;
    const double sqecc = ecc*ecc;
    const double sqrt3 = sqrt(3.0);

    // with respect to the Simulator, our time is already in years
    // but the indexing of spacecraft is the same
    
    double alpha = twopi*t + kappa;
    double beta = twopi*(craft-1)/3.0 + lambda;

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

void MontanaEccentric::putp(Vector &p, int craft, double t) {
    if (t != cachetime[craft]) settime(craft,t);

    for(int i=0;i<3;i++)
        p[i] = cachep[craft][i];
}

// --- NoisyLISA class -----------------------------------------------------

// for simplicity, the error on the arms is defined in seconds

NoisyLISA::NoisyLISA(LISA *clean,double starm,double sdarm) {
    cleanlisa = clean;
    
    // to estimate size of noisebuffer, take an armlength at time zero,
    // and add 10%; need to convert from years to seconds (we use the factor from lisasim-tdi.cpp)
    // ah, this might not work for OriginalLISA and strange geometries

    double lighttime = 1.10 * (cleanlisa->armlength(1,0.0) / 3.17098E-8);

    // create InterpolateNoise objects for the arm determination noises
    // we need triple retardations

    double pbtarm = 3.0 * lighttime;

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

// convert the error from seconds to years

double NoisyLISA::armlength(int arm, double t) {
    if(arm > 0)
        return(cleanlisa->armlength(arm,t) + (*uperror[arm])[t/3.17098E-8] * 3.17098E-8);
    else
        return(cleanlisa->armlength(arm,t) + (*downerror[-arm])[t/3.17098E-8] * 3.17098E-8);
}
