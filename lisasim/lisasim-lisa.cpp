#include "lisasim.h"

// --- generic LISA class --------------------------------------------------------------

// This is the generic version of "armlength", which takes differences between positions

double LISA::armlength(int arm, double t) {
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
        
    return(sqrt(diff.dotproduct(diff)));
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
