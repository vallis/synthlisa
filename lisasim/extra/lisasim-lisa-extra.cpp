#include "lisasim-lisa.h"
#include "lisasim-noise.h"
#include "lisasim-lisa-extra.h"

using namespace std;
#include <iostream>

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
