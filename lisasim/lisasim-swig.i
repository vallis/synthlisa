/* File : lisasim-swig.i */
%module lisaswig
%{
#include "lisasim.h"
%}

class LISA;

class OriginalLISA : public LISA {
    public:

        // accept the armlength in seconds

	OriginalLISA(double arm1,double arm2,double arm3);
};

class ModifiedLISA : public OriginalLISA {
    public:

        // accept the armlength in seconds
    
	ModifiedLISA(double arm1,double arm2,double arm3);
};

class CircularRotating : public LISA {
    public:
	
	// three arguments: eta0, xi0, 2<->3 switch (1.0 or -1.0) 
    
        CircularRotating(double eta0,double xi0,double sw);
};

class NoisyLISA : public LISA {
    public:
        NoisyLISA(LISA *clean,double starm,double sdarm);
        ~NoisyLISA();
};

class TDInoise {
    public:

        // claser is a correlation e-folding time

        TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser);
        
        // the two lisa pointers indicate which LISA is used to determine the TDI times
        // and which to determine the physical laser delays
        
        TDInoise(LISA *mylisa, LISA *physlisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser);
        
        ~TDInoise();
        
        void reset();

        // leave these here so we can show the cancellation of laser noise

        double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
        double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);

        double X(double t);
        double Y(double t);
        double Z(double t);

        double alpha(double t);
        double beta(double t);
        double gamma(double t);

        double zeta(double t);

        double P(double t);
        
        double E(double t);

        double Xm(double t);
};

extern void printnoise(char *filename,TDInoise *mynoise,int samples,double samplingtime,char *observables);

// The "%new" syntax here does not seem to be working. thisown is set to 0.

%new extern TDInoise *stdnoise(LISA *mylisa);

