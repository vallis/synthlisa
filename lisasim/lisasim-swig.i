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
    
    virtual double armlength(int arm, double t);
};

class ModifiedLISA : public OriginalLISA {
public:

    // accept the armlength in seconds
    
    ModifiedLISA(double arm1,double arm2,double arm3);
        
    double armlength(int arm, double t);
};

class CircularRotating : public LISA {
public:
	
    // three arguments: eta0, xi0, 2<->3 switch (1.0 or -1.0) 
    
    CircularRotating(double eta0,double xi0,double sw);
    
    double armlength(int arm, double t);
    
    double genarmlength(int arms, double t);
};

class NoisyLISA : public LISA {
public:
    NoisyLISA(LISA *clean,double starm,double sdarm);
    ~NoisyLISA();
	
    double armlength(int arm, double t);
};

class InterpolateNoise {
public:
    InterpolateNoise(double st, double pbt, double sd, double ex);
    ~InterpolateNoise();

    void reset();
        
    double inoise(double time);
};

class ExpGaussNoise {
public:
    
    ExpGaussNoise(double samplinginterval, double lapseinterval, double foldingtime, double spectraldensity);
    ~ExpGaussNoise();
    
    void reset();

    double enoise(double time);
};

class Wave;

class SimpleBinary : public Wave {
public:
    SimpleBinary(double freq, double initphi, double inc, double amp, double d, double a, double p);
};

class InterpolateMemory : public Wave {
public:
    InterpolateMemory(double *hpa, double *hca, long samples, double samplingtime, double lookback, double d, double a, double p);
};

class TDI {
public:
    TDI() {};
    virtual ~TDI() {};

    virtual void reset() {};

    virtual double X(double t);
    virtual double Y(double t);
    virtual double Z(double t);
    
    virtual double alpha(double t);
    virtual double beta(double t);
    virtual double gamma(double t);
    
    virtual double zeta(double t);

    virtual double P(double t);
    virtual double E(double t);
    virtual double U(double t);
    
    virtual double Xm(double t);

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);
};

class TDInoise : public TDI {
public:
    TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser);
    TDInoise(LISA *mylisa, LISA *physlisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser, double claser);
    ~TDInoise();
    
    void reset();

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
    double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);
};

class TDIsignal : public TDI {
public:
    TDIsignal(LISA *mylisa, Wave *mywave);
    TDIsignal(LISA *mylisa, LISA *physlisa, Wave *mywave);
            
    double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
};

%newobject stdnoise;
extern TDInoise *stdnoise(LISA *mylisa);

%newobject stdlisa;
extern LISA *stdlisa();

extern void printtdi(char *filename,TDI *mytdi,int samples,double samplingtime,char *observables);

%include numpy.i

%apply double* IN_1D_DOUBLE { double *array };
extern void settdi(double *array,TDI *mytdi,int samples,double samplingtime,char *observables);

extern void printn(CircularRotating *mylisa,int arm,double t);
