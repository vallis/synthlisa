/* File : lisasim-swig.i */
%module lisaswig
//%include typemaps.i
%{
#include "lisasim.h"
//#include "arrayobject.h"
//#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)
//#define PyArray_CONTIGUOUS(m) (ISCONTIGUOUS(m) ? Py_INCREF(m), m : \
//(PyArrayObject *)(PyArray_ContiguousFromObject((PyObject *)(m), (m)->descr->type_num, 0,0)))
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
        TDI(LISA *mylisa, Wave *mywave);
    
        double X(double t);
        double Y(double t);
        double Z(double t);
    
        double alpha(double t);
        double beta(double t);
        double gamma(double t);
    
        double zeta(double t);

	double P(double t);
	double E(double t);
	double U(double t);
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

        double U(double t);

        double Xm(double t);
};

%newobject stdnoise;
extern TDInoise *stdnoise(LISA *mylisa);

extern void printnoise(char *filename,TDInoise *mynoise,int samples,double samplingtime,char *observables);
extern void printsignal(char *filename,TDI *mysignal,int samples,double samplingtime,char *observables);

%include numpy.i

%apply double* IN_1D_DOUBLE { double *array };
extern void setnoise(double *array, TDInoise *mynoise,int samples,double samplingtime,char *observables);
extern void setsignal(double *array, TDI *mysignal,int samples,double samplingtime,char *observables);
