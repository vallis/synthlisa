/* File : lisasim-swig.i */

%module lisaswig
%{
#include "lisasim.h"
%}

%include lisasim-typemaps.i

/* -------- LISA objects -------- */

class LISA;

class OriginalLISA : public LISA {
public:
    
    // accept the armlength in seconds

    OriginalLISA(double arm1,double arm2,double arm3);

    virtual double armlength(int arm, double t);

    virtual double armlengthbaseline(int arm, double t);
    virtual double armlengthaccurate(int arm, double t);

};

class ModifiedLISA : public OriginalLISA {
public:

    // accept the armlength in seconds
    
    ModifiedLISA(double arm1,double arm2,double arm3);
        
    double armlength(int arm, double t);

    double genarmlength(int arms, double t);
};

class CircularRotating : public LISA {
public:
	
    // three arguments: eta0, xi0, 2<->3 switch (1.0 or -1.0) 
    
    CircularRotating(double eta0,double xi0,double sw);
    CircularRotating(double myL,double eta0,double xi0,double sw);
    
    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
    
    double genarmlength(int arms, double t);
};

class EccentricInclined : public LISA {
public:
    EccentricInclined(double kappa0,double lambda0);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
    
    double genarmlength(int arms, double t);
};

class NoisyLISA : public LISA {
public:
    NoisyLISA(LISA *clean,double starm,double sdarm);
    ~NoisyLISA();
	
    double armlength(int arm, double t);
};

/* -------- Noise objects -------- */

class Noise;

%apply double *NUMPY_ARRAY_DOUBLE { double *noisebuf };

class InterpolateNoise : public Noise {
public:
    InterpolateNoise(double sampletime,double prebuffer,double density,double exponent);
    InterpolateNoise(double *noisebuf,long samples,double sampletime,double prebuffer,double norm);

    ~InterpolateNoise();

    void reset();
        
    virtual double inoise(double time);
};

class InterpolateNoiseBetter : public InterpolateNoise {
 public:

    InterpolateNoiseBetter(double sampletime,double prebuffer,double density,double exponent,int swindow);
    InterpolateNoiseBetter(double *noisebuf,long samples,double sampletime,double prebuffer,double norm,int swindow);

    double inoise(double time);
};

/* -------- Wave objects -------- */

class Wave;

class SimpleBinary : public Wave {
public:
    SimpleBinary(double freq, double initphi, double inc, double amp, double d, double a, double p);
};

%apply double *NUMPY_ARRAY_DOUBLE { double *hpa };
%apply double *NUMPY_ARRAY_DOUBLE { double *hca };

class InterpolateMemory : public Wave {
public:
    InterpolateMemory(double *hpa, double *hca, long samples, double samplingtime, double lookback, double d, double a, double p);
};

/* -------- TDI objects -------- */

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
    virtual double Ym(double t);
    virtual double Zm(double t);

    virtual double X1(double t);

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

%apply double PYTHON_SEQUENCE_DOUBLE[ANY] {double stproof[6], double sdproof[6], double stshot[6], double sdshot[6], double stlaser[6], double sdlaser[6], double claser[6]}

%apply Noise *PYTHON_SEQUENCE_NOISE[ANY] {Noise *proofnoise[6], Noise *shotnoise[6], Noise *lasernoise[6]}

class TDInoise : public TDI {
public:
    TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot, double stlaser, double sdlaser);

    TDInoise(LISA *mylisa, double stproof[6], double sdproof[6], double stshot[6], double sdshot[6], double stlaser[6], double sdlaser[6]);

    TDInoise(LISA *mylisa, Noise *proofnoise[6],Noise *shotnoise[6],Noise *lasernoise[6]);

    ~TDInoise();

    void setphlisa(LISA *mylisa);
    
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

/* -------- Helper functions -------- */

/* LISA */

%newobject stdlisa;
extern LISA *stdlisa();

/* Noise */

%newobject stdproofnoise;
extern Noise *stdproofnoise(LISA *lisa,double stproof,double sdproof);

%newobject stdopticalnoise;
extern Noise *stdopticalnoise(LISA *lisa,double stshot,double sdshot);

%newobject stdlasernoise;
extern Noise *stdlasernoise(LISA *lisa,double stlaser,double sdlaser);

%newobject newstdlasernoise;
extern Noise *newstdlasernoise(LISA *lisa,double stlaser,double sdlaser,int window);

%newobject stdnoise;
extern TDInoise *stdnoise(LISA *mylisa);

/* TDI interface (defined in lisasim-swig.i) */

extern void printtdi(char *filename,TDI *mytdi,int samples,double samplingtime,char *observables);

%apply double *NUMPY_ARRAY_DOUBLE { double *array };
extern void settdi(double *array,TDI *mytdi,int samples,double samplingtime,char *observables);

%apply double *NUMPY_ARRAY_DOUBLE { double *aa };
%apply double *NUMPY_ARRAY_DOUBLE { double *ab };
%apply double *NUMPY_ARRAY_DOUBLE { double *ag };
%apply double *NUMPY_ARRAY_DOUBLE { double *ax };
extern void setabg(double *aa, double *ab, double *ag, TDI *mytdi,int samples,double samplingtime);
extern void setabgx(double *aa, double *ab, double *ag, double *ax, TDI *mytdi,int samples,double samplingtime);

extern double retardation(LISA *mylisa,int ret1,int ret2,int ret3,int ret4,int ret5,int ret6,int ret7,int ret8,double t);

extern void printp(LISA *lisa,double t);
extern void printn(LISA *lisa,double t);
