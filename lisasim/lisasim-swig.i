/* File : lisasim-swig.i */

%module(directors="1") lisaswig
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
};

class ModifiedLISA : public OriginalLISA {
public:

    // accept the armlength in seconds
    
    ModifiedLISA(double arm1,double arm2,double arm3);
        
    double armlength(int arm, double t);
};

class CircularRotating : public LISA {
 public:
	
    CircularRotating(double eta0=0.0,double xi0=0.0,double sw=0.0,double t0=0.0);
    CircularRotating(double myL,double eta0,double xi0,double sw,double t0);
    
    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
    
    double genarmlength(int arms, double t);
};

%feature("director") AnyLISA;
class AnyLISA : public LISA {
 public:
    AnyLISA(LISA *clean);
    virtual ~AnyLISA();

    virtual void reset();

    virtual double armlength(int arm, double t);

    virtual double armlengthbaseline(int arm, double t);
    virtual double armlengthaccurate(int arm, double t);
};

%apply PyObject* PYTHONFUNC { PyObject *func };

class PyLISA : public LISA {
  public:
    PyLISA(LISA *base,PyObject *func);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

class EccentricInclined : public LISA {
 public:

    EccentricInclined(double eta0 = 0.0,double xi0 = 0.0,double sw = 1.0,double t0=0.0);
    
    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
    
    double genarmlength(int arm, double t);
};

class NoisyLISA : public LISA {
public:
    NoisyLISA(LISA *clean,double starm,double sdarm);

    // nontrivial destructor should appear here

    ~NoisyLISA(); 

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

class NominalLISA : public LISA {
 public:
    double Lnom, emod, cmod, toff;

    NominalLISA(double eta0,double xi0,double sw,double t0);
    ~NominalLISA();

    void setparameters(double l,double cm,double em,double toff);
    void setparameters(double cm,double em,double toff);
    void setparameters3(double l,double cm,double em);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

%apply double PYTHON_SEQUENCE_DOUBLE[ANY] {double dl[6], double dldt[6]}

class LinearLISA : public LISA {
 public:
    LinearLISA(double eta0,double xi0,double sw,double t0);
    ~LinearLISA();

    void settimeoffset(double toff);
    void setparameters(double dl[6],double dldt[6]);
    
    double armlengtherror(int arm, double t);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

class CacheLISA : public LISA {
public:
    CacheLISA(LISA *basic);
    ~CacheLISA(); 

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

class MeasureLISA : public LISA {
 public:
    MeasureLISA(LISA *clean,double starm,double sdarm,int swindow = 1);
    ~MeasureLISA();
        
    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

/* -------- Noise objects -------- */

class Noise {
 public:
    virtual void reset() {};

    virtual double noise(double time) = 0;

    virtual double noise(double timebase,double timecorr);
};

%apply double *NUMPY_ARRAY_DOUBLE { double *noisebuf };

class InterpolateNoise : public Noise {
 public:
    InterpolateNoise(double sampletime,double prebuffer,double density,double exponent,int swindow = 1);

    InterpolateNoise(double *noisebuf,long samples,double sampletime,double prebuffer,double density, double exponent = 0.0, int swindow = 1);

    ~InterpolateNoise();
    
    void reset();
    double noise(double time);
    double noise(double timebase,double timecorr);

    void setinterp(int window);
};

/* -------- Wave objects -------- */

class Wave;

class SimpleBinary : public Wave {
public:
    SimpleBinary(double freq, double initphi, double inc, double amp, double d, double a, double p);
};

class SimpleMonochromatic : public Wave {
public:
    SimpleMonochromatic(double freq, double phi, double gamma, double amp, double d, double a, double p);
};

%apply double *NUMPY_ARRAY_DOUBLE { double *hpa };
%apply double *NUMPY_ARRAY_DOUBLE { double *hca };

class InterpolateMemory : public Wave {
public:
    InterpolateMemory(double *hpa, double *hca, long samples, double samplingtime, double lookback, double d, double a, double p);
};

%apply PyObject* PYTHONFUNC { PyObject *hpf, PyObject *hcf };

class PyWave : public Wave {
  public:
    PyWave(PyObject *hpf, PyObject *hcf, double d, double a, double p);
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

    virtual double alpham(double t);
    virtual double betam(double t);
    virtual double gammam(double t);

    virtual double alpha1(double t);
    virtual double alpha2(double t);
    virtual double alpha3(double t);
    
    virtual double zeta(double t);

    virtual double zeta1(double t);
    virtual double zeta2(double t);
    virtual double zeta3(double t);

    virtual double P(double t);
    virtual double E(double t);
    virtual double U(double t);
    
    virtual double Xm(double t);
    virtual double Ym(double t);
    virtual double Zm(double t);

    virtual double Xmlock1(double t);
    virtual double Xmlock2(double t);
    virtual double Xmlock3(double t);

    virtual double X1(double t);
    virtual double X2(double t);
    virtual double X3(double t);

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

    // just when are virtual destructors needed?
    virtual ~TDInoise();

    void setphlisa(LISA *mylisa);

    void lock(int master);
    
    void reset();

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

class TDIaccurate : public TDInoise {
 public:
    TDIaccurate(LISA *mylisa, Noise *proofnoise[6],Noise *shotnoise[6],Noise *lasernoise[6]);

    ~TDIaccurate();

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

class TDIsignal : public TDI {
public:
    TDIsignal(LISA *mylisa, Wave *mywave);
            
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

/* TDInoise */

%newobject stdnoise;
extern TDInoise *stdnoise(LISA *mylisa);

/* debugging */

%apply double *NUMPY_ARRAY_DOUBLE { double *array };
extern void settdi(double *array,TDI *mytdi,int samples,double samplingtime,char *observables);

%apply Noise *PYTHON_SEQUENCE_NOISE[ANY] {Noise *proofnoise[6], Noise *shotnoise[6], Noise *lasernoise[6]}

extern double retardation(LISA *mylisa,int ret1,int ret2,int ret3,int ret4,int ret5,int ret6,int ret7,int ret8,double t);
