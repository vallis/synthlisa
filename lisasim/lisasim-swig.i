/* File : lisasim-swig.i */

// Directors are not needed (only for AnyLISA, deprecated)
// %module(directors="1") lisaswig

%module lisaswig
%{
#include "lisasim.h"
%}

%include lisasim-typemaps.i

%define initsave(theclass)
%feature("addtofunc") theclass::theclass {
        self.initargs = args
}
%enddef

/* -------- LISA objects -------- */

%nodefault LISA;

class LISA {
 public:
    virtual void putp(Vector &outvector, int arm, double t);
    virtual void putn(Vector &outvector, int arm, double t);
};

class OriginalLISA : public LISA {
public:
    
    // accept the armlength in seconds

    OriginalLISA(double arm1,double arm2,double arm3);

    virtual void putn(Vector &outvector, int arm, double t);
    virtual void putp(Vector &outvector, int craft, double t);

    virtual double armlength(int arm, double t);
};

class ModifiedLISA : public OriginalLISA {
public:

    // accept the armlength in seconds
    
    ModifiedLISA(double arm1,double arm2,double arm3);

    void putp(Vector &outvector, int craft, double t);
        
    double armlength(int arm, double t);
};

class CircularRotating : public LISA {
 public:
	
    CircularRotating(double eta0=0.0,double xi0=0.0,double sw=0.0,double t0=0.0);
    CircularRotating(double myL,double eta0,double xi0,double sw,double t0);

    void putp(Vector &outvector, int craft, double t);
    
    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
    
    double genarmlength(int arms, double t);
};

initsave(PyLISA)

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

    void putp(Vector &outvector,int craft,double t);
    
    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
    
    double genarmlength(int arm, double t);
};

initsave(NoisyLISA)

class NoisyLISA : public LISA {
public:
    NoisyLISA(LISA *clean,double starm,double sdarm);

    // nontrivial destructor should appear here

    ~NoisyLISA(); 

    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

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

    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

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

    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

initsave(CacheLISA)

class CacheLISA : public LISA {
public:
    CacheLISA(LISA *basic);
    ~CacheLISA(); 

    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

initsave(MeasureLISA)

class MeasureLISA : public LISA {
 public:
    MeasureLISA(LISA *clean,double starm,double sdarm,int swindow = 1);
    ~MeasureLISA();
        
    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

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

/* Who gets ownership of the Numpy arrays? Should I worry about this
   and fix it with %feature("shadow")? Maybe so. */

initsave(InterpolateNoise)

// should add machinery to infer its own length...

%apply double *NUMPY_ARRAY_DOUBLE { double *noisebuf };

class InterpolateNoise : public Noise {
 public:
    InterpolateNoise(double sampletime,double prebuffer,double density,double exponent,int swindow = 1);

    InterpolateNoise(double *noisebuf,long samples,double sampletime,double prebuffer,double density, double exponent, int swindow = 1);

    ~InterpolateNoise();
    
    void reset();
    double noise(double time);
    double noise(double timebase,double timecorr);

    void setinterp(int window);
};

/* -------- Wave objects -------- */

class WaveObject;

%nodefault Wave;

class Wave : public WaveObject {
 public:
    Vector k;
    Tensor pp, pc;

    void putwave(Tensor &outtensor, double t);

    virtual double hp(double t);
    virtual double hc(double t);  
};

class SimpleBinary : public Wave {
 public:
    SimpleBinary(double freq, double initphi, double inc, double amp, double d, double a, double p);

    double hp(double t);
    double hc(double t);
};

class SimpleMonochromatic : public Wave {
public:
    SimpleMonochromatic(double freq, double phi, double gamma, double amp, double d, double a, double p);

    double hp(double t);
    double hc(double t);
};

initsave(NoiseWave)

%apply double *NUMPY_ARRAY_DOUBLE { double *hpa };
%apply double *NUMPY_ARRAY_DOUBLE { double *hca };

class NoiseWave : public Wave {
    public:
        // with Noise objects, given directly
	NoiseWave(Noise *noisehp, Noise *noisehc, double d, double a, double p);

        // allocates its own Noise objects
	NoiseWave(double sampletime, double prebuffer, double density, double exponent, int swindow, double d, double a, double p);
        
        // from sampled buffers (using filters and interpolation...)
	NoiseWave(double *hpa, double *hca, long samples, double sampletime, double prebuffer, double density, double exponent, int swindow, double d, double a, double p);

	~NoiseWave();

	double hp(double t);
	double hc(double t);
};

%newobject SampledWave;
extern NoiseWave *SampledWave(double *hpa, double *hca, long samples, double sampletime, double prebuffer, double density, double exponent, int swindow, double d, double a, double p);

// this class is now redundant, since it can be obtained as a special form of NoiseWave...

initsave(InterpolateMemory)

%apply double *NUMPY_ARRAY_DOUBLE { double *hpa };
%apply double *NUMPY_ARRAY_DOUBLE { double *hca };

class InterpolateMemory : public Wave {
public:
    InterpolateMemory(double *hpa, double *hca, long samples, double samplingtime, double lookback, double d, double a, double p);

    double hp(double t);
    double hc(double t);
};

initsave(PyWave)

%apply PyObject* PYTHONFUNC { PyObject *hpf, PyObject *hcf };

class PyWave : public Wave {
  public:
    PyWave(PyObject *hpf, PyObject *hcf, double d, double a, double p);

    double hp(double t);
    double hc(double t);
};

initsave(WaveArray)

class WaveArray : public WaveObject {
	public:
		WaveArray(Wave **WaveSeq, int WaveNum);
		~WaveArray();

		Wave *firstwave();
		Wave *nextwave();
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

/* We're holding on to the constructor args so that the LISA object
   won't get destroyed if they fall out of scope: we may still need
   them! */

initsave(TDIquantize)

class TDIquantize : public TDI {
 public:
    TDIquantize(TDI *bt,double qlev,int qbits,int qsat);
    ~TDIquantize();

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

%apply double PYTHON_SEQUENCE_DOUBLE[ANY] {double stproof[6], double sdproof[6], double stshot[6], double sdshot[6], double stlaser[6], double sdlaser[6], double claser[6]}

%apply Noise *PYTHON_SEQUENCE_NOISE[ANY] {Noise *proofnoise[6], Noise *shotnoise[6], Noise *lasernoise[6]}

/* We're holding on to the constructor args so that the LISA/Noise
   objects won't get destroyed if they fall out of scope: we may still
   need them! */

initsave(TDInoise)

/* Same for physical LISA objects */

%feature("addtofunc") TDInoise::setphlisa {
        self.phlisa = args
}

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

/* We're holding on to the constructor args so that the LISA/Noise
   objects won't get destroyed if they fall out of scope: we may still
   need them! */

initsave(TDIaccurate)

class TDIaccurate : public TDInoise {
 public:
    TDIaccurate(LISA *mylisa, Noise *proofnoise[6],Noise *shotnoise[6],Noise *lasernoise[6]);
    ~TDIaccurate();

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);
};

/* We're holding on to the constructor args so that the LISA/Wave
   objects won't get destroyed if they fall out of scope: we may still
   need them for TDInoise! */

initsave(TDIsignal)

class TDIsignal : public TDI {
public:
    TDIsignal(LISA *mylisa, WaveObject *mywave);
            
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
