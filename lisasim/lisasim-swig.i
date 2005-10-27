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

%feature("docstring") LISA "
LISA is the base class for all LISA geometries. All derived LISA
classes must define putn(l,t), putp(i,t), and armlength(l,t). LISA,
and most of its derived classes, are defined in lisasim-lisa.cpp"

%feature("docstring") LISA::putn "
LISA.putn(l,t) -> (nlx,nly,nlz) returns a 3-tuple with the SSB
components (normalized) of arm l (1,2,3,-1,-2,-3) at time t (seconds)."

%feature("docstring") LISA::putp "
LISA.putp(i,t) -> (pix,piy,piz) returns a 3-tuple with the SSB
coordinates of spacecraft i (1,2,3) at time t (seconds)."

%feature("docstring") LISA::armlength "
LISA.armlength(l,t) the armlength (seconds) of LISA link l
(1,2,3,-1,-2,-3) at time t (seconds)."

%nodefault LISA;

class LISA {
  public:
    virtual void putp(Vector &outvector, int arm, double t);
    virtual void putn(Vector &outvector, int arm, double t);

    virtual double armlength(int arm, double t);
};

%feature("docstring") OriginalLISA "
OriginalLISA(L1,L2,L3) returns a static LISA object with armlengths
given by L1, L2, L3 (in seconds). If omitted, the Li take the standard
value Lstd.

The positions of the spacecraft are chosen so that 1) the baricenter
is in the SSB origin; 2) piz = 0 for all spacecraft; 3) the spacecraft
sequence 1 -> 2 -> 3 -> 1 traces a clockwise path, as seen from above;
p1y = 0."

class OriginalLISA : public LISA {
  public:
    OriginalLISA(double arm1 = Lstd,double arm2 = Lstd,double arm3 = Lstd);

    virtual void putn(Vector &outvector, int arm, double t);
    virtual void putp(Vector &outvector, int craft, double t);

    virtual double armlength(int arm, double t);
};

%feature("docstring") ModifiedLISA "
ModifiedLISA(L1,L2,L3) returns a stationary LISA object set equal to
OriginalLISA(L1,L2,L3) at time 0, but rotating around its baricenter
(at the SSB origin) with angular speed 2*pi/yr. L1, L2, L3 are given
in seconds; if omitted, they take the standard value Lstd.

Because of the rotation, the armlengths depend on the
direction of light propagation."

class ModifiedLISA : public OriginalLISA {
  public:
    ModifiedLISA(double arm1 = Lstd,double arm2 = Lstd,double arm3 = Lstd);

    void putp(Vector &outvector, int craft, double t);
        
    double armlength(int arm, double t);
};

%feature("docstring") CircularRotating "
CircularRotating(eta0=0,xi0=0,sw=1,t0=0) and
CircularRotating(myL,eta0,xi0,sw,t0) return rigid LISA objects where
the armlengths are equal and constant before aberration is taken into
account. The parameters eta0 and xi0 are the true anomaly of the LISA
guiding center and the initial phase of the LISA array at t=t0; sw<0
will swap spacecraft 2 and 3, so that the spacecraft sequence 1 -> 2
-> 3 -> 1 goes ccw as seen from above. If given (in which case all
parameters must be specified), myL is the common armlength; if not
given, it is set to Lstd.

CircularRotating() is legal and replaces stdlisa(). The LISA
baricenter moves on a circular orbit in the ecliptic, while the
constellation rotates around the guiding center with the same angular
velocity. The orbits follow Krolak et al, PRD 70, 022003 (2004)."

class CircularRotating : public LISA {
  public:
    CircularRotating(double eta0=0.0, double xi0=0.0, double sw=0.0, double t0=0.0);
    CircularRotating(double myL, double eta0, double xi0, double sw, double t0);

    void putp(Vector &outvector, int craft, double t);
    
    double armlength(int arm, double t);

    // FIX: not clear to me that I want to expose (or even need) the following

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
    
    double genarmlength(int arms, double t);
};

%feature("docstring") EccentricInclined "
EccentricInclined(eta0=0,xi0=0,sw=1,t0=0) returns a realistic LISA
geometry modeled up to second order in the eccentricity, following
Cornish and Rubbo, PRD 67, 022001 (2003), but with the approximate
parametrization of CircularRotating (eta0 and xi0 true anomaly of
baricenter and array phase at t0=0; sw<0 swaps spacecraft)." 

class EccentricInclined : public LISA {
 public:
    EccentricInclined(double eta0=0.0, double xi0=0.0, double sw=1.0, double t0=0.0);
    // DEV: do I want an enhanced version that can set myL, as CircularRotating?

    void putp(Vector &outvector,int craft,double t);
    
    double armlength(int arm, double t);

    // FIX: not clear to me that I want to expose (or even need) the following

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
    
    double genarmlength(int arm, double t);
};

%feature("docstring") PyLISA "
PyLISA(baseLISA,Func) returns a LISA object with physical geometry
(used for light propagation and GW and noise responses) given by the
LISA object baseLISA, but with nominal armlengths (as used for the TDI
delays) given by the Python function Func(i,t) (i=1,2,3,-1,-2,-3, time
in seconds).

PyLISA attempts to provide a more general (if less efficient)
mechanism for the nominal-armlength computations previously performed
with NoisyLISA, NominalLISA, LinearLISA, MeasureLISA (all
experimental and now removed from the main distribution)."

initsave(PyLISA)

%apply PyObject* PYTHONFUNC { PyObject *func };

class PyLISA : public LISA {
  public:
    LISA *baseLISA;

    PyLISA(LISA *base,PyObject *func);

    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector,int craft,double t);

    double armlength(int arm, double t);

    // FIX: not clear to me that I want to expose (or even need) the following

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

%feature("docstring") CacheLISA "
CacheLISA(baseLISA) works by routing all putp-putn calls to the
baseLISA object, passed on construction. It will however interpose a
layer of its own making for retard() calls, effectively caching
chained retardations for the most recently accessed time.

CacheLISA, defined in lisasim-retard.h, might improve performance for
complicated (or multiple) TDI-variable evaluations, especially when
performed in conjunction with PyLISA."

initsave(CacheLISA)

class CacheLISA : public LISA {
  public:
    CacheLISA(LISA *base);
    ~CacheLISA();

    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

    double armlength(int arm, double t);

    // FIX: not clear to me that I want to expose (or even need) the following

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

// FIX: Who gets ownership of the Numpy arrays? Should I worry about
//      this and fix it with %feature("shadow")? Maybe so.
//      I wrote to Scott Ransom about it, but never got a response.

%feature("docstring") InterpolateNoise "
InterpolateNoise(gentime,prebuf,psd,exp,swin=1) and
InterpolateNoise(array,samples,sampletime,prebuf,norm,filter=0,window=1)
return a Noise object from either filtered and interpolated
pseudorandom whitenoise (first form) or from a filtered and
interpolated array (second form).

For the first form:

- gentime is the time spacing (seconds) of the pseudorandom samples at
  generation ;

- prebuf (seconds??) sets the minimum interval of buffering of the
  pseudorandom sequences (after the noise has been requested at time
  t, the earliest noise value guaranteed to be available will be at
  time (t-prebuf)

- psd and exp set the one-sided PSD of the filtered pseudorandom noise
  to psd*(f/Hz)^exp Hz^-1 (currently exp = -2, 0, 2 are implemented);

- wind (> 1) sets the semiwidth of the data window used in Lagrange
  interpolation (1 yields linear interpolation)."

initsave(InterpolateNoise)

// DEV: Should add machinery to let the typemaps infer the length of
// the array and pass it on

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
    SimpleBinary(double freq, double initphi, double inc, double amp, double elat, double elon, double p);

    double hp(double t);
    double hc(double t);
};

class SimpleMonochromatic : public Wave {
public:
    SimpleMonochromatic(double freq, double phi, double gamma, double amp, double elat, double elon, double p);

    double hp(double t);
    double hc(double t);
};

initsave(NoiseWave)

%apply double *NUMPY_ARRAY_DOUBLE { double *hpa };
%apply double *NUMPY_ARRAY_DOUBLE { double *hca };

class NoiseWave : public Wave {
    public:
        // with Noise objects, given directly
	NoiseWave(Noise *noisehp, Noise *noisehc, double elat, double elon, double p);

        // allocates its own Noise objects
	NoiseWave(double sampletime, double prebuffer, double density, double exponent, int swindow, double elat, double elon, double p);
        
        // from sampled buffers (using filters and interpolation...)
	NoiseWave(double *hpa, double *hca, long samples, double sampletime, double prebuffer, double density, double exponent, int swindow, double elat, double elon, double p);

	~NoiseWave();

	double hp(double t);
	double hc(double t);
};

%newobject SampledWave;
extern NoiseWave *SampledWave(double *hpa, double *hca, long samples, double sampletime, double prebuffer, double density, double exponent, int swindow, double elat, double elon, double p);

// this class is now redundant, since it can be obtained as a special form of NoiseWave...

initsave(InterpolateMemory)

%apply double *NUMPY_ARRAY_DOUBLE { double *hpa };
%apply double *NUMPY_ARRAY_DOUBLE { double *hca };

class InterpolateMemory : public Wave {
public:
    InterpolateMemory(double *hpa, double *hca, long samples, double samplingtime, double lookback, double elat, double elon, double p);

    double hp(double t);
    double hc(double t);
};

initsave(PyWave)

%apply PyObject* PYTHONFUNC { PyObject *hpf, PyObject *hcf };

class PyWave : public Wave {
  public:
    PyWave(PyObject *hpf, PyObject *hcf, double elat, double elon, double p);

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
