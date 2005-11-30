/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

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

%define initdoc(theclass)
%feature("docstring") theclass::theclass "
Constructor. See above.
"
%feature("docstring") theclass::~theclass "
Destructor. See above.
"
%enddef

%define exceptionhandle(thefunction,theexception,theerror)
%exception thefunction {
	try {
		$action
	} catch (theexception &e) {
		PyErr_SetString(theerror,"");
		return NULL;
	}
};
%enddef

/* -------- LISA objects -------- */


%feature("docstring") LISA "
LISA is the base class for all LISA geometries. All derived LISA
classes must define putn(l,t), putp(i,t), and armlength(l,t). LISA,
and most of its derived classes, are defined in lisasim-lisa.cpp"

%feature("docstring") LISA::putn "
LISA.putn(l,t) -> (nlx,nly,nlz) returns a 3-tuple with the SSB
components (normalized) of arm l (1,2,3,-1,-2,-3) at time t [s]."

%feature("docstring") LISA::putp "
LISA.putp(i,t) -> (pix,piy,piz) returns a 3-tuple with the SSB
coordinates of spacecraft i (1,2,3) at time t [s]."

%feature("docstring") LISA::armlength "
LISA.armlength(l,t) the armlength [s] of LISA link l (1,2,3,-1,-2,-3)
at time t [s]."

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

initdoc(OriginalLISA)

class OriginalLISA : public LISA {
  public:
    OriginalLISA(double arm1 = Lstd,double arm2 = Lstd,double arm3 = Lstd);
};


%feature("docstring") ModifiedLISA "
ModifiedLISA(L1,L2,L3) returns a stationary LISA object set equal to
OriginalLISA(L1,L2,L3) at time 0, but rotating around its baricenter
(at the SSB origin) with angular speed 2*pi/yr. L1, L2, L3 are given
in seconds; if omitted, they take the standard value Lstd.

Because of the rotation, the armlengths depend on the
direction of light propagation."

initdoc(ModifiedLISA)

class ModifiedLISA : public OriginalLISA {
  public:
    ModifiedLISA(double arm1 = Lstd,double arm2 = Lstd,double arm3 = Lstd);
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

initdoc(CircularRotating)

class CircularRotating : public LISA {
  public:
    CircularRotating(double eta0=0.0, double xi0=0.0, double sw=0.0, double t0=0.0);
    CircularRotating(double myL, double eta0, double xi0, double sw, double t0);
};


%feature("docstring") EccentricInclined "
EccentricInclined(eta0=0,xi0=0,sw=1,t0=0) returns a realistic LISA
geometry modeled up to second order in the eccentricity, following
Cornish and Rubbo, PRD 67, 022001 (2003), but with the approximate
parametrization of CircularRotating (eta0 and xi0 true anomaly of
baricenter and array phase at t0=0; sw<0 swaps spacecraft)." 

initdoc(EccentricInclined)

class EccentricInclined : public LISA {
 public:
    EccentricInclined(double eta0=0.0, double xi0=0.0, double sw=1.0, double t0=0.0);
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

initdoc(PyLISA)

initsave(PyLISA)

%apply PyObject* PYTHONFUNC { PyObject *func };

class PyLISA : public LISA {
  public:
    LISA *baseLISA;

    PyLISA(LISA *base,PyObject *func);
};


%feature("docstring") CacheLISA "
CacheLISA(baseLISA)
returns a LISA object that works by routing all putp-putn calls to the
baseLISA object, passed on construction. It will however interpose a
layer of its own making for retard() calls, effectively caching
chained retardations for the most recently accessed time.

CacheLISA, defined in lisasim-retard.h, might improve performance for
complicated (or multiple) TDI-variable evaluations, especially when
performed in conjunction with PyLISA."

initdoc(CacheLISA)

initsave(CacheLISA)

class CacheLISA : public LISA {
  public:
    CacheLISA(LISA *base);
};


%feature("docstring") SampledLISA "
SampledLISA(p1,p2,p3,deltat,prebuffer,interp)
returns a LISA object that takes the positions of its spacecraft from
the 2-dimensional Numeric arrays p1, p2, p3; each of these consists of
three columns that give the SSB coordinates of corresponding LISA
spacecraft.

The array data is understood to be spaced by intervals deltat [s], and
it is offset so that (for instance)

p1x(t) = p1[(t - prebuffer)/deltat,0],
p1y(t) = p1[(t - prebuffer)/deltat,1],
p1z(t) = p1[(t - prebuffer)/deltat,2].
  
Last, interp (> 1) sets the semiwidth of the data window used in
Lagrange interpolation of the positions.

Note 1: SampledLISA makes copies of the positions arrays, so these can
be safely destroyed after calling the SampledLISA constructor; on the
other hand, modifying the arrays passed to the constructor will have no
effect on the positions returned by SampledLISA.putp().

Note 2: currently armlength are computed explicitly by SampledLISA by
solving the backward light propagation equation."

initdoc(SampledLISA)

/* SampledLISA makes copies of the positions arrays, so initsave is not
   needed */

class SampledLISA : public LISA {
 public:
    SampledLISA(double *numarray,long length,double *numarray,long length,double *numarray,long length,double deltat,double prebuffer,int interplen = 1);
};


/* -------- Signal/Noise objects -------- */

exceptionhandle(SignalSource::__getitem__,ExceptionOutOfBounds,PyExc_IndexError)

%nodefault SignalSource;
class SignalSource {
 public:
	virtual void reset(unsigned long seed = 0);

	%extend {
		double __getitem__(long pos) {
			return (*self)[pos];
		};
	};
};

class WhiteNoiseSource : public SignalSource {
 public:
	WhiteNoiseSource(long len,unsigned long seed = 0,double norm = 1.0);	
};

class SampledSignalSource : public SignalSource {
 public:
	SampledSignalSource(double *numarray,long length,double norm = 1.0);
};

%nodefault Filter;
class Filter;

class NoFilter : public Filter {
 public:
	NoFilter();
};

class IntFilter : public Filter {
 public:
    IntFilter(double a = 0.9999);
};

class DiffFilter : public Filter {
 public:
	DiffFilter();
};

class FIRFilter: public Filter {
 public:
 	FIRFilter(double *doublearray,int doublenum);
};

class IIRFilter : public Filter {
 public:
 	IIRFilter(double *doublearray,int doublenum,double *doublearray,int doublenum);
};

class SignalFilter : public SignalSource {
  public:
	SignalFilter(long length,SignalSource *src,Filter *flt,double norm = 1.0);
};

exceptionhandle(Signal::value,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(Signal::__call__,ExceptionOutOfBounds,PyExc_IndexError)

%nodefault Signal;
class Signal {
 public:
    virtual void reset(unsigned long seed = 0);

	virtual double value(double time);
	virtual double value(double timebase,double timecorr);
	
	%extend {
		double __call__(double time) {
			return self->value(time);
		};
	};
};

typedef Signal Noise;

%nodefault Interpolator;
class Interpolator;

class NearestInterpolator : public Interpolator {
 public:
	NearestInterpolator();
};

class LinearInterpolator : public Interpolator {
 public:
	LinearInterpolator();
};

class LinearExtrapolator : public Interpolator {
 public:
	LinearExtrapolator();
};

class LagrangeInterpolator : public Interpolator {
 public:
    LagrangeInterpolator(int semiwin);
};

extern Interpolator *getInterpolator(int interplen);

class InterpolatedSignal : public Signal {
 public:
	InterpolatedSignal(SignalSource *src,Interpolator *interp,
					   double deltat,double prebuffer = 0.0,double norm = 1.0);
	
	void setinterp(Interpolator *inte);
};


%feature("docstring") PowerLawNoise "
PowerLawNoise(deltat,prebuffer,psd,exponent,interplen = 1,seed = 0)
"

initdoc(PowerLawNoise)

class PowerLawNoise : public Signal {
 public:
	PowerLawNoise(double deltat,double prebuffer,double psd,double exponent,
				  int interplen = 1,unsigned long seed = 0);
};


%feature("docstring") SampledSignal "
SampledSignal(numarray,deltat,prebuffer,norm = 1.0,filter = 0,interplen = 1)
"

initdoc(SampledSignal)

initsave(SampledSignal)

class SampledSignal : public Signal {
 public:
	SampledSignal(double* numarray,long length,double deltat,double prebuffer,
				  double norm = 1.0,Filter *filter = 0,int interplen = 1);
};


/* -------- Wave objects -------- */


%feature("docstring") Wave "
Wave is the base class for all GW objects; all Wave objects represent
plane waves incoming from a definite position in the sky (as
characterized by its SSB ecliptic latitude and longitude) with a
definite polarization angle.

- Wave provides putk(), which returns the GW direction of propagation,
  and putwave(t), which returns the gravitational strain tensor at
  time t at the SSB.

- All derived Wave objects must define hp(t) and hc(t), which describe
  the time dependence of the two independent GW polarizations.

Wave, and most of its derived classes, are defined in lisasim-wave.cpp"

%feature("docstring") Wave::putk "
Wave.putk() -> (kx,ky,kz) returns a 3-tuple with the SSB components
(normalized) of GW propagation vector."

%feature("docstring") Wave::putwave "
Wave.putwave(t) -> ((hxx,hxy,hxz),(hyx,hyy,hyz),(hzx,hzy,hzz)) returns
a nested 3-tuple with the gravitational strain tensor at time t [s] at
the SSB."

%feature("docstring") Wave::hp "
Wave.hp(t) returns the hp polarization of the GW Wave at time t [s] at
the SSB."

%feature("docstring") Wave::hc "
Wave.hc(t) returns the hc polarization of the GW Wave at time t [s] at
the SSB."

%feature("docstring") Wave::putep "
Wave.putep(elat,elon,pol) returns the basic ep polarization tensor
for a plane Wave object at ecliptic latitude elat, longitude elon,
and polarization pol. This is a class (static) method.
"

%feature("docstring") Wave::putec "
Wave.putec(elat,elon,pol) returns the basic ec polarization tensor
for a plane Wave object at ecliptic latitude elat, longitude elon,
and polarization pol. This is a class (static) method.
"

class WaveObject;

%nodefault Wave;

class Wave : public WaveObject {
 public:
    void putk(Vector &outvector);
    void putwave(Tensor &outtensor, double t);

    virtual double hp(double t);
    virtual double hc(double t);
    
    static void putep(Tensor &outtensor,double b,double l,double p);
    static void putec(Tensor &outtensor,double b,double l,double p);
};


%feature("docstring") SimpleBinary "
SimpleBinary(f,phi0,inc,amp,elat,elon,pol) returns a Wave object that
models the waveform emitted by a simple monochromatic binary with the
following parameters:

- the wave has frequency f [Hz] (so the binary has orbital frequency
  f/2) and initial phase phi0;

- the hp and hc polarizations have a relative phase shift pi/2 and
  have an amplitude amp distributed consistently with a binary
  inclination angle inc:
  hp = amp * (1 + cos(inc)^2) * cos(2*pi*f*t + phi0),
  hc = amp * 2*cos(inc) *       sin(2*pi*f*t + phi0);

- the wave is incoming from sky position (elat,elon), with
  polarization pol."

initdoc(SimpleBinary)

class SimpleBinary : public Wave {
 public:
    SimpleBinary(double freq, double phi0, double inc, double amp, double elat, double elon, double pol);
};


%feature("docstring") SimpleMonochromatic "
SimpleMonochromatic(f,phi,gamma,amp,elat,elon,pol) returns a Wave
object that implements a simple sinusoidal wave with the following
parameters:

- the wave has frequency f [Hz];

- the hp and hc polarizations have a relative phase shift phi and
  amplitude amp distributed according to the angle gamma:
  hp = amp * sin(gamma) * sin(2*pi*f*t + phi),
  hc = amp * cos(gamma) * sin(2*pi*f*t);

- the wave is incoming from sky position (elat,elon), with
  polarization pol."

initdoc(SimpleMonochromatic)

class SimpleMonochromatic : public Wave {
 public:
    SimpleMonochromatic(double freq, double phi, double gamma, double amp, double elat, double elon, double pol);
};


%feature("docstring") GaussianPulse "
GaussianPulse(t0,efold,gamma,amp,elat,elon,pol) returns a Wave object
that implements a Gaussian wavepulse with the following parameters:

- the pulse is centered at SSB baricentric time t0 [seconds];

- the pulse has e-folding time efold [seconds];

- the pulse has amplitude amp, distributed as {ap,ac} = 
  amp*{sin(gamma),cos(gamma)} between the two polarizations
  (the same convention as SimpleMonochromatic)

- the pulse is incoming from sky position (elat,elon), with
  polarization pol.

The amplitude of the pulse is cut to zero at 10 efolding times from
the central time (this is set by GaussianPulse::sigma_cutoff in
lisasim-wave.cpp)."

initdoc(GaussianPulse)

class GaussianPulse : public Wave {
 public:
    GaussianPulse(double time, double decay, double gamma, double amp, double b, double l, double p);
};


%feature("docstring") SineGaussian "
SineGaussian(t0,efold,f,phi0,gamma,amp,elat,elon,pol) returns a Wave object
that implements a Sine-Gaussian wavepulse with the following parameters:

- the pulse is centered at SSB baricentric time t0 [seconds];

- the pulse has e-folding time efold [seconds];

- the pulse is modulated by a sinusoid of frequency f, centered at time
  t0, with relative phase phi0 between the two polarizations;

- the pulse has amplitude amp, distributed as {ap,ac} = 
  amp*{sin(gamma),cos(gamma)} between the two polarizations
  (the same convention as SimpleMonochromatic)

- the pulse is incoming from sky position (elat,elon), with
  polarization pol.

The amplitude of the pulse is cut to zero at 10 efolding times from
the central time (this is set by SineGaussian::sigma_cutoff in
lisasim-wave.cpp)."

initdoc(SineGaussian)

class SineGaussian : public Wave {
 public:
    SineGaussian(double time, double decay, double freq, double phase0, double gamma, double amp, double b, double l, double p);
};


%feature("docstring") NoiseWave "
NoiseWave(hpnoise,hcnoise,elat,elon,pol),
NoiseWave(deltat,prebuffer,psd,filter,interp,elat,elon,pol)
return Wave objects that represent a noise-like plane GW incoming from
sky position (elat,elon) with polarization pol.

With the first form, hpnoise and hcnoise are Noise objects (previously
instantiated, e.g., with InterpolateNoise) that provide the GW
polarizations directly from the methods hpnoise.noise(t) and
hpnoise.noise(t)

With the second form, hp(t) and hc(t) are fed from two independent
pseudorandom noise streams (implemented with InterpolateNoise), where:

- deltat [s] is the time spacing of the pseudorandom samples at
  generation;

- prebuffer [s ???] sets the minimum buffering interval for the
  pseudorandom sequences, which are stored in ring buffers (i.e.,
  after pseudorandom noise has been requested and computed at time t,
  the earliest value guaranteed to be available will be at time t -
  prebuffer);

- psd and filter set the one-sided PSD of the pseudorandom noise
  stream to psd*(f/Hz)^filter Hz^-1 (currently only filter = -2, 0, 2
  is implemented);

- interp (> 1) sets the semiwidth of the data window used in Lagrange
  interpolation (1 yields linear interpolation).

The parameter prebuffer should generally be set higher for a NoiseWave
object than it would be for an InterpolateNoise object, because hp(t)
and hc(t) are defined in terms of SSB time, so the polarization values
needed to build (say) X at t = 0 may require accessing hp and hc at
the SSB t = -500 s, or earlier; on the other hand, when Synthetic LISA
builds new InterpolateNoise objects it fills the ring buffer with
values in the interval [-prebuffer,0].

NoiseWave can also be created from arrays of previously
generated waveform samples, using the SampledWave syntax (see help for
SampledWave):

NoiseWave(hparray,hcarray,len,deltat,prebuffer,norm,filter,interp,
          elat,elon,pol)

TODO: details to be addressed: does the choice of interpolation affect
the effective prebuffering interval? ???"

initdoc(NoiseWave)

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
	NoiseWave(double *hpa, double *hca, long samples, double sampletime, double prebuffer, double norm, Filter *filter, int swindow, double b, double l, double p);

    ~NoiseWave();

    double hp(double t);
    double hc(double t);
};


/* SampledWave is implemented as a Python factory function for
   NoiseWave; this is the best way (better than a C++ factory
   function) to assure that garbage collection does not disturb its
   arguments. */

// ??? It would be interesting not to have to give the length of the arrays...


%pythoncode %{
def SampledWave(hparray,hcarray,len,deltat,prebuffer,norm,filter,interp,elat,elon,pol):
    """Returns a Wave object that represent a sampled plane GW incoming from
sky position (elat,elon) with polarization pol, where:

- hparray and hcarray are 1-D Numeric arrays of length len containing
  time series for the hp and hc polarizations;

- deltat [s] is the time spacing of the time series;

- prebuffer [s ???] sets the indexing offset of the arrays with
  respect to the SSB time, so that (temporarily forgetting
  interpolation) hp(t) = hparray[(t - prebuffer)/deltat].
  This padding is needed because hp(t) and hc(t) are defined in terms
  of SSB time: the polarization values needed to build (say) X at t =
  0 may require accessing hp and hc at the SSB t = -500 s, or earlier;

- norm and filter set the normalization and filtering of the sampled
  arrays. The parameter filter must be a SynthLISA Filter object, or
  0 if only normalization is required;

- interp (> 1) sets the semiwidth of the data window used in Lagrange
  interpolation (1 yields linear interpolation).

SampledWave is represented internally using NoiseWave.

TODO: check if the prebuffering formula above is exact or if there is
another displacement by one or so."""

    return NoiseWave(hparray,hcarray,len,deltat,prebuffer,norm,filter,interp,elat,elon,pol)
%}


/* Provide backward compatibility to InterpolateMemory (essentially
   the same as SampledWave, without possibility of normalizing,
   filtering, or changing the interpolation scheme); note that the
   prebuffering convention might be slightly different. */

%pythoncode %{
def InterpolateMemory(hparray,hcarray,len,deltat,prebuffer,elat,elon,pol):
    """For backward compatibility, equivalent to 
NoiseWave(hparray,hcarray,len,deltat,prebuffer,1.0,0,1,elat,elon,pol)."""

    return NoiseWave(hparray,hcarray,len,deltat,prebuffer,1.0,0,1,elat,elon,pol)
%}

%feature("docstring") PyWave "
PyWave(hpfunc,hcfunc,elat,elon,pol)
returns a Wave object that represents a generic plane GW incoming from
ecliptic latitude elat and longitude elon, with polarization pol. The
parameters hpfunc(t) and hcfunc(t) are Python functions that must
return the hp and hc polarizations at SSB time t.

While not too efficient, PyWave may be the simplest way to extend the
Synthetic-LISA built-in Wave objects."

initdoc(PyWave)

initsave(PyWave)

%apply PyObject* PYTHONFUNC { PyObject *hpf, PyObject *hcf };

class PyWave : public Wave {
 public:
    PyWave(PyObject *hpf, PyObject *hcf, double elat, double elon, double p);

    double hp(double t);
    double hc(double t);
};


%feature("docstring") WaveArray "
WaveArray(waves[]) returns a WaveArray object that represents the
linear superposition of an array of Wave objects (each of which can be
of a different type, have different parameters, or have different sky
position and polarization angle).

Note that if WaveArray is initialized with a reference to a list,
successive modifications to the list will not be mirrored in the
contents of the WaveArray, and could result in single Wave objects
being deallocated when they shouldn't. To avoid the possibility of
this, WaveArray should be created with a tuple."

%feature("docstring") WaveArray::firstwave "
WaveArray.firstwave() returns the first Wave object in a WaveArray."

%feature("docstring") WaveArray::nextwave "
WaveArray.nextwave() can be used repeatedly after
WaveArray.firstwave() to return all the Wave objects in a
WaveArray. After the last object, it will return None."

initdoc(WaveArray)

initsave(WaveArray)

%exception WaveArray::WaveArray {
	try {
		$action
	} catch (ExceptionOutOfBounds &e) {
		PyErr_SetString(PyExc_ValueError,"Wrong arguments.");
		return NULL;
	}
};

class WaveArray : public WaveObject {
 public:
    WaveArray(Wave **WaveSeq, int WaveNum);
    ~WaveArray();

    Wave *firstwave();
    Wave *nextwave();
};


/* -------- TDI objects -------- */

exceptionhandle(TDI::X,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::Y,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::Z,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::alpha,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::beta,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::gamma,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::zeta,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::P,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::E,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::U,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::y,ExceptionOutOfBounds,PyExc_IndexError)

class TDI {
public:
    TDI() {};
    virtual ~TDI() {};

    virtual void reset() {};

    // ??? should I keep all of these?

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

/* Note that if TDInoise is initialized with a reference to a list,
   successive modifications to the list will not be mirrored in the
   contents of the TDInoise, and could result in unwanted behavior.
   To avoid the possibility of this, TDInoise should be created with a
   tuple. */

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

// ??? Helper functions used to be here...


