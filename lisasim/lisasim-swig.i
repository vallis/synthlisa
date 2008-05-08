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

%pythoncode %{
import numpy.oldnumeric as Numeric

import math
import sys
import random
import time
%}

/* ------- random numbers ------- */

%pythoncode %{
globalseed = 0

cseeds = []

def setglobalseed(seed = 0):
    globalseed = seed
    random.seed(globalseed)
   
def getglobalseed():
    if globalseed == 0:
        setglobalseed(int(time.time()))

    return globalseed
    
def getcseed():
    if globalseed == 0:
        setglobalseed(int(time.time()))

    while True:
        ret = random.randint(0,2**32)
    
        if not ret in cseeds:
            cseeds.append(ret)
            return ret
%}

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

%feature("docstring") LISA::putv "
LISA.putv(i,t) -> (vix,viy,viz) returns a 3-tuple with the SSB
coordinate speed of spacecraft i (1,2,3) at time t [s], given
in units of the speed of light."

%feature("docstring") LISA::armlength "
LISA.armlength(l,t) the armlength [s] of LISA link l (1,2,3,-1,-2,-3)
for laser pulse reception at time t [s]."
    
%feature("docstring") LISA::dotarmlength "
LISA.armlength(l,t) the instantaneous rate of change of the armlength
of LISA link l (1,2,3,-1,-2,-3) for laser pulse reception a time t [s],
given in units of the speed of light."
    
%feature("docstring") LISA::reset "
LISA.reset() resets any underlying pseudo-random or ring-buffer
elements used by the LISA object."

%nodefault LISA;

class LISA {
  public:
    virtual void putp(Vector &outvector, int craft, double t);
    virtual void putn(Vector &outvector, int arm, double t);
    virtual void putv(Vector &outvector, int craft, double t);
    
    virtual double armlength(int arm, double t);
    virtual double dotarmlength(int arm, double t);
    
    virtual void reset();
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

initsave(OriginalLISA)

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

initsave(ModifiedLISA)

class ModifiedLISA : public OriginalLISA {
  public:
    ModifiedLISA(double arm1 = Lstd,double arm2 = Lstd,double arm3 = Lstd);
};

%nodefault ApproxLISA;

class ApproxLISA {
    public:
      virtual double geteta0();
      virtual double getxi0();
      virtual double getsw();
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

initsave(CircularRotating)

class CircularRotating : public LISA, public ApproxLISA {
  public:
    CircularRotating(double eta0=0.0, double xi0=0.0, double sw=0.0, double t0=0.0);
    CircularRotating(double myL,double e0,double x0,double sw,double t0);
};


%feature("docstring") EccentricInclined "
EccentricInclined(eta0=0,xi0=0,sw=1,t0=0) returns a realistic LISA
geometry modeled up to second order in the eccentricity, following
Cornish and Rubbo, PRD 67, 022001 (2003), but with the approximate
parametrization of CircularRotating (eta0 and xi0 true anomaly of
baricenter and array phase at t0=0; sw<0 swaps spacecraft)." 

initdoc(EccentricInclined)

initsave(EccentricInclined)

class EccentricInclined : public LISA, public ApproxLISA {
 public:
    EccentricInclined(double eta0=0.0, double xi0=0.0, double sw=1.0, double t0=0.0);

    double genarmlength(int arm,double t);
};


class ZeroLISA : public LISA {
  public:
    ZeroLISA();
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

%feature("docstring") AllPyLISA "
AllPyLISA(putpfunc,armlengthfunc = 0)

returns a LISA object that takes the positions of spacecraft from
the Python function putpfunc(craft,time), where craft = 1,2,3, and
time is given in seconds; and that takes the armlengths from the
Python function armlengthfunc(link,time), where link = 1,2,3,-1,-2,-3.
If armlengthfunc is not given, the armlengths will be determined from
the s/c positions using the exact light-propagation equation (solved
in the base LISA class).

Note: simulations that use AllPyLISA can be speeded up somewhat by
enclosing the AllPyLISA object in a CacheLengthLISA object."

initdoc(AllPyLISA)

initsave(AllPyLISA)

%apply PyObject* PYTHONFUNC { PyObject *sfunc, PyObject *afunc };

class AllPyLISA : public LISA {
  public:
    AllPyLISA(PyObject *sfunc,PyObject *afunc = 0);
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


%feature("docstring") CacheLengthLISA "
CacheLengthLISA(baseLISA,bufferlength,deltat,interplen = 1)
returns a LISA object that caches and interpolates armlengths found by
solving the light-propagation equation for the spacecraft positions
returned by baseLISA.putp(). The light-propagation equation is solved
every deltat seconds, and results remain available in a time window of
duration bufferlength*deltat. Last, interplen is the semiwidth of the
interpolation kernel (with 0 nearest-neighbor interpolation and 1 linear
interpolation).

Note: the current implementation does not support changing the physical
LISA of baseLISA on the fly, and is untested for different nominal and
physical LISAs."

initdoc(CacheLengthLISA)

initsave(CacheLengthLISA)

class CacheLengthLISA : public LISA {
 public:
    CacheLengthLISA(LISA *lisa,long length,double deltat,int interplen = 4);
    ~CacheLengthLISA();
};

extern double retardation(LISA *lisa,int ret1,int ret2,int ret3,int ret4,int ret5,int ret6,int ret7,int ret8,double t);

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


%feature("docstring") WhiteNoiseSource::setglobalseed "
WhiteNoiseSource.setglobalseed(seed) sets the global seed that will be
used to initialize the next pseudorandom-number generator to be created
or reset. The argument should be an unsigned long integer. The global
seed is increased by one after each creation or initialization. This is
a class (static) method. "

%feature("docstring") WhiteNoiseSource::getglobalseed "
WhiteNoiseSource.getglobalseed() returns the global seed that will be
used to initialize the next pseudorandom-number generator to be created
or reset. The seed is an unsigned long integer. The global seed is
increased by one after each creation or initialization. This is a class
(static) method. "

class WhiteNoiseSource : public SignalSource {
 public:
    WhiteNoiseSource(long len,unsigned long seed = 0,double norm = 1.0);    

    static void setglobalseed(unsigned long seed = 0);
    static unsigned long getglobalseed();
};

initsave(SampledSignalSource)

class SampledSignalSource : public SignalSource {
 public:
    SampledSignalSource(double *numarray,long length,double norm = 1.0);
};

initsave(SampledSignalSource)

class FileSignalSource : public SignalSource {
 public:
    FileSignalSource(char *filename,long bufferlen,long prebuffer,int endianness = -1,double norm = 1.0);
    ~FileSignalSource();
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

class BandIntFilter : public Filter {
 public:
    BandIntFilter(double deltat,double flow,double fhi);
};

class FIRFilter: public Filter {
 public:
    FIRFilter(double *doublearray,int doublenum);
};

class IIRFilter : public Filter {
 public:
    IIRFilter(double *doublearray,int doublenum,double *doublearray,int doublenum);
};

initsave(SignalFilter)

class SignalFilter : public SignalSource {
  public:
    SignalFilter(long length,SignalSource *src,Filter *flt);
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


class DotLagrangeInterpolator : public Interpolator {
 public:
    DotLagrangeInterpolator(int semiwin);
};


class NewLagrangeInterpolator : public Interpolator {
 public:
    NewLagrangeInterpolator(int semiwin);
};

class NoSignal : public Signal {
 public:
    NoSignal();
};

initsave(SumSignal)

class SumSignal : public Signal {
 public:
    SumSignal(Signal *s1,Signal *s2);
};

initsave(InterpolatedSignal)

class InterpolatedSignal : public Signal {
 public:
    InterpolatedSignal(SignalSource *src,Interpolator *interp,
                       double deltat,double prebuffer = 0.0,double norm = 1.0);
    ~InterpolatedSignal();
    
    void setinterp(Interpolator *inte);
};

%pythoncode %{
def getInterpolator(interplen=1):
    if interplen == 0:
        ret = NearestInterpolator()
    elif interplen == -1:
        ret = LinearExtrapolator()
    elif interplen == 1:
        ret = LinearInterpolator()
    elif interplen > 1:
        ret = LagrangeInterpolator(interplen)
    elif interplen < -1:
        ret = DotLagrangeInterpolator(-interplen)
    else:
        raise NotImplementedError, "getInterpolator: undefined interpolator length %s (lisasim-swig.i)." % interplen

    return ret

def getDerivativeInterpolator(interplen=2):
    if interplen > 1:
        return DotLagrangeInterpolator(interplen)
    else: 
        raise NotImplementedError, "getDerivativeInterpolator: undefined interpolator length %s (lisasim-swig.i)." % interplen
%}

%pythoncode %{
def PowerLawNoise(deltat,prebuffer,psd,exponent,interplen=1,seed=0):
    nyquistf = 0.5 / deltat

    if exponent == 0:
        filter = NoFilter()
        normalize = math.sqrt(psd) * math.sqrt(nyquistf)
    elif exponent == 2:
        filter = DiffFilter()
        normalize = math.sqrt(psd) * math.sqrt(nyquistf) / (2.00 * math.pi * deltat)
    elif exponent == -2:
        filter = IntFilter()
        normalize = math.sqrt(psd) * math.sqrt(nyquistf) * (2.00 * math.pi * deltat)
    elif exponent == -4:
        filter = IIRFilter([1,2,1],[0,0.9999*2,-0.9999])
        normalize = math.sqrt(psd) * math.sqrt(nyquistf) * (math.pi * deltat)**2
    else:
        raise NotImplementedError, "PowerLawNoise: undefined PowerLaw exponent %s (lisasim-swig.i)." % exponent

    if seed == 0:
        seed = getcseed()

    whitenoise = WhiteNoiseSource(int(prebuffer/deltat+32),seed)
    filterednoise = SignalFilter(int(prebuffer/deltat+32),whitenoise,filter)

    interp = getInterpolator(interplen)

    noise = InterpolatedSignal(filterednoise,interp,deltat,prebuffer,normalize)

    noise.xmltype = 'PowerLawNoise'
    noise.xmlargs = [deltat,prebuffer,psd,exponent,interplen,seed]

    return noise
%}

%pythoncode %{
def SampledSignal(array,deltat,buffer = 136.0,norm = 1.0,filter = None,interplen = 1,timeoffset = 0.0,endianness = -1,readbuffer=2**20):
    interp = getInterpolator(interplen)

    if isinstance(array,Numeric.ndarray):
        samplednoise = SampledSignalSource(array,norm)
    elif isinstance(array,str):
        readbuffer = max(int(buffer/deltat),int(readbuffer))

        # the endianness parameter supports the MATLAB convention,
        # as well as the XSIL descriptors

        if endianness == 'n':
            endianness = -1
        elif endianness == 'b' or endianness == 'BigEndian':
            endianness = 0
        elif endianness == 'l' or endianness == 'LittleEndian':
            endianness = 1
        
        samplednoise = FileSignalSource(array,readbuffer,int(buffer/deltat),endianness,norm)
    else:
        raise NotImplementedError, "SampledSignal: need Numeric array or filename as first argument (lisasim-swig.i)."

    if not filter:
        filteredsamples = 0
        interpolatednoise = InterpolatedSignal(samplednoise,interp,deltat,-timeoffset)
    else:
        filteredsamples = SignalFilter(int(buffer/deltat),samplednoise,filter)
        interpolatednoise = InterpolatedSignal(filteredsamples,interp,deltat,-timeoffset)

    return interpolatednoise
%}

%feature("docstring") CachedSignal "
CachedSignal(Signal,bufferlen,deltat,interplen = 4)
"

initsave(CachedSignal)

class CachedSignal : public Signal {
 public:
    CachedSignal(Signal *s,long length,double deltat,int interplen = 4);
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

initsave(SimpleBinary)

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

    wave = NoiseWave(hparray,hcarray,len,deltat,prebuffer,norm,filter,interp,elat,elon,pol)

    wave.xmltype = 'SampledWave'
    wave.xmlargs = [hparray,hcarray,len,deltat,prebuffer,norm,filter,interp,elat,elon,pol]

    return wave
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

exceptionhandle(WaveArray::WaveArray,ExceptionOutOfBounds,PyExc_ValueError)

class WaveArray : public WaveObject {
 public:
    WaveArray(Wave **WaveSeq, int WaveNum);
    ~WaveArray();

    Wave *firstwave();
    Wave *nextwave();
};


/* -------- TDI objects -------- */

exceptionhandle(TDI::alpham,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::betam,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::gammam,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::zetam,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::alpha1,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::alpha2,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::alpha3,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::zeta1,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::zeta2,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::zeta3,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::P,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::E,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::U,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::Xm,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::Ym,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::Zm,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::Xmlock1,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::Xmlock2,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::Xmlock3,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::X1,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::X2,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::X3,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::y,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::z,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::y123,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::y231,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::y312,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::y321,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::y132,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::y213,ExceptionOutOfBounds,PyExc_IndexError)

exceptionhandle(TDI::z123,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::z231,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::z312,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::z321,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::z132,ExceptionOutOfBounds,PyExc_IndexError)
exceptionhandle(TDI::z213,ExceptionOutOfBounds,PyExc_IndexError)

%nodefault TDIobject;
class TDIobject : public Signal {
 public:
    TDIobject(TDI *t);
    virtual ~TDIobject();
    
    virtual double value(double t) = 0;
};

// %nodefault TDIobjectpnt;
class TDIobjectpnt  : public TDIobject {
 public:
    TDIobjectpnt(TDI *t,double (TDI::*o)(double t));
    ~TDIobjectpnt();

    double value(double t);
};

%nodefault timeobject;
class timeobject : public Signal {};

exceptionhandle(fastgetobsc,ExceptionKeyboardInterrupt,PyExc_KeyboardInterrupt)

extern void fastgetobs(double *numarray,long length,long samples,double stime,Signal **thesignals,int signals,double inittime);
extern void fastgetobsc(double *numarray,long length,long samples,double stime,Signal **thesignals,int signals,double inittime);

%newobject TDI::alpham();
%newobject TDI::betam();
%newobject TDI::gammam();

%newobject TDI::zetam();

%newobject TDI::alpha1();
%newobject TDI::alpha2();
%newobject TDI::alpha3();

%newobject TDI::zeta1();
%newobject TDI::zeta2();
%newobject TDI::zeta3();

%newobject TDI::P();
%newobject TDI::E();
%newobject TDI::U();

%newobject TDI::Xm();
%newobject TDI::Ym();
%newobject TDI::Zm();

%newobject TDI::Xmlock1();
%newobject TDI::Xmlock2();
%newobject TDI::Xmlock3();

%newobject TDI::X1();
%newobject TDI::X2();
%newobject TDI::X3();

%newobject TDI::y123();
%newobject TDI::y231();
%newobject TDI::y312();

%newobject TDI::y213();
%newobject TDI::y321();
%newobject TDI::y132();

%newobject TDI::z123();
%newobject TDI::z231();
%newobject TDI::z312();

%newobject TDI::z213();
%newobject TDI::z321();
%newobject TDI::z132();

class TDI {
 public:
    TDI() {};
    virtual ~TDI() {};

    virtual void reset() {};
    
    virtual double alpham(double t);
    TDIobject *alpham();
    virtual double betam(double t);
    TDIobject *betam();
    virtual double gammam(double t);
    TDIobject *gammam();

    virtual double zetam(double t);
    TDIobject *zetam();

    virtual double alpha1(double t);
    TDIobject *alpha1();
    virtual double alpha2(double t);
    TDIobject *alpha2();
    virtual double alpha3(double t);
    TDIobject *alpha3();

    virtual double zeta1(double t);
    TDIobject *zeta1();
    virtual double zeta2(double t);
    TDIobject *zeta2();
    virtual double zeta3(double t);
    TDIobject *zeta3();

    // P, E, U still have non-signed delays
    
    virtual double P(double t);
    TDIobject *P();
    virtual double E(double t);
    TDIobject *E();
    virtual double U(double t);
    TDIobject *U();
    
    virtual double Xm(double t);
    TDIobject *Xm();
    virtual double Ym(double t);
    TDIobject *Ym();
    virtual double Zm(double t);
    TDIobject *Zm();

    virtual double Xmlock1(double t);
    TDIobject *Xmlock1();
    virtual double Xmlock2(double t);
    TDIobject *Xmlock2();
    virtual double Xmlock3(double t);
    TDIobject *Xmlock3();
    
    virtual double X1(double t);
    TDIobject *X1();
    virtual double X2(double t);
    TDIobject *X2();
    virtual double X3(double t);
    TDIobject *X3();

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);

    virtual double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);
    virtual double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t);

    double y123(double t);
    TDIobject *y123();
    double y231(double t);
    TDIobject *y231();
    double y312(double t);
    TDIobject *y312();

    double y321(double t);
    TDIobject *y321();
    double y132(double t);
    TDIobject *y132();
    double y213(double t);
    TDIobject *y213();
    
    double z123(double t);
    TDIobject *z123();
    double z231(double t);
    TDIobject *z231();
    double z312(double t);
    TDIobject *z312();
    
    double z321(double t);
    TDIobject *z321();
    double z132(double t);
    TDIobject *z132();
    double z213(double t);
    TDIobject *z213();
    
    double time(double t);
    timeobject *time();
    double t(double t);
    timeobject *t();
};

initsave(SampledTDI)

%apply Noise *PYTHON_SEQUENCE_NOISE[ANY] {Noise *yijk[6], Noise *zijk[6]}

class SampledTDI : public TDI {
 public:
    SampledTDI(LISA *lisa,Noise *yijk[6],Noise *zijk[6]);
    ~SampledTDI();
};

class SampledTDIaccurate : public SampledTDI {
 public:
    SampledTDIaccurate(LISA *lisa,Noise *yijk[6],Noise *zijk[6]);
    ~SampledTDIaccurate();
};

/* We're holding on to the constructor args so that the LISA object
   won't get destroyed if they fall out of scope: we may still need
   them! */

initsave(TDIquantize)

class TDIquantize : public TDI {
 public:
    TDIquantize(TDI *bt,double qlev,int qbits,int qsat);
    ~TDIquantize();
};

%feature("docstring") TDInoise "
TDInoise(lisa,PMnoise[6],      SHnoise[6],      LSnoise[6])
TDInoise(lisa,PMdt,   PMpsd,   SHdt,   SHpsd,   LSdt,   LSpsd)

all return TDI objects that implement the inter- and intra-spacecraft
phase measurements using the standard noise transfer functions for
proof-mass noise, shot noise, and laser noise. All the TDI observables
defined in the base TDI class are available from TDInoise.

- In the first form of the constructor, Synthetic LISA Noise/Signal
  objects (e.g., as created with PowerLawNoise or SampledSignal) must be
  provided for all 18 fundamental noises, in the orders [according to
  the PRD 71, 022001 (2005) naming]

  [pn_1,pn_1*,pn_2,pn_2*,pn_3,pn_3*]
  [y^sh_{12},y^sh_{21},y^sh_{23},y^sh_{32},y^sh_{31},y^sh_{13}]
  [C_1,C_1*,C_2,C_2*,C_3,C_3*]
  
- In the second form of the constructor, the noise object are created
  internally as pseudorandom noise objects (equivalent to PowerLawNoise)
  with sampling times PMdt, SHdt, and LSdt, and with approximate
  one-sided power spectral densities

  PMpsd*(f/Hz)^-2 Hz^-1
  SHpsd*(f/Hz)^2  Hz^-1
  LSpsd           Hz^-1

Note: resetting the TDInoise object will reset all the component noise
objects. 

TODO: might want to allow different interpolation widths for the
self-built pseudorandom noise objects."

initdoc(TDInoise)

initsave(TDInoise)

// not so sure about this...

%feature("addtofunc") TDInoise::setphlisa {
        self.phlisa = args
}

// utility Python functions needed by TDInoise constructor

%pythoncode %{
    def lighttime(lisa):
        # to estimate size of noisebuffer, take maximum armlength at time zero,
        # and add 10% for uplink-downlink uncertainty and flexing

        return 1.10 * max(lisa.armlength(1,0.0),
                          lisa.armlength(2,0.0),
                          lisa.armlength(3,0.0))

    def stdproofnoise(lisa,stproof,sdproof,interp=1,seed=0):
        # we need quadruple retardations for the V's appearing in the z's
        # (octuple for 2nd-gen TDI); we add two sampling times to allow linear
        # interpolation for large sampling times
    
        pbtproof = 8.0 * lighttime(lisa) + 2.0*stproof
    
        return PowerLawNoise(stproof,pbtproof,sdproof,-2.0,interp,seed)

    def stdproofnoisepink(lisa,stproof,sdproof,sf0proof,interp=1,seed=0):
        # we need quadruple retardations for the V's appearing in the z's
        # (octuple for 2nd-gen TDI); we add two sampling times to allow linear
        # interpolation for large sampling times
    
        pbtproof = 8.0 * lighttime(lisa) + 2.0*stproof
    
        return PowerLawNoise(stproof,pbtproof,sdproof * sf0proof**2,-4.0,interp,seed)

    def stdopticalnoise(lisa,stshot,sdshot,interp=1,seed=0):
        # we need only triple retardations for the shot's appearing in the y's
        # (septuple for 2nd-gen TDI); we add two sampling times to allow linear
        # interpolation for large sampling times
    
        pbtshot = 7.0 * lighttime(lisa) + 2.0*stshot
    
        return PowerLawNoise(stshot,pbtshot,sdshot,2.0,interp,seed)

    def stdlasernoise(lisa,stlaser,sdlaser,interp=1,seed=0):
        pbtlaser = 8.0 * lighttime(lisa) + 2.0*stlaser

        return PowerLawNoise(stlaser,pbtlaser,sdlaser,0.0,interp,seed)
%}

%feature("pythonprepend") TDInoise::TDInoise %{
        self.lisa = args[0]
        args = args[1:]
        
        # if no parameters are passed, used default value
        # if only one parameter is passed, use it as stime
        
        if len(args) == 0:
            args = (1.0, 2.5e-48, 1.0, 1.8e-37, 1.0, 1.1e-26)
        elif len(args) == 1:
            args = (args[0], 2.5e-48, args[0], 1.8e-37, args[0], 1.1e-26)
        
        # proof-mass: use deltat-psd model if we find a number, or assume
        # six Noise objects otherwise
        
        if type(args[0]) in (int,float):
            self.pm = [stdproofnoise(self.lisa,args[0],args[1]) for i in range(6)]
            args = args[2:]
        else:
            self.pm = args[0]
            args = args[1:]
        
        # shot noise: same story
        
        if type(args[0]) in (int,float):
            self.pd = [stdopticalnoise(self.lisa,args[0],args[1]) for i in range(6)]
            args = args[2:]
        else:
            self.pd = args[0]
            args = args[1:]
        
        # laser noise may not be passed: if so, set it to zero
        
        if len(args) > 0:
            if type(args[0]) in (int,float):
                self.c = [stdlasernoise(self.lisa,args[0],args[1]) for i in range(6)]
                args = args[2:]
            else:
                self.c = args[0]
                args = args[1:]    
        else:
            self.c = [NoSignal() for i in range(6)]
        
        args = (self.lisa,self.pm,self.pd,self.c)
%}

%apply Noise *PYTHON_SEQUENCE_NOISE[ANY] {Noise *proofnoise[6], Noise *shotnoise[6], Noise *lasernoise[6]}

%apply double PYTHON_SEQUENCE_DOUBLE[ANY] {double stproof[6], double sdproof[6], double stshot[6], double sdshot[6], double stlaser[6], double sdlaser[6], double claser[6]}

class TDInoise : public TDI {
 public:
    TDInoise(LISA *mylisa, Noise *proofnoise[6], Noise *shotnoise[6], Noise *lasernoise[6]);

    virtual ~TDInoise();

    void setphlisa(LISA *mylisa);

    void lock(int master);
    
    void reset(unsigned long seed = 0);
};        

/* We're holding on to the constructor args so that the LISA/Noise
   objects won't get destroyed if they fall out of scope: we may still
   need them! */

initsave(TDIaccurate)

class TDIaccurate : public TDInoise {
 public:
    TDIaccurate(LISA *mylisa, Noise *proofnoise[6],Noise *shotnoise[6],Noise *lasernoise[6]);
    ~TDIaccurate();
};

/* We're holding on to the constructor args so that the LISA/Wave
   objects won't get destroyed if they fall out of scope: we may still
   need them for TDInoise! */

initsave(TDIdoppler)

class TDIdoppler : public TDInoise {
 public:
    TDIdoppler(LISA *mylisa, Noise *proofnoise[6],Noise *shotnoise[6],Noise *lasernoise[6]);
    ~TDIdoppler();
};

/* We're holding on to the constructor args so that the LISA/Wave
   objects won't get destroyed if they fall out of scope: we may still
   need them for TDInoise! */

%apply double PYTHON_SEQUENCE_DOUBLE[ANY] {double laserfreqs[6]}

initsave(TDIcarrier)

class TDIcarrier : public TDInoise {
 public:
    TDIcarrier(LISA *mylisa,double laserfreqs[6]);
    ~TDIcarrier();

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

    double Phi(int slink,double t);
};
