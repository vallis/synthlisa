// FIX: the following have all rather specific uses, and may be
//      replaced by the suitable use of PyLISA. They are moved to this
//      private C++ and interface file.

// NOTE: this interface file is incomplete and won't compile without changes

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

class Filter;

class FIR : public Filter {
 public:
 	FIR(double *doublearray, int doublenum); 
 	~FIR();
};

class IIR : public Filter {
 public:
 	IIR(double *doublearray, int doublenum, double *doublearray, int doublenum);
    ~IIR();
};

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
InterpolateNoise(deltat,prebuffer,psd,filter,interp=1) and
InterpolateNoise(array,length,deltat,prebuffer,norm,filter=0,interp=1)
return a Noise object from either filtered and interpolated
pseudorandom whitenoise (first form) or from a filtered and
interpolated array (second form).

For the first form:

- deltat is the time spacing (seconds) of the pseudorandom samples at
  generation;

- prebuffer (seconds??) sets the minimum interval of buffering of the
  pseudorandom sequences (after the noise has been requested at time
  t, the earliest noise value guaranteed to be available will be at
  time (t-prebuf);

- psd and filter set the one-sided PSD of the filtered pseudorandom noise
  to psd*(f/Hz)^filter Hz^-1 (currently exp = -2, 0, 2 are implemented);

- wind (> 1) sets the semiwidth of the data window used in Lagrange
  interpolation (1 yields linear interpolation).

For the second form:

...
"

initdoc(InterpolateNoise)

initsave(InterpolateNoise)

// DEV: Should add machinery to let the typemaps infer the length of
// the array and pass it on

%apply double *NUMPY_ARRAY_DOUBLE { double *noisebuf };

class InterpolateNoise : public Noise {
 public:
    InterpolateNoise(double sampletime,double prebuffer,double density,double exponent,int swindow = 1);

    InterpolateNoise(double *noisebuf,long samples,double sampletime,double prebuffer,double density, double exponent, int swindow = 1);

    virtual ~InterpolateNoise();
    
    void reset();
    
    double noise(double time);
    double noise(double timebase,double timecorr);

    void setinterp(int window);
};


class InterpolateNoise2 : public Noise {
 public:
	InterpolateNoise2(double deltat,double prebuffer,double psd,double filterexp,int interplen = 1);
	InterpolateNoise2(double deltat,double prebuffer,double psd,Filter *myfilter,int interplen = 1);

	InterpolateNoise2(double *numarray,long length,double deltat,double prebuffer,double norm,double filterexp = 0.0,int interplen = 1);
	InterpolateNoise2(double *numarray,long length,double deltat,double prebuffer,double norm,Filter *myfilter,int interplen = 1);

	virtual ~InterpolateNoise2();

    void reset();

    double noise(double time);
    double noise(double timebase,double timecorr);
};


