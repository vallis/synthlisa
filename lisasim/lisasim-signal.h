/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */
 
#ifndef _LISASIM_SIGNAL_H_
#define _LISASIM_SIGNAL_H_

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

class RingBuffer {
 private:
    double *data;
    long length;

 public:
	RingBuffer(long len);
	~RingBuffer();
	
	void reset();
	
	inline double& operator[](long pos);
};

inline double& RingBuffer::operator[](long pos) {
	if(pos < 0) pos += length * (abs(pos)/length + 1);

	return data[pos % length];
}


/* SignalSource interface: can return [pos], may be reset to beginning,
   will throw exception (ExceptionOutOfBounds) at EOF or for stale 
   access */

class SignalSource {
 public:
	virtual ~SignalSource() {};

	virtual void reset(unsigned long seed = 0) {};
	virtual double operator[](long pos) = 0;
};


/* BufferedSignalSource interface: inherit from it, redefine
   getvalue(pos), which will be called only once to get each new sample */

class BufferedSignalSource : public SignalSource {
 private:
	RingBuffer buffer;
	long length;

	long current;

 public:
	BufferedSignalSource(long len);
	virtual ~BufferedSignalSource() {}; // ??? would "= 0" do here?

	virtual double getvalue(long pos) = 0;

	virtual void reset(unsigned long seed = 0); // ??? redefining default
	virtual double operator[](long pos);
};


/* WhiteNoiseSource uses the GSL random number routines; perhaps there
   is something more efficient. */

#include "GSL/gsl_rng.h"

class WhiteNoiseSource : public BufferedSignalSource {
 private:	
    gsl_rng *randgen;

    int cacheset;
    double cacherand;

	double normalize;

    void seedrandgen(unsigned long seed);
    double deviate();

    static unsigned long globalseed;

 public:
	WhiteNoiseSource(long len,unsigned long seed = 0,double norm = 1.0);
	~WhiteNoiseSource();	

	double getvalue(long pos);
		
	void reset(unsigned long seed = 0);  // ??? redefining default

    static void setglobalseed(unsigned long seed = 0);
    static unsigned long getglobalseed();
};


class SampledSignalSource : public SignalSource {
 private:
	double *data;
	long length, warn;

	double normalize;

 public:
	SampledSignalSource(double *darray,long len,double norm = 1.0);

	double operator[](long pos);
};


class FileSignalSource : public BufferedSignalSource {
 private:
    FILE *file;
    int convert;
    
    double *filebuffer;
    long length, initpos;

    double normalize;

    void loadbuffer();
 
 public:
    FileSignalSource(char *filename,long bufferlen,long prebuffer,int endianness = -1,double norm = 1.0);
    ~FileSignalSource();
    
	double getvalue(long pos);
	
	void reset(unsigned long seed = 0);  // ??? redefining default
};


/* Filters! It might be tempting to inline the getvalue() call, but it
   would probably do little good since they are called as virtual methods
   from a base pointer. */

class Filter {
 public:
    virtual ~Filter() {};
 
    virtual double getvalue(SignalSource &x,SignalSource &y,long pos) = 0;
};


class NoFilter : public Filter {
 public:
    double getvalue(SignalSource &x,SignalSource &y,long pos);
};


class IntFilter : public Filter {
 private:
    double alpha;

 public:
    IntFilter(double a = 0.9999);

    double getvalue(SignalSource &x,SignalSource &y,long pos);
};


class DiffFilter : public Filter {
 public:
    double getvalue(SignalSource &x,SignalSource &y,long pos);
};


class BandIntFilter : public Filter {
 private:
    double alpha0, alpha1, beta1;
    
 public:
    BandIntFilter(double deltat,double flow,double fhi);

    double getvalue(SignalSource &x,SignalSource &y,long pos);
};


class FIRFilter : public Filter {
 private:
 	double *a;
 	int length;

 public:
 	FIRFilter(double *aarray,int length);
 	~FIRFilter();

    double getvalue(SignalSource &x,SignalSource &y,long pos);
};


class IIRFilter : public Filter {
 private:
 	double *a, *b;
 	int lengtha, lengthb;

 public:
 	IIRFilter(double *aarray,int lena,double *barray,int lenb);
 	~IIRFilter();
 
    double getvalue(SignalSource &x,SignalSource &y,long pos);
};


class SignalFilter : public BufferedSignalSource {
 private:
	SignalSource *source;
	Filter *filter;

 public:
	SignalFilter(long length,SignalSource *src,Filter *flt);

	double getvalue(long pos);
	
	void reset(unsigned long seed = 0);  // ??? redefining default
};


/* Interface for Signal: value(time) and value(timebase,timecorr). Also
   reset(). */

class Signal {
 public:
	virtual ~Signal() {};

    virtual void reset(unsigned long seed = 0) {};  // ??? redefining default

	virtual double value(double time) = 0;

	// standard implementation if we are not more precise than this
	virtual double value(double timebase,double timecorr) {
		return value(timebase + timecorr);
	};

	// for backward compatibility

	virtual double operator[](double time) { return value(time); };
	virtual double operator()(double timebase,double timecorr) { return value(timebase,timecorr); }
	
    virtual double noise(double time) { return value(time); };    
	virtual double noise(double timebase,double timecorr) { return value(timebase,timecorr); }
};


// With this typedef, Signal = Noise for backward compatibility

typedef Signal Noise;


// Interpolators!

class Interpolator {
 public:
	virtual ~Interpolator() {};

    virtual double getvalue(SignalSource &y,long ind,double dind) = 0;
};


class NearestInterpolator : public Interpolator {
 public:
	double getvalue(SignalSource &y,long ind,double dind);
};


class LinearInterpolator : public Interpolator {
 public:
	double getvalue(SignalSource &y,long ind,double dind);
};


class LinearExtrapolator : public Interpolator {
 public:
	double getvalue(SignalSource &y,long ind,double dind);
};


class LagrangeInterpolator : public Interpolator {
 private:
    int window, semiwindow;

    double *xa,*ya;
    double *c,*d;

    double polint(double x);

 public:
    LagrangeInterpolator(int semiwin);
    virtual ~LagrangeInterpolator();

    double getvalue(SignalSource &y,long ind,double dind);
};

class DotLagrangeInterpolator : public Interpolator {
 private:
    int window, semiwindow;

    double *xa,*ya,*yd;

    double dpolint(double x);

 public:
    DotLagrangeInterpolator(int semiwin);
    virtual ~DotLagrangeInterpolator();

    double getvalue(SignalSource &y,long ind,double dind);
};

class NewLagrangeInterpolator : public Interpolator {
 private:
    int window;
    double semiwindow;

    double *xa;
    double *ya;
    
    double *c,*d;

    double polint(double x);

 public:
    NewLagrangeInterpolator(int semiwin);
    virtual ~NewLagrangeInterpolator();

	double getvalue(SignalSource &y,long ind,double dind);
};


// get one of the above by choosing its length (-1 for Extrapolator)

Interpolator *getInterpolator(int interplen);
Interpolator *getDerivativeInterpolator(int interplen);

// --- InterpolatedSignal ---

class NoSignal : public Signal {
 public:
    NoSignal() {};
    
    void reset(unsigned long seed = 0) {};
    
    double value(double time) { return 0.0; };
    double value(double timebase,double timecorr) { return 0.0; }
};


class SumSignal : public Signal {
 private:
    Signal *signal1, *signal2;

 public:
    SumSignal(Signal *s1,Signal *s2) : signal1(s1), signal2(s2) {};
    
    void reset(unsigned long seed = 0) {
        signal1->reset();
        signal2->reset();        
    };

    double value(double time) {
        return signal1->value(time) + signal2->value(time);
    };

    double value(double timebase,double timecorr) {
        return signal1->value(timebase,timecorr) +
               signal2->value(timebase,timecorr);
    };
};


class InterpolatedSignal : public Signal {
 private:
	SignalSource *source;
	Interpolator *interp;
	
	double samplingtime, prebuffertime, normalize;
	
 public:
	InterpolatedSignal(SignalSource *src,Interpolator *inte,
					   double deltat,double prebuffer = 0.0,double norm = 1.0);
    ~InterpolatedSignal() {};

    void reset(unsigned long seed = 0);  // ??? redefining default

	double value(double time);
	double value(double timebase,double timecorr);
	
	void setinterp(Interpolator *inte);
};


/* ??? Since PowerLawNoise and SampledSignal really exist as separate
   classes to garbage collect gracefully their components (and to be
   convenient), it might be wise to implement them directly in Python, or
   also in Python). This would save one indirection (currently inlined,
   although the inline is probably not realized by the compiler because
   the classes are virtual. */

// --- PowerLawNoise ---

class PowerLawNoise : public Signal {
 private:
	WhiteNoiseSource *whitenoise;
	
	Filter *filter;
	SignalFilter *filterednoise;

	Interpolator *interp;
	InterpolatedSignal *interpolatednoise;

 public:
	PowerLawNoise(double deltat,double prebuffer,
		double psd,double exponent,int interplen = 1,unsigned long seed = 0);
	~PowerLawNoise();

	void reset(unsigned long seed = 0);  // ??? redefining default

	double value(double time);
	double value(double timebase,double timecorr);
};

inline double PowerLawNoise::value(double time) {
	return interpolatednoise->value(time);
}

inline double PowerLawNoise::value(double timebase,double timecorr) {
	return interpolatednoise->value(timebase,timecorr);
}


// --- SampledSignal ---

class SampledSignal : public Signal {
 private:
	SampledSignalSource *samplednoise;

	SignalFilter *filteredsamples;

	Interpolator *interp;
	InterpolatedSignal *interpolatednoise;

 public:
	SampledSignal(double *narray,long length,double deltat,double prebuffer,
		double norm = 1.0,Filter *filter = 0,int interplen = 1);
	~SampledSignal();

    // nothing to reset...

	double value(double time);
	double value(double timebase,double timecorr);
};

inline double SampledSignal::value(double time) {
	return interpolatednoise->value(time);
}

inline double SampledSignal::value(double timebase,double timecorr) {
	return interpolatednoise->value(timebase,timecorr);
}


// --- CachedSignal (uses ResampledSignalSource) ---

class ResampledSignalSource : public BufferedSignalSource {
 private:
	double deltat, prebuffer;

	Signal *signal;

 public:
	ResampledSignalSource(long len,double dt,double pbt,Signal *s)
		: BufferedSignalSource(len), deltat(dt), prebuffer(pbt), signal(s) {};

	void reset(unsigned long seed = 0);  // ??? redefining default

	double getvalue(long pos);
};

class CachedSignal : public Signal {
 private:
	ResampledSignalSource *resample;
	Interpolator *interp;
	InterpolatedSignal *interpsignal;

 public:
    CachedSignal(Signal *s,long length,double deltat,int interplen = 4);
	~CachedSignal();

    void reset(unsigned long seed = 0);  // ??? redefining default

	double value(double time);
	double value(double timebase,double timecorr);
};

inline double CachedSignal::value(double time) {
	return interpsignal->value(time);
}

inline double CachedSignal::value(double timebase,double timecorr) {
	return interpsignal->value(timebase,timecorr);
}

#endif /* _LISASIM_SIGNAL_H_ */











