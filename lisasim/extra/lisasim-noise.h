/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */
 
#ifndef _LISASIM_NOISE_H_
#define _LISASIM_NOISE_H_

#include <iostream>
using namespace std;

#include <limits.h>
#include "GSL/gsl_rng.h"

// the class structure used by the noise classes is somewhat complex
// but it is very efficient and very extendable

// here we go:

// at the top sits InterpolateNoise, which is the main user interface
// depending on how it is called, it will select different origins for
// the noise (pseudorandom time series, memory buffers, files), and
// different interpolation schemes

// the interpolation scheme is a function object derived from the base
// class "Interpolator"

// the noise (before interpolation, at a constant cadence) is made in the object
// MakeNoise, which applies a noise filter to the time series returned
// by a noise getter; MakeNoise uses to internal ring buffers, implemented
// by the class BufferNoise

// the noise getters are derived from GetNoise;
// currently two getters are implemented, WhiteNoise for pseudorandom noise
// and SampledNoise for noise read from memory

//        BufferNoise  GetNoise  Filter
//                  \      |     /
//                   \     |    /
//                    \    |   /  
//                     MakeNoise  Interpolator
//                            |        |
//                         InterpolateNoise


// ??? Ah, but the negative indexes should also wrap around, otherwise
//     the transient of IIR filters will be wrong!

class BufferNoise {
 private:
    double *data;
    long maxbuffer;

	double zero;

 public:
    BufferNoise(long l) : maxbuffer(l), zero(0.0) {
		if ( !(data = new double[maxbuffer]) ) {
			cout << "RingBufferNoise::RingBufferNoise: out of memory" << endl;
			abort();
		}

		for(int i=0;i<maxbuffer;i++)
			data[i] = 0.0;
	};
	
	~BufferNoise() {
		delete [] data;
	};
	
	void reset() {
		for(int i=0;i<maxbuffer;i++)
			data[i] = 0.0;
	};
	
	double& operator[](long pos) {
		if(pos < 0) {
			if(pos <= -maxbuffer)
				return zero;
			else
				return data[(pos + maxbuffer) % maxbuffer];
		} else {
			return data[pos % maxbuffer];
		}
    };
};

class GetNoise {
 public:
    virtual void reset() {};
    virtual void operator()(BufferNoise &x,long pos) = 0;
};

class WhiteNoise : public GetNoise {
 private:
    gsl_rng *randgen;

    int cacheset;
    double cacherand;

    void seedrandgen(unsigned long seed);
    double deviate();

 public:
    WhiteNoise(unsigned long seed = 0);
    virtual ~WhiteNoise();

    void reset(unsigned long seed = 0) {
		seedrandgen(seed);
    };

    void operator()(BufferNoise &x, long pos) {
		x[pos] = deviate();
    };
};

// this version of SampledNoise gets a buffer allocated externally
// and works from there

class SampledNoise : public GetNoise {
 private:
    double *data;
    long maxbuffer;

 public:
    SampledNoise(double *d,long m) : data(d), maxbuffer(m) {};

    void operator()(BufferNoise &x, long pos) {
		if (pos < 0 || pos > maxbuffer) {
			cout << "SampledNoise::(): index (" << pos << ") out of range" << endl;
			abort();
		} else {
			x[pos] = data[pos];
		}
    }
};

// Filters... should not know about the length of the buffer
// this means that they should be passed intelligent objects that can
// do their own indexing...

class Filter {
 public:
    // need a virtual constructor, too?
    virtual ~Filter() {};
 
 	virtual int length() { return 0; };
 
    virtual void operator()(BufferNoise &x,BufferNoise &y,long pos) = 0;
};

class NoFilter : public Filter {
 public:
    void operator()(BufferNoise &x,BufferNoise &y,long pos) {
		y[pos] = x[pos];
    };
};

// these filters should be causal; they pose a problem only if the
// length of the filter is longer than the bufferlength

class IntFilter : public Filter {
 private:
    double alpha;

 public:
    IntFilter(double a = 0.9999) : alpha(a) {};

	int length() {
		return 2;
	};

    void operator()(BufferNoise &x,BufferNoise &y,long pos) {
		y[pos] = alpha * y[pos - 1] + x[pos];
    };
};

class DiffFilter : public Filter {
 public:
 	int length() {
 		return 2;
 	};
 
    void operator()(BufferNoise &x,BufferNoise &y,long pos) {
		y[pos] = x[pos] - x[pos-1];
    };
};

class FIR : public Filter {
 private:
 	double *a;
 	int len;

 public:
 	FIR(double *aarray,int length) : len(length) {
 		a = new double[len];
 		
 		for(int i=0;i<len;i++) a[i] = aarray[i];
 	};
 
 	~FIR() {
 		delete [] a;
 	};
 
 	int length() {
 		return len;
 	};
 
 	void operator()(BufferNoise &x,BufferNoise &y,long pos) {
		double acc = 0.0;
		
		for(int i=0;i<len;i++) acc += a[i] * x[pos-i];

		y[pos] = acc;
 	};
};

class IIR : public Filter {
 private:
 	double *a, *b;
 	int lena, lenb;

 public:
 	IIR(double *aarray,int lengtha,double *barray,int lengthb)
 		: lena(lengtha), lenb(lengthb) {
 		a = new double[lena];
 		b = new double[lenb];
 		
 		for(int i=0;i<lena;i++) a[i] = aarray[i];
 		for(int j=0;j<lenb;j++) b[j] = barray[j];
 	};

 	~IIR() {
 		delete [] b;
 		delete [] a;
 	};
 
 	int length() {
		return lena > lenb ? lena : lenb;
 	};
 
 	void operator()(BufferNoise &x,BufferNoise &y,long pos) {
		double acc = 0.0;
		
		for(int i=0;i<lena;i++) acc += a[i] * x[pos-i];
		for(int j=1;j<lenb;j++) acc += b[j] * y[pos-j];

		y[pos] = acc;
 	};
};

// the main fixed-timestep noise class

class MakeNoise {
 public:
    virtual void reset() = 0;

    virtual double operator[](long pos) = 0;
};

class FilterMakeNoise : public MakeNoise {
 private:
    BufferNoise *x, *y;
    long maxbuffer;
    long earliest, latest;
	long filterlength;

    GetNoise *get;
    Filter *filter;

 public:
    FilterMakeNoise(GetNoise *g,Filter *f,long l)
	: maxbuffer(l), latest(-1), get(g), filter(f) {
		filterlength = f->length();

		if(filterlength + 1 >= maxbuffer) {
			cout << "FilterMakeNoise::FilterMakeNoise(): the filter is longer "
				 << "than the Ring Buffer!"
				 << " [" << __FILE__ << ":" << __LINE__ << "]" << endl;

			abort();
		}

		earliest = -filterlength;

		x = new BufferNoise(maxbuffer);
		y = new BufferNoise(maxbuffer);
    };

    virtual ~FilterMakeNoise() {
		delete y;
		delete x;
    }

    void reset() {
		get->reset();

		x->reset();
		y->reset();

		earliest = -filterlength;
		latest = -1;
    };

	// it would be better to throw an exception!

    double operator[](long pos) {
		if(pos < earliest) {
			cout << "FilterMakeNoise::[]: trying to access noise element ("
				 << pos << ") before oldest kept"
				 << " [" << __FILE__ << ":" << __LINE__ << "]" << endl;
			
			abort();
		} else if (pos > latest) {
			for(int i=latest+1;i<=pos;i++) {
				(*get)(*x,i);
				(*filter)(*x,*y,i);
			}
	
			latest = pos;
			
			if(latest - earliest >= maxbuffer)
				earliest = latest - maxbuffer + 1;
	
			return (*y)[pos];
		} else {
			return (*y)[pos];
		}
    };
};


// interpolators

class Interpolator {
 public:
  	virtual int length() { return 1; };
 
    virtual double operator()(MakeNoise &y,long ind,double dind) = 0;
};

class NearestInterpolator : public Interpolator {
 public:
  	int length() {
  		return 1;
  	};

    double operator()(MakeNoise &y,long ind,double dind) {
		return (dind < 0.5 ? y[ind] : y[ind+1]);
    };
};

// 0 < dind < 1; the desired sample would be at ind + dind

class LinearInterpolator : public Interpolator {
 public:
  	int length() {
  		return 2;
  	};

    double operator()(MakeNoise &y,long ind,double dind) {
		return (1.0 - dind) * y[ind] + dind * y[ind+1];
    };
};

// this is to use only "old" values, with (implicitly) 1 < dind < 2

class LinearExtrapolator : public Interpolator {
 public:
  	int length() {
  		return 2;
  	};

    double operator()(MakeNoise &y,long ind,double dind) {
		return (-dind) * y[ind-1] + (1.0 + dind) * y[ind];
    };
};

class LagrangeInterpolator : public Interpolator {
 private:
    int window, semiwindow;

    double *xa,*ya;
    double *c,*d;

    double polint(double x);

 public:
    LagrangeInterpolator(int sw);

    virtual ~LagrangeInterpolator();

  	int length() {
  		return window;
  	};

    double operator()(MakeNoise &y,long ind,double dind) {
 		for(int i=0;i<semiwindow;i++) {
	    	ya[semiwindow-i] = y[ind-i];
	    	ya[semiwindow+i+1] = y[ind+i+1];
		}

		return polint(semiwindow+dind);
    };
};

// we leave "Noise" as a generic interface which defines only [] and "noise"
// [] is a pure virtual function, hence it MUST be implemented by derived classes

class Noise {
 public:
    virtual void reset() {};

	virtual double operator[](double time) = 0;

    virtual double noise(double time) {
		return (*this)[time];
    };
    
    virtual double noise(double timebase,double timecorr) {
		return noise(timebase + timecorr);
    };
};

// dimensioned interpolated noise

class InterpolateNoise : public Noise {
 private:
    double samplingtime, nyquistf;

    double prebuffertime;
    double maxtime;

    double timewindow;
    double lasttime;
  
    GetNoise *getnoise;

    Filter *filter;

    Interpolator *interp;

    void setfilter(double ex);
    void setnorm(double sd, double ex);
    void setnormsampled(double sd, double ex);

 protected:
    MakeNoise *thenoise;
    double normalize;

 public:
    InterpolateNoise(double sampletime,double prebuffer,double density,double exponent,int swindow = 1);

    InterpolateNoise(double *noisebuf,long samples,double sampletime,double prebuffer,double density,double exponent = 0.0,int swindow = 1);

    virtual ~InterpolateNoise();
    
    void reset();

    double operator[](double time);

    double noise(double time) {
		return (*this)[time];
    };

    double noise(double timebase,double timecorr);

    // expose the setinterp method so it becomes possible to change
    // the interpolation window dynamically for a noise object
    void setinterp(int window);
};


class InterpolateNoise2 : public Noise {
 private:
    double samplingtime, nyquistf;
    double prebuffertime;

    Filter *filter; int ownfilter;
    Interpolator *interp; int owninterp;

    double normalize;
  
    GetNoise *getnoise;
    MakeNoise *thenoise;

    double maxtime, lasttime, timewindow;

	long efflength;
	

	Filter *getfilter(double filterexp);
	Interpolator *getinterp(int interplen);

	double getnormrandom(double psd,double filterexp);
	double getnormsampled(double norm,double filterexp);

	void settime();

	void initrandom();
	void initsampled(double *narray,long length);

 public:
	InterpolateNoise2(double deltat,double prebuffer,double psd,double filterexp,int interplen = 1);
	InterpolateNoise2(double deltat,double prebuffer,double psd,Filter *myfilter,int interplen = 1);
	InterpolateNoise2(double deltat,double prebuffer,double psd,double filterexp,Interpolator *myinterp);
	InterpolateNoise2(double deltat,double prebuffer,double psd,Filter *myfilter,Interpolator *myinterp);

	InterpolateNoise2(double *narray,long length,double deltat,double prebuffer,double norm,double filterexp = 0.0,int interplen = 1);
	InterpolateNoise2(double *narray,long length,double deltat,double prebuffer,double norm,Filter *myfilter,int interplen = 1);
	InterpolateNoise2(double *narray,long length,double deltat,double prebuffer,double norm,double filterexp,Interpolator *myinterp);
	InterpolateNoise2(double *narray,long length,double deltat,double prebuffer,double norm,Filter *myfilter,Interpolator *myinterp);

	virtual ~InterpolateNoise2();

    void reset();

    double operator[](double time);

    double noise(double time) {
		return (*this)[time];
    };

    double noise(double timebase,double timecorr);

    void setinterp(Interpolator *myinterp);
};

#endif /* _LISASIM_NOISE_H_ */
