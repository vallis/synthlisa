#ifndef _LISASIM_NOISE_H_
#define _LISASIM_NOISE_H_

#include <iostream>
using namespace std;

#include <limits.h>
#include "gsl_rng.h"

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


class BufferNoise {
 private:
    double *data;
    long maxbuffer;

    double zero;

 public:
    BufferNoise(long l) : maxbuffer(l) {
	if ( !(data = new double[maxbuffer]) ) {
	    cout << "RingBufferNoise::RingBufferNoise: out of memory" << endl;
	    abort();
	}

	for(int i=0;i<maxbuffer;i++)
	    data[i] = 0.0;
    };

    ~BufferNoise() {
	delete data;
    };

    void reset() {
	for(int i=0;i<maxbuffer;i++)
	    data[i] = 0.0;
    };

    double& operator[](long pos) {
	if(pos < 0) {
	    zero = 0.0;
	    return zero;
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
	if (pos < 0 && pos > maxbuffer) {
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
    virtual void operator()(BufferNoise &x,BufferNoise &y,long pos) = 0;
};

class NoFilter : public Filter {
 public:
    void operator()(BufferNoise &x,BufferNoise &y,long pos) {
	y[pos] = x[pos];
    };
};

class IntFilter : public Filter {
 private:
    double alpha;

 public:
    IntFilter(double a = 0.9999) : alpha(a) {};

    void operator()(BufferNoise &x,BufferNoise &y,long pos) {
	y[pos] = alpha * y[pos - 1] + x[pos];
    };
};

class DiffFilter : public Filter {
 public:
    void operator()(BufferNoise &x,BufferNoise &y,long pos) {
	y[pos] = x[pos] - x[pos-1];
    };
};

// the main fixed-timestep noise class

class MakeNoise {
 private:
    BufferNoise *x, *y;
    long maxbuffer;
    long earliest, latest;

    GetNoise *get;
    Filter *filter;

 public:
    MakeNoise(GetNoise *g,Filter *f,long l)
	: maxbuffer(l), earliest(-1), latest(-1), get(g), filter(f) {
	x = new BufferNoise(maxbuffer);
	y = new BufferNoise(maxbuffer);
    };

    ~MakeNoise() {
	delete y;
	delete x;
    }

    void reset() {
	get->reset();

	x->reset();
	y->reset();

	earliest = -1;
	latest = -1;
    };

    double operator[](long pos) {
	if(pos < earliest) {
	    cout << "MakeNoise::[] trying to access noise element (" << pos << ") before oldest kept." << endl;
	    abort();
	} else if (pos > latest) {
	    for(int i=latest+1;i<=pos;i++) {
		(*get)(*x,i);
		(*filter)(*x,*y,i);
	    }

	    latest = pos;
	    if(latest - earliest > maxbuffer)
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
    virtual double operator()(MakeNoise &y,long ind,double dind) = 0;
};

class NearestInterpolator : public Interpolator {
 public:
    double operator()(MakeNoise &y,long ind,double dind) {
	return (dind < 0.5 ? y[ind] : y[ind+1]);
    };
};

// 0 < dind < 1; the desired sample would be at ind + dind

class LinearInterpolator : public Interpolator {
 public:
    double operator()(MakeNoise &y,long ind,double dind) {
	return (1.0 - dind) * y[ind] + dind * y[ind+1];
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
};

// dimensioned interpolated noise

class InterpolateNoise : public Noise {
 private:
    double samplingtime;
    double nyquistf;

    double prebuffertime;
    double maxtime;

    double normalize;
  
    GetNoise *getnoise;
    Filter *thefilter;

    MakeNoise *thenoise;
    Interpolator *interp;

    void setfilter(double ex);
    void setnorm(double sd, double ex);
    void setinterp(int window);

 public:
    InterpolateNoise(double sampletime,double prebuffer,double density,double exponent, int swindow = 1);

    InterpolateNoise(double *noisebuf,long samples,double sampletime,double prebuffer,double density, double exponent = 0.0, int swindow = 1);

    virtual ~InterpolateNoise();
    
    void reset() {
	thenoise->reset();
    };

    double operator[](double time);

    double noise(double time) {
	return (*this)[time];
    };
};


#endif /* _LISASIM_NOISE_H_ */
