#include "lisasim-noise.h"

#include <iostream>
using namespace std;

#include <cmath>
#include <stdlib.h>

#include <float.h>

// from contrib-source/GSL-1.4
#include "gsl_rng.h"

// for seedrandgen
#include <sys/time.h>

// --- SampledNoise ----------------------------------------------------

SampledNoise::SampledNoise(double *nb, long samples) {
    noisebuffer = nb;

    maxsamples = samples;
}

double SampledNoise::operator[](long pos) {
    if (pos > maxsamples) {
      cout << "SampledNoise::[] accessing index larger than available" << endl;
      abort();
    }

    return noisebuffer[pos];
}


// --- "Ring" classes --------------------------------------------------

//  default constructor for RingNoise

RingNoise::RingNoise(long bs) {
    randgen = gsl_rng_alloc(gsl_rng_ranlux);
    seedrandgen();

    buffersize = bs;
    
    bufferx = new double[bs];
    buffery = new double[bs];
    
    for(int i=0; i<bs; i++) {
        bufferx[i] = 0.00;
        buffery[i] = 0.00;
    }
    
    earliest = -1;
    latest = -1;
}

RingNoise::~RingNoise() {
    delete [] buffery;
    delete [] bufferx;

    gsl_rng_free(randgen);
}

void RingNoise::reset() {
    seedrandgen();

    for(int i=0; i<buffersize; i++) {
        bufferx[i] = 0.00;
        buffery[i] = 0.00;
    }
    
    earliest = -1;
    latest = -1;    
}

void RingNoise::updatebuffer(long pos) {
    for(int i=latest+1; i<=pos; i++) {
        bufferx[i % buffersize] = deviate();
        buffery[i % buffersize] = filter(i);
    }

    latest = pos;    

    if (latest - earliest >= buffersize)
        earliest = latest - buffersize + 1;
}

double RingNoise::operator[](long pos) {
    if (pos > latest) {
        updatebuffer(pos);
    } if(pos < earliest) {
        cout << "RingNoise::[] trying to access element before oldest kept." << endl;
        abort();
    } else {
        return buffery[pos % buffersize];
    }
}

void RingNoise::seedrandgen() {
    struct timeval tv;

    gettimeofday(&tv,0);
  
    // printf("seconds is: %ld, microseconds is: %ld\n",tv.tv_sec,tv.tv_usec);
    // printf("seed is: %ld\n",tv.tv_sec+tv.tv_usec);
  
    gsl_rng_set(randgen,tv.tv_sec+tv.tv_usec); 
    cacheset = 0;
}

// Box-Muller transform to get Gaussian deviate from uniform random number
// adapted from GSL 1.4 randist/gauss.c

double RingNoise::deviate() {
  double x, y, r2;

  // this will never return a cached 0.0, but what's the chance?

  if (cacheset == 0) {
      do {
	  x = -1.0 + 2.0 * gsl_rng_uniform(randgen);
	  y = -1.0 + 2.0 * gsl_rng_uniform(randgen);
	  
	  r2 = x * x + y * y;
      } while (r2 > 1.0 || r2 == 0);

      double root = sqrt (-2.0 * log (r2) / r2);

      cacheset = 1;
      cacherand = x * root;

      return y * root;
  } else {
      cacheset = 0;
      return cacherand;
  }
}

// the default is to apply no filter (yielding pure uncorrelated Gaussian white noise)

double RingNoise::filter(long pos) {
    return bufferx[pos % buffersize];
}

// derived classes: differentiating and integrating filters

DiffNoise::DiffNoise(long bs) : RingNoise(bs) {};

double DiffNoise::filter(long pos) {
    return (bufferx[pos % buffersize] - bufferx[(pos + buffersize - 1) % buffersize]);
}

IntNoise::IntNoise(long bs, double ic) : RingNoise(bs) {
    intconst = ic;
};

double IntNoise::filter(long pos) {
    return (intconst * buffery[(pos + buffersize - 1) % buffersize] + bufferx[pos % buffersize]);
}

// --- dimensioned-noise classes --------------------------------------------------

// constructor based on SampledNoise: need pointer to noise buffer, number of samples, sampling time,
// requested prebuffer time (probably a multiple of some nominal armlength), and the one-sided power spectral density

InterpolateNoise::InterpolateNoise(double *noisebuf,long samples,double st,double pbt,double norm) {
    samplingtime = st;
    nyquistf = 1.00 / (2*st);

    // compute the prebuffer offset

    prebuffertime = long(floor(1e-10 + pbt / samplingtime));

    // the corresponding maximum sample time that can be accessed is
    
    maxitime = samples - prebuffertime - 1;    

    buffernoise = new SampledNoise(noisebuf,samples);
    normalize = norm;
}

// constructor based on RingNoise: need to pass the sampling time (in secs), the requested prebuffer time
// (probably a multiple of some nominal armlength), and the one-sided power spectral density,
// given as sd*(f/Hz)^ex, where sd is in Hz^-1, and ex is either 0.00, -2.00, or 2.00.
// The corresponding ringnoise object is then allocated automatically

InterpolateNoise::InterpolateNoise(double st, double pbt, double sd, double ex) {
    samplingtime = st;
    nyquistf = 1.00 / (2*st);

    // compute the prebuffer offset

    prebuffertime = long(floor(1e-10 + pbt / samplingtime));

    // the corresponding maximum sample time that can be accessed is
    
    maxitime = LONG_MAX * 1.00 - prebuffertime - 1;    

    // set the RingNoise buffer as a power of two larger than the offset

    long rnbuffer = 2;
    while (prebuffertime > rnbuffer) rnbuffer *= 2;

    // create the RingNoise object corresponding to the requested digital filter
    // do not have 'switch' for floats; manage with if
    
    if (ex == 0.00) {
        buffernoise = new RingNoise(rnbuffer);
        normalize = sqrt(sd) * sqrt(nyquistf);
    } else if (ex == 2.00) {
        buffernoise = new DiffNoise(rnbuffer);
        normalize = sqrt(sd) * sqrt(nyquistf) / (2.00 * M_PI * samplingtime);
    } else if (ex == -2.00) {
        buffernoise = new IntNoise(rnbuffer,0.9999);
        normalize = sqrt(sd) * sqrt(nyquistf) * (2.00 * M_PI * samplingtime);
    } else {
        cout << "InterpolateNoise::InterpolateNoise() spectral shape f^" << ex << " not implemented!" << endl;
        abort();
    }
}

InterpolateNoise::~InterpolateNoise() {
    delete buffernoise;
}

void InterpolateNoise::reset() {
    buffernoise->reset();
}

double InterpolateNoise::inoise(double time) {
    double ctime = time / samplingtime;
    double itime = floor(ctime);

    long index;
    
    if (itime > maxitime) {
        cout << "InterpolateNoise::inoise() accessing index larger than maximum allowed" << endl;
        abort();
    } else {
        index = long(itime) + prebuffertime;
    }

    double interp = (ctime - itime) * (*buffernoise)[index+1] + (1.00 - (ctime - itime)) * (*buffernoise)[index];

    return( normalize * interp );
}

double InterpolateNoise::operator[](double time) {
    return inoise(time);
}

// --- windowed-sinc interpolation ------------------------------------------------

InterpolateNoiseBetter::InterpolateNoiseBetter(double sampletime,double prebuffer,double density,double exponent,int swindow) : InterpolateNoise(sampletime,prebuffer,density,exponent) {
    semiwindow = swindow;
    semiconst = M_PI / swindow;
}

InterpolateNoiseBetter::InterpolateNoiseBetter(double *noisebuf,long samples,double sampletime,double prebuffer,double norm,int swindow) : InterpolateNoise(noisebuf,samples,sampletime,prebuffer,norm) {
    semiwindow = swindow;
    semiconst = M_PI / swindow;
}

inline double InterpolateNoiseBetter::sinc(double time) {
    if(time==0.0)
	return 1.0;
    else
	return ( sin(M_PI * time)/(M_PI * time) );
}

inline double InterpolateNoiseBetter::window(double time) {
    return ( 0.54 + 0.46 * cos(time * semiconst) );
}

double InterpolateNoiseBetter::inoise(double time) {
    double ctime = time / samplingtime;
    double itime = floor(ctime);

    long index;
    
    if (itime > maxitime) {
        cout << "InterpolateNoise::inoise() accessing index larger than maximum allowed" << endl;
        abort();
    } else {
        index = long(itime) + prebuffertime;
    }

    // or should I implement this with a tolerance?

    if (itime == ctime) {
	return ( normalize * (*buffernoise)[index] );
    } else {
	double interp = 0.0;
	double t;

	// here we work with integer bin abscissae

	for(int i=0;i<semiwindow;i++) {
	    t = ctime - itime + i;
	    interp += sinc(t) * window(t) * (*buffernoise)[index-i];
	    
	    t = -(1.00 - (ctime - itime)) - i;
	    interp += sinc(t) * window(t) * (*buffernoise)[index+1+i];
	}
	
	return ( normalize * interp );
    }
}

double InterpolateNoiseBetter::operator[](double time) {
    return inoise(time);
}
