/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#include "lisasim-noise.h"

#include <iostream>
using namespace std;

#include <cmath>
#include <stdlib.h>
#include <limits.h>

#include <float.h>

// from contrib-source/GSL-1.4

#include "GSL/gsl_rng.h"

// for seedrandgen

#include <sys/time.h>


WhiteNoise::WhiteNoise(unsigned long seed) {
    randgen = gsl_rng_alloc(gsl_rng_ranlux);

    seedrandgen(seed);
}

WhiteNoise::~WhiteNoise() {
    gsl_rng_free(randgen);
}

void WhiteNoise::seedrandgen(unsigned long seed) {
    if (seed == 0) {
		struct timeval tv;
	
		gettimeofday(&tv,0);
	  
		gsl_rng_set(randgen,tv.tv_sec+tv.tv_usec);
    } else {
		gsl_rng_set(randgen,seed);
    }

    cacheset = 0;
    cacherand = 0;
}

// Box-Muller transform to get Gaussian deviate from uniform random number
// adapted from GSL 1.4 randist/gauss.c

double WhiteNoise::deviate() {
  double x, y, r2;

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

// constructor based on WhiteNoise. Need to pass:
// - the sampling time (in secs)
// - the requested prebuffer time (in secs)
// - the one-sided power spectral density, given as sd*(f/Hz)^ex,
//   where sd is in Hz^-1, and ex is either 0.00, -2.00, or 2.00.
// - the noise exponent ex

// modified from Numerical Recipes

LagrangeInterpolator::LagrangeInterpolator(int sw)
    : window(2*sw), semiwindow(sw) {
	xa = new double[window+1];
	ya = new double[window+1];
		
	c = new double[window+1];
	d = new double[window+1];
		
	for(int i=1;i<=window;i++) {
		xa[i] = 1.0*i;
		ya[i] = 0.0;
	}    
};

LagrangeInterpolator::~LagrangeInterpolator() {
    delete [] d;
    delete [] c;
	
    delete [] ya;
    delete [] xa;
};

double LagrangeInterpolator::polint(double x) {
    int n = window;
    int i,m,ns=1;
    double den,dif,dift,ho,hp,w;

    double res,dres; // dres is the error estimate

    dif=fabs(x-xa[1]);

    for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
	
		c[i]=ya[i];
		d[i]=ya[i];
    }

    res=ya[ns--];

    for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			den=ho-hp;
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
	
		res += (dres=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }

    return res;
}

// InterpolateNoise2

Filter* InterpolateNoise2::getfilter(double filterexp) {
    if (filterexp == 0.00)
		return new NoFilter();
    else if (filterexp == 2.00)
		return new DiffFilter();
    else if (filterexp == -2.00)
		return new IntFilter();
    else {
		cout << "InterpolateNoise::InterpolateNoise: noise spectral shape f^"
			 << filterexp << " not implemented. Defaulting to no filtering." << endl;

		return new NoFilter();
    }
}

Interpolator* InterpolateNoise2::getinterp(int interplen) {
    if (interplen == 0)
    	return new NearestInterpolator();
    else if (interplen == -1)
		return new LinearExtrapolator();
    else if (interplen == 1)
		return new LinearInterpolator();
    else
		return new LagrangeInterpolator(interplen);
}

double InterpolateNoise2::getnormrandom(double psd,double filterexp) {
    if (filterexp == 0.00)
        return sqrt(psd) * sqrt(nyquistf);
    else if (filterexp == 2.00)
        return sqrt(psd) * sqrt(nyquistf) / (2.00 * M_PI * samplingtime);
    else if (filterexp == -2.00)
        return sqrt(psd) * sqrt(nyquistf) * (2.00 * M_PI * samplingtime);
    else {
        cout << "InterpolateNoise::InterpolateNoise: noise spectral shape f^"
        	 << filterexp << " not implemented. Defaulting to no white normalization." << endl;

        return sqrt(psd) * sqrt(nyquistf);
    }
}

double InterpolateNoise2::getnormsampled(double norm,double filterexp) {
    if (filterexp == 0.00)
        return norm;                // we're just normalizing
    else if (filterexp == 2.00)
        return norm / samplingtime; // we're differentiating
    else if (filterexp == -2.00)
        return norm * samplingtime; // we're integrating
    else {
        cout << "InterpolateNoise::InterpolateNoise: noise spectral shape f^"
        	 << filterexp << " not implemented. Defaulting to no filtering." << endl;
        
        return norm;
    }
}

void InterpolateNoise2::settime() {
    maxtime = samplingtime * (efflength - 1) - prebuffertime;
    
    lasttime = 0.0;

	timewindow = prebuffertime - samplingtime * (interp->length());
}

void InterpolateNoise2::initrandom() {
	getnoise = new WhiteNoise();
    thenoise = new FilterMakeNoise(getnoise,filter,long(prebuffertime/samplingtime));
    
    efflength = LONG_MAX;
    settime();
}

void InterpolateNoise2::initsampled(double *narray,long length) {
    getnoise = new SampledNoise(narray,length);
    thenoise = new FilterMakeNoise(getnoise,filter,length);

	efflength = length;
	settime();
}

InterpolateNoise2::InterpolateNoise2(double deltat,double prebuffer,double psd,double filterexp,int interplen)
	: samplingtime(deltat), nyquistf(0.5 / deltat), prebuffertime(prebuffer) {

	filter = getfilter(filterexp); ownfilter = 1;
	interp = getinterp(interplen); ownfilter = 1;

	normalize = getnormrandom(psd,filterexp);

	initrandom();
}

InterpolateNoise2::InterpolateNoise2(double deltat,double prebuffer,double psd,Filter *myfilter,int interplen)
	: samplingtime(deltat), nyquistf(0.5 / deltat), prebuffertime(prebuffer) {

	filter = myfilter; ownfilter = 0;
	interp = getinterp(interplen); ownfilter = 1;

	normalize = getnormrandom(psd,0.0);

	initrandom();
}

InterpolateNoise2::InterpolateNoise2(double deltat,double prebuffer,double psd,double filterexp,Interpolator *myinterp)
	: samplingtime(deltat), nyquistf(0.5 / deltat), prebuffertime(prebuffer) {

	filter = getfilter(filterexp); ownfilter = 1;
	interp = myinterp; ownfilter = 0;

	normalize = getnormrandom(psd,filterexp);

	initrandom();
}

InterpolateNoise2::InterpolateNoise2(double deltat,double prebuffer,double psd,Filter *myfilter,Interpolator *myinterp)
	: samplingtime(deltat), nyquistf(0.5 / deltat), prebuffertime(prebuffer) {

	filter = myfilter; ownfilter = 0;
	interp = myinterp; ownfilter = 0;

	normalize = getnormrandom(psd,0.0);

	initrandom();
}

InterpolateNoise2::InterpolateNoise2(double *narray,long length,double deltat,double prebuffer,double norm,double filterexp,int interplen)
	: samplingtime(deltat), nyquistf(0.5 / deltat), prebuffertime(prebuffer) {

	filter = getfilter(filterexp); ownfilter = 1;
	interp = getinterp(interplen); ownfilter = 1;

	normalize = getnormsampled(norm,filterexp);

	initsampled(narray,length);
}

InterpolateNoise2::InterpolateNoise2(double *narray,long length,double deltat,double prebuffer,double norm,Filter *myfilter,int interplen)
	: samplingtime(deltat), nyquistf(0.5 / deltat), prebuffertime(prebuffer) {

	filter = myfilter; ownfilter = 0;
	interp = getinterp(interplen); ownfilter = 1;

	normalize = getnormsampled(norm,0.0);

	initsampled(narray,length);
}

InterpolateNoise2::InterpolateNoise2(double *narray,long length,double deltat,double prebuffer,double norm,double filterexp,Interpolator *myinterp)
	: samplingtime(deltat), nyquistf(0.5 / deltat), prebuffertime(prebuffer) {

	filter = getfilter(filterexp); ownfilter = 1;
	interp = myinterp; ownfilter = 0;

	normalize = getnormsampled(norm,filterexp);

	initsampled(narray,length);
}


InterpolateNoise2::InterpolateNoise2(double *narray,long length,double deltat,double prebuffer,double norm,Filter *myfilter,Interpolator *myinterp)
	: samplingtime(deltat), nyquistf(0.5 / deltat), prebuffertime(prebuffer) {

	filter = myfilter; ownfilter = 0;
	interp = myinterp; ownfilter = 0;

	normalize = getnormsampled(norm,0.0);

	initsampled(narray,length);
}

InterpolateNoise2::~InterpolateNoise2() {
    delete thenoise;
    delete getnoise;
    if (owninterp) delete interp;
    if (ownfilter) delete filter;
}

void InterpolateNoise2::reset() {
    thenoise->reset();

    lasttime = 0.0;
};

// it would be better to throw an exception!

double InterpolateNoise2::operator[](double time) {
	if (normalize == 0.0) return 0.0;

    if (time > maxtime) {
        cout << "InterpolateNoise::[]: time requested (" << time <<
	    		") too large, maxtime is " << maxtime <<
	    		" [" << __FILE__ << ":" << __LINE__ << "]" << endl;
	    		
        abort();
    } else if (lasttime - time > timewindow) {
		cout << "InterpolateNoise::[]: time requested (" << time <<
	    		") too old, last was " << lasttime <<
	    		", prebuffer is " << prebuffertime <<
	    		", timewindow is " << timewindow <<
	    		" [" << __FILE__ << ":" << __LINE__ << "]" << endl;
	    		
		abort();
    } else {
    	if (time > lasttime) lasttime = time;

		double rind = (time + prebuffertime) / samplingtime;
		double dind = rind - floor(rind);
	
		// the floor is to make sure we can handle negative indices
		// correctly
		long ind = long(floor(rind));
	
		return normalize * (*interp)(*thenoise,ind,dind);
    }
}

double InterpolateNoise2::noise(double timebase,double timecorr) {
	if (normalize == 0.0) return 0.0;

    double time = timebase + timecorr;

    if (time > maxtime) {
        cout << "InterpolateNoise::[]: time requested (" << time << ") too large" <<
	    "[" << __FILE__ << ":" << __LINE__ << "]" << endl;
	    
        abort();
    } else if (lasttime - time > timewindow) {
		cout << "InterpolateNoise::[]: time requested (" << time << 
	    ") too old, last was " << lasttime << ", pbt is " << prebuffertime <<
	    "[" << __FILE__ << ":" << __LINE__ << "]" << endl;
	    
		abort();
    } else {
	    if (time > lasttime) lasttime = time;

		double rindbase = (timebase + prebuffertime) / samplingtime;
		double dindbase = rindbase - floor(rindbase);
	
		double rindcorr = timecorr / samplingtime;
		double dindcorr = rindcorr - floor(rindcorr);
	
		double rind = rindbase + rindcorr + floor(dindbase + dindcorr);
		double dind = dindbase + dindcorr - floor(dindbase + dindcorr);
	
		// the floor is to make sure we can handle negative indices
		// correctly
		long ind = long(floor(rind));
	
		return normalize * (*interp)(*thenoise,ind,dind);
    }
}

// ??? It should be verified whether this works or not!

void InterpolateNoise2::setinterp(Interpolator *myinterp) {
	if (owninterp) delete interp;
	
	interp = myinterp; owninterp = 1;

	double newtimewindow = prebuffertime - samplingtime * (interp->length());

	if (newtimewindow > timewindow)
		lasttime = lasttime + (newtimewindow - timewindow);

	timewindow = newtimewindow;
}


InterpolateNoise::InterpolateNoise(double st, double pbt, double sd, double ex, int win) {
    samplingtime = st;
    nyquistf = 0.5 / st;

    prebuffertime = pbt;
    maxtime = samplingtime * (LONG_MAX - 1) - prebuffertime;    
    lasttime = 0.0;

    setfilter(ex);
    setnorm(sd,ex);

    interp = 0;
    timewindow = 0.0;
    setinterp(win);

    // what about we set this to exactly what requested?
    long buffersize = long(prebuffertime/samplingtime);

    // long buffersize = 2;
    // while (buffersize < long(prebuffertime/samplingtime))
    //     buffersize *= 2;

    // WhiteNoise object

    getnoise = new WhiteNoise();

    thenoise = new FilterMakeNoise(getnoise,filter,buffersize);
}

// constructor based on SampledNoise. Need to pass:
// - pointer to noise buffer
// - number of samples
// - sampling time
// - requested prebuffer time
// - normalization factor
// - the noise exponent ex, which determines the filtering

InterpolateNoise::InterpolateNoise(double *nb,long sl,double st,double pbt,double sd,double ex,int win) {
    samplingtime = st;
    nyquistf = 0.5 / st;

    prebuffertime = pbt;
    maxtime = samplingtime * (sl - 1) - prebuffertime;
    lasttime = 0.0;

    setfilter(ex);
    setnormsampled(sd,ex);

    interp = 0;
    timewindow = 0.0;
    setinterp(win);

    // SampledNoise object

    getnoise = new SampledNoise(nb,sl);

    long &buffersize = sl;
    thenoise = new FilterMakeNoise(getnoise,filter,buffersize);
}

InterpolateNoise::~InterpolateNoise() {
    delete interp;
    delete thenoise;
    delete filter;
    delete getnoise;
}

void InterpolateNoise::reset() {
    thenoise->reset();
    lasttime = 0.0;
};

void InterpolateNoise::setfilter(double ex) {
    if (ex == 0.00) {
		filter = new NoFilter();
    } else if (ex == 2.00) {
		filter = new DiffFilter();
    } else if (ex == -2.00) {
		filter = new IntFilter();
    } else {
		cout << "InterpolateNoise::InterpolateNoise: noise spectral shape f^"
		     << ex << " not implemented. Defaulting to no filtering." << endl;
		filter = new NoFilter();
    }
}

void InterpolateNoise::setnorm(double sd, double ex) {
    if (ex == 0.00) {
        normalize = sqrt(sd) * sqrt(nyquistf);
    } else if (ex == 2.00) {
        normalize = sqrt(sd) * sqrt(nyquistf) / (2.00 * M_PI * samplingtime);
    } else if (ex == -2.00) {
        normalize = sqrt(sd) * sqrt(nyquistf) * (2.00 * M_PI * samplingtime);
    } else {
        cout << "InterpolateNoise::InterpolateNoise: noise spectral shape f^"
        	 << ex << " not implemented. Defaulting to no filtering." << endl;
        normalize = sqrt(sd) * sqrt(nyquistf);
    }
}

void InterpolateNoise::setnormsampled(double sd, double ex) {
    if (ex == 0.00) {
        normalize = sd;                // we're just normalizing
    } else if (ex == 2.00) {
        normalize = sd / samplingtime; // we're differentiating
    } else if (ex == -2.00) {
        normalize = sd * samplingtime; // we're integrating
    } else {
        cout << "InterpolateNoise::InterpolateNoise: noise spectral shape f^"
        	 << ex << " not implemented. Defaulting to no filtering." << endl;
        normalize = sqrt(sd) * sqrt(nyquistf);
    }
}

void InterpolateNoise::setinterp(int window) {
    if (interp != 0) delete interp;

    if (window == 0) {
		interp = new NearestInterpolator();
    } else if (window == -1) {
		interp = new LinearExtrapolator();
    } else if (window == 1) {
		interp = new LinearInterpolator();
    } else {
		interp = new LagrangeInterpolator(window);
    }

    // Here we set the window of interpolation times that will be
    // accessible going back from the last requested time.  We need to
    // be careful if we're changing the interpolation scheme
    // dynamically.

    // Intuitively, the -1 extrapolator should be equivalent to a 2 interpolator
    if (window == -1) window = 2;

    double tw = prebuffertime - 2.0 * window * samplingtime;

    if (timewindow > 0.0 && tw > timewindow) {
		lasttime = lasttime + (tw - timewindow);
    }

    timewindow = tw;
}

// should update the paper to say that the earliest accessible time
// actually depends on the interpolation scheme!

double InterpolateNoise::operator[](double time) {
	if (normalize == 0.0) return 0.0;

    if (time > maxtime) {
        cout << "InterpolateNoise::[]: time requested (" << time <<
	    		") too large [" << __FILE__ << ":" << __LINE__ << "]" << endl;
	    		
        abort();
    } else if (lasttime - time > timewindow) {
		cout << "InterpolateNoise::[]: time requested (" << time <<
	    		") too old, last was " << lasttime << ", pbt is " << prebuffertime <<
	    		" [" << __FILE__ << ":" << __LINE__ << "]" << endl;
	    		
		abort();
    } else {
    	if (time > lasttime) lasttime = time;

		double rind = (time + prebuffertime) / samplingtime;
		double dind = rind - floor(rind);
	
		// the floor is to make sure we can handle negative indices
		// correctly
		long ind = long(floor(rind));
	
		return normalize * (*interp)(*thenoise,ind,dind);
    }
}

double InterpolateNoise::noise(double timebase,double timecorr) {
	if (normalize == 0.0) return 0.0;

    double time = timebase + timecorr;

    if (time > maxtime) {
        cout << "InterpolateNoise::[]: time requested (" << time << ") too large" <<
	    "[" << __FILE__ << ":" << __LINE__ << "]" << endl;
	    
        abort();
    } else if (lasttime - time > timewindow) {
		cout << "InterpolateNoise::[]: time requested (" << time << 
	    ") too old, last was " << lasttime << ", pbt is " << prebuffertime <<
	    "[" << __FILE__ << ":" << __LINE__ << "]" << endl;
	    
		abort();
    } else {
	    if (time > lasttime) lasttime = time;

		double rindbase = (timebase + prebuffertime) / samplingtime;
		double dindbase = rindbase - floor(rindbase);
	
		double rindcorr = timecorr / samplingtime;
		double dindcorr = rindcorr - floor(rindcorr);
	
		double rind = rindbase + rindcorr + floor(dindbase + dindcorr);
		double dind = dindbase + dindcorr - floor(dindbase + dindcorr);
	
		// the floor is to make sure we can handle negative indices
		// correctly
		long ind = long(floor(rind));
	
		return normalize * (*interp)(*thenoise,ind,dind);
    }
}
