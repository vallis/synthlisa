/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#include "lisasim-signal.h"
#include "lisasim-except.h"

#include <iostream>
#include <cmath>

#include <sys/time.h>

// --- RingBuffer ---

RingBuffer::RingBuffer(long len)
	: data(new double[len]), length(len) {

	reset();
}

RingBuffer::~RingBuffer() {
	delete [] data;
}

void RingBuffer::reset() {
	for(int i=0;i<length;i++) data[i] = 0.0;
}


// --- BufferedSignalSource ---

BufferedSignalSource::BufferedSignalSource(long len)
	: buffer(len), length(len), current(-1) {}

void BufferedSignalSource::reset(unsigned long seed) {
	buffer.reset();
	
	current = -1;
}

double BufferedSignalSource::operator[](long pos) {
	if (pos <= current - length) {
		std::cerr << "BufferedSignalSource::operator[](long): stale sample access at "
		          << pos << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;

		ExceptionOutOfBounds e;
		throw e;
	} else if (pos > current) {
		for(int i=current+1;i<=pos;i++) {
			buffer[i] = getvalue(i);
		}
		
		current = pos;

		return buffer[pos];
	} else {
		return buffer[pos];
	}
}


// --- WhiteNoiseSource ---

WhiteNoiseSource::WhiteNoiseSource(long len,unsigned long seed,double norm) 
	: BufferedSignalSource(len), normalize(norm) {

	randgen = gsl_rng_alloc(gsl_rng_taus2);
    // randgen = gsl_rng_alloc(gsl_rng_ranlux);

    seedrandgen(seed);
}

WhiteNoiseSource::~WhiteNoiseSource() {
    gsl_rng_free(randgen);
}

unsigned long WhiteNoiseSource::globalseed = 0;

void WhiteNoiseSource::setglobalseed(unsigned long seed) {
	if(seed == 0) {
		struct timeval tv;
	
		gettimeofday(&tv,0);	

		globalseed = tv.tv_sec+tv.tv_usec;
	} else {
		globalseed = seed;
	}
}

unsigned long WhiteNoiseSource::getglobalseed() {
	if(globalseed == 0) setglobalseed(0);
	
	return globalseed;
}

void WhiteNoiseSource::seedrandgen(unsigned long seed) {
    if (seed == 0) {
		gsl_rng_set(randgen,globalseed);
		
		// should replace this simple seed enumeration
		// with something a bit more sophisticated
		globalseed++;
    } else {
		gsl_rng_set(randgen,seed);
    }

    cacheset = 0;
    cacherand = 0;
}

void WhiteNoiseSource::reset(unsigned long seed) {
	seedrandgen(seed);
	
	BufferedSignalSource::reset(seed);
}

/* Box-Muller transform to get Gaussian deviate from uniform random,
   number, adapted from GSL 1.4 randist/gauss.c

   Since pos is not actually used, this will give sensible results only
   within a buffered environment */

double WhiteNoiseSource::getvalue(long pos) {
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

      return normalize * y * root;
  } else {
      cacheset = 0;
      return normalize * cacherand;
  }
}


// --- ResampledSignalSource ---

// the idea is to take a Signal (continuous) and to use it to feed a
// SignalSource (discrete), with the eventual purpose of using an
// InterpolatedSignal on top of that, to cache it

double ResampledSignalSource::getvalue(long pos) {
    return signal->value(pos*deltat - prebuffer);
}

void ResampledSignalSource::reset(unsigned long seed) {
	signal->reset(seed);

	BufferedSignalSource::reset(seed);
}


// --- FileSignalSource ---

// reads in batches from a file, owns the data buffer array

/* Endianness conversion: a bit rough, but should work */
/* From http://www.codeproject.com/cpp/endianness.asp */

#define BIGENDIAN      0
#define LITTLEENDIAN   1

static int testbyteorder()
{
   short int word = 0x0001;
   char *byte = (char *) &word;
   return(byte[0] ? LITTLEENDIAN : BIGENDIAN);
}

/* This will work as long as doubles are 8 bytes */

static void convertendianness(double *val) {
    unsigned char cval[8];
    
    cval[0] = ((unsigned char *)val)[7];
    cval[1] = ((unsigned char *)val)[6];
    cval[2] = ((unsigned char *)val)[5];
    cval[3] = ((unsigned char *)val)[4];
    cval[4] = ((unsigned char *)val)[3];    
    cval[5] = ((unsigned char *)val)[2];
    cval[6] = ((unsigned char *)val)[1];
    cval[7] = ((unsigned char *)val)[0];
    
    *val = *((double *)cval);
}

// pass 0 for bigendian, 1 for littleendian, -1 for native

FileSignalSource::FileSignalSource(char *filename,long bufferlen,long prebuffer,int endianness,double norm)
    : BufferedSignalSource(prebuffer) {
    file = fopen(filename,"rb");
    if(file == NULL) {
		std::cerr << "FileSignalSource::FileSignalSource(...): cannot open file "
		          << filename << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;
	
		ExceptionFileError e;
		throw e;
    }
    
    length = bufferlen;    
    filebuffer = new double[length];
    
    if(endianness == -1 || endianness == testbyteorder()) {
        convert = 0;
    } else {
        convert = 1;
    }
    
    initpos = 0;
    loadbuffer();
    
    normalize = norm;
}

FileSignalSource::~FileSignalSource() {
    delete [] filebuffer;

    fclose(file);
}

void FileSignalSource::loadbuffer() {
    long actuallyread = fread(filebuffer,sizeof(double),length,file);

    // convert endianness if necessary

    if(convert == 1) {
        for(int i=0;i<actuallyread;i++)
            convertendianness(&filebuffer[i]);
    }

    // pad with zeros if we haven't read enough

    if(actuallyread < length) {
        for(int i=actuallyread;i<length;i++)
            filebuffer[i] = 0.0;
    } 
}

double FileSignalSource::getvalue(long pos) {
    while(pos - initpos >= length) {
        initpos += length;
        loadbuffer();
    }
    
    return normalize * filebuffer[pos - initpos];
}

void FileSignalSource::reset(unsigned long seed ) {
    rewind(file);
    
    initpos = 0;
    loadbuffer();
}


// --- SampledSignalSource ---

// does not own the data array!

SampledSignalSource::SampledSignalSource(double *darray,long len,double norm)
	: data(darray), length(len), warn(0), normalize(norm) {}

// will pad with zeros on the negative index side
// NEW BEHAVIOR: will pad with zeros outside the sampled length, but will
// warn once

double SampledSignalSource::operator[](long pos) {
	if(pos < 0) {
		/* if(!(warn & 1)) {
			std::cerr << "SampledSignalSource::operator[](long): padding with zeros at pos "
		          << pos << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;

			warn |= 1;
		} */

		return 0.0;
	} else if(pos >= length) {
		/* if(!(warn & 2)) {
			std::cerr << "SampledSignalSource::operator[](long): padding with zeros at pos "
		          << pos << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;

			warn |= 2;
		} */
			
		return 0.0;
	} else {
		return normalize * data[pos];
	}
}


// --- Filters ---

double NoFilter::getvalue(SignalSource &x,SignalSource &y,long pos) {
	return x[pos];
}


IntFilter::IntFilter(double a)
	: alpha(a) {}

double IntFilter::getvalue(SignalSource &x,SignalSource &y,long pos) {
	return alpha * y[pos - 1] + x[pos];
}


double DiffFilter::getvalue(SignalSource &x,SignalSource &y,long pos) {
	return x[pos] - x[pos-1];
}


/* Note that FIR makes a copy of the coefficient array. Normally filters
   are defined so that a[0] = 1 */

FIRFilter::FIRFilter(double *aarray,int len)
	: a(new double[len]), length(len) {
		
	for(int i=0;i<length;i++)
		a[i] = aarray[i];
}

FIRFilter::~FIRFilter() {
 	delete [] a;	
}

double FIRFilter::getvalue(SignalSource &x,SignalSource &y,long pos) {
	double acc = 0.0;
	
	for(int i=0;i<length;i++)
		acc += a[i] * x[pos-i];

	return acc;
}


/* Note that IIR makes copies of the coefficient arrays. Normally filters
   are defined so that a[0] = 1 and b[0] = 0; b[0] is not used, anyway */

IIRFilter::IIRFilter(double *aarray,int lena,double *barray,int lenb)
	: a(new double[lena]), b(new double[lenb]), lengtha(lena), lengthb(lenb) {

	for(int i=0;i<lena;i++)
		a[i] = aarray[i];

	for(int j=0;j<lenb;j++)
		b[j] = barray[j];		
}

IIRFilter::~IIRFilter() {
	delete [] b;
	delete [] a;
}

double IIRFilter::getvalue(SignalSource &x,SignalSource &y,long pos) {
	double acc = 0.0;
	
	for(int i=0;i<lengtha;i++)
		acc += a[i] * x[pos-i];

	for(int j=1;j<lengthb;j++)
		acc += b[j] * y[pos-j];

	return acc;
}


// --- SignalFilter ---

SignalFilter::SignalFilter(long len,SignalSource *src,Filter *flt)
	: BufferedSignalSource(len), source(src), filter(flt) {}

void SignalFilter::reset(unsigned long seed) {
	source->reset(seed);
	
	BufferedSignalSource::reset(seed);
}

double SignalFilter::getvalue(long pos) {
	return filter->getvalue(*source,*this,pos);
}

// --- Interpolators ---

double NearestInterpolator::getvalue(SignalSource &y,long ind,double dind) {
	return (dind < 0.5 ? y[ind] : y[ind+1]);
}


// 0 < dind < 1; the desired sample would be at ind + dind

double LinearInterpolator::getvalue(SignalSource &y,long ind,double dind) {
	return (1.0 - dind) * y[ind] + dind * y[ind+1];
}


// this is to use only "old" values, with (implicitly) 1 < dind < 2

double LinearExtrapolator::getvalue(SignalSource &y,long ind,double dind) {
	return (-dind) * y[ind-1] + (1.0 + dind) * y[ind];
}


// --- DotLagrangeInterpolator ---

double DotLagrangeInterpolator::getvalue(SignalSource &y,long ind,double dind) {
	// ya[1] .. ya[2*semiwindow]
	// xa[1] .. xa[2*semiwindow] = 1 .. 2*semiwindow

	for(int i=0;i<semiwindow;i++) {
		ya[semiwindow-i] = y[ind-i];
		ya[semiwindow+i+1] = y[ind+i+1];
	}

	for(int i=1;i<=window;i++) {
		yd[i] = 0.0;
	};

	// x = semiwindow + dind is between ya[semiwindow] and ya[semiwindow+1]

	return dpolint(semiwindow+dind);
}

double DotLagrangeInterpolator::dpolint(double x) {
	for(int gen=1;gen<=window-1;gen++) {
		for(int i=1;i<=window-gen;i++) {
			double &xb = xa[i+gen];
			double &xt = xa[i];

			double invden = 1.0 / (xt - xb);

			yd[i] = (ya[i] - ya[i+1] + (x - xb) * yd[i] + (xt - x) * yd[i+1]) * invden;
			ya[i] = ((x - xb) * ya[i] + (xt - x) * ya[i+1]) * invden;
		}
	}

	return yd[1];
}

DotLagrangeInterpolator::DotLagrangeInterpolator(int semiwin)
    : window(2*semiwin), semiwindow(semiwin),
      xa(new double[2*semiwin+1]), ya(new double[2*semiwin+1]),
      yd(new double[2*semiwin+1]) {
                
	for(int i=1;i<=window;i++) {
		xa[i] = 1.0*i;

		ya[i] = 0.0;
		yd[i] = 0.0;
	}    
}

DotLagrangeInterpolator::~DotLagrangeInterpolator() {
	delete [] yd;
    delete [] ya;
    delete [] xa;
}

double LagrangeInterpolator::getvalue(SignalSource &y,long ind,double dind) {
	for(int i=0;i<semiwindow;i++) {
		ya[semiwindow-i] = y[ind-i];
		ya[semiwindow+i+1] = y[ind+i+1];
	}

	return polint(semiwindow+dind);
}

LagrangeInterpolator::LagrangeInterpolator(int semiwin)
    : window(2*semiwin), semiwindow(semiwin),
      xa(new double[2*semiwin+1]), ya(new double[2*semiwin+1]),
      c(new double[2*semiwin+1]),  d(new double[2*semiwin+1]) {
                
	for(int i=1;i<=window;i++) {
		xa[i] = 1.0*i;
		ya[i] = 0.0;
	}    
}

LagrangeInterpolator::~LagrangeInterpolator() {
    delete [] d;
    delete [] c;
        
    delete [] ya;
    delete [] xa;
}

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


NewLagrangeInterpolator::NewLagrangeInterpolator(int semiwin)
    : window(2*semiwin), semiwindow(1.0*semiwin),
      xa(new double[2*semiwin+1]), ya(new double[2*semiwin+1]),
      c(new double[2*semiwin+1]),  d(new double[2*semiwin+1]) {
		
	for(int i=1;i<=window;i++) {
		xa[i] = 1.0*i;
		ya[i] = -1.0/xa[i];
	}    
}

NewLagrangeInterpolator::~NewLagrangeInterpolator() {
    delete [] d;
    delete [] c;
	
    delete [] ya;
    delete [] xa;
}

double NewLagrangeInterpolator::getvalue(SignalSource &y,long ind,double dind) {
	int bind = ind - window/2;

	for(int i=window;i>0;i--) {
		c[i] = d[i] = y[bind + i];
	}

	/* was:	for(int i=0;i<semiwindow;i++) {
    			ya[semiwindow-i] = y[ind-i];
    			ya[semiwindow+i+1] = y[ind+i+1];
			}
	   with ya[] copied into c[] and d[] in polint */

	return polint(semiwindow+dind);
}

double NewLagrangeInterpolator::polint(double x) {
    int n = window, ns = 1;

    double dif, mindif = fabs(x-xa[1]);

    for (int i=2;i<=n;i++) {
		dif = fabs(x-xa[i]);
		
		if(dif < mindif) {
			ns = i;
			mindif = dif;
		}
	}

    double den, res = c[ns--];

    for (int m=1;m<n;m++) {
		for (int i=1;i<=n-m;i++) {		
			den = ya[m] * (c[i+1] - d[i]);

			c[i] = (xa[i]   - x) * den;
			d[i] = (xa[i+m] - x) * den;
		}
	
		// the summand here is the error estimate
		res += (2*ns < (n-m) ? c[ns+1] : d[ns--]);
    }

    return res;
}


// getInterpolator

Interpolator *getInterpolator(int interplen) {
	if (interplen == 0)
    	return new NearestInterpolator();
	else if (interplen == -1)
		return new LinearExtrapolator();
	else if (interplen == 1)
		return new LinearInterpolator();
	else if (interplen > 0)
		return new LagrangeInterpolator(interplen);
	else if (interplen < -1)
	    return getDerivativeInterpolator(-interplen);
	else {
		std::cerr << "getInterpolator(...): undefined interpolator length "
		          << interplen << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;
	
		ExceptionUndefined e;
		throw e;
	}
}


Interpolator *getDerivativeInterpolator(int interplen) {
	if (interplen > 0)
		return new DotLagrangeInterpolator(interplen);
	else {
		std::cerr << "getDerivativeInterpolator(...): undefined interpolator length "
		          << interplen << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;
	
		ExceptionUndefined e;
		throw e;
	}
}


// InterpolatedSignal

InterpolatedSignal::InterpolatedSignal(SignalSource *src,Interpolator *inte,
									   double deltat,double prebuffer,double norm)
	: source(src), interp(inte),
	  samplingtime(deltat), prebuffertime(prebuffer), normalize(norm) {}

void InterpolatedSignal::reset(unsigned long seed) {
	source->reset(seed);
}

double InterpolatedSignal::value(double time) {
	if (normalize == 0.0) return 0.0;

	try {
		double ireal, iint, ifrac;

		ireal = (time + prebuffertime) / samplingtime;
		iint  = floor(ireal);
		ifrac = ireal - iint;

		return normalize * interp->getvalue(*source,long(iint),ifrac);
	} catch (ExceptionOutOfBounds &e) {
		std::cerr << "InterpolateSignal::value(double) : OutOfBounds while accessing "
		          << time << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;
		
		throw e;
	}
}

// do everything much more carefully (and slowly) to avoid rounding error
// there may be a more efficient way to do this

double InterpolatedSignal::value(double timebase,double timecorr) {
	if (normalize == 0.0) return 0.0;

	try {
		double irealb, iintb, ifracb;
		double irealc, iintc, ifracc;
		double ifrac;

		// irealb = (timebase + prebuffertime) / samplingtime;
		irealb = timebase / samplingtime;
		iintb  = floor(irealb);
		ifracb = irealb - iintb;

		// irealc = timecorr / samplingtime;
		irealc = (timecorr + prebuffertime) / samplingtime;
		iintc  = floor(irealc);
		ifracc = irealc - iintc;

		ifrac = ifracb + ifracc;

		if (ifrac >= 1.0) {
			return normalize * interp->getvalue(*source,long(iintb+iintc)+1,ifrac-1.0);
		} else {
			return normalize * interp->getvalue(*source,long(iintb+iintc),ifrac);
		}
	} catch (ExceptionOutOfBounds &e) {
		std::cerr << "InterpolateSignal::value(double,double): OutOfBounds while accessing "
		          << "(" << timebase << "," << timecorr << ")"
				  << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;
		
		throw e;
	}
}

void InterpolatedSignal::setinterp(Interpolator *inte) {
	interp = inte;
}


// PowerLawNoise

/* ??? I wonder what the 32 is doing in "prebuffer/deltat+32". Does it imply
   a maximum length for interpolators? Also in SampledSignal ??? */

PowerLawNoise::PowerLawNoise(double deltat,double prebuffer,
				double psd,double exponent,int interplen,unsigned long seed) {

	double nyquistf = 0.5 / deltat;
	double normalize;
	
	if (exponent == 0.00) {
		filter = new NoFilter();
		normalize = sqrt(psd) * sqrt(nyquistf);
	} else if (exponent == 2.00) {
		filter = new DiffFilter();
		normalize = sqrt(psd) * sqrt(nyquistf) / (2.00 * M_PI * deltat);
	} else if (exponent == -2.00) {
		filter = new IntFilter();
        normalize = sqrt(psd) * sqrt(nyquistf) * (2.00 * M_PI * deltat);
	} else {
		std::cerr << "PowerLawNoise::PowerLawNoise(...): undefined PowerLaw exponent "
		          << exponent << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;

		ExceptionUndefined e;
		throw e;
	}

	whitenoise = new WhiteNoiseSource(long(prebuffer/deltat+32),seed);
	filterednoise = new SignalFilter(long(prebuffer/deltat+32),whitenoise,filter);

	try {
		interp = getInterpolator(interplen);	
	} catch (ExceptionUndefined &e) {
		std::cerr << "PowerLawNoise::PowerLawNoise(...): undefined interpolator length "
		          << interplen << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;

		throw e;		
	}

	interpolatednoise = new InterpolatedSignal(filterednoise,interp,deltat,prebuffer,normalize);
}

PowerLawNoise::~PowerLawNoise() {
	delete interpolatednoise;
	delete interp;
	delete filterednoise;
	delete whitenoise;
	delete filter;
}

void PowerLawNoise::reset(unsigned long seed) {
	interpolatednoise->reset(seed);
}


// SampledSignal

SampledSignal::SampledSignal(double *narray,long length,double deltat,double prebuffer,
	double norm,Filter *filter,int interplen) {
	
	try {
		interp = getInterpolator(interplen);	
	} catch (ExceptionUndefined &e) {
		std::cerr << "SampledSignal::SampledSignal(...): undefined interpolator length "
		          << interplen << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;

		throw e;
	}

	/* if (interplen > prebuffer/deltat) {
		std::cerr << "WARNING: SampledSignal::SampledSignal(...): for t = 0, interpolator (semiwin=" 
				  << interplen << ") will stray beyond prebuffer, yielding zeros." << std::endl;
	
	} */

	samplednoise = new SampledSignalSource(narray,length,norm);

	if (!filter) {
		filteredsamples = 0;
		interpolatednoise = new InterpolatedSignal(samplednoise,interp,deltat,prebuffer);		
	} else {
		filteredsamples = new SignalFilter(long(160.0/deltat+32.0),samplednoise,filter);
		interpolatednoise = new InterpolatedSignal(filteredsamples,interp,deltat,prebuffer);
	}
}

SampledSignal::~SampledSignal() {
	delete interpolatednoise;
	delete filteredsamples;
	delete samplednoise;
	delete interp;
}

// Derive this from InterpolatedSignal, or simply contain it?

CachedSignal::CachedSignal(Signal *signal,long length,double deltat,int interplen) {
	try {
		interp = getInterpolator(interplen);	
	} catch (ExceptionUndefined &e) {
		std::cerr << "CachedSignal::CachedSignal(...): undefined interpolator length "
				  << interplen << " [" << __FILE__ << ":" << __LINE__ << "]." << std::endl;
 
		throw e;
	}

	double prebuffer = interplen * deltat;

	resample = new ResampledSignalSource(length,deltat,prebuffer,signal);
	interpsignal = new InterpolatedSignal(resample,interp,deltat,prebuffer);
}

CachedSignal::~CachedSignal() {
	delete interpsignal;
	delete resample;
}

void CachedSignal::reset(unsigned long seed) {
	// InterpolatedSignal will reset also the SignalSource

	interpsignal->reset();
}











