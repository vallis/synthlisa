#include "lisasim.h"
#include "lisasim-noise.h"

#include "nr.h"
#include <cmath>
#include <stdlib.h>

// --- uniform random number generator --------------------------------------------------

// from numerical recipes is C++:
// "Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added
// safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
// values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
// successive deviates in a sequence. RNMX should approximate the largest floating value that is
// less than 1.

// modified by changing idum to a global variable

int idum = 0;

DP ran1()
{
        const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
        const int NDIV=(1+(IM-1)/NTAB);
        const DP EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
        static int iy=0;
        static Vec_INT iv(NTAB);
        int j,k;
        DP temp;

        if (idum <= 0 || !iy) {
                if (-idum < 1) idum=1;
                else idum = -idum;
                for (j=NTAB+7;j>=0;j--) {
                        k=idum/IQ;
                        idum=IA*(idum-k*IQ)-IR*k;
                        if (idum < 0) idum += IM;
                        if (j < NTAB) iv[j] = idum;
                }
                iy=iv[0];
        }
        k=idum/IQ;
        idum=IA*(idum-k*IQ)-IR*k;
        if (idum < 0) idum += IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j] = idum;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}

// from numerical recipes in C++:
// Returns a normally distributed deviate with zero mean and unit variance,
// using ran1(idum) as the source of uniform deviates.

DP gasdev()
{
        static int iset=0;
        static DP gset;
        DP fac,rsq,v1,v2;

        if (iset == 0) {
                do {
                        v1=2.0*ran1()-1.0;
                        v2=2.0*ran1()-1.0;
                        rsq=v1*v1+v2*v2;
                } while (rsq >= 1.0 || rsq == 0.0);
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

// --- "Ring" classes --------------------------------------------------

//  default constructor for RingNoise
//  initialize ran1 if needed; should do with clock or with seed provided from main

RingNoise::RingNoise(long bs) {
    if (idum == 0) idum = -1;

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

double RingNoise::deviate() {
    return gasdev();
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

// constructor: need to pass the sampling time (in secs), the requested prebuffer time
// (probably a multiple of some nominal armlength), and the one-sided power spectral density,
// given as sd*(f/Hz)^ex, where sd is in Hz^-1, and ex is either 0.00, -2.00, or 2.00.
// The corresponding ringnoise object is then allocated automatically

InterpolateNoise::InterpolateNoise(double st, double pbt, double sd, double ex) {
    samplingtime = st;
    nyquistf = 1.00 / (2*st);

    // compute the prebuffer offset

    prebuffertime = long(floor(pbt / samplingtime)) + 1;

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

double InterpolateNoise::inoise(double time) {
    double ctime = time / samplingtime;
    double itime = floor(ctime);

    long index;
    
    if (itime > maxitime) {
        cout << "InterpolateNoise::inoise() accessing index larger than LONG_MAX" << endl;
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
