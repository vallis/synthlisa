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

// --- correlated-Gaussian classes ------------------------------------------------

ExpGaussNoise::ExpGaussNoise(double samplinginterval, double lapseinterval, double foldingtime, double spectraldensity) {
    if (idum == 0) idum = -1;

    samplingtime = samplinginterval;
    lapsetime = lapseinterval;

    // the buffersize should be set from the sampling time and the number
    // of in-between samples needed in TDI; say 32?

    buffersize = 32 * long(floor(lapsetime / samplingtime) + 1);

    // e-folding time

    lambda = -1.00 * foldingtime;

    // need a sqrt of the Nyquist frequency here? probably yes

    normalize = sqrt(spectraldensity) *  sqrt(1.00 / (2*samplingtime));

    // now allocate the buffer of pointers

    ptbuffer = new (GaussSample *)[buffersize];
    
    for(int i=0;i<buffersize;i++)
        ptbuffer[i] = new GaussSample;
        
    bufferlevel = buffersize - 1;
    
    // need to create first and last sample, at -lapsetime
    
    first = ptbuffer[bufferlevel];
    bufferlevel = bufferlevel - 1;
    last = first;
    
    first->time = -lapsetime;
    first->value = gasdev();
    
    first->prev = 0;
    first->next = 0;
}

ExpGaussNoise::~ExpGaussNoise() {
    // ah, need to be careful here; first release the buffer, then the list

    GaussSample *nextone;

    while(first) {
        nextone = first->next;
        delete first;
        first = nextone;
    }

    while(bufferlevel >= 0) {
        delete ptbuffer[bufferlevel];
        bufferlevel = bufferlevel - 1;
    }
        
    delete ptbuffer;
}

double ExpGaussNoise::operator[](double time) {
    GaussSample *current = first;

    while(current->time <= time) {
        if(current->time == time) {
            // we already have the sample
    
            return(normalize * current->value);
        }
        
        current = current->next;
        if(!current) break;
    };
    
    // we'll need a new sample
    
    GaussSample *newone = ptbuffer[bufferlevel];
    ptbuffer[bufferlevel] = 0;
    
    bufferlevel = bufferlevel - 1;
    if(bufferlevel < 0) {
        cout << "ExpGaussNoise::out of GaussSample instances in buffer" << endl;
        abort();  
    }
    
    // now, let's check if we're at the end of the chain, at the beginning, or between samples
    
    if(!current) {
        // we're at the end of the chain!
        // create a new sample using the "last" rule
        // remember that lambda is negative

        double r = exp(lambda * (time - last->time));

        newone->time = time;  
        newone->value = sqrt(1 - r*r) * gasdev() + r * (last->value);
    
        // adjust the links in the chain
    
        newone->prev = last;
        newone->next = 0;
        
        last->next = newone;
        last = newone;
        
        // since we added at the end, purge from the beginning

        while(time - first->time > lapsetime) {
            bufferlevel = bufferlevel + 1;
            ptbuffer[bufferlevel] = first;
            
            first = first->next;
            first->prev = 0;
        }
    } else {
        GaussSample *previous = current->prev;

        if(!previous) {
            // we're at the beginning of the chain!

            double r = exp(lambda * (current->time - time));

            newone->time = time;  
            newone->value = sqrt(1 - r*r) * gasdev() + r * (current->value);

            newone->prev = 0;
            newone->next = last;

            first->prev = newone;
            first = newone;
        } else {
            // we're in between! "current" holds the larger time, "previous" the smaller
            // create a new sample using the "bridge" rule
            // remember that lambda is negative
        
            double r1 = exp(lambda * (time - previous->time));
            double r2 = exp(lambda * (current->time - time));

            newone->time = time;
            newone->value = sqrt( (1 - r1*r1) * (1 - r2*r2) / (1 - r1*r1*r2*r2) ) * gasdev() + 
                                  (1 - r2*r2) * r1 * (previous->value) / (1 - r1*r1*r2*r2) +
                                  (1 - r1*r1) * r2 * (current->value)  / (1 - r1*r1*r2*r2);

            // adjust the links in the chain

            newone->prev = previous;
            newone->next = current;
        
            previous->next = newone;
            current->prev = newone;
        }
    }
    
    return(normalize * newone->value);
}

