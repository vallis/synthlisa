#ifndef _LISASIM_NOISE_H_
#define _LISASIM_NOISE_H_

#include "gsl_rng.h"

class BufferNoise {
 public:
    virtual void reset() {};

    virtual double operator[](long pos) = 0;
};

class SampledNoise : public BufferNoise {
 private:
    double *noisebuffer;

    long maxsamples;

 public:
    SampledNoise(double *nb, long samples);

    double operator[](long pos);
};

// current implementation limits the total number of total samples to MAX_LONG

class RingNoise : public BufferNoise {
 private:
    gsl_rng *randgen;

    int cacheset;
    double cacherand;

 public:
    long buffersize;
    long earliest, latest;
    
    double *bufferx, *buffery;
    
    void updatebuffer(long pos);

    RingNoise(long bs);
    virtual ~RingNoise();
	
    void reset();
        
    double operator[](long pos);

    void seedrandgen();
    
    virtual double deviate();
    virtual double filter(long pos);
};

// derived classes: differentiating and integrating filters

class DiffNoise : public RingNoise {
    public:
        DiffNoise(long bs);

        double filter(long pos);
};

class IntNoise : public RingNoise {
    private:
        double intconst;

    public:
        IntNoise(long bs, double ic);
        
        double filter(long pos);
};

// the following are the dimensioned noises that are used in the TDInoise class
// we make them descendants of a generic Noise class that defines only the [] override
// 

class Noise {
    public:
        Noise() {};
        virtual ~Noise() {};

        virtual void reset() {};

        virtual double operator[](double time) = 0;
};

// dimensioned interpolated noise

class InterpolateNoise : public Noise {
    protected:
        BufferNoise *buffernoise;

        // use time in seconds

        double samplingtime;
        double nyquistf;

        long prebuffertime;
        double maxitime;
        
        // sqrt spectral density in Hz^-1/2
        
        double normalize;

    public:
        InterpolateNoise(double sampletime,double prebuffer,double density,double exponent);
	InterpolateNoise(double *noisebuf,long samples,double sampletime,double prebuffer,double norm);
        ~InterpolateNoise();

        void reset();
        
        virtual double inoise(double time);
        virtual double operator[](double time);
};

class InterpolateNoiseBetter : public InterpolateNoise {
 private:
    double sinc(double time);
    double window(double time);

    int semiwindow;
    double semiconst;

 public:
    InterpolateNoiseBetter(double sampletime,double prebuffer,double density,double exponent,int swindow);
    InterpolateNoiseBetter(double *noisebuf,long samples,double sampletime,double prebuffer,double norm,int swindow);

    double inoise(double time);
    double operator[](double time);
};

const int maxwindow = 100;

class InterpolateNoiseLagrange : public InterpolateNoise {
 private:
    double polint(double x, int n);

    int semiwindow;

    double xa[2*maxwindow+1],ya[2*maxwindow+1];
    double c[2*maxwindow+1],d[2*maxwindow+1];

 public:
    InterpolateNoiseLagrange(double sampletime,double prebuffer,double density,double exponent,int swindow);
    InterpolateNoiseLagrange(double *noisebuf,long samples,double sampletime,double prebuffer,double norm,int swindow);

    double lagnoise(double time,int semiorder);

    double inoise(double time);
    double operator[](double time);

    double pnoise(double time);
};

#endif /* _LISASIM_NOISE_H_ */
