#ifndef _LISASIM_NOISE_H_
#define _LISASIM_NOISE_H_

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

// global seed for random number generator

extern int idum;

// current implementation limits the total number of total samples to MAX_LONG

class RingNoise : public BufferNoise {
    public:
        long buffersize;
        long earliest, latest;
    
        double *bufferx, *buffery;
    
        void updatebuffer(long pos);

        RingNoise(long bs);
        virtual ~RingNoise();
        
        void reset();
        
        double operator[](long pos);
    
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


// correlated-Gaussian class

class GaussSample {
    public:
        GaussSample *prev, *next;
        
        double time;
        double value;
};

class ExpGaussNoise : public Noise {
    private:
        GaussSample **ptbuffer;
        int buffersize;
        int bufferlevel;
    
        GaussSample *first, *last;
        
        double samplingtime;
        double lapsetime;
        
        double lambda;
        double normalize;
    
    public:
        ExpGaussNoise(double samplinginterval, double lapseinterval, double foldingtime, double spectraldensity);
        ~ExpGaussNoise();
    
        void reset();
    
        double enoise(double time);
        double operator[](double time);
};

#endif /* _LISASIM_NOISE_H_ */
