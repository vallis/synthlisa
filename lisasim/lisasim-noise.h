#ifndef _LISASIM_NOISE_H_
#define _LISASIM_NOISE_H_

// global seed for random number generator

extern int idum;

// current implementation limits the total number of total samples to MAX_LONG

class RingNoise {
    public:
    
        long buffersize;
        long earliest, latest;
    
        double *bufferx, *buffery;
    
        RingNoise(long bs);
        virtual ~RingNoise();
        
        void reset();
        
        void updatebuffer(long pos);
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
    public:

        double intconst;
        
        IntNoise(long bs, double ic);
        
        double filter(long pos);
};

// dimensioned interpolated noise

class InterpolateNoise {
    public:

        RingNoise *buffernoise;

        // use time in seconds

        double samplingtime;
        double nyquistf;

        long prebuffertime;
        double maxitime;
        
        // sqrt spectral density in Hz^-1/2
        
        double normalize;

        InterpolateNoise(double st, double pbt, double sd, double ex);
        ~InterpolateNoise();

        void reset();
        
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

class ExpGaussNoise {
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
