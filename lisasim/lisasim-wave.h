#ifndef _LISASIM_WAVE_H_
#define _LISASIM_WAVE_H_

#include "lisasim-tens.h"

class Wave {
    public:
        // position in the sky

        double dec, asc, pol;

        // k vector

        Vector k;
	double kArray[3];

        // polarization tensors
    
        Tensor pp, pc;
	double ppArray[9], pcArray[9];

        Wave(double d, double a, double p);
        virtual ~Wave() {};

        virtual double hp(double t) = 0;
        virtual double hc(double t) = 0;  

        void putwave(Tensor &h, double t);
        void putwave(double **h, double t);
};

class SimpleBinary : public Wave {
    private:
        // frequency, initial phase

        double f, phi0;

        // inclination, polarization amplitudes

        double i, a, ap, ac;

    public:
        SimpleBinary(double freq, double initphi, double inc, double amp, double d, double a, double p);

        double hp(double t);
        double hc(double t);
};

class SimpleMonochromatic : public Wave {
    private:
        // frequency

        double f;

        // polarization angles using John's convention

        double gm, ph, ap, ac;

    public:
        SimpleMonochromatic(double freq, double phi, double gamma, double amp, double d, double a, double p);

        double hp(double t);
        double hc(double t);
};

class InterpolateMemory : public Wave {
    private:
        double *hpbuffer, *hcbuffer;

        long maxsamples;
        double sampletime;
        double lkback;

    public:
        InterpolateMemory(double *hpa, double *hca, long samples, double samplingtime, double lookback, double d, double a, double p);

        double hp(double t);
        double hc(double t);
};

#endif /* _LISASIM_WAVE_H_ */
