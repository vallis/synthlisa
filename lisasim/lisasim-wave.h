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

        // default constructor

        Wave(double d, double a, double p);

        // default destructor; do I need this?

        virtual ~Wave() {};

        virtual double hp(double t) = 0;
        virtual double hc(double t) = 0;  

        void putwave(Tensor &h, double t);
        void putwave(double **h, double t);
};

class SimpleBinary : public Wave {
    public:

        // frequency, initial phase

        double f, phi0;

        // inclination, polarization amplitudes

        double i, a, ap, ac;

        SimpleBinary(double freq, double initphi, double inc, double amp, double d, double a, double p);

        double hp(double t);
        double hc(double t);
};

class SimpleMonochromatic : public Wave {
    public:

        // frequency

        double f;

        // polarization angles using John's convention

        double gm, ph, ap, ac;

        SimpleMonochromatic(double freq, double phi, double gamma, double amp, double d, double a, double p);

        double hp(double t);
        double hc(double t);
};

#endif /* _LISASIM_WAVE_H_ */
