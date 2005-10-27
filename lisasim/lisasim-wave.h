#ifndef _LISASIM_WAVE_H_
#define _LISASIM_WAVE_H_

#include "lisasim-tens.h"
#include "lisasim-noise.h"

class Wave;

class WaveObject {
	public:
		// need a virtual constructor, too?
		virtual ~WaveObject() {};
		
		virtual Wave *firstwave() = 0;
		virtual Wave *nextwave() = 0;
};

class WaveArray : public WaveObject {
	private:
		Wave **wavearray;
		int wavenum;
		int wavecurrent;

	public:
		WaveArray(Wave **warray, int wnum);
		~WaveArray();

		Wave *firstwave(), *nextwave();
};

class Wave : public WaveObject {
    public:
        // position in the sky
        // beta is the SSB ecliptic latitude
        // lambda is the SSB ecliptic longitude from the vernal point

        double beta, lambda, pol;

        // k vector

        Vector k;
	double kArray[3];

        // polarization tensors
    
        Tensor pp, pc;
	double ppArray[9], pcArray[9];

        Wave(double b, double l, double p);

	Wave *firstwave() { return this; }
	Wave *nextwave()  { return 0; }

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
        SimpleBinary(double freq, double initphi, double inc, double amp, double b, double l, double p);

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
        SimpleMonochromatic(double freq, double phi, double gamma, double amp, double b, double l, double p);

        double hp(double t);
        double hc(double t);
};

class GaussianPulse : public Wave {
    private:
       double t0, dc; // offset time, decay
       double gm, a;  // polarization and amplitude
       double ap, ac; // polarization amplitudes

	   static const double sigma_cutoff; // defined in lisasim-wave.cpp

    public:
	   GaussianPulse(double time, double decay, double gamma, double amp, double b, double l, double p);

       double hp(double t);
       double hc(double t);
};

class NoiseWave : public Wave {
    private:
        Noise *np, *nc;

	// set this to one if we are allocating noise objects

	int allocated;

    public:
	NoiseWave(Noise *noisehp, Noise *noisehc, double b, double l, double p);
	NoiseWave(double sampletime, double prebuffer, double density, double exponent, int swindow, double b, double l, double p);
	NoiseWave(double *hpa, double *hca, long samples, double sampletime, double prebuffer, double density, double exponent, int swindow, double b, double l, double p);

	~NoiseWave();

	double hp(double t) { return np->noise(t); };
	double hc(double t) { return nc->noise(t); };
};

class InterpolateMemory : public Wave {
    private:
        double *hpbuffer, *hcbuffer;

        long maxsamples;
        double sampletime;
        double lkback;

    public:
        InterpolateMemory(double *hpa, double *hca, long samples, double samplingtime, double lookback, double b, double l, double p);

        double hp(double t);
        double hc(double t);
};

#include <Python.h>

class PyWave : public Wave {
 private:
    PyObject *hpfunc, *hcfunc;

 public:
    PyWave(PyObject *hpf, PyObject *hcf, double b, double l, double p)
	: Wave(b,l,p), hpfunc(hpf), hcfunc(hcf) {};
    virtual ~PyWave() {};

    double hp(double t) {
	PyObject *arglist, *result;

	double dres = 0.0;
   
	arglist = Py_BuildValue("(d)",t);             // Build argument list
	result = PyEval_CallObject(hpfunc,arglist);  // Call Python
	Py_DECREF(arglist);                           // Trash arglist
	if (result) dres = PyFloat_AsDouble(result);  // If no errors, return double
	Py_XDECREF(result);                           // Trash result
	return dres;
    }

    double hc(double t) {
	PyObject *arglist, *result;

	double dres = 0.0;
   
	arglist = Py_BuildValue("(d)",t);             // Build argument list
	result = PyEval_CallObject(hcfunc,arglist);  // Call Python
	Py_DECREF(arglist);                           // Trash arglist
	if (result) dres = PyFloat_AsDouble(result);  // If no errors, return double
	Py_XDECREF(result);                           // Trash result
	return dres;
    }
};

#endif /* _LISASIM_WAVE_H_ */
