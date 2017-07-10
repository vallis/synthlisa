/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#ifndef _LISASIM_WAVE_H_
#define _LISASIM_WAVE_H_

#include "lisasim-tens.h"
#include "lisasim-signal.h"

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


// --- Wave ---

class Wave : public WaveObject {
 public:
    // position in the sky
    // beta is the SSB ecliptic latitude
    // lambda is the SSB ecliptic longitude from the vernal point

    double beta, lambda, pol;

    // k vector

    Vector k;

    // polarization tensors

    Tensor pp, pc;

    Wave(double b, double l, double p);

    Wave *firstwave() { return this; }
    Wave *nextwave()  { return 0; }

    virtual int inscope(double t) { return 1; }

    virtual double hp(double t) = 0;
    virtual double hc(double t) = 0;

    void putk(Vector &k);
    void putwave(Tensor &h, double t);

    static void putep(Tensor &h,double b,double l,double p);
    static void putec(Tensor &h,double b,double l,double p);
};


// --- SimpleBinary ---

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


// --- GalacticBinary ---

class GalacticBinary : public Wave {
 private:
	// frequency, initial phase

	double f, fdot, fddot, eps, phi0;

	// inclination, polarization amplitudes

	double i, a, ap, ac;

 public:
	GalacticBinary(double freq, double freqdot, double b, double l, double amp, double inc, double p, double initphi, double fddot = 0.0, double epsilon = 0.0);

	double hp(double t);
	double hc(double t);
};


// --- SimpleMonochromatic ---

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


// --- GaussianPulse ---

class GaussianPulse : public Wave {
 private:
    double t0, dc; // offset time, decay
    double gm, a;  // polarization and amplitude
    double ap, ac; // polarization amplitudes

    static const double sigma_cutoff; // defined in lisasim-wave.cpp

 public:
    GaussianPulse(double time, double decay, double gamma, double amp, double b, double l, double p);

    int inscope(double t);

    double hp(double t);
    double hc(double t);
};


// --- SineGaussian ---

class SineGaussian : public Wave {
 private:
    double t0, dc;  // offset time, decay
	double f, phi0; // base frequency and phase
    double gm, a;   // polarization and amplitude
    double ap, ac;  // polarization amplitudes

    static const double sigma_cutoff; // defined in lisasim-wave.cpp

 public:
    SineGaussian(double time, double decay, double freq, double phase0, double gamma, double amp, double b, double l, double p);

    int inscope(double t);

    double hp(double t);
    double hc(double t);
};


// --- NoiseWave ---

class NoiseWave : public Wave {
 private:
    Noise *np, *nc;

	// set this to one if we are allocating noise objects

	int allocated;

 public:
	NoiseWave(Noise *noisehp, Noise *noisehc, double b, double l, double p);
	NoiseWave(double sampletime, double prebuffer, double density, double exponent, int swindow, double b, double l, double p);
	NoiseWave(double *hpa, double *hca, long samples, double sampletime, double prebuffer, double norm, Filter *filter, int swindow, double b, double l, double p);

	~NoiseWave();

	double hp(double t) { return np->noise(t); };
	double hc(double t) { return nc->noise(t); };
};


// --- SampledWave ---

/* This is really a frontend for NoiseWave, so we're just providing
   this factory function. The Python interface has its own factory,
   written in Python, which deals with garbage collection issues,
   etc. Perhaps there's a better way to do this in C++... */

NoiseWave *SampledWave(double *hpa, double *hca, long samples, double sampletime, double prebuffer, double density, Filter *filter, int swindow, double d, double a, double p);

// --- PyWave ---

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
