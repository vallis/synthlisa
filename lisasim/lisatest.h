#ifndef _LISATEST_H_
#define _LISATEST_H_

// All taken from the standard Montana LISAconstants.h

// Physical constants: c (m/s), G (mks), AU (m), yr (s), Msun (kg), kpc (m), Gpc (m)

#define c 2.99792458e8
#define G 6.67259e-11
#define AU 1.49597870660e11
#define year 3.15581498e7
#define Msun 1.9889e30
#define kpc 3.0856675807e19
#define Gpc 3.0856675807e25

// Power of 2 that gives the total number of sample data points
#define nsource 6.0

// All of these taken from the standard Montana NewtonianParameters.h

// Observation time (seconds) (this does not include padding)
#define ObsTime (1.0*year)

// Source location: polar and azimuthal angles in heliocentric-ecliptic coordinates

// #define theta 0.917311
// #define phi 2.97378
#define theta 1.57
#define phi 0.00

// Basic source parameters: freq (Hz), inclination angle, polarization angle, initial phase

#define fgw (2.0/1028.76)
#define inc 1.53921

// #define psi 2.46531
#define psi 0.00

#define phase 0.0

// parameters that establish the amplitude of the signal: rest masses and luminosity distance

#define M1 (0.5*Msun)
#define M2 (0.033*Msun)
#define r (0.1*kpc)

#endif /*  _LISATEST_H_ */
