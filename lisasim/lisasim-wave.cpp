#include "lisasim-tens.h"
#include "lisasim-wave.h"

Wave::Wave(double d, double a, double p) {
    dec = d;
    asc = a;
    pol = p;

    kArray[0] = -cos(asc)*cos(dec);
    kArray[1] = -sin(asc)*cos(dec);
    kArray[2] = -sin(dec);

    k[0] = kArray[0];
    k[1] = kArray[1];
    k[2] = kArray[2];

    Tensor stdpp, stdpc;
    
    stdpp[0][0] = 1.0; stdpp[1][1] = -1.0;
    stdpc[0][1] = 1.0; stdpc[1][0] =  1.0;

    Tensor A, At;
    
    A.seteuler(dec,asc,pol);
    At = A; 
    At.seteuler(dec,asc,pol);
    At.settranspose(); 
    
    Tensor tmp;
    
    tmp.setproduct(stdpp,At);
    pp.setproduct(A,tmp);
    
    tmp.setproduct(stdpc,At);
    pc.setproduct(A,tmp);      

    for(int i=0;i<3;i++) {
      for(int j=0;j<3;j++) {
	ppArray[3*i+j] = pp[i][j];
	pcArray[3*i+j] = pc[i][j];
      }
    }
}

void Wave::putwave(Tensor &h, double t) {
  double hp_temp;
  double hc_temp; 
  hp_temp = hp(t);
  hc_temp = hc(t);
  
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      h[i][j] = hp_temp * pp[i][j] + hc_temp * pc[i][j];
    }        
  }
}

void Wave::putwave(double **h, double t) {
  double hp_temp;
  double hc_temp; 
  int i,j;
  hp_temp = hp(t);
  hc_temp = hc(t);
  
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      h[i][j] = hp_temp * ppArray[3*i+j] + hc_temp * pcArray[3*i+j];
    }        
  }
}

// full constructor for SimpleBinary; takes frequency in Hertz
// d and a are (notwithstanding their name, which should be changed) heliocentric ecliptic latitude and longitude

SimpleBinary::SimpleBinary(double freq, double initphi, double inc, double amp, double d, double a, double p) : Wave(d,a,p) {
    f = freq;

    phi0 = initphi;

    i = inc;
    a = amp;
    
    ap = a * (1.0 + cos(i)*cos(i));
    ac = a * (2.0 * cos(i));
}

double SimpleBinary::hp(double t) {
    const double twopi = 2.0*M_PI;

    return ap * cos(twopi*f*t + phi0);
}

double SimpleBinary::hc(double t) {
    const double twopi = 2.0*M_PI;
    
    return ac * sin(twopi*f*t + phi0);
}            


// --- SimpleMonochromatic wave class --------------------------------------------------

// originally written to compare with John's fortran code

SimpleMonochromatic::SimpleMonochromatic(double freq, double phi, double gamma, double amp, double d, double a, double p) : Wave(d,a,p) {
    f = freq;

    ph = phi;
    gm = gamma;
    
    ap = amp*sin(gm);
    ac = amp*cos(gm);
}

double SimpleMonochromatic::hp(double t) {
    const double twopi = 2.0*M_PI;

    return ap * sin(twopi*f*t + ph);
}

double SimpleMonochromatic::hc(double t) {
    const double twopi = 2.0*M_PI;
    
    return ac * sin(twopi*f*t);
}

// --- InterpolateMemory wave class --------------------------------------------------

InterpolateMemory::InterpolateMemory(double *hpa, double *hca, long samples, double samplingtime, double lookback, double d, double a, double p) : Wave(d,a,p) {
  hpbuffer = hpa;
  hcbuffer = hca;

  maxsamples = samples;
  sampletime = samplingtime;

  lkback = lookback;
}

double InterpolateMemory::hp(double time) {
    double ctime = (time + lkback) / sampletime;
    double itime = floor(ctime);
    long index = long(itime);
    
    if (index > maxsamples) {
      cout << "InterpolateMemory::hp() accessing index larger than available" << endl;
      abort();
    } else if (index < 0) {
      cout << "InterpolateMemory::lookback buffer insufficient" << endl;
      abort();
    }

    return ( (ctime - itime) * hpbuffer[index+1] + (1.00 - (ctime - itime)) * hpbuffer[index] );
}

double InterpolateMemory::hc(double time) {
    double ctime = (time + lkback) / sampletime;
    double itime = floor(ctime);
    long index = long(itime);
    
    if (index > maxsamples) {
      cout << "InterpolateMemory::hp() accessing index larger than available" << endl;
      abort();
    } else if (index < 0) {
      cout << "InterpolateMemory::lookback buffer insufficient" << endl;
      abort();
    }

    return ( (ctime - itime) * hcbuffer[index+1] + (1.00 - (ctime - itime)) * hcbuffer[index] );
}
