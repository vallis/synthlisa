#include <iostream>
#include <math.h>
#include <unistd.h>
#include <fstream>
#include "lisasim.h"

// a rather longwinded way to unwind the phase...

double unwind(double a, double b, int index) {
    double diff; 
    static double last[3] = {0.0,0.0,0.0};
    
    diff = a - b;
    
    if (diff > 0.0)
        diff -= (2.0*M_PI)*floor(diff/(2.0*M_PI));
    else
        diff += (2.0*M_PI)*floor(-diff/(2.0*M_PI));

    while (diff - last[index] < -M_PI)
        diff += 2.0*M_PI;

    while (diff - last[index] > M_PI)
        diff -= 2.0*M_PI;

    last[index] = diff;

    return(diff);
}

// our main

int main(int argc, char **argv) {
    double fsec = 1E-3;
    double dec = 0.0, asc = 0.0, beta = 0.0, lambda = 0.0, phi0 = 0.0, psi = 0.0, inc = M_PI/2.0;
    int samples = 1024;
    double srate = 1.0;
    int diff = 0, mono = 0, equat = 0, eclip = 0, fast=0, SkipCalc=0;
    double deltat = 0.0;
    char *filename = "tdi.txt";
    TDIsignal *mytdi=0;
    TDIfast *mytdifast=0;
    

    opterr = 0;

    int c;
    
    while ((c = getopt(argc, argv, "+f:d:a:b:l:0:p:i:g:#:r:MBPFhS")) != -1)
        switch (c) {
            case 'f':	// frequency of source
                fsec = atof(optarg);
                break;
            case 'd':	// equatorial declination
                dec = atof(optarg);
                equat = 1;
                break;
            case 'a':   // equatorial right ascension
                asc = atof(optarg);
                equat = 1;
                break;
            case 'b': 	// ecliptic latitude (beta)
                beta = atof(optarg);
                eclip = 1;
                break;
            case 'l': 	// ecliptic longitude (lambda)
                lambda = atof(optarg);
                eclip = 1;
                break;
            case '0':	// initial phase
                phi0 = atof(optarg);
                break;
            case 'p':	// source polarization convention
                psi = atof(optarg);
                break;
            case 'i':	// inclination angle (binaries)
                inc = atof(optarg);
                break;
            case 'g':	// source polarization (monochromatic)
                inc = atof(optarg);
                break;
            case '#':	// number of samples
                samples = atoi(optarg);
                break;
            case 'r':	// sampling interval in seconds
                srate = atof(optarg);
                break;
            case 'P':	// toggle print amplitude + phase
                diff = 1;
                deltat = 0.01 / fsec;
                break;
            case 'M':	// toggle monochromatic source
                mono = 1;
                break;
            case 'B':	// toggle binary
                mono = 0;
                break;
	    case 'F': // Use TDIfast
	        fast = 1;
		break;
	    case 'S': // Skip Calculation (Setup only.)
	        SkipCalc = 1;
		break;
            case 'h':	// help
	      fprintf(stderr, "Usage: %s [options] [outputfile]\n", argv[0]);
		fprintf(stderr, "where options are:\n");
		fprintf(stderr, "\t-f freq(Hz)\n");
		fprintf(stderr, "\t[-d eqdec(rad) -a eqrasc(rad) | -b eclat -l eclon]\n");
		fprintf(stderr, "\t-0 phi0(rad)\n");
		fprintf(stderr, "\t-p psi(rad)\n");
		fprintf(stderr, "\t-i inc(rad)\n");
		fprintf(stderr, "\t-g gamma(rad)\n");
		fprintf(stderr, "\t-# nsamples\n");
		fprintf(stderr, "\t-r samplerate(sec)");
		fprintf(stderr, "The following flags effect processing and output:\n");
		fprintf(stderr, "\t-P\tThe default behavior is to print X, Y, Z.  This option will print amplitude/phase pairs instead.\n");
		fprintf(stderr, "\t-M\tChange the source from simple binary to simple monochromatic.\n");
		fprintf(stderr, "\t-F\tUse TDIfast caching of orbits.  (Note that -P is not supported with this option.)\n");
		fprintf(stderr, "\t-S\tSkip actual calculation -- Setup only\n");
            case '?':	// unrecognized option
                fprintf(stderr, "Unknown option `-%c'.\n", optopt);
		fprintf(stderr, "Use -h to see available options.");
                return 1;
                break;
            default:
                abort();
        }

    // get filename if it is present

    if (optind < argc)
        filename = argv[optind];

    // plot only five cycles if rate is set to 0

    if (srate == 0.0)
        srate = 5.0 / fsec / samples;

    // convert equatorial to ecliptic

    if (equat == 1) {
        if (eclip == 1) {
            fprintf(stderr, "Please choose either ecliptic or equatorial source coordinates.\n");
            return 1;
            }
            
        double eps = 0.408931; // ecliptic inclination
        
        beta = asin(sin(dec)*cos(eps)-cos(dec)*sin(eps)*sin(asc));
        
        if (cos(beta) == 0.00) { // we're at either pole
            lambda = 0.00;
        } else {
            lambda = acos(cos(asc)*cos(dec)/cos(beta));
            if (asin((sin(dec)-sin(beta)*cos(eps))/(cos(beta)*sin(eps))) < 0.0)
                lambda = 2*M_PI-lambda;
        }
    }

    // Output parameters

    cout << "lisasim " << endl;
    cout << "  f = " << fsec << "; lambda = " << lambda << "; beta = " << beta;
    cout << "; phi0 = " << phi0 << "; psi = " << psi << "; inc/gamma = " << inc << endl;
    cout << "  computing " << samples << " samples with a " << srate << " s sampling interval" << endl;


    if (diff == 1) {
      if (fast == 1) {
	cout << "Finite-difference amplitude and phase (-P) is not supported with orbit caching (-F). Aborting." << endl;
	return(1);
      } else {
        cout << "  doing finite-difference amplitude and phase" << endl;
      }
    }

    if (mono == 1)
        cout << "  using simple monochromatic source" << endl;

    cout << "  writing output to " << filename << endl;

    Wave *mywave;
    
    if (mono==1)
        mywave = new SimpleMonochromatic(fsec, phi0, inc, 1.0, beta, lambda, psi);
    else
        mywave = new SimpleBinary(fsec, phi0, inc, 1.0, beta, lambda, psi);

    //    CircularRotating mylisa(0.0,0.0);
    CircularRotating mylisa(1.50*M_PI,0.0,1.0);

    if (fast == 0) {
      mytdi = new TDIsignal(&mylisa,mywave);
    } else {
      mytdifast = new TDIfast(&mylisa,mywave, srate, samples);
      mytdifast->CacheX();
      mytdifast->CacheY();
      mytdifast->CacheZ();
      
      cout << "X, Y, Z initialized" << endl;
    }

    if (!SkipCalc) {
      ofstream fout(filename);
      
      for(int i=0;i<samples;i++) {
	double time = i * srate;
	
	if (diff==1) {
	  double x = mytdi->X(time);
	  double y = mytdi->Y(time);
	  double z = mytdi->Z(time);
	  
	  double dx = (mytdi->X(time+deltat) - x) / deltat / (2*M_PI*fsec);
	  double dy = (mytdi->Y(time+deltat) - y) / deltat / (2*M_PI*fsec);
	  double dz = (mytdi->Z(time+deltat) - z) / deltat / (2*M_PI*fsec);            
	  
	  double cphase = 2*M_PI*fsec*time;
	  
	  fout << time << " " <<
	    sqrt(x*x + dx*dx) << " " << unwind(atan2(-dx,x),cphase,0) << " " <<
	    sqrt(y*y + dy*dy) << " " << unwind(atan2(-dy,y),cphase,1) << " " <<
	    sqrt(z*z + dz*dz) << " " << unwind(atan2(-dz,z),cphase,2) << " " <<  endl;
	} else { 
	  if (fast == 1){
	    fout << time << " " << mytdifast->Xfast(i) << " " << mytdifast->Yfast(i) << " " << mytdifast->Zfast(i) << endl;        
	  } else {
	    fout << time << " " << mytdi->X(time) << " " << mytdi->Y(time) << " " << mytdi->Z(time) << endl;        
	  }
	}  
      }
    }
    cout << "  done!" << endl;
    
    delete mywave;
    if(!fast){
      delete mytdi;
    } else {
      delete mytdifast;
    }
}

