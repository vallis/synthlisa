#include <iostream.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include "lisasim.h"
#include "lisasim-noise.h"
#include "lisasim-tdinoise.h"

// outputs noise on X, Y, Z observables

int main(int argc, char **argv) {
    // take as only input the number of noise samples to output

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " [samples/year] " << endl;
        exit(0);
    }

    int samples = atoi(argv[1]);
    
    // arbitrary upper limit
    
    if (samples < 0 || samples > 10e6) {
        cout << "Number of samples unacceptable!" << endl;
        exit(0);
    }

    CircularRotating mylisa(0.0,3.0*M_PI/2.0,-1.0);

    // use this is the nominal time interval of the noises

    double dstep = 1.0;

    TDInoise mytdinoise(&mylisa,dstep,2.5e-48,dstep,1.8e-37);

    // initialize random-number-generator seed
    // here we're actually passing a long (on PPC); should work as long
    // is the same as int

    idum = -time(0);

    ofstream xout("noise-X.txt");
//    ofstream yout("noise-Y.txt");
//    ofstream zout("noise-Z.txt");

    ofstream aout("noise-alpha.txt");

    ofstream pout("noise-P.txt");
    ofstream eout("noise-E.txt");
    
    // set this to a fraction of dstep to oversample by interpolation
    // (see the results on the spectra!)

    double sstep = 1.0;

    for(int i=0;i<samples;i++) {
        xout << mytdinoise.X(sstep * i) << endl;
//        yout << mytdinoise.Y(sstep * i) << endl;
//        zout << mytdinoise.Z(sstep * i) << endl;

        aout << mytdinoise.alpha(sstep * i) << endl;

        pout << mytdinoise.P(sstep * i) << endl;
        eout << mytdinoise.E(sstep * i) << endl;
    }
}
