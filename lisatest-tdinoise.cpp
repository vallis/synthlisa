#include <iostream.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include "lisasim.h"
#include "lisasim-lisa.h"
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
    
    if (samples < 0 || samples > 10e8) {
        cout << "Number of samples unacceptable!" << endl;
        exit(0);
    }

    OriginalLISA mylisa(16.6782,16.6782,16.6782);
    ModifiedLISA mymlisa(16.6782,16.6782,16.6782);
                                        
    // use this is the nominal time interval of the noises

    double dstep = 1.0;

    TDInoise mytdinoise(&mylisa,dstep,2.5e-48,dstep,1.8e-37,dstep,1.1e-26,1.0e-6);
    TDInoise mymtdinoise(&mymlisa,dstep,2.5e-48,dstep,1.8e-37,dstep,1.1e-26,1.0e-6);
    
    // initialize random-number-generator seed
    // here we're actually passing a long (on PPC); should work as long
    // as "long" is the same as "int"

    idum = -time(0);

    // set this to a fraction of dstep to oversample by interpolation
    // (see the results on the spectra!)
    // it seems this does not work so well... maybe the normalization of the noise is wrong

    double sstep = 1.0;

    ofstream xout("noise-X.txt");
    ofstream xmout("noise-Xm.txt");
    ofstream xmcout("noise-Xmc.txt");

//  ofstream yout("noise-Y.txt");
//  ofstream zout("noise-Z.txt");

//  ofstream aout("noise-alpha.txt");

//  ofstream pout("noise-P.txt");
//  ofstream eout("noise-E.txt");
    
    ofstream y21out("noise-y21.txt");
    
    for(int i=0;i<samples;i++) {
 	xout << mytdinoise.X(sstep * i)   << endl;
        xmout << mymtdinoise.X(sstep * i)   << endl;
        xmcout << mymtdinoise.Xm(sstep * i)   << endl;        
        y21out << mytdinoise.y(3, 2, 1, 0, 0, 0, sstep*i) << endl;

// 	xmout  << mytdinoise.Xm(sstep * i)  << endl;
// 	xmcout << mymtdinoise.Xm(sstep * i) << endl;
        
//      yout << mytdinoise.Y(sstep * i) << endl;
//      zout << mytdinoise.Z(sstep * i) << endl;

//      aout << mytdinoise.alpha(sstep * i) << endl;

//      pout << mytdinoise.P(sstep * i) << endl;
//      eout << mytdinoise.E(sstep * i) << endl;

//      y21out << mytdinoise.y(3, 2, 1, 0, 0, 0, sstep*i) << endl;
    }
   
/*  MontanaEccentric thlisa(0.0,0.0);
    
    TDInoise thtdinoise(&thlisa,dstep,2.5e-48,dstep,1.8e-37,dstep,1.1e-26,0.1);
    TDInoise thtdinoise2(&thlisa,dstep,2.5e-48,dstep,1.8e-37,dstep,1.1e-26,0.001);
    TDInoise thtdinoise3(&thlisa,dstep,2.5e-48,dstep,1.8e-37,dstep,1.1e-26,0.00001);

    ofstream xunout("noise-X-unequal.txt");
    ofstream xunout2("noise-X-unequal-2.txt");
    ofstream xunout3("noise-X-unequal-3.txt");
        
    ofstream y21unout("noise-y21-unequal.txt");
    ofstream y21unout2("noise-y21-unequal-2.txt");
    ofstream y21unout3("noise-y21-unequal-3.txt");

    for(int i=0;i<samples;i++) {
        xunout << thtdinoise.X(sstep * i) << endl;
        y21unout << thtdinoise.y(3, 2, 1, 0, 0, 0, sstep*i) << endl;
        
        xunout2 << thtdinoise2.X(sstep * i) << endl;
        y21unout2 << thtdinoise2.y(3, 2, 1, 0, 0, 0, sstep*i) << endl;
        
        xunout3 << thtdinoise3.X(sstep * i) << endl;
        y21unout3 << thtdinoise3.y(3, 2, 1, 0, 0, 0, sstep*i) << endl;        
    }
*/
}
