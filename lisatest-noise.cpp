#include <iostream.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include "lisasim.h"
#include "lisasim-noise.h"

// outputs simple noise from RingNoise classes
// and dimensioned noise from InterpolateNoise

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

    // use arbitrary buffer length of 256

    // here we're actually passing a long (on PPC); should work as long
    // is the same as int

    idum = -time(0);

    // RingNoise myrnoise(256);
    // DiffNoise mydnoise(256);
    // IntNoise myinoise(256, 0.9999);

    // this is the nominal time interval of the noises

    double dstep = 1.0;

    InterpolateNoise lasernoise(dstep,256.0,1.1e-26,0.0);
    InterpolateNoise proofnoise(dstep,256.0,2.5e-48,-2.0);
    InterpolateNoise shotnoise(dstep,256.0,1.8e-37,2.0);

    // ExpGaussNoise expnoise(dstep,20.00,1.00/0.01,1.1e-26);

    // set this to a fraction of dstep to oversample by interpolation
    // (see the results on the spectra!)

    double sstep = 0.1;

    // ofstream rout("noise-white.txt");
    // ofstream dout("noise-diff.txt");
    // ofstream iout("noise-int.txt");    

    ofstream lout("noise-laser.txt");
    ofstream pout("noise-proof.txt");
    ofstream sout("noise-shot.txt");

    // ofstream eout("noise-exp.txt");

    for(int i=0;i<samples;i++) {
    //    rout << myrnoise[i] << endl;
    //    dout << mydnoise[i] << endl;
    //    iout << myinoise[i] << endl;
        
        lout << lasernoise[sstep * i] << endl;
        pout << proofnoise[sstep * i] << endl;
        sout << shotnoise[sstep * i] << endl;

    //    eout << expnoise[sstep * i] << endl;
    } 

    // out-of-order test
/*
    double *doublebuffer = new double[samples];

    ExpGaussNoise expnoise2(dstep,20.00,1.00/0.01,1.1e-26);
    
    for(int i=0;i<samples;i++) {
        if(i % 2 == 0) {
            doublebuffer[i+1] = expnoise2[sstep * (i+1)];
        } else {
            doublebuffer[i-1] = expnoise2[sstep * (i-1)];
        }
    }
    
    ofstream eout2("noise-exp2.txt");    

    for(int i=0;i<samples;i++)
        eout2 << doublebuffer[i] << endl;

    delete doublebuffer;

    // now increase the correlation

    ExpGaussNoise expnoise3(dstep,20.00,1.00/0.5,1.1e-26);    

    ofstream eout3("noise-exp3.txt");

    for(int i=0;i<samples;i++)
        eout3 << expnoise3[sstep * i] << endl;

    // again
        
    ExpGaussNoise expnoise4(dstep,20.00,1.00/1.0,1.1e-26);    

    ofstream eout4("noise-exp4.txt");

    for(int i=0;i<samples;i++)
        eout4 << expnoise4[sstep * i] << endl;
        
    // again
        
    ExpGaussNoise expnoise5(dstep,20.00,1.00/2.0,1.1e-26);    

    ofstream eout5("noise-exp5.txt");

    for(int i=0;i<samples;i++)
        eout5 << expnoise5[sstep * i] << endl;
*/
}
