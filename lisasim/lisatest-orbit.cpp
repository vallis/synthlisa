#include <iostream.h>
#include <math.h>
#include <unistd.h>
#include "lisasim.h"

// compares the CircularRotating and EccentricInclined LISA models

int main(int argc, char **argv) {
    // take as only input the number of positions to sample throughout the year

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

    // need to pass eta0 and xi0; we need to set xi0 and invert crafts 2 and 3 to get the same orbits

    CircularRotating mylisa(0.0,3.0*M_PI/2.0,-1.0);

    // need to pass the initial azimuthal position of the guiding center and the initial orientation of the spacecraft

    EccentricInclined thlisa(0.0,0.0);

    ofstream dout("distances.txt");
    ofstream myout("circular.txt");
    ofstream thout("eccentric.txt");    

    // remember that the LISA functions take times as fractions of a year

    Vector spc1, spc2, spc3;

    for(int i=0;i<samples;i++) {
        double time = (1.0*i)/samples;
        
//        mylisa.putp(spc1,1,time); mylisa.putp(spc2,2,time); mylisa.putp(spc3,3,time);

//        myout << time << " " << spc1[0] << " " << spc1[1] << " " << spc1[2];
//        myout         << " " << spc2[0] << " " << spc2[1] << " " << spc2[2];
//        myout         << " " << spc3[0] << " " << spc3[1] << " " << spc3[2] << endl;

//        thlisa.putp(spc1,1,time); thlisa.putp(spc2,2,time); thlisa.putp(spc3,3,time);

//        thout << time << " " << spc1[0] << " " << spc1[1] << " " << spc1[2];
//        thout         << " " << spc2[0] << " " << spc2[1] << " " << spc2[2];
//        thout         << " " << spc3[0] << " " << spc3[1] << " " << spc3[2] << endl;
        
        dout << thlisa.armlength(1,time) << " " << thlisa.armlength(1,time) << " " << thlisa.armlength(1,time) << endl;
    }
}
