#include <iostream.h>
#include <math.h>
#include <unistd.h>
#include "lisasim.h"
#include "lisatest.h"

// this program compares the TDI response to an incoming "Newtonian" signal,
// as computed by the LISA Simulator and by Synthetic LISA

int main(int argc, char **argv) {
    TDI *prova = new TDI(new CircularRotating(0.0,3.0*M_PI/2.0,-1.0),new SimpleBinary(0.00194545,0,1.60238,1.88392e-22,0.000796327,0,0));
    
    cout << prova->X(0) << endl;

    // --- definition of the waveform

    // find orbital separation
    
    double forb = 0.5*fgw;
    double R = pow( G*(M1+M2)*pow(1.0/(2.0*M_PI*forb),2) , 1.0/3.0);
    
    // find overall amplitude
    
    double A = 2.0*M1*M2/(r*R) * pow(G/(c*c),2);
    
    // I don't have a minus in my definition of inclination as in Newtonian.c, so I pass pi - inc [cos(pi-inc) = -cos(inc)]
    // need to adjust my frequency for the different definition of year in the Montana and JPL codes
    // need to convert my ecliptic latitude to the Montana latitude (which is Pi/2 - ...)
    // need to use minus my polarization angle
    
    cout << fgw*1.000702365550482 << " " << phase << " " << (M_PI-inc) << " " << A << " " << 0.5*M_PI - theta << " " << phi << " " << -psi << endl;
    
    SimpleBinary mybinary(fgw*1.000702365550482, phase, (M_PI-inc), A, 0.5*M_PI - theta, phi, -psi);
 
    // --- definition of LISA and TDI

    CircularRotating mylisa(0.0,3.0*M_PI/2.0,-1.0);    
    // MontanaEccentric mylisa(0.0,0.0);
    
    TDI mytdi(&mylisa,&mybinary);
 
    // --- write to file
          
    ofstream myout("newtoniantdi.txt");
  
    long samples = (long)pow(2.0,nsource);
    
    // set timestep (in seconds) for 1 year of observations
     
    double Dt = ObsTime / samples;

    // use this to write fewer samples at the same sampling rate

    long Ntot = samples;    
    // long Ntot = samples/8;

    // compute and write the TDI signals
 
    cout << "Computing " << Ntot << " samples" << endl;
 
    // now, the TDI functions are programmed to take an argument in seconds, which they then convert into years
    // but I need to adjust the conversion factor to the sidereal year used in the Montana code
 
    for(int i=0;i<Ntot;i++) {
        double time = 0.9992978747981357*(i*Dt);
        
        myout << mytdi.X(time) << " " << mytdi.Y(time) << " " << mytdi.Z(time) << endl;
                
        if (i == (int)floor(Ntot/5.))    printf("  20 percent done\n");
        if (i == (int)floor(2.*Ntot/5.)) printf("  40 percent done\n");
        if (i == (int)floor(3.*Ntot/5.)) printf("  60 percent done\n");
        if (i == (int)floor(4.*Ntot/5.)) printf("  80 percent done\n");
        if (i == Ntot-1)                 printf(" 100 percent done\n");
    }
}
