#include <iostream.h>
#include "lisasim.h"
#include <time.h>
#include <math.h>

TDInoise *stdnoise(LISA *mylisa) {
    return new TDInoise(mylisa,1.0,2.5e-48,1.0,1.8e-37,1.0,1.1e-26,1.0e-6);
}

void printnoise(char *filename,TDInoise *mynoise,int samples,double samplingtime,char *observables) {
    ofstream outfile(filename);    

    mynoise->reset();
    
    long timebeg = (unsigned long)clock();
    
    for(int i=0;i<samples;i++) {
        double t = i * samplingtime;

        if(i==8192) {
            double sest = ((samples - 8192.0)/8192.0)*(1.0*((unsigned long)clock() - timebeg))/CLOCKS_PER_SEC;        
            cout << "Estimating " << floor(sest/60.0) << "m" << floor(remainder(sest,60.0)) << "s (CPU) to completion" << endl;
        }

        char *obs = observables;

        while(obs[0]) {
            switch(obs[0]) {
                case 't':
                    outfile << t;
                    break;
                case 'X':
                    if(obs[1] == 'm') {
                        outfile << mynoise->Xm(t);
                        obs++;
                    } else {
                        outfile << mynoise->X(t);
                    }
                
                    break;
                case 'Y':
                    outfile << mynoise->Y(t);
                    break;
                case 'Z':
                    outfile << mynoise->Z(t);
                    break;
                case 'a':
                    outfile << mynoise->alpha(t);
                    break;
                case 'b':
                    outfile << mynoise->beta(t);
                    break;
                case 'g':
                    outfile << mynoise->gamma(t);
                    break;
                case 'z':
                    outfile << mynoise->zeta(t);
                    break;
                case 'P':
                    outfile << mynoise->P(t);
                    break;
                case 'E':
                    outfile << mynoise->E(t);
                    break;
                case 'y':
                    outfile << mynoise->y(3,-2,1,0,0,0,t);
                    break;
                case 'w':
                    outfile << mynoise->z(3,-2,1,0,0,0,0,t);
                    break;
                default:
                    break;
            }
        
            outfile << " ";
            obs++;
        }
    
        outfile << endl;
    }
}
