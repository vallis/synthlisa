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
    
    struct timeval tv;
    
    gettimeofday(&tv,0);
    double timebeg = 1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec;
    
    for(int i=0;i<samples;i++) {
        double t = i * samplingtime;

        if(i % 16384 == 0 && i != 0) {
            gettimeofday(&tv,0);
            double speed = (1.0*i) / (1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec - timebeg);
            double sest = (samples - i) / speed;

            if(speed > 0.0)
                cout << "\rEstimating " << floor(sest/60.0) << "m" << floor(sest-60.0*floor(sest/60.0)) << "s to completion [" << floor(speed) << " (multi)samples/s]                    ";
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
                case 'U':
                    outfile << mynoise->U(t);
                    break;
                case 'y':
                    outfile << mynoise->y(3,1,2,0,0,0,t);
                    break;
                case 'w':
                    outfile << mynoise->z(3,1,2,0,0,0,0,t);
                    break;
                default:
                    break;
            }
        
            outfile << " ";
            obs++;
        }
    
        outfile << endl;
    }

    gettimeofday(&tv,0);
    double lapse = (1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec - timebeg);
    double speed = samples / lapse;
    cout << "\rCompleted in " << floor(lapse/60.0) << "m" << floor(lapse-60.0*floor(lapse/60.0)) << "s [" << floor(speed) << " (multi)samples/s]                      " << endl;
}

void printsignal(char *filename,TDI *mysignal,int samples,double samplingtime,char *observables) {
    ofstream outfile(filename);    

    // Can't reset signals at the moment, but no need to do it

    // mysignal->reset();
    
    struct timeval tv;
    
    gettimeofday(&tv,0);
    double timebeg = 1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec;
    
    (unsigned long)time(0);
    
    for(int i=0;i<samples;i++) {
        double t = i * samplingtime;

        if(i % 16384 == 0 && i != 0) {
            gettimeofday(&tv,0);
            double speed = (1.0*i) / (1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec - timebeg);
            double sest = (samples - i) / speed;

            if(speed > 0.0)
                cout << "\rEstimating " << floor(sest/60.0) << "m" << floor(sest-60.0*floor(sest/60.0)) << "s to completion [" << floor(speed) << " (multi)samples/s]                    ";
        }

        char *obs = observables;

        while(obs[0]) {
            switch(obs[0]) {
                case 't':
                    outfile << t;
                    break;
                case 'X':
                    outfile << mysignal->X(t);
                    break;
                case 'Y':
                    outfile << mysignal->Y(t);
                    break;
                case 'Z':
                    outfile << mysignal->Z(t);
                    break;
                case 'a':
                    outfile << mysignal->alpha(t);
                    break;
                case 'b':
                    outfile << mysignal->beta(t);
                    break;
                case 'g':
                    outfile << mysignal->gamma(t);
                    break;
                case 'z':
                    outfile << mysignal->zeta(t);
                    break;
                case 'P':
                    outfile << mysignal->P(t);
                    break;
                case 'E':
                    outfile << mysignal->E(t);
                    break;
                case 'U':
                    outfile << mysignal->U(t);
                    break;
                default:
                    break;
            }
        
            outfile << " ";
            obs++;
        }
    
        outfile << endl;
    }

    gettimeofday(&tv,0);
    double lapse = (1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec - timebeg);
    double speed = samples / lapse;
    cout << "\rCompleted in " << floor(lapse/60.0) << "m" << floor(lapse-60.0*floor(lapse/60.0)) << "s [" << floor(speed) << " (multi)samples/s]                   " << endl;
}
