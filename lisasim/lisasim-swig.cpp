#include <iostream>
#include <fstream>
using namespace std;

#include "lisasim.h"

#include <time.h>
#include <sys/time.h>
#include <math.h>

void printp(LISA *lisa,double t)
{
    Vector p;

    for(int i=1;i<4;i++) {
	lisa->putp(p,i,t);
	printf("%d(%f): %f %f %f\n",i,t,p[0],p[1],p[2]);
    }
}

void printn(LISA *lisa,double t)
{
    Vector n;

    for(int i=1;i<4;i++) {
	lisa->putn(n,i,t);
	printf("%d(%f): %f %f %f\n",i,t,n[0],n[1],n[2]);
    }
}

void printw(Wave *wave,double t)
{
    Tensor h;

    wave->putwave(h,t);
    printf("(%f): %f %f %f %f %f %f %f %f %f\n",t,h[0][0],h[0][1],h[0][2],h[1][0],h[1][1],h[1][2],h[2][0],h[2][1],h[2][2]);
}

extern void retardone(LISA *lisa,int ret,double t,double *retardedtime,double *totalretardbaseline,double *totalretardaccurate);

double retardation(LISA *lisa,int ret1,int ret2,int ret3,int ret4,int ret5,int ret6,int ret7,int ret8,double t) {
    double retardedtime = t;

    double totalretardbaseline = 0.0;
    double totalretardaccurate = 0.0;

    if(ret8 != 0) retardone(lisa,ret8,t,&retardedtime,&totalretardbaseline,&totalretardaccurate);
    if(ret7 != 0) retardone(lisa,ret7,t,&retardedtime,&totalretardbaseline,&totalretardaccurate);
    if(ret6 != 0) retardone(lisa,ret6,t,&retardedtime,&totalretardbaseline,&totalretardaccurate);
    if(ret5 != 0) retardone(lisa,ret5,t,&retardedtime,&totalretardbaseline,&totalretardaccurate);
    if(ret4 != 0) retardone(lisa,ret4,t,&retardedtime,&totalretardbaseline,&totalretardaccurate);
    if(ret3 != 0) retardone(lisa,ret3,t,&retardedtime,&totalretardbaseline,&totalretardaccurate);
    if(ret2 != 0) retardone(lisa,ret2,t,&retardedtime,&totalretardbaseline,&totalretardaccurate);
    if(ret1 != 0) retardone(lisa,ret1,t,&retardedtime,&totalretardbaseline,&totalretardaccurate);

    return totalretardbaseline + totalretardaccurate;
}

double observable(TDI *mytdi,char **obs,double t) {
    switch((*obs)[0]) {
    case 't':
	return t;
	break;
    case 'X':
	if((*obs)[1] == 'm') {
	    (*obs)++;
	    return mytdi->Xm(t);
	} else if((*obs)[1] == '1') {
	    (*obs)++;
	    return mytdi->X1(t);
	} else if((*obs)[1] == '2') {
	    (*obs)++;
	    return mytdi->X2(t);
	} else if((*obs)[1] == '3') {
	    (*obs)++;
	    return mytdi->X3(t);
	} else
	    return mytdi->X(t);
	break;
    case 'Y':
	if((*obs)[1] == 'm') {
	    (*obs)++;
	    return mytdi->Ym(t);
	} else
	    return mytdi->Y(t);          
	break;
    case 'Z':
	if((*obs)[1] == 'm') {
	    (*obs)++;
	    return mytdi->Zm(t);
	} else
	    return mytdi->Z(t);           
	break;
    case 'a':
	if((*obs)[1] == 'm') {
	    (*obs)++;
	    return mytdi->alpham(t);
	} else if((*obs)[1] == '1') {
	    (*obs)++;
	    return mytdi->alpha1(t);
	} else if((*obs)[1] == '2') {
	    (*obs)++;
	    return mytdi->alpha2(t);
	} else if((*obs)[1] == '3') {
	    (*obs)++;
	    return mytdi->alpha3(t);
	} else
	    return mytdi->alpha(t);
	break;
    case 'b':
	if((*obs)[1] == 'm') {
	    (*obs)++;
	    return mytdi->betam(t);
	} else
	    return mytdi->beta(t);          
	break;
    case 'g':
	if((*obs)[1] == 'm') {
	    (*obs)++;
	    return mytdi->gammam(t);
	} else
	    return mytdi->gamma(t);          
	break;
    case 'z':
	return mytdi->zeta(t);
	break;
    case 'P':
	return mytdi->P(t);
	break;
    case 'E':
	return mytdi->E(t);
	break;
    case 'U':
	return mytdi->U(t);
	break;
    case 'y':
	return mytdi->y(3,1,2,0,0,0,t);
	break;
    case 'w':
	return mytdi->z(3,1,2,0,0,0,0,t);
	break;
    default:
	break;
    }

    return 0.0;
}

void printtdi(char *filename,TDI *mytdi,int samples,double samplingtime,char *observables) {
    ofstream outfile(filename);    

    mytdi->reset();
    
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
                cout << "\rEstimating " << floor(sest/60.0) << "m" << floor(sest-60.0*floor(sest/60.0)) << "s to completion [" << floor(speed) << " (multi)samples/s]                    " << flush;
        }

        char *obs = observables;

        while(obs[0]) {
	    outfile << observable(mytdi,&obs,t) << " ";
            obs++;
        }
    
        outfile << endl;
    }

    gettimeofday(&tv,0);
    double lapse = (1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec - timebeg);
    double speed = samples / lapse;
    cout << "\rCompleted in " << floor(lapse/60.0) << "m" << floor(lapse-60.0*floor(lapse/60.0)) << "s [" << floor(speed) << " (multi)samples/s]                      " << endl;
}

void settdi(double *array, TDI *mytdi,int samples,double samplingtime,char *observables) {
    mytdi->reset();
    
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
                cout << "\rEstimating " << floor(sest/60.0) << "m" << floor(sest-60.0*floor(sest/60.0)) << "s to completion [" << floor(speed) << " (multi)samples/s]                    " << flush;
        }

        char *obs = observables;

	// for the moment, only the first observable is returned

	array[i] = observable(mytdi,&obs,t);
    }

    gettimeofday(&tv,0);
    double lapse = (1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec - timebeg);
    double speed = samples / lapse;
    cout << "\rCompleted in " << floor(lapse/60.0) << "m" << floor(lapse-60.0*floor(lapse/60.0)) << "s [" << floor(speed) << " (multi)samples/s]                      " << endl;
}

void setabg(double *aa, double *ab, double *ag, TDI *mytdi,int samples,double samplingtime) {
    mytdi->reset();
    
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
                cout << "\rEstimating " << floor(sest/60.0) << "m" << floor(sest-60.0*floor(sest/60.0)) << "s to completion [" << floor(speed) << " (multi)samples/s]                    " << flush;
        }

	aa[i] = mytdi->alpha(t);
	ab[i] = mytdi->beta(t);
	ag[i] = mytdi->gamma(t);
    }

    gettimeofday(&tv,0);
    double lapse = (1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec - timebeg);
    double speed = samples / lapse;
    cout << "\rCompleted in " << floor(lapse/60.0) << "m" << floor(lapse-60.0*floor(lapse/60.0)) << "s [" << floor(speed) << " (multi)samples/s]                      " << endl;
}

void setabgx(double *aa, double *ab, double *ag, double *ax, TDI *mytdi,int samples,double samplingtime) {
    mytdi->reset();
    
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
                cout << "\rEstimating " << floor(sest/60.0) << "m" << floor(sest-60.0*floor(sest/60.0)) << "s to completion [" << floor(speed) << " (multi)samples/s]                    " << flush;
        }

	aa[i] = mytdi->alpha(t);
	ab[i] = mytdi->beta(t);
	ag[i] = mytdi->gamma(t);
	ax[i] = mytdi->X(t);
    }

    gettimeofday(&tv,0);
    double lapse = (1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec - timebeg);
    double speed = samples / lapse;
    cout << "\rCompleted in " << floor(lapse/60.0) << "m" << floor(lapse-60.0*floor(lapse/60.0)) << "s [" << floor(speed) << " (multi)samples/s]                      " << endl;
}
