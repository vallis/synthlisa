#ifndef _LISASIM_LISA_H_
#define _LISASIM_LISA_H_

#include "lisasim-tens.h"

class LISA {
    public:

        LISA() {};

        // default destructor: do I need this?

        virtual ~LISA() {};

        virtual void putn(Vector &n, int arm, double t) = 0;
        virtual void putp(Vector &p, int craft, double t) = 0;

        virtual void putn(double n[], int arm, double t) {
            Vector temp;
            putn(temp, arm, t);
            n[0] = temp[0];
            n[1] = temp[1];
            n[2] = temp[2];
	}
        
        virtual void putp(double p[], int craft, double t) {
            Vector temp;
            putp(temp, craft, t);
            p[0] = temp[0];
            p[1] = temp[1];
            p[2] = temp[2];
	}
	  
        virtual double armlength(int arm, double t);
};

// use the following two classes only for TDInoise

class OriginalLISA : public LISA {
    public:
    
    double L[4];

    // accept the armlength in seconds

    OriginalLISA(double armlength[]);

    // we need to define these, but we leave them empty for the moment

    void putn(Vector &n, int arm, double t) {}
    void putp(Vector &p, int craft, double t) {}

    double armlength(int arm, double t);
};

class ModifiedLISA : public LISA {
    public:
    
    double L[4], Lp[4];

    // accept the armlength in seconds

    ModifiedLISA(double armlength[],double armlengthp[]);

    // we need to define these, but we leave them empty for the moment

    void putn(Vector &n, int arm, double t) {}
    void putp(Vector &p, int craft, double t) {}

    double armlength(int arm, double t);
};

// The following classes model the geometry of LISA and can be used also for "signal" TDI

class CircularRotating : public LISA {
    public:
    
        double scriptl;
        double R;
        double L;

        double eta0;
        double xi0;

        // Trick: we use 1-3 indexing for LISA positions and vectors, so we need to allocate 4

        Vector initn[4];
        Vector initp[4];
        
        double rotationtime;
        Vector center;
        Tensor rotation;

        CircularRotating(double eta0 = 0.0,double xi0 = 0.0,double sw = 1.0);
        
        void settime(double t);

        void putn(Vector &n,int arm,double t);

        void putp(Vector &p,int craft,double t);
        
        double armlength(int arm, double t);
};

class MontanaEccentric : public LISA {
    public:
    
        // orbital radius of the guiding center (yrs)
        // LISA simulator has Rgc = 1.49597870660e11 m

        static const double Rgc = 0.0000158233;

        // mean arm length of the LISA detector (yrs)
        // LISA simulator has L = 5.0e9 m

        static const double L = 5.288624035993024E-7;

        // eccentricity of the LISA spacecraft orbits
        // LISA simulator has L/(2.0*sqrt(3.0)*Rgc)
        // with L = 5.0e9, equal to 0.00964837

        static const double ecc = 0.00964839;
        
        // initial azimuthal position of the guiding center
    
        double kappa;
        
        // initial orientation of the spacecraft
        
        double lambda;
        
        // caching the positions of spacecraft
        
        Vector cachep[4];
        double cachetime[4];
        
        // public methods (perhaps not all should be)
        
        MontanaEccentric(double kappa0 = 0.0,double lambda0 = 0.0);

        void settime(int craft,double t);
        
        void putn(Vector &n,int arm,double t);

        void putp(Vector &p,int craft,double t);
};

#endif /* _LISASIM_LISA_H_ */
