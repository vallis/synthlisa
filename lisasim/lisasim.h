#include "nr.h"
#include "nrutil.h"

class Tensor : public Mat_IO_DP {
    public:

    Tensor();
    
    void seteuler(double d, double a, double p);
    
    void setproduct(Tensor &fac1, Tensor &fac2);
                                
    void settranspose();
};

class Vector : public Vec_IO_DP {
    public:
    
    Vector();
    
    void setproduct(Tensor &mat, Vector &vec);

    double dotproduct(Vector &vec);
};

class Wave {
    public:

        // position in the sky

        double dec, asc, pol;

        // k vector

        Vector k;
	double kArray[3];

        // polarization tensors
    
        Tensor pp, pc;
	double ppArray[9], pcArray[9];

        // default constructor

        Wave(double d, double a, double p);

        // default destructor; do I need this?

        virtual ~Wave() {};

        virtual double hp(double t) = 0;
        virtual double hc(double t) = 0;  

        void putwave(Tensor &h, double t);
        void putwave(double **h, double t);
};

class SimpleBinary : public Wave {
    public:

        // frequency, initial phase

        double f, phi0;

        // inclination, polarization amplitudes

        double i, a, ap, ac;

        SimpleBinary(double freq, double initphi, double inc, double amp, double d, double a, double p);

        double hp(double t);
        double hc(double t);
};

class SimpleMonochromatic : public Wave {
    public:

        // frequency

        double f;

        // polarization angles using John's convention

        double gm, ph, ap, ac;

        SimpleMonochromatic(double freq, double phi, double gamma, double amp, double d, double a, double p);

        double hp(double t);
        double hc(double t);
};

class LISA {
    public:

        LISA() {};

        // default destructor: do I need this?

        virtual ~LISA() {};

        virtual void putn(Vector &n, int arm, double t) = 0;
        virtual void putp(Vector &p, int craft, double t) = 0;

        virtual void putn(double n[], int arm, double t){
	  Vector temp;
	  putn(temp, arm, t);
	  n[0] = temp[0];
	  n[1] = temp[1];
	  n[2] = temp[2];
	}
        virtual void putp(double p[], int craft, double t){
	  Vector temp;
	  putp(temp, craft, t);
	  p[0] = temp[0];
	  p[1] = temp[1];
	  p[2] = temp[2];
	}
	  
        virtual double armlength(int arm, double t);
};

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

class TDI {
    public:
        LISA *lisa;
        Wave *wave;

        TDI(LISA *mylisa, Wave *mywave) {
            lisa = mylisa;
            wave = mywave;
        }


        double M(double t);
        double N(double t);
        double O(double t);
    
        double X(double t);
        double Y(double t);
        double Z(double t);
    
        double alpha(double t);
        double beta(double t);
        double gamma(double t);
    
        double zeta(double t);

    private:
	double psi(int arm, double t, double twave);
//        double psi(int arm, double t);
    
        double retard(int craft, double t);
    
        double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);




};

class TDIfast : public TDI {
 public:

        TDIfast(LISA *mylisa, Wave *mywave, double srate, long samples);
        virtual ~TDIfast() ;
	/*
        double M(double t);
        double N(double t);
        double O(double t);
	*/
        double Xfast(int t);
        double Yfast(int t);
        double Zfast(int t);

	int CacheX(void);
        int CacheY(void);
	int CacheZ(void);
	
	/*
	  double alpha(double t);
	  double beta(double t);
	  double gamma(double t);
	  
	  double zeta(double t);
	*/
 private:
	int samples;
	double srate; // Stored in years

	double **cwave;

	//	double psi(Vector lisan, double twave);
	double psi(double lisan[], double twave);
	//        double psi(int arm, double t, double twave);
	//        double psi(int arm, double t);
    
        double retardfast(int craft, int tIndex);
	
        double y(int send, int link, int recv, int ret1, int ret2, int ret3, int tIndex);

        int yc(int send, int link, int recv, int ret1, int ret2, int ret3);
	
	
	double **Storedarmlength;
	bool *Checkarmlength;
	double  ***Storedputn;
	bool *Checkputn;
	double ***Storedputp;
	bool *Checkputp;
	
	double **StoredRetardedTime;
	bool ***CheckRetardedTime;
	//Vector ******StoredputnRet;
	double ***StoredputnRet;
	bool ****CheckputnRet;

	
};


