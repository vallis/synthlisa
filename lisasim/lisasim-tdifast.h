#ifndef _LISASIM_TDIFAST_H_
#define _LISASIM_TDIFAST_H_

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

#endif /* _LISASIM_TDIFAST_H_ */
