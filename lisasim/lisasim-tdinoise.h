class TDInoise {
    public:
        LISA *lisa;

        InterpolateNoise *pm[4], *pms[4];
        
        // I label shot noises by sending and receiving spacecraft, not by link and receiving
        
        InterpolateNoise *shot[4][4];

        TDInoise(LISA *mylisa, double stproof, double sdproof, double stshot, double sdshot);
        ~TDInoise();

        double X(double t);
        double Y(double t);
        double Z(double t);

        double alpha(double t);
        double beta(double t);
        double gamma(double t);

        double zeta(double t);

        double P(double t);
        
        double E(double t);
    
    private:
        double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
        double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, double t);
};
