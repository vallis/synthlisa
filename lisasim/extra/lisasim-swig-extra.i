// FIX: the following have all rather specific uses, and may be
//      replaced by the suitable use of PyLISA. They are moved to this
//      private C++ and interface file.

// NOTE: this interface file is incomplete and won't compile without changes

initsave(NoisyLISA)

class NoisyLISA : public LISA {
public:
    NoisyLISA(LISA *clean,double starm,double sdarm);

    // nontrivial destructor should appear here

    ~NoisyLISA(); 

    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

class NominalLISA : public LISA {
 public:
    double Lnom, emod, cmod, toff;

    NominalLISA(double eta0,double xi0,double sw,double t0);
    ~NominalLISA();

    void setparameters(double l,double cm,double em,double toff);
    void setparameters(double cm,double em,double toff);
    void setparameters3(double l,double cm,double em);

    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

%apply double PYTHON_SEQUENCE_DOUBLE[ANY] {double dl[6], double dldt[6]}

class LinearLISA : public LISA {
 public:
    LinearLISA(double eta0,double xi0,double sw,double t0);
    ~LinearLISA();

    void settimeoffset(double toff);
    void setparameters(double dl[6],double dldt[6]);
    
    double armlengtherror(int arm, double t);

    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};

initsave(MeasureLISA)

class MeasureLISA : public LISA {
 public:
    MeasureLISA(LISA *clean,double starm,double sdarm,int swindow = 1);
    ~MeasureLISA();
        
    void putn(Vector &outvector, int arm, double t);
    void putp(Vector &outvector, int craft, double t);

    double armlength(int arm, double t);

    double armlengthbaseline(int arm, double t);
    double armlengthaccurate(int arm, double t);
};
