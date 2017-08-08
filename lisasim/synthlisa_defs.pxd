cdef extern from "lisasim-tens.h":
    cdef cppclass Vector:
        Vector()
        double& operator[](int)

    cdef cppclass Tensor:
        Tensor()
        double* operator[](int)

cdef extern from "lisasim-lisa.h":
    cdef cppclass LISA:
        void putp(Vector&, int, double)
        void putn(Vector&, int, double)

        double armlength(int, double)
        double dotarmlength(int, double)

    cdef cppclass EccentricInclined(LISA):
        EccentricInclined(double myL,double eta0,double xi0,double sw,double t0) except+

    cdef cppclass CircularRotating(LISA): 
        CircularRotating(double myL,double eta0,double xi0,double sw,double t0) except+

cdef extern from "lisasim-wave.h":
    cdef cppclass WaveObject:
        pass

    cdef cppclass Wave(WaveObject):
        double beta
        double lam "lambda"
        double pol

        Vector k

        Tensor pp
        Tensor pc

        double hp(double)
        double hc(double)

        void putk(Vector&)
        void putwave(Tensor&, double)

    cdef cppclass SimpleBinary(Wave):
        SimpleBinary(double, double, double, double, double, double, double)

    cdef cppclass GalacticBinary(Wave):
        GalacticBinary(double, double, double, double, double, double, double, double, double, double)

cdef extern from "lisasim-tdi.h":
    cdef cppclass TDI:
        void reset()

        double Xm(double)
        double Ym(double)
        double Zm(double)        

        double X1(double)
        double X2(double)
        double X3(double)

        double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t)
        double z(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, int ret8, double t)

cdef extern from "lisasim-tdisignal.h":
    cdef cppclass TDIsignal(TDI):
        TDIsignal(LISA*, WaveObject*)
