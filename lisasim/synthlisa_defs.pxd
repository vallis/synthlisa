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
    cdef cppclass Wave:
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
