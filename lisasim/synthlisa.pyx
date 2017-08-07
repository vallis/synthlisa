import numpy
cimport numpy

speedoflight = 2.99792458e8     # m/s
year2sec = 365.25636 * 86400    # s

standard_armlength = 2.5e9      # m

cdef extern from "lisasim-tens.h":
    cdef cppclass Vector:
        Vector()
        double& operator[](int)

cdef extern from "lisasim-lisa.h":
    # base class; renamed to avoid collision with Cython object defined below
    cdef cppclass cppLISA "LISA":
        void putp(Vector&, int, double)
        void putn(Vector&, int, double)

        double armlength(int, double)
        double dotarmlength(int, double)

    cdef cppclass EccentricInclined(cppLISA):
        EccentricInclined(double myL,double eta0,double xi0,double sw,double t0) except+

    cdef cppclass CircularRotating(cppLISA): 
        CircularRotating(double myL,double eta0,double xi0,double sw,double t0) except+

cdef class LISA:
    cdef cppLISA *lisa

    # standard armlength to 2.5
    def __cinit__(self, armlength=standard_armlength/speedoflight, eccentric=True):
        if eccentric:
            # extra parameters are eta0, xi0, sw, t0; where eta0 and xi0 are
            # the true anomaly of the LISA guiding center and the initial phase
            # of the LISA array at t=t0; sw<0 will swap spacecraft 2 and 3,
            # so that the spacecraft sequence 1 -> 2 -> 3 -> 1 goes ccw
            # as seen from above.

            self.lisa = new EccentricInclined(armlength, 0.0, 0.0, 1, 0.0)
        else:
            self.lisa = new CircularRotating(armlength, 0.0, 0.0, 1, 0.0)

    def p(self, int arm, double t):
        cdef Vector v        
        self.lisa.putp(v, arm, t)
        
        return numpy.array([v[0],v[1],v[2]],'d')

    def n(self, int arm, double t):
        cdef Vector v
        self.lisa.putn(v, arm, t)
        
        return numpy.array([v[0],v[1],v[2]],'d')

    def armlength(self, int arm, double t):
        return self.lisa.armlength(arm, t)

    def dotarmlength(self, int arm, double t):
        return self.lisa.dotarmlength(arm, t)
