import numpy
cimport numpy

cdef extern from "lisasim-tens.h":
    cdef cppclass CPPVector "Vector":
        pass

cdef extern from "lisasim-lisa.h":
    # proof of principle
    cdef cppclass CPPOriginalLISA "OriginalLISA":
        LISA() except +
        LISA(float, float, float) except +

        void putn(Vector, int, double)
        void putp(Vector, int, double)

        double armlength(int, double)
        double dotarmlength(int, double)
