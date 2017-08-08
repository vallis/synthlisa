import numpy
cimport numpy

cimport synthlisa_defs as defs


speedoflight = 2.99792458e8     # m/s
year2sec = 365.25636 * 86400    # s

standard_armlength = 2.5e9      # m


cdef class LISA:
    cdef defs.LISA *lisa

    # standard armlength to 2.5
    def __cinit__(self, armlength=standard_armlength/speedoflight, eccentric=True):
        if eccentric:
            # extra parameters are eta0, xi0, sw, t0; where eta0 and xi0 are
            # the true anomaly of the LISA guiding center and the initial phase
            # of the LISA array at t=t0; sw<0 will swap spacecraft 2 and 3,
            # so that the spacecraft sequence 1 -> 2 -> 3 -> 1 goes ccw
            # as seen from above.

            self.lisa = new defs.EccentricInclined(armlength, 0.0, 0.0, 1, 0.0)
        else:
            self.lisa = new defs.CircularRotating(armlength, 0.0, 0.0, 1, 0.0)

    def __dealloc__(self):
        del self.lisa

    def p(self, int arm, double t):
        cdef defs.Vector v        
        self.lisa.putp(v, arm, t)
        
        return numpy.array([v[0],v[1],v[2]],'d')

    def n(self, int arm, double t):
        cdef defs.Vector v
        self.lisa.putn(v, arm, t)
        
        return numpy.array([v[0],v[1],v[2]],'d')

    def armlength(self, int arm, double t):
        return self.lisa.armlength(arm, t)

    def dotarmlength(self, int arm, double t):
        return self.lisa.dotarmlength(arm, t)


cdef class Wave:
    cdef defs.Wave *wave

    @property
    def beta(self):
        return self.wave.beta

    @property
    def lam(self):
        return self.wave.lam

    @property
    def pol(self):
        return self.wave.pol

    @property
    def k(self):
        return numpy.array([self.wave.k[0],self.wave.k[1],self.wave.k[2]],'d')

    @property
    def pp(self):
        return numpy.array([[self.wave.pp[0][0],self.wave.pp[0][1],self.wave.pp[0][2]],
                            [self.wave.pp[1][0],self.wave.pp[1][1],self.wave.pp[1][2]],
                            [self.wave.pp[2][0],self.wave.pp[2][1],self.wave.pp[2][2]]],'d')

    @property
    def pc(self):
        return numpy.array([[self.wave.pc[0][0],self.wave.pc[0][1],self.wave.pc[0][2]],
                            [self.wave.pc[1][0],self.wave.pc[1][1],self.wave.pc[1][2]],
                            [self.wave.pc[2][0],self.wave.pc[2][1],self.wave.pc[2][2]]],'d')

    def hp(self, t):
        return self.wave.hp(t)

    def hc(self, t):
        return self.wave.hc(t)

    def htens(self, t):
        cdef defs.Tensor h
        self.wave.putwave(h, t)

        return numpy.array([[h[0][0],h[0][1],h[0][2]],
                            [h[1][0],h[1][1],h[1][2]],
                            [h[2][0],h[2][1],h[2][2]]],'d')

# TODO: may want to rationalize the argument order and names

cdef class SimpleBinary(Wave):
    def __cinit__(self, freq, initphi, inc, amp, b, l, p):
        self.wave = new defs.SimpleBinary(freq, initphi, inc, amp, b, l, p)

    def __dealloc__(self):
        del self.wave

cdef class GalacticBinary(Wave):
    def __cinit__(self, freq, freqdot, b, l, amp, inc, p, initphi, fddot=0, epsilon=0):
        self.wave = new defs.GalacticBinary(freq, freqdot, b, l, amp, inc, p, initphi, fddot, epsilon)

    def __dealloc__(self):
        del self.wave


cdef class TDI:
    cdef defs.TDI *tdi

    def reset(self):
        self.tdi.reset()

    def Xm(self, double t):
        return self.tdi.Xm(t)

    def Ym(self, double t):
        return self.tdi.Ym(t)

    def Zm(self, double t):
        return self.tdi.Zm(t)

    def X1(self, double t):
        return self.tdi.X1(t)

    def X2(self, double t):
        return self.tdi.X2(t)

    def X3(self, double t):
        return self.tdi.X3(t)

    def y(self, int send, int link, int recv, double t, int ret1 = 0, int ret2 = 0, int ret3 = 0, int ret4 = 0, int ret5 = 0, int ret6 = 0, int ret7 = 0):
        return self.tdi.y(send, link, recv, ret1, ret2, ret3, ret4, ret5, ret6, ret7, t)

    def z(self, int send, int link, int recv, double t, int ret1 = 0, int ret2 = 0, int ret3 = 0, int ret4 = 0, int ret5 = 0, int ret6 = 0, int ret7 = 0, int ret8 = 0):
        return self.tdi.z(send, link, recv, ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, t)

cdef class TDIsignal(TDI):
    def __cinit__(self, LISA lisa, Wave wave):
        self.tdi = new defs.TDIsignal(lisa.lisa, wave.wave)

    def __dealloc__(self):
        del self.tdi
