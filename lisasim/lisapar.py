import mpi

import synthlisa
import numpy

import time
import pickle

# tetrad must be given as (lisa, srcfunc, parameters, obs)

class LISApar:
    """LISApar() instantiates a class that provides access to
    the synthLISA parallel-computing methods (so far, only
    getobsp); the LISApar constructor calls MPI_Init, and
    the size of the CPU group under which python is being
    executed, as well as the rank of each CPU, are available
    as lp.size and lp.rank, where lp is the LISApar
    instance. MPI_Finalize is called upon destruction of the
    LISApar instance. You should probably create only one
    LISApar instance."""

    def __init__(self):
        """Constructor. See above."""
        self.rank, self.size = mpi.init()

    def __del__(self):
        """Destructor. See above."""
        mpi.finalize()

    def getobsp(self,snum,stime,tetrad,zerotime=0.0,debug=0):

        """
        
        LISApar.getobsp(length,deltat,tetrad,zerotime=0.0)
        is the parallel-computing equivalent of getobs and
        getobsc, and it is used to compute the TDI responses
        of large sets of Wave objects. It must be called
        from an instance of LISApar, with the following
        parameters:
        
        - length is the total length of the TDI-observable
          arrays that will be returned;
        
        - deltat is the cadence of the time series;
        
        - zerotime is the initial time for the time series;
        
        - tetrad is a tuple (lisa,wavefactory,parameterlist,
          observables) of four elements:

          * lisa is an instance of a LISA class, which
            should be the same for every CPU taking part in
            the computation;

          * wavefactory is a Python function taking any
            number of parameters, and returning an instance of
            a synthLISA Wave object; the function must be
            defined for every CPU taking part in the
            computation;

          * parameterlist is a list of source parameters (or
            of parameter n-tuples, if wavefactory takes more
            than one parameter), which will be distributed
            among the CPUs, and passed to the Wave Factory to
            construct synthLISA Wave objects; the parameter
            sets need to be defined only on the root CPU, but
            it won't hurt to define them everywhere. They can
            contain any Python types (they are pickled before
            distribution), but not synthLISA objects;

          * observables is a list or tuple of TDI
            observables, which must be given as unbound
            methods, such as synthlisa.TDI.X1 or
            synthlisa.TDI.time.
        
        The distribution of the parameter sets among the
        CPUs tries to balance the load of the computation.
        If the number of sources is not divisible by the
        number of CPUs, it will assign a smaller number of
        sources to the root CPU, and the same number of
        sources to all other CPUs."""

        # accept four levels (0-4) of debugging info
        
        inittime = time.time()

        myrank = self.rank
        size = self.size
    
        try:
            (lisa, srcfunc, parameters, obs) = tetrad
        except:
            if myrank == 0:
                print "LISApar.getobsp(...): third parameter must be a 4-tuple containing a",
                print "LISA instance, a Wave factory, an array of parameters for the factory,",
                print "and a set of TDI observables given as class methods (such as synthlisa.TDI.X)."
            raise IndexError
        
        if type(parameters) not in (list,tuple,numpy.ndarray):
            if myrank == 0:
                print "LISApar.getobsp(...): needs a list of parameters to feed to the factory!"
            raise IndexError
                
        if size == 1:
            if myrank == 0:
                print "LISApar.getobsp(...): must be run with more than one cpu!"
            raise NotImplementedError
    
        if size > len(parameters):
            if myrank == 0:
                print "LISApar.getobsp(...): needs to run with more sources than cpus!"
            raise IndexError
    
        # root may get zero processors
    
        blocksize, remain = divmod(len(parameters),size)

        if remain > 0:
            blockadd, remain = divmod(remain,size-1)
            blocksize = blocksize + blockadd

        if myrank == 0 and debug > 2:
            print "Standard block: ", blocksize,
            print "; root block: ", len(parameters) - blocksize * (size-1)
    
        if myrank == 0:
            if debug > 3:
                print "Preparing for parallel execution..."
    
            for cpu in range(1,size):
                blockstart, blockend = (cpu-1)*blocksize, cpu*blocksize
    
                serial_pars = pickle.dumps(parameters[blockstart:blockend])         
                len_pars = len(serial_pars)
    
                mpi.isend(len_pars,1,mpi.MPI_INT,cpu,0,mpi.MPI_COMM_WORLD)
                mpi.isend(serial_pars,len_pars,mpi.MPI_CHAR,cpu,1,mpi.MPI_COMM_WORLD)
    
            mypars = parameters[blockend:]
        else:
            len_pars = mpi.recv(1,mpi.MPI_INT,0,0,mpi.MPI_COMM_WORLD)
            serial_pars = mpi.recv(len_pars,mpi.MPI_CHAR,0,1,mpi.MPI_COMM_WORLD)
    
            mypars = pickle.loads(serial_pars)
    
        if debug > 2:
            print "CPU ", myrank, " received ", len(mypars), " source parameters ", mypars
        
        try:
            if type(mypars[0]) in (list,tuple,numpy.ndarray):
                sources = map(lambda x: srcfunc(*x),mypars)
            else:
                sources = map(srcfunc,mypars)
        
            if len(filter(lambda x: not isinstance(x,synthlisa.Wave),sources)) > 0:
                raise TypeError
        except:
            if myrank == 0:
                print "LISApar.getobsp(...): srcfunc must return a synthlisa.Wave when applied",
                print "to each element of the parameter list"
            raise TypeError
    
        if debug > 3:
            print "CPU ", myrank, " created sources ", sources
    
        wavearray = synthlisa.WaveArray(sources)
    
        if not isinstance(lisa,synthlisa.LISA):
            if myrank == 0:
                print "LISApar.getobsp(...): lisa must be an instance of synthlisa.LISA."
            raise TypeError
    
        tdisignal = synthlisa.TDIsignal(lisa,wavearray)
    
        # is it possible to permanently bind an unbound method?
        # yes, by doing bound_obs = obs.__get__(tdisignal)
        # but it's not clear this will yield a faster call
    
        if type(obs) == list or type(obs) == tuple:
            multobs = len(obs)
    
            array = numpy.zeros((snum,multobs),dtype='d')
            for i in numpy.arange(0,snum):
                for j in range(0,multobs):
                    array[i,j] = obs[j](tdisignal,zerotime+i*stime)
        else:
            multobs = 1
    
            array = numpy.zeros(snum,dtype='d')
            for i in numpy.arange(0,snum):
                array[i] = obs(tdisignal,zerotime+i*stime)
    
        sumresults = mpi.reduce(array,snum*multobs,mpi.MPI_DOUBLE,mpi.MPI_SUM,0,mpi.MPI_COMM_WORLD)
    
        if myrank == 0 and debug > 0:
            currenttime = time.time() - inittime
    
            vel = snum/currenttime
            print "Completed in %d s [%d (multi)samples/s]." % (int(currenttime),int(vel))
    
        if myrank == 0:
            if multobs == 1:
                return sumresults
            else:
                return sumresults.reshape(snum,multobs)
        else:
            return None
