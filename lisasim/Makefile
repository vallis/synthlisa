#OPTIONS = -Wall -g 
#OPTIONS = -Wall -O3 -pg -g 
OPTIONS = -Wall -O3 

lisasim: lisasim-lisa.o lisasim-tdi.o lisasim-tens.o lisasim-wave.o lisasim-main.o lisasim-tdifast.o
	g++ $(OPTIONS) lisasim-tdi.o lisasim-lisa.o lisasim-wave.o lisasim-tens.o lisasim-main.o lisasim-tdifast.o -o lisasim

%.o: %.cpp lisasim.h
	g++ $(OPTIONS) -c $<

clean:
	rm -rf *.o lisasim

.PHONEY: clean

SWIG= swig/bin/swig

PYTHON_INCLUDE= -DHAVE_CONFIG_H -I/usr/include/python2.2 -I/usr/lib/python2.2/config
PYTHON_LIB    = /usr/lib/python2.2/config
SO=		.so

#LDSHARED=	gcc -bundle -undefined suppress -flat_namespace
LDSHARED=	g++ -bundle -bundle_loader /usr/bin/python

CFILES= lisasim-lisa.cpp lisasim-noise.cpp lisasim-tdi.cpp lisasim-tdinoise.cpp lisasim-tens.cpp lisasim-wave.cpp

OFILES= $(CFILES:.cpp=.o)

lisasim-swig_wrap.cxx: lisasim-swig.i
	$(SWIG) -c++ -python lisasim-swig.i

lisasim-swig_wrap.o: lisasim-swig_wrap.cxx
	g++ -c lisasim-swig_wrap.cxx -I$(PYTHON_INCLUDE)

swig: $(OFILES) lisasim-swig.o lisasim-swig_wrap.o
	$(LDSHARED) $(OFILES) lisasim-swig.o lisasim-swig_wrap.o -lstdc++ -L$(PYTHON_LIB) -o _lisaswig$(SO)
	cp lisaswig.py lisaswig.pyc _lisaswig.so python/.

lisatest-tdi: $(OFILES) lisatest-tdi.o
	g++ $(OFILES) lisatest-tdi.o -o lisatest-tdi
