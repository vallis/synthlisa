#OPTIONS = -Wall -g 
#OPTIONS = -Wall -O3 -pg -g 
OPTIONS = -Wall -O3 

CFILES= lisasim-lisa.cpp lisasim-noise.cpp lisasim-tdi.cpp lisasim-tdisignal.cpp lisasim-tdinoise.cpp lisasim-tens.cpp lisasim-wave.cpp
HFILES= $(CFILES:.cpp=.h)

OFILES= $(CFILES:.cpp=.o)

lisasim: $(OFILES) lisasim-main.o lisasim-tdifast.o
	g++ $(OPTIONS) $(OFILES) lisasim-main.o lisasim-tdifast.o -o lisasim

%.o: %.cpp lisasim.h
	g++ $(OPTIONS) -c $<

lisasim.h: $(HFILES)

clean:
	rm -rf *.o lisasim
	rm -rf lisasim-swig_wrap.cxx lisaswig.py lisaswig.pyc _lisaswig.so
	rm -rf literate/*.htex literate/*.ctex literate/*.aux literate/*.log literate/*.pdf literate/*.ps

.PHONEY: clean

# --- testing ---

lisatest-tdi: $(OFILES) lisatest-tdi.o
	g++ $(OFILES) lisatest-tdi.o -o lisatest-tdi

# --- swig ---

SWIG= swig/bin/swig

PYTHON_INCLUDE=	-DHAVE_CONFIG_H -I/usr/include/python2.2 -I/usr/lib/python2.2/config
PYTHON_LIB=	/usr/lib/python2.2/config
SO=	.so

#LDSHARED=	gcc -bundle -undefined suppress -flat_namespace
LDSHARED=	g++ -bundle -bundle_loader /usr/bin/python

lisasim-swig_wrap.cxx: lisasim-swig.i
	$(SWIG) -c++ -python lisasim-swig.i

lisasim-swig_wrap.o: lisasim-swig_wrap.cxx
	g++ -c lisasim-swig_wrap.cxx -I$(PYTHON_INCLUDE)

swig: $(OFILES) lisasim-swig.o lisasim-swig_wrap.o
	$(LDSHARED) $(OFILES) lisasim-swig.o lisasim-swig_wrap.o -lstdc++ -L$(PYTHON_LIB) -o _lisaswig$(SO)
	cp lisaswig.py _lisaswig.so python/.

# --- literate programming ---

TEXFILES= $(CFILES:%.cpp=literate/%.ctex) $(HFILES:%.h=literate/%.htex)

literate/%.ctex: %.cpp
	literate/comments.py $< > $@

literate/%.htex: %.h
	literate/comments.py $< > $@

literate/program.pdf: $(TEXFILES) literate/program.tex
	cd literate; pdflatex program.tex; pdflatex program.tex

doc: literate/program.pdf
