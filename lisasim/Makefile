# main lisasim Makefile
# Michele 12/18/2003

all: lisaswig lisasim lisatest doc

clean: lisaswig-clean lisasim-clean lisatest-clean objects-clean doc-clean

# standard compilation options
OPTIONS = -Wall -O3

# profiling compilation options
# OPTIONS = -Wall -O3 -pg -g 

# debugging compilation options
# OPTIONS = -Wall -g 

CURDIR=$(shell pwd)

# detect platform and setup directories
PLATFORM=$(shell uname -s)-$(shell uname -r)

BUILD=../build/$(PLATFORM)
BINS=../bin/$(PLATFORM)
LIBS=../lib/$(PLATFORM)

# make directories if they are not there
$(BUILD):
	mkdir -p $(BUILD)

$(LIBS):
	mkdir -p $(LIBS)

$(BINS):
	mkdir -p $(BINS)

# --- compile basic lisasim objects ---

# note no lisasim-main, lisasim-swig, lisasim-tdifast, or lisatest files
CFILES= lisasim-lisa.cpp lisasim-noise.cpp lisasim-tdi.cpp lisasim-tdisignal.cpp lisasim-tdinoise.cpp lisasim-tens.cpp lisasim-wave.cpp

# each of these cpp files has a header
HFILES= $(CFILES:.cpp=.h)

# we build the object files in the build directory
OFILES= $(foreach cfile,$(CFILES),$(BUILD)/$(cfile:.cpp=.o))

# if you include lisasim.h (which includes all the headers), then you should depend on them
lisasim.h: $(HFILES)

GSL_DIR=../contrib-source/GSL-1.4
GSL_INCLUDE=-I$(GSL_DIR)

# what do we need from GSL?

GSL_OFILES= $(BUILD)/rng.o $(BUILD)/ranlux.o $(BUILD)/mt.o

$(GSL_OFILES):
	make -C $(GSL_DIR) libs

# basic compilation step
$(BUILD)/%.o: %.cpp lisasim.h
	g++ $(OPTIONS) $(GSL_INCLUDE) -c $< -o $@

# clean objects
objects-clean:
	-rm -f $(OFILES)

# --- build lisaswig ---

# load some Makefile directives taken from the swig example directory

include ../swig/$(PLATFORM)/Makefile.python

lisaswig: $(BUILD) $(LIBS) $(LIBS)/_lisaswig$(PYTHON_SO)

SWIG= ../swig/$(PLATFORM)/bin/swig -w402

# build the swig wrapper from the swig interface
# this is really a multiple target operation; it also builds $(LIBS)/lisasim.py

$(BUILD)/lisasim-swig_wrap.cxx $(LIBS)/lisasim.py: lisasim-swig.i lisasim-typemaps.i lisasim.h
	$(SWIG) -c++ -python lisasim-swig.i
	mv lisasim-swig_wrap.cxx $(BUILD)
	mv lisaswig.py $(LIBS)

NUMERIC_INCLUDE=-I../numeric/$(PLATFORM)/include/$(shell ls ../numeric/$(PLATFORM)/include | head -1)/Numeric
NUMARRAY_INCLUDE=-I../numarray/$(PLATFORM)/include/$(shell ls ../numarray/$(PLATFORM)/include | head -1)/numarray

# compile the swig wrapper in the build directory
# need -Wno-long-double for OS X? It is found in CFLAGS within swig Makefile
# we are including arrayobject.h from the Numeric distribution
# we are including lisasim.h from .
$(BUILD)/lisasim-swig_wrap.o: $(BUILD)/lisasim-swig_wrap.cxx
	g++ $(OPTIONS) -c $(BUILD)/lisasim-swig_wrap.cxx -I. $(NUMERIC_INCLUDE) $(PYTHON_INCLUDE) $(GSL_INCLUDE) -o $(BUILD)/lisasim-swig_wrap.o

# assemble everything in the lisaswig library
$(LIBS)/_lisaswig$(PYTHON_SO): $(OFILES) $(BUILD)/lisasim-swig.o $(BUILD)/lisasim-swig_wrap.o $(GSL_OFILES)
	$(CXXSHARED) $(OFILES) $(GSL_OFILES) $(BUILD)/lisasim-swig.o $(BUILD)/lisasim-swig_wrap.o $(IOBJS) $(PYTHON_DLNK) $(CPP_DLLIBS) -o $(LIBS)/_lisaswig$(PYTHON_SO)

lisaswig-clean:
	-rm -f $(BUILD)/lisasim-swig_wrap.cxx
	-rm -f $(BUILD)/lisasim-swig.o $(BUILD)/lisasim-swig_wrap.o
	-rm -f $(LIBS)/lisaswig.py $(LIBS)/lisaswig.pyc $(LIBS)/_lisaswig$(PYTHON_SO)

# --- build lisasim tool, including tdifast ---

lisasim: $(BINS)/lisasim

$(BINS)/lisasim: $(BINS) $(BUILD) $(OFILES) $(BUILD)/lisasim-main.o $(BUILD)/lisasim-tdifast.o
	g++ $(OPTIONS) $(OFILES) $(GSL_OFILES) $(BUILD)/lisasim-main.o $(BUILD)/lisasim-tdifast.o -o $(BINS)/lisasim

lisasim-clean:
	-rm -f $(BUILD)/lisasim-main.o $(BUILD)/lisasim-tdifast.o
	-rm -f $(BINS)/lisasim

# --- build lisatest ---

lisatest: $(BINS)/lisatest-tdi

$(BINS)/lisatest-tdi: $(OFILES) $(BUILD)/lisatest-tdi.o
	g++ $(OPTIONS) $(OFILES) $(GSL_OFILES) $(BUILD)/lisatest-tdi.o -o $(BINS)/lisatest-tdi

lisatest-clean:
	-rm -f $(BUILD)/lisatest-tdi.o
	-rm -f $(BINS)/lisatest-tdi

# --- build package docs from literate programming ---

DOCDIR=../doc
DOCBUILD=../doc/build

TEXFILES= $(CFILES:%.cpp=$(DOCBUILD)/%.ctex) $(HFILES:%.h=$(DOCBUILD)/%.htex)

doc: $(DOCDIR)/program.pdf

$(DOCDIR):
	mkdir -p $(DOCDIR)

$(DOCBUILD):
	mkdir -p $(DOCBUILD)

# might use this step to set the path to python in comments.py
$(BINS)/comments.py: ../bin/comments.py
	cp -p ../bin/comments.py $(BINS)

$(DOCBUILD)/%.ctex: $(BINS)/comments.py $(DOCBUILD) %.cpp
	$(BINS)/comments.py $< > $@

$(DOCBUILD)/%.htex: $(BINS)/comments.py $(DOCBUILD) %.h
	$(BINS)/comments.py $< > $@

$(DOCDIR)/program.pdf: $(DOCDIR) program.tex $(TEXFILES)
	cp -p program.tex $(DOCBUILD)
	cd $(DOCBUILD); pdflatex --interaction batchmode program.tex; pdflatex --interaction batchmode program.tex; mv program.pdf ..; cd $(CURDIR)

doc-clean:
	-rm -rf $(DOCBUILD)/*.htex $(DOCBUILD)/*.ctex $(DOCBUILD)/*.aux $(DOCBUILD)/*.log $(DOCDIR)/*.pdf

.PHONEY: clean
