# main lisasim Makefile
# Michele 12/18/2003

all: lisaswig lisasim doc

clean: lisaswig-clean lisasim-clean objects-clean doc-clean

# standard compilation options
OPTIONS = -Wall -O3

# profiling compilation options
# OPTIONS = -Wall -O3 -pg -g 

# debugging compilation options
# OPTIONS = -Wall -g 

CURDIR=$(shell pwd)

# detect platform and setup directories
PLATFORM=$(shell uname)

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

# basic compilation step
$(BUILD)/%.o: %.cpp lisasim.h
	g++ $(OPTIONS) -c $< -o $@

# clean objects
objects-clean:
	rm $(OFILES)

# --- build lisaswig ---

# load some Makefile directives taken from the swig example directory

include ../swig/$(PLATFORM)/Makefile.python

lisaswig: $(LIBS)/_lisaswig$(PYTHON_SO)

SWIG= ../swig/$(PLATFORM)/bin/swig

# build the swig wrapper from the swig interface
# this is really a multiple target operation; it also builds $(LIBS)/lisasim.py

$(BUILD)/lisasim-swig_wrap.cxx $(LIBS)/lisasim.py: lisasim-swig.i $(LIBS)
	$(SWIG) -c++ -python lisasim-swig.i
	mv lisasim-swig_wrap.cxx $(BUILD)
	mv lisaswig.py $(LIBS)

# compile the swig wrapper in the build directory
# need -Wno-long-double for OS X? It is found in CFLAGS within swig Makefile
# will need to include also arrayobject.h, presumably from the Numeric distribution (currently from .)
$(BUILD)/lisasim-swig_wrap.o: $(BUILD)/lisasim-swig_wrap.cxx
	g++ $(OPTIONS) -c $(BUILD)/lisasim-swig_wrap.cxx -I. -I$(PYTHON_INCLUDE) -o $(BUILD)/lisasim-swig_wrap.o

# assemble everything in the lisaswig library
$(LIBS)/_lisaswig$(PYTHON_SO): $(OFILES) $(BUILD)/lisasim-swig.o $(BUILD)/lisasim-swig_wrap.o
	$(CXXSHARED) $(OFILES) $(BUILD)/lisasim-swig.o $(BUILD)/lisasim-swig_wrap.o $(IOBJS) $(PYTHON_DLNK) $(CPP_DLLIBS) -o $(LIBS)/_lisaswig$(PYTHON_SO)

lisaswig-clean:
	rm $(BUILD)/lisasim-swig.o $(BUILD)/lisasim-swig_wrap.o
	rm $(LIBS)/lisaswig.py $(LIBS)/lisaswig.pyc $(LIBS)/_lisaswig$(PYTHON_SO)

# --- build lisasim tool, including tdifast ---

lisasim: $(BINS) $(BUILD) $(OFILES) $(BUILD)/lisasim-main.o $(BUILD)/lisasim-tdifast.o
	g++ $(OPTIONS) $(OFILES) $(BUILD)/lisasim-main.o $(BUILD)/lisasim-tdifast.o -o $(BINS)/lisasim

lisasim-clean:
	rm $(BUILD)/lisasim-main.o $(BUILD)/lisasim-tdifast.o $(BINS)/lisasim

.PHONEY: clean

# --- testing ---

lisatest-tdi: $(OFILES) lisatest-tdi.o
	g++ $(OFILES) lisatest-tdi.o -o lisatest-tdi

# --- package docs from literate programming ---

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

$(DOCBUILD)/%.ctex: $(DOCBUILD) %.cpp
	$(BINS)/comments.py $< > $@

$(DOCBUILD)/%.htex: %.h
	$(BINS)/comments.py $< > $@

$(DOCDIR)/program.pdf: $(DOCDIR) program.tex $(TEXFILES)
	cp -p program.tex $(DOCBUILD)
	cd $(DOCBUILD)
	pdflatex program.tex; pdflatex program.tex
	mv $(DOCBUILD)/program.pdf $(DOCDIR)
	cd $(CURDIR)

doc-clean:
	rm -rf $(DOCBUILD)/*.htex $(DOCBUILD)/*.ctex $(DOCBUILD)/*.aux $(DOCBUILD)/*.log $(DOCDIR)/*.pdf
