# detect platform, set directories

PLATFORM=$(shell uname -s)-$(shell uname -r)

# get current directory

CURDIR=$(shell pwd)

# get the top directory

TOPDIR=$(shell cd ..; pwd)
PLADIR=$(TOPDIR)/$(PLATFORM)

# build and lib directories

BUILD=$(CURDIR)/build/$(PLATFORM)
LIBS=$(PLADIR)/lib

# standard compilation options

OPTIONS = -Wall -O3

# these would be profiling and debugging compilation options

# OPTIONS = -Wall -O3 -pg -g 
# OPTIONS = -Wall -g 

# so, what do you want to build today?

all: lisaswig lisautils

# the copying of lisautils.py
lisautils: $(LIBS)/lisautils.py

$(LIBS)/lisautils.py: lisautils.py
	cp lisautils.py $(LIBS)/.

# get ready to compile basic lisasim objects
CFILES= lisasim-lisa.cpp lisasim-noise.cpp lisasim-tdi.cpp lisasim-tdisignal.cpp lisasim-tdinoise.cpp lisasim-tens.cpp lisasim-wave.cpp

# each of these cpp files has a header
HFILES= $(CFILES:.cpp=.h)

# lisasim.h depends on all the headers
lisasim.h: $(HFILES)

# we build the object files in the build directory
OFILES= $(foreach cfile,$(CFILES),$(BUILD)/$(cfile:.cpp=.o))

# basic compilation step
$(BUILD)/%.o: %.cpp $(HFILES)
	g++ $(OPTIONS) $(GSL_INCLUDE) -c $< -o $@

# make the directories if they are not there
$(BUILD):
	mkdir -p $(BUILD)
$(LIBS):
	mkdir -p $(LIBS)

# load some Makefile directives taken from the swig example directory
include $(PLADIR)/include/Makefile.python

# this is what we have to build
lisaswig: $(BUILD) $(LIBS) $(LIBS)/_lisaswig$(PYTHON_SO) $(LIBS)/lisaswig.py

# where can we find swig?
SWIG= $(PLADIR)/bin/swig -w402

# build the swig wrapper from the swig interface
# this is a multiple target operation; it also builds $(LIBS)/lisaswig.py
$(BUILD)/lisasim-swig_wrap.cxx $(LIBS)/lisaswig.py: lisasim-swig.i lisasim-typemaps.i lisasim.h
	$(SWIG) -c++ -python lisasim-swig.i
	mv lisasim-swig_wrap.cxx $(BUILD)
	mv lisaswig.py $(LIBS)

# where can we find Numeric includes?
NUMERIC_INCLUDE=-I$(PLADIR)/include/$(shell ls $(PLADIR)/include | grep python2 | head -1)/Numeric

# ah, the GSL libraries!
GSL_DIR=$(TOPDIR)/contrib-source/GSL-1.4
GSL_INCLUDE=-I$(GSL_DIR)
GSL_BUILD=$(GSL_DIR)/build/$(PLATFORM)

# what do we need from GSL?
GSL_OFILES= $(GSL_BUILD)/rng.o $(GSL_BUILD)/ranlux.o $(GSL_BUILD)/mt.o

$(GSL_OFILES):
	make -C $(GSL_DIR) libs

# compile the swig wrapper in the build directory
$(BUILD)/lisasim-swig_wrap.o: $(BUILD)/lisasim-swig_wrap.cxx
	g++ $(OPTIONS) -c $(BUILD)/lisasim-swig_wrap.cxx -I. $(NUMERIC_INCLUDE) $(PYTHON_INCLUDE) $(GSL_INCLUDE) -o $(BUILD)/lisasim-swig_wrap.o

# assemble everything in the lisaswig library
$(LIBS)/_lisaswig$(PYTHON_SO): $(OFILES) $(BUILD)/lisasim-swig.o $(BUILD)/lisasim-swig_wrap.o $(GSL_OFILES)
	$(CXXSHARED) $(OFILES) $(GSL_OFILES) $(BUILD)/lisasim-swig.o $(BUILD)/lisasim-swig_wrap.o $(IOBJS) $(PYTHON_DLNK) $(CPP_DLLIBS) -o $(LIBS)/_lisaswig$(PYTHON_SO)
