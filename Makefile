# detect platform, set directories

MAKE=make

PLATFORM=$(shell uname -s)-$(shell uname -r)

# get current directory

CURDIR=$(shell pwd)

# get the top directory

TOPDIR=$(shell pwd)
PLADIR=$(TOPDIR)/$(PLATFORM)

all: contribs lisaswig

# create contribs

contribs:
	cd contrib-source; $(MAKE)

# create Synthetic LISA

lisaswig:
	cd lisasim; $(MAKE) all

