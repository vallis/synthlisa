# detect platform, set directories

PLATFORM=$(shell uname -s)-$(shell uname -r)

# get current directory

CURDIR=$(shell pwd)

# get the top directory

TOPDIR=$(shell pwd)
PLADIR=$(TOPDIR)/$(PLATFORM)

all: contribs lisaswig

# create contribs

contribs:
	cd contrib-source; make

# create Synthetic LISA

lisaswig:
	cd lisasim; make all

