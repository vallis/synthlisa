==============================
Synthetic LISA, v. 1.0 (alpha)
==============================

by Michele Vallisneri and John Armstrong
Jet Propulsion Laboratory, Caltech
4800 Oak Grove Dr., Pasadena CA 91109

contact: vallis@vallis.org

---------------------------------------------
In this file:

0. Quickstart
1. Introduction to Synthetic LISA
2. License terms
3. System requirements
4. Installation of source package
   4a. Compilation notes for Mac OS X
   4b. Compilation notes for Linux RedHat 9.0
   4c. Compilation notes for SunOS 5.8
5. Installation of binary package
6. Running the examples
7. Synthetic LISA usage
8. User involvement
9. Asking for help 
10.To do list
---------------------------------------------

=============
0. Quickstart
=============

- Unpack ("tar zxvf") the distribution, do "make all" (if you are
  compiling from sources).

- Run "source setdir.sh" or "source setdir.csh" (depending on your
  shell; you'll need to do this every time you want to work with
  Synthetic LISA).

- Go to the example directory, run a few of the Python (*.py) scripts
  ("python scriptname.py") and plot their results ("gnuplot -persist
  scriptname.plt").

- Read the tutorial in Section 7, the documentation (synthlisa.pdf),
  and the user discussions and feedback on the Synthetic LISA
  distribution website (OCF website, TBD).

This was quick, what did you expect? For more detail, read on.

===============
1. Introduction
===============

Synthetic LISA, developed by Michele Vallisneri and John Armstrong at
the Jet Propulsion Laboratory under the auspices of the LISA Mission
Science Office, is a C++/Python package to simulate the LISA science
process at the level of scientific and technical requirements. It
generates synthetic time series of the LISA fundamental noises, as
filtered through all the TDI observables; it provides a streamlined
module to compute the TDI responses to gravitational waves, according
to a full model of TDI (including the motion of the LISA array, and
the temporal and directional dependence of the armlengths).

================
2. License terms
================

Please see license.pdf or license.doc. I would appreciate it if you
could return me a signed license, if you have not done so yet.

======================
3. System requirements
======================

The Synthetic LISA library was developed for UNIX platforms that
provide the GNU C/C++ compilers (gcc), GNU make, and the standard GNU
utilities (tar, grep, sed, etc.) In addition, the Python interface
requires Python >= v2.2, and the example scripts use gnuplot (and of
course X11) to plot spectra. If you download a binary distribution for
your system, you will not need the compilers, but you will still need
Python and gnuplot.

All this software is rather standard for modern scientific
workstations; ask you system administrator about it, or download it
from the sites

- www.gnu.org
- www.python.org
- www.gnuplot.info

The Synthetic LISA Python interface was built using SWIG
(www.swig.org), and it relies on the Python package Numeric
(numpy.sourceforge.net) to manage large arrays. A recent release of
Numeric is included in the binary and source distributions; SWIG is
needed only to compile Synthetic LISA from sources, and it is included
in the source distribution. (If you already have Numeric and/or SWIG
in your system, Synthetic LISA might or might not be able to work with
them, depending on their version; if in doubt, use the versions
provided).

=================================
4. Installation of source package
=================================

To compile Synthetic LISA from sources, you will need to unpack the
tar.gz archive in a directory of your choice,

> cd mydirectory
> tar zxvf synthlisa-source-version.tar.gz

Then run

> make all

this will first create the Python packages SWIG and Numeric from the
tar archives contained in contrib-source (don't worry about compiler
warnings), and then compile the Synthetic LISA library.

In the unlikely event (ehm) that you should get compilation errors
(probably due to differences between your compiler or standard C++
libraries and those used in developing Synthetic LISA), please see
section 9 ("Asking for help") and we'll try to sort them out. Even
better, if you are able to sort them out, we'd love to hear about it,
to improve future versions.

After compilation is done, set the Python path as described in the
next section, and you're ready to go.

---------------------------------
4a.Compilation notes for Mac OS X
---------------------------------

This software was developed under Mac OS X, so we thought we'd give a
few pointers to other Mac addicts. We're giving you a binary
distribution, but it's always good to be able to recompile from
scratch.

First, a couple of bugs that you need to fix (you'll need
administrator privileges: then "sudo" your commands to get access to
protected files). These are based on the standard distribution of
Development Tools from Apple (on some computers, obtained by clicking
the "pkg" within the "Applications/Installers" directory), so other
distributions, such as fink, might not have them.

- In Jaguar, there is a bug in
  /usr/include/gcc/darwin/3.1/g++-v3/cmath that will undo the
  definition of isnan. To fix it, comment out the line #undef
  isnan. In Panther, the file to change is
  /usr/include/gcc/darwin/3.3/c++/cmath.

- In Jaguar, malloc.h does not exist (and packages will go looking for
  it). "touch" it in your compilation directory to create a dummy,
  empty version. Panther seems to have the file.

- In Jaguar, there is a bug in /usr/lib/python2.2/config/Makefile,
  where the architecture flag "LDFLAGS= -arch i386 -arch ppc" should
  be replaced with "LDFLAGS= -arch ppc". This bug might prevent the
  compilation of some Python libraries (although it does not seem to
  interfere with the current versions of Swig and Numeric).

If you're wondering where to find the required software, X11 is
distributed by Apple for Jaguar, and it is standard in Panther; Python
gets installed with Apple's Development tools; and you can get gnuplot
(and many other things) from fink.sourceforge.net, or from
http://www.lee-phillips.org/info/Macintosh/gnuplot.html, and probably
in other places, as well.

------------------------------
4b.Compilation notes for Linux
------------------------------

The package compiles fine under RedHat 9.0, which can be installed
(and usually is) with X, Python and gnuplot. Other recent versions of
Linux should also have no problems, and might be compatible with our
binaries.

One compilation caution: it seems that setdir.sh and setdir.csh are
broken by colored "ls" output. Just do

alias ls="ls --color=none" (bash)
alias ls="ls --color=none" (csh)

before sourcing setdir.(c)sh, or modify setdir.sh and setdir.csh
accordingly.

----------------------------------
4c.Compilation notes for SunOS 5.8
----------------------------------

A similar procedure might apply to other close versions of SunOS.

First of all, the default "make" and "tar" on Suns is usually not
GNU. You should change the "MAKE=make" and "GNUTAR=tar" lines in the
synthlisa/Makefile, synthlisa/contrib-source/Makefile, and
synthlisa/lisasim/Makefile to reflect the location of your GNU make
and tar.

Second, on Suns grep does not seem to have the "-A" option, and I have
not been able to find a fix so far. Replace the "GREPA=grep -A 1" line
in synthlisa/contrib-source/Makefile with "GREPA=grep"; this, however,
will break one of the Makefiles. To correct the problem, do first

make contribs

and then edit the file SunOS-5.8/include/Makefile.python, replacing
the "CPP_DLLIBS..." line with

CPP_DLLIBS = -L/usr/local/lib -lstdc++ -lgcc

Then you can do

make lisaswig

and you should be all set.

=================================
5. Installation of binary package
=================================

If your platform is among those for which a binary distribution is
provided and you wish to use it, you will need to unpack the tar.gz
archive in a directory of your choice,

> cd mydirectory
> tar zxvf synthlisa-binary-platform-version.tar.gz

[do not type ">", which denotes the command prompt, and replace
"mydirectory", "platform", and "version" as appropriate]. Before each
Synthetic LISA session, you will need to set the Python path by going
to the main "synthlisa" directory, and running

> source setdir.csh (if you use csh or tcsh)

or

> source setdir.sh (if you use bash or sh)

You can check that the installation works correctly by entering the
Python interactive shell,

> python

and then, at the Python prompt, typing

python> from lisaswig import *

[again, "python>" denotes the Python prompt]. If Python does not
complain, you're all set.

If the library cannot be found, it might be that the identification
string for your operating system, given by the shell command

> echo `uname -s`-`uname -r`

(for instance, "Linux-2.4.20-24.9") is slightly different from the
library directory, of similar name, within the main synthlisa
directory. You can try renaming that directory to the identification
string for your system, but if you keep having problems, it might be
that the binary distribution is incompatible with your system. In that
case, you should compile Synthetic LISA from sources (and then perhaps
make the result available to the community; see section 8).

=======================
6. Running the examples
=======================

The best way to learn about Synthetic LISA is to study the example
scripts, found in the directory "examples". The files ending in ".py"
are the Python scripts that run the simulations; these will produce
output (usually noise or signal spectra) as ASCII ".txt" files, in the
"data" directory; the files ending in ".plt" are gnuplot scripts that
will plot the results; and the files ending in ".eps", in the "eps"
directory, are the resulting postscript plots (view them with
ghostview or a similar utility).

Run the Python scripts by typing

> python example.py

[replace "example" as needed]; produce the plots by typing

> gnuplot -persist example.plt

["-persist" will keep the graph in a window until you type "q"].

If the simulations run too slowly on your system, you may want to
reduce the value of the variable "samples" found near the top of each
script. If they run too fast (and you want better detail), increase
"samples".

=======================
7. Synthetic LISA usage
=======================

A general description of the formulation, implementation, and usage of
SyntheticLISA can be found in the document "synthlisa.pdf" included
with the package.

Study the example scripts to learn about the basic syntax of a
Synthetic LISA session. A good sequence of examples to examine is the following:

- test-proofnoise.py
- test-tdiequal-X.py
- test-tdiequal.py
- test-tdibadmass.py

You will see that after some standard
library-loading incantations,

> from lisaswig import *
> from Numeric import *
> from lisautils import *

the script goes into creating a hierarchical structure of objects: for
instance, to get time series of a TDI noise we first need to create a
LISA geometry object,

> originallisa = OriginalLISA(16.6782,16.6782,16.6782)

where the "16.6782" are the LISA armlengths in seconds for a
stationary LISA; then we use the LISA geometry object to create a
TDInoise object,

> originalTDI = TDInoise(originallisa,
>                        1.0, 2.5e-48, # proof-mass noise parameters
>                        1.0, 1.8e-37, # optical-path noise parameters
>                        1.0, 1.1e-26) # laser frequency noise parameters

[the "#" denote comments, and the two parameters of each noise type
are a sampling time, in seconds, and a spectral density, as described
in synthlisa.pdf]. Next, we can get (say) 16384 values of the TDI
observable X, at times starting from 0, and separated by 0.1 seconds:

> noiseX = getobs(16384,0.1,originalTDI.X)

Now "noiseX" is a Numeric array of length 16384, which can be written
to disk as an ASCII file,

> writearray('myXnoise.txt',noiseX)

or which can be written as a binary file (of "double" floating-point
numbers, with the byte ordering of the platform where Synthetic LISA
is being run),

> writebinary('myXnoise.bin',noiseX)

You can also use the function "spect" to compute a triangle-windowed,
averaged spectrum of "noiseX":

> specX = spect(noiseX,0.1,32)

[we need to tell "spect" the sampling timestep 0.1 (s), and the number
of periods (here 32) that the signal should be divided into to compute
the averaged spectrum]. The result is a 2D Numeric array, where each
element is a pair of frequency (in seconds), and noise density (in
1/Hertz).

===================
8. User involvement
===================

Although we regard the current version of Synthetic LISA as a very
usable and useful product in its current version, we are aware that
its eventual impact is strongly dependent on the involvement of users.

There are many aspects that users can help (and help themselves) with:

- Pointing out bugs so that we can fix them

- Requesting new features and/or improvements

- Providing binary distributions and/or compilation tricks and fixes
  for additional systems

To this purpose, we have set up the Synthetic LISA webpage on the Open
Channel Foundation website (www.openchannelfoundation.org) with a user
forum, bug report forms, and other tools for interaction. Be sure to
use them!

==================
9. Asking for help
==================

Please refer to the forthcoming Synthetic LISA distribution site at
the Open Channel Foundation (www.openchannelfoundation.org) to report
bugs, request features, ask for explanations, and generally interact
with other Synthetic LISA designers and with other users.

=============
10.To do list
=============

A partial list of planned improvements and upgrades:

- Improved documentation

- Inclusion of USO noise

- More Wave models

- Caching of the LISA and Noise classes for improved speed

- Support for the Numarray (the successor to Numeric)

- Support for the LISA Mock Data Format

Throughout these changes, we shall strive to keep the Synthetic LISA
interface backward compatible, so that your scripts and applications
will continue to be useful as you updgrade SynthLISA.
