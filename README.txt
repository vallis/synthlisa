====================================
Synthetic LISA, v. 1.2.3, 2005-07-22
====================================

by M. Vallisneri and J. W. Armstrong
Jet Propulsion Laboratory, Caltech
4800 Oak Grove Dr., Pasadena CA 91109

contact: Michele Vallisneri, vallis@vallis.org

---------------------------------------------
In this file:

0. Quickstart
1. Introduction to Synthetic LISA
2. License terms
3. System requirements
4. Installation of source package
   4a. Compilation notes for Mac OS X
   4b. Compilation notes for Linux
   4c. Compilation notes for Windows (cygwin)
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

- I assume you have Python >= 2.3, GNU gcc, GNU Make, and Gnuplot (if
  you want to plot). If you use Cygwin, you also need the Cygwin SWIG
  package.

- Unpack the Synthetic LISA distribution (synthLISA-1.2.1.tar.gz) and cd
  to it; download recent versions of Numeric and SWIG from
  www.vallis.org/syntheticlisa into the unpacked Synthetic LISA
  distribution directory.

- Run "./default-install.sh" ("./cygwin-install.sh" if you use
  Cygwin). This will install Numeric, SWIG, and Synthetic LISA
  locally, into "lib", "bin", "include", and "share" within the
  Synthetic LISA distribution directory. If you don't get any fatal
  errors, you're set. Otherwise, I'm afraid you'll have to read all
  the installation instructions below.

- If you're running csh/tcsh, do "source bin/synthlisa-setdir.csh"; if you're
  running bash, do "source bin/synthlisa-setdir.sh". This will set the correct
  path for Numeric and Synthetic LISA. [The path needs to be set each
  time you open a new shell session and want to use Synthetic LISA.]

- Now you can go to the "examples" directory, run a few of the Python
  (*.py) scripts ("python scriptname.py") and plot their results
  ("gnuplot -persist scriptname.plt").

- Read the tutorial in Section 7, the documentation (synthlisa.pdf),
  and the user discussions and feedback on the Synthetic LISA
  distribution website (which is being set up...)

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
could send me a signed copy of the license.

======================
3. System requirements
======================

The Synthetic LISA library was developed for (and is known to work on)
UNIX platforms that include the GNU gcc compilers, GNU make, and
Python >= v2.3. If you are running Windows (any version), you can set
up a reasonable Unix platform using Cygwin (www.cygwin.com). The
Synthetic LISA scripts included in the directory "examples" are
accompanied by Gnuplot scripts to plot their results, although any
plotting package can certainly be substituted with little work. This
said, it should be possible to install Synthetic LISA on any platform
with a working installation of Python and of a standards-compliant
C/C++ compiler.

Synthetic LISA relies on the Python package Numerical Python
("Numeric") to manage large arrays. Numeric can be downloaded from
numeric.scipy.org; Synthetic LISA is known to work with versions >=
23.1.

In addition, if you want to modify the Synthetic LISA code you will
need SWIG to regenerate Python wrappers. SWIG can be downloaded from
www.swig.org; Synthetic LISA is known to compile with versions >=
1.3.20.

If you need the GNU tools, Python, or Gnuplot, please see

- www.gnu.org
- www.python.org
- www.gnuplot.info

=================================
4. Installation of Synthetic LISA
=================================

To install and use Synthetic LISA from the source distribution
(synthLISA-1.2.1.tar.gz) you need a working installation of Python and
of a Python-interoperable C/C++ compiler, preferably gcc. The setup of
these is not discussed in this document.

--------------------------------
4a.Installing Numeric (required)
--------------------------------

Numeric might already be installed on your system, as you can check by
doing

> python
> import Numeric

No response is a good response. If Numeric is not available, you'll
get "ImportError: No module named Numeric".

You might also be able to install Numeric system-wide using a standard
distribution tool for your system, such as rpm, apt-get, or Fink. In
either case, skip the instructions to follow, although you might still
want to reinstall Numeric if the current version is old or incomplete.

Download the most recent platform-independent (Numeric-**.*.tar.gz)
version of Numeric from numeric.scipy.org, or directly from
SourceForce.

http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=1351

Unpack the distribution in a directory of your choice (you might have
a different version number).

> tar ztvf Numeric-23.8.tar.gz

At this point, you have two choices: if you have administrator access
to your workstation, you can install Numeric to the system-wide Python
"site-packages" directory; if not, or if you prefer not to modify the
system-wide directories, you can install locally to your home
directory, or to another directory.

Numeric is installed system-wide by cd'ing to the unpacked directory
and running

> python setup.py install

To get sufficient privileges for the install you may need to work as
root (use "su"), or you may just need to prefix the "python setup.py
install" command by "sudo"; the latter is the correct procedure on OS
X.

Numeric is installed locally by running

> python setup.py install --prefix=[userdir]

where you should replace [userdir] by the directory where you want to
install: it could be "/home/yourname", or
"/home/yourname/numeric". The setup process will create directories
"lib" and "include" within [userdir].

All done, but if you encounter any problems with this process (I hope
not!)  you should refer to the documentation at
numeric.scipy.org. More generally, a description of the python
"setup.py" (distutils) installation process can be found at

http://www.python.org/doc/current/inst/inst.html

--------------------------------------------------------------------
4b.Installing SWIG (required to modify and recompile Synthetic LISA)
--------------------------------------------------------------------

SWIG might already be installed on your system, as you can check by
trying

> swig -version

If you get an error, you probably don't have swig.

You might also be able to install SWIG system-wide using a standard
distribution tool for your system, such as rpm, apt-get, Fink, or
Cygwin's setup.exe (if you are using Cygwin, *definitely* install SWIG
using setup.exe; I have not been able to compile it from scratch.)  In
either case, skip the instructions to follow, although you might still
want to reinstall SWIG if the current version is old or incomplete.

Download the most recent platform-independent (swig-*.*.*.tar.gz)
version of SWIG from www.swig.org, or directly from

http://www.swig.org/download.html

Unpack it in a directory of your choice (you might have a different
version number).

> tar ztvf swig-1.3.24.tar.gz

As for Numeric, you have two choices: system-wide installation (files
will go to /usr/local), or local installation.

System-wide installation is achieved by cd'ing to the unpacked SWIG
directory and running

> ./configure
> make
> make install

after securing the required privileges. Local installation is achieved
by running

> ./configure --prefix=[userdir]
> make
> make install

where you should replace [userdir] by the directory where you want to
install (it could "/home/yourname", or "/home/yourname/numeric"). The
setup process will create directories "bin" and "share" within
[userdir].

All done, but if you encounter any problems with this step (I sure
hope not!) you should refer to the documentation at www.swig.org.

----------------------------
4c.Installing Synthetic LISA
----------------------------

If you're reading this file, you have already unpacked the Synthetic
LISA distribution. Again, you have the choice of installing
system-wide or locally. In addition, if you have installed Numeric
locally, you'll need to tell the Synthetic LISA installer where it can
be found.

So we have four cases:

- System-wide Synthetic LISA with system-wide Numeric: run

> python setup.py install

- System-wide Synthetic LISA with local Numeric installed in
  [numericdir]: run

> python setup.py install --with-numeric=[numericdir]

- Local Synthetic LISA with system-wide Numeric: run

> python setup.py install --prefix=[synthlisadir]

- Local Synthetic LISA with local Numeric installed in [numericdir}: run

> python setup.py install --prefix=[synthlisadir] --with-numeric=[numericdir]

In the last two options, you should replace [synthlisadir] by the
directory where you want to install. This could be "/home/yourname",
or it might be convenient to install to the very directory where you
have unpacked Synthetic LISA, using "--prefix=`pwd`").

In the unlikely event (ehm) that you should get compilation errors
(probably due to differences between your compiler or standard C++
libraries and those used in developing Synthetic LISA), please see
section 9 ("Asking for help") and I'll try to sort them out. Even
better, if you are able to sort them out, I'd love to hear about it,
to improve future distributions.

After compilation is done, set the Python path as described in the
next section, and you're ready to go. If you need to recompile
Synthetic LISA after modifying it, see section 4e.

--------------------------
4d.Setting the Python path
--------------------------

To use Numeric and Synthetic LISA, the Python interpreter needs to
know where they can be found. If you have performed system-wide
installs, this will be the case automatically. But if you have
performed local installs, you will need to change the value of the
system variable PYTHONPATH.

For your convenience, the Synthetic LISA installs creates scripts to
do this, which can be found in the "bin" directory of the Synthetic
LISA installation. If you're running csh/tcsh, do "source
bin/synthlisa-setdir.csh"; if you're running bash, do "source
bin/synthlisa-setdir.sh".

The PYTHONPATH needs to be set each time you open a new shell session
and want to work with Numeric or Synthetic LISA. If you want this job
to be done automatically, you can call the correct synthlisa-setdir
script to your .cshrc, .profile, or .login file.

---------------------------------
4e.Recompiling and Synthetic LISA
---------------------------------

All the Synthetic LISA files that you may want to modify are in the
subdirectory "lisasim". Unless you are changing only the C/C++ LISA
sources (and not the headers, or the SWIG interface
"lisasim/lisasim-swig.i"), you will need to have installed SWIG (as
described in section 4b, or by other means).

If SWIG was installed system-wide, and you can access it by simply
typing "swig", then you can just use the "python setup.py" commands
given above.

Otherwise, you will need to tell "setup.py" where to find the swig
executable, by adding "--with-swig=[swigdir]/bin/swig" to the "python
setup.py" commands given above. Yes, you must pass the fully qualified
swig executable name, not just the swig directory (I have my
reasons...).

If you do not want to install yet, but just build, replace "install"
with "build" in the "python setup.py" commands.

========================
5. System-specific notes
========================

-------------------------------------
5a.Mac OS X
-------------------------------------

First of all, at this time LISA is known to work correctly with
Panther (OS X 10.3) and Tiger (OS X 10.4).

Panther has a recent enough version of Python; gcc is included in the
Apple Development Tools (on your OS X CD, or already installed); you
can get gnuplot from fink.sourceforge.net, or from

http://www.lee-phillips.org/info/Macintosh/gnuplot.html

and probably in other places, as well; if X11 is not installed, you can get it from

http://www.apple.com/downloads/macosx/apple/x11formacosx.html

--------
5b.Linux
--------

The package compiles fine under RedHat Enterprise Linux WS 4, which
can be installed (and usually is) with X, Python and gnuplot. There
should be no problems also with RedHat Enterprise Linux ES 4 and
Fedora Core 3, which are very similar.

The packages compiles fine also under Debian Sarge (which should
become 3.1).

Other recent versions of Linux should also have no problems, as long
as Python >= 2.3. (For instance, RedHat 9.0 has Python 2.2.2, which
needs to be upgraded to work with Synthetic LISA.)

------------------------------------------
4c. Compilation notes for Windows (cygwin)
------------------------------------------

Synthetic LISA is known to work correctly with Cygwin 1.5.16-1, under
Windows ME (but I don't see why other flavors of Windows should have
more problems.) You will need a reasonably complete installation of
Cygwin, including certain standard utilities (tar, sharutils,
binutils, autoconf), the GCC compilers, Python, and SWIG (you should
use the version obtained through the Cygwin setup.exe rather than
trying to compile your own).

I wouldn't hope to get a stellar performance from your processor under
Windows/Cygwin. But it works...

Also, do yourself a favor:  escape the strictures of the Windows shell
by installing rxvt (through Cygwin's setup.exe), and configuring it as
per

http://c2.com/cgi/wiki?BetterCygwinTerminal

Under rxvt, it's even possible to run Cygwin's emacs (but first do "export
TERM=rxvt").

Good luck, and if you are more of a Windows wizard than I am, perhaps
you can give me some suggestions, or even compile Synthetic LISA under
ActiveState Python and some other compiler than GCC...

=======================
6. Running the examples
=======================

The best way to learn about Synthetic LISA is to study the example
scripts, found in the subdirectory "examples" in the Synthetic LISA
distribution directory. The files ending in ".py" are the Python
scripts that run the simulations; these will produce output (usually
noise or signal spectra) as ASCII ".txt" files, in the "data"
directory; the files ending in ".plt" are gnuplot scripts that will
plot the results; while the files ending in ".eps", in the "eps"
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
SyntheticLISA can be found in the documents "synthlisa.pdf" and
"manual.pdf" included with the package.

Study the example scripts to learn about the basic syntax of a
Synthetic LISA session. A good sequence of examples to examine is the
following:

- test-proofnoise.py
- test-tdiequal-X.py
- test-tdiequal.py
- test-tdibadmass.py

You will see that after some standard library-loading incantation,

> from synthlisa import *
> from Numeric import *

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

> noiseX = getobsc(16384,0.1,originalTDI.X)

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

Although I regard the current version of Synthetic LISA as a very
usable and useful product in its current version, I am aware that
its eventual impact is strongly dependent on the involvement of users.

There are many aspects that users can help (and help themselves) with:

- Pointing out bugs so that I can fix them

- Requesting new features and/or improvements

- Providing compilation tricks and fixes for additional systems

To this purpose, a Synthetic LISA webpage is being set up on the Open
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

- Support for the LISA Mock Data Format

Throughout these changes, I shall strive to keep the Synthetic LISA
interface backward compatible, so that your scripts and applications
will continue to be useful as you updgrade Synthetic LISA.
