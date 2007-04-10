====================================
Synthetic LISA, v. 1.3.2, 2006-08-28
====================================

by M. Vallisneri and J. W. Armstrong
Jet Propulsion Laboratory, Caltech
4800 Oak Grove Dr., Pasadena CA 91109

contact: Michele Vallisneri, vallis@vallis.org

-------------------------------------------------
In this file:

-1.Changes, changes
 0. Quickstart
 1. Introduction to Synthetic LISA
 2. License terms
 3. System requirements
 4. Installation of Synthetic LISA
    4a. Installing numpy (required)
    4b. Installing SWIG (required to modify and
                         recompile Synthetic LISA)
    4c.Installing pyRXP (required)
    4d.Installing PyX (optional)
    4e.Installing matplotlib (optional)
    4f.Installing Synthetic LISA
    4g.Setting the Python path
    4h.Recompiling and Synthetic LISA
 5. System-specific notes
    5a. Mac OS X
    5b. Linux
    5c. Cygwin
 6. Running the examples
 7. Synthetic LISA usage
 8. User involvement
 9. Asking for help 
10. To do list
-------------------------------------------------

====================
-1. Changes, changes
====================

1.3.2 (2006-08-28)
------------------

- With this release, synthlisa is ported to numpy, the new standard array
  package for Python. This helps the interoperability of LISA with many
  modern Python packages. Changes to existing synthlisa scripts should be
  minimal. At worst, if you were doing "import Numeric", replace it with
  "import numpy.oldnumeric as Numeric". If you were doing "from Numeric
  import *", now do "from numpy import *".
 
- The standard getobs and getobsc functions used to collect arrays of TDI
  observables or noises have become much faster, thanks to their
  reimplementation in C++. This makes synthlisa 40% faster altogether. [But
  the new functions will be used only for "stock" observables and noises,
  and not for pseudoobservables written as lambda functions, which I have
  been known to do occasionally. This should hardly concern you...]
  
- The plugin architecture for MLDC sources has been improved somewhat. Now
  all XML definitions for the new sources are in the contrib directory.

=============
0. Quickstart
=============

- I assume you have Python >= 2.3, GNU gcc, GNU Make, and gnuplot (if you
  want to plot; an alternative is PyX, see 4d below, or matplotlib, see
  4e blow). If you use Cygwin, you also need the Cygwin SWIG package.

- Unpack the Synthetic LISA distribution (synthLISA-1.3.2.tar.gz) and
  cd to it.

- Run "./default-install.py". This will install numpy,
  SWIG, and Synthetic LISA locally in the Synthetic LISA
  distribution directory "synthLISA-1.3.2". If you don't get
  any fatal errors, you're set. Otherwise, I'm afraid you'll
  have to read all the installation instructions below.
  [Note: this quickstart procedure is currently untested on Cygwin.]

- Set the Python path by running "source bin/synthlisa-setdir.sh" (if
  your shell is bash) or "source bin/synthlisa-setdir.csh" (if your
  shell is csh or tcsh). You will need to do this at the beginning of
  every synthLISA session in order for Python to find the correct
  libraries.

- Now you can go to the "examples" directory, run a few of the Python
  (*.py) scripts ("python scriptname.py") and plot their results
  ("gnuplot -persist scriptname.plt").
  [Note: after the XML-ization of synthLISA output, the examples
  directory has become scarcely populated, as the example scripts
  are converted. Stay tuned.]

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
to a full model of TDI (including the motion of the LISA array, and the
temporal and directional dependence of the armlengths).

================
2. License terms
================

Synthetic LISA is currently licensed under the Caltech public domain
license (doc/license.pdf and doc/license.doc). If you have downloaded
synthLISA anonymously from the www.vallis.org webpage (i.e, you have
not inserted your name, institution, and e-mail), please e-mail me this
information so that I can keep track of the synthLISA user base.

Later, I hope to switch Synthetic LISA to the less restrictive MIT
license.

======================
3. System requirements
======================

The Synthetic LISA library was developed for (and is known to work on)
UNIX platforms that include the GNU gcc compilers, GNU make, and Python
>= v2.3. If you are running Windows (any version), you can set up a
reasonable Unix platform using Cygwin (www.cygwin.com). The Synthetic
LISA scripts included in the directory "examples" are accompanied by
Gnuplot scripts to plot their results, although any plotting package
can certainly be substituted with little work. This said, it should be
possible to install Synthetic LISA on any platform with a working
installation of Python and of a standards-compliant C/C++ compiler.

Synthetic LISA relies on a number of Python packages: the latest
version of these known to interact correctly with synthLISA is included
in the synthLISA distribution. (This does not mean that newer versions
won't work, just that I haven't tested them.) These packages include

- numpy (currently v. 1.0b1), see http://numpy.scipy.org

- pyRXP (v. 1.07), see http://www.reportlab.org/pyrxp.html

- SWIG (v. 1.3.29), see http://www.swig.org

- PyX (v. 0.8.1), see http://pyx.sourceforge.net

Of these, numpy is needed for synthLISA to deal with arrays in
Python; pyRXP is needed for XML input; SWIG is needed if you want to
modify the Synthetic LISA code, to regenerate Python wrappers; and PyX
is a stand-alone plotting package used for some of the example scripts.

If you need the GNU tools, Python, or Gnuplot, please see

- www.gnu.org
- www.python.org
- www.gnuplot.info

=================================
4. Installation of Synthetic LISA
=================================

To install and use Synthetic LISA from the source distribution
(synthLISA-1.3.2.tar.gz) you need a working installation of Python and
of a Python-interoperable C/C++ compiler, preferably gcc. The setup of
these is not discussed in this document.

The default-install.py script (see the quickstart in section 0 above)
will perform a complete installation of Synthetic LISA and all
attendant packages in the synthLISA unpacking (a.k.a. "distribution")
directory. The single Python packages (or the single synthlisa package)
can also be installed by running default-install.py with one or more of
the arguments "numpy", "pyRXP", "SWIG", "PyX", and "synthLISA".

The instructions below pertain to the case where you want to install
synthLISA and/or one or more packages *globally* in your system (i.e.,
in the main Python library location, shared by all users, which will
require that you have administrator privileges), or you want to
customize your installation in other ways.

------------------------------
4a.Installing numpy (required)
------------------------------

Numpy might already be installed on your system, as you can check by
doing

> python
> import numpy

No response is a good response. If numpy is not available, you'll
get "ImportError: No module named numpy".

You might also be able to install numpy system-wide using a standard
distribution tool for your system, such as rpm, apt-get, or Fink. In
either case, skip the instructions to follow, although you might still
want to reinstall numpy if the current version is old or incomplete.

Download the most recent platform-independent (numpy-*.tar.gz) version of
numpy from numpy.scipy.org. Unpack the distribution in a directory of your
choice (you might have a different version number).

> tar ztvf numpy-1.0b1.tar.gz

At this point, you have two choices: if you have administrator access
to your workstation, you can install numpy to the system-wide Python
"site-packages" directory; if not, or if you prefer not to modify the
system-wide directories, you can install locally to your home
directory, or to another directory.

numpy is installed system-wide by cd'ing to the unpacked directory
and running

> python setup.py install

To get sufficient privileges for the install you may need to work as
root (use "su"), or you may just need to prefix the "python setup.py
install" command by "sudo"; the latter is the correct procedure on OS
X.

numpy is installed locally by running

> python setup.py install --prefix=[userdir]

where you should replace [userdir] by the directory where you want to
install: it could be "/home/yourname", or "/home/yourname/numpy". The
setup process will create directories "lib" and "include" within
[userdir].

All done, but if you encounter any problems with this process (I hope
not!) you should refer to the documentation at numpy.scipy.org. More
generally, a description of the python "setup.py" (distutils)
installation process can be found at

http://www.python.org/doc/current/inst/inst.html

If numpy gives you problems related to the BLAS and LAPACK matrix
libraries, it may be that the versions of these libraries on your system
are not compatible with the use that numpy makes of them. While numpy
will use these libraries (for speed) if it finds them, it does not really need them. To compile numpy without the system BLAS and LAPACK, set the environment variables ATLAS, BLAS, and LAPACK to None before the python setup.py install step. In bash,

> export ATLAS=None; export BLAS=None; export LAPACK=None

and in tcsh,

> setenv ATLAS None; setenv BLAS None; setenv LAPACK None


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

> tar ztvf swig-1.3.27.tar.gz

As for numpy, you have two choices: system-wide installation (files
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

-------------------
4c.Installing pyRXP
-------------------

Get the latest pyRXP from http://www.reportlab.org/pyrxp.html, and
unpack it in a directory of your choice.

For system-wide installation, cd into the unpacking directory and issue
the command

> python setup.py install

For a local installation, cd into the unpacking directory and issue the
command

> python setup.py install --prefix=[userdir]

-----------------
4d.Installing PyX
-----------------

PyX is a very nice alternative to Gnuplot, better suited to creating
high-quality PDF and EPS images than to interactive viewing. PyX is a
Python module which requires requires Python 2.1 or newer, and a
working TeX installation that includes Type1 fonts.  Try to start
'python', 'tex' and 'kpsewhich cmr10.pfb' (the later should issue a
full path of the requested Type1 font) to see if you have
everything. You can get PyX at the URL

http://sourceforge.net/project/showfiles.php?group_id=45430

For system-wide installation, untar the PyX archive, cd into the
distribution directory, and issue the command

python setup.py install

For a local installation, cd into the PyX distribution directory and
issue the command

> python setup.py install --prefix=. --root=[userdir]

where you should replace [userdir] with the directory used to install
numeric. The setup process will create the directory "etc" (and "lib",
if needed) within [userdir]. Then the synthlisa-setdir.csh and
synthlisa-setdir.sh scripts will automatically setup the Python path
correctly to access PyX.

4e.Installing matplotlib
------------------------

Use matplotlib to get MATLAB-style plotting in Python. But I warn you,
the installation is a bit involved. Get it from
http://matplotlib.sourceforge.net and then keep into consideration all
of the following:

- matplotlib will happily use numpy, which is now standard with synthLISA;
  the Python path to include numpy needs to be set before running
  matplotlib's setup.py, probably using bin/synthlisa-setdir.[csh|sh];

- matplotlib needs freetype, zlib, and libpng. If these libraries are
  not already installed on the system, they must be downloaded from

  http://freetype.sourceforge.net
  http://www.zlib.net
  http://libpng.sourceforge.net

  and installed. Further, you need to make sure that matplotlib can
  reach them by adding the installation basepath to the variable
  basedir (for the right platform) in matplotlib's setup.py file.
  For instance, if the packages were installed in /opt/lib
  on linux, the 'linux' line in the basedir definition must be
  changed to
  
  'linux'  : ['/usr/local', '/usr', '/opt'],

- it seems that OS X 10.4 already has freetype and zlib, which are included
  with the X server, or at least with the X development kit. However, it's
  necessary to tell matplotlib that they can be found in the base directory
  '/usr/X11R6'; on the contrary, libpng must be downloaded and installed;

- the best way to run matplotlib-calling code is from ipython -pylab; see

  http://ipython.scipy.org/moin/

- [NOTE: this is fixed in matplotlib 0.87.6] on OS X, there are some
  problems with system fonts, resulting in broken EPS output; the
  workaround is to put the matplotlib-provided "Vera" fonts first in the
  'font.sans-serif', 'font.serif', 'font.monospace' fields in matplotlibrc
  (copy it from the matplotlib distribution directory into ~/.matplotlib);
  at runtime, you can still do something like

  rcParams['font.sans-serif'] = ['Bitstream Vera Sans']
  rcParams['font.serif'] = ['Bitstream Vera Serif']
  rcParams['font.monospace'] = ['Bitstream Vera Sans Mono']

- otherwise, you can also use the matplotlib LaTeX output by setting
  text.usetex to True in matplotlibrc (or at runtime in rcParams); this
  will be slower (and it will require a working installation of LaTeX
  including dvipng), but probably better-looking; some older versions of
  libpng require text.dvipnghack in matplotlibrc to be set to True to
  avoid jagging in the Tk display (but the postscript is fine, anyway).

- update: the most reliable way to install matplotlib is from the
  Python "egg" available from http://matplotlib.sourceforge.net. See
  http://peak.telecommunity.com/DevCenter/EasyInstall for a discussion 
  of easy_install and Python eggs. Presumably the dependencies given above
  are linked statically, but I'm not sure. 

----------------------------
4f.Installing Synthetic LISA
----------------------------

If you're reading this file, you have already unpacked the Synthetic
LISA distribution. Again, you have the choice of installing
system-wide or locally. In addition, if you have installed numpy
locally, you'll need to tell the Synthetic LISA installer where it can
be found.

So we have four cases:

- System-wide Synthetic LISA with system-wide numpy: run

> python setup.py install

- System-wide Synthetic LISA with local numpy installed in
  [numericdir]: run

> python setup.py install --with-numeric=[numericdir]

- Local Synthetic LISA with system-wide numpy: run

> python setup.py install --prefix=[synthlisadir]

- Local Synthetic LISA with local numpy installed in [numericdir}: run

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
4g.Setting the Python path
--------------------------

To use numpy and Synthetic LISA, the Python interpreter needs to
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
and want to work with numpy or Synthetic LISA. If you want this job
to be done automatically, you can call the correct synthlisa-setdir
script to your .cshrc, .profile, or .login file.

---------------------------------
4h.Recompiling and Synthetic LISA
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

Note: if you have installed locally with default-install.py, you can
use bin/recompile-here.sh to recompile and reinstall synthLISA.

========================
5. System-specific notes
========================

-------------------------------------
5a.Mac OS X
-------------------------------------

First of all, at this time LISA is known to work correctly with
Panther (OS X 10.3) and Tiger (OS X 10.4).

Panther and Tiger have recent enough versions of Python. To make your
life easier, install the enhancements suggested at the URL

http://www.vallis.org/blogspace/osxtricks/macpython.html

Most important, enabling readline is vital for interactive Python
sessions. Try the following:

python `python -c "import pimp; print pimp.__file__"` -i readline

As for the other tools, gcc is included in the Apple Development Tools
(on your OS X CD, or already installed); you can get gnuplot from
fink.sourceforge.net, or from

http://www.lee-phillips.org/info/Macintosh/gnuplot.html

and probably in other places, as well; if X11 is not installed, you can
get it from

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

Other recent versions of Linux should also have no problems, as long as
the Python version >= 2.3. (For instance, RedHat 9.0 has Python 2.2.2,
which needs to be upgraded to work with Synthetic LISA.)

------------------------------------------
4c. Compilation notes for Windows (cygwin)
------------------------------------------

The following notes apply to the installation under cygwin of lisatools,
which includes synthlisa. The patches mentioned below are already included
in the lisatools distribution.

Cygwin platform: on Windows XP, using Cygwin 1.5.24-2, I installed its
Python (2.4.3-1), its Subversion (1.4.2-1), openssh (4.5p1-1), SWIG
(1.3.29-2), its gcc (3.4.4-3, but probably not needed since mingw is
standard), patchutils (0.2.31-1) and make (3.81-1).

To compile numpy, it is necessary to disable
numarray, by commenting the line

# config.add_subpackage('numarray')

It is also necessary to patch lisaXML/io-C/ezxml.c (which is also
replicated inside the LISA Simulator):

#define EZXML_NOMMAP 

I wouldn't hope to get a stellar performance from your processor under
Windows/Cygwin. But it works...

Also, do yourself a favor: escape the strictures of the Windows shell
by installing rxvt (through Cygwin's setup.exe), and configuring it as
per

http://c2.com/cgi/wiki?BetterCygwinTerminal

Under rxvt, it's even possible to run Cygwin's emacs (but first do "export
TERM=rxvt").

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
> from numpy import *

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

Now "noiseX" is a numpy array of length 16384, which can be written
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
the averaged spectrum]. The result is a 2D numpy array, where each
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

To this purpose, a SourceForge webpage is being set up with a user
forum, bug report forms, and other tools for interaction. Be sure to
use them!

==================
9. Asking for help
==================

Please refer to www.vallis.org/syntheticlisa to report bugs, request
features, ask for explanations, and generally interact with other
Synthetic LISA designers and with other users.

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
