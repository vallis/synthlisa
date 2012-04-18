# Synthetic LISA #

This is Michele Vallisneri's _Synthetic LISA_, now on GitHub with a cleaned-up version 2.

I thought I would celebrate the 2011 demise of the U.S. LISA project and the 2012 non-selection of the European eLISA by sharing this code as widely as possibly with posterity. 

## What is this? ##

_Synthetic LISA_, developed by Michele Vallisneri and John Armstrong at the Jet Propulsion Laboratory under the auspices of the LISA Mission Science Office, is a C++/Python package to simulate the LISA science process at the level of scientific and technical requirements. Synthetic LISA:

* generates synthetic time series of the LISA fundamental noises, as filtered through all the TDI observables;
* provides a streamlined module to compute the TDI responses to gravitational waves, according to a full model of TDI, including the motion of the LISA array, and the temporal and directional dependence of the armlengths;
* was a central component of [lisatools](http://lisatools.googlecode.com), the software pipeline used to generate the [Mock LISA Data Challenge](http://astrogravs.nasa.gov/docs/mldc) datasets;
* can be used for a variety of LISA-like mission designs.

See the homepage at [www.vallis.org](http://www.vallis.org/syntheticlisa).

## Requirements ##

* A a working installation of Python (2.6.X or 2.7.X) and of a Python-interoperable C/C++ compiler, preferably gcc (>= 4.0).
* [numpy](http://numpy.scipy.org). If you have [pip](http://www.pip-installer.org), `pip install numpy` will do.
* [pyRXP](http://www.reportlab.com/software/opensource/pyrxp). You can do `pip install http://svn.reportlab.com/svn/public/pyRXP/trunk/src`.
* [SWIG](http://www.swig.org). But you'll need it only if you need to modify the Python API.

## Installation ##

If you have all the requirements,

    python setup.py install

will do, which will install to the default location for your Python setup. To install to a specific directory, use

    python setup.py install --prefix=$INSTALLDIR
    
where `$INSTALLDIR` could be `$HOME`, `/usr/local`, `$VIRTUAL_ENV` (if you use [virtualenv](http://www.virtualenv.org)), etc.

## Usage ##

A general description of the formulation, implementation, and usage of
_Synthetic LISA_ can be found in [`doc/synthlisa.pdf`](https://github.com/vallis/synthlisa/blob/master/doc/synthlisa.pdf) and
[`doc/manual.pdf`](https://github.com/vallis/synthlisa/blob/master/doc/manual.pdf) included with the package. Those documents however are somewhat outdated. See also [`doc/history.txt`](https://github.com/vallis/synthlisa/blob/master/doc/history.txt) for a list of changes.

The briefest summary is that _Synthetic LISA_ implements (at the C++ level) a number of objects that describe the LISA orbits, noises, and TDI observables, as well as gravitational-wave sources. The user writes Python scripts to create and connect these objects, and to generate synthetic data form them. The Python API is documented reasonably well in the docstrings found in [`lisasim/lisasim-swig.i`](https://github.com/vallis/synthlisa/blob/master/lisasim/lisasim-swig.i), which are accessible with Python's `help`. The example scripts in [`examples`](https://github.com/vallis/synthlisa/tree/master/examples) (and especially [`examples/manual-examples`](https://github.com/vallis/synthlisa/tree/master/examples/manual-examples)) are a good place to start.

## Credit ##

If you use _Synthetic LISA_ in your work, please cite M. Vallisneri, "Synthetic LISA: Simulating Time Delay Interferometry in a Model LISA," Phys. Rev. D 71, 022001 (2005), which you can find in [PRD](http://dx.doi.org/10.1103/PhysRevD.71.022001) and as [gr-qc/0407102](http://www.arxiv.org/abs/gr-qc/0407102).

## License ##

_Synthetic LISA_ is licensed under the Caltech [public domain license](https://github.com/vallis/synthlisa/blob/master/doc/license.txt).
