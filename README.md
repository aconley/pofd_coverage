pofd_coverage
===========

This is a package to empirically determine the appropriate likelihood
normalization for a P(D) analysis for a broken power law model, or (in
2D) a broken power law coupled with a log normal color ratio
distribution.  It is a companion to
[pofd_affine](https://github.com/aconley314/pofd_affine).  The basic
idea is to implement (essentially) the same number counts model, but
only fit for one parameter -- the total number of sources -- so that
it is trivial to simulate input images and to map out the likelhood
curve.  This can be used to determine the right normalization in order
to obtain the correct statistical coverage for a P(D) analysis.
Unfortunately, the correct expression is not known from theory.

Additionally, a python library (analyze_pofd_coverage) is provided
to process the results of pofd_coverage runs, which should be
accessed through the script analyze_coverage.py.

### Installation

Installation is via the standard UNIX `configure` and
`make`. pofd_affine depends on the following packages:
* [FFTW](http://www.fftw.org/).  Generating a 
  [FFTW wisdom](http://www.fftw.org/fftw-wisdom.1.html)
  file will speed up the code signficantly.
* [cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/)
* The [GNU scientific library](http://www.gnu.org/software/gsl/)
* [HDF5](http://www.hdfgroup.org/HDF5/)

It may be necessary to tell configure where to look for these
libraries -- see `configure --help`.

The python library (analyze_pofd_coverage) depends on the python
packages depends on having a recent version of python (tested on
3.3 and above) as well as:
* [numpy](http://numpy.scipy.org/)
* [scipy](http://numpy.scipy.org/)
* [astropy](http://www.astropy.org/)
* [h5py](http://www.h5py.org/)

This library can be installed by going into the analyze_pofd_coverage
directory and executing the standard python setup script:

        python setup.py install

which provides (and installs) a analyze_coverage.py script, some
description of which can be accessed via

        analyze_coverage.py --help


### Documentation

There isn't much.  Doxygen style documentation can be
generated using

	make docs

but that isn't really all that helpful.   All of the command line
routines (e.g., pofd_coverage_runSim) come with some documentation
invoked by providing the --help option:

        pofd_coverage_runSim --help

### References
* The original P(D) paper is [Scheuer (1957)](http://dx.doi.org/10.1017/S0305004100032825),
* A more recent and easier to read discussion is
  [Patanchon et al. (2009)](http://dx.doi.org/10.1088/0004-637X/707/2/1750),
  which gives results from BLAST.
* The most recent P(D) results from SPIRE can be found at
  [Glenn et al. (2010)](http://dx.doi.org/10.1111/j.1365-2966.2010.17781.x).
