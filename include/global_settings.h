#ifndef __global_settings__
#define __global_settings__

#include<vector>
#include<string>
#include<cmath>
#include<limits>
#include<utility>

/*!
\mainpage pofd_coverage

This is the source documentation for pofd_coverage, a package
which tests one and two-dimensional P(D) fitting against simulated
maps.  It is related to pofd_affine, but:
  - Can only fit for one variable, the number of sources per square degree.
    This makes it simpler to test the code.
  - Deals purely with simulated data, including methods for generating
    such data.
  - Uses a slightly simpler (brute force) implementation of the R computation.
  - The beam is strictly Gaussian, although it can be filtered.

There are three intended purposes for this code:
  - Because the one-parameter fit is fast and simple, it is possible to
    establish with much more certainty that this code works correctly.
    It can then be compared with pofd_affine for debugging purposes.
  - It can be used to test out modifications more easily (e.g. floating point
    FFTs, aliasing, the effects of clustering, matched filtering) on simulated
    data.
  - In addition to finding the best fit, it is capable of mapping out
    likelihood curves, and thus can be used to empirically determine the
    likelihood correction factors to match the observed scatter (that is --
    to get the statistical coverage right).

Because this code is so strongly related to pofd_affine, the documentation
is a bit sparser; some things may be better described there.  That's also
where most of the unit tests live.
*/


/*!
  \brief Global convenience variables
*/
namespace pofd_coverage {
  const char version[] = "0.3.2"; //Version number

  const double n_zero_pad = 7.5; //!< Zero padding size in sigma
  const double pi = 3.141592653589793238462643383279502884197; //!< \f$\pi\f$
  const double two_pi = 2.0*pi; //!< \f$2 \pi \f$
  const double sqrt_pi_over_two = sqrt(pi / 2.0); //!< \f$\sqrt{\pi/2}\f$
  const double isqrt_two_pi = 1.0 / sqrt(two_pi);
  const double logfac = log2(10.0); //!< Conversion base 2, base 10
  const double ilogfac = 1.0 / logfac; //!< Inverse conversion factor
  const double smallval= exp2(-100); //!< Log of small number (base 2)
  const double smalllogval=-100; //!< Log of small number
  const double log2toe = log(2); //!< Multiply by this to go from log2 to ln
  const double rhofac = 4 * std::log(2) * 3600.0 * 3600.0;
  const double fwhmfac = 1.0 / sqrt(8 * std::log(2)); //!< Conv FWHM to sigma
  const double qnan = std::numeric_limits<double>::quiet_NaN();

  /*! \brief Maximum transform size.  Make sure it's a power of 2*/
  const unsigned int nmax = 16777216;

  //Powers of 2
  const int npow2 = 24; //!< Number of powers of 2 in pow2
  const int pow2[npow2+1] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
			     2048, 4096, 8192, 16384, 32768, 65536, 131072,
			     262144, 524288, 1048576, 2097152, 4194304,
			     8388608, 16777216}; //!< Powers of 2
}

// Typedefs
typedef std::pair<double, double> dblpair;
typedef std::pair<bool, bool> blpair;

#endif
