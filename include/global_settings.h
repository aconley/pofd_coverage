#ifndef __global_settings__
#define __global_settings__

#include<vector>
#include<string>
#include<cmath>
#include<fftw3.h>

// FFT precision type
#if USEFFT32BIT

#define FFTFLOAT float
#define FITSIOIMG FLOAT_IMG
#define TFITSFLOAT TFLOAT
#define FFTWMALLOC fftwf_malloc
#define FFTWFREE fftwf_free
#define FFTWCOMPLEX fftwf_complex
#define FFTWPLAN fftwf_plan
#define FFTWDFTR2C1D fftwf_plan_dft_r2c_1d
#define FFTWDFTC2R1D fftwf_plan_dft_c2r_1d
#define FFTWDFTR2C2D fftwf_plan_dft_r2c_2d
#define FFTWDFTC2R2D fftwf_plan_dft_c2r_2d
#define FFTWDESTROYPLAN fftwf_destroy_plan
#define FFTWIMPORTWIS fftwf_import_wisdom_from_file
#define FFTWEXECUTE fftwf_execute
#define FFTWEXECUTEDFTR2C fftwf_execute_dft_r2c

#else

#define FFTFLOAT double
#define FITSIOIMG DOUBLE_IMG
#define TFITSFLOAT TDOUBLE
#define FFTWMALLOC fftw_malloc
#define FFTWFREE fftw_free
#define FFTWCOMPLEX fftw_complex
#define FFTWPLAN fftw_plan
#define FFTWDFTR2C1D fftw_plan_dft_r2c_1d
#define FFTWDFTC2R1D fftw_plan_dft_c2r_1d
#define FFTWDFTR2C2D fftw_plan_dft_r2c_2d
#define FFTWDFTC2R2D fftw_plan_dft_c2r_2d
#define FFTWDESTROYPLAN fftw_destroy_plan
#define FFTWIMPORTWIS fftw_import_wisdom_from_file
#define FFTWEXECUTE fftw_execute
#define FFTWEXECUTEDFTR2C fftw_execute_dft_r2c

#endif

/*!
  \brief Global convenience variables
*/
namespace pofd_coverage {
  const char version[] = "0.2.3"; //Version number

  const double n_sigma_shift = 8.0; //!< Shift amount
  const double n_sigma_pad = 10.0; //!< Noise padding size in sigma
  const double n_sigma_shift2d = 4.0; //!< Shift amount
  const double n_sigma_pad2d = 6.0; //!< Noise padding size in sigma, 2D
  const double pi = 3.141592653589793238462643383279502884197; //!< \f$\pi\f$
  const double two_pi = 2.0*pi; //!< \f$2 \pi \f$
  const double sqrt_pi_over_two = sqrt(pi/2.0); //!< \f$\sqrt{\pi/2}\f$
  const double isqrt_two_pi = 1.0/sqrt(two_pi); 
  const double logfac = std::log(10.0); //!< Conversion base e, base 10
  const double ilogfac = 1.0/logfac; //!< Inverse conversion factor
  const FFTFLOAT smallval=exp2(-100); //!< Log of small number (base 2)
  const FFTFLOAT smalllogval=-100; //!< Log of small number
  const double log2toe = log(2); //!< Multiply by this to go from log2 to ln
  const double rhofac = 4*std::log(2)*3600.0*3600.0;
  const double fwhmfac = 1.0/sqrt(8*std::log(2)); //!< Conv FWHM to sigma

  /*! \brief Maximum transform size.  Make sure it's a power of 2*/
  const unsigned int nmax = 16777216;

  //Powers of 2
  const int npow2 = 24; //!< Number of powers of 2 in pow2
  const int pow2[npow2+1] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 
			     2048, 4096, 8192, 16384, 32768, 65536, 131072, 
			     262144, 524288, 1048576, 2097152, 4194304,
			     8388608, 16777216}; //!< Powers of 2
}
#endif
