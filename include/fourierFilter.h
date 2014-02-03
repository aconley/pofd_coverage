//matchedFilter.h

#ifndef __matchedFilter__
#define __matchedFilter__

#include<fftw3.h>

#include "../include/global_settings.h"

/*!
  \brief A class to apply Fourier space based filters.

  Two types are supported (which can be combined):
    1) A High pass filter with an apodized Gaussian edge
    2) A confusion matched filter.

  The confusion matched filter is based on Chapin et al. 2011, MNRAS, 
  411, 505, Appendix A

  Construction of the actual filter is deferred until the first
  time it is applied.

  You set this up with the pixel scale and filter parameters of the data you
  plan to apply it to at construction.  It can't then be changed,
  although the physical extent of the image it is applied to can be.
*/
class fourierFilter {
 private:
  mutable bool initialized; //!< Has the filter been set up
  bool doHipass; //!< Do Highpass filtering
  bool doMatched; //!< Do matched filtering

  mutable unsigned int nx; //!< Size of filter region, x
  mutable unsigned int ny; //!< Size of filter region, y
  double pixscale; //!< Pixel scale, in arcseconds

  // Hipass filter stuff
  double filtscale; //!< Filter cutoff radius in arcsec
  double qfactor; //!< Gaussian sigma for apodization in terms of filtscale

  // Matched filter stuff
  double fwhm; //!< FWHM in arcsec
  double sig_inst; //!< Instrument noise in Jy
  double sig_conf; //!< Confusion noise in Jy

  // FFTW stuff
  unsigned int fftw_plan_style; //!< What type of FFTW planning to use
  mutable unsigned int nyhalf; //!< ny / 2 + 1
  mutable fftw_plan plan; //!< Forward transform plan
  mutable fftw_plan plan_inv; //!< Backwards transform
  mutable fftw_complex* filt_fft; //!< conjugated FFT of filter
  mutable fftw_complex* map_fft; //!< Working variable for FFTed map

  void setup_beam(double* const bm) const; //!< Helper for setup_matched

  void setup_plans(double* const rl, fftw_complex* const im) const;
  void setup_matched() const; //!< Sets up the matched filter
  void setup_hipass() const; //!< Sets up the hipass filter
  bool setup(unsigned int NX, unsigned int NY) const; //!< Sets up filters

  double meanSub(double* const data) const; //!< Subtracts mean from data

 public:

  /*! \brief Matched filtering only constructor */
  explicit fourierFilter(double pixsize, double fwhm, double sigi, double sigc, 
			 bool quickfft=false);
  /*! \brief Hipass filtering only constructor */
  explicit fourierFilter(double pixsize, double fscale, double q=0.1, 
			 bool quickfft=false);
  /*! \brief Both types of filtering constructor */
  explicit fourierFilter(double pixsize, double fwhm, double sigi, 
			 double sigc, double fscale, double q=0.1, 
			 bool quickfft=false);
  ~fourierFilter();
  
  bool isMatched() const { return doMatched; }
  bool isHipass() const { return doHipass; }
  bool isBoth() const { return doMatched && doHipass; }
  double getFiltScale() const { return filtscale; }
  double getQFactor() const { return qfactor; }
  double getFWHM() const { return fwhm; }
  double getSigInst() const { return sig_inst; }
  double getSigConf() const { return sig_conf; }

  /*! \brief Apply filtering*/
  void filter(unsigned int n1, unsigned int n2, double pixsize,
	      double* const data) const;
};

#endif
