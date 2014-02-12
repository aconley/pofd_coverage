//positionGenerator.h

#ifndef __positionGenerator__
#define __positionGenerator__

#include<string>
#include<vector>
#include<utility>

#include<fftw3.h>

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_interp.h>  

#include<hdf5.h>

#include "../include/ran.h"

/*!
  \brief Holds power spectrum

  Uses log spline interpolation internally
*/
class powerSpectrum {
 private:
  unsigned int nk; //!< Number of k values
  double mink; //!< Minimum value of k
  double maxk; //!< Maximum value of k
  double *logk; //!< log k
  double *logpk; //!< log P(k)

  gsl_interp_accel *acc; //!< Spline lookup accelerator
  gsl_spline *spline_logpk; //!< Spline of P(k)

  void init(const std::vector<double>&, const std::vector<double>&);

 public:
  powerSpectrum(const std::string&); //!< Constructor from file
  powerSpectrum(const std::vector<double>& k, 
		const std::vector<double>& p_k); //!< Constructor from vecs
  ~powerSpectrum();

  double getPk(double) const; //!< Get P_k for specified k
  double getLogPk(double) const; //!< Get Log P_k for specified k

  void writeToHDF5Handle(hid_t objid) const;
};

// In principle, we could also have a uniform position generator,
// maybe with a base class, etc.  But the uniform one is so trivial
// it's just silly.

/*!
  \brief Generates positions on the sky obeying a power spectrum
*/
class positionGeneratorClustered {
 private:
  unsigned int nx; //!< x dimension to generate over
  unsigned int ny; //!< y dimension to generate over
  unsigned int nyhalf; //!< ny / 2 + 1
  double pixsize; //!< Pixel size in arcsec (square pixels assumed)

  // Power spectrum
  powerSpectrum powspec;

  //Internal storage -- only allocated when needed
  double *scl; //!< Power spectrum scaling array (row major 2D nx by ny)
  double *probarr; //!< Normalized probability array (row major 2D nx by ny)
  fftw_complex* probarr_trans; //!< Fourier transformed prob array

  /*! \brief Use bilinear interpolation on prob image*/
  double interpolate(double, double) const;

  // FFTW stuff
  fftw_plan plan;     //!< Holds forward transformation plan
  fftw_plan plan_inv; //!< Holds inverse transformation plan

 public:
  positionGeneratorClustered(unsigned int, unsigned int, double,
			     const std::string&);
  ~positionGeneratorClustered();

  void generate(ran&); //!< generate from power spectrum
 
  /*!\brief Get position of single source, no interpolation*/
  std::pair<unsigned int, unsigned int> getPosition(ran&) const;

  /*!\brief Get position of single source, with oversampling and interpolation*/
  std::pair<unsigned int, unsigned int> getPosition(ran&, unsigned int) const;

  int writeProbToFits(const std::string&) const; //!< Write prob image to FITS file

  void writeToHDF5Handle(hid_t objid) const;

};

#endif
