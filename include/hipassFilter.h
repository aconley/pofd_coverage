//hipassFilter.h

#ifndef __hipassFilter__
#define __hipassFilter__

#include<fftw3.h>

/*!
  \brief A high-pass filter with an apodized gaussian edge
*/
class hipassFilter {
 private:
  double filtscale; //!< Filter cutoff radius in arcseconds
  double qfactor; //!< Gaussian sigma in terms of filtscale

  unsigned int nx; //!< Current maximum capacity data array
  unsigned int ny; //!< Current maximum capacity data array
  unsigned int nyhalf; //!< ny / 2 + 1
  fftw_complex* transdata; //!< Transformed data array
  
  unsigned int nxplan; //!< Current plan size, x
  unsigned int nyplan; //!< Current plan size, y
  fftw_plan plan; //!< Forward transform plan from input data to transdata
  fftw_plan plan_inv; //!< Backwards transform from transdata back to input data

  double meanSub(unsigned int, double* const) const; //!< Subtracts mean from input

 public:
  hipassFilter(double f, double q=0.1); //!< Constructor
  ~hipassFilter();
  
  double getFiltScale() const { return filtscale; }
  double getQFactor() const { return qfactor; }
  void setFiltScale(double val) { filtscale = val; }
  void setQFactor(double val) { qfactor = val; }

  /*! \brief Apply filtering*/
  void filter(double pixscale, unsigned int NX, unsigned int NY, 
	      double* const data);
};

#endif
