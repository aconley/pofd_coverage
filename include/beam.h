//beam.h

//1-band Gaussian PSF 

#ifndef __beam__
#define __beam__

/*!
  \brief Represents PSF parameters for single Gaussian beam

  Note that the beam is stored discretely.

  \ingroup Beams
*/
class beam {
 private :
  double fwhm; //!< FWHM of beam, in arcsec
  double rhosq; //!< Convenience variable \f$\rho = 4 \log\left(2\right) 3600^2/FWHM^2)\f$
  
  /*!\brief Inner beam generator, no filtering*/
  void getRawBeam(unsigned int n, double pixsize, double* const bm) const;

  /*!\brief Inner beam generator, no filtering, with oversampling*/
  void getRawBeam(unsigned int n, double pixsize, unsigned int oversamp,
		  double* const bm) const;

 public :
  beam(double FWHM=10.0); //!< Constructor with FWHM

  void setFWHM(double); //!< Set the FWHM

  double getFWHM() const { return fwhm; } //!< Get the FWHM
  double getRhoSq() const { return rhosq; } //!< Get \f$\rho^2\f$

  double getEffectiveArea() const; //<! Get effective area of beam in sq deg
  double getEffectiveAreaSq() const; //!< Get effective area of beam^2 in sq deg

  /*!\brief Get factorized beam, no filtering*/
  void getBeamFac(unsigned int, double, double* const) const; 

  /*!\brief Get 2D beam*/
  void getBeam(unsigned int n, double pixsize, double* const, double) const;

  /*!\brief Get 2D beam with oversampling*/
  void getBeam(unsigned int n, double pixsize, unsigned int oversamp,
	       double* const, double) const;

  /*!\brief Get histogrammed 2D beam*/
  void getBeamHist(unsigned int, double, unsigned int,
		   unsigned int&, unsigned int* const, double* const,
		   unsigned int&, unsigned int* const, double* const,
		   double, bool=false) const;		       

  /*!\brief Get histogrammed 2D beam with oversampling*/
  void getBeamHist(unsigned int, double, unsigned int, unsigned int,
		   unsigned int&, unsigned int* const, double* const,
		   unsigned int&, unsigned int* const, double* const,
		   double, bool=false) const;		       
};

#endif
