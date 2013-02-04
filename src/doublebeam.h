//beam.h

//2-band Gaussian PSF 

#ifndef __doublebeam__
#define __doublebeam__

/*!
  \brief Represents PSF parameters for two Gaussian beams.

  Note that the beams are stored discretely.

  \ingroup Beams
*/
class doublebeam {
 private :
  double fwhm1; //!< FWHM of beam, in arcsec, band 1
  double fwhm2; //!< FWHM of beam, in arcsec, band 2
  double rhosq1; //!< Convenience variable \f$\rho = 4 \log\left(2\right) 3600^2/FWHM^2)\f$, band 1
  double rhosq2; //!< Convenience variable \f$\rho = 4 \log\left(2\right) 3600^2/FWHM^2)\f$, band 2
  
  void getBeamFac(unsigned int, unsigned int, double, double* const) const;
  void getBeam(unsigned int, unsigned int, double, double* const) const;

 public :
  doublebeam(double FWHM1=10.0, double FWHM2=15.0); //!< Constructor with FWHM

  void setFWHM(double, double); //!< Set the FWHM values

  //Because of the way these are used, it is convenient to just define
  // two versions of each accessor rather than returning both always, or
  // using a selection argument to chose which one. 
  double getFWHM1() const {return fwhm1;} //!< Get the FWHM value, band 1
  double getRhoSq1() const {return rhosq1;} //!< Get the $\f\rho^2\f$ value, band 1
  double getFWHM2() const {return fwhm2;} //!< Get the FWHM value, band 2
  double getRhoSq2() const {return rhosq2;} //!< Get the $\f\rho^2\f$ value, band 2

  /*! \brief Get effective areas of beam in sq deg, band 1*/
  double getEffectiveArea1() const; 
  /*! \brief Get effective areas of beam^2 in sq deg, band 1*/
  double getEffectiveAreaSq1() const;
  /*! \brief Get effective areas of beam in sq deg, band 2*/
  double getEffectiveArea2() const; 
  /*! \brief Get effective areas of beam^2 in sq deg, band 2*/
  double getEffectiveAreaSq2() const;

  /*!\brief Get factorized beam, band 1*/
  void getBeamFac1(unsigned int, double, double* const) const; 
  /*!\brief Get factorized beam, band 2*/
  void getBeamFac2(unsigned int, double, double* const) const; 

  /*!\brief Get 2D beam, band 1*/
  void getBeam1(unsigned int, double, double* const) const;
  /*!\brief Get 2D beam, band 2*/
  void getBeam2(unsigned int, double, double* const) const;


  /*!\brief Get 2D histogrammed beams */
  void getBeamHist(unsigned int, double, unsigned int,
		    unsigned int&, unsigned int* const,
		    double* const, double* const, bool=false) const;

};

#endif
