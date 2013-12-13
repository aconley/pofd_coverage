//beam.h

//1-band Gaussian PSF 

#ifndef __beam__
#define __beam__

#include<utility>

#include "../include/hipassFilter.h"

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
  void getBeam(unsigned int n, double pixsize, double* const, 
	       hipassFilter* const=NULL) const;

  /*!\brief Get 2D beam with oversampling*/
  void getBeam(unsigned int n, double pixsize, unsigned int oversamp,
	       double* const, hipassFilter* const=NULL) const;
};

/*!
  \brief Holds histogrammed beam

  \ingroup beams
*/

class beamHist {
 private:
  bool has_data; //!< Have we been filled?

  bool inverse; //!< Are we holding the inverse beam instead of the beam?
  unsigned nbins; //!< Number of bins

  double fwhm; //!< FWHM of beam we are storing
  double pixsize; //!< Pixsize of beam sampling
  double eff_area; //!< Effective area of beam in deg^2

  bool n_pos; //!< Number of positive beam histogram elements filled
  unsigned int* wt_pos; //!< Positive weights
  double *bm_pos; //!< Beam elements in each bin.

  bool n_neg; //!< Number of negative beam histogram elements filled
  unsigned int* wt_neg; //!< Negative weights
  double *bm_neg; //!< Beam elements in each bin.
 public:
  beamHist(unsigned int NBINS); //!< Constructor
  ~beamHist(); //!< Destructor
  
  bool hasData() const { return has_data; }
  bool isInverse() const { return inverse; }
  unsigned int getNbins() const { return nbins; }
  double getFWHM() const { return fwhm; }
  double getPixsize() const { return pixsize; }
  double getEffectiveArea() const { return eff_area; }
  unsigned int getNPos() const { return n_pos; }
  unsigned int getNNeg() const { return n_neg; }

  // This is very bad, but also very, very useful for speed purposes
  const unsigned int* getWtPos() const { return wt_pos; }
  const double* getBmPos() const { return bm_pos; }
  const unsigned int* getWtNeg() const { return wt_neg; }
  const double* getBmNeg() const { return bm_neg; }
  
  // Min/max values
  std::pair<double, double> getMinMaxPos() const;
  std::pair<double, double> getMinMaxNeg() const;

  // Fill from beam
  void fill(const beam& bm, double nfwhm, double pixsize,
	    hipassFilter* const filt=NULL, double filtscale=360.0,
	    bool inv=false, unsigned int oversamp=1);
};

#endif
