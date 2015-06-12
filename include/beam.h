//beam.h

//1-band Gaussian PSF 

#ifndef __beam__
#define __beam__

#include "../include/global_settings.h"
#include "../include/ran.h"
#include "../include/fourierFilter.h"

/*!
  \brief Represents PSF parameters for single Gaussian beam

  Note that the beam is stored discretely.

  \ingroup Beams
*/
class beam {
 private :
  double fwhm; //!< FWHM of beam, in arcsec
  double rhosq; //!< Convenience variable \f$\rho = 4 \log\left(2\right) 3600^2/FWHM^2)\f$
  double noise; //!< Fractional beam noise

  mutable ran *rangen; //!< Random number generator, if needed for fills

  /*! \brief Inner function to add noise */
  void addNoise(unsigned int, double, double* const) const;

  /*!\brief Inner beam generator, no filtering, arbitrary extent*/
  void getRawBeam(unsigned int n1, unsigned int n2, double pixsize, 
                  double* const bm) const;

  /*!\brief Inner beam generator, no filtering, with oversampling, arb extent*/
  void getRawBeam(unsigned int n1, unsigned int n2, double pixsize, 
                  unsigned int oversamp, double* const bm) const;

 public :
  beam(double FWHM=10.0) noexcept; //!< Constructor with FWHM
  ~beam(); //!< Desctructor

  void setFWHM(double) noexcept; //!< Set the FWHM

  double getFWHM() const noexcept { return fwhm; } //!< Get the FWHM
  double getRhoSq() const noexcept { return rhosq; } //!< Get \f$\rho^2\f$

  void setNoise(double NOISE) noexcept { noise = NOISE; } //!< Set the beam noise
  double getNoise() const noexcept { return noise; } //!< Get beam noise

  double getEffectiveArea() const noexcept; //<! Get effective area of beam in sq deg
  double getEffectiveAreaSq() const noexcept; //!< Get effective area of beam^2 in sq deg

  /*!\brief Get factorized beam, no filtering*/
  void getBeamFac(unsigned int n, double pixsize, double* const fac) const; 

  /*!\brief Get 2D beam, square*/
  void getBeam(unsigned int n, double pixsize, double* const, 
               const fourierFilter* const=nullptr) const;

  /*!\brief Get 2D beam, arbitrary size*/
  void getBeam(unsigned int n1, unsigned int n2, double pixsize, double* const, 
               const fourierFilter* const=nullptr) const;
  
  /*!\brief Get 2D beam with oversampling, square*/
  void getBeam(unsigned int n, double pixsize, unsigned int oversamp,
               double* const, const fourierFilter* const=nullptr) const;

  /*!\brief Get 2D beam with oversampling, arbitrary size*/
  void getBeam(unsigned int n1, unsigned int n2, double pixsize, 
               unsigned int oversamp, double* const, 
               const fourierFilter* const=nullptr) const;

  /*!\brief Write the beam to a FITS file*/
  void writeToFits(const std::string& outfile, double pixsize, 
                   double nfwhm=3.5, unsigned int oversamp=1,
                   const fourierFilter* const=nullptr, bool inverse=false) const;
};

/*!
  \brief Holds histogrammed beam

  Note that the absolute value of the negative beam is held.

  \ingroup beams
*/

class beamHist {
 private:
  bool has_data; //!< Have we been filled?

  bool inverse; //!< Are we holding the inverse beam instead of the beam?
  unsigned nbins; //!< Number of bins

  double fwhm; //!< FWHM of beam we are storing
  double nfwhm; //!< Number of FWHM out we go
  double nfwhmkeep; //!< Number of FWHM we actually keep (after filtering)
  double pixsize; //!< Pixsize of beam sampling
  double eff_area; //!< Effective area of beam in deg^2
  unsigned int oversamp; //!< Oversampling factor

  bool isHipass; //!< Was hipass filtering applied during fill?
  double filtscale; //!< High-pass filtering scale, in arcsec
  double qfactor; //!< High-pass filtering apodization

  bool isMatched; //!< Was matched filtering applied during the fill?
  double matched_fwhm; //!< FWHM of matched filter (doesn't have to match this beam)
  double matched_sigi; //!< Instrument sigma of matched filter
  double matched_sigc; //!< Confusion sigma of matched filter

  unsigned int n_pos; //!< Number of positive beam histogram elements filled
  unsigned int* wt_pos; //!< Positive weights
  double *bm_pos; //!< Beam elements in each bin.

  unsigned int n_neg; //!< Number of negative beam histogram elements filled
  unsigned int* wt_neg; //!< Negative weights
  double *bm_neg; //!< Beam elements in each bin.

  dblpair minmax_pos; //!< Min/max value of pos beam (not inverse beam!)
  dblpair minmax_neg; //!< Min/max value of |neg| beam (not inverse beam!)

 public:
  beamHist(unsigned int NBINS); //!< Constructor
  ~beamHist(); //!< Destructor
  
  bool hasData() const { return has_data; }
  bool isInverse() const { return inverse; }
  unsigned int getNbins() const { return nbins; }
  double getFWHM() const { return fwhm; }
  double getNFWHM() const { return nfwhm; }
  double getNFWHMKeep() const { return nfwhmkeep; }
  double getPixsize() const { return pixsize; }
  unsigned int getOversamp() const { return oversamp; }
  double getEffectiveArea() const { return eff_area; }
  bool hasPos() const { return n_pos > 0; }
  bool hasNeg() const { return n_neg > 0; }
  unsigned int getNPos() const { return n_pos; }
  unsigned int getNNeg() const { return n_neg; }

  // This is very bad, but also very, very useful for speed purposes
  const unsigned int* getWtPos() const { return wt_pos; }
  const double* getBmPos() const { return bm_pos; }
  const unsigned int* getWtNeg() const { return wt_neg; }
  const double* getBmNeg() const { return bm_neg; }
  
  bool isHipassFiltered() const { return isHipass; }
  double getFiltScale() const { return filtscale; } //!< Get filtering scale
  double getFiltQFactor() const { return qfactor; } //!< Get filtering scale
  bool isMatchFiltered() const { return isMatched; }
  double getFiltFWHM() const { return matched_fwhm; }
  double getFiltSigInst() const { return matched_sigi; }
  double getFiltSigConf() const { return matched_sigc; }

  // Min/max values of real beam (not inverse)
  dblpair getMinMaxPos() const {return minmax_pos;} //!< Get min/max pos beam
  dblpair getMinMaxNeg() const {return minmax_neg;} //!< Get min/max neg beam

  /*!\brief Fill from beam, symmetric number of fwhm version*/
  void fill(const beam& bm, double nfwhm, double pixsize,
            bool inv=false, unsigned int oversamp=1,
            const fourierFilter* const filt=nullptr, double num_fwhm_keep=0.0);

  /*!\brief Fill from beam, arbitrary size version */
  void fill(const beam& bm, unsigned int n1, unsigned int n2, double pixsize,
            bool inv=false, unsigned int oversamp=1,
            const fourierFilter* const filt=nullptr, double num_fwhm_keep=0.0);

  /*! \brief Write out as FITS file*/
  void writeToFits(const std::string&) const;

};

#endif
