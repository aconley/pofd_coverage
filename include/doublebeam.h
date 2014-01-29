//beam.h

//2-band Gaussian PSF 

#ifndef __doublebeam__
#define __doublebeam__

#include "../include/global_settings.h"
#include "../include/hipassFilter.h"

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
  
  /*!\brief Inner beam generator, no filtering*/
  void getRawBeam(unsigned int band, unsigned int n, double pixsize, 
		  double* const bm) const;

  /*!\brief Inner beam generator, no filtering, with oversampling*/
  void getRawBeam(unsigned int band, unsigned int n, double pixsize, 
		  unsigned int oversamp, double* const bm) const;

 public :
  doublebeam(double FWHM1=10.0, double FWHM2=15.0); //!< Constructor with FWHM

  void setFWHM(double, double); //!< Set the FWHM values

  dblpair getFWHM() const { return std::make_pair(fwhm1, fwhm2); }
  dblpair getRhoSq() const {return std::make_pair(rhosq1, rhosq2);}

  /*! \brief Get effective areas of beam in sq deg*/
  dblpair getEffectiveArea() const; 
  /*! \brief Get effective areas of beam^2 in sq deg*/
  dblpair getEffectiveAreaSq() const;

  /*! \brief Get factorized beam*/
  void getBeamFac(unsigned int band, unsigned int n, 
		  double pixsize, double* const fac) const;

  /*!\brief Get 2D beam*/
  void getBeam(unsigned int band, unsigned int n, double pixsize, 
	       double* const, hipassFilter* const=NULL) const;
  /*!\brief Get 2D beam with oversampling*/
  void getBeam(unsigned int band, unsigned int n, double pixsize, 
	       unsigned int oversamp, double* const, 
	       hipassFilter* const=NULL) const;

  /*!\brief Write the beams to a FITS file*/
  void writeToFits(const std::string& outfile, double pixsize, 
		   double nfwhm=3.5, unsigned int oversamp=1,
		   hipassFilter* const=NULL, bool inverse=false) const;
};

/*!
  \brief Holds histogrammed double beam

  Note that the absolute value of the negative beam is held.

  There are 4 sign components to keep track of here, based on the
  relative signs of the two beams.  We keep them in the order
  pp, pn, np, nn.  The histograms are also two dimensional, but
  of course we don't keep empty bins.

  \ingroup beams
*/

class doublebeamHist {
 private:
  bool has_data; //!< Have we been filled?

  bool inverse; //!< Are we holding the inverse beam instead of the beam?
  unsigned nbins; //!< Number of bins

  double fwhm1; //!< FWHM of beam we are storing, band 1
  double fwhm2; //!< FWHM of beam we are storing, band 2
  double nfwhm; //!< Number of FWHM out we go (using the larger of fwhm1/2)
  double pixsize; //!< Pixsize of beam sampling
  double eff_area1; //!< Effective area of beam in deg^2, band 1
  double eff_area2; //!< Effective area of beam in deg^2, band 2
  unsigned int oversamp; //!< Oversampling factor

  bool keep_filt; //!< Keep filt allocated, or dealloc between calls
  double filtscale; //!< Filtering scale, in arcsec
  hipassFilter* filt; //!< Hi pass filter, if any is being applied

  unsigned int n[4]; //!< Number of histogram elements filled in pp, pn, np, nn
  unsigned int* wt[4]; //!< Weights
  double *bm1[4]; //!< Beam elements in each bin, band 1
  double *bm2[4]; //!< Beam elements in each bin, band 2

  dblpair minmax1[4]; //!< Min/max beam (not inverse beam!) in each sign component, band 1
  dblpair minmax2[4]; //!< Min/max beam (not inverse beam!) in each sign component, band 2
 public:

  doublebeamHist(unsigned int NBINS, double FILTSCALE=0.0,
		 bool KEEP_FILT_INMEM=true); //!< Constructor
  ~doublebeamHist(); //!< Destructor
  
  bool hasData() const { return has_data; }
  bool isInverse() const { return inverse; }
  unsigned int getNbins() const { return nbins; }
  dblpair getFWHM() const { return std::make_pair(fwhm1, fwhm2); }
  double getNFWHM() const { return nfwhm; }
  double getPixsize() const { return pixsize; }
  unsigned int getOversamp() const { return oversamp; }
  dblpair getEffectiveArea() const { return std::make_pair(eff_area1, eff_area2); }
  bool hasSign(unsigned int idx) const { return n[idx] > 0; }
  unsigned int getN(unsigned int idx) const { return n[idx]; }

  // This is very bad, but also very, very useful for speed purposes
  const unsigned int* getWt(unsigned int idx) const { return wt[idx]; }
  const double* getBm1(unsigned int idx) const { return bm1[idx]; }
  const double* getBm2(unsigned int idx) const { return bm2[idx]; }
  
  bool isFiltered() const { return filtscale>0; } //!< Is the beam filtered
  double getFiltScale() const { return filtscale; } //!< Get filtering scale

  // Min/max values
  /*!\brief Get min/max band 1 (non-inverse beam)*/
  dblpair getMinMax1(unsigned int i) const {return minmax1[i];} 
  /*!\brief Get min/max band 2 (non-inverse beam)*/
  dblpair getMinMax2(unsigned int i) const {return minmax2[i];}

  /*!\brief Fill from beam*/
  void fill(const doublebeam& bm, double nfwhm, double pixsize,
	    bool inv=false, unsigned int oversamp=1);

  /*! \brief Write out as FITS file*/
  void writeToFits(const std::string&) const;
};

#endif
