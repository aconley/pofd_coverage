//simImage.h

#ifndef __simImage__
#define __simImage__

#include<string>
#include<hdf5.h>

#include "../include/ran.h"
#include "../include/numberCounts.h"
#include "../include/positionGenerator.h"
#include "../include/fourierFilter.h"

/*!
  \brief Data class for creating and holding simulated
  image with Gaussian beams
 */

class simImage {
 private:
  unsigned int n1; //!< Number of pixels in dimension 1
  unsigned int n2; //!< Number of pixels in dimension 2
  unsigned int oversample; //!< Oversampling in generation
  double pixsize; //!< Size of pixels in arcsec.  Assumed square
  unsigned int ngen1; //!< Number of generated pixels in dimension 1
  unsigned int ngen2; //!< Number of generated pixels in dimension 2
  double pixsize_gen; //!< Size of pixels in internal mode

  double fwhm; //!< Beam FWHM in arcsec
  double sigi; //!< Raw instrumental noise, before any smoothing/filtering
  double esmooth; //!< Additional Gaussian smoothing FWHM in arcsec

  bool isHipass; //!< Was hipass filtering applied during fill?
  double filtscale; //!< High-pass filtering scale, in arcsec
  double qfactor; //!< High-pass filtering apodization

  bool isMatched; //!< Was matched filtering applied during the fill?
  double matched_fwhm; //!< FWHM of matched filter (doesn't have to match this beam)
  double matched_sigi; //!< Instrument sigma of matched filter
  double matched_sigc; //!< Confusion sigma of matched filter

  // sigi_final is only computed when needed, since it isn't
  // totally free to do if there is filtering
  mutable bool sigi_final_computed; //!< Has sig_final been computed
  mutable unsigned int sigi_final_ntrials; //!< Number of trials we did to determine sigi_final
  mutable double sigi_final; //!< Instrumental noise after smoothing/filtering.

  bool is_binned; //!< Has data been binned
  unsigned int nbins; //!< Number of bins
  unsigned int bin_sparcity; //!< Sampling rate of binning (1 means fully sampled)
  double bincent0; //!< Center of bin 0
  double bindelta; //!< Delta in bins
  unsigned int* binval; //!< Number of elements in each bin

  mutable ran rangen; //!< Random number generator

  // Generator for positions if not uniform
  bool use_clustered_pos; //!< Use clustered positions (rather than uniform)
  positionGeneratorClustered *posgen; //!< Position generator if not uniform

  bool is_full; //!< A realization has been generated

  double *data; //!< Data

  //Working variables
  mutable double *work; //!< Holds temporary array during convolution
  mutable double *gen_image; //!< Used to generate raw image if oversampling

  //Working arrays for holding beam
  unsigned int ngauss; //!< Number of elements in gauss
  double *gauss; //!< Holds convolution array

  //Additional smoothing arrays
  unsigned int ngauss_add; //!< Number of additional smoothing elements
  double *gauss_add; //!< Additional smoothing array

  /*! \brief Downsample array */
  void downSample(unsigned int, unsigned int, double* const,
                  unsigned int, unsigned int, double* const);

  /*! \brief Does model have valid params */
  bool isValid() const;

  void convolveInner(unsigned int, const double* const,
                     unsigned int, unsigned int, double* const,
                     double* const) const; //!< Inner convolution function

  /*! \brief Do convolution with beams in place*/
  void convolveWithBeamInPlace(unsigned int, unsigned int, double* const);
  /*! \brief Do convolution with beams, also possibly downsampling */
  void convolveWithBeam(unsigned int, unsigned int, double* const,
                        unsigned int, unsigned int, double* const);
  /*! \brief Do convolution with extra smoothing bit */
  void convolveWithAdd();

 public:
  /*!\brief Constructor */
  simImage(unsigned int N1, unsigned int N2, double PIXSIZE,
           double FWHM, double SIGI, double ESMOOTH=0.0,
           unsigned int OVERSAMPLE=1, unsigned int NBINS=1000,
           const std::string& powerfile="");
  ~simImage(); //!< Destructor

  /*! \brief Set random number generator seed */
  void setSeed(unsigned long long int seed) const { rangen.set_seed(seed); }

  /*! \brief Generate realization of model */
  void realize(const numberCounts& model, double n0,
               bool meansub=false, const fourierFilter* const filt=nullptr,
               bool bin=false, unsigned int sparsebin=1);

  bool isClustered() const { return use_clustered_pos; } //!< Are we using clustered positions?

  bool isBinned() const { return is_binned; } //!< Is data binned?
  void applyBinning(unsigned int=1); //!< Takes an unbinned image and bins it
  unsigned int binSparcity() const { return bin_sparcity; } //!< Bin sampling rate
  unsigned int getNBins() const { return nbins; } //!< Get number of bins
  double getBinCent0() const { return bincent0; } //!< Get bin 0 center
  double getBinDelta() const { return bindelta; } //!< Get bin size

  bool isHipassFiltered() const { return isHipass; }
  double getFiltScale() const { return filtscale; } //!< Get filtering scale
  double getFiltQFactor() const { return qfactor; } //!< Get filtering scale
  bool isMatchFiltered() const { return isMatched; }
  double getFiltFWHM() const { return matched_fwhm; }
  double getFiltSigInst() const { return matched_sigi; }
  double getFiltSigConf() const { return matched_sigc; }

  unsigned int isSmoothed() const { return esmooth > 0.0; }
  double getEsmooth() const;

  /*! \brief Get Raw instrument noise */
  double getInstNoise() const { return sigi; }

  /*! \brief Returns noise level estimate for image after smoothing or
    filtering */
  double getFinalNoise(unsigned int ntrials=3,
                       const fourierFilter* const filt=nullptr) const;

  double meanSubtract(); //!< Subtract mean from image
  void getMinMax(double& min, double& max) const; //!< Get minima and maxima of data
  double getMean() const; //!< Get the mean of the image
  void getMeanAndVar(double&, double&) const; //!< Get mean and variance

  double getBeamFWHM() const { return fwhm; } //!< Get FWHM (arcsec)
  double getBeamSum() const; //!< Returns sum of beam (i.e., area in pixels)
  double getBeamSumSq() const; //!< Returns sum of beam squared

  unsigned int getN1() const { return n1; } //!< Gets number of pixels in band1
  unsigned int getN2() const { return n2; } //!< Gets number of pixels in band2
  bool isOversampled() const { return oversample > 1; }
  unsigned int getOversampling() const { return oversample; } //!< Gets amount of oversampling
  double getPixSize() const { return pixsize; } //!< Gets pixel size (arcsec)
  const double* getData() const { return data; } //!< Gets data pointer
  const unsigned int* getBinnedData() const { return binval; } //!< Gets binned data

  double getArea() const; //!< Area of simulated image in sq deg

  /*! \brief Unchecked data access */
  double getData(unsigned int i, unsigned int j) const {return data[i*n2+j];}

  int writeToFits(const std::string& file) const; //!< Write as fits file

  int writeProbImageToFits(const std::string& file) const; //!< Write prob image as fits
  void writePositionGeneratorToHDF5Handle(hid_t obj_id) const;
};

#endif
