//simImageDouble.h

#ifndef __simImageDouble__
#define __simImageDouble__

#include<string>
#include<utility>

#include "../include/numberCountsDouble.h"
#include "../include/ran.h"
#include "../include/positionGenerator.h"
#include "../include/fourierFilter.h"

/*!
  \brief Data class for creating and holding double (two-band) simulated
  images with Gaussian beams
 */

class simImageDouble {
 private:
  unsigned int n1; //!< Number of pixels along dimension 1.  Same in both bands
  unsigned int n2; //!< Number of pixels along dimension 2.  Same in both bands
  double pixsize; //!< Size of pixels in arcsec.  Assumed square.

  unsigned int oversample; //!< Oversampling in generation
  unsigned int ngen1; //!< Number of generated pixels in dimension 1
  unsigned int ngen2; //!< Number of generated pixels in dimension 2
  double pixsize_gen; //!< Size of pixels in internal mode

  double fwhm1; //!< Beam FWHM in arcsec, band 1
  double fwhm2; //!< Beam FWHM in arcsec, band 2
  double sigi1; //!< Instrumental (white) noise in band 1
  double sigi2; //!< Instrumental (white) noise in band 2
  double esmooth1; //!< Additional Gaussian smoothing FWHM in arcsec, band 1
  double esmooth2; //!< Additional Gaussian smoothing FWHM in arcsec, band 2
  mutable fourierFilter* filt; //!< High pass filter

  // The sigi_final are only computed when needed, since they isn't 
  // totally free to do if there is filtering
  mutable bool sigi_final_computed; //!< Has sig_final been computed
  mutable unsigned int sigi_final_ntrials; //!< Number of trials we did to determine sigi_final
  mutable double sigi_final1; //!< Instrumental noise after smoothing/filtering, band 1
  mutable double sigi_final2; //!< Instrumental noise after smoothing/filtering, band 2

  bool is_binned; //!< Has data been binned
  unsigned int nbins; //!< Number of bins along each dimension
  unsigned int bin_sparcity; //!< Sampling rate of binning (1 means fully sampled)
  double bincent01, bincent02; //!< Center of bin 0 along dimension 1,2
  double bindelta1, bindelta2; //!< Delta in bins along each dimension
  unsigned int* binval; //!< Number of elements in each bin (row-major order)

  mutable ran rangen; //!< Random number generator

  // Generator for positions if not uniform
  bool use_clustered_pos; //!< Use clustered positions (rather than uniform)
  positionGeneratorClustered *posgen; //!< Position generator if not uniform

  bool is_full; //!< A realization has been generated

  double *data1; //!< Data in band 1, final
  double *data2; //!< Data in band 2, final
  
  //Working variables
  mutable double *work; //!< Temporary array to hold convolution product
  mutable double *gen_1; //!< Holds generated image, band 1 if oversample > 1
  mutable double *gen_2; //!< Holds generated image, band 2 if oversample > 1

  //Working arrays for holding beam
  unsigned int ngauss1; //!< Number of elements in gauss1
  double *gauss1; //!< Holds convolution array, band 1
  unsigned int ngauss2; //!< Number of elements in gauss2
  double *gauss2; //!< Holds convolution array, band 2

  //Additional smoothing arrays
  unsigned int ngauss_add1; //!< Number of additional smoothing elements, band1
  double *gauss_add1; //!< Additional smoothing array
  unsigned int ngauss_add2; //!< Number of additional smoothing elements, band1
  double *gauss_add2; //!< Additional smoothing array

  /*! \brief Inner convolution function */
  void convolveInner(double*, unsigned int, const double* const);

  /*! \brief Downsample arrays to final resolution */
  void downSample(unsigned int, unsigned int, double* const,
		  unsigned int, unsigned int, double* const); 

  double getFinalNoiseHelper(unsigned int ntrials, double* const data, 
			     double sigi, double fwhm, double esmooth, 
			     unsigned int ngauss_add,
			     const double* const gauss_add,
			     fourierFilter* const filt) const;

  //Convolution stuff
  void convolveInner(unsigned int, const double* const,
		     unsigned int, unsigned int, double* const,
		     double* const) const; //!< Inner convolution function

  /*! \brief Do convolution with beams */
  void convolveWithBeam();
  /*! \brief Do convolution with extra smoothing bit */
  void convolveWithAdd();

  /*! \brief Does model have valid params */
  bool isValid() const;

  /*! \brief Fits writer helper */
  int writeFits(const std::string& outfile, unsigned int band) const;

 public:

  simImageDouble(unsigned int N1, unsigned int N2, double PIXSIZE, 
		 double FWHM1, double FWHM2, double SIGI1,
		 double SIGI2, double ESMOOTH1=0.0, double ESMOOTH2=0.0, 
		 double FILTERSCALE=0.0, unsigned int OVERSAMPLE =1, 
		 unsigned int NBINS=1000, const std::string& powerspec="",
		 bool quickfft=false);
  ~simImageDouble(); //!< Destructor

  /*! \brief Set random number generator seed */
  void setSeed(unsigned long long int seed) const { rangen.set_seed(seed); }
  
  /*! Generate realization of model */
  void realize(const numberCountsDouble& model, double n0,
	       bool meansub=false, bool bin=false, unsigned int sparsebin=1);

  bool isClustered() const { return use_clustered_pos; } //!< Are we using clustered positions?

  bool isBinned() const { return is_binned; } //!< Is data binned?
  void applyBinning(unsigned int=1); //!< Takes an unbinned image and bins it
  unsigned int binSparcity() const { return bin_sparcity; } //!< Bin sampling rate
  unsigned int getNBins() const { return nbins; } //!< Get number of bins
  double getBinCent01() const { return bincent01; }
  double getBinDelta1() const { return bindelta1; }
  double getBinCent02() const { return bincent02; }
  double getBinDelta2() const { return bindelta2; }

  bool isFiltered() const { return filt != NULL; } //!< Is the image filtered
  double getFiltScale() const; //!< Get image filter scale

  bool isSmoothed() const { return (esmooth1 > 0.0) || (esmooth2 > 0.0); }
  std::pair<double, double> getEsmooth() const {
    return std::make_pair(esmooth1, esmooth2); }

  /*! \brief Get raw instrument noise */
  std::pair<double,double> getInstNoise() const {
    return std::make_pair(sigi1, sigi2); }

  /*! \brief Returns noise level estimate after smoothing or filtering*/
  std::pair<double,double> getFinalNoise(unsigned int ntrials=3) const; 

  std::pair<double, double> meanSubtract(); //!< Subtract off means
  std::pair<double, double> getMean() const; //!< Get mean in each band
  /*! \brief Get mean and variance of image in each band */
  void getMeanAndVar(double&, double&, double&, double&) const;
  /*! \brief Get minimum and maximum of image in each band */
  void getMinMax(double& min1, double& max1, double& min2, 
		 double& max2) const; //!< Get minima and maxima in both bands

  /*! \brief Returns sum of beams (i.e., area in pixels) */
  std::pair<double,double> getBeamSum() const; 
  /*! \brief Returns sum of beams squared (i.e., area in pixels of bm^2) */
  std::pair<double,double> getBeamSumSq() const; 

  unsigned int getN1() const { return n1; } //!< Get image extent, band 1
  unsigned int getN2() const { return n2; } //!< Get image extent, band 1
  bool isOversampled() const { return oversample > 1; }
  double getOversampling() const { return oversample; } //!< Gets amount of oversampling

  double getPixSize() const { return pixsize; } //!< Get pixel size (arcsec)
  const double* const getData1() const { return data1; } //!< Band 1 data access
  const double* const getData2() const { return data2; } //!< Band 2 data access
  const unsigned int* const getBinnedData() const { return binval; } //!< Binned data direct access

  double getArea() const; //!< Area of images in sq deg

  /*! \brief Unchecked data access, band 1 (unbinned) */
  double getData1(unsigned int i, unsigned int j) const {return data1[i*n2+j];}
  /*! \brief Unchecked data access, band 2 (unbinned) */
  double getData2(unsigned int i, unsigned int j) const {return data2[i*n2+j];}

  int writeToFits(const std::string& outfile1,
		  const std::string& outfile2) const; //!< Write as fits files

};

#endif
