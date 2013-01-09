//simImage.h

#ifndef __simImage__
#define __simImage__

#include<string>

#include<ran.h>

#include<numberCounts.h>

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
  double bm_pixarea; //!< Beam area in pixels
  double sigi; //!< Instrumental noise
  double esmooth; //!< Additional Gaussian smoothing FWHM in arcsec

  bool is_binned; //!< Has data been binned
  unsigned int nbins; //!< Number of bins
  double bincent0; //!< Center of bin 0
  double bindelta; //!< Delta in bins
  unsigned int* binval; //!< Number of elements in each bin

  mutable ran rangen; //!< Random number generator

  bool is_full; //!< A realization has been generated

  double *data; //!< Data
  
  //Working variables
  double *work; //!< Holds temporary array during convolution
  double *work2; //!< Target of final convolution

  //Working arrays for holding beam
  unsigned int ngauss; //!< Number of elements in gauss
  double *gauss; //!< Holds convolution array

  //Additional smoothing arrays
  unsigned int ngauss_add; //!< Number of additional smoothing elements
  double *gauss_add; //!< Additional smoothing array
  bool smooth_applied; //!< Was esmooth applied to this fill

  /*! \brief Downsample array */
  void downSample(unsigned int, unsigned int, double* const,
		  unsigned int, unsigned int, double* const); 
  void convolveInner(unsigned int, const double* const,
		     unsigned int, unsigned int, double* const,
		     double* const); //!< Inner convolution function

  /*! \brief Does model have valid params */
  bool isValid() const;

  /*! \brief Do convolution with beams */
  void convolveWithBeam(unsigned int, unsigned int, double* const,
			unsigned int, unsigned int, double* const);
  /*! \brief Do convolution with extra smoothing bit */
  void convolveWithAdd();

 public:

  simImage(unsigned int, unsigned int, double,
	   double, double, double=0.0, unsigned int=1, 
	   unsigned int=1000); //!< Constructor with parameters
  ~simImage(); //!< Destructor

  /*! \brief Set random number generator seed */
  void setSeed(unsigned long long int seed) const { rangen.set_seed(seed); }
  
  void realize(const numberCounts&, double, bool=false,bool=false,
	       bool=false); //!< Generate realization of model

  bool isBinned() const { return is_binned; } //!< Is data binned?
  void applyBinning(); //!< Takes an unbinned image and bins it
  unsigned int getNBins() const { return nbins; } //!< Get number of bins
  double getBinCent0() const { return bincent0; } //!< Get bin 0 center
  double getBinDelta() const { return bindelta; } //!< Get bin size

  /*! \brief Noise in current simulated image */
  double getNoise() const;

  double getSmoothedNoiseEstimate() const; //!< Returns noise level estimate for smoothed image

  double meanSubtract(); //!< Subtract mean from image
  void getMinMax(double& min, double& max) const; //!< Get minima and maxima of data
  double getMean() const; //!< Get the mean of the image
  void getMeanAndVar(double&, double&) const; //!< Get mean and variance

  double getBeamSum() const; //!< Returns sum of beam (i.e., area in pixels)
  double getBeamSumSq() const; //!< Returns sum of beam squared

  unsigned int getN1() const { return n1; } //!< Gets number of pixels in band1
  unsigned int getN2() const { return n2; } //!< Gets number of pixels in band2
  double getOversampling() const { return oversample; } //!< Gets amount of oversampling
  double getPixSize() const { return pixsize; } //!< Gets pixel size (arcsec)
  const double* const getData() const { return data; } //!< Gets data pointer
  const unsigned int* const getBinnedData() const { return binval; } //!< Gets binned data

  double getArea() const; //!< Area of simulated image in sq deg

  /*! \brief Unchecked data access */
  double getData(unsigned int i, unsigned int j) const {return data[i*n2+j];}

  int writeToFits(const std::string& file) const; //!< Write as fits file

};

#endif
