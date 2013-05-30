//binneddata

#ifndef __PDDouble__
#define __PDDouble__

#include<string>
#include<ostream>

#include "simImageDouble.h"
#include "global_settings.h"

/*!
  \brief Class to hold P(D) in 2 dimensions.  Supports interpolation.

  The interpolation is simple bilinear.  Only monotonic
  flux value grids are allowed.

  By default, the log of the P(D) is stored (and interpolated on),
  since this is what you want for the likelihood calculation.
  The user can specify the real-space P(D).
  The memory usage is something like std::vector -- each object
  has a capacity and a current size, with size < capacity.  This
  tries to avoid resizing if the memory is already allocated.
 */
class PDDouble {
 private:
  static const double lowsigval; //!< Constant used in edgeFix

  unsigned int n1; //!< Current size along dimension 1
  unsigned int n2; //!< Current size along dimension 2
  unsigned int capacity; //!< Current capacity

  double getLogLikeBinned(const simImageDouble&) const;
  double getLogLikeUnbinned(const simImageDouble&) const;

 public:
  PDDouble(unsigned int N1=0, FFTFLOAT MINFLUX1=0.0, FFTFLOAT DFLUX1=0.0,
	   unsigned int N2=0, FFTFLOAT MINFLUX2=0.0, FFTFLOAT DFLUX2=0.0,
	   bool LOG=true); //!< Constructor
  ~PDDouble(); //!< Destructor

  //Public for efficient filling -- bad form, but speed matters here
  bool logflat; //!< True if log( P(D1,D2) ) is stored instead of P(D1,D2)
  FFTFLOAT minflux1; //!< Minimum flux along axis1
  FFTFLOAT dflux1; //!< Flux step along axis1
  FFTFLOAT minflux2; //!< Minimum flux along axis2
  FFTFLOAT dflux2; //!< Flux step along axis2
  FFTFLOAT* pd_; //!< Actual P(D), stored as row-major array

  void shrink(); //!< Shrink memory requirements to user size
  void strict_resize(unsigned int,unsigned int); //!< Resize data arrays, forcing actual resizing in all cases
  void resize(unsigned int,unsigned int); //!< Resize data arrays

  double getTotal() const; //!< Get sum of entries
  double getIntegral() const; //!< Get integral of entries
  void normalize(); //!< Normalize integral of model
  
  bool isLog() const { return logflat; } //!< Is log(P(D)) stored?

  void applyLog(bool=false); //!< Logify if not already
  void deLog(); //!< de-logify if not already

  void edgeFix(bool donorm=true); //!< Apply Gaussian replacement to bottom edges

  void getMeans(double&, double&, bool donorm=true) const; //!< Get mean along each axis

  void getMeansAndVars(double&, double&, double&, double&, bool donorm=true) const; //!< Get means and variances along each axis
  
  PDDouble& operator=(const PDDouble&); //!< Copy

  /*! \brief Fill contents from row major array */
  void fill(unsigned int, FFTFLOAT, FFTFLOAT,
	    unsigned int, FFTFLOAT, FFTFLOAT,
	    const FFTFLOAT* const, bool LOG=true); 

  /*! \brief Get flux value (band 1) corresponding to index */
  FFTFLOAT getFluxVal1(unsigned int i) const { return minflux1+static_cast<FFTFLOAT>(i)*dflux1; }
  /*! \brief Get flux value (band 2) corresponding to index */
  FFTFLOAT getFluxVal2(unsigned int i) const { return minflux2+static_cast<FFTFLOAT>(i)*dflux2; }
  /*! \brief Get PD value corresponding to indices */
  FFTFLOAT getPDVal(unsigned int i, 
		  unsigned int j) const { return pd_[n2*i+j]; }
  /*! \brief Get PD value corresponding to flux values (using interpolation) */
  FFTFLOAT getPDVal(double,double,bool=false) const; //!< Interpolation

  /*! \brief PD element access */
  const FFTFLOAT* operator[](unsigned int i) const { return pd_+n2*i; }
  /*! \brief Number of elements in PD, band 1 */
  unsigned int getDim1() const { return n1; }
  /*! \brief Number of elements in PD, band 2 */
  unsigned int getDim2() const { return n2; }

  /*! \brief Get Log likelihood of data set*/
  double getLogLike(const simImageDouble&) const;

  std::ostream& writeToStream(std::ostream& os) const; //!< Write summary

  int writeToFits(const std::string& file) const; //!< Write as fits file
  
};

/*! \brief Output to stream operator */
std::ostream& operator<<(std::ostream& os, const PDDouble&);

#endif
