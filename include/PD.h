//binneddata

#ifndef __PD__
#define __PD__

#include<string>
#include<ostream>

#include<simImage.h>

/*!
  \brief Class to hold P(D).  Supports interpolation.

  Only monotonic flux value grids are allowed.

  By default, the log of the P(D) is stored (and interpolated on),
  since this is what you want for the likelihood calculation.
  The user can specify the real-space P(D).
  The memory usage is something like std::vector -- each object
  has a capacity and a current size, with size < capacity.  This
  tries to avoid resizing if the memory is already allocated.
 */
class PD {
 private:
  static const double lowsigval; //!< Constant used in edgeFix

  unsigned int n; //!< Current size
  unsigned int capacity; //!< Current capacity

  double getLogLikeBinned(const simImage&) const;
  double getLogLikeUnbinned(const simImage&) const;

 public:
  PD(unsigned int N=0, double MINFLUX=0.0, double DFLUX=0.0,
     bool LOG=true); //!< Constructor
  ~PD(); //!< Destructor

  //Public for efficient filling -- bad form, but speed matters here
  bool logflat; //!< True if log( P(D1,D2) ) is stored instead of P(D1,D2)
  double minflux; //!< Minimum flux
  double dflux; //!< Flux step along axis
  double* pd_; //!< Actual P(D)

  void shrink(); //!< Shrink memory requirements to user size
  void strict_resize(unsigned int); //!< Resize data arrays, forcing actual resizing in all cases
  void resize(unsigned int); //!< Resize data arrays

  double getTotal() const; //!< Get sum of entries
  double getIntegral() const; //!< Get integral of entries
  void normalize(); //!< Normalize integral of model
  
  bool isLog() const { return logflat; } //!< Is log(P(D)) stored?

  void applyLog(bool=false); //!< Logify if not already
  void deLog(); //!< de-logify if not already

  void edgeFix(bool donorm=true); //!< Apply Gaussian replacement to bottom edges

  void getMean(double&, bool donorm=true) const; //!< Get mean 

  void getMeanAndVar(double&, double&, bool donorm=true) const; //!< Get mean and variance 
  
  PD& operator=(const PD&); //!< Copy

  /*! \brief Fill contents from array*/
  void fill(unsigned int, double, double,
	    const double* const, bool LOG=true); 

  double getFluxVal(unsigned int i) const { return minflux+static_cast<double>(i)*dflux; }
  double getPDVal(unsigned int i) const { return pd_[i]; }
  double getPDVal(double) const; //!< Interpolation

  const double operator[](unsigned int i) const { return pd_[i]; }
  unsigned int getDim() const { return n; }

  /*! \brief Get Log likelihood of data set*/
  double getLogLike(const simImage&) const;

  std::ostream& writeToStream(std::ostream& os) const; //!< Write summary

  int writeToFits(const std::string& file) const; //!< Write as fits file
  
};

std::ostream& operator<<(std::ostream& os, const PD&);

#endif
