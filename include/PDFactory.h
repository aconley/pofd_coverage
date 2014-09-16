#ifndef __pdfactory__
#define __pdfactory__

#include<string>

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_errno.h>
#include<fftw3.h>

#include "../include/global_settings.h"
#include "../include/numberCounts.h"
#include "../include/PD.h"

/*!
  Always call initPD before using getPD for a given model.
 */

class PDFactory {
 private :
  bool rinitialized; //!< R is filled
  bool initialized; //!< forward transformed R is filled

  unsigned int currsize; //!< Current variable sizes 
  double sigma; //!< Current supported instrumental \f$\sigma\f$
  double max_n0; //!< Current maximum supported model \f$N_0\f$
  double base_n0; //!< Model base \f$N_0\f$
  double mn; //!< Expected mean for max n0 supported
  double varnoi; //!< Expected variance for max_n0 without instrument noise
  double sg; //!< Expected sigma, inc instrument noise, max n0 model

  //Working variables for transformation
  fftw_plan plan, plan_inv; //!< Hold plans
  bool plans_valid; //!< Are the current plans valid
  double* rvals; //!< Working space for R computation
  bool rdflux; //!< Is rvals R or R*dflux?
  fftw_complex *rtrans; //!< Holds forward transformed base R
  fftw_complex* pval; //!< Working variable holding p = exp( stuff )
  double* pofd; //!< Internal P(D) variable.

  double dflux; //!< Flux size step of last computation
  bool doshift; //!< Apply shifting
  double shift; //!< Shift amount

  unsigned int fftw_plan_style; //!< FFTW plan flags
  bool has_wisdom; //!< Has FFTW wisdom 
  std::string wisdom_file; //!< Name of wisdom file

  bool verbose; //!< Outputs information about each iter

  // Working variable for main R fill
  double *RFlux; //!< Holds R flux values for fill
  double minflux_R; //!< Minimum flux in RFlux; it wraps, so good to keep track fo
  
  void init(); //!< Initializes memory
  bool resize(unsigned int); //!< Sets internal storage to specified size
  
  /*! \brief Sets RFlux, with wrapping */
  void initRFlux(unsigned int n, double minflux, double maxflux);

  /*! \brief Figure out non-zero range of R */
  dblpair getMinMaxR(const numberCounts& model, const beamHist& bm) const;

  /*! \brief Get mean and variance from R integrals*/
  dblpair getRMoments(unsigned int n, const numberCounts&, const beamHist&,
		      dblpair range);

  /*! \brief Moves P(D) over to output variable inside getPD */
  void unwrapPD(double n0, unsigned int n, PD& pd) const;

#ifdef TIMING
  mutable std::clock_t RTime, p0Time, fftTime, posTime, copyTime, normTime;
  mutable std::clock_t meanTime, logTime, starttime;
#endif

 public :
  PDFactory(); //!< Default constructor
  PDFactory(const std::string&); //!< Constructor with wisdom file
  ~PDFactory(); //!< Destructor

  void setVerbose() { verbose = true; } //!< Sets verbose mode
  void unsetVerbose() { verbose = false; } //!< Unset verbose mode

  double getMaxN0() const { return max_n0; }
  unsigned int getCurrSize() const { return currsize; }
  double getSigma() const { return sigma; }
  double getBaseN0() const { return base_n0; }
  
  /*! \brief Adds wisdom file*/
  bool addWisdom(const std::string& filename);

  /*! \brief Initializes R*/
  void initR(unsigned int n, double minflux, double maxflux, 
	     const numberCounts& model, const beamHist& bm,
	     bool muldflux=false);

  /*! \brief Initializes P(D) by computing R */
  void initPD(unsigned int n, double inst_sigma, double maxflux, 
	      double maxn0, const numberCounts&, const beamHist&);

  /*! \brief Gets P(D) of specified transform size */
  void getPD(double, PD&, bool setLog=true);

  /*! \brief Write out current R to text file*/
  void writeRToFile(const std::string&) const;

  /*! \brief Write out current R to HDF5 file*/
  void writeRToHDF5(const std::string&) const;

#ifdef TIMING
  void resetTime();
  void summarizeTime(unsigned int=0) const;
#endif

};

#endif
