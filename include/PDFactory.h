#ifndef __pdfactory__
#define __pdfactory__

#include<string>

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_errno.h>

#include "global_settings.h"
#include "fftstuff.h"
#include "numberCounts.h"
#include "PD.h"

/*!
  Always call initPD before using getPD for a given model.
 */

class PDFactory {
 private :
  bool initialized; //!< forward transformed R is filled

  unsigned int currsize; //!< Current memory allocation
  FFTFLOAT sigma; //!< Current supported instrumental \f$\sigma\f$
  FFTFLOAT max_n0; //!< Current maximum supported model \f$N_0\f$
  FFTFLOAT base_n0; //!< Model base \f$N_0\f$
  FFTFLOAT mn; //!< Expected mean, base model
  FFTFLOAT sg; //!< Expected sigma, inc instrument noise, base model

  //Working variables for transformation
  unsigned int plan_size; //!< Size of plans
  FFTWPLAN plan, plan_inv; //!< Hold plans
  bool plans_valid; //!< Are the current plans valid
  FFTFLOAT* rvals; //!< Working space for R computation
  FFTWCOMPLEX *rtrans; //!< Holds forward transformed base R
  FFTWCOMPLEX* pval; //!< Working variable holding p = exp( stuff )
  FFTFLOAT* pofd; //!< Internal P(D) variable.  
  bool isRTransAllocated; //!< Is rtrans allocated

  FFTFLOAT dflux; //!< Flux size step of last computation
  unsigned int maxidx; //!< Max non-zero index in Rs
  bool doshift; //!< Apply shifting
  FFTFLOAT shift; //!< Shift amount

  unsigned int fftw_plan_style; //!< FFTW plan flags
  bool has_wisdom; //!< Has FFTW wisdom 
  std::string wisdom_file; //!< Name of wisdom file

  bool verbose; //!< Outputs information about each iter

  //Working variable for main R fill
  double *RFlux; //!< Holds R flux values for fill

  void init(); //!< Initializes memory
  bool resize(unsigned int); //!< Sets transform size arrays
  void strict_resize(unsigned int); //!< Sets transform size arrays
  
  //*! \brief Initializes R*/
  void initR(unsigned int n, FFTFLOAT maxflux, const numberCounts& model,
	     const beam& bm, double, double, unsigned int);

#ifdef TIMING
  std::clock_t RTime, p0Time, fftTime, posTime, copyTime, normTime, edgeTime;
  std::clock_t meanTime, logTime;
#endif

 public :
  PDFactory(); //!< Default constructor
  PDFactory(const std::string&); //!< Constructor with wisdom file
  ~PDFactory(); //!< Destructor

  void setVerbose() { verbose = true; } //!< Sets verbose mode
  void unsetVerbose() { verbose = false; } //!< Unset verbose mode

  unsigned int getPlanSize() const { return plan_size; }
  FFTFLOAT getMaxN0() const { return max_n0; }

  /*! \brief Adds wisdom file*/
  bool addWisdom(const std::string& filename);

  /*! \brief Initializes P(D) by computing R */
  void initPD(unsigned int, FFTFLOAT, FFTFLOAT, FFTFLOAT,
	      const numberCounts&, const beam&, 
	      double, double, unsigned int);

  /*! \brief Gets P(D) of specified transform size */
  void getPD(double, PD&, bool setLog=true, 
	     bool edgeFix=false);

  /*! \brief Write out current R to text file*/
  void writeRToFile(const std::string&) const;

  /*! \brief Get first n integrals of R*/
  void getRIntegrals(unsigned int, FFTFLOAT, const numberCounts&,
		     const beam&, double, double, unsigned int,
		     std::vector<FFTFLOAT>& vec, 
		     unsigned int nmom=7);

#ifdef TIMING
  void resetTime();
  void summarizeTime(unsigned int=0) const;
#endif

};

#endif
