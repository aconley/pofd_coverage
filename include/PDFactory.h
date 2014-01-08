#ifndef __pdfactory__
#define __pdfactory__

#include<string>

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_errno.h>
#include<fftw3.h>

#include "../include/numberCounts.h"
#include "../include/PD.h"

/*!
  Always call initPD before using getPD for a given model.
 */

class PDFactory {
 private :
  bool initialized; //!< forward transformed R is filled

  unsigned int currsize; //!< Current memory allocation
  double sigma; //!< Current supported instrumental \f$\sigma\f$
  double max_n0; //!< Current maximum supported model \f$N_0\f$
  double base_n0; //!< Model base \f$N_0\f$
  double mn; //!< Expected mean, base model
  double sg; //!< Expected sigma, inc instrument noise, base model

  //Working variables for transformation
  unsigned int plan_size; //!< Size of plans
  fftw_plan plan, plan_inv; //!< Hold plans
  bool plans_valid; //!< Are the current plans valid
  double* rvals; //!< Working space for R computation
  fftw_complex *rtrans; //!< Holds forward transformed base R
  fftw_complex* pval; //!< Working variable holding p = exp( stuff )
  double* pofd; //!< Internal P(D) variable.  
  bool isRTransAllocated; //!< Is rtrans allocated

  double dflux; //!< Flux size step of last computation
  double minflux_R; //!< Minimum flux in RFlux
  bool doshift; //!< Apply shifting
  double shift; //!< Shift amount

  unsigned int fftw_plan_style; //!< FFTW plan flags
  bool has_wisdom; //!< Has FFTW wisdom 
  std::string wisdom_file; //!< Name of wisdom file

  bool verbose; //!< Outputs information about each iter

  //Working variable for main R fill
  double *RFlux; //!< Holds R flux values for fill

  void init(); //!< Initializes memory
  bool resize(unsigned int); //!< Sets transform size arrays
  void strict_resize(unsigned int); //!< Sets transform size arrays
  
  /*! \brief Sets RFlux, with wrapping */
  void initRFlux(unsigned int n, double minflux, double maxflux);

  /*! \brief Initializes R*/
  void initR(unsigned int n, double minflux, double maxflux, 
	     const numberCounts& model, const beamHist& bm);

  /*! \brief Get mean and variance from R integrals*/
  std::pair<double, double> getRMoments(unsigned int n, const numberCounts&, 
					const beamHist&);

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
  double getMaxN0() const { return max_n0; }
  unsigned int getCurrSize() const { return currsize; }
  double getSigma() const { return sigma; }
  double getBaseN0() const { return base_n0; }
  
  /*! \brief Adds wisdom file*/
  bool addWisdom(const std::string& filename);

  /*! \brief Initializes P(D) by computing R */
  void initPD(unsigned int n, double inst_sigma, double maxflux, 
	      double maxn0, const numberCounts&, const beamHist&);

  /*! \brief Gets P(D) of specified transform size */
  void getPD(double, PD&, bool setLog=true, 
	     bool edgeFix=false);

  /*! \brief Write out current R to text file*/
  void writeRToFile(const std::string&) const;

#ifdef TIMING
  void resetTime();
  void summarizeTime(unsigned int=0) const;
#endif

};

#endif
