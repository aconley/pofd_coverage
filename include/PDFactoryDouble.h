#ifndef __pdfactorydouble__
#define __pdfactorydouble__

#include<string>

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_errno.h>
#include<fftw3.h>

#include "../include/numberCountsDouble.h"
#include "../include/PDDouble.h"

/*!
  \brief Class for computing P(D) values from a set of parameters,
  two band version

  Only square FFTs are supported.
  
  Always call initPD before using getPD for a given model.
  Unlike the one band version, doesn't interpolate on R.  Two dimensional
  interpolation was either inaccurate or not faster than directly computing
  R in tests.
 */

class PDFactoryDouble {
 private :
  //Constants controlling edge integrations
  static const double lowEdgeRMult; //!< How far down to go from model limits in edge integrals
  static const bool use_edge_log_x; //!< Integrate Rx in x log space
  static const bool use_edge_log_y; //!< Integrate Ry in y log space

  bool initialized; //!< forward transformed R is filled

  unsigned int currsize; //!< Current memory allocation
  double sigma1; //!< Current supported instrumental \f$\sigma\f$, band 1
  double sigma2; //!< Current supported instrumental \f$\sigma\f$, band 2
  double max_n0; //!< Current maximum supported model \f$N_0\f$
  double base_n0; //!< Model base \f$N_0\f$
  double mn1; //!< Expected mean, band 1
  double mn2; //!< Expected mean, band 2
  double var_noi1; //!< Expected variance without instrumental noise, band 1
  double var_noi2; //!< Expected variance without instrumental noise, band 2
  double sg1; //!< Expected sigma (inc instrument noise), band 1
  double sg2; //!< Expected sigma (inc instrument noise), band 2

  unsigned int plan_size; //!< Size of currently computed plans
  fftw_plan plan;     //!< Holds forward transformation plan
  fftw_plan plan_inv; //!< Holds inverse transformation plan

  //Working variables for transformation
  bool rvars_allocated; //!< Are R variables (rest of this block) allocated
  bool is_r_set; //!< Has R been set to something?
  unsigned int rsize; //!< Size of R fill along each dimension
  double *RFlux1; //!< Holds R flux values for fill
  double *RFlux2; //!< Holds R flux values for fill
  double minfluxR_1; //!< Minimum flux in RFlux1
  double minfluxR_2; //!< Minimum flux in RFlux2
  double* rvals; //!< Working space for R computation, row major order
  fftw_complex *rtrans; //!< Holds FFTed rvals 
  fftw_complex* pval; //!< Working variable holding p = exp( stuff )
  double* pofd; //!< Internal P(D) variable.

  void allocateRvars(); //!< Allocates R variables if needed
  void freeRvars(); //!< Free R variables

  //Edge variables
  bool edgevars_allocated; //!< Are Edge variables (this block) allocated
  unsigned int nedge; //!< Number of edge integral steps
  double* REdgeFlux1; //!< Holds flux for R edge integration
  double* REdgeFlux2; //!< Holds flux for R edge integration
  double* REdgeWork; //!< Holds R in edge integration
    
  void allocateEdgevars(); //!< Allocate Edge variables
  void freeEdgevars(); //!< Free edge variables

  double dflux1; //!< Flux size step of last computation, band 1
  double dflux2; //!< Flux size step of last computation, band 2
  bool doshift1; //!< Apply shifting, band 2
  bool doshift2; //!< Apply shifting, band 2
  double shift1; //!< Shift amount, band 2
  double shift2; //!< Shift amount, band 2
  
  unsigned int fftw_plan_style; //!< FFTW plan flags
  bool has_wisdom; //!< Has FFTW wisdom 
  std::string wisdom_file; //!< Name of wisdom file

  bool verbose; //!< Outputs information about each iter

  void init(unsigned int); //!< Initializes memory
  bool resize(unsigned int); //!< Sets transform size arrays
  void strict_resize(unsigned int); //!< Sets transform size arrays
  void free(); //!< Frees memory

  /*! \brief Sets RFlux1, Rflux2, with wrapping */
  void initRFlux(unsigned int n, double minflux1, double maxflux1,
		 double minflux2, double maxflux2);

  /*! \brief Initializes R */
  void initR(unsigned int n, double minflux1, double maxflux1, 
	     double minflux2, double maxflux2, const numberCountsDouble& model,
	     const doublebeamHist& bm, bool setEdge=true);

  /*! \brief Figure out min/max flux densities where R is non-zero */
  std::pair<dblpair, dblpair> getMinMaxR(const numberCountsDouble&,
					 const doublebeamHist&) const;

  //*! \brief Compute means and variances from R */
  std::pair<dblpair, dblpair> getRMoments(unsigned int n, 
					  const numberCountsDouble&,
					  const doublebeamHist&,
					  const std::pair<dblpair, dblpair>&,
					  bool setEdge=true);

  /*! \brief Moves P(D) over to output variable inside getPD */
  void unwrapPD(unsigned int n, PDDouble& pd) const;

#ifdef TIMING
  std::clock_t RTime, p0Time, fftTime, posTime, copyTime;
  std::clock_t normTime, meanTime, logTime;
#endif

 public :

  PDFactoryDouble(unsigned int nedge=256); //!< Default constructor
  PDFactoryDouble(const std::string&, unsigned int nedge=256); //!< Constructor with wisdom file
  ~PDFactoryDouble(); //!< Destructor

  void setVerbose() { verbose = true; } //!< Sets verbose mode
  void unsetVerbose() { verbose = false; } //!< Unset verbose mode

  double getMaxN0() const { return max_n0; }

  /*! \brief Get size of last FFT */
  unsigned int getPlanSize() const { return plan_size; }

  /*! \brief Returns edge integration length */
  unsigned int getNEdge() const { return nedge; }
  void setNEdge(unsigned int); //!< Set nedge

  /*! \brief Adds FFTW wisdom file*/
  bool addWisdom(const std::string& filename);

  /*! \brief Initializes P(D) by computing R and forward transforming it*/
  void initPD(unsigned int n, double inst_sigma1, double inst_sigma2, 
	      double maxflux1, double maxflux2, double maxn0,
	      const numberCountsDouble& model, const doublebeamHist& bm,
	      bool setEdge=true);

  /*! \brief Gets P(D) with specified noise levels */
  void getPD(double, PDDouble&, bool setLog=true);

  void writeRToFile(const std::string&) const;
  
#ifdef TIMING
  void resetTime(); //!< Reset timing information
  void summarizeTime(unsigned int=0) const; //!< Summarize timing information
#endif
};

#endif
