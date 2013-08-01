//numberCounts.h

#ifndef __numberCounts__
#define __numberCounts__

#include<string>
#include<ostream>

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_integration.h>

#include "../include/global_settings.h"
#include "../include/beam.h"

/*!
  \brief Spline based law number counts
*/

class numberCounts {
 private :

  static const double ftol; //!< General floating point tolerance for equality

  // Internally we represent the number counts as spline in log/log space
  double base_n0; //!< Total number of sources in base model
  unsigned int nknots; //!< Number of flux density knots
  double* knotpos; //!< Positions of knots, length nknots
  double* logknotpos; //!< Log2 of knot positions, length nknots
  double* knotvals; //!< Values of differential number counts at knotpos
  double* logknotvals; //!< Log2 values of differential number counts at knotpos

  gsl_interp_accel *acc; //!< Spline lookup accelerator
  gsl_spline *splinelog; //!< Spline in log/log space

  gsl_integration_workspace *gsl_work; //!< Integration workspace for QAG
  /*! \brief Integrates flux^power over number counts */
  double splineInt(double power) const; 
  void **varr; //!< Internal evil casting array for integration

  // Stuff for generating sources
  unsigned int gen_ninterp; //!< Number of generated interpolation sources
  double *gen_interp_flux; //!< Flux densities in interpolation
  double *gen_interp_cumsum; //!< Cumulative probability value
  gsl_interp *gen_interp; //!< Linear interpolant function
  gsl_interp_accel *gen_interp_acc; //!< Accelerator
  
  double base_flux; //!< Flux per area for base model
  double base_fluxsq; //!< Flux squared per area for base model

  //Working variables for histogrammed beam
  mutable unsigned int nbm; //!< Number of elements in inverse beam
  mutable unsigned int* bm_wts; //!< Histogrammed inverse beam weights
  mutable double* inv_bm; //!< Histogrammed inverse beam

 public:
  explicit numberCounts(const std::string&, unsigned int=2000);  //!< Constructor with model file
  ~numberCounts();

  bool isValid() const; //!< See if model params are valid
  
  /*! \brief Get number of sources per area in base model*/
  double getBaseN0() const;

  /*! \brief Get flux density per unit area for base model*/
  double getBaseFluxPerArea() const;

  /*! \brief Get flux density squared per unit area for base model*/
  double getBaseFluxSqPerArea() const;

  /*!\brief Get differential number counts for base model*/
  double getdNdS(double) const;

  /*!\brief Get number of source responses for base model */
  double getR(double, const beam&, double, double,
	      unsigned int) const;
  /*!\brief Get number of source responses for base model */
  double getR(double, const beam&, double, double,
	      unsigned int, unsigned int) const;

  /*!\brief Get number of source responses for base model, general case, array*/
  void getR(unsigned int, double, double, const beam&, 
	    double, double, unsigned int, double*) const;
  /*!\brief Get number of source responses for base model, general case, array*/
  void getR(unsigned int, double, double, const beam&, 
	    double, double, unsigned int, unsigned int, double*) const;
  
  /*! \brief Generate a source flux from model */
  double genSource(double val) const;

  /*! \brief Get minimum knot position*/
  double getMinKnotPosition() const { return knotpos[0]; }

  /*! \brief Get maximum knot position*/
  double getMaxKnotPosition() const { return knotpos[nknots-1]; }

  /*! \brief Get number of knots */
  unsigned int getNKnots() const { return nknots; }

  /*! \brief Get Knot position */
  double getKnotPosition(unsigned int idx) const { return knotpos[idx]; }

  /*! \brief Get Log10 differential number counts knot value */
  double getLog10KnotValue(unsigned int idx) const { 
    return logknotvals[idx] * pofd_coverage::ilogfac; }

  bool writeToStream(std::ostream& os) const; //!< Output to stream

};

/*! Stream output operator */
std::ostream& operator<<(std::ostream& os, const numberCounts& b);

#endif
