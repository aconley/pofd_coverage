//numberCounts.h

#ifndef __numberCounts__
#define __numberCounts__

#include<string>
#include<ostream>

#include "../include/global_settings.h"
#include "../include/beam.h"

/*!
  \brief Broken power law number counts
*/

class numberCounts {
 private :

  static const double ftol; //!< General floating point tolerance for equality

  // Internally we represent the number counts as a_i S^-gamma_i
  double base_n0; //!< Total number of sources in base model
  unsigned int nknots; //!< Number of flux density knots
  double* knotpos; //!< Positions of knots, length nknots
  double* knotvals; //!< Values of differential number counts at knotpos
  double* logknotvals; //!< Log values of differential number counts at knotpos

  //Internal information
  void initMParams(); //< Initialize internal arrays based on model
  double* a; //!< Model a parameter for N_0 = base_n0.  Length nknots-1
  double* gamma; //!< Model gammas.  Length nknots-1
  bool *gamone; //!< Is gamma 1?  len nknots-1
  double *fk; //!< Cumulative source numbers, len nknots-1
  double *omg; //!< 1.0 - gamma, len knots-1
  double *iomg; //!< 1.0 / (1.0 - gamma), len nknots-1
  double *powarr; //!< Internal convenience array, len nknots-1

  double base_meanflux; //!< Mean flux per area for base model
  double base_meanfluxsq; //!< Mean flux squared per area for base model

  //Working variables for histogrammed beam
  mutable unsigned int nbm; //!< Number of elements in inverse beam
  mutable unsigned int* bm_wts; //!< Histogrammed inverse beam weights
  mutable double* inv_bm; //!< Histogrammed inverse beam

 public:
  explicit numberCounts(const std::string&);  //!< Constructor with model file
  ~numberCounts();

  bool isValid() const; //!< See if model params are valid
  
  /*! \brief Get number of sources per area in base model*/
  double getBaseN0() const;

  /*! \brief Get Mean Flux per unit area for base model*/
  double getMeanFluxPerArea() const;

  /*! \brief Get Mean Flux squared per unit area for base model*/
  double getMeanFluxSqPerArea() const;

  /*!\brief Get differential number counts for base model*/
  double getdNdS(double) const;

  /*!\brief Get number of source responses for base model */
  double getR(double, const beam&, double, double,
	      unsigned int) const;

  /*!\brief Get number of source responses for base model, general case, array*/
  void getR(unsigned int, double, double, const beam&, 
	    double, double, unsigned int, double*) const;
  
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
