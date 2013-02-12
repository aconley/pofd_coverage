//numberCountsDouble.h

#ifndef __numberCountsDouble__
#define __numberCountsDouble__

#include<string>
#include<fstream>
#include<utility>

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_integration.h>

#include<doublebeam.h>
#include<global_settings.h>

/*!
  \brief Broken power law times log normal number counts
  \ingroup Models

  The full expression is
    \f[
       \frac{dN}{dS_1\, dS_2} = \frac{n_1\left(S_1 \right)}{S_1} 
            \mathrm{L} \left( \frac{S_2}{S_1}; \mu\left(S_1\right),
            \sigma\left(S_1\right) \right)
    \f] 
    where \f$n_1\f$ is just as in the one-dimensional case represented
    by numberCounts.  The division by \f$S_1\f$ may seem odd, but assures 
    that the distribution marginalized over \f$S_2\f$ is just \f$ n_1 \f$.
    \f$\sigma\f$ and \f$\mu\f$ are also strored as splines. \f$L\f$
    is just the Log-Normal form:
    \f[
      \mathrm{L}\left( x ; \mu, \sigma \right) =
       \frac{1}{\sqrt{2 \pi} \sigma} \frac{1}{x}
       \exp\left[ \frac{ - \left(\log x - \mu\right)^2 }{ \sigma^2 } \right] .
    \f]
    Note that \f$\sigma\f$ and \f$\mu\f$ are not the mean and square root
    of the variance of the actual distribution, but rather the mean
    and square root of the variance of the log quantities.  These
    are related by
    \f[
       \left< S_2/S_1 \right> = \exp \left[ \mu\left(S_1\right) +
             \frac{1}{2}\sigma^2\left(S_1\right) \right]
    \f]
    and
    \f[
       Var[ S_2/S_1 ] = \exp \left( 2 \mu \left( S_1 \right) +
           \sigma^2\left(S_1\right) \right) 
          \left[ \exp \left( \sigma^2\left(S_1\right) \right) - 1 \right]
    \f].
    Note that the number counts require that \f$S_2/S_1 > 0\f$.
    We also explicitly require that \f$S_1 > 0\f$.  Both are justified
    on physical grounds.
*/
class numberCountsDouble {

 private:

  static const double ftol; //!< General floating point tolerance for equality

  double base_n0; //!< Total number of sources in base model

  //Band 1
  unsigned int nknots; //!< Number of knot positions for counts in band one
  double *knotpos; //!< Locations of knots, band 1
  double* knotvals; //!< Values of differential number counts at knotpos, 
  double *logknotvals; //!< Ln values of differential number counts at knots, band 1

  //LogNormal for band 2
  unsigned int nsigma; //!< Number of knot positions in sigma
  double *sigmapos; //!< Positions of sigma knots
  double *sigmavals; //!< Sigma values at positions of knots
  gsl_interp *sigmainterp; //!< Sigma interpolator
  gsl_interp_accel *accsigma; //!< Sigma accelerator
  unsigned int noffset; //!< Number of knot positions in offset
  double *offsetpos; //!< Positions of offset knots
  double *offsetvals; //!< Offset values at positions of knots
  gsl_interp *offsetinterp; //!< Offset interpolator
  gsl_interp_accel *accoffset; //!< Offset accelerator

  //Internal information, band 1
  void initM1Params(); //< Initialize internal arrays based on model, band 1
  double* a; //!< Model a parameter for N_0 = base_n0.  Length nknots-1
  double* gamma; //!< Model gammas.  Length nknots-1
  bool *gamone; //!< Is gamma 1?  len nknots-1
  double *fk; //!< Cumulative source numbers, len nknots-1
  double *omg; //!< 1.0 - gamma, len knots-1
  double *iomg; //!< 1.0 / (1.0 - gamma), len nknots-1
  double *powarr; //!< Internal convenience array, len nknots-1

  double base_meanflux1; //!< Mean flux per area for base model, band 1
  double base_meanfluxsq1; //!< Mean flux squared per area for base model, b1
  double base_meanflux2; //!< Mean flux per area for base model, band 2
  double base_meanfluxsq2; //!< Mean flux squared per area for base model, b2

  //Working variables for histogrammed beam
  mutable unsigned int nbm; //!< Number of elements in beam working arrays 
  mutable unsigned int* bm_wts; //!< Histogrammed beam weights
  mutable double* inv_bm1; //!< Histogrammed inverse beam, band 1
  mutable double* inv_bm2; //!< Histogrammed inverse beam, band 2

  //Workspace
  gsl_integration_workspace *gsl_work; //!< Integration workspace for QAG

  //Convenience functions for computations without input checking
  double getSigmaInner(double) const; //!< Inner sigma computation
  double getOffsetInner(double) const; //!< Inner offset computation
  /*! \brief Inner number counts computation */
  double getNumberCountsInner(double, double) const; 

    /*! \brief Integrate powers of fluxes over number counts */
  double powerInt(double alpha, double beta) const;
  static const unsigned int nvarr; //!< Number of elements in varr
  void **varr; //!< Internal evil casting array for splineInt

  bool isValidLoaded() const; //!< See if model params are valid
 public:
  
  numberCountsDouble(const std::string&);
  ~numberCountsDouble();

  /*! \brief Get number of sources per area in base model*/
  double getBaseN0() const;
  bool isValid() const; //!< See if model params are valid

  /*! \brief Get Mean Flux per unit area, band 1 */
  double getMeanFluxPerArea1() const;
  /*! \brief Get Mean Flux per unit area, band 2 */
  double getMeanFluxPerArea2() const;
  /*! \brief Get Mean Flux squared per unit area for base model, band 1*/
  double getMeanFluxSqPerArea1() const;
  /*! \brief Get Mean Flux squared per unit area for base model, band 2*/
  double getMeanFluxSqPerArea2() const;

  /*! \brief Crude estimate of maximum flux from model */
  std::pair<double, double> getMaxFluxEstimate() const;

  /*!\brief Get differential number counts for base model*/
  double getdNdS(double, double) const;

  /*!\brief Evaluate sigma*/
  double getSigma(double) const;

  /*!\brief Evaluate offset*/
  double getOffset(double) const;

  /*! \brief Get number of knots in 1D model, sigma, and offset*/
  unsigned int getNTotalKnots() const { return nknots + nsigma + noffset; }
  /*! \brief Get number of knots */
  unsigned int getNKnots() const { return nknots; }
  /*! \brief Get Knot position */
  double getKnotPosition(unsigned int idx) const { return knotpos[idx]; }
  /*! \brief Get number of sigma knots */
  unsigned int getNSigmaKnots() const { return nsigma; }
  /*! \brief Get Sigma Knot position */
  double getSigmaKnotPosition(unsigned int idx) const { return sigmapos[idx]; }
  /*! \brief Get number of offset knots */
  unsigned int getNOffsetKnots() const { return noffset; }
  /*! \brief Get Offset Knot position */
  double getOffsetKnotPosition(unsigned int idx) const { return offsetpos[idx]; }
  /*! \brief Get Log10 differential number counts knot value */
  double getLog10KnotValue(unsigned int idx) const { 
    return logknotvals[idx] * pofd_coverage::ilogfac; }
  /*! \brief Get sigma knot value */
  double getSigmaKnotValue(unsigned int idx) const { return sigmavals[idx]; }
  /*! \brief Get offset knot value */
  double getOffsetKnotValue(unsigned int idx) const { return offsetvals[idx]; }

  /*! \brief Get number of source responses, single value version */
  double getR(double, double, const doublebeam&, double, double,
	      unsigned int) const;
  
  /*! \brief Get number of source responses, array version*/
  void getR(unsigned int, const double*, unsigned int, const double*, 
	    const doublebeam&, double, double, unsigned int,
	    double*) const;

  /*! \brief Generate a source flux from model */
  std::pair<double, double> genSource(double udev, double gdev) const;

  bool writeToStream(std::ostream& os) const; //!< Output to stream
};

/*! \brief Write to stream */
std::ostream& operator<<(std::ostream& os, const numberCountsDouble& b);

#endif
