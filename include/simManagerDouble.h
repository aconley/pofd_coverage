//simManagerDouble

#ifndef __simManagerDouble__
#define __simManagerDouble__

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_math.h>
#include<gsl/gsl_min.h>

#include "simImageDouble.h"
#include "doublebeam.h"
#include "PDDouble.h"
#include "PDFactoryDouble.h"
#include "numberCountsDouble.h"

/*!
  \brief Class for managing simulations of double maps and mapping
  out their likelihoods as a function of n0

  This works in two steps.  First, the maximum likelihood is found.
  Then (optionally) the likelihood is mapped out around that peak.
*/
class simManagerDouble {

 private:

  unsigned int nsims; //!< Number of simulations to do
  double n0initrange; //!< Initial range for likelihood peak finding step

  //Stuff for likelihood map
  bool do_map_like; //!< Create the likelihood map?
  unsigned int nlike; //!< Number of likelihoods to compute
  double n0rangefrac; //!< Fractional range in n0 to cover

  unsigned int fftsize; //!< Size of FFT (in each dim)

  //Holds result of likelihood maximization
  double *bestn0; //!< N0 that maximized the likelihood
  double *bestlike; //!< 1D array holding best likelihood value
  
  //Holds results of likelihood map

  double **likearr; //!< 2D Array of computed log likelihoods (nsims x nlike)
  double *min_n0; //!< Minimum n0 for each sim (nsims elements)
  double *delta_n0; //!< Step in n0 for each sim (nsims elements)

  //Model params
  double n0; //!< Number of sources per sq deg in input model
  double sig_i1; //!< Instrument noise, band 1
  double sig_i2; //!< Instrument noise, band 2
  double sig_i1_sm; //!< Smoothed instrument noise, band 1
  double sig_i2_sm; //!< Smoothed instrument noise, band 2

  //Stuff for doing individual sims
  double fwhm1, fwhm2; //!< Beam sizes
  double pixsize; //!< Pixel size of image and beam
  doublebeam bm; //!< Beam
  mutable simImageDouble simim; //!< Simulated image
  bool use_binning; //!< Work on binned images in likelihood
  mutable PDDouble pd; //!< Holds computed P(D)
  mutable PDFactoryDouble pdfac; //!< Computes P(D)
  mutable numberCountsDouble model; //!< Model variable

  bool do_extra_smooth; //!< Apply additional smoothing?
  double esmooth1, esmooth2; //!< Additional smoothing amount

  //Stuff for GSL minimization call
  void **varr; //!< Internal evil casting array for minimization
  gsl_min_fminimizer *s; //!< Function minimizer
  gsl_function F; //!< Interface to minimizer

#ifdef TIMING
  //Timing variables
  std::clock_t initTime, getTime, getLikeTime;
#endif

 public:
  simManagerDouble(const std::string& MODELFILE,
		   unsigned int NSIMS=200, double N0INITRANGE=0.3,
		   bool MAPLIKE=true, unsigned int NLIKE=500, 
		   double N0RANGEFRAC=0.1, unsigned int FFTSIZE=4096, 
		   unsigned int N1=720, unsigned int N2=720, 
		   double PIXSIZE=5, double FWHM1=15, double FWHM2=20, 
		   double SIGI1=0.004, double SIGI2=0.006, double N0=2.63e3, 
		   double ESMOOTH1=0, double ESMOOTH2=0, 
		   unsigned int OVERSAMPLE=1,
		   bool USEBIN=false, unsigned int NBINS=1000);
  ~simManagerDouble();

  void setSeed(unsigned long long int seed) { simim.setSeed(seed); }

  void addWisdom(const std::string& filename) { pdfac.addWisdom(filename); }
  void setNEdge(unsigned int edglen) { pdfac.setNEdge(edglen); }

  unsigned int getN1() const { return simim.getN1(); }
  unsigned int getN2() const { return simim.getN2(); }
  double getArea() const { return simim.getArea(); }
  double getFWHM1() const { return bm.getFWHM1(); }
  double getFWHM2() const { return bm.getFWHM2(); }

#ifdef TIMING
  void resetTime();
  void summarizeTime() const;
#endif

  void doSims(bool verbose); //!< Run a set of simulations

  int writeToFits(const std::string& file) const; //!< Write out results

};

#endif
