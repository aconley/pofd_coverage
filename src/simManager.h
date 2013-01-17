//simManager

#ifndef __simManager__
#define __simManager__

#ifdef TIMING
#include<ctime>
#endif

#include<string>

#include<gsl/gsl_math.h>
#include<gsl/gsl_min.h>

#include<simImage.h>
#include<beam.h>
#include<PD.h>
#include<PDFactory.h>
#include<numberCounts.h>

/*!
  \brief Class for managing simulations of single band maps and mapping
  out their likelihoods as a function of n0

  This works in two steps.  First, the maximum likelihood is found.
  Then (optionally) the likelihood is mapped out around that peak.
*/
class simManager {

 private:

  unsigned int nsims; //!< Number of simulations to do
  double n0initrange; //!< Initial range for likelihood peak finding step
  
  //Stuff for likelihood map
  bool do_map_like; //!< Create the likelihood map?
  unsigned int nlike; //!< Number of likelihoods to compute
  double n0rangefrac; //!< Fractional range in n0 to cover

  unsigned int fftsize; //!< Size of FFT

  //Holds result of likelihood maximization
  double *bestn0; //!< N0 that maximized the likelihood
  double *bestlike; //!< 1D array holding best likelihood value
  
  //Holds results of likelihood map
  double **likearr; //!< 2D Array of computed log likelihoods (nsims x nlike)
  double *min_n0; //!< Minimum n0 for each sim (nsims elements)
  double *delta_n0; //!< Step in n0 for each sim (nsims elements)

  //Model params
  double n0; //!< Number of sources per sq deg in input model
  std::string modelfile; //!< File we got model from
  double sig_i; //!< Instrument noise
  double sig_i_sm; //!< Smoothed instrument noise; is sig_i if no smoothing

  //Stuff for doing individual sims
  double fwhm; //!< FWHM of image beam
  double pixsize; //!< Pixel size of image and beam
  beam bm; //!< Beam
  mutable simImage simim; //!< Simulated image
  unsigned int oversample; //!< Amount of oversampling when generating image
  bool use_binning; //!< Work on binned images in likelihood
  mutable PD pd; //!< Holds computed P(D)
  mutable PDFactory pdfac; //!< Computes P(D)
  mutable numberCounts model; //!< Model variable

  //Additional smoothing
  bool do_extra_smooth; //!< Extra smoothing step
  double esmooth; //!< Additional smoothing

  //Stuff for GSL minimization call
  void **varr; //!< Internal evil casting array for minimization
  gsl_min_fminimizer *s; //!< Function minimizer
  gsl_function F; //!< Interface to minimizer

#ifdef TIMING
  //Timing variables
  std::clock_t initTime, getTime, getLikeTime;
#endif

 public:
  simManager(const std::string& MODELFILE,
	     unsigned int NSIMS=1000, double N0INITRANGE=0.3, 
	     bool MAPLIKE=true, unsigned int NLIKE=2000, 
	     double N0RANGEFRAC=0.1, unsigned int FFTSIZE=262144, 
	     unsigned int N1=720, unsigned int N2=720,
	     double PIXSIZE=5.0, double FWHM=15.0, double SIGI=0.005, 
	     double N0=2.6e3, double ESMOOTH=0.0,
	     unsigned int OVERSAMPLE=1, bool USEBIN=false,
	     unsigned int NBINS=1000);
  ~simManager();

  void setSeed(unsigned long long int seed) { simim.setSeed(seed); }

  void addWisdom(const std::string& filename) { pdfac.addWisdom(filename); }

  unsigned int getN1() const { return simim.getN1(); }
  unsigned int getN2() const { return simim.getN2(); }
  double getArea() const { return simim.getArea(); }
  double getFWHM() const { return bm.getFWHM(); }

#ifdef TIMING
  void resetTime();
  void summarizeTime() const;
#endif

  void doSims(bool verbose); //!< Run a set of simulations

  int writeToFits(const std::string& file) const; //!< Write out results

};

#endif
