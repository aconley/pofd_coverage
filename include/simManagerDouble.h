//simManagerDouble

#ifndef __simManagerDouble__
#define __simManagerDouble__

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_math.h>
#include<gsl/gsl_min.h>

#include "../include/global_settings.h"
#include "../include/simImageDouble.h"
#include "../include/doublebeam.h"
#include "../include/PDDouble.h"
#include "../include/PDFactoryDouble.h"
#include "../include/numberCountsDouble.h"

/*!
  \brief Class for managing simulations of double maps and mapping
  out their likelihoods as a function of n0

  This works in two steps.  First, the maximum likelihood is found.
  Then (optionally) the likelihood is mapped out around that peak.
*/
class simManagerDouble {

 private:

  static const unsigned int nnoisetrials; //!< Number of trials to carry out when computing filtered noise level

  unsigned int nsims; //!< Number of simulations to do
  double n0initrange; //!< Initial range for likelihood peak finding step

  //Stuff for likelihood map
  bool do_map_like; //!< Create the likelihood map?
  unsigned int nlike; //!< Number of likelihoods to compute
  double n0rangefrac; //!< Fractional range in n0 to cover
  unsigned int like_sparcity; //!< Sparcity of likelihood computation

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

  //Stuff for doing individual sims
  doublebeam bm; //!< Beam
  doublebeamHist inv_bmhist; //!< Histogrammed inverse beam
  simImageDouble *simim; //!< Simulated image
  bool use_binning; //!< Work on binned images in likelihood
  mutable PDDouble pd; //!< Holds computed P(D)
  mutable PDFactoryDouble pdfac; //!< Computes P(D)
  mutable numberCountsDouble model; //!< Model variable

  // Filtering.  If only filt1 is set, use that for both.
  fourierFilter *filt1; //!< Band 1 (and maybe 2) fourier space filter
  fourierFilter *filt2; //!< Band 2 fourier space filter

  double esmooth1, esmooth2; //!< Additional smoothing amount

  dblpair sigi_final; //!< Final noise in each band.

  //Stuff for GSL minimization call
  void **varr; //!< Internal evil casting array for minimization
  gsl_min_fminimizer *s; //!< Function minimizer
  gsl_function F; //!< Interface to minimizer

#ifdef TIMING
  //Timing variables
  mutable std::clock_t initTime, getTime, getLikeTime, starttime;
#endif

 public:
  explicit simManagerDouble(const std::string& MODELFILE,
			    unsigned int NSIMS=200, double N0INITRANGE=0.3,
			    bool MAPLIKE=true, unsigned int NLIKE=401,
			    double N0RANGEFRAC=0.1, unsigned int FFTSIZE=4096,
			    unsigned int N1=720, unsigned int N2=720,
			    double PIXSIZE=5.0, double FWHM1=15,
			    double FWHM2=20, double NFWHM=15.0,
			    double SIMFWHM1=pofd_coverage::qnan,
			    double SIMFWHM2=pofd_coverage::qnan,
			    double SIGI1=0.004, double SIGI2=0.006,
			    bool SINGLEFILT=true, double FILTSCALE=0.0,
			    double QFACTOR=0.2, bool MATCHED=false,
			    double FILTFWHM=0.0, double SIGMI1=0.0,
			    double SIGMI2=0.0, double SIGC1=0.006,
			    double SIGC2=0.006, unsigned int NBEAMBINS=150,
			    double N0=2.63e3, double ESMOOTH1=0,
			    double ESMOOTH2=0, unsigned int OVERSAMPLE=1,
			    const std::string& POWERSPECFILE="",
			    unsigned int SPARCITY=1, bool USEBIN=false,
			    unsigned int NBINS=1000);
  ~simManagerDouble();

  void setSeed(unsigned long long int seed) { simim->setSeed(seed); }

  void addWisdom(const std::string& filename) { pdfac.addWisdom(filename); }
  void setNEdge(unsigned int edglen) { pdfac.setNEdge(edglen); }

  unsigned int getN1() const { return simim->getN1(); }
  unsigned int getN2() const { return simim->getN2(); }
  double getArea() const { return simim->getArea(); }
  std::pair<double, double> getFWHM() const { return bm.getFWHM(); }

#ifdef TIMING
  void resetTime();
  void summarizeTime() const;
#endif

  void doSims(bool verbose); //!< Run a set of simulations

  int writeToFits(const std::string& file) const; //!< Write out results
  void writeToHDF5(const std::string& file) const; //!< Write to HDF5
};

#endif
