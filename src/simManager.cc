#include<cstring>
#include<cmath>
#include<sstream>
#include<limits>

#include<gsl/gsl_errno.h>
#include<fitsio.h>
#include<hdf5.h>

#include "../include/simManager.h"
#include "../include/global_settings.h"
#include "../include/pofdExcept.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

const unsigned int simManager::nnoisetrials = 9;

//This is the function we call to find the best fitting n0
// It can't be in simManager because the GSL doesn't like
// being given class methods
/*!
  \param[in] x Value of n0
  \param[in] params Messy array to hold everything to compute likelihood
*/
static double minfunc(double x, void* params) {
  //Basically in the order we need them
  //This is very, very ugly, but seems to be the only way
  // to interface class members to the GSL
  //PDFactory::initPD should have been called first

  void** vptr = static_cast<void**>(params);
  PDFactory *pdfac = static_cast<PDFactory*>(vptr[0]);
  PD *pd = static_cast<PD*>(vptr[1]);
  simImage *im = static_cast<simImage*>(vptr[2]);
  unsigned int sparcity = *static_cast<unsigned int*>(vptr[3]);

  if (x > pdfac->getMaxN0()) {
    std::stringstream errmsg;
    errmsg << "N0 out of prepared range; x value: " << x << " max N0: "
	   << pdfac->getMaxN0();
    throw pofdExcept("", "minfunc", errmsg.str(), 1);
  }

  pdfac->getPD(x, *pd, true);
  double loglike = pd->getLogLike(*im, sparcity);

  return -loglike; //Remember -- we want to minimize this

}

//////////////////////////////////////////////


/*!
  \param[in] MODELFILE Name of file containing base model
  \param[in] NSIMS Number of simulations to do
  \param[in] N0INITRANGE Initial values to maximize likelihood are
              n0*(1-n0initrange) to n0*(1+n0initrange)
  \param[in] MAPLIKE Do make likelihood map?
  \param[in] NLIKE   Number of likelihoods in likelihood map
  \param[in] N0RANGEFRAC Likelihood map covers n0*(1+n0rangefrac) to
                         n0*(1-n0rangefrac)
  \param[in] FFTSIZE Number of elements to use in FFT
  \param[in] N1 Number of pixels along first axis in simulated image
  \param[in] N2 Number of pixels along second axis in simulated image
  \param[in] PIXSIZE Pixel size, in arcsec
  \param[in] FWHM Fwhm of beam, in arcsec
  \param[in] NFWHM Number of FWHM out to go on beam.  This is the number
              after filtering (if filtering is applied)
  \param[in] SIMFWHM FWHM of beam used in constructing simulated maps.
              If this is NaN (the default), the same FWHM is used
              as is used for the P(D) computations.  Allowing these to
              be different is a way of doing systematics testing.
  \param[in] SIGI Instrument noise (without smoothing or filtering) in Jy
  \param[in] FILTSCALE High-pass filtering scale, in arcsec.  If 0, no
              high-pass filtering is applied.
  \param[in] QFACTOR High-pass Gaussian filtering apodization sigma as
              a fraction of FILTSCALE.
  \param[in] MATCHED Apply matched filtering using the FWHM of the beam,
               the instrument noise (SIGI), and SIGC
  \param[in] SIGMI The matched filter instrument noise, if matched filtering
              is used.  If zero, defaults to SIGI
  \param[in] SIGMC The matched filter confusion noise, if matched filtering
              is used
  \param[in] NBEAMBINS Number of bins to use in beam histogram; def 100
  \param[in] N0 Simulated number of sources per sq deg.
  \param[in] ESMOOTH Amount of extra Gaussian smoothing to apply
  \param[in] OVERSAMPLE Amount of oversampling to use when simulating image
  \param[in] POWERSPECFILE File containing power spectrum for on-sky source
              distribution.  If not provided, uniform sampling is used.
  \param[in] SPARCITY Sparcity of data sampling for likelihood computation.
             0 or 1 means fully sample, using all data.
  \param[in] USEBIN Bin the data in the likelihood calculation
  \param[in] NBINS Number of bins in binned data
 */
simManager::simManager(const std::string& MODELFILE,
		       unsigned int NSIMS, double N0INITRANGE,
		       bool MAPLIKE, unsigned int NLIKE,
		       double N0RANGEFRAC, unsigned int FFTSIZE,
		       unsigned int N1, unsigned int N2,
		       double PIXSIZE, double FWHM, double NFWHM,
		       double SIMFWHM, double SIGI, double FILTSCALE,
           double QFACTOR, bool MATCHED, double SIGMI, double SIGMC,
		       unsigned int NBEAMBINS, double N0,
		       double ESMOOTH, unsigned int OVERSAMPLE,
		       const std::string& POWERSPECFILE,
		       unsigned int SPARCITY,
		       bool USEBIN, unsigned int NBINS) :
  nsims(NSIMS), n0initrange(N0INITRANGE), do_map_like(MAPLIKE),
  nlike(NLIKE), n0rangefrac(N0RANGEFRAC), like_sparcity(SPARCITY),
  fftsize(FFTSIZE), n0(N0), inv_bmhist(NBEAMBINS), simim(nullptr),
  use_binning(USEBIN), model(MODELFILE), filt(nullptr), esmooth(ESMOOTH) {

#ifdef TIMING
  initTime = getTime = getLikeTime = 0;
#endif

  // Setup simim
  if (std::isnan(SIMFWHM) || SIMFWHM <= 0.0) SIMFWHM = FWHM;
  simim = new simImage(N1, N2, PIXSIZE, SIMFWHM, SIGI, ESMOOTH,
    OVERSAMPLE, NBINS, POWERSPECFILE);

  if (esmooth > 0)
    bm.setFWHM(std::sqrt(FWHM * FWHM + esmooth * esmooth));
  else
    bm.setFWHM(FWHM);

  if (nsims > 0) {
    bestn0 = new double[nsims];
    bestlike = new double[nsims];
  } else {
    bestn0 = nullptr;
    bestlike = nullptr;
  }

  if (do_map_like && (nsims > 0)) {
    likearr = new double*[nsims];
    for (unsigned int i = 0; i < nsims; ++i) likearr[i] = new double[nlike];
    min_n0 = new double[nsims];
    delta_n0 = new double[nsims];
  } else {
    likearr = nullptr;
    min_n0 = nullptr;
    delta_n0 = nullptr;
  }

  // Set up filter if needed
  if (MATCHED && SIGMI == 0) SIGMI = SIGI;
  if (FILTSCALE > 0) {
    if (MATCHED) // Hipass and matched
      filt = new fourierFilter(PIXSIZE, FWHM, SIGMI, SIGMC, FILTSCALE,
			       QFACTOR, false, true);
    else // hipass only
      filt = new fourierFilter(PIXSIZE, FILTSCALE, QFACTOR, false, true);
  } else if (MATCHED) // Matched only
    filt = new fourierFilter(PIXSIZE, FWHM, SIGMI, SIGMC, false, true);

  // Set up the histogrammed beam
  //  If we are filtering take a large piece
  if (filt != nullptr)
    inv_bmhist.fill(bm, N1, N2, PIXSIZE, true, OVERSAMPLE, filt, NFWHM);
  else
    inv_bmhist.fill(bm, NFWHM, PIXSIZE, true, OVERSAMPLE, filt, NFWHM);

  // Setup GSL interface
  varr = new void*[4];
  s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

  sigi_final = pofd_coverage::qnan;
}

simManager::~simManager() {
  if (simim != nullptr) delete simim;
  if (bestn0 != nullptr) delete[] bestn0;
  if (bestlike != nullptr) delete[] bestlike;
  delete[] varr;
  gsl_min_fminimizer_free(s);
  if (likearr != nullptr) {
    for (unsigned int i = 0; i < nsims; ++i) delete[] likearr[i];
    delete[] likearr;
  }
  if (min_n0 != nullptr) delete[] min_n0;
  if (delta_n0 != nullptr) delete[] delta_n0;
  if (filt != nullptr) delete filt;
}

/*!
  \param[in] verbose Output status messages as it runs; basically, what
                     simulation are we on.
*/
void simManager::doSims(bool verbose=false) {
  //Two steps; first: find the best fit n0.  Then -- maybe --
  // map out the likelihood
  const unsigned int max_iter = 100; // Maximum number of iters for min
  const unsigned int max_expiter = 10; //!< Maximum number of bracket steps
  const double reltol = 0.001; // Relative tolerance in n0: 0.1%
  double maxflux; // Params to pass to minimizer

  //Turn off gsl error handler during call
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

  // Should only be computed the first time
  sigi_final = simim->getFinalNoise(nnoisetrials, filt);

  //Now, set up the parameters to pass to the minimizer.  Ugly!
  varr[0] = static_cast<void*>(&pdfac);
  varr[1] = static_cast<void*>(&pd);
  varr[2] = static_cast<void*>(simim);
  varr[3] = static_cast<void*>(&like_sparcity);

  void *params;
  params = static_cast<void*>(varr);
  F.function = &minfunc;
  F.params = params;

  //First guess at range
  double init_b = n0 * (1.0 + n0initrange);
  double max_n0ratio = 1.0 + n0initrange;

  //We are going to reuse R for all compuations of all models.
  // That is, we call pdfactory.initPD once and re-use it forever.
  // But that means making an informed guess about the maximum flux
  // to ask for.
  if (verbose)
    std::cout << "Initializing P(D)" << std::endl;
#ifdef TIMING
  starttime = std::clock();
#endif
  maxflux = 1.05 * max_n0ratio * model.getMaxKnotPosition();
  pdfac.initPD(fftsize, sigi_final, maxflux, init_b, model,
	       inv_bmhist);
#ifdef TIMING
  initTime += std::clock() - starttime;
#endif

  //Main loop
  // 1) Make a simulated image
  // 2) Find the best fitting n0
  // 3) Optionally, map out the likelihood around that
  for (unsigned int i = 0; i < nsims; ++i) {

    if (verbose)
      std::cout << "Doing simulation " << i + 1 << " of " << nsims
		<< std::endl;

    //Make simulated image (mean subtracted, possibly with filtering)
    simim->realize(model, n0, true, filt, use_binning, like_sparcity);

    //Now set up the minimization; note this involves calling minfunc
    // so we can't do it until all the arguments are ready
    //We do this explicitly ourselves so that we can iterate outwards
    // on init_a, init_b as needed
    double curr_initrange, next_initrange;
    double a, b, fa, fm, fb;
    fm = minfunc(n0, params);

    //We move in logarithmically with a maximum number of steps
    unsigned int exp_iter = 0;
    next_initrange = n0initrange;
    do {
      ++exp_iter;
      curr_initrange = next_initrange;
      if (curr_initrange >= 1.0) //Shouldn't be possible, but...
	throw pofdExcept("simManager", "doSims",
			 "Logic error generating bracketing range", 1);

      //Try with current init range
      a = n0 * (1.0 - curr_initrange);
      b = n0 * (1.0 + curr_initrange);
      fa = minfunc(a, params);
      fb = minfunc(b, params);

      //Now, update current init range.  We step a quarter way in log
      // space between here and 1.0.
      next_initrange = exp2(0.75 * log2(curr_initrange));

    } while ((exp_iter < max_expiter) && ((fa <= fm) || (fb <= fm)));

    if ((fa <= fm) || (fb <= fm)) {
      std::stringstream errstr;
      errstr << "Unable to bracket minimum in -log likelihood."
	     << std::endl;
      errstr << "Final values:" << std::endl;
      errstr << "\tf(" << a << ") = " << fa << std::endl;
      errstr << "\tf(" << n0 << ") = " << fm << std::endl;
      errstr << "\tf(" << b << ") = " << fb << std::endl;
      errstr << "n0 range: " << curr_initrange << std::endl;
      throw pofdExcept("simManager", "doSims", errstr.str(), 2);
    }

    int status;
    status = gsl_min_fminimizer_set_with_values(s, &F, n0, fm,
						a, fa, b, fb);
    if (status != GSL_SUCCESS) {
      std::stringstream errstr;
      errstr << "GSL minimizer setup failed with code: " << status
	     << std::endl;
      errstr << "GSL error message: " << gsl_strerror(status);
      throw pofdExcept("simManager", "doSims", errstr.str(), 3);
    }

    //Now that we have bracketed, do the actual minimization
    unsigned int iter;
    double currmax, currmin; //Current limits
    iter = 0;
    do {
      ++iter;

      //Do a minimization step
      status = gsl_min_fminimizer_iterate(s);
      if (status != GSL_SUCCESS) break;

      //Get current limits
      currmin = gsl_min_fminimizer_x_lower(s);
      currmax = gsl_min_fminimizer_x_upper(s);

      //Convergence test
      status = gsl_min_test_interval(currmin, currmax, 0.0, reltol);
    } while (status == GSL_CONTINUE && iter < max_iter);

    if (status != GSL_SUCCESS) {
      std::stringstream errstr;
      errstr << "GSL minimizer failed to converge with code: " << status
	     << std::endl;
      errstr << "GSL error message: " << gsl_strerror(status);
      throw pofdExcept("simManager", "doSims",
		       errstr.str(), 4);
    }

    bestn0[i] = gsl_min_fminimizer_x_minimum(s);
    bestlike[i] = -gsl_min_fminimizer_f_minimum(s); //From -Log like to log like

    //Now -- optional likelihood map
    if (do_map_like) {
      //Compute n0 array
      if (nlike > 1) {
	double max_n0;
	min_n0[i] = bestn0[i] * (1.0 - n0rangefrac);
	max_n0 = bestn0[i] * (1.0 + n0rangefrac);
	delta_n0[i] = (max_n0 - min_n0[i]) /
	  (static_cast<double>(nlike) - 1.0);
      } else {
	min_n0[i] = n0;
	delta_n0[i] = 0.0;
      }

      //Compute P(D) and like
#ifdef TIMING
      std::clock_t starttime;
#endif
      double curr_n0;
      for (unsigned int j = 0; j < nlike; ++j) {
	curr_n0 = min_n0[i] + static_cast<double>(j) * delta_n0[i];

#ifdef TIMING
	starttime = std::clock();
#endif
	pdfac.getPD(curr_n0, pd, true);
#ifdef TIMING
	getTime += std::clock() - starttime;
#endif

	//Get like
#ifdef TIMING
	starttime = std::clock();
#endif
	likearr[i][j] =  pd.getLogLike(*simim, like_sparcity);
#ifdef TIMING
	getLikeTime += std::clock() - starttime;
#endif
      }
    }
  }
#ifdef TIMING
  summarizeTime();
  resetTime();
#endif

  //Restore the error handler
  gsl_set_error_handler(old_handler);

}

#ifdef TIMING
void simManager::resetTime() {
  initTime = getTime = getLikeTime = 0;
  pdfac.resetTime();
}

void simManager::summarizeTime() const {
  std::cout << "initPD time: " << 1.0*initTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "getPD time: " << 1.0*getTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "getLogLike time: " << 1.0*getLikeTime/CLOCKS_PER_SEC
	    << std::endl;
  std::cout << "Within PDFactory: " << std::endl;
  pdfac.summarizeTime(2);
}
#endif

//Output
int simManager::writeToFits(const std::string& outputfile) const {
  //Make the fits file
  int status = 0;
  fitsfile *fp;

  fits_create_file(&fp, outputfile.c_str(), &status);
  if (status) {
    fits_report_error(stderr,status);
    return status;
  }

  //We are going to do this as a big ol table
  //The contents of the table depend if we mapped out the likelihood or
  // just maximized
  //A tricky/irritating bit is constructing the tform argument
  // for the LOGLIKE column
  if (do_map_like) {
    char* ttype[] = {"BEST_N0", "BEST_LOGLIKE", "MIN_N0",
		     "DELTA_N0", "LOGLIKE"};
    char* tform[5];
    char* tform0 = "1D";
    tform[0] = tform0; tform[1] = tform0; tform[2] = tform0; tform[3] = tform0;
    int ndigits = static_cast<int>(log10(nlike) + 1.99999999999);
    char* tform2 = new char[ndigits + 5];
    sprintf(tform2, "%uE", nlike); //E is the format code for 32 bit float
    tform[4] = tform2;
    fits_create_tbl(fp, BINARY_TBL, 0, 5, ttype, tform, nullptr, "RESULTS",
		    &status );
    delete[] tform2;
  } else {
    char* ttype[] = {"BEST_N0","BEST_LOGLIKE"};
    char* tform[2];
    char* tform0 = "1D";
    tform[0] = tform0; tform[1] = tform0;
    fits_create_tbl(fp, BINARY_TBL, 0, 2, ttype, tform, nullptr,
		    "RESULTS", &status);
  }

  //Write header stuff
  unsigned int utmp = 0;
  int itmp = 0;
  double dtmp = 0.0;
  fits_write_key(fp, TSTRING, const_cast<char*>("MODEL"),
		 const_cast<char*>("Spline"),
		 const_cast<char*>("Model type"),
		 &status );
  dtmp = n0; //Must copy to temporary for const type handling
  fits_write_key(fp, TDOUBLE, const_cast<char*>("N0"), &dtmp,
		 const_cast<char*>("Number of sources per sq deg"),
		 &status);

  //Base model parameters
  dtmp = model.getBaseN0();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BASEN0"), &dtmp,
		 const_cast<char*>("Base number of sources per sq deg"),
		 &status);
  utmp = model.getNKnots();
  fits_write_key(fp, TUINT, const_cast<char*>("NKNOTS"), &utmp,
		 const_cast<char*>("Number of model knots"),
		 &status);

  //Sim params
  dtmp = inv_bmhist.getFWHM();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM"), &dtmp,
		 const_cast<char*>("Beam fwhm [arcsec]"),
		 &status);
  dtmp = inv_bmhist.getNFWHMKeep();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("NFWHM"), &dtmp,
		 const_cast<char*>("Number of FWHM kept"),
		 &status);
  utmp = inv_bmhist.getNbins();
  fits_write_key(fp, TUINT, const_cast<char*>("NBMBINS"), &utmp,
		 const_cast<char*>("Number of Beam hist bins"),
		 &status);
  dtmp = simim->getBeamFWHM();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIMFWHM"), &dtmp,
     const_cast<char*>("Simulated Image Beam FWHM [arcsec]"),
     &status);
  dtmp = simim->getBeamSum();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BMAREA"), &dtmp,
		 const_cast<char*>("Simulated Beam area [pixels]"),
		 &status);
  dtmp = simim->getBeamSumSq();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BMAREASQ"), &dtmp,
		 const_cast<char*>("Simulated Beam squared area [pixels]"),
		 &status);

  if (simim->isSmoothed()) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &itmp,
		   const_cast<char*>("Additional smoothing applied"),
		   &status);
    dtmp = simim->getEsmooth();
    fits_write_key(fp, TDOUBLE, const_cast<char*>("ESMOOTH"), &dtmp,
		   const_cast<char*>("Extra smoothing fwhm [arcsec]"),
		   &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &itmp,
		   const_cast<char*>("Additional smoothing applied"),
		   &status);
  }
  dtmp = simim->getInstNoise();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGI"), &dtmp,
		 const_cast<char*>("Raw instrument noise"),
		 &status);
  dtmp = sigi_final;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGIFNL"), &dtmp,
		 const_cast<char*>("Final instrument noise"),
		 &status);
  if (simim->isOversampled()) {
    utmp = simim->getOversampling();
    fits_write_key(fp, TUINT, const_cast<char*>("OVERSMPL"), &utmp,
		   const_cast<char*>("Beam Oversampling factor"),
		   &status);
  }

  itmp = static_cast<int>(simim->isClustered());
  fits_write_key(fp, TLOGICAL, const_cast<char*>("CLUSTPOS"), &itmp,
		 const_cast<char*>("Use clustered positions"), &status);

  if (use_binning) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("USEBIN"), &itmp,
		   const_cast<char*>("Use binned likelihood"),
		   &status);
    utmp = simim->getNBins();
    fits_write_key(fp, TUINT, const_cast<char*>("NBINS"), &utmp,
		   const_cast<char*>("Number of bins in Likelihood"),
		   &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("USEBIN"), &itmp,
		   const_cast<char*>("Use binned likelihood"),
		   &status);
  }

  dtmp = simim->getPixSize();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("PIXSIZE"), &dtmp,
		 const_cast<char*>("Simulation pixel size [arcsec]"),
		 &status);
  dtmp = simim->getArea();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("AREA"), &dtmp,
		 const_cast<char*>("Simulation area size [sq deg]"),
		 &status);

  utmp = nsims;
  fits_write_key(fp, TUINT, const_cast<char*>("NSIMS"), &utmp,
		 const_cast<char*>("Number of simulations"),
		 &status);

  utmp = fftsize;
  fits_write_key(fp, TUINT, const_cast<char*>("FFTSIZE"), &utmp,
		 const_cast<char*>("Size of FFT transformation"),
		 &status);
  dtmp = n0initrange;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("N0INIRNG"), &dtmp,
		 const_cast<char*>("N0 range fraction for minimization"),
		 &status);

  itmp = simim->isHipassFiltered();
  fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT"), &itmp,
		 const_cast<char*>("Has hipass filtering been applied?"),
		 &status);
  if (itmp) {
    dtmp = simim->getFiltScale();
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTSCL"), &dtmp,
		 const_cast<char*>("Hipass filtering scale [arcsec]"),
		 &status);
    dtmp = simim->getFiltQFactor();
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTQ"), &dtmp,
		 const_cast<char*>("Hipass filtering apodization"),
		 &status);
  }
  itmp = simim->isMatchFiltered();
  fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT"), &itmp,
		 const_cast<char*>("Has matched filtering been applied?"),
		 &status);
  if (itmp) {
    dtmp = simim->getFiltFWHM();
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTFWHM"), &dtmp,
		 const_cast<char*>("Matched filtering fwhm [arcsec]"),
		 &status);
    dtmp = simim->getFiltSigInst();
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTSIGI"), &dtmp,
		 const_cast<char*>("Matched filtering sig_i [Jy]"),
		 &status);
    dtmp = simim->getFiltSigConf();
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTSIGC"), &dtmp,
		 const_cast<char*>("Matched filtering sig_c [Jy]"),
		 &status);
  }

  utmp = getN1();
  fits_write_key(fp, TUINT, const_cast<char*>("N1"), &utmp,
		 const_cast<char*>("Image extent, dimension 1"),
		 &status);
  utmp = getN2();
  fits_write_key(fp, TUINT, const_cast<char*>("N2"), &utmp,
		 const_cast<char*>("Image extent, dimension 2"),
		 &status);

  if (like_sparcity > 1) {
    utmp = like_sparcity;
    fits_write_key(fp, TUINT, const_cast<char*>("LIKESPAR"), &utmp,
		   const_cast<char*>("Like sampling sparcity"),
		   &status);
  }

  if (do_map_like) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MAPLIKE"), &itmp,
		   const_cast<char*>("Do map out likelihood"),
		   &status);
    utmp = nlike;
    fits_write_key(fp, TUINT, const_cast<char*>("NLIKE"), &utmp,
		   const_cast<char*>("Number of likelihoods"),
		   &status);
    dtmp = n0rangefrac;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("N0RANGE"), &dtmp,
		   const_cast<char*>("N0 range fraction"),
		   &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MAPLIKE"), &itmp,
		   const_cast<char*>("Do map out likelihood"),
		   &status);
  }

  fits_write_key(fp, TSTRING, const_cast<char*>("VERSION"),
		 const_cast<char*>(pofd_coverage::version),
		 const_cast<char*>("pofd_coverage version"),
		 &status);

  fits_write_history(fp,
		     const_cast<char*>("Simulation results from pofd_coverage"),
		     &status);
  fits_write_date(fp, &status);

  //Now write out the data.  We write the actual array of likelihoods
  // (if present) as the logLikes minus bestlike so that we can use
  // floats rather than doubles
  fits_insert_rows(fp, 0, nsims, &status);
  if (do_map_like) {
    float *like_row = new float[nlike];
    double *rowptr;
    double val, blike;
    for (unsigned int i = 0; i < nsims; ++i) {
      val = bestn0[i];
      fits_write_col(fp, TDOUBLE, 1, i+1, 1, 1, &val, &status);
      blike = bestlike[i];
      fits_write_col(fp, TDOUBLE, 2, i+1, 1, 1, &blike, &status);
      val = min_n0[i];
      fits_write_col(fp, TDOUBLE, 3, i+1, 1, 1, &val, &status);
      val = delta_n0[i];
      fits_write_col(fp, TDOUBLE, 4, i+1, 1, 1, &val, &status);
      rowptr = likearr[i];
      for (unsigned int j = 0; j < nlike; ++j)
	like_row[j] = static_cast<float>(rowptr[j] - blike);
      fits_write_col(fp, TFLOAT, 5, i+1, 1, nlike, like_row, &status);
    }
    delete[] like_row;
  } else
    for (unsigned int i = 0; i < nsims; ++i) {
      double val;
      val = bestn0[i];
      fits_write_col(fp,TDOUBLE,1,i+1,1,1,&val,&status);
      val = bestlike[i];
      fits_write_col(fp,TDOUBLE,2,i+1,1,1,&val,&status);
    }

  //Add the model information to another extension
  char* mttype[] = {"KNOTPOS","LOG10KNOTVAL"};
  char* mtform[] = {"1D", "1D"};
  fits_create_tbl(fp, BINARY_TBL, 0, 2, mttype, mtform, nullptr,
		  "BASEMODEL", &status);
  //Base model parameters, write to this extension header as well
  dtmp = model.getBaseN0();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BASEN0"), &dtmp,
		 const_cast<char*>("Base number of sources per sq deg"),
		 &status);
  utmp = model.getNKnots();
  fits_write_key(fp, TUINT, const_cast<char*>("NKNOTS"), &utmp,
		 const_cast<char*>("Number of model knots"),
		 &status);
  fits_insert_rows(fp, 0, utmp, &status);
  for (unsigned int i = 0; i < utmp; ++i) {
    double val;
    val = model.getKnotPosition(i);
    fits_write_col(fp,TDOUBLE,1,i+1,1,1,&val,&status);
    val = model.getLog10KnotValue(i);
    fits_write_col(fp,TDOUBLE,2,i+1,1,1,&val,&status);
  }

  //Close up and go home
  fits_close_file(fp,&status);
  if (status) {
    fits_report_error(stderr,status);
    return status;
  }
  return status;
}

/*!
  \param[in] outputfile Name of file to write to as HDF5
*/
void simManager::writeToHDF5(const std::string& outputfile) const {

  hid_t file_id;
  file_id = H5Fcreate(outputfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);
  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw pofdExcept("simManager", "writeToHDF5",
		     "Failed to open HDF5 file to write", 1);
  }

  // Properties (e.g., metadata)
  hsize_t adims = 1;
  hid_t mems_id, att_id, dat_id;
  hbool_t bl;
  unsigned int utmp;
  double dtmp;

  // Model info (as a group)
  hid_t group_id;
  group_id = H5Gcreate(file_id, "Model", H5P_DEFAULT, H5P_DEFAULT,
		      H5P_DEFAULT);
  if (H5Iget_ref(group_id) < 0)
    throw pofdExcept("simManager", "writeToHDF5",
		     "Failed to create HDF5 model group", 2);
  model.writeToHDF5Handle(group_id);
  H5Gclose(group_id);

  // Beam group
  group_id = H5Gcreate(file_id, "Beam", H5P_DEFAULT, H5P_DEFAULT,
		       H5P_DEFAULT);
  if (H5Iget_ref(group_id) < 0)
    throw pofdExcept("simManager", "writeToHDF5",
		     "Failed to create HDF5 beam group", 3);
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);

  dtmp = inv_bmhist.getFWHM();
  att_id = H5Acreate2(group_id, "FWHM", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
  H5Aclose(att_id);

  dtmp = inv_bmhist.getNFWHMKeep();
  att_id = H5Acreate2(group_id, "NFWHMKept", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
  H5Aclose(att_id);

  utmp = inv_bmhist.getNbins();
  att_id = H5Acreate2(group_id, "NBeamBins", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &utmp);
  H5Aclose(att_id);

  dtmp = simim->getBeamFWHM();
  att_id = H5Acreate2(group_id, "SimBeamFWHM", H5T_NATIVE_DOUBLE,
          mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
  H5Aclose(att_id);

  dtmp = simim->getBeamSum();
  att_id = H5Acreate2(group_id, "BeamPixArea", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
  H5Aclose(att_id);

  dtmp = simim->getBeamSumSq();
  att_id = H5Acreate2(group_id, "BeamSqPixArea", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
  H5Aclose(att_id);

  if (simim->isOversampled()) {
    utmp = simim->getOversampling();
    att_id = H5Acreate2(group_id, "Oversampling", H5T_NATIVE_UINT,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_UINT, &utmp);
    H5Aclose(att_id);
  }

  H5Sclose(mems_id);
  H5Gclose(group_id);

  // Filter group
  group_id = H5Gcreate(file_id, "Filter", H5P_DEFAULT, H5P_DEFAULT,
		       H5P_DEFAULT);
  if (H5Iget_ref(group_id) < 0)
    throw pofdExcept("simManager", "writeToHDF5",
		     "Failed to create HDF5 filter group", 4);
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  bl = static_cast<hbool_t>(simim->isHipassFiltered());
  att_id = H5Acreate2(group_id, "IsHighPassFiltered", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (bl) {
    dtmp = simim->getFiltScale();
    att_id = H5Acreate2(group_id, "HighPassFiltScale", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
    H5Aclose(att_id);

    dtmp = simim->getFiltQFactor();
    att_id = H5Acreate2(group_id, "HighPassQFactor", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
    H5Aclose(att_id);
  }
  bl = static_cast<hbool_t>(simim->isMatchFiltered());
  att_id = H5Acreate2(group_id, "IsMatchFiltered", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (bl) {
    dtmp = simim->getFiltFWHM();
    att_id = H5Acreate2(group_id, "MatchedFiltFWHM", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
    H5Aclose(att_id);

    dtmp = simim->getFiltSigInst();
    att_id = H5Acreate2(group_id, "MatchedFiltSigInst", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
    H5Aclose(att_id);

    dtmp = simim->getFiltSigConf();
    att_id = H5Acreate2(group_id, "MatchedFiltSigConf", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
    H5Aclose(att_id);
  }
  H5Sclose(mems_id);
  H5Gclose(group_id);

  // Sim group
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);

  group_id = H5Gcreate(file_id, "Simulations", H5P_DEFAULT, H5P_DEFAULT,
		      H5P_DEFAULT);
  if (H5Iget_ref(group_id) < 0)
    throw pofdExcept("simManager", "writeToHDF5",
		     "Failed to create HDF5 simulation group", 5);

  att_id = H5Acreate2(group_id, "N0", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &n0);
  H5Aclose(att_id);


  bl = static_cast<hbool_t>(simim->isSmoothed());
  att_id = H5Acreate2(group_id, "Smoothed", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (bl) {
    dtmp = simim->getEsmooth();
    att_id = H5Acreate2(group_id, "ExtraSmoothing", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
    H5Aclose(att_id);
  }

  dtmp = simim->getInstNoise();
  att_id = H5Acreate2(group_id, "SigmaInst", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
  H5Aclose(att_id);

  dtmp = sigi_final;
  att_id = H5Acreate2(group_id, "SigmaInstFinal", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
  H5Aclose(att_id);

  if (simim->isOversampled()) {
    utmp = simim->getOversampling();
    att_id = H5Acreate2(group_id, "Oversampling", H5T_NATIVE_UINT,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_UINT, &utmp);
    H5Aclose(att_id);
  }

  // Write information about positions
  simim->writePositionGeneratorToHDF5Handle(group_id);

  bl = static_cast<hbool_t>(use_binning);
  att_id = H5Acreate2(group_id, "DataBinned", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (use_binning) {
    utmp = simim->getNBins();
    att_id = H5Acreate2(group_id, "NDataBins", H5T_NATIVE_UINT,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_UINT, &utmp);
    H5Aclose(att_id);
  }

  utmp = getN1();
  att_id = H5Acreate2(group_id, "N1", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &utmp);
  H5Aclose(att_id);

  utmp = getN2();
  att_id = H5Acreate2(group_id, "N2", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &utmp);
  H5Aclose(att_id);

  dtmp = simim->getPixSize();
  att_id = H5Acreate2(group_id, "PixelSize", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
  H5Aclose(att_id);

  dtmp = simim->getArea();
  att_id = H5Acreate2(group_id, "SimAreaSqDeg", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dtmp);
  H5Aclose(att_id);

  att_id = H5Acreate2(group_id, "NSims", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nsims);
  H5Aclose(att_id);

  att_id = H5Acreate2(group_id, "FFTSize", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &fftsize);
  H5Aclose(att_id);

  att_id = H5Acreate2(group_id, "N0InitRange", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &n0initrange);
  H5Aclose(att_id);

  if (like_sparcity > 1) {
    att_id = H5Acreate2(group_id, "LikelihoodSparcity", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_UINT, &like_sparcity);
    H5Aclose(att_id);
  }

  bl = static_cast<hbool_t>(do_map_like);
  att_id = H5Acreate2(group_id, "DoMapLike", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (do_map_like) {
    att_id = H5Acreate2(group_id, "NLike", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_UINT, &nlike);
    H5Aclose(att_id);

    att_id = H5Acreate2(group_id, "N0LikeRangeFraction", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &n0rangefrac);
    H5Aclose(att_id);
  }

  H5Sclose(mems_id); // The rest is not scalar

  // Write the data
  if (nsims > 0) {
    adims = nsims;
    mems_id = H5Screate_simple(1, &adims, nullptr);

    dat_id = H5Dcreate2(group_id, "BestN0", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
	     H5P_DEFAULT, bestn0);
    H5Dclose(dat_id);

    dat_id = H5Dcreate2(group_id, "BestLike", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
	     H5P_DEFAULT, bestlike);
    H5Dclose(dat_id);

    if (do_map_like) {
      dat_id = H5Dcreate2(group_id, "MinN0", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
	       H5P_DEFAULT, min_n0);
      H5Dclose(dat_id);
      dat_id = H5Dcreate2(group_id, "DeltaN0", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
	       H5P_DEFAULT, delta_n0);
      H5Dclose(dat_id);
    }

    H5Sclose(mems_id);

    // Multi-D.  We will write the likelihood as float after
    //  subtracting of bestlike to save space
    if (do_map_like) {
      hsize_t cdims[2] = {nsims, nlike};
      float *wrk, *wrowptr;
      wrk = new float[nsims * nlike];
      double cmin, *rowptr;
      for (unsigned int i = 0; i < nsims; ++i) {
	cmin = min_n0[i];
	rowptr = likearr[i];
	wrowptr = wrk + i * nlike;
	for (unsigned int j = 0; j < nlike; ++j)
	  wrowptr[j] = static_cast<float>(rowptr[j] - cmin);
      }
      // We further write with compression
      hid_t plist;
      mems_id = H5Screate_simple(2, cdims, nullptr);
      plist = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_chunk(plist, 2, cdims);
      H5Pset_deflate(plist, 6); // Gzip, compression 6
      dat_id = H5Dcreate2(group_id, "LogLikelihood", H5T_NATIVE_FLOAT,
			  mems_id, H5P_DEFAULT, plist, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, wrk);
      H5Dclose(dat_id);
      H5Pclose(plist);
      H5Sclose(mems_id);
      delete[] wrk;
    }
  }

  H5Gclose(group_id);  // Simulation group

  H5Fclose(file_id); // File
}
