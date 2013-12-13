#include<cstring>
#include<cmath>
#include<sstream>

#include<gsl/gsl_errno.h>
#include<fitsio.h>

#include "../include/simManager.h"
#include "../include/global_settings.h"
#include "../include/pofdExcept.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

const unsigned int simManager::nbeambins = 80; 
const double simManager::nfwhm = 3.0; 

//This is the function we call to find the best fitting n0
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

  if (x > pdfac->getMaxN0())
    throw pofdExcept("", "minfunc", "N0 out of prepared range", 1);

  pdfac->getPD(x, *pd, true, true);
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
  \param[in] FILTSCALE Filtering scale, in arcsec.  If 0, no filtering is
              applied.
  \param[in] SIGI Instrument noise (without smoothing) in Jy
  \param[in] N0 Simulated number of sources per sq deg.
  \param[in] ESMOOTH Amount of extra smoothing to apply
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
		       double PIXSIZE, double FWHM, double FILTSCALE,
		       double SIGI, double N0, double ESMOOTH,
		       unsigned int OVERSAMPLE, 
		       const std::string& POWERSPECFILE,
		       unsigned int SPARCITY,
		       bool USEBIN, unsigned int NBINS) :
  nsims(NSIMS), n0initrange(N0INITRANGE), do_map_like(MAPLIKE),
  nlike(NLIKE), n0rangefrac(N0RANGEFRAC), like_sparcity(SPARCITY),
  fftsize(FFTSIZE), n0(N0), sig_i(SIGI), sig_i_sm(SIGI), fwhm(FWHM), 
  pixsize(PIXSIZE), inv_bmhist(nbeambins), filtscale(FILTSCALE), filt(NULL),
  simim(N1, N2, PIXSIZE, FWHM, SIGI, ESMOOTH, OVERSAMPLE, 
	NBINS, POWERSPECFILE), 
  oversample(OVERSAMPLE), use_binning(USEBIN), model(MODELFILE), 
  esmooth(ESMOOTH) {

#ifdef TIMING
  initTime = getTime = getLikeTime = 0;
#endif

  if (ESMOOTH > 0.0) {
    do_extra_smooth = true;
    bm.setFWHM(std::sqrt(fwhm*fwhm+esmooth*esmooth));
    sig_i_sm = simim.getSmoothedNoiseEstimate();
  } else {
    do_extra_smooth = false;
    bm.setFWHM(fwhm);
  }

  if (nsims > 0) {
    bestn0 = new double[nsims];
    bestlike = new double[nsims];
  } else {
    bestn0 = NULL;
    bestlike = NULL;
  }
  
  if (do_map_like && (nsims > 0)) {
    likearr = new double*[nsims];
    for (unsigned int i = 0; i < nsims; ++i) likearr[i] = new double[nlike];
    min_n0 = new double[nsims];
    delta_n0 = new double[nsims];
  } else {
    likearr = NULL;
    min_n0 = NULL;
    delta_n0 = NULL;
  }

  // Set up the histogrammed beam
  if (filtscale > 0)
    filt = new hipassFilter(filtscale / PIXSIZE);
  inv_bmhist.fill(bm, nfwhm, PIXSIZE, filt, filtscale, true, oversample);

  varr = new void*[4];
  s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

}

simManager::~simManager() {
  if (bestn0 != NULL) delete[] bestn0;
  if (bestlike != NULL) delete[] bestlike;
  delete[] varr;
  gsl_min_fminimizer_free(s);
  if (likearr != NULL) {
    for (unsigned int i = 0; i < nsims; ++i) delete[] likearr[i];
    delete[] likearr;
  }
  if (min_n0 != NULL) delete[] min_n0;
  if (delta_n0 != NULL) delete[] delta_n0;
  if (filt != NULL) delete filt;
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
  double sigval, maxflux; // Params to pass to minimizer
  
  //Turn off gsl error handler during call
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

  if (do_extra_smooth) sigval = sig_i_sm; else sigval = sig_i;    

  //Now, set up the parameters to pass to the minimizer.  Ugly!
  varr[0] = static_cast<void*>(&pdfac);
  varr[1] = static_cast<void*>(&pd);
  varr[2] = static_cast<void*>(&simim);
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
#ifdef TIMING
  starttime = std::clock();
#endif
  maxflux = max_n0ratio * model.getMaxKnotPosition();
  pdfac.initPD(fftsize, sigval, maxflux, init_b, model,
	       inv_bmhist);
#ifdef TIMING
  initTime += std::clock() - starttime;
#endif 

  //Main loop
  // 1) Make a simulated image
  // 2) Find the best fitting n0
  // 3) Optionally, map out the likelihood around that
  for (unsigned int i = 0; i < nsims; ++i) {

    if (verbose) {
      std::cout << "Doing simulation " << i + 1 << " of " << nsims 
		<< std::endl;
    }

    //Make simulated image (mean subtracted, possibly with filtering)
    simim.realize(model, n0, filt, do_extra_smooth, true, use_binning,
		  like_sparcity);

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
    } while ( (exp_iter < max_expiter) && ( (fa <= fm) || (fb <= fm) ) );
    
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
	pdfac.getPD(curr_n0, pd, true, true);
#ifdef TIMING
	getTime += std::clock() - starttime;
#endif

	//Get like
#ifdef TIMING
	starttime = std::clock();
#endif
	likearr[i][j] =  pd.getLogLike(simim, like_sparcity);
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
    fits_create_tbl(fp, BINARY_TBL, 0, 5, ttype, tform, NULL, "RESULTS", 
		    &status );
    delete[] tform2;
  } else {
    char* ttype[] = {"BEST_N0","BEST_LOGLIKE"};
    char* tform[2];
    char* tform0 = "1D";
    tform[0] = tform0; tform[1] = tform0;
    fits_create_tbl(fp, BINARY_TBL, 0, 2, ttype, tform, NULL, 
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
  dtmp = fwhm;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM"), &dtmp, 
		 const_cast<char*>("Beam fwhm [arcsec]"), 
		 &status);
  dtmp = simim.getBeamSum();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BMAREA"), &dtmp, 
		 const_cast<char*>("Beam area [pixels]"), 
		 &status);
  dtmp = simim.getBeamSumSq();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BMAREASQ"), &dtmp, 
		 const_cast<char*>("Beam squared area [pixels]"), 
		 &status);


  if (do_extra_smooth) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &itmp,
		   const_cast<char*>("Additional smoothing applied"), 
		   &status);
    dtmp = esmooth;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("ESMOOTH"), &dtmp, 
		   const_cast<char*>("Extra smoothing fwhm [arcsec]"), 
		   &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &itmp,
		   const_cast<char*>("Additional smoothing applied"), 
		   &status);
  }
  dtmp = sig_i;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGI"), &dtmp, 
		 const_cast<char*>("Instrument noise"), 
		 &status);
  if (do_extra_smooth) {
    dtmp = sig_i_sm;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGISM"), &dtmp, 
		   const_cast<char*>("Smoothed instrument noise"), 
		   &status);
  }

  if (oversample > 1) {
    utmp = oversample;
    fits_write_key(fp, TUINT, const_cast<char*>("OVERSMPL"), &utmp, 
		   const_cast<char*>("Oversampling factor"), 
		   &status);
  }

  itmp = static_cast<int>(simim.isClustered());
  fits_write_key(fp, TLOGICAL, const_cast<char*>("CLUSTPOS"), &itmp,
		 const_cast<char*>("Use clustered positions"), &status);

  if (use_binning) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("USEBIN"), &itmp,
		   const_cast<char*>("Use binned likelihood"), 
		   &status);
    utmp = simim.getNBins();
    fits_write_key(fp, TUINT, const_cast<char*>("NBINS"), &utmp,
		   const_cast<char*>("Number of bins in Likelihood"), 
		   &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("USEBIN"), &itmp,
		   const_cast<char*>("Use binned likelihood"), 
		   &status);
  }

  dtmp = simim.getPixSize();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("PIXSIZE"), &dtmp, 
		 const_cast<char*>("Simulation pixel size [arcsec]"), 
		 &status);
  dtmp = simim.getArea();
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
  fits_create_tbl(fp, BINARY_TBL, 0, 2, mttype, mtform, NULL, 
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
