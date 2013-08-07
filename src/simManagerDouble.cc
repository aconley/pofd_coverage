#include<cstring>
#include<iostream>
#include<sstream>

#include<gsl/gsl_errno.h>
#include<fitsio.h>

#include "../include/simManagerDouble.h"
#include "../include/global_settings.h"
#include "../include/pofdExcept.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

//This is the function we call to find the best fitting n0
/*!
  \param[in] x Value of n0
  \param[in] params Messy array to hold everything to compute likelihood
 */
static double minfunc(double x, void* params) {
  //Basically in the order we need them
  //This is very, very ugly, but seems to be the only way
  // to interface class members to the GSL
  void** vptr = static_cast<void**>(params);

  PDFactoryDouble *pdfac = static_cast<PDFactoryDouble*>(vptr[0]);
  PDDouble *pd = static_cast<PDDouble*>(vptr[1]);
  simImageDouble *im = static_cast<simImageDouble*>(vptr[2]);

  if (x > pdfac->getMaxN0())
    throw pofdExcept("","minfunc","N0 out of prepared range",1);

  pdfac->getPD(x, *pd, true, true);
  double loglike = pd->getLogLike(*im);

  return -loglike; //Remember -- we want to minimize this
}

/////////////////////////////////////////////

/*!
  \param[in] MODELFILE Name of file containing base model
  \param[in] NSIMS Number of simulations to do
  \param[in] N0INITRANGE Initial values to maximize likelihood are 
              n0*(1-n0initrange) to n0*(1+n0initrange)
  \param[in] MAPLIKE Do make likelihood map?
  \param[in] NLIKE   Number of likelihoods in likelihood map
  \param[in] N0RANGEFRAC Likelihood map covers n0*(1+n0rangefrac) to
                         n0*(1-n0rangefrac)
  \param[in] FFTSIZE Number of elements along each axis to use in FFT
  \param[in] N1 Number of pixels along first axis in simulated image
  \param[in] N2 Number of pixels along second axis in simulated image
  \param[in] PIXSIZE Pixel size, in arcsec
  \param[in] FWHM1 Fwhm of band 1 beam, in arcsec
  \param[in] FWHM2 Fwhm of band 2 beam, in arcsec
  \param[in] SIGI1 Instrument noise (without smoothing) in Jy, band 1
  \param[in] SIGI2 Instrument noise (without smoothing) in Jy, band 2
  \param[in] SIGRNG Instrument noise range (max - min in sigma units)
  \param[in] N0 Simulated number of sources per sq deg.
  \param[in] ESMOOTH Amount of extra smoothing to apply
  \param[in] OVERSAMPLE Oversampling of simulated image
  \param[in] USEBIN Bin the data in the likelihood calculation
  \param[in] NBINS Number of bins in binned data
 */
simManagerDouble::simManagerDouble(const std::string& MODELFILE,
				   unsigned int NSIMS, double N0INITRANGE,
				   bool MAPLIKE, unsigned int NLIKE, 
				   double N0RANGEFRAC, unsigned int FFTSIZE,
				   unsigned int N1, unsigned int N2, 
				   double PIXSIZE, double FWHM1, double FWHM2, 
				   double SIGI1, double SIGI2, double SIGRNG,
				   double N0, double ESMOOTH1, double ESMOOTH2,
				   unsigned int OVERSAMPLE,
				   bool USEBIN, unsigned int NBINS) :
  nsims(NSIMS), n0initrange(N0INITRANGE), do_map_like(MAPLIKE),
  nlike(NLIKE), n0rangefrac(N0RANGEFRAC), fftsize(FFTSIZE),
  n0(N0), sig_i1(SIGI1), sig_i2(SIGI2), sigrng(SIGRNG), sig_i1_sm(SIGI1), 
  sig_i2_sm(SIGI2), fwhm1(FWHM1), fwhm2(FWHM2), pixsize(PIXSIZE),
  simim(N1, N2, PIXSIZE, FWHM1, FWHM2, SIGI1, SIGI2, SIGRNG, ESMOOTH1, 
	ESMOOTH2, OVERSAMPLE, NBINS), 
  use_binning(USEBIN), model(MODELFILE), 
  esmooth1(ESMOOTH1), esmooth2(ESMOOTH2) {

#ifdef TIMING
  initTime = getTime = getLikeTime = 0;
#endif

  if (esmooth1 > 0.0 || esmooth2 > 0.0) do_extra_smooth = true;
  else do_extra_smooth = false;
  if (do_extra_smooth) {
    bm.setFWHM(std::sqrt(fwhm1*fwhm1+esmooth1*esmooth1),
	       std::sqrt(fwhm2*fwhm2+esmooth2*esmooth2));
    std::pair<double,double> pr;
    pr = simim.getSmoothedNoiseEstimate();
    sig_i1_sm = pr.first;
    sig_i2_sm = pr.second;
  } else {
    bm.setFWHM(fwhm1, fwhm2);
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

  varr = new void*[3];
  s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

}

simManagerDouble::~simManagerDouble() {
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
}

/*!
  \param[in] verbose Output status messages as it runs; basically, what
                     simulation are we on.
*/
void simManagerDouble::doSims(bool verbose=false) {
  //Two steps; first: find the best fit n0.  Then -- maybe --
  // map out the likelihood
  const float nfwhm = 3.0; //!< How far out to go on beam
  const unsigned int nbeambins = 80; // Number of beam histogram bins
  const unsigned int max_iter = 100; //!< Maximum number of iters for min
  const unsigned int max_expiter = 20; //!< Maximum number of bracket steps
  const double reltol = 0.001; //!< Relative tolerance in n0: 0.1%
  double sigval1, sigval2; // Params to pass to minimizer
  std::pair<double, double> maxflux;

  //Turn off gsl error handler during call
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

  if (do_extra_smooth) {
    sigval1 = sig_i1_sm; 
    sigval2 = sig_i2_sm; 
  } else {
    sigval1 = sig_i1;
    sigval2 = sig_i2;
  }

  //Now, set up the parameters to pass to the minimizer.  Ugly!
  varr[0] = static_cast<void*>(&pdfac);
  varr[1] = static_cast<void*>(&pd);
  varr[2] = static_cast<void*>(&simim);

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
  maxflux = model.getMaxFluxEstimate();
  maxflux.first *= max_n0ratio;
  maxflux.second *= max_n0ratio;
#ifdef TIMING
  starttime = std::clock();
#endif
  pdfac.initPD(fftsize, sigval1, sigval2, maxflux.first, maxflux.second, 
	       init_b, model, bm, pixsize, nfwhm, nbeambins, true);
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
    //Make simulated image, don't mean sub until we have mn estimate
    // so we can estimate shifts
    simim.realize(model, n0, do_extra_smooth, true, use_binning);

    //Now set up the minimization; note this involves calling minfunc
    // so we can't do it until all the arguments are ready
    //First -- we have to bracket the minimum
    double curr_initrange, next_initrange;
    double a, b, fa, fm, fb;
    fm = minfunc(n0, params);

    //We move out logarithmically with a maximum number of steps
    unsigned int exp_iter = 0;
    next_initrange = n0initrange;
    do {
      ++exp_iter;
      curr_initrange = next_initrange;
      if (curr_initrange >= 1.0) //Shouldn't be possible, but...
	throw pofdExcept("simManagerDouble", "doSims",
			 "Logic error generating bracketing range", 1);

      //Try with current init range
      a = n0 * (1.0 - curr_initrange);
      b = n0 * (1.0 + curr_initrange);
      fa = minfunc(a, params);
      fb = minfunc(b, params);   

      //Now, update current init range.  We step a quarter way in log
      // space between here and 1.0. 
      next_initrange = exp2(0.75 * log2(curr_initrange)); 
    } while ((exp_iter < max_expiter) && ( (fa <= fm) || (fb <= fm) ));
    
    if ((fa <= fm) || (fb <= fm)) {
      std::stringstream errstr;
      errstr << "Unable to bracket minimum in -log likelihood."
	     << std::endl;
      errstr << "Final values:" << std::endl;
      errstr << "\tf(" << a << ") = " << fa << std::endl;
      errstr << "\tf(" << n0 << ") = " << fm << std::endl;
      errstr << "\tf(" << b << ") = " << fb << std::endl;
      errstr << "n0 range: " << curr_initrange << std::endl;
      throw pofdExcept("simManagerDouble", "doSims", errstr.str(), 2);
    }

    int status;
    status = gsl_min_fminimizer_set_with_values(s, &F, n0, fm,
						a, fa, b, fb);
    if (status != GSL_SUCCESS) {
      std::stringstream errstr;
      errstr << "GSL minimizer setup failed with code: " << status
	     << std::endl;
      errstr << "GSL error message: " << gsl_strerror(status);
      throw pofdExcept("simManagerDouble", "doSims", errstr.str(), 3);
    }

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
      throw pofdExcept("simManagerDouble", "doSims", 
		       errstr.str(), 2);
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
	//Note that for 1 likelihood we just hit n0, not the best fit
	//This is because we already recoreded the likelihood there,
	// so it wouldn't add anything
	min_n0[i] = n0;
	delta_n0[i] = 0.0;
      }

#ifdef TIMING
      std::clock_t starttime;
#endif
      //Compute P(D) and like
      double curr_n0;
      for (unsigned int j = 0; j < nlike; ++j) {
	curr_n0 = min_n0[i] + static_cast<double>(j)*delta_n0[i];
	
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
	likearr[i][j] =  pd.getLogLike(simim);
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
void simManagerDouble::resetTime() {
  initTime = getTime = getLikeTime = 0;
  pdfac.resetTime();
}

void simManagerDouble::summarizeTime() const {
  std::cout << "initPD time: " << 1.0*initTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "getPD time: " << 1.0*getTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "getLogLike time: " << 1.0*getLikeTime/CLOCKS_PER_SEC 
	    << std::endl;
  std::cout << "Within PDFactoryDouble: " << std::endl;
  pdfac.summarizeTime(2);
}
#endif


//Output
int simManagerDouble::writeToFits(const std::string& outputfile) const {
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
    char* tform2 = new char[ndigits+5];
    sprintf(tform2, "%uE", nlike); //E is the format code for 32 bit float
    tform[4] = tform2;
    fits_create_tbl(fp, BINARY_TBL, 0, 5, ttype, tform, NULL, NULL, &status);
    delete[] tform2;
  } else {
    char* ttype[] = {"BEST_N0", "BEST_LOGLIKE"};
    char* tform[2];
    char* tform0 = "1D";
    tform[0] = tform0; tform[1] = tform0;
    fits_create_tbl(fp, BINARY_TBL, 0, 2, ttype, tform, NULL, NULL, &status);
  }

  //Write header stuff
  unsigned int utmp;
  double dtmp;
  fits_write_key(fp, TSTRING, const_cast<char*>("MODEL"),
		 const_cast<char*>("Spline-Log Normal"), 
		 const_cast<char*>("Model type"),
		 &status);
  dtmp = n0; //Must copy to temporary for const type handling
  fits_write_key(fp, TDOUBLE, const_cast<char*>("N0"), &dtmp, 
		 const_cast<char*>("Number of sources per sq deg"), 
		 &status);

  //Base model parameters
  dtmp = model.getBaseN0();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BASEN0"), &dtmp, 
		 const_cast<char*>("Base number of sources per sq deg"), 
		 &status);

  //Sim params
  dtmp = bm.getFWHM1();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM1"), &dtmp, 
		 const_cast<char*>("Beam fwhm, band 1 [arcsec]"), 
		 &status);
  dtmp = bm.getFWHM2();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM2"), &dtmp, 
		 const_cast<char*>("Beam fwhm, band 2 [arcsec]"), 
		 &status);
  std::pair<double,double> dpr;
  dpr = simim.getBeamSum();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BMAREA1"), &dpr.first, 
		 const_cast<char*>("Beam area, band 1 [pixels]"), 
		 &status);
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BMAREA2"), &dpr.second, 
		 const_cast<char*>("Beam area, band 2 [pixels]"), 
		 &status);
  dpr = simim.getBeamSumSq();
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BMARESQ1"), &dpr.first, 
		 const_cast<char*>("Beam squared area, band 1 [pixels]"), 
		 &status);
  fits_write_key(fp, TDOUBLE, const_cast<char*>("BMARESQ2"), &dpr.second, 
		 const_cast<char*>("Beam squared area, band 2 [pixels]"), 
		 &status);

  if (do_extra_smooth) {
    int tmpi = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &tmpi,
		   const_cast<char*>("Additional smoothing applied"), 
		   &status);
    dtmp = esmooth1;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("ESMOOTH1"), &dtmp, 
		   const_cast<char*>("Extra smoothing fwhm, band 1 [arcsec]"), 
		   &status);
    fits_write_key(fp, TDOUBLE, const_cast<char*>("ESMOOTH2"), &dtmp, 
		   const_cast<char*>("Extra smoothing fwhm, band 2 [arcsec]"), 
		   &status);
  } else {
    int tmpi = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &tmpi,
		   const_cast<char*>("Additional smoothing applied"), 
		   &status);
  }
  
  dtmp = sig_i1;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGI_1"), &dtmp, 
		 const_cast<char*>("Instrument noise, band 1"), 
		 &status);
  dtmp = sig_i2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGI_2"), &dtmp, 
		 const_cast<char*>("Instrument noise, band 1"), 
		 &status);
  dtmp = sigrng;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGRNG"), &dtmp, 
		 const_cast<char*>("Instrument noise range"), 
		 &status);
  if (do_extra_smooth) {
    dtmp = sig_i1_sm;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGISM1"), &dtmp, 
		   const_cast<char*>("Smoothed instrument noise, band 1"), 
		   &status);
    dtmp = sig_i2_sm;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGISM2"), &dtmp, 
		   const_cast<char*>("Smoothed instrument noise, band 2"), 
		   &status);
  }
  
  if (use_binning) {
    int tmpi = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("USEBIN"), &tmpi,
		   const_cast<char*>("Use binned likelihood"), 
		   &status);
    utmp = simim.getNBins();
    fits_write_key(fp, TUINT, const_cast<char*>("NBINS"), &utmp,
		   const_cast<char*>("Number of bins in Likelihood"), 
		   &status);
  } else {
    int tmpi = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("USEBIN"), &tmpi,
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
  
  if (do_map_like) {
    int tmpi = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MAPLIKE"), &tmpi,
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
    int tmpi = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MAPLIKE"), &tmpi,
		   const_cast<char*>("Do map out likelihood"), 
		   &status);
  }
  
  
  fits_write_key(fp, TSTRING, const_cast<char*>("VERSION"),
		 const_cast<char*>(pofd_coverage::version), 
		 const_cast<char*>("pofd_coverage version"),
		 &status);
  
  fits_write_history(fp, 
		     const_cast<char*>("Simulation results from pofd_delta"),
		     &status);
  fits_write_date(fp, &status);
  
  //Now write out the data We write the actual array of likelihoods
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
		 const_cast<char*>("Number of 1D model knots"), 
		 &status);
  utmp = model.getNSigmaKnots();
  fits_write_key(fp, TUINT, const_cast<char*>("NSIGKNOTS"), &utmp, 
		 const_cast<char*>("Number of model sigma knots"), 
		 &status);
  utmp = model.getNOffsetKnots();
  fits_write_key(fp, TUINT, const_cast<char*>("NOFFKNOTS"), &utmp, 
		 const_cast<char*>("Number of model offset knots"), 
		 &status);
  utmp = model.getNTotalKnots();
  fits_write_key(fp, TUINT, const_cast<char*>("NTOTKNOT"), &utmp, 
		 const_cast<char*>("Number of total knots"), 
		 &status);
  fits_insert_rows(fp, 0, utmp, &status);
  for (unsigned int i = 0; i < model.getNKnots(); ++i) {
    double val;
    val = model.getKnotPosition(i);
    fits_write_col(fp, TDOUBLE, 1, i+1, 1, 1, &val, &status);
    val = model.getLog10KnotValue(i);
    fits_write_col(fp, TDOUBLE, 2, i+1, 1, 1, &val, &status);
  }
  int idxoff;
  idxoff = model.getNKnots();
  for (unsigned int i = 0; i < model.getNSigmaKnots(); ++i) {
    double val;
    val = model.getSigmaKnotPosition(i);
    fits_write_col(fp, TDOUBLE, 1, i+idxoff+1, 1, 1, &val, &status);
    val = model.getSigmaKnotValue(i);
    fits_write_col(fp, TDOUBLE, 2, i+idxoff+1, 1, 1, &val, &status);
  }
  idxoff += model.getNSigmaKnots();
  for (unsigned int i = 0; i < model.getNOffsetKnots(); ++i) {
    double val;
    val = model.getOffsetKnotPosition(i);
    fits_write_col(fp, TDOUBLE, 1, i+idxoff+1, 1, 1, &val, &status);
    val = model.getOffsetKnotValue(i);
    fits_write_col(fp, TDOUBLE, 2, i+idxoff+1, 1, 1, &val, &status);
  }

  //Close up and go home
  fits_close_file(fp, &status);
  if (status) {
    fits_report_error(stderr, status);
    return status;
  }
  return status;
}
