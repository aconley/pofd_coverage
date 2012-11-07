#include<sstream>
#include<cmath>
#include<cstring>
#include<limits>

#include<global_settings.h>
#include<PDFactory.h>
#include<pofdExcept.h>

/*!
 */
PDFactory::PDFactory() { init(); }

/*
  \param[in] wisfile Wisdom file filename
 */
PDFactory::PDFactory(const std::string& wisfile) {
  init();
  addWisdom(wisfile);
}

PDFactory::~PDFactory() {
  if (RFlux != NULL) fftw_free(RFlux);
  if (rvals != NULL) fftw_free(rvals);
  if (rtrans != NULL) fftw_free(rtrans);
  if (pofd != NULL) fftw_free(pofd);
  if (pval != NULL) fftw_free(pval);
  if (plan != NULL) fftw_destroy_plan(plan); 
  if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
}

void PDFactory::init() {
  lastfftlen = 0;
  currsize = 0;

#ifdef TIMING
  resetTime();
#endif

  RFlux = NULL;
  rvals = NULL;
  rtrans = NULL;
  isRTransAllocated = false;
  pofd = NULL;
  pval = NULL;

  dflux = 0.0;

  plan = plan_inv = NULL;
  plans_valid = false;

  verbose = false;
  has_wisdom = false;
  fftw_plan_style = FFTW_ESTIMATE;

  max_sigma = 0.0;
  initialized = false;
}

#ifdef TIMING
void PDFactory::resetTime() {
  RTime = p0Time = fftTime = posTime = copyTime = normTime = 0;
  edgeTime = meanTime = logTime = 0;
}

void PDFactory::summarizeTime(unsigned int nindent) const {
  std::string prestring(nindent,' ');
  prestring = "  ";
    
  std::cout << "R time: " << prestring 
	    << 1.0*RTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "p0 time: " << prestring 
	    << 1.0*p0Time/CLOCKS_PER_SEC << std::endl;
  std::cout << "fft time: " << prestring 
	    << 1.0*fftTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "pos time: " << prestring 
	    << 1.0*posTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "copy time: " << prestring 
	    << 1.0*copyTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "norm time: " << prestring 
	    << 1.0*normTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "edge time: " << prestring 
	    << 1.0*edgeTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "mean time: " << prestring 
	    << 1.0*meanTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "log time: " << prestring 
	    << 1.0*logTime/CLOCKS_PER_SEC << std::endl;
}
#endif


/*
  Doesn't force resize if there is already enough room available

  \param[in] NSIZE new size
  \returns True if a resize was needed
 */
bool PDFactory::resize(unsigned int NSIZE) {
  if (NSIZE > currsize) {
    strict_resize(NSIZE);
    return true;
  } else return false;
}

/*!
  Forces a resize

  \param[in] NSIZE new size (must be > 0)
 */
//I don't check for 0 NSIZE because that should never happen
// in practice, and it isn't worth the idiot-proofing
// rtrans always gets nulled in here, then refilled when you
//  call initPD
void PDFactory::strict_resize(unsigned int NSIZE) {
  if (NSIZE == currsize) return;

  if (RFlux != NULL) fftw_free(RFlux);
  if (rvals != NULL) fftw_free(rvals);
  if (rtrans != NULL) fftw_free(rtrans);
  rtrans = NULL;
  isRTransAllocated = false;
  if (pval != NULL) fftw_free(pval);
  if (pofd != NULL) fftw_free(pofd);

  RFlux = (double*) fftw_malloc(sizeof(double)*NSIZE);
  rvals = (double*) fftw_malloc(sizeof(double)*NSIZE);
  pofd = (double*) fftw_malloc(sizeof(double)*NSIZE);
  unsigned int fsize = NSIZE/2+1;
  pval = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fsize);
  plans_valid = false;

  currsize = NSIZE;
  initialized = false;
}

/*!
  \param[in] filename Name of wisdom file
*/
bool PDFactory::addWisdom(const std::string& filename) {
  FILE *fp = NULL;
  fp = fopen( filename.c_str(), "r" );
  if (!fp) {
    std::stringstream str;
    str << "Error opening wisdom file: " << filename;
    throw pofdExcept("PDFactory","addWisdom",str.str(),1);
  }
  if (fftw_import_wisdom_from_file(fp) == 0) {
    std::stringstream str;
    str << "Error reading wisdom file: " << wisdom_file;
    throw pofdExcept("PDFactory","addWisdom",str.str(),1);
  }
  fclose(fp);
  fftw_plan_style = FFTW_WISDOM_ONLY;
  has_wisdom = true;
  wisdom_file = filename;
  if (plan != NULL) {
    fftw_destroy_plan(plan); 
    plan = NULL;
  }
  if (plan_inv != NULL) {
    fftw_destroy_plan(plan_inv);
    plan_inv = NULL;
  }
  plans_valid = false;
  initialized = false;

  return true;
}

/*!
  Compute integrals of R

  \param[in] n       Size of transform 
  \param[in] maxflux Maximum flux generated in R
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] beam    Beam fwhm
  \param[out] vec    Returns moments
  \param[in] nmom    Number of moments

*/
void PDFactory::getRIntegrals(unsigned int n, double maxflux,
			      const numberCounts& model, const beam& bm, 
			      std::vector<double>& vec, 
			      unsigned int nmom) {
  //We are totally going to screw things up, so mark as unclean
  initialized = false;  
  if (n == 0)
    throw pofdExcept("PDFactory","getRIntegrals",
		     "n must be positive", 1);
  if (nmom == 0)
    throw pofdExcept("PDFactory","getRIntegrals",
		     "nmom must be positive", 2);

  vec.resize(nmom);
  initR(n, maxflux, model, bm);

  double dx = RFlux[1] - RFlux[0];

  //We use the trap rule
  //Do zeroth integral
  double mom;
  mom = 0.5 * rvals[0];
  for (unsigned int i = 1; i < n-1; ++i)
    mom += rvals[i];
  mom += 0.5 * rvals[n-1];
  vec[0] = mom * dx;

  //Need a working array for the flux product
  double *wrkarr = new double[n];
  for (unsigned int i = 0; i < n; ++i)
    wrkarr[i] = 1.0;

  for (unsigned int i = 1; i < nmom; ++i) {
    for (unsigned int j = 0; j < n; ++j)
      wrkarr[j] *= RFlux[j];  //so for i=1 it is Rflux, for i=2 Rflux^2, etc.
    mom = 0.5 * wrkarr[0] * rvals[0];
    for (unsigned int j = 1; j < n-1; ++j)
      mom += wrkarr[j] * rvals[j];
    mom += 0.5 * wrkarr[n-1] * rvals[n-1];
    vec[i] = mom * dx;
  }
}


/*!
  Prepare R.
 
  \param[in] n       Size of transform 
  \param[in] maxflux Maximum flux generated in R
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] beam    Beam fwhm

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux will end up being -less- than maxflux in practice
  by about the mean flux + 10 sigma.
 */
void PDFactory::initR(unsigned int n, double maxflux, 
		      const numberCounts& model, const beam& bm ) {

  //Make sure we have enough room
  resize(n);

  double inm1 = 1.0 / static_cast<double>(n-1);
  dflux = maxflux * inm1;

  //Fill in flux values (note minflux is 0)
  RFlux[0] = 0.0;
  for (unsigned int i = 1; i < n; ++i)
    RFlux[i] = static_cast<double>(i)*dflux;

  //Now fill in R.  
#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif
  model.getR(n, 0.0, maxflux, bm, rvals);
#ifdef TIMING
  RTime += std::clock() - starttime;
#endif
}


/*!
  Gets ready for P(D) computation by preparing R
 
  \param[in] n       Size of transform 
  \param[in] sigma   Maximum allowed sigma
  \param[in] maxflux Maximum flux generated in R
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] beam    Beam fwhm

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux will end up being -less- than maxflux in practice
  by about the mean flux + 10 sigma.
 */
void PDFactory::initPD(unsigned int n, double sigma,
		       double maxflux, const numberCounts& model,
		       const beam& bm ) {

  //Fill R
  initR(n, maxflux, model, bm);

  //Use R to estimate the mean and standard deviation of
  // the resulting P(D) -> <D> = \int x R dx, Var[D] = \int x^2 R dx
  mn = rvals[1]; //Noting that RFlux[0] = 0
  for (unsigned int i = 2; i < n-1; ++i)
    mn += rvals[i]*static_cast<double>(i);
  mn += 0.5*rvals[n-1]*static_cast<double>(n-1);
  mn *= dflux*dflux; //Twice -- one for the step size, one for RFlux step

  double varR = rvals[1]; //Again, using the fact that RFlux[0] = 0
  for (unsigned int i = 2; i < n-1; ++i)
    varR += rvals[i]*static_cast<double>(i)*static_cast<double>(i);
  varR += 0.5*rvals[n-1]*static_cast<double>(n-1)*
    static_cast<double>(n-1);
  varR *= dflux*dflux*dflux;
  double sigR = sqrt(varR + sigma*sigma);

  //Multiply R by dflux factor to represent the actual
  // number of sources in each bin
  for (unsigned int i = 0; i < n; ++i)
    rvals[i] *= dflux;

  //Make the plans, or keep the old ones if possible
  //Note that the forward transform dumps into out_part,
  // but the backwards one comes from pval.  The idea
  // is that out_part holds the working bit.  These are
  // convolved together into pval.  This is inefficient
  // if there is only one sign present, but the special
  // case doesn't seem worth the effort
  //If we resized, we must make the new plans because the
  // addresses changed
  //We will have to use the advanced interfact to point
  // specifically at the rvals subindex we are using on the
  // forward plan, but the backwards plan is fine
  int intn = static_cast<int>(n);
  if ( (!plans_valid) || (lastfftlen != n) ) {
    if (plan != NULL) fftw_destroy_plan(plan);
    plan = fftw_plan_dft_r2c_1d(intn, rvals, rtrans,
				fftw_plan_style);
    if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
    plan_inv = fftw_plan_dft_c2r_1d(intn, pval, pofd,
				    fftw_plan_style);
    if (plan == NULL) {
      std::stringstream str;
      str << "Plan creation failed for forward transform of size: " << 
	n << std::endl;
      throw pofdExcept("PDFactory","initPD",str.str(),32);
    }
    if (plan_inv == NULL) {
      std::stringstream str;
      str << "Plan creation failed for inverse transform of size: " << 
	n << std::endl;
      throw pofdExcept("PDFactory","initPD",str.str(),64);
    }
    plans_valid = true;
  }

  //Decide if we will shift and pad, and if so by how much
  //Only do shift if the noise is larger than one actual step size
  // Otherwise we can't represent it well.
  bool dopad = (sigma > dflux);
  doshift = ( dopad && ( mn < pofd_delta::n_sigma_shift*sigR) );
  if ( doshift ) shift = pofd_delta::n_sigma_shift*sigR - mn; else
    shift=0.0;

  if (verbose) {
    std::cout << " Initial mean estimate: " << mn << std::endl;
    std::cout << " Initial stdev estimate: " << sigR << std::endl;
    if (doshift)
      std::cout << " Additional shift applied: " << shift << std::endl;
  }

  //Make sure that maxflux is large enough that we don't get
  // bad aliasing wrap from the top around into the lower P(D) values.
  if (maxflux <= pofd_delta::n_sigma_pad*sigR)
    throw pofdExcept("PDFactory","initPD","Top wrap problem",
		     128);

  //The other side of the equation is that we want to zero-pad the
  // top, and later discard that stuff.
  // The idea is as follows:
  // the 'target mean' of the calculated P(D) will lie at mn+shift.
  // We assume that anything within n_sigma_pad2d*sg
  // is 'contaminated'.  That means that, if n_sigma_pad*sg >
  // mn+shift, there will be some wrapping around the bottom of the P(D)
  // to contaminate the top by an amount n_sigma_pad*sg - (mn+shift).
  // We therefore zero pad and discard anything above
  // maxflux - (n_sigma_pad*sg - (mn+shift))
  if (dopad) {
    double contam = pofd_delta::n_sigma_pad*sigR - (mn+shift);
    if (contam < 0) maxidx = n; else {
      double topflux = maxflux - contam;
      if (topflux < 0)
	throw pofdExcept("PDFactory","initPD","Padding problem",
			 256);
      maxidx = static_cast< unsigned int>( topflux/dflux );
      if (maxidx > n)
	throw pofdExcept("PDFactory","initPD",
			 "Padding problem",
			 512);
      //Actual padding
      for (unsigned int i = maxidx; i < n; ++i)
	rvals[i] = 0.0;
    }
  } else maxidx = n;

  //Allocate memory if needed; this is a way of not allocating
  // these until we run into a beam that needs them
  if (! isRTransAllocated ) {
    if (rtrans != NULL) fftw_free(rtrans);
    unsigned int fsize = currsize/2+1;
    rtrans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fsize);
    isRTransAllocated = true;
  }

  //Compute forward transform of this r value, store in rtrans
  //Have to use argument version, since the address of rtrans can move
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute_dft_r2c(plan,rvals,rtrans); 
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif
  
  lastfftlen = n;
  max_sigma = sigma;
  initialized = true;
}

/*!
  Calculates P(D) for a dataset using direct fills
 
  \param[in] sigma Instrument noise sigma
  \param[out] pd Holds P(D) on output
  \param[in] setLog If true, pd is log(P(D) on output; convenient
              for likelihood evaluation.
  \param[in] edgeFix  Apply a fix to the lower edges to minimize wrapping
                      effects using a Gaussian to each row/col

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux often won't quite match the target values.

  You must call initPD first, or bad things will probably happen.
*/
void PDFactory::getPD( double sigma, PD& pd, bool setLog, 
		       bool edgeFix) {

  // The basic idea is to compute the P(D) from the previously filled
  // R values, adding in noise and all that fun stuff, filling pd
  // for output

  if (! initialized )
    throw pofdExcept("PDFactory","getPD",
		     "Must call initPD first",1);
  if (sigma > max_sigma) {
    std::stringstream errstr("");
    errstr << "Sigma value " << sigma
	   << " larger than maximum prepared value " << max_sigma
	   << std::endl;
    errstr << "initPD should have been called with at least " << sigma;
    throw pofdExcept("PDFactory","getPD",errstr.str(),2);
  }

  //Output array from 2D FFT is n/2+1
  unsigned int n = lastfftlen;
  unsigned int ncplx = n/2 + 1;
      
  //Calculate p(omega) = exp( r(omega) - r(0) ),
  // convolving together all the bits into pval, which
  // is what we will transform back into pofd.
  // The forward transform of R are stored in rtrans in the
  //  usual sign order.
  // There are some complications because of shifts and all that.
  // The frequencies are:
  //  f = i/dflux*n  
  // We work in w instead of f (2 pi f)
  // and actually compute 
  //  exp( r(omega) - r(0) - i*shift*omega - 1/2 sigma^2 omega^2 )
  if (verbose) std::cout << "  Computing p(w)" << std::endl;

#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif

  double r0, expfac, rval, ival;
  r0 = rtrans[0][0]; //r[0] is pure real
  double iflux = pofd_delta::two_pi / (n * dflux);

  if (doshift) {
    //This should be the most common case,
    // and corresponds to having some noise
    double sigfac = 0.5*sigma*sigma;
    double w;
    for (unsigned int idx = 1; idx < ncplx; ++idx) {
      w    = iflux * static_cast<double>(idx);
      rval = rtrans[idx][0] - r0 - sigfac*w*w;
      ival = rtrans[idx][1] - shift*w;
      expfac = exp( rval );
      pval[idx][0] = expfac*cos(ival);
      pval[idx][1] = expfac*sin(ival);
    } 
  } else {
    //No shift, sigma must be zero
    for (unsigned int idx = 1; idx < ncplx; ++idx) {
      expfac = exp(rtrans[idx][0]-r0);
      ival = rtrans[idx][1];
      pval[idx][0] = expfac*cos(ival);
      pval[idx][1] = expfac*sin(ival);
    }
  }
  //p(0) is special
  pval[0][0] = 1.0;
  pval[0][1] = 0.0;

#ifdef TIMING
  p0Time += std::clock() - starttime;
#endif

  //Transform back
  if (verbose) std::cout << " Reverse transform" << std::endl;
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute(plan_inv); //overwrites pofd with reverse transform of pval
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif

  //Enforce positivity
#ifdef TIMING
  starttime = std::clock();
#endif
  for (unsigned int idx = 0; idx < n; ++idx)
    if (pofd[idx] < 0) pofd[idx] = 0.0;
#ifdef TIMING
  posTime += std::clock() - starttime;
#endif

  //Copy into output variable
#ifdef TIMING
  starttime = std::clock();
#endif
  pd.resize(maxidx);
  for (unsigned int i = 0; i < maxidx; ++i)
    pd.pd_[i] = pofd[i];
  pd.logflat = false;
  pd.minflux = 0.0; pd.dflux = dflux;
#ifdef TIMING
  copyTime += std::clock() - starttime;
#endif

  //Normalize
#ifdef TIMING
  starttime = std::clock();
#endif
  pd.normalize();
#ifdef TIMING
  normTime += std::clock() - starttime;
#endif

  //Fix up the edge
#ifdef TIMING
  starttime = std::clock();
#endif
  if (edgeFix) pd.edgeFix();
#ifdef TIMING
  edgeTime += std::clock() - starttime;
#endif

  //Now mean subtract flux axis
#ifdef TIMING
  starttime = std::clock();
#endif
  double tmn; //True mean
  pd.getMean(tmn,false);
  if ( std::isinf(tmn) || std::isnan(tmn) ) {
    std::stringstream str;
    str << "Un-shift amounts not finite: " << tmn << " " << std::endl;
    str << "At length: " << n << " with noise: " << sigma;
    throw pofdExcept("PDFactory","getPD",str.str(),8);
  }
  if (verbose) {
    std::cerr << " Expected mean: " << shift+mn 
	      << " Realized mean: " << tmn << std::endl;
  }
  pd.minflux = -tmn;
#ifdef TIMING
  meanTime += std::clock() - starttime;
#endif

  //Turn PD to log for more efficient log computation of likelihood
#ifdef TIMING
  starttime = std::clock();
#endif
  if (setLog) pd.applyLog(false);
#ifdef TIMING
  logTime += std::clock() - starttime;
#endif
}
 
