#include<sstream>
#include<cmath>
#include<cstring>
#include<fstream>
#include<limits>

#include<hdf5.h>

#include "../include/global_settings.h"
#include "../include/PDFactory.h"
#include "../include/pofdExcept.h"

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
  currsize = 0;

#ifdef TIMING
  resetTime();
#endif

  RFlux = NULL;
  rvals = NULL;
  rtrans = NULL;
  pofd = NULL;
  pval = NULL;

  dflux = 0.0;
  minflux_R = 0.0;

  plan = plan_inv = NULL;
  plans_valid = false;

  verbose = false;
  has_wisdom = false;
  fftw_plan_style = FFTW_MEASURE;

  mn = sg = varnoi = std::numeric_limits<double>::quiet_NaN();

  sigma = std::numeric_limits<double>::quiet_NaN();
  max_n0 = std::numeric_limits<double>::quiet_NaN();
  rinitialized = false;
  initialized = false;
  rdflux = false;
}

#ifdef TIMING
void PDFactory::resetTime() {
  RTime = p0Time = fftTime = posTime = copyTime = normTime = 0;
  meanTime = logTime = 0;
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
  std::cout << "mean time: " << prestring 
	    << 1.0*meanTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "log time: " << prestring 
	    << 1.0*logTime/CLOCKS_PER_SEC << std::endl;
}
#endif


/*
  \param[in] NSIZE new size
  \returns True if there was actually a resize
*/
//I don't check for 0 NSIZE because that should never happen
// in practice, and it isn't worth the idiot-proofing
// rtrans always gets nulled in here, then refilled when you
//  call initPD
bool PDFactory::resize(unsigned int NSIZE) {
  if (NSIZE == currsize) return false;

  if (RFlux != NULL) fftw_free(RFlux);
  RFlux = NULL;
  if (rvals != NULL) fftw_free(rvals);
  rvals = NULL;
  if (rtrans != NULL) fftw_free(rtrans);
  rtrans = NULL;
  if (pval != NULL) fftw_free(pval);
  pval = NULL;
  if (pofd != NULL) fftw_free(pofd);
  pofd = NULL;

  void *alc;
  alc = fftw_malloc(sizeof(double) * NSIZE);
  if (alc == NULL)
    throw pofdExcept("PDFactory", "strict_resize", 
		     "Failed to alloc RFlux", 1);
  RFlux = (double*) alc;
  alc = fftw_malloc(sizeof(double) * NSIZE);
  if (alc == NULL)
    throw pofdExcept("PDFactory", "strict_resize", 
		     "Failed to alloc rvals", 2);
  rvals = (double*) alc;
  alc = fftw_malloc(sizeof(double)*NSIZE);
  if (alc == NULL)
    throw pofdExcept("PDFactory", "strict_resize", 
		     "Failed to alloc pofd", 3);
  pofd = (double*) alc;
  unsigned int fsize = NSIZE / 2 + 1;
  alc = fftw_malloc(sizeof(fftw_complex)*fsize);
  if (alc == NULL)
    throw pofdExcept("PDFactory", "strict_resize", 
		     "Failed to alloc pval", 4);
  pval = (fftw_complex*) alc;
  plans_valid = false;

  currsize = NSIZE;
  rinitialized = false;
  initialized = false;
  return true;
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
    str << "Error reading wisdom file: " << filename;
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
  \param[in] n Number of elements
  \param[in] minflux Minimum flux to use in R
  \param[in] maxflux Maximum flux to use in R

  Sets up RFlux.  Rflux may not quite get to minflux/maxflux
  if minflux is negative.  
*/
void PDFactory::initRFlux(unsigned int n, double minflux, double maxflux) {
  // Make sure there is room
  resize(n);

  if (n == 0)
    throw pofdExcept("PDFactory", "initRFlux", "Invalid (0) n", 1);
  if (n == 1) {
    dflux = 0.1;
    RFlux[0] = minflux;
    minflux_R = minflux;
    return;
  }

  if (maxflux < minflux) std::swap(minflux, maxflux);
  double inm1 = 1.0 / static_cast<double>(n-1);
  dflux = (maxflux - minflux) * inm1;
  if (minflux >= 0.0) {
    RFlux[0] = minflux;
    for (unsigned int i = 1; i < n; ++i)
      RFlux[i] = static_cast<double>(i) * dflux + minflux;
    minflux_R = minflux;
  } else {
    // Here the complication is that we would really like to have
    // RFlux = 0 included in the array.  Also, we are wrapping 
    // negative fluxes around to the top of the array.
    // We do this by tweaking minflux slightly
    dflux = maxflux / (n - floor(-minflux / dflux) - 2.0);

    // Figure out what index we go up to with positive fills
    unsigned int maxpos = static_cast<unsigned int>(maxflux / dflux);
    RFlux[0] = 0.0;
    for (unsigned int i = 1; i < maxpos + 1; ++i) // Pos Rflux
      RFlux[i] = static_cast<double>(i) * dflux;

    // Figure out new minimum flux
    double wrapval = - static_cast<double>(n) * dflux;
    for (unsigned int i = maxpos + 1; i < n; ++i) // Wrapped neg Rflux
      RFlux[i] = static_cast<double>(i) * dflux + wrapval;
    minflux_R = RFlux[maxpos + 1];
  }
  initialized = false; // If there was a P(D), it's no longer valid
  rinitialized = false;
}

/*!
  \param[in] n        Size of transform 
  \param[in] minflux  Minimum flux to use in R
  \param[in] maxflux  Maximum flux to use in R
  \param[in] model    number counts model to use for fill.  Params must be set
  \param[in] bm       Histogrammed inverse beam
  \param[in] muldflux Multiply R by dflux

  dflux is also set, as well as Rflux.

  This is called by initPD, so you don't need to call this to compute
  the P(D).  The reason to use this is if you just want R, but not the P(D).

  The computed R is for the base model.
*/
void PDFactory::initR(unsigned int n, double minflux, double maxflux, 
		      const numberCounts& model, const beamHist& bm,
		      bool muldflux) {

  //Fill in Rflux values, set dflux.  Also resizes.
  initRFlux(n, minflux, maxflux);

  //Now fill in R.  
#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif
  model.getR(n, RFlux, bm, rvals);
  if (muldflux) {
    for (unsigned int i = 0; i < n; ++i) rvals[i] *= dflux;
    rdflux = true;
  } else rdflux = false;
#ifdef TIMING
  RTime += std::clock() - starttime;
#endif
  initialized = false; // If there was a P(D), it's no longer valid
  rinitialized = true;
}


/*!
  \param[in] n0 Value of n0 used
  \param[in] n Number of elements
  \param[out] pd Holds P(D) on output, normalized, mean subtracted,
                 and with positivity enforced.
  
  This does the unwrapping of the internal P(D) into pd.
*/
// This should only ever be called by getPD, so we don't really
// check the inputs
void PDFactory::unwrapPD(double n0, unsigned int n, PD& pd) const {
  
  // Our acceptance testing is a bit complicated.
  // If the minimum point is more than nsig1 away from the expected mean (0)
  //  then we just accept it.  If it is more than nsig2 away, we make sure
  //  that the min/max ratio of the P(D) along that axis is more than
  //  maxminratio.  If it is less than nsig2, we just flat out reject.
  const double nsig1 = 4.0; 
  const double nsig2 = 2.0;
  const double maxminratio = 1e5;
  
  // First, Enforce positivity
#ifdef TIMING
  starttime = std::clock();
#endif
  for (unsigned int idx = 0; idx < n; ++idx)
    if (pofd[idx] < 0) pofd[idx] = 0.0;
#ifdef TIMING
  posTime += std::clock() - starttime;
#endif

  // This is slightly tricky -- we may (probably do) have wrapping issues
  // so we need to un-wrap.  We do this by scanning from the top to
  // find the minimum in the P(D), then split the array there.
  // We also find the maximum, but don't keep track of its index
#ifdef TIMING
  starttime = std::clock();
#endif
  int mdx = static_cast<int>(n - 1);
  double cval, minval, maxval;
  minval = maxval = pofd[mdx];
  for (int i = n - 2; i >= 0; --i) {
    cval = pofd[i];
    if (cval < minval) {
      minval = cval;
      mdx = i;
    } else if (cval > maxval) maxval = cval;
  }
  unsigned int minidx = static_cast<unsigned>(mdx);

  // Make sure this is sane!
  //  Recall that varnoi was computed for max_n0, not this one
  double curr_sigma = sqrt(n0 * varnoi / max_n0 + sigma * sigma);
  double fwrap_plus = RFlux[minidx]; // Wrap in pos flux
  double fwrap_minus = static_cast<double>(n - minidx) * dflux; // Abs neg wrap
  double cs1, cs2;
  cs1 = nsig1 * curr_sigma;
  cs2 = nsig2 * curr_sigma;
  if ((fwrap_plus > cs1) || (fwrap_minus > cs1)) {
    // Worth further investigation
    if (fwrap_plus < cs2) {
      std::stringstream errstr;
      errstr << "Top wrapping problem; wrapping point at "
	     << fwrap_plus << " which is only " << fwrap_plus / curr_sigma
	     << " sigma away from expected (0) mean with sigma "
	     << curr_sigma << " at n0: " << n0;
      throw pofdExcept("PDFactory", "unwrapPD", errstr.str(), 1);
    }
    if (fwrap_minus < cs2) {
      std::stringstream errstr;
      errstr << "Bottom wrapping problem; wrapping point at "
	     << -fwrap_minus << " which is only " << fwrap_minus / curr_sigma
	     << " sigma away from expected (0) mean, with sigma "
	     << curr_sigma << " at n0: " << n0;
      throw pofdExcept("PDFactory", "unwrapPD", errstr.str(), 2);
    } 
    // Min/max ratio test
    if (maxval / minval < maxminratio) {
      std::stringstream errstr;
      errstr << "Wrapping problem with wrapping fluxes: "
	     << fwrap_plus << " and " << -fwrap_minus << " with min/max ratio: "
	     << maxval / minval << " and sigma: " << curr_sigma
	     << " with n0: " << n0;
      throw pofdExcept("PDFactory", "unwrapPD", errstr.str(), 3);
    }
  }

  // Copy over.  Things above minidx in pofd go into the bottom of
  // pd.pd_, then the stuff below that in pofd goes above that in pd.pd_
  // in the same order.
  pd.resize(n);
  double *ptr_curr, *ptr_out; // convenience vars
  ptr_curr = pofd + minidx;
  ptr_out = pd.pd_;
  //for (unsigned int i = 0; i < n - minidx; ++i)
  //  ptr_out[i] = ptr_curr[i];
  std::memcpy(ptr_out, ptr_curr, (n - minidx) * sizeof(double));
  ptr_curr = pofd;
  ptr_out = pd.pd_ + n - minidx;
  //for (unsigned int i = 0; i < minidx; ++i)
  //  ptr_out[i] = ptr_curr[i];
  std::memcpy(ptr_out, ptr_curr, minidx * sizeof(double));

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

  //Now mean subtract flux axis
#ifdef TIMING
  starttime = std::clock();
#endif
  double tmn; //True mean
  pd.getMean(tmn, false);
  if (std::isinf(tmn) || std::isnan(tmn)) {
    std::stringstream str;
    str << "Un-shift amounts not finite: " << tmn << " " << std::endl;
    str << "At length: " << n << " with noise: " << sigma;
    throw pofdExcept("PDFactory", "unwrapPD", str.str(), 3);
  }
  pd.minflux = -tmn;
#ifdef TIMING
  meanTime += std::clock() - starttime;
#endif

}

/*!
  \param[in] model Number counts model
  \param[in] bm Inverse beam histogram
  \returns Estimate of min/max flux values where R is non-zero
*/
dblpair PDFactory::getMinMaxR(const numberCounts& model,
			      const beamHist& bm) const {

  double minFRnonzero, maxFRnonzero;
  double maxknot = model.getMaxKnotPosition();
  maxFRnonzero = 1.01 * maxknot * bm.getMinMaxPos().second; 
  if (bm.hasNeg()) minFRnonzero = -1.01 * maxknot * bm.getMinMaxNeg().second;
  else minFRnonzero = 0.0;
  
  return std::make_pair(minFRnonzero, maxFRnonzero);
}


/*!
  Compute mean and variance from R

  \param[in] n       Size of transform 
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] bm      Histogrammed beam
  \param[in] range   A pair giving the range to compute R over.  Ideally
                      is the minimum and maximum range over which R is non-zero
  \returns A pair of the mean and variance
		      
  Will resize and change Rflux		      
*/
dblpair PDFactory::getRMoments(unsigned int n, const numberCounts& model, 
			       const beamHist& bm, dblpair range) {

  // Recall that the formulae for the mean and central 2nd moment
  // (not including instrumental noise) are
  //  <x> = \int x R dx
  //  <(x - <x>)^2> = \int x^2 R dx

  //We are totally going to screw things up, so mark as unclean
  initialized = rinitialized = false;  
  if (n == 0)
    throw pofdExcept("PDFactory", "getRMoments",
		     "n must be positive", 1);

  //Set R, Rflux, and dflux. 
  initR(n, range.first, range.second, model, bm); 

  //We use the trap rule. 
  double mean, var, cf, prod;
  cf = RFlux[0];
  prod = 0.5 * cf * rvals[0];
  mean = prod;
  var = cf * prod;
  for (unsigned int i = 1; i < n; ++i) {
    cf = RFlux[i];
    prod = cf * rvals[i];
    mean += prod;
    var += cf * prod;
  }
  prod = 0.5 * cf * rvals[n-1];
  mean += prod;
  var += cf * prod;

  // Include step size
  mean *= dflux;
  var *= dflux;

  return std::make_pair(mean, var);
}

/*!
  Gets ready for P(D) computation by preparing R and computing its forward
   transform
 
  \param[in] n       Size of transform 
  \param[in] inst_sigma Instrument noise, in Jy
  \param[in] maxflux Desired maximum image flux in Jy
  \param[in] maxn0   Maximum n0 supported (number of sources per area)
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] bm      Histogrammed beam

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux will end up being -less- than maxflux in practice
  by about the mean flux + 10 sigma.
 */
void PDFactory::initPD(unsigned int n, double inst_sigma, double maxflux, 
		       double maxn0, const numberCounts& model,
		       const beamHist& bm) {

  if (!model.isValid())
    throw pofdExcept("PDFactory", "initPD", "model not valid", 1);
  if (maxn0 <= 0.0)
    throw pofdExcept("PDFactory", "initPD", "invalid (non-positive) n0", 2);

  //This will cause R wrapping problems, so check maxn0 relative to
  // the model base n0 value
  if (maxn0 < model.getBaseN0()) {
    std::stringstream errstr;
    errstr << "maxn0 (" << maxn0 << ") must be greater than model base N0 ("
	   << base_n0 << ")";
    throw pofdExcept("PDFactory", "initPD", errstr.str(), 3);
  }

  //Allocate/resize internal arrays
  resize(n);
  if (rtrans == NULL) {
    //This has to be done before we plan 
    unsigned int fsize = n / 2 + 1;
    void* alc;
    alc = fftw_malloc(sizeof(fftw_complex) * fsize);
    if (alc == NULL)
      throw pofdExcept("PDFactory", "initPD", "Failed to alloc rtrans", 4);
    rtrans = (fftw_complex*) alc;
  }
  
  //Make the plans, or keep the old ones if possible
  //Do this before initializing rvals, rtrans, pofd, etc., 
  // since this will mess with those values.
  //Note that the forward transform dumps into rtrans
  // but the backwards one comes from pval.  The idea
  // is that we can re-use the forward transform, updating
  // pval as we change n0, and including instrument noise and
  // the shift.  That is, rtrans doesn't change after we call
  // initPD, but pvals will change each time we call getPD
  //If we resized, we must make the new plans because the
  // addresses changed
  int intn = static_cast<int>(n);
  if (!plans_valid) {
    if (plan != NULL) fftw_destroy_plan(plan);
    plan = fftw_plan_dft_r2c_1d(intn, rvals, rtrans,
				fftw_plan_style);
    if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
    plan_inv = fftw_plan_dft_c2r_1d(intn, pval, pofd,
				    fftw_plan_style);
    if (plan == NULL) {
      std::stringstream str;
      str << "Plan creation failed for forward transform of size: " << n;
      if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			  << " that size";
      throw pofdExcept("PDFactory", "initPD", str.str(), 5);
    }
    if (plan_inv == NULL) {
      std::stringstream str;
      str << "Plan creation failed for inverse transform of size: " << n;
      if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			  << " that size";
      throw pofdExcept("PDFactory", "initPD", str.str(), 6);
    }
    plans_valid = true;
  }

  base_n0 = model.getBaseN0();
  double n0ratio = maxn0 / base_n0;

  //Estimate the mean and standard deviation of the resulting P(D) using R.
  //Since the usage model for pofd_coverage is to call initPD once and then
  // re-use it a bunch of times, we can afford to compute R twice to get
  // a better estimate.  Start by getting the non-zero range of R
  if (!bm.hasPos())
    throw pofdExcept("PDFactory", "initPD", 
		     "Code assumes positive beam is present", 7);
  dblpair rangeR = getMinMaxR(model, bm);

  // Estimate the mean model flux and sigma.
  // We do this by computing R and using it's moments
  // Note we compute these for the maximum n0
  dblpair mom; // Mean, var
  mom = getRMoments(n, model, bm, rangeR);
  mn = n0ratio * mom.first; // mean goes as n0 -- so this is for maxn0 model
  varnoi = n0ratio * mom.second;
  sg = sqrt(varnoi + inst_sigma * inst_sigma); // for maxn0

  // Now compute the range we will ask for R over in the actual
  // computation.  Ensure zero padding at the top.
  double maxflux_R = maxflux + pofd_coverage::n_zero_pad * sg;

  //Now prepare final R.  Note this uses the base model, even
  // though we computed the maximum R value based on the maximum n0 value
  // We go from the smallest non-zero R to the value we just found.
  // Note that this automatically pads R with zeros as needed!
  // We actually ask for R * dflux, because that's what we want to transform
  initR(n, rangeR.first, maxflux_R, model, bm, true);

  //Decide if we will shift, and if so by how much
  // The idea is to shift the mean to zero -- but we only
  // do the shift if the sigma is larger than one actual step size
  // because otherwise we can't represent it well.
  doshift = (sg > dflux) && (fabs(mn) > dflux);
  if (doshift) shift = - mn; 
  else shift=0.0;

  if (verbose) {
    std::cout << " Initial mean estimate: " << mn << std::endl;
    std::cout << " Initial stdev estimate: " << sg << std::endl;
    if (doshift)
      std::cout << " Additional shift applied: " << shift << std::endl;
    else std::cout << " Not applying additional shift" << std::endl;
  }

  //Compute forward transform of this r value, store in rtrans
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute(plan); 
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif
  
  max_n0 = maxn0;
  sigma = inst_sigma;
  initialized = true;
}

/*!
  Calculates P(D) for a dataset based on already computed R (by initPD)
 
  \param[in] n0 Number of sources in current computation
  \param[out] pd Holds P(D) on output
  \param[in] setLog If true, pd is log(P(D) on output; convenient
              for likelihood evaluation.

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux often won't quite match the target values.

  You must call initPD first.
*/
void PDFactory::getPD(double n0, PD& pd, bool setLog) {

  // The basic idea is to compute the P(D) from the previously filled
  // R values, adding in noise and all that fun stuff, filling pd
  // for output

  if (!initialized)
    throw pofdExcept("PDFactory", "getPD",
		     "Must call initPD first", 1);
  if (n0 > max_n0) {
    std::stringstream errstr("");
    errstr << "N_0 " << n0
	   << " larger than maximum prepared value " << max_n0
	   << std::endl;
    errstr << "initPD should have been called with at least " << n0;
    throw pofdExcept("PDFactory", "getPD", errstr.str(), 2);
  }
  double n0ratio = n0 / base_n0;
  
  //Output array from 2D FFT is n/2+1
  unsigned int n = currsize;
  unsigned int ncplx = n/2 + 1;
      
  //Calculate p(omega) = exp(r(omega) - r(0)),
  // convolving together all the bits into pval, which
  // is what we will transform back into pofd.
  // The forward transform of R are stored in rtrans in the
  //  usual sign order, but for the base model.
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

  double r0 = n0ratio * rtrans[0][0]; //r[0] is pure real
  double iflux = pofd_coverage::two_pi / (n * dflux);

  if (doshift) {
    //This should be the most common case.  Note that the shift
    // was computed based on max_n0, not the current or base n0
    // so we have to correct for that
    double curr_shift = n0 * shift / max_n0;
    double sigfac = 0.5 * sigma * sigma;
    double w, expfac, rval, ival;
    for (unsigned int idx = 1; idx < ncplx; ++idx) {
      w = iflux * static_cast<double>(idx);
      rval = n0ratio * rtrans[idx][0] - r0 - sigfac * w * w;
      ival = n0ratio * rtrans[idx][1] - curr_shift * w;
      expfac = exp(rval);
      pval[idx][0] = expfac * cos(ival);
      pval[idx][1] = expfac * sin(ival);
    } 
  } else {
    double expfac, ival;
    if (sigma > 0.0) {
      double w, rval;
      double sigfac = 0.5 * sigma * sigma;
      for (unsigned int idx = 1; idx < ncplx; ++idx) {
	w = iflux * static_cast<double>(idx);
	rval = n0ratio * rtrans[idx][0] - r0 - sigfac * w * w;
	ival = n0ratio * rtrans[idx][1];
	expfac = exp(rval);
	pval[idx][0] = expfac * cos(ival);
	pval[idx][1] = expfac * sin(ival);
      }
    } else {
      // No instrument sigma, simple
      for (unsigned int idx = 1; idx < ncplx; ++idx) {
	expfac = exp(n0ratio * rtrans[idx][0] - r0);
	ival = n0ratio * rtrans[idx][1];
	pval[idx][0] = expfac * cos(ival);
	pval[idx][1] = expfac * sin(ival);
      }
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

  // Copy into output variable, also normalizing, mean subtracting, 
  // making positive
  unwrapPD(n0, n, pd);

  //Turn PD to log for more efficient log computation of likelihood
#ifdef TIMING
  starttime = std::clock();
#endif
  if (setLog) pd.applyLog(false);
#ifdef TIMING
  logTime += std::clock() - starttime;
#endif
}
 
/*!
  Writes out current R to a text file
 
  \param[in] filename File to write to

  You must call initPD or initR first, or bad things will probably happen.
*/
void PDFactory::writeRToFile(const std::string& filename) const {
  if (!rinitialized )
    throw pofdExcept("PDFactory", "writeRToFile",
		     "Must call initPD or initR first", 1);

  std::ofstream ofs(filename.c_str());
  if (!ofs)
    throw pofdExcept("PDFactory", "writeRToFile",
		     "Couldn't open output file", 2);

  ofs << currsize << std::endl;
  if (rdflux) {
    double idflux = 1.0 / dflux;
    for (unsigned int i = 0; i < currsize; ++i)
      ofs << RFlux[i] << " " << rvals[i] * idflux << std::endl;
  } else {
    for (unsigned int i = 0; i < currsize; ++i)
      ofs << RFlux[i] << " " << rvals[i] << std::endl;
  }

  ofs.close();
}

/*!
  Writes out current R to a HDF5 file
 
  \param[in] filename File to write to

  You must call initPD or initR first, or bad things will probably happen.
*/
void PDFactory::writeRToHDF5(const std::string& filename) const {
  if (!rinitialized )
    throw pofdExcept("PDFactory", "writeRToHDF5",
		     "Must call initPD or initR first", 1);

  hid_t file_id;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);

  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw pofdExcept("PDFactory", "writeToHDF5",
		     "Failed to open HDF5 file to write", 2);
  }

  // Write it as one dataset -- Rflux, R. 
  hsize_t adims;
  hid_t mems_id, att_id, dat_id;
  
  // Properties
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(file_id, "dflux", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux);
  H5Aclose(att_id);
  att_id = H5Acreate2(file_id, "N0", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &base_n0);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  // Rflux
  adims = currsize;
  mems_id = H5Screate_simple(1, &adims, NULL);
  dat_id = H5Dcreate2(file_id, "RFlux", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, RFlux);
  H5Dclose(dat_id);

  // R -- which we may need to copy to remove the dflux
  dat_id = H5Dcreate2(file_id, "R", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (rdflux) {
    double* tmp = new double[currsize];
    double idflux = 1.0 / dflux;
    for (unsigned int i = 0; i < currsize; ++i) tmp[i] = rvals[i] * idflux;
    H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,  H5P_DEFAULT, tmp);
    delete[] tmp;
  } else
    H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,  H5P_DEFAULT, rvals);
  H5Dclose(dat_id);
  H5Sclose(mems_id);

  // Done
  H5Fclose(file_id);
}
