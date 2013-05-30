#include<sstream>
#include<cmath>
#include<cstring>
#include<fstream>
#include<limits>

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
  if (RFlux != NULL) FFTWFREE(RFlux);
  if (rvals != NULL) FFTWFREE(rvals);
  if (rtrans != NULL) FFTWFREE(rtrans);
  if (pofd != NULL) FFTWFREE(pofd);
  if (pval != NULL) FFTWFREE(pval);
  if (plan != NULL) FFTWDESTROYPLAN(plan); 
  if (plan_inv != NULL) FFTWDESTROYPLAN(plan_inv);
}

void PDFactory::init() {
  plan_size = 0;
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
  fftw_plan_style = FFTW_MEASURE;

  sigma = 0.0;
  max_n0 = 0.0;
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

  if (RFlux != NULL) FFTWFREE(RFlux);
  if (rvals != NULL) FFTWFREE(rvals);
  if (rtrans != NULL) FFTWFREE(rtrans);
  rtrans = NULL;
  isRTransAllocated = false;
  if (pval != NULL) FFTWFREE(pval);
  if (pofd != NULL) FFTWFREE(pofd);

  RFlux = (double*) FFTWMALLOC(sizeof(double) * NSIZE);
  rvals = (FFTFLOAT*) FFTWMALLOC(sizeof(FFTFLOAT) * NSIZE);
  pofd = (FFTFLOAT*) FFTWMALLOC(sizeof(FFTFLOAT)*NSIZE);
  unsigned int fsize = NSIZE/2+1;
  pval = (FFTWCOMPLEX*) FFTWMALLOC(sizeof(FFTWCOMPLEX)*fsize);
  plans_valid = false;

  currsize = NSIZE;
  initialized = false;
}

/*!
  \param[in] filename Name of wisdom file
*/
bool PDFactory::addWisdom(const std::string& filename) {
  FILE *fp = NULL;
  fp = fopen(filename.c_str(), "r");
  if (!fp) {
    std::stringstream str;
    str << "Error opening wisdom file: " << filename;
    throw pofdExcept("PDFactory","addWisdom",str.str(),1);
  }
  if (FFTWIMPORTWIS(fp) == 0) {
    std::stringstream str;
    str << "Error reading wisdom file: " << wisdom_file;
    throw pofdExcept("PDFactory","addWisdom",str.str(),1);
  }
  fclose(fp);
  fftw_plan_style = FFTW_WISDOM_ONLY;
  has_wisdom = true;
  wisdom_file = filename;
  if (plan != NULL) {
    FFTWDESTROYPLAN(plan); 
    plan = NULL;
  }
  if (plan_inv != NULL) {
    FFTWDESTROYPLAN(plan_inv);
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
  \param[in] pixsize Pixel size in arcseconds
  \param[in] nfwhm   Number of FWHM out to go in beam
  \param[in] nbins   Number of bins used in histogrammed beam
  \param[out] vec    Returns moments
  \param[in] nmom    Number of moments.  Note we start with the 0th one,
                      so if you want the 2nd moment, nmom must be at least 3.

*/
void PDFactory::getRIntegrals(unsigned int n, FFTFLOAT maxflux,
			      const numberCounts& model, const beam& bm, 
			      double pixsize, double nfwhm, unsigned int nbins,
			      std::vector<FFTFLOAT>& vec, unsigned int nmom) {
  //We are totally going to screw things up, so mark as unclean
  initialized = false;  
  if (n == 0)
    throw pofdExcept("PDFactory","getRIntegrals",
		     "n must be positive", 1);
  if (nmom == 0)
    throw pofdExcept("PDFactory","getRIntegrals",
		     "nmom must be positive", 2);

  vec.resize(nmom);
  initR(n, maxflux, model, bm, pixsize, nfwhm, nbins); //Set R and dflux

  //We use the trap rule.  Note R already is multiplied by dflux
  //Do zeroth integral
  double mom;
  mom = 0.5 * rvals[0];
  for (unsigned int i = 1; i < n-1; ++i) mom += rvals[i];
  mom += 0.5 * rvals[n-1];
  vec[0] = static_cast<FFTFLOAT>(mom);

  //Need a working array for the flux product.  The idea is to steadily
  // accumulate Rflux^n * R in wrkarr.
  double *wrkarr = new double[n];
  for (unsigned int i = 0; i < n; ++i)
    wrkarr[i] = rvals[i];

  for (unsigned int i = 1; i < nmom; ++i) {
    for (unsigned int j = 0; j < n; ++j)
      wrkarr[j] *= RFlux[j];  //so for i=1 it is Rflux*R, for i=2 Rflux^2*R
    mom = 0.5 * wrkarr[0];
    for (unsigned int j = 1; j < n-1; ++j) mom += wrkarr[j];
    mom += 0.5 * wrkarr[n-1];
    vec[i] = static_cast<FFTFLOAT>(mom);
  }
}


/*!
  Prepare R.
 
  \param[in] n       Size of transform 
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] beam    Beam fwhm
  \param[in] pixsize Pixel size in arcseconds
  \param[in] nfwhm   Number of FWHM out to go in beam
  \param[in] nbins   Number of bins used in histogrammed beam

  The R value that is set is actually R * dflux for convenience.
  dflux is also set.

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux will often end up a bit off from the desired value.
*/
void PDFactory::initR(unsigned int n, FFTFLOAT maxflux, 
		      const numberCounts& model, const beam& bm,
		      double pixsize, double nfwhm, 
		      unsigned int nbins) {

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
  model.getR(n, 0.0, maxflux, bm, pixsize, nfwhm, nbins, rvals);
  for (unsigned int i = 1; i < n; ++i) rvals[i] *= dflux;
#ifdef TIMING
  RTime += std::clock() - starttime;
#endif
}


/*!
  Gets ready for P(D) computation by preparing R and computing its forward
   transform
 
  \param[in] n       Size of transform 
  \param[in] inst_sigma Instrument noise, in Jy
  \param[in] maxflux Desired maximum image flux in Jy
  \param[in] maxn0   Maximum n0 supported (number of sources per area)
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] beam    Beam fwhm
  \param[in] pixsize Pixel size in arcseconds
  \param[in] nfwhm   Number of FWHM out to go in beam
  \param[in] nbins   Number of bins used in histogrammed beam

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux will end up being -less- than maxflux in practice
  by about the mean flux + 10 sigma.
 */
void PDFactory::initPD(unsigned int n, FFTFLOAT inst_sigma, 
		       FFTFLOAT maxflux, FFTFLOAT maxn0, 
		       const numberCounts& model,
		       const beam& bm, double pixsize, double nfwhm,
		       unsigned int nbins) {

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
  if (! isRTransAllocated ) {
    if (rtrans != NULL) FFTWFREE(rtrans);
    unsigned int fsize = n / 2 + 1;
    rtrans = (FFTWCOMPLEX*) FFTWMALLOC(sizeof(FFTWCOMPLEX)*fsize);
    isRTransAllocated = true;
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
  if ( (!plans_valid) || (plan_size != n) ) {
    if (plan != NULL) FFTWDESTROYPLAN(plan);
    plan = FFTWDFTR2C1D(intn, rvals, rtrans, fftw_plan_style);
    if (plan_inv != NULL) FFTWDESTROYPLAN(plan_inv);
    plan_inv = FFTWDFTC2R1D(intn, pval, pofd, fftw_plan_style);
    if (plan == NULL) {
      std::stringstream str;
      str << "Plan creation failed for forward transform of size: " << n;
      if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			  << " that size";
      throw pofdExcept("PDFactory","initPD",str.str(),4);
    }
    if (plan_inv == NULL) {
      std::stringstream str;
      str << "Plan creation failed for inverse transform of size: " << n;
      if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			  << " that size";
      throw pofdExcept("PDFactory","initPD",str.str(),5);
    }
    plans_valid = true;
  }

  base_n0 = model.getBaseN0();
  FFTFLOAT n0ratio = maxn0 / base_n0;

  //We will iterate to try to get the maximum right.  Given R,
  // we can compute the mean and standard deviation of the resulting P(D).
  // Thanks to shifting (to avoid aliasing and include mean subtraction),
  // what we compute R out to to achieve this is slightly tricky.
  //Since the usage model for pofd_coverage is to compute R once and then
  // re-use it a bunch of times, we can afford to compute R twice to get
  // a better estimate

  //Estimate the mean model flux and sigma crudely
  double maxflux_R, s_ave, est_shift, var;
  s_ave = n0ratio * model.getBaseFluxPerArea();
  mn =  s_ave * bm.getEffectiveArea();
  var = model.getBaseFluxSqPerArea() * bm.getEffectiveArea() - s_ave*s_ave;
  if (var <= 0) var = 0.0;
  sg = sqrt(n0ratio * n0ratio * var + inst_sigma*inst_sigma);
  est_shift = mn + pofd_coverage::n_sigma_shift * sg;
  maxflux_R = maxflux + est_shift;

  //Compute R integrals to update estimates for shift
  std::vector<FFTFLOAT> mom(3); //0th, 1st, 2nd moment
  getRIntegrals(n, maxflux_R, model, bm, pixsize, nfwhm, nbins, mom, 3);

  mn = n0ratio * mom[1];
  //Note the variance goes up as n0ratio, not the sigma
  sg = sqrt(n0ratio * mom[2] + inst_sigma * inst_sigma);
  est_shift = mn + pofd_coverage::n_sigma_shift * sg;
  maxflux_R = maxflux + est_shift;

  //Now prepare final R.  Note this uses the base model, even
  // though we computed the maximum R value based on the maximum n0 value
  // The returned value is R * dflux, and dflux is set
  initR(n, maxflux_R, model, bm, pixsize, nfwhm, nbins);

  //Decide if we will shift and pad, and if so by how much
  //Only do shift if the noise is larger than one actual step size
  // Otherwise we can't represent it well.
  bool dopad = (inst_sigma > dflux);
  doshift = ( dopad && ( mn < pofd_coverage::n_sigma_shift * sg) );
  if ( doshift ) shift = pofd_coverage::n_sigma_shift*sg - mn; else
    shift=0.0;

  if (verbose) {
    std::cout << " Initial mean estimate: " << mn << std::endl;
    std::cout << " Initial stdev estimate: " << sg << std::endl;
    if (doshift)
      std::cout << " Additional shift applied: " << shift << std::endl;
    else std::cout << " Not applying additional shift" << std::endl;
  }

  //Make sure that maxflux is large enough that we don't get
  // bad aliasing wrap from the top around into the lower P(D) values.
  if (maxflux_R <= pofd_coverage::n_sigma_pad * sg)
    throw pofdExcept("PDFactory", "initPD", "Top wrap problem", 4);

  //The other side of the equation is that we want to zero-pad the
  // top, and later discard that stuff.
  // The idea is as follows:
  // the 'target mean' of the calculated P(D) will lie at mn+shift.
  // We assume that anything within n_sigma_pad*sg
  // is 'contaminated'.  That means that, if n_sigma_pad*sg >
  // mn+shift, there will be some wrapping around the bottom of the P(D)
  // to contaminate the top by an amount n_sigma_pad*sg - (mn+shift).
  // We therefore zero pad and discard anything above
  // maxflux - (n_sigma_pad*sg - (mn+shift))
  if (dopad) {
    double contam = pofd_coverage::n_sigma_pad*sg - (mn+shift);
    if (contam < 0) maxidx = n; else {
      double topflux = maxflux_R - contam;
      if (topflux < 0)
	throw pofdExcept("PDFactory", "initPD", "Padding problem", 6);
      maxidx = static_cast< unsigned int>(topflux/dflux);
      if (maxidx > n)
	throw pofdExcept("PDFactory","initPD", "Padding problem", 7);
      //Actual padding
      for (unsigned int i = maxidx; i < n; ++i)
	rvals[i] = 0.0;
    }
  } else maxidx = n;

  //Compute forward transform of this r value, store in rtrans
#ifdef TIMING
  starttime = std::clock();
#endif
  FFTWEXECUTE(plan); 
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif
  
  max_n0 = maxn0;
  sigma = inst_sigma;
  plan_size = n;
  initialized = true;
}

/*!
  Calculates P(D) for a dataset using direct fills
 
  \param[in] n0 Number of sources in current computation
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
void PDFactory::getPD(double n0, PD& pd, bool setLog, bool edgeFix) {

  // The basic idea is to compute the P(D) from the previously filled
  // R values, adding in noise and all that fun stuff, filling pd
  // for output

  if (!initialized)
    throw pofdExcept("PDFactory","getPD",
		     "Must call initPD first",1);
  if (n0 > max_n0) {
    std::stringstream errstr("");
    errstr << "N_0 " << n0
	   << " larger than maximum prepared value " << max_n0
	   << std::endl;
    errstr << "initPD should have been called with at least " << n0;
    throw pofdExcept("PDFactory","getPD",errstr.str(),2);
  }
  double n0ratio = n0 / base_n0;
  
  //Output array from 2D FFT is n/2+1
  unsigned int n = plan_size;
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
    //This should be the most common case,
    // and corresponds to having some noise
    double sigfac = 0.5*sigma*sigma;
#pragma omp parallel
    {
      double w, expfac, rval, ival;
#pragma omp for
      for (unsigned int idx = 1; idx < ncplx; ++idx) {
	w = iflux * static_cast<double>(idx);
	rval = n0ratio * rtrans[idx][0] - r0 - sigfac*w*w;
	ival = n0ratio * rtrans[idx][1] - shift*w;
	expfac = exp(rval);
	pval[idx][0] = static_cast<FFTFLOAT>(expfac * cos(ival));
	pval[idx][1] = static_cast<FFTFLOAT>(expfac * sin(ival));
      } 
    }
  } else {
    //No shift, sigma must be zero
#pragma omp parallel
    {
      double expfac, ival;
#pragma omp for
      for (unsigned int idx = 1; idx < ncplx; ++idx) {
	expfac = exp(n0ratio * rtrans[idx][0] - r0);
	ival = n0ratio * rtrans[idx][1];
	pval[idx][0] = static_cast<FFTFLOAT>(expfac*cos(ival));
	pval[idx][1] = static_cast<FFTFLOAT>(expfac*sin(ival));
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
  FFTWEXECUTE(plan_inv); //overwrites pofd with reverse transform of pval
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
  if (verbose) std::cout << " Expected mean: " << shift+mn 
			 << " Realized mean: " << tmn << std::endl;
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
 
/*!
  Writes out current R
 
  \param[in] filename File to write to

  You must call initPD first, or bad things will probably happen.
*/
void PDFactory::writeRToFile(const std::string& filename) const {
  if (! initialized )
    throw pofdExcept("PDFactory","writeRToFile",
		     "Must call initPD first",1);

  std::ofstream ofs(filename.c_str());
  if (!ofs)
    throw pofdExcept("PDFactory","writeRToFile",
		     "Couldn't open output file",2);

  //Recall after initPD we are storing R * dflux
  ofs << plan_size << std::endl;
  for (unsigned int i = 0; i < plan_size; ++i)
    ofs << RFlux[i] << " " << rvals[i] / dflux << std::endl;

  ofs.close();
}
