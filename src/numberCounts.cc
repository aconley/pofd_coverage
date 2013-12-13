#include<limits>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

#include "../include/pofdExcept.h"
#include "../include/numberCounts.h"
#include "../include/utility.h"

const double numberCounts::ftol = 1e-4;

//Function to pass to GSL integrator
static double evalPowfNKnotsSpline(double, void*); //!< Evaluates f^pow dN/dS

/*!
  \param[in] modelfile Name of file to read base model from
*/
numberCounts::numberCounts(const std::string& modelfile,
			   unsigned int NINTERP) : gen_ninterp(NINTERP) {

  //We have to read in the input file
  std::ifstream ifs(modelfile.c_str());
  if (!ifs)
    throw pofdExcept("numberCounts","numberCounts",
		     "Unable to open input model file",1);
  
  //Do the read into temporary vectors, then copy
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double currval;
  std::vector<double> kp, kv;
  const unsigned int nreq = 2;
  while (!ifs.eof()) {
    std::getline(ifs,line);
    if (line[0] == '#' || line[0] == '%') continue; //Skip comments
    
    //Parse into words, stipping spaces
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#' || words[0][0] == '%') continue; //Comment line
    if (words.size() < nreq) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currval;
    kp.push_back(currval);
    str.str(words[1]); str.clear(); str >> currval;
    kv.push_back(currval);
  }
  ifs.close();
  nknots = kp.size();
  if (nknots < 2)
    throw pofdExcept("numberCounts","numberCounts",
		     "Need two valid knots from modelfile",2);
  knotpos = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    knotpos[i] = kp[i];
  logknotpos = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    logknotpos[i] = log2(knotpos[i]);
  
  //These are read in as log_10, so we must convert to log_2
  logknotvals = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = kv[i] * pofd_coverage::logfac;

  knotvals = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotvals[i] = exp2(logknotvals[i]);

  acc = gsl_interp_accel_alloc();
  splinelog = gsl_spline_alloc(gsl_interp_cspline,
			       static_cast<size_t>(nknots));
  gsl_spline_init(splinelog, logknotpos, logknotvals, 
		  static_cast<size_t>(nknots));

  // Don't set histogrammed beam currently
  nbm = 0;
  bm_wts_pos = bm_wts_neg = NULL;
  inv_bm_pos = inv_bm_neg = NULL;

  // Need to set base_n0, base_flux, etc.
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
  base_n0 = 1.0; //Fake value to allow isValid to pass
  base_n0 = splineInt(0.0);
  base_flux = splineInt(1.0);
  base_fluxsq = splineInt(2.0);

  // Set up variables for generating sources
  // Note we integrate -down- in flux because the number counts
  // are expected to drop rapidly to lower values
  // We use log spaced points rather than linearly spaced to concentrate
  // more of the tabulated points at low fluxes where there should be
  // more sources
  if (gen_ninterp < 2)
    throw pofdExcept("numberCounts","numberCounts",
		     "Number of interpolation generation points must be >2", 3);
  // Two ingredients -- the cumulative normalized probablity (gen_interp_cumsum)
  // and the corresponding flux densities.
  // First, set up the flux densities -- note they are the log2 values!
  gen_interp_flux = new double[gen_ninterp];
  double lmaxf = logknotpos[nknots-1];
  double lminf = logknotpos[0];
  double dlogf = (lmaxf - lminf) / static_cast<double>(gen_ninterp-1);
  gen_interp_flux[0] = lmaxf;
  for (unsigned int i = 1; i < gen_ninterp-1; ++i)
    gen_interp_flux[i] = lmaxf - static_cast<double>(i) * dlogf;
  gen_interp_flux[gen_ninterp-1] = lminf;

  // Now integrate down the number counts
  // Get cumulative distribution function using trapezoidal rule, at first
  // un-normalized.
  gen_interp_cumsum = new double[gen_ninterp];
  double cumsum, currf, prevf, df, currcnts, prevcnts;
  prevf = exp2(gen_interp_flux[0]);
  prevcnts = 0.0;
  cumsum = gen_interp_cumsum[0] = 0.0; //First bit has 0 counts
  for (unsigned int i = 1; i < gen_ninterp; ++i) {
    currf = exp2(gen_interp_flux[i]);
    df = prevf - currf; // Sign flip because we are going down
    currcnts = exp2(gsl_spline_eval(splinelog, gen_interp_flux[i], acc));
    cumsum += 0.5 * df * (currcnts + prevcnts);
    gen_interp_cumsum[i] = cumsum;
    prevf = currf;
    prevcnts = currcnts;
  }
  //Normalize to the range 0 to 1; can skip the first entry (which is 0)
  double norm = 1.0 / cumsum;
  for (unsigned int i = 1; i < gen_ninterp-1; ++i) gen_interp_cumsum[i] *= norm;
  gen_interp_cumsum[gen_ninterp-1] = 1.0;

  // Now stick in the interpolation object.  Use linear interpolation
  // to avoid the probability oscillating negative or something
  gen_interp = gsl_interp_alloc(gsl_interp_linear,
				static_cast<size_t>(gen_ninterp));
  gen_interp_acc = gsl_interp_accel_alloc();
  gsl_interp_init(gen_interp, gen_interp_cumsum, gen_interp_flux,
		  static_cast<size_t>(gen_ninterp));
}

numberCounts::~numberCounts() {
  delete[] knotpos;
  delete[] logknotpos;
  delete[] knotvals;
  delete[] logknotvals;
  gsl_interp_accel_free(acc);
  gsl_spline_free(splinelog);
  gsl_interp_accel_free(gen_interp_acc);
  gsl_interp_free(gen_interp);
  gsl_integration_workspace_free(gsl_work);
  delete[] gen_interp_cumsum;
  delete[] gen_interp_flux;
  delete[] varr;
  if (bm_wts_pos != NULL) delete[] bm_wts_pos;
  if (inv_bm_pos != NULL) delete[] inv_bm_pos;
  if (bm_wts_neg != NULL) delete[] bm_wts_neg;
  if (inv_bm_neg != NULL) delete[] inv_bm_neg;
}

bool numberCounts::isValid() const {
  if (std::isnan(base_n0)) return false;
  if (std::isinf(base_n0)) return false;
  if (base_n0 <= 0) return false;
  return true;
}

/*!
  Computes
  \f[
   \int dS\, S^{\alpha} \frac{dN}{dS}
  \f]

  \param[in] alpha  Power of flux
  \returns Integral
 */
double numberCounts::splineInt(double alpha) const {
  if (nknots < 2) return std::numeric_limits<double>::quiet_NaN();
  if (! isValid()) return std::numeric_limits<double>::quiet_NaN();
  double result, error;
  void *params;
  gsl_function F;
  double minknot = knotpos[0];
  double maxknot = knotpos[nknots-1];

  varr[0] = static_cast<void*>(&alpha);
  varr[1] = static_cast<void*>(splinelog);
  varr[2] = static_cast<void*>(acc);
  varr[3] = static_cast<void*>(&minknot);
  varr[4] = static_cast<void*>(&maxknot);
  params = static_cast<void*>(varr);

  F.function = &evalPowfNKnotsSpline;
  F.params = params;

  gsl_integration_qag(&F, minknot, maxknot, 0.0, 1e-5, 1000,
  		      GSL_INTEG_GAUSS41, gsl_work, &result, &error); 
  return result;
}



/*!
  \returns The base number of sources per area
 */
double numberCounts::getBaseN0() const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  return base_n0;
}

/*!
  \returns Flux density per square degree for base model
 */
double numberCounts::getBaseFluxPerArea() const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  return base_flux;
}

/*!
  \returns Flux density^2 per square degree for base model
 */
double numberCounts::getBaseFluxSqPerArea() const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  return base_fluxsq;
}

/*!
  \param[in] flux Flux density number counts are desired at

  \returns Differential number counts at flux density flux for base model
 */

double numberCounts::getdNdS(double flux) const {
  if (flux < knotpos[0]) return 0.0;
  if (flux >= knotpos[nknots-1]) return 0.0;
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  return exp2(gsl_spline_eval(splinelog, log2(flux), acc)); 
}

/*!
  \param[in] x The value R is desired for
  \param[in] bm The beam
  \param[in] pixsize The pixel size, in arcseconds
  \param[in] nfwhm The number of beam fwhm to use in the computation
  \param[in] nbins The number of bins in the beam histogramming
  \param[in] filterparam The filtering param.


  \returns R(x) computed for the base input model.
*/
double numberCounts::getR(double x, const beam& bm,
			  double pixsize, double nfwhm,
			  unsigned int nbins, double filterparam) const {
  return getR(x, bm, pixsize, nfwhm, nbins, 1, filterparam);
}

/*!
  \param[in] x The value R is desired for
  \param[in] bm The beam
  \param[in] pixsize The pixel size, in arcseconds
  \param[in] nfwhm The number of beam fwhm to use in the computation
  \param[in] nbins The number of bins in the beam histogramming
  \param[in] oversamp Oversampling amount for beam
  \param[in] filterparam The filtering param.

  \returns R(x) computed for the base input model.
*/
double numberCounts::getR(double x, const beam& bm,
			  double pixsize, double nfwhm,
			  unsigned int nbins,
			  unsigned int oversamp,
			  double filterparam) const {

  //This can be done much more efficiently (see the old pofd_mcmc
  // code), but the idea here is to only really compute R once
  // and then reuse it again and again so we prize simplicity over
  // efficiency.

  if (pixsize <= 0.0)
    throw pofdExcept("numberCounts", "getR", "Invalid (non-positive) pixsize",
		     1);
  if (nfwhm <= 0.0)
    throw pofdExcept("numberCounts", "getR", "Invalid (non-positive) nfwhm",
		     2);
  if (nbins == 0)
    throw pofdExcept("numberCounts", "getR", "Invalid (zero) nbins",
		     3);
  if (!isValid())
    throw pofdExcept("numberCounts", "getR", "Invalid model",
		     4);

  double s_min = knotpos[0];
  double s_max = knotpos[nknots-1];
  if (x >= s_max) return 0.0;
  if (x <= 0.0) return 0.0;

  if (nbm < nbins) {
    //Must expand
    if (bm_wts_pos == NULL) delete[] bm_wts_pos;
    if (inv_bm_pos == NULL) delete[] inv_bm_pos;
    bm_wts_pos = new unsigned int[nbins];
    inv_bm_pos = new double[nbins];
    if (bm_wts_neg == NULL) delete[] bm_wts_neg;
    if (inv_bm_neg == NULL) delete[] inv_bm_neg;
    bm_wts_neg = new unsigned int[nbins];
    inv_bm_neg = new double[nbins];
    nbm = nbins;
  }

  //Get number of pixels out we will go
  unsigned int npix = static_cast<unsigned int>(nfwhm * bm.getFWHM()/pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;
  unsigned int nnonzero_pos, nnonzero_neg;

  //Load inverse histogrammed beam
  bm.getBeamHist(npix, pixsize, nbins, oversamp, nnonzero_pos, 
		 bm_wts_pos, inv_bm_pos, nnonzero_neg, bm_wts_neg, 
		 inv_bm_neg, filterparam, true);

  //And now the actual computation
  double cval, cR, R, ibm;
  R = 0.0;
  for (unsigned int i = 0; i < nnonzero_pos; ++i) {
    ibm = inv_bm_pos[i]; //1 / eta
    cval = x * ibm;
    if (cval < s_min) continue;
    if (cval >= s_max) continue;
    cR = exp2(gsl_spline_eval(splinelog, log2(cval), acc));
    R += bm_wts_pos[i] * cR * ibm;
  }
  for (unsigned int i = 0; i < nnonzero_neg; ++i) {
    ibm = inv_bm_neg[i];
    cval = x * ibm;
    if (cval < s_min) continue;
    if (cval >= s_max) continue;
    cR = exp2(gsl_spline_eval(splinelog, log2(cval), acc));
    R += bm_wts_neg[i] * cR * ibm;
  }

  double prefac = pixsize / 3600.0;
  return prefac * prefac * R;
}


/*!\brief Get number of source responses, vector version 

  \param[in] n The number of values to compute R for
  \param[in] minflux The minimum value to compute R for
  \param[in] maxflux The maximum value to compute R for
  \param[in] bm The beam
  \param[in] pixsize The pixel size, in arcseconds
  \param[in] nfwhm The number of beam fwhm to use in the computation
  \param[in] nbins The number of bins in the beam histogramming
  \param[out] R The returned values of R, of length n.  Must
                be pre-allocated by the caller
  \param[in] filterparam The filter model param.

  R is computed for the base model		 
*/
void numberCounts::getR(unsigned int n, double minflux,
			double maxflux, const beam& bm, 
			double pixsize, double nfwhm,
			unsigned int nbins, 
			double* R, double filterparam) const {
  getR(n, minflux, maxflux, bm, pixsize, nfwhm, nbins, 1, R, filterparam);
}

/*!\brief Get number of source responses, vector version 

  \param[in] n The number of values to compute R for
  \param[in] minflux The minimum value to compute R for
  \param[in] maxflux The maximum value to compute R for
  \param[in] bm The beam
  \param[in] pixsize The pixel size, in arcseconds
  \param[in] nfwhm The number of beam fwhm to use in the computation
  \param[in] nbins The number of bins in the beam histogramming
  \param[in] oversamp The oversampling in the beam
  \param[out] R The returned values of R, of length n.  Must
                be pre-allocated by the caller
  \param[in] filterparam The filter model param.		

  R is computed for the base model		 
*/

void numberCounts::getR(unsigned int n, double minflux,
			double maxflux, const beam& bm, 
			double pixsize, double nfwhm,
			unsigned int nbins, unsigned int oversamp,
			double* R, double filterparam) const {

  //As for the scalar version of getR above, this could be done much more
  // efficiently but we aim for simplicity rather than efficiency
  // in this implementation.
  
  if (n == 0) return;
  if (pixsize <= 0.0)
    throw pofdExcept("numberCounts", "getR", "Invalid (non-positive) pixsize",
		     1);
  if (nfwhm <= 0.0)
    throw pofdExcept("numberCounts", "getR", "Invalid (non-positive) nfwhm",
		     2);
  if (nbins == 0)
    throw pofdExcept("numberCounts", "getR", "Invalid (zero) nbins",
		     3);
  

  double s_min = knotpos[0];
  double s_max = knotpos[nknots-1]; 
  double dflux = (maxflux - minflux) / static_cast<int>(n - 1);
  if (n == 1) dflux = 0.0;

  //Space for the inverse beam
  if (nbm < nbins) {
    //Must expand
    if (bm_wts_pos == NULL) delete[] bm_wts_pos;
    if (inv_bm_pos == NULL) delete[] inv_bm_pos;
    bm_wts_pos = new unsigned int[nbins];
    inv_bm_pos = new double[nbins];
    if (bm_wts_neg == NULL) delete[] bm_wts_neg;
    if (inv_bm_neg == NULL) delete[] inv_bm_neg;
    bm_wts_neg = new unsigned int[nbins];
    inv_bm_neg = new double[nbins];
    nbm = nbins;
  }
  unsigned int npix = static_cast<unsigned int>(nfwhm * bm.getFWHM()/pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;
  unsigned int nnonzero_pos, nnonzero_neg;
  bm.getBeamHist(npix, pixsize, nbins, oversamp, 
		 nnonzero_pos, bm_wts_pos, inv_bm_pos, 
		 nnonzero_neg, bm_wts_neg, inv_bm_neg,
		 filterparam, true);
  
  double prefac = pixsize / 3600.0;
  prefac = prefac * prefac;

  //And now the actual computation.  We loop over each flux
  double cflux, cval, cR, ibm, workR;
  
  for (unsigned int i = 0; i < n; ++i) {
    cflux = minflux + static_cast<double>(i) * dflux;
    if (cflux <= 0.0 || cflux >= s_max) {
      //R is zero for this x
      R[i] = 0.0;
      continue;
    }
    workR = 0.0;
    for (unsigned int j = 0; j < nnonzero_pos; ++j) {
      ibm = inv_bm_pos[j]; //1 / eta
      cval = cflux * ibm;
      if (cval < s_min) continue;
      if (cval >= s_max) continue;
      cR = exp2(gsl_spline_eval(splinelog, log2(cval), acc));
      workR += bm_wts_pos[j] * cR * ibm;
    }
    for (unsigned int j = 0; j < nnonzero_neg; ++j) {
      ibm = inv_bm_neg[j]; 
      cval = cflux * ibm;
      if (cval < s_min) continue;
      if (cval >= s_max) continue;
      cR = exp2(gsl_spline_eval(splinelog, log2(cval), acc));
      workR += bm_wts_neg[j] * cR * ibm;
    }
    R[i] = prefac * workR;
  }
}

/*!
  \param[in] udev Uniform deviate [0,1).  

  udev is not checked for validity in the interests of speed, 
  so screwy things will happen if you provide an invalid one.
 */
double numberCounts::genSource(double udev) const {
  // Note that we can ignore curr_n0 here, since changing n_0
  // just changes the number of sources, not their distribution
  // in flux density.  Recall that gen_interp_flux is the log2 flux.
  return exp2(gsl_interp_eval(gen_interp, gen_interp_cumsum,
			      gen_interp_flux, udev, gen_interp_acc));
  
}

bool numberCounts::writeToStream(std::ostream& os) const {
  os << "Base model parameters: " << std::endl;
  os << " " << std::left << std::setw(13) << "#Flux knot" << "  "
     << std::setw(13) << "Log10 Knot value" << std::endl;
  for (unsigned int i = 0; i < nknots; ++i)
    os << " " << std::left << std::setw(13) << knotpos[i] << "  "
       << std::setw(13) << logknotvals[i] * pofd_coverage::ilogfac 
       << std::endl; 
  return true;
}


std::ostream& operator<<(std::ostream& os, const numberCounts& b) {
  b.writeToStream(os);
  return os;
}

static double evalPowfNKnotsSpline(double x, void* params) {
  //Params are:
  // parmas[0] power
  // params[1] spline (log2)
  // params[2] accelerator
  // params[3] minknot
  // params[4] maxknot
  //But this really has to be an array of pointers to void to work
  void** vptr = static_cast<void**>(params);

  double power = *static_cast<double*>(vptr[0]);
  gsl_spline* spl = static_cast<gsl_spline*>(vptr[1]);
  gsl_interp_accel* acc = static_cast<gsl_interp_accel*>(vptr[2]);
  double minknot = *static_cast<double*>(vptr[3]);
  double maxknot = *static_cast<double*>(vptr[4]);
  
  if (x < minknot || x >= maxknot) return 0.0;

  double splval = exp2(gsl_spline_eval(spl, log2(x), acc));

  if (fabs(power) < 1e-2) return splval;
  if (fabs(power-1.0) < 1e-2) return x * splval;
  if (fabs(power-2.0) < 1e-2) return x * x * splval;
  return pow(x, power) * splval;
}
