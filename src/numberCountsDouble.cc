#include<cmath>
#include<iomanip>
#include<sstream> 
#include<fstream>
#include<limits>
#include<cstring>

#include "../include/utility.h"
#include "../include/numberCountsDouble.h"
#include "../include/global_settings.h"
#include "../include/pofdExcept.h"

/* Given two flux densities, tell which sign component they match */
unsigned int signComp(double x1, double x2) {
  if (x1 >= 0) {
    if (x2 >= 0) return 0;
    else return 1;
  } else {
    if (x2 >= 0) return 2;
    else return 3;
  }
}


//Function to pass to GSL integrator
/*! \brief Evaluates flux1^power1 * exp(const1*mu + const2*sigma^2) dN/dS1 */
static double evalPowfNDoubleLogNormal(double, void*); 

const unsigned int numberCountsDouble::nvarr = 17;
const double numberCountsDouble::ftol = 1e-4;

/*!
  \param[in] modelfile Name of file to read base model from
*/
numberCountsDouble::numberCountsDouble(const std::string& modelfile,
                                       unsigned int NINTERP) :
  gen_ninterp(NINTERP) {
  //Read in file
  unsigned int nk, ns, no; //Number of knots, sigmas, offsets
  std::vector<double> wvec1, wvec2;
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double currval;

  std::ifstream initfs(modelfile.c_str());
  if (!initfs) {
    initfs.close();
    std::stringstream errmsg;
    errmsg << "Unable to open file:" << modelfile << std::endl;
    throw pofdExcept("numberCountsDouble", "numberCountsDouble",
                       errmsg.str(), 1);
  }

  //Read in number of knots in band1, sigmas, offsets
  bool has_nknots = false;
  while (!initfs.eof()) {
    std::getline(initfs,line);
    if (line[0] == '#' || line[0] == '%') continue; //Comment
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words.size() < 3) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> nk;
    str.str(words[1]); str.clear(); str >> ns;
    str.str(words[2]); str.clear(); str >> no;
    has_nknots = true;
    break;
  }
  if (!has_nknots) {
    initfs.close();
    throw pofdExcept("numberCountsDouble", "numberCountsDouble",
                     "Unable to find number of knots line", 2);
  }
  if (nk < 2) {
    initfs.close();
    throw pofdExcept("numberCountsDouble", "numberCountsDouble",
                       "Need at least 2 band 1 knots", 3);
  }
  if (ns < 1) {
    initfs.close();
    throw pofdExcept("numberCountsDouble", "numberCountsDouble",
                       "Need at least one sigma color model knot", 4);

  }
  if (no < 1) {
    initfs.close();
    throw pofdExcept("numberCountsDouble", "numberCountsDouble",
                       "Need at least one offset color model knot", 5);
  }
  
  //Read in values
  while (!initfs.eof()) {
    std::getline(initfs,line);
    if (line[0] == '#' || line[0] == '%') continue; //Comment
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#' || words[0][0] == '%') continue; //Comment line
    if (words.size() < 2) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currval;
    wvec1.push_back(currval);
    str.str(words[1]); str.clear(); str >> currval;
    wvec2.push_back(currval); 
  }
  initfs.close();

  unsigned int ntot = nk+ns+no;
  if (wvec1.size() != ntot) {
    std::stringstream errstr;
    errstr << "Expected " << ntot << " values from " << modelfile << ", got: " 
           << wvec1.size();
    throw pofdExcept("numberCountsDouble", "numberCountsDouble",
                       errstr.str(), 6);
  }

  //Set up band 1
  nknots = nk;
  knotpos = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotpos[i] = wvec1[i];
  logknotpos = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    logknotpos[i] = log2(knotpos[i]);

  //These are read in as log_10, so we must convert to log_2
  logknotvals = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = wvec2[i] * pofd_coverage::logfac;

  //Get non-log values
  knotvals = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotvals[i] = exp2(logknotvals[i]);

  acc = gsl_interp_accel_alloc();
  splinelog = gsl_spline_alloc(gsl_interp_cspline,
                               static_cast<size_t>(nknots));
  gsl_spline_init(splinelog, logknotpos, logknotvals, 
                  static_cast<size_t>(nknots));

  //Now set up band 2 info
  // First, the sigma spline
  nsigma = ns;
  sigmapos = new double[nsigma];
  for (unsigned int i = 0; i < nsigma; ++i) sigmapos[i] = wvec1[i + nk];
  sigmavals = new double[nsigma];
  for (unsigned int i = 0; i < nsigma; ++i) sigmavals[i] = wvec2[i + nk];

  sigmainterp = nullptr;
  accsigma = nullptr;
  if (nsigma == 2)
      sigmainterp = gsl_interp_alloc(gsl_interp_linear,
                                     static_cast<size_t>(nsigma));
  else if (nsigma > 2)
    sigmainterp = gsl_interp_alloc(gsl_interp_cspline,
                                   static_cast<size_t>(nsigma));
  if (nsigma > 1) {
    gsl_interp_init(sigmainterp, sigmapos, sigmavals,
                    static_cast<size_t>(nsigma));
    accsigma = gsl_interp_accel_alloc();
  }

  // Offset spline
  noffset = no;
  offsetpos = new double[noffset];
  for (unsigned int i = 0; i < noffset; ++i) offsetpos[i] = wvec1[i+nk+ns];
  offsetvals = new double[noffset];
  for (unsigned int i = 0; i < noffset; ++i) offsetvals[i] = wvec2[i+nk+ns];
  offsetinterp = nullptr;
  accoffset = nullptr;
  if (noffset == 2)
      offsetinterp = gsl_interp_alloc(gsl_interp_linear,
                                      static_cast<size_t>(noffset));
  else if (noffset > 2)
    offsetinterp = gsl_interp_alloc(gsl_interp_cspline,
                                    static_cast<size_t>(noffset));
  if (noffset > 1) {
    gsl_interp_init(offsetinterp, offsetpos, offsetvals,
                    static_cast<size_t>(noffset));
    accoffset = gsl_interp_accel_alloc();
  }

  //Make sure what we read makes sense
  if (!isValidLoaded())
    throw pofdExcept("numberCountsDouble", "numberCountsDouble",
                     "Invalid base model parameters", 7);

  // Need to set base_n0, base_flux, etc.
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[nvarr];
  base_n0 = 1.0; //Fake value to allow isValid to pass
  base_n0 = powerInt(0.0, 0.0);
  base_flux1 = powerInt(1.0, 0.0);
  base_fluxsq1 = powerInt(2.0, 0.0);
  base_flux2 = powerInt(0.0, 1.0);
  base_fluxsq2 = powerInt(0.0, 2.0);

  // Set up variables for generating sources in band 1
  // Note we integrate -down- in flux because the number counts
  // are expected to drop rapidly to lower values
  // We use log spaced points rather than linearly spaced to concentrate
  // more of the tabulated points at low fluxes where there should be
  // more sources
  if (gen_ninterp < 2)
    throw pofdExcept("numberCountsDouble","numberCountsDouble",
                     "Number of interpolation generation points must be >2", 8);
  // Two ingredients -- the cumulative normalized probablity (gen_interp_cumsum)
  // and the corresponding flux densities.
  // First, set up the flux densities -- note they are the log2 values!
  gen_interp_flux = new double[gen_ninterp];
  double lmaxf = logknotpos[nknots - 1];
  double lminf = logknotpos[0];
  double dlogf = (lmaxf - lminf) / static_cast<double>(gen_ninterp - 1);
  gen_interp_flux[0] = lmaxf;
  for (unsigned int i = 1; i < gen_ninterp-1; ++i)
    gen_interp_flux[i] = lmaxf - static_cast<double>(i) * dlogf;
  gen_interp_flux[gen_ninterp - 1] = lminf;

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

numberCountsDouble::~numberCountsDouble() {
  if (acc != nullptr) gsl_interp_accel_free(acc);
  if (splinelog != nullptr) gsl_spline_free(splinelog);
  if (accsigma != nullptr) gsl_interp_accel_free(accsigma);
  if (sigmainterp != nullptr) gsl_interp_free(sigmainterp);
  if (accoffset != nullptr) gsl_interp_accel_free(accoffset);
  if (offsetinterp != nullptr) gsl_interp_free(offsetinterp);
  gsl_integration_workspace_free(gsl_work);
  delete[] varr;

  delete[] knotpos;
  delete[] logknotpos;
  delete[] knotvals;
  delete[] logknotvals;
}

bool numberCountsDouble::isValidLoaded() const {
  for (unsigned int i = 0; i < nknots; ++i)
    if (std::isnan(knotpos[i])) return false;
  if (knotpos[0] <= 0.0 ) return false;
  for (unsigned int i = 1; i < nknots; ++i)
    if (knotpos[i] <= knotpos[i-1]) return false;
  for (unsigned int i = 0; i < nknots; ++i)
    if (std::isnan(logknotvals[i])) return false;
  
  for (unsigned int i = 0; i < nsigma; ++i)
    if (std::isnan(sigmapos[i])) return false;
  if (sigmapos[0] <= 0.0) return false;
  for (unsigned int i = 1; i < nsigma; ++i)
    if (sigmapos[i] <= sigmapos[i-1]) return false;
  for (unsigned int i = 0; i < nsigma; ++i)
    if (sigmavals[i] <= 0.0) return false;
  for (unsigned int i = 0; i < nsigma; ++i)
    if (std::isnan(sigmavals[i])) return false;

  for (unsigned int i = 0; i < noffset; ++i)
    if (std::isnan(offsetpos[i])) return false;
  if (offsetpos[0] <= 0.0) return false;
  for (unsigned int i = 1; i < noffset; ++i)
    if (offsetpos[i] <= offsetpos[i-1]) return false;
  for (unsigned int i = 0; i < noffset; ++i)
    if (std::isnan(offsetvals[i])) return false;
  return true;
}

/*!
  \returns True if the model parameters are valid
 */
bool numberCountsDouble::isValid() const {
  //Already covered by isValidLoaded.
  return true;
}

/*!
  \param[in] f1 Band 1 flux value to evaluate sigma for
  \returns Value of band 2 color model sigma at f1

  Assumes validity already checked
*/
double numberCountsDouble::getSigmaInner(double f1) const {
  if (nsigma == 1) return sigmavals[0];
  if (f1 <= sigmapos[0]) return sigmavals[0];
  if (f1 >= sigmapos[nsigma-1]) return sigmavals[nsigma - 1];
  return gsl_interp_eval(sigmainterp, sigmapos, sigmavals, 
                         f1, accsigma );
}

/*!
  \param[in] f1 Band 1 flux value to evaluate offset for
  \returns Value of band 2 color model offset at f1

  Assumes validity already checked
*/
double numberCountsDouble::getOffsetInner(double f1) const {
  if (noffset == 1) return offsetvals[0];
  if (f1 <= offsetpos[0]) return offsetvals[0];
  if (f1 >= offsetpos[noffset-1]) return offsetvals[noffset - 1];
  return gsl_interp_eval(offsetinterp, offsetpos, offsetvals,
                         f1, accoffset );
}

/*!
  \param[in] f1 Band 1 flux value to evaluate sigma for
  \returns Value of band 2 color model sigma at f1

  Does validity checks on model state.
*/
double numberCountsDouble::getSigma(double f1) const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (nsigma < 1 ) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  return getSigmaInner(f1);
}

/*!
  \param[in] f1 Band 1 flux value to evaluate offset for
  \returns Value of band 2 color model offset at f1

  Does validity checks on model state.
*/
double numberCountsDouble::getOffset(double f1) const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (noffset < 1) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  return getOffsetInner(f1);
}

/*!
  \param[in] f1 Flux density in band 1
  \param[in] f2 Flux density in band 2
  \returns Number counts at f1, f2

  Assumes model validity already checked
*/
double numberCountsDouble::
getNumberCountsInner(double f1, double f2) const {
  const double normfac = 1.0/sqrt(2 * M_PI);
  if (f1 < knotpos[0] || f1 >= knotpos[nknots - 1] || f2 <= 0.0) 
    return 0.0; //Out of range

  //This is the n_1 bit
  double cnts = exp2(gsl_spline_eval(splinelog, log2(f1), acc));

  //Counts in band 2, Log Normal in f2/f1, multiply them onto n_1
  double if1 = 1.0 / f1;
  double isigma = 1.0 / getSigmaInner(f1);
  double tfac = (log(f2 * if1) - getOffsetInner(f1)) * isigma;
  cnts *= normfac * isigma * exp(-0.5 * tfac * tfac) / f2; //yes, it's 1/f2 here
  return cnts;
}

/*!
  \param[in] f1 Flux density in band 1
  \param[in] f2 Flux density in band 2
  \returns Number counts at f1, f2

  Does validity checks on input.
*/
double numberCountsDouble::getdNdS(double f1, double f2) 
  const {
  if ((nknots < 2) || (nsigma < 1) || (noffset < 1))
    return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f2) || std::isinf(f2)) 
    return std::numeric_limits<double>::quiet_NaN();
  return getNumberCountsInner(f1, f2);
}

/*!
  \param[in] f1 Flux density in band 1
  \returns Band 1 differential number counts at f1

  The input is checked.
 */
double numberCountsDouble::getBand1dNdS(double f1) const {
  if ( (nknots < 2) || (nsigma < 1) || (noffset < 1) )
    return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  if (f1 < knotpos[0] || f1 >= knotpos[nknots-1])
    return 0.0; //Out of range

  return exp2(gsl_spline_eval(splinelog, log2(f1), acc));
}

/*!
  \returns The base number of sources per area
 */
//base_n0 is set by initM1Params
double numberCountsDouble::getBaseN0() const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  return base_n0;
}

//Currently not set -- needs to be added
double numberCountsDouble::getBaseFluxPerArea1() const {return base_flux1;}
double numberCountsDouble::getBaseFluxPerArea2() const {return base_flux2;}
double numberCountsDouble::getBaseFluxSqPerArea1() const {
  return base_fluxsq1; 
}
double numberCountsDouble::getBaseFluxSqPerArea2() const {
  return base_fluxsq2; 
}

dblpair numberCountsDouble::getMaxFluxEstimate() const {
  if (!isValid())
    throw pofdExcept("numberCountsDouble", "getMaxFluxEstimate",
                     "Invalid model", 1);
  //Easy in band one. 
  double mflux1 = knotpos[nknots-1];

  //In band 2 take some number of sigma out (in real space) using the
  // mean/var relationships for a log normal
  double sig = getSigmaInner(mflux1);
  double off = getOffsetInner(mflux1);
  double var = sig*sig;
  double mnratio = exp(off + 0.5 * var);
  double varratio = exp(2*off + var) * (exp(var) - 1.0);
  double mflux2 = mflux1 * (mnratio + 3.0 * sqrt(varratio));
  return std::make_pair(mflux1, mflux2);
}


/*!
  Compute
  \f[
    \int dS_1 \int dS_2 \, S_1^{\alpha} S_2^{\beta} \frac{dN}{dS_1 dS_2} =
      \int dS_1 \, S_1^{\alpha+\beta} n_1\left(S_1\right) \exp \left[
       \beta \mu \left(S_1\right) + \frac{1}{2} \beta^2 \sigma^2
       \left( S_1 \right) \right]
  \f]
  where \f$n_1\left(S_1\right)\f$ is the number counts in band 1.
  This simplification takes advantage of the rather special form
  of the number counts model, which was basically designed to keep
  this simple and avoid actually having to integrate over \f$S_2\f$.  

  Note that to compute the average flux of each object, you have to 
  divide by the total number (i.e., with \f$\alpha = \beta = 0\f$).
  If you want the mean flux per area, however, you don't divide.

  \param[in] alpha   Power of flux in band 1
  \param[in] beta    Power of flux in band 2
  \returns Integral

  The function evaluation is done in evalPowfNDoubleLogNormal, 
  which is the thing that is passed to the GSL integration routines.
 */
double numberCountsDouble::powerInt(double alpha, double beta) const {
  if (nknots < 2) return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  double result, error;
  void *params;
  gsl_function F;
  double minknot = knotpos[0];
  double maxknot = knotpos[nknots-1];
  unsigned int noff = noffset;
  unsigned int nsig = nsigma;
  
  //What we actually pass to evalPowfNDoubleLogNormal is
  // power = alpha+beta
  // const1 = beta
  // const2 = 1/2 beta^2
  //So it evaluates 
  // S_1^power1 n_1(S_1) exp( const1*mu(S_1) + const2*sigma^2(S_1) )
  double power = alpha + beta;
  double const1 = beta;
  double const2 = 0.5 * beta * beta;

  //There are a -ton- of other things to set though, so that
  // evalfN knows what to do in detail (minima, maxima, etc.)

  //Stuff we always need
  varr[0] = static_cast<void*>(&power);
  varr[1] = static_cast<void*>(&const1);
  varr[2] = static_cast<void*>(&const2);
  varr[3] = static_cast<void*>(splinelog);
  varr[4] = static_cast<void*>(acc);
  varr[5] = static_cast<void*>(&minknot);
  varr[6] = static_cast<void*>(&maxknot);
  if (const1 != 0) {
    //We will be evaluating the mu term, so need to pass these
    varr[7]  = static_cast<void*>(offsetinterp);
    varr[8]  = static_cast<void*>(accoffset);
    varr[9]  = static_cast<void*>(&noff);
    varr[10] = static_cast<void*>(offsetpos);
    varr[11] = static_cast<void*>(offsetvals);
  }
  if (const2 != 0) {
    //We will be evaluating the sigma term, so need to pass these
    varr[12] = static_cast<void*>(sigmainterp);
    varr[13] = static_cast<void*>(accsigma);
    varr[14] = static_cast<void*>(&nsig);
    varr[15] = static_cast<void*>(sigmapos);
    varr[16] = static_cast<void*>(sigmavals);
  }

  params = static_cast<void*>(varr);

  F.function = &evalPowfNDoubleLogNormal;
  F.params = params;

  gsl_integration_qag(&F, minknot, maxknot, 0, 1e-5, 1000,
                      GSL_INTEG_GAUSS41, gsl_work, &result, &error); 
  return result;
}

/*!
  \param[in] x1   Source response, band 1
  \param[in] x2   Source response, band 2
  \param[in] bm   The histogrammed inverse beam

  \returns R(x1, x2) computed for the base input model.
*/
double numberCountsDouble::getR(double x1, double x2, 
                                const doublebeamHist& bm) const {

  if (!isValid())
    throw pofdExcept("numberCountsDouble", "getR", 
                     "Invalid model", 1);
  if (!bm.hasData())
    throw pofdExcept("numberCountsDouble", "getR", 
                     "Beam histogram not filled", 2);
  if (!bm.isInverse())
    throw pofdExcept("numberCountsDouble", "getR", 
                     "Beam histogram not inverse", 3);

  //Do actual R computation
  unsigned int curr_n;
  double ieta1, ieta2, retval, ax1, ax2, cts;
  const unsigned int* wtptr;
  const double* ibmptr1;
  const double* ibmptr2;
  retval = 0.0;
  unsigned int sgn = signComp(x1, x2);
  curr_n = bm.getN(sgn);
  if (curr_n == 0) return 0;  // None of that sign component, so there
  wtptr = bm.getWt(sgn);
  ibmptr1 = bm.getBm1(sgn);
  ibmptr2 = bm.getBm2(sgn);
  ax1 = fabs(x1);
  ax2 = fabs(x2);
  for (unsigned int i = 0; i < curr_n; ++i) {
    ieta1 = ibmptr1[i];
    ieta2 = ibmptr2[i];
    cts = getNumberCountsInner(ax1 * ieta1, ax2 * ieta2);
    if (cts > 0) retval += wtptr[i] * ieta1 * ieta2 * cts;
  } 

  double prefac = bm.getPixsize() / 3600.0;
  return prefac * prefac * retval;
}

/*!
  \param[in] n1   Number of source responses, band 1
  \param[in] x1   Source response band 1, length n1
  \param[in] n2   Number of source responses, band 2
  \param[in] x2   Source response band 2, length n2
  \param[in] bm   The histogrammed inverse beam
  \param[out] R   R value, preallocated by caller of dimension n1*n2.  
*/
void numberCountsDouble::getR(unsigned int n1, const double* const x1,
                              unsigned int n2, const double* const x2,
                              const doublebeamHist& bm, double* const R) const {
  
  if (!isValid())
    throw pofdExcept("numberCountsDouble", "getR", 
                     "Invalid model", 1);
  if (!bm.hasData())
    throw pofdExcept("numberCountsDouble", "getR", 
                     "Beam histogram not filled", 2);
  if (!bm.isInverse())
    throw pofdExcept("numberCountsDouble", "getR", 
                     "Beam histogram not inverse", 3);

  //Do R computation
  // It is possible to do this computation much more efficiently
  // by tabulating and saving information.  However, for this code,
  // where this will be a very sub-dominant cost because we re-use it, we go
  // for brute simplicity
  std::memset(R, 0, n1 * n2 * sizeof(double));

  // Load beam sign components into local arrays to avoid
  // inner loop function call overhead
  unsigned int nsgn[4];
  for (unsigned int i = 0; i < 4; ++i) nsgn[i] = bm.getN(i);
  const unsigned int *wtptr4[4];
  for (unsigned int i = 0; i < 4; ++i) 
    if (nsgn[i] > 0) wtptr4[i] = bm.getWt(i); else wtptr4[i] = nullptr;
  const double *ibmptr14[4];
  for (unsigned int i = 0; i < 4; ++i) 
    if (nsgn[i] > 0) ibmptr14[i] = bm.getBm1(i); else ibmptr14[i] = nullptr;
  const double *ibmptr24[4];
  for (unsigned int i = 0; i < 4; ++i) 
    if (nsgn[i] > 0) ibmptr24[i] = bm.getBm2(i); else ibmptr24[i] = nullptr;

  unsigned int sgn, sgn1, curr_n;
  double ieta1, ieta2, cx1, cx2, ax1, ax2, Rsum, cts;
  const unsigned int* wtptr;
  const double* ibmptr1;
  const double* ibmptr2;
  double *rowptr; // Points into rows of R
  bool hasNegX1 = (nsgn[2] > 0) || (nsgn[3] > 0);
  for (unsigned int j = 0; j < n1; ++j) {
    cx1 = x1[j];
    if ((cx1 <= 0) && (!hasNegX1)) continue; //R will be zero for this x1
    rowptr = R + j * n2;
    ax1 = fabs(cx1);
    sgn1 = cx1 >= 0 ? 0 : 2; // Beam sign component when combined below
    for (unsigned int k = 0; k < n2; ++k) {
      cx2 = x2[k];
      sgn = sgn1 + (cx2 >= 0 ? 0 : 1); // Beam sign component
      curr_n = nsgn[sgn]; // Number of beam elements for matching component
      if (curr_n > 0) {
        ax2 = fabs(cx2);
        wtptr = wtptr4[sgn];
        ibmptr1 = ibmptr14[sgn];
        ibmptr2 = ibmptr24[sgn];
        Rsum = 0.0;
        for (unsigned int i = 0; i < curr_n; ++i) {
          ieta1 = ibmptr1[i];
          ieta2 = ibmptr2[i];
          cts = getNumberCountsInner(ax1 * ieta1, ax2 * ieta2);
          if (cts > 0) Rsum += wtptr[i] * ieta1 * ieta2 * cts;
        }
        rowptr[k] += Rsum;
      }
    }
  }

  double prefac = bm.getPixsize() / 3600.0;
  prefac = prefac * prefac;
  for (unsigned int i = 0; i < n1 * n2; ++i) R[i] *= prefac;
}

/*!
  \param[in] udev Uniform deviate [0,1).  
  \param[in] gdev Gaussian deviate with mean 0 and variance 1.
  \returns A pair of fluxes, one for each band, drawn from the model

  udev is not checked for validity in the interests of speed, 
  so screwy things will happen if you provide an invalid one.
  The model is also not checked for validity.
*/
dblpair numberCountsDouble::genSource(double udev, double gdev) const {

  double f1, f2of1;

  //We first generate a flux from band 1.  This is easy because
  // we pre-tabulated what we needed
  f1 = exp2(gsl_interp_eval(gen_interp, gen_interp_cumsum,
                            gen_interp_flux, udev, gen_interp_acc));

  //Now we need the flux in band 2.  This is also easy because
  // a log normal distribution is just one where the log follows
  // a gaussian distribution -- that is, we generate a Gaussian
  // variable with the appropriate mean and sigma, then exponentiate it.
  //Fortunately, one of the inputs is supposed to be a Gaussian deviate
  //Keep in mind that the log normal distribution is in f2/f1, not f2 alone.
  f2of1 = exp(getSigmaInner(f1) * gdev + getOffsetInner(f1));

  return std::make_pair(f1, f2of1 * f1);
}


bool numberCountsDouble::writeToStream(std::ostream& os) const {
  os << "Model parameters: " << std::endl;
  os << " " << std::left << std::setw(13) << "#Flux knot" << "  "
     << std::setw(13) << "Knot value" << std::endl;
  //Convert to log10 for output
  for (unsigned int i = 0; i < nknots; ++i)
    os << " " << std::left << std::setw(13) << knotpos[i] << "  "
       << std::setw(13) << pofd_coverage::ilogfac * logknotvals[i] << std::endl; 
  
  os << " " << std::left << std::setw(13) << "#Sigma knot" << "  "
     << std::setw(13) << "Sigma value" << std::endl;
  for (unsigned int i = 0; i < nsigma; ++i)
    os << " " << std::left << std::setw(13) << sigmapos[i] << "  "
       << std::setw(13) << sigmavals[i] << std::endl; 
  
  os << " " << std::left << std::setw(13) << "#Offset knot" << "  "
     << std::setw(13) << "Offset value" << std::endl;
  for (unsigned int i = 0; i < noffset; ++i)
    os << " " << std::left << std::setw(13) << offsetpos[i] << "  "
       << std::setw(13) << offsetvals[i] << std::endl; 
  return true;
}


/*!
  \param[in] obj_id HDF5 object ID to write to
*/
void numberCountsDouble::writeToHDF5Handle(hid_t obj_id) const {
  if (H5Iget_ref(obj_id) < 0)
    throw pofdExcept("numberCountsDouble", "writeToHDF5Handle",
                     "Given non-open objid to write to", 1);

  if (!isValid())
    throw pofdExcept("numberCountsDouble", "writeToHDF5Handle",
                     "Asked to write invalid model", 2);

  hsize_t adims;
  hid_t mems_id, att_id, dat_id;

  // Single item attributes
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);

  const char modeltype[] = "numberCountsDouble";
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, strlen(modeltype)); 
  att_id = H5Acreate1(obj_id, "ModelType", datatype,
                      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, datatype, modeltype);
  H5Aclose(att_id);

  att_id = H5Acreate2(obj_id, "BaseN0", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &base_n0);
  H5Aclose(att_id);  

  att_id = H5Acreate2(obj_id, "NKnots", H5T_NATIVE_UINT,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nknots);
  H5Aclose(att_id);  

  att_id = H5Acreate2(obj_id, "NSigma", H5T_NATIVE_UINT,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nsigma);
  H5Aclose(att_id);  

  att_id = H5Acreate2(obj_id, "NOffset", H5T_NATIVE_UINT,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &noffset);
  H5Aclose(att_id);  

  H5Sclose(mems_id);
  
  // Knot positions and values as data
  adims = nknots;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  dat_id = H5Dcreate2(obj_id, "KnotPositions", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
           H5P_DEFAULT, knotpos);
  H5Dclose(dat_id);

  // Copy over knot values to temp variable to make log10
  double *ktmp = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    ktmp[i] = getLog10KnotValue(i);
  dat_id = H5Dcreate2(obj_id, "Log10KnotValues", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
           H5P_DEFAULT, ktmp);
  H5Dclose(dat_id);

  delete[] ktmp;
  H5Sclose(mems_id);

  //Sigma knots
  adims = nsigma;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  dat_id = H5Dcreate2(obj_id, "SigmaKnotPositions", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
           H5P_DEFAULT, sigmapos);
  H5Dclose(dat_id);

  dat_id = H5Dcreate2(obj_id, "SigmaKnotValues", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
           H5P_DEFAULT, sigmavals);
  H5Dclose(dat_id);
  H5Sclose(mems_id);

  //Offset
  adims = noffset;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  dat_id = H5Dcreate2(obj_id, "OffsetKnotPositions", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
           H5P_DEFAULT, offsetpos);
  H5Dclose(dat_id);

  dat_id = H5Dcreate2(obj_id, "OffsetKnotValues", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
           H5P_DEFAULT, offsetvals);
  H5Dclose(dat_id);
  H5Sclose(mems_id);

}

/*!
  Internal function for use in model integrations.
 */
static double evalPowfNDoubleLogNormal(double s1, void* params) {
  //Model is: ( s1^power * exp( const1*mu + const2*sigma^2 ) ) * n1(f)
  // where f is the flux in the first band and n1 is the number
  // counts in band 1 -- see powerInt.  s1 is the flux in band 1.
  //Params are:
  // params[0]  power
  // params[1]  const1
  // params[2]  const2
  // params[3]  knot interpolant (log)
  // params[4]  knot accelerator
  // params[5]  minknot
  // params[6]  maxnot
  // params[7]  offset interpolant (log)
  // params[8]  offset accelerator
  // params[9]  noffsets
  // params[10] offset positions
  // params[11] offset values
  // params[12] sigma interpolant (log)
  // params[13] sigma accelerator
  // params[14] nsigmas
  // params[15] sigma positions
  // params[16] sigma values
  //But this really has to be an array of pointers to void to work
  void** vptr = static_cast<void**>(params);

  //First get min/max knot in band 1 for quick return if we are outside that
  double minknot = *static_cast<double*>(vptr[5]);
  double maxknot = *static_cast<double*>(vptr[6]);
  if (s1 < minknot || s1 >= maxknot) return 0.0;

  //Get coeffs
  double power  = *static_cast<double*>(vptr[0]);
  double const1 = *static_cast<double*>(vptr[1]);
  double const2 = *static_cast<double*>(vptr[2]);

  //Construct thing we multiply n1 counts by
  double prefac;
  //Construct s1^power part
  if (fabs(power) < 1e-6)
    prefac = 1.0;
  else if (fabs(power - 1.0) < 1e-6) 
    prefac = s1; 
  else if (fabs(power - 2.0) < 1e-6) 
    prefac = s1 * s1;
  else prefac = pow(s1, power);

  //Now exponential part
  if (const1 != 0 || const2 != 0) {
    double expbit; //Part to exponentiate
    if (const1 != 0) {
      //Evaluate offset at s1
      unsigned int noffsets = *static_cast<unsigned int*>(vptr[9]);
      double *offsetpos = static_cast<double*>(vptr[10]);
      double *offsetval = static_cast<double*>(vptr[11]);
      double mu;
      if (noffsets == 1) mu = offsetval[0];
      else if (s1 <= offsetpos[0]) mu = offsetval[0];
      else if (s1 >= offsetpos[noffsets-1]) mu = offsetval[noffsets-1];
      else {
        gsl_interp* ospl = static_cast<gsl_interp*>(vptr[7]);
        gsl_interp_accel* oacc = static_cast<gsl_interp_accel*>(vptr[8]);
        mu = gsl_interp_eval(ospl, offsetpos, offsetval, s1, oacc);
      }
      expbit = const1 * mu;
    } else expbit = 0.0;

    if (const2 != 0) {
      //Get sigma bit -> const2 * sigma^2
      double sigma;
      unsigned int nsigmas = *static_cast<unsigned int*>(vptr[14]);
      double *sigmapos = static_cast<double*>(vptr[15]);
      double *sigmaval = static_cast<double*>(vptr[16]);
      if (nsigmas == 1) sigma = sigmaval[0];
      else if (s1 <= sigmapos[0]) sigma = sigmaval[0];
      else if (s1 >= sigmapos[nsigmas-1]) sigma = sigmaval[nsigmas-1];
      else {
        gsl_interp* sspl = static_cast<gsl_interp*>(vptr[12]);
        gsl_interp_accel* sacc = static_cast<gsl_interp_accel*>(vptr[13]);
        sigma = gsl_interp_eval(sspl, sigmapos, sigmaval, s1, sacc);
      }
      expbit += const2 * sigma * sigma;
    } 

    prefac *= exp(expbit);

  } //Otherwise exp(expbit) is just 1

  //Now multiply in n(band1)
  //Now multiply in n(band1)
  gsl_spline* spl = static_cast<gsl_spline*>(vptr[3]);
  gsl_interp_accel* acc = static_cast<gsl_interp_accel*>(vptr[4]);
  double splval = exp2(gsl_spline_eval(spl,log2(s1),acc));
  return prefac * splval;
}

std::ostream& operator<<(std::ostream& os, const numberCountsDouble& b) {
  b.writeToStream(os);
  return os;
}
