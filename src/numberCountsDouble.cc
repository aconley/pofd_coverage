#include<cmath>
#include<iomanip>
#include<sstream> 
#include<fstream>

#include<utility.h>
#include<numberCountsDouble.h>
#include<global_settings.h>
#include<pofdExcept.h>

//Function to pass to GSL integrator
/*! \brief Evaluates flux1^power1 * exp(const1*mu + const2*sigma^2) dN/dS1 */
static double evalPowfNDoubleLogNormal(double,void*); 

const unsigned int numberCountsDouble::nvarr = 17;
const double numberCountsDouble::ftol = 1e-4;

/*!
  \param[in] modelfile Name of file to read base model from
*/
numberCountsDouble::numberCountsDouble(const std::string& modelfile) {
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
    throw pofdExcept("numberCountsDouble","numberCountsDouble",
		       errmsg.str(),1);
  }

  //Read in number of knots in band1, sigmas, offsets
  bool has_nknots = false;
  while (!initfs.eof()) {
    std::getline(initfs,line);
    if (line[0] == '#') continue; //Comment
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
    throw pofdExcept("numberCountsDouble","numberCountsDouble",
		     "Unable to find number of knots line",2);
  }
  if ( nk < 2 ) {
    initfs.close();
    throw pofdExcept("numberCountsDouble","numberCountsDouble",
		       "Need at least 2 band 1 knots",3);
  }
  if ( ns < 1 ) {
    initfs.close();
    throw pofdExcept("numberCountsDouble","numberCountsDouble",
		       "Need at least one sigma color model knot",4);

  }
  if ( no < 1 ) {
    initfs.close();
    throw pofdExcept("numberCountsDouble","numberCountsDouble",
		       "Need at least one offset color model knot",5);
  }
  
  //Read in values
  while (!initfs.eof()) {
    std::getline(initfs,line);
    if (line[0] == '#') continue; //Comment
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
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
    errstr << "Expected " << ntot << " values, got: " 
	   << wvec1.size();
    throw pofdExcept("numberCountsDouble", "numberCountsDouble",
		       errstr.str(), 6);
  }

  //Set up band 1
  nknots = nk;
  knotpos = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotpos[i] = wvec1[i];

  //These are read in as log_10, so we must convert to log_e
  logknotvals = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = wvec2[i] * pofd_coverage::logfac;

  //Get non-log values
  knotvals = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotvals[i] = exp(logknotvals[i]);

  //Set up internal convenience variables for band 1
  a = gamma = omg = iomg = fk = powarr = NULL;
  gamone = NULL;
  initM1Params();

  //Now set up band 2 info
  // First, the sigma spline
  nsigma = ns;
  sigmapos = new double[nsigma];
  for (unsigned int i = 0; i < nsigma; ++i) sigmapos[i] = wvec1[i+nk];
  sigmavals = new double[nsigma];
  for (unsigned int i = 0; i < nsigma; ++i) sigmavals[i] = wvec2[i+nk];
  if (nsigma > 2)
    sigmainterp = gsl_interp_alloc(gsl_interp_cspline,
				   static_cast<size_t>(nsigma));
  else
    sigmainterp = gsl_interp_alloc(gsl_interp_linear,
				   static_cast<size_t>(nsigma));
  if (nsigma > 1)
    gsl_interp_init(sigmainterp, sigmapos, sigmavals,
		    static_cast<size_t>(nsigma));
  accsigma = gsl_interp_accel_alloc();

  // Offset spline
  noffset = no;
  offsetpos = new double[noffset];
  for (unsigned int i = 0; i < noffset; ++i) offsetpos[i] = wvec1[i+nk+ns];
  offsetvals = new double[noffset];
  for (unsigned int i = 0; i < noffset; ++i) offsetvals[i] = wvec2[i+nk+ns];
  if (noffset > 2)
    offsetinterp = gsl_interp_alloc(gsl_interp_cspline,
				    static_cast<size_t>(noffset));
  else
    offsetinterp = gsl_interp_alloc(gsl_interp_linear,
				    static_cast<size_t>(noffset));
  if (noffset > 1)
    gsl_interp_init(offsetinterp, offsetpos, offsetvals,
		    static_cast<size_t>(noffset));
  accoffset = gsl_interp_accel_alloc();
  
  //Make sure what we read makes sense
  if (!isValidLoaded())
    throw pofdExcept("numberCountsDouble", "numberCountsDouble",
		     "Invalid base model parameters",7);

  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[nvarr];

  //Compute band 2 mean flux and mean flux^2.  Band 1 was taken
  // care of in initM1Params
  base_meanflux2 = powerInt(0.0, 1.0) / base_n0;
  base_meanfluxsq2 = powerInt(0.0, 2.0) / base_n0;

  nbm = 0;
  bm_wts = NULL;
  inv_bm1 = inv_bm2 = NULL;
}

numberCountsDouble::~numberCountsDouble() {
  gsl_interp_accel_free(accsigma);
  gsl_interp_accel_free(accoffset);
  if (sigmainterp != NULL) gsl_interp_free(sigmainterp);
  if (offsetinterp != NULL) gsl_interp_free(offsetinterp);
  gsl_integration_workspace_free(gsl_work);
  delete[] varr;

  delete[] knotpos;
  delete[] knotvals;
  delete[] logknotvals;
  delete[] gamma;
  delete[] a;
  delete[] omg;
  delete[] iomg;
  delete[] gamone;
  delete[] fk;
  delete[] powarr;

  if (bm_wts != NULL) delete[] bm_wts;
  if (inv_bm1 != NULL) delete[] inv_bm1;
  if (inv_bm2 != NULL) delete[] inv_bm2;
}


void numberCountsDouble::initM1Params() {
  
  //Set up gamma and a
  if (gamma != NULL) delete[] gamma;
  gamma = new double[nknots-1];
  for (unsigned int i = 0; i < nknots-1; ++i)
    gamma[i] = - (logknotvals[i+1] - logknotvals[i]) / 
      log(knotpos[i+1]/knotpos[i]);

  if (a != NULL) delete[] a;
  a = new double[nknots-1];
  for (unsigned int i = 0; i < nknots-1; ++i)
    a[i] = knotvals[i] * pow(knotpos[i], gamma[i]);

  if (gamone != NULL) delete[] gamone;
  gamone = new bool[nknots-1];
  if (omg != NULL) delete[] omg;
  omg = new double[nknots-1];
  if (iomg != NULL) delete[] iomg;
  iomg = new double[nknots-1];
  double val;
  for (unsigned int i = 0; i < nknots-1; ++i) {
    omg[i] = val = 1.0 - gamma[i];
    if (fabs(val) < ftol) gamone[i] = true; else {
      iomg[i] = 1.0 / val;
      gamone[i] = false;
    }
  }

  if (fk != NULL) delete[] fk;
  if (powarr != NULL) delete[] powarr;
  fk = new double[nknots-1];
  powarr = new double[nknots-1];

  //Compute the total number of sources and the fk, powarr arrays
  double m, tmp;
  for (unsigned int i = 0; i < nknots-1; ++i) {
    if (gamone[i]) {
      tmp = log(knotpos[i]);
      m = a[i] * (log(knotpos[i+1]) - tmp);
    } else {
      tmp = pow(knotpos[i], omg[i]);
      m = a[i] * iomg[i] * (pow(knotpos[i+1], omg[i]) - tmp);
    }
    powarr[i] = tmp;
    if (i == 0) fk[i] = m;
    else fk[i] = m + fk[i-1];
  }
  base_n0 = fk[nknots-2];

  //Compute the mean flux per area for the base model, band 1
  base_meanflux1 = 0.0;
  double tmg;
  for (unsigned int i = 0; i < nknots-1; ++i) {
    tmg = 2.0 - gamma[i];
    if (fabs(tmg) < ftol)
      base_meanflux1 += a[i] * log(knotpos[i+1]/knotpos[i]);
    else
      base_meanflux1 += 
	a[i] * (pow(knotpos[i+1], tmg) - pow(knotpos[i], tmg)) / tmg;
  }

  //And the mean flux^2 per area, band 1
  base_meanfluxsq1 = 0.0;
  for (unsigned int i = 0; i < nknots-1; ++i) {
    tmg = 3.0 - gamma[i];
    if (fabs(tmg) < ftol)
      base_meanfluxsq1 += a[i] * log(knotpos[i+1]/knotpos[i]);
    else
      base_meanfluxsq1 += 
	a[i] * (pow(knotpos[i+1], tmg) - pow(knotpos[i], tmg)) / tmg;
  }
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
  if (f1 >= sigmapos[nsigma-1]) return sigmavals[nsigma-1];
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
  if (f1 >= offsetpos[noffset-1]) return offsetvals[noffset-1];
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
  const double normfac = 1.0/sqrt(2*M_PI);
  if (f1 < knotpos[0] || f1 >= knotpos[nknots-1] || f2 <= 0.0) 
    return 0.0; //Out of range

  //This is the n_1 bit
  unsigned int loc;
  loc = utility::binary_search_lte(f1, knotpos, nknots);
  double cnts = a[loc] * pow(f1, -gamma[loc]);

  //Counts in band 2, Log Normal in f2/f1, multiply them onto n_1
  double if1 = 1.0/f1;
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
  if ( (nknots < 2) || (nsigma < 1) || (noffset < 1) )
    return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if ( std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  if ( std::isnan(f2) || std::isinf(f2)) 
    return std::numeric_limits<double>::quiet_NaN();
  return getNumberCountsInner(f1,f2);
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
double numberCountsDouble::getMeanFluxPerArea1() const {return base_meanflux1;}
double numberCountsDouble::getMeanFluxPerArea2() const {return base_meanflux2;}
double numberCountsDouble::getMeanFluxSqPerArea1() const {
  return base_meanfluxsq1; 
}
double numberCountsDouble::getMeanFluxSqPerArea2() const {
  return base_meanfluxsq2; 
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
  unsigned int n = nknots;
  unsigned int noff = noffset;
  unsigned int nsig = nsigma;
  
  //What we actually pass to evalPowfNDoubleLogNormal is
  // power = alpha+beta
  // const1 = beta
  // const2 = 1/2 beta^2
  //So it evaluates 
  // S_1^power1 n_1(S_1) exp( const1*mu(S_1) + const2*sigma^2(S_2) )
  double power = alpha+beta;
  double const1 = beta;
  double const2 = 0.5*beta*beta;

  //There are a -ton- of other things to set though, so that
  // evalfN knows what to do in detail (minima, maxima, etc.)

  //Stuff we always need
  varr[0] = static_cast<void*>(&power);
  varr[1] = static_cast<void*>(&const1);
  varr[2] = static_cast<void*>(&const2);
  varr[3] = static_cast<void*>(&n);
  varr[4] = static_cast<void*>(knotpos);
  varr[5] = static_cast<void*>(a);
  varr[6] = static_cast<void*>(gamma);
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
  \param[in] bm   The beam
  \param[in] pixsize The pixel size, in arcseconds, in both bands
  \param[in] nfwhm The number of beam fwhm to use in the computation
  \param[in] nbins The number of bins in the beam histogramming

  \returns R(x1, x2) computed for the base input model.
*/
double numberCountsDouble::getR(double x1,double x2, const doublebeam& bm,
				double pixsize, double nfwhm,
				unsigned int nbins) const {

  if (pixsize <= 0.0)
    throw pofdExcept("numberCountsDouble", 
		     "getR", "Invalid (non-positive) pixsize", 1);
  if (nfwhm <= 0.0)
    throw pofdExcept("numberCountsDouble", "getR", 
		     "Invalid (non-positive) nfwhm", 2);
  if (nbins == 0)
    throw pofdExcept("numberCountsDouble", "getR", 
		     "Invalid (zero) nbins", 3);
  if (!isValid())
    throw pofdExcept("numberCountsDouble", "getR", 
		     "Invalid model", 4);

  double s1_max = knotpos[nknots-1];
  if (x1 >= s1_max) return 0.0;
  if (x1 <= 0.0 || x2 <= 0.0) return 0.0;

  //Storage for beam pixels
  unsigned int nbm_needed = nbins * nbins;
  if (nbm < nbm_needed) {
    //Must expand
    if (bm_wts == NULL) delete[] bm_wts;
    if (inv_bm1 == NULL) delete[] inv_bm1;
    if (inv_bm2 == NULL) delete[] inv_bm2;
    nbm = nbm_needed;
    bm_wts = new unsigned int[nbm];
    inv_bm1 = new double[nbm];
    inv_bm2 = new double[nbm];
  }
  //Get number of pixels out we will go.  This has to be the same
  // in each beam
  double maxfwhm = bm.getFWHM1();
  if (bm.getFWHM2() > maxfwhm) maxfwhm = bm.getFWHM2();
  unsigned int npix = static_cast<unsigned int>(nfwhm * maxfwhm/pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;

  //Load inverse histogrammed beams into bm_wts, inv_bm?
  unsigned int nnonzero;
  bm.getBeamHist(npix, pixsize, nbins, nnonzero, bm_wts, 
		 inv_bm1, inv_bm2, true);

  //Do actual R computation
  double ieta1, ieta2, retval;
  retval = 0.0;
  for (unsigned int i = 0; i < nnonzero; ++i) {
      ieta1 = inv_bm1[i];
      ieta2 = inv_bm2[i];
      retval += bm_wts[i] * ieta1 * ieta2 * 
	getNumberCountsInner(x1 * ieta1, x2 * ieta2);
    } 

  double prefac = pixsize / 3600.0;
  return prefac * prefac * retval;
}

/*!
  \param[in] n1   Number of source responses, band 1
  \param[in] x1   Source response band 1, length n1
  \param[in] n2   Number of source responses, band 2
  \param[in] x2   Source response band 2, length n2
  \param[in] bm   The beam
  \param[in] pixsize The pixel size, in arcseconds, in both bands
  \param[in] nfwhm The number of beam fwhm to use in the computation
  \param[in] nbins The number of bins in the beam histogramming
  \param[out] R   R value dimension n1*n2.  
*/
void numberCountsDouble::getR(unsigned int n1, const double* const x1,
			      unsigned int n2, const double* const x2,
			      const doublebeam& bm, double pixsize,
			      double nfwhm, unsigned int nbins,
			      double* R) const {
  
  if (pixsize <= 0.0)
    throw pofdExcept("numberCountsDouble", 
		     "getR", "Invalid (non-positive) pixsize", 1);
  if (nfwhm <= 0.0)
    throw pofdExcept("numberCountsDouble", "getR", 
		     "Invalid (non-positive) nfwhm", 2);
  if (nbins == 0)
    throw pofdExcept("numberCountsDouble", "getR", 
		     "Invalid (zero) nbins", 3);
  if (!isValid())
    throw pofdExcept("numberCountsDouble", "getR", 
		     "Invalid model", 4);

  double s1_max = knotpos[nknots-1];

  //Storage for beam pixels
  unsigned int nbm_needed = nbins * nbins;
  if (nbm < nbm_needed) {
    //Must expand
    if (bm_wts == NULL) delete[] bm_wts;
    if (inv_bm1 == NULL) delete[] inv_bm1;
    if (inv_bm2 == NULL) delete[] inv_bm2;
    nbm = nbm_needed;
    bm_wts = new unsigned int[nbm];
    inv_bm1 = new double[nbm];
    inv_bm2 = new double[nbm];
  }
  //Get number of pixels out we will go.  This has to be the same
  // in each beam
  double maxfwhm = bm.getFWHM1();
  if (bm.getFWHM2() > maxfwhm) maxfwhm = bm.getFWHM2();
  unsigned int npix = static_cast<unsigned int>(nfwhm * maxfwhm/pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;

  //Load inverse histogrammed beams into bm_wts, inv_bm?
  unsigned int nnonzero;
  bm.getBeamHist(npix, pixsize, nbins, nnonzero, bm_wts, 
		 inv_bm1, inv_bm2, true);

  //Do R computation
  // It is possible to do this computation much more efficiently
  // by tabulating and saving information.  However, for this code,
  // where this will likely be a very sub-dominant cost, we go
  // for brute simplicity
  double ieta1, ieta2, cx1, cx2, Rsum;
  double *rowptr;
  double prefac = pixsize / 3600.0;
  prefac = prefac * prefac;
  for (unsigned int j = 0; j < n1; ++j) {
    cx1 = x1[j];
    rowptr = R + j * n2;
    if (cx1 <= 0 || cx1 >= s1_max)
      for (unsigned int k = 0; k < n2; ++k)
	rowptr[k] = 0.0;
    else 
      for (unsigned int k = 0; k < n2; ++k) {
	cx2 = x2[k];
	if (cx2 <= 0.0) rowptr[k] = 0.0; else {
	  //x1, x2 are in bounds, must do full sum over beams
	  Rsum = 0.0;
	  for (unsigned int i = 0; i < nnonzero; ++i) {
	    ieta1 = inv_bm1[i];
	    ieta2 = inv_bm2[i];
	    Rsum += bm_wts[i] * ieta1 * ieta2 * 
	      getNumberCountsInner(cx1 * ieta1, cx2 * ieta2);
	  }
	  rowptr[k] = prefac * Rsum;
	}
      }
  }
}

/*!
  \param[in] udev Uniform deviate [0,1).  
  \param[in] gdev Gaussian deviate with mean 0 and variance 1.
  \returns A pair of fluxes, one for each band, drawn from the model

  udev is not checked for validity in the interests of speed, 
  so screwy things will happen if you provide an invalid one.
  The model is also not checked for validity.
*/
std::pair<double, double> 
numberCountsDouble::genSource(double udev, double gdev) const {

  double f1, f2of1;

  //We first generate a flux from band 1.  This is fairly easy
  // for a power law model.
  // We want F[k-1] <= udev < F[k], so that S_k <= A < S_{k+1}
  //  where A is the value we are trying to generate
  double prod = udev * base_n0;
  if (prod < fk[0]) {
    //Between the 0th and 1st knot, special case
    if (gamone[0])
      f1 = knotpos[0] * exp(prod / a[0]);
    else
      f1 = pow(omg[0] * prod / a[0] + powarr[0], iomg[0]);
  } else {
    unsigned int km1 = utility::binary_search_lte(prod, fk, nknots-1);
  
    double delta = prod - fk[km1];
    unsigned int k = km1 + 1;
    if (fabs(delta / prod) < ftol) {
      //Close enough! A = S_k
      f1 = knotpos[k];
    } else {
      delta /= a[k];
      if (gamone[k]) 
	f1 = knotpos[k] * exp(delta);
      else 
        f1 = pow(omg[k] * delta + powarr[k], iomg[k]);
    }
  }

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
  // params[3]  nknots in band 1 model
  // params[4]  knot positions in band 1 model
  // params[5]  a variable in band 1 model
  // params[6]  gamma variable in band 1 model
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
  unsigned int nknots = *static_cast<unsigned int*>(vptr[3]);
  double *knotpos = static_cast<double*>(vptr[4]);
  double minknot = knotpos[0];
  double maxknot = knotpos[nknots-1];

  if (s1 < minknot || s1 >= maxknot) return 0.0;

  //Get coeffs
  double power  = *static_cast<double*>(vptr[0]);
  double const1 = *static_cast<double*>(vptr[1]);
  double const2 = *static_cast<double*>(vptr[2]);

  //Construct thing we multiply n1 counts by
  double prefac;
  //Construct s1^power part
  if (power == 0) 
    prefac = 1.0; 
  else {
    if (fabs(power-1.0) < 1e-6) prefac = s1; 
    else if (fabs(power-2.0) < 1e-6) prefac = s1*s1;
    else prefac = pow(s1,power);
  }


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
	mu = gsl_interp_eval( ospl, offsetpos, offsetval, s1, oacc );
      }
      expbit = const1 * mu;
    } else expbit = 0.0;

    if (const2 != 0) {
      //Get sigma bit -> const2* sigma^2
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
	sigma = gsl_interp_eval( sspl, sigmapos, sigmaval, s1, sacc );
      }
      expbit += const2 * sigma * sigma;
    } 

    prefac *= exp(expbit);

  } //Otherwise exp(expbit) is just 1

  //Now multiply in n(band1)
  unsigned int loc;
  loc = utility::binary_search_lte(s1, knotpos, nknots);
  double *a = static_cast<double*>(vptr[5]);
  double *gamma = static_cast<double*>(vptr[6]);
  double n1cnts = a[loc] * pow(s1, -gamma[loc]);

  return prefac * n1cnts;
}

std::ostream& operator<<(std::ostream& os, const numberCountsDouble& b) {
  b.writeToStream(os);
  return os;
}
