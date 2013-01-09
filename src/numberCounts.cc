#include<limits>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

#include<pofdExcept.h>
#include<global_settings.h>
#include<numberCounts.h>
#include<utility.h>

const double numberCounts::ftol = 1e-4;

numberCounts::numberCounts(const std::string& modelfile) {

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
    if (line[0] == '#') continue; //Skip comments
    
    //Parse into words, stipping spaces
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
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
  knotvals = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    knotvals[i] = kv[i];

  a = gamma = omg = iomg = fk = powarr = NULL;
  gamone = NULL;
  initMParams();

  //Don't set histogrammed beam currently
  nwrk = 0;
  wrk_wts = NULL;
  wrk_bm = NULL;

}


numberCounts::numberCounts(unsigned int NKNOTS, const double* const KNOTPOS,
			   const double* const KNOTVAL) {

  //Input checks
  if (NKNOTS < 2)
    throw pofdExcept("numberCounts","numberCounts",
		     "nknots must be >= 2",1);

  for (unsigned int i = 0; i < NKNOTS; ++i)
    if (KNOTPOS[i] <= 0.0)
      throw pofdExcept("numberCounts","numberCounts",
		       "knot positions must all be positive",2);

  for (unsigned int i = 1; i < NKNOTS; ++i)
    if (KNOTPOS[i] <= KNOTPOS[i-1])
      throw pofdExcept("numberCounts","numberCounts",
		       "knot positions must be monotonically increasing",3);

  for (unsigned int i = 0; i < NKNOTS; ++i)
    if (KNOTVAL[i] <= 0.0)
      throw pofdExcept("numberCounts","numberCounts",
		       "knot values must all be positive",4);

  //Now construct internal model parameters
  nknots = NKNOTS;
  knotpos = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotpos[i] = KNOTPOS[i];
  knotvals = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotvals[i] = KNOTVAL[i];

  a = gamma = omg = iomg = fk = powarr = NULL;
  gamone = NULL;
  initMParams();

  nwrk = 0;
  wrk_wts = NULL;
  wrk_bm = NULL;
}

void numberCounts::initMParams() {
  
  //Set up gamma and a
  if (gamma != NULL) delete[] gamma;
  gamma = new double[nknots-1];
  for (unsigned int i = 0; i < nknots-1; ++i)
    gamma[i] = - log(knotvals[i+1]/knotvals[i]) / 
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
    omg[i] = val = gamma[i] - 1.0;
    if (fabs(val) < ftol) {
      gamone[i] = true; 
      iomg[i] = 1.0 / val;
    } else gamone[i] = false;
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
  base_n0 = fk[nknots - 2];

  //Compute the mean flux per area for the base model
  base_meanflux = 0.0;
  double tmg;
  for (unsigned int i = 0; i < nknots-1; ++i) {
    tmg = 2.0 - gamma[i];
    if (fabs(tmg) < ftol)
      base_meanflux += a[i] * log(knotpos[i+1]/knotpos[i]);
    else
      base_meanflux += 
	a[i] * (pow(knotpos[i+1], tmg) - pow(knotpos[i], tmg)) / tmg;
  }

  //And the mean flux^2 per area
  base_meanfluxsq = 0.0;
  for (unsigned int i = 0; i < nknots-1; ++i) {
    tmg = 3.0 - gamma[i];
    if (fabs(tmg) < ftol)
      base_meanfluxsq += a[i] * log(knotpos[i+1]/knotpos[i]);
    else
      base_meanfluxsq += 
	a[i] * (pow(knotpos[i+1], tmg) - pow(knotpos[i], tmg)) / tmg;
  }
}

numberCounts::~numberCounts() {
  delete[] knotpos;
  delete[] knotvals;
  delete[] gamma;
  delete[] a;
  delete[] omg;
  delete[] iomg;
  delete[] gamone;
  delete[] fk;
  delete[] powarr;
  if (wrk_wts != NULL) delete[] wrk_wts;
  if (wrk_bm != NULL) delete[] wrk_bm;
}

bool numberCounts::isValid() const {
  if (std::isnan(base_n0)) return false;
  if (std::isinf(base_n0)) return false;
  if (base_n0 <= 0) return false;
  return true;
}

/*!
  \returns The base number of sources per area
 */
double numberCounts::getBaseN0() const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  return base_n0;
}

/*!
  \returns Mean flux per square degree for base model
 */
double numberCounts::getMeanFluxPerArea() const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  return base_meanflux;
}

/*!
  \returns Mean flux^2 per square degree for base model
 */
double numberCounts::getMeanFluxSqPerArea() const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  return base_meanfluxsq;
}

/*!
  \param[in] flux Flux density number counts are desired at

  \returns Differential number counts at flux density flux for base model
 */

double numberCounts::getdNdS(double flux) const {
  if (flux < knotpos[0]) return 0.0;
  if (flux >= knotpos[nknots-1]) return 0.0;
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  unsigned int loc;
  loc = utility::binary_search_lte(flux, knotpos, nknots);
  return a[loc] * pow(flux, -gamma[loc]);
}

/*!
  \param[in] x The value R is desired for
  \param[in] bm The beam
  \param[in] pixsize The pixel size, in arcseconds
  \param[in] nfwhm The number of beam fwhm to use in the computation
  \param[in] nbins The number of bins in the beam histogramming

  \returns R(x) computed for the base input model.
*/
double numberCounts::getR(double x, const beam& bm,
			  double pixsize, double nfwhm,
			  unsigned int nbins) const {

  //This could be done much more efficiently (see the old pofd_mcmc
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

  double s_min = knotpos[0];
  double s_max = knotpos[nknots-1];
  if (x >= s_max) return 0.0;
  if (x <= 0.0) return 0.0;

  //Get the inverse binned beam histogram
  if (nwrk < nbins) {
    //Must expand
    if (wrk_wts == NULL) delete[] wrk_wts;
    if (wrk_bm == NULL) delete[] wrk_bm;
    wrk_wts = new unsigned int[nbins];
    wrk_bm = new double[nbins];
    nwrk = nbins;
  }
  unsigned int npix = static_cast<unsigned int>(nfwhm * bm.getFWHM()/pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;
  unsigned int nnonzero;
  bm.getBeamHist(npix, pixsize, nbins, nnonzero, wrk_wts, wrk_bm, true);

  //And now the actual computation
  double cval, cR, R, ibm;
  unsigned int loc;
  R = 0.0;
  for (unsigned int i = 0; i < nnonzero; ++i) {
    ibm = wrk_bm[i]; //1 / eta
    cval = x * ibm;
    if (cval < s_min) continue;
    if (cval >= s_max) continue;
    loc = utility::binary_search_lte(cval, knotpos, nknots);
    cR = a[loc] * pow(cval, -gamma[loc]);
    R += wrk_wts[i] * cR * ibm;
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

  R is computed for the base model		 
*/

void numberCounts::getR(unsigned int n, double minflux,
			double maxflux, const beam& bm, 
			double pixsize, double nfwhm,
			unsigned int nbins, double* R) const {

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

  //Get the inverse binned beam histogram
  if (nwrk < nbins) {
    //Must expand
    if (wrk_wts == NULL) delete[] wrk_wts;
    if (wrk_bm == NULL) delete[] wrk_bm;
    wrk_wts = new unsigned int[nbins];
    wrk_bm = new double[nbins];
    nwrk = nbins;
  }
  unsigned int npix = static_cast<unsigned int>(nfwhm * bm.getFWHM()/pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;
  unsigned int nnonzero;
  bm.getBeamHist(npix, pixsize, nbins, nnonzero, wrk_wts, wrk_bm, true);

  double prefac = pixsize / 3600.0;
  prefac = prefac * prefac;

  //And now the actual computation.  We loop over each flux
  double cflux, cval, cR, ibm, workR;
  unsigned int loc;
  for (unsigned int i = 0; i < n; ++i) {
    cflux = minflux + static_cast<double>(i) * dflux;
    if (cflux <= 0.0 || cflux >= s_max) {
      //R is zero for this x
      R[i] = 0.0;
      continue;
    }
    workR = 0.0;
    for (unsigned int j = 0; j < nnonzero; ++j) {
      ibm = wrk_bm[j]; //1 / eta
      cval = cflux * ibm;
      if (cval < s_min) continue;
      if (cval >= s_max) continue;
      loc = utility::binary_search_lte(cval, knotpos, nknots);
      cR = a[loc] * pow(cval, -gamma[loc]);
      workR += wrk_wts[j] * cR * ibm;
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
  // in flux density
  //  We want F[k-1] <= udev < F[k], so that S_k <= A < S_{k+1}
  // where A is the value we are trying to generate
  double prod = udev * base_n0;

  if (prod < fk[0]) {
    //Between the 0th and 1st knot, special case
    if (gamone[0])
      return knotpos[0] * exp(prod / a[0]);
    else
      return pow(omg[0] * prod / a[0] + powarr[0], iomg[0]);
  }

  unsigned int km1 = utility::binary_search_lte(prod, fk, nknots-1);
  
  double delta = prod - fk[km1];
  unsigned int k = km1 + 1;
  if (fabs(delta / prod) < ftol) {
    //Close enough! A = S_k
    return knotpos[k];
  } 
  
  delta /= a[k];
  if (gamone[k]) 
    return knotpos[k] * exp(delta);
  else 
    return pow(omg[k] * delta + powarr[k], iomg[k]);
}

bool numberCounts::writeToStream(std::ostream& os) const {
  os << "Base model parameters: " << std::endl;
  os << " " << std::left << std::setw(13) << "#Flux knot" << "  "
     << std::setw(13) << "Log10 Knot value" << std::endl;
  for (unsigned int i = 0; i < nknots; ++i)
    os << " " << std::left << std::setw(13) << knotpos[i] << "  "
       << std::setw(13) << log10( a[i] * pow(knotpos[i], -gamma[i])) 
       << std::endl; 
  return true;
}


std::ostream& operator<<(std::ostream& os, const numberCounts& b) {
  b.writeToStream(os);
  return os;
}
