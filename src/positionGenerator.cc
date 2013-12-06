#include<fstream>
#include<sstream>
#include<limits>


#include "../include/positionGenerator.h"
#include "../include/pofdExcept.h"
#include "../include/utility.h"

////////////////// Powerspectrum ////////////////
/*!
  \param[in] filename File to read k, P(k) from
*/
powerSpectrum::powerSpectrum(const std::string& filename) {
  // Just read the file, pass it to init
  std::ifstream ifs(filename.c_str());
  if (!ifs)
    throw pofdExcept("numberCounts","numberCounts",
		     "Unable to open input model file",1);

  //Do the read into temporary vectors, then copy
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double currval;
  std::vector<double> kin, pkin;
  const unsigned int nreq = 2;
  while (!ifs.eof()) {
    std::getline(ifs, line);
    if (line[0] == '#') continue; //Skip comments
    
    //Parse into words, stipping spaces
    utility::stringwords(line, words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
    if (words.size() < nreq) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currval;
    kin.push_back(currval);
    str.str(words[1]); str.clear(); str >> currval;
    pkin.push_back(currval);
  }
  ifs.close();

  init(kin, pkin);
}

/*!
  \param[in] k k vector in inverse arcmin.  Assumed sorted
  \param[in] pk P(k)
*/
powerSpectrum::powerSpectrum(const std::vector<double>& k,
			     const std::vector<double>& pk) {
  init(k, pk);
}

powerSpectrum::~powerSpectrum() {
  delete[] logk;
  delete[] logpk;
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline_logpk);
}


/*!
  \param[in] k k vector in inverse arcmin.  Assumed sorted
  \param[in] pk P(k)
*/
void powerSpectrum::init(const std::vector<double>& k,
			 const std::vector<double>& pk) {
  nk = k.size();
  if (nk == 0)
    throw pofdExcept("powerSpectrum", "powerSpectrum",
		     "No elements in k", 1);
  if (nk != pk.size())
    throw pofdExcept("powerSpectrum", "powerSpectrum",
		     "P(k) not same length as k", 2);
  if (nk < 3)
    throw pofdExcept("powerSpectrum", "powerSpectrum",
		     "Need at least 3 elements in k, P(k)", 3);

  logk = new double[nk];
  mink = k[0]; // We assume k is sorted
  maxk = k[nk-1];
  if (mink <= 0.0) throw pofdExcept("powerSpectrum", "powerSpectrum",
				    "Invalid (non-positive) k", 4);

  for (unsigned int i = 0; i < nk; ++i)
    logk[i] = log(k[i]);

  logpk = new double[nk];
  for (unsigned int i = 0; i < nk; ++i)
    if (k[i] <= 0.0)
      throw pofdExcept("powerSpectrum", "powerSpectrum",
		       "Invalid (non-positive) P(k)", 5);
  for (unsigned int i = 0; i < nk; ++i)
    logpk[i] = log(pk[i]);

  acc = gsl_interp_accel_alloc();
  spline_logpk = gsl_spline_alloc(gsl_interp_cspline,
				  static_cast<size_t>(nk));
  gsl_spline_init(spline_logpk, logk, logpk, static_cast<size_t>(nk));
}



/*!
  \param[in] k k vector in inverse arcmin.  Assumed sorted
  \param[in] pk P(k)
*/
void powerSpectrum::init(const std::vector<double>& k,
			 const std::vector<double>& pk) {
  nk = k.size();
  if (nk == 0)
    throw pofdExcept("powerSpectrum", "powerSpectrum",
		     "No elements in k", 1);
  if (nk != pk.size())
    throw pofdExcept("powerSpectrum", "powerSpectrum",
		     "P(k) not same length as k", 2);
  if (nk < 3)
    throw pofdExcept("powerSpectrum", "powerSpectrum",
		     "Need at least 3 elements in k, P(k)", 3);

  logk = new double[nk];
  mink = k[0]; // We assume k is sorted
  maxk = k[nk-1];
  if (mink <= 0.0) throw pofdExcept("powerSpectrum", "powerSpectrum",
				    "Invalid (non-positive) k", 4);

  for (unsigned int i = 0; i < nk; ++i)
    logk[i] = log(k[i]);

  logpk = new double[nk];
  for (unsigned int i = 0; i < nk; ++i)
    if (k[i] <= 0.0)
      throw pofdExcept("powerSpectrum", "powerSpectrum",
		       "Invalid (non-positive) P(k)", 5);
  for (unsigned int i = 0; i < nk; ++i)
    logpk[i] = log(pk[i]);

  acc = gsl_interp_accel_alloc();
  spline_logpk = gsl_spline_alloc(gsl_interp_cspline,
				  static_cast<size_t>(nk));
  gsl_spline_init(spline_logpk, logk, logpk, static_cast<size_t>(nk));
}


double positionGenerator::getPk(double k) const {
  // Taking advantage of the fact that mink must be positive
  if (k < mink || k > maxk) return 0;
  return exp(gsl_spline_eval(spline_logpk, log(k), acc));
}

double positionGenerator::getLogPk(double k) const {
  if (k < mink || k > maxk) std::numeric_limits<double>::quiet_NaN();
  return gsl_spline_eval(spline_logpk, log(k), acc);
}

/////////////////// positionGeneratorClustered ///////////////////
positionGeneratorClustered(unsigned int NX, unsigned int NY, double pixsize,
			   const std::string& powerfile) {

  if (NX == 0)
    throw pofdExcept("positionGeneratorClustered", 
		     "positionGeneratorClustered", 
		     "invalid (0) x dimension", 1);
  if (NY == 0)
    throw pofdExcept("positionGeneratorClustered", 
		     "positionGeneratorClustered", 
		     "invalid (0) y dimension", 1);

  nx = NX;
  ny = NY;
  probarr = NULL;
  probarr_trans = NULL;

  // Now, generate k
  powerSpectum pspec(powerfile);
  
  unsigned int nkx = nx / 2 + 1;
  double kxfac = tail_invk / nx;
  unsigned int nky = ny / 2 + 1;
  double kyfac = tail_invk / ny;
  double kxfacsq = kxfac * kxfac;
  double kyfacsq = kyfac * kyfac;

  // fill in k; k first stores the k values, then
  // stores the p(k) values
  unsigned int xtmp, ytmp;
  double kxval, kyval;
  k = new double[nx * ny];
  for (unsigned int i = 0; i < nkx; ++i) {
    kxval = kxfacsq * static_cast<double>(i * i);
    for (unsigned int j = 0; j < nky; ++j) {
      kyval = kyfacsq * static_cast<double>(j * j);
      k[ny * i + j] = sqrt(kxval + kyval);
    }
    for (unsigned int j = nky; j < ny; ++j) {
      ytmp = ny - j;
      kyval = kyfacsq * static_cast<double>(ytmp * ytmp);
      k[ny * i + j] = sqrt(kxval + kyval);
    }
  }
  for (unsigned int i = nkx; i < nx; ++i) {
    xtmp = nx - i;
    kxval = kxfacsq * static_cast<double>(xtmp * xtmp);
    for (unsigned int j = 0; j < nky; ++j) {
      kyval = kyfacsq * static_cast<double>(j * j);
      k[ny * i + j] = sqrt(kxval + kyval);
    }
    for (unsigned int j = nky; j < ny; ++j) {
      ytmp = ny - j;
      kyval = kyfacsq * static_cast<double>(ytmp * ytmp);
      k[ny * i + j] = sqrt(kxval + kyval);
    }
  }

  // Apply the power spectrum
  for (unsigned int i = 0; i < nx * ny; ++j)
    k[i] = powspec.getPk(k[i]);
  k[0] = 0.0; //Mean 0
}

positionGeneratorClustered::~positionGeneratorClustered() {
  delete[] k;
  if (probarr != NULL) fftw_free(probarr);
  if (probarr_trans != NULL) fftw_free(probarr_trans);
  if (plan != NULL) fftw_destroy_plan(plan); 
  if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
}

/*!
  \param[in] rangen Random number generator
 */
void positionGeneratorClustered::generate(ran& rangen) {

  if (probarr == NULL) {
    probarr = (*double) fftw_malloc(sizeof(double) * nx * ny);
    unsigned int nyhalf = ny / 2 + 1;
    probarr_trans = (fftw_complex*) 
      fftw_malloc(sizeof(fftw_complex) * nx * nyhalf)
    unsigned intx = static_cast<int>(nx);
    unsigned inty = static_cast<int>(ny);
    plan = fftw_plan_dft_r2c_2d(intx, inty, probarr, probarr_trans,
				FFTW_ESTIMATE);
    if (plan == NULL) {
      std::stringstream str;
      str << "Plan creation failed for forward transform of size: " << 
	  intx << " by " << inty;
      throw pofdExcept("positionGeneratorClustered", "generate",
		       str.str(), 1);
    }
    plan_inv = fftw_plan_dft_c2r_2d(intx, inty, probarr_trans, probarr,
				    FFTW_ESTIMATE);
    if (plan_inv == NULL) {
      std::stringstream str;
      str << "Inverse plan creation failed for forward transform of size: " << 
	  intx << " by " << inty;
      throw pofdExcept("positionGeneratorClustered", "generate",
		       str.str(), 2);
    }
  }

  // Fill
  for (unsigned int i = 0; i < nx * ny; ++i)
    probarr[i] = rangen.doub();
  
  // Transform
  fftw_execute(plan);

  // Multiply by k
  double val;
  for (unsigned int i = 0; i < nx; ++i)
    for (unsigned int j = 0; j <= nyhalf; ++j) {
      val = k[i * ny + j];
      probarr_trans[i * nyhalf + j][0] *= val;
      probarr_trans[i * nyhalf + j][1] *= val;
    }

  // Transform back
  fftw_execute(plan_inv);

  // Renormalize
  double minval, maxval, currval;
  minval = maxval = probarr[0];
  for (unsigned int i = 0; i < nx * ny; ++i) {
    currval = probarr[i * ny + j];
    if (currval > maxval) maxval = currval;
    if (currval < minval) minval = currval;
  }
  
  double rng = maxval - minval;
  if (rng == 0)
    for (unsigned int i = 0; i < nx * ny; ++i)
      probarr[i] -= minval - 1.0; //Make it all one
  else {
    double irng = 1.0 / rng;
    for (unsigned int i = 0; i < nx * ny; ++i)
      probarr[i] = (probarr[i] - minval) * irng; 
  }
    
}

std::pair<unsigned int> 
positinGeneratorClustered::getPosition(ran& rangen) const {
  unsigned int idx1, idx2;
  double zval;

  do {
    idx1 = static_cast<unsigned int>(rangen.doub() * nx);
    idx2 = static_cast<unsigned int>(rangen.doub() * ny);
    z = rangen.doub();
  } while (z > probarr[idx1 * ny + idx2]);

  return std::make_pair(idx1, idx2);
}
 
