//beam.cc

#include<sstream>

#include "../include/beam.h"
#include "../include/global_settings.h"
#include "../include/pofdExcept.h"

/*!
  \param[in] FWHM     FWHM of beam, in arcsec
*/
beam::beam( double FWHM ) { setFWHM(FWHM); }

/*!
  \param[in] fwhm_    New FWHM of beam, in arcsec
*/
void beam::setFWHM(double fwhm_) {
  fwhm = fwhm_;
  rhosq = pofd_coverage::rhofac / (fwhm*fwhm);
}

double beam::getEffectiveArea() const {
  return pofd_coverage::pi/rhosq;
}

double beam::getEffectiveAreaSq() const {
  return 0.5*pofd_coverage::pi/rhosq;
}

/*!
  \params[in] n  Number of pixels along each dimension.  Should be odd
  \params[in] pixsize Size of pixels in arcsec
  \params[out] fac Returns beam factor.  Must be pre-allocated by
                   caller and be of length n

  The beam is center normalized.	     
 */
void beam::getBeamFac(unsigned int n, double pixsize,
		      double* const fac) const {

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("beam", "getBeamFac", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("beam", "getBeamFac", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("beam", "getBeamFac", errstr.str(), 3);
  }
  if (fac == NULL)
    throw pofdExcept("beam", "getBeamFac", "fac is not allocated", 4);

  if (n == 1) {
    fac[0] = 1.0;
    return;
  }

  double sig = fwhm * pofd_coverage::fwhmfac / pixsize; //Gaussian sigma in pix
  double expfac = - 0.5 / (sig * sig);
  int ni = static_cast<int>(n);
  int no2 = ni / 2;
  double dist;
  for (int i = 0; i < ni; ++i) {
    dist = fabs(i - no2);
    fac[i] = exp(expfac * dist * dist);
  }
}


/*!
  \params[in] n  Number of pixels along each dimension.  Should be odd
  \params[in] pixsize Size of pixels in arcsec
  \params[out] bm Returns beam in row-major order.  Must be pre-allocated by
                   caller and be of length n * n

  The beam is center normalized.	     
 */
void beam::getBeam(unsigned int n, double pixsize,
		   double* const bm) const {

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("beam", "getBeam", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("beam", "getBeam", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("beam", "getBeam", errstr.str(), 3);
  }
  if (bm == NULL)
    throw pofdExcept("beam", "getBeam", "bm is not allocated", 4);

  if (n == 1) {
    bm[0] = 1.0;
    return;
  }

  //Make factorized beam, then multiply
  double* fac = new double[n];
  double sig = fwhm * pofd_coverage::fwhmfac / pixsize; //Gaussian sigma in pix
  double expfac = - 0.5 / (sig * sig);
  int ni = static_cast<int>(n);
  int no2 = ni / 2;
  double dist;
  for (int i = 0; i < ni; ++i) {
    dist = fabs(i - no2);
    fac[i] = exp(expfac * dist * dist);
  }

  double val;
  double* rowptr;
  for (unsigned int i = 0; i < n; ++i) {
    val = fac[i];
    rowptr = bm + i * n;
    for (unsigned int j = 0; j < n; ++j)
      rowptr[j] = val * fac[j];
  }

  delete[] fac;
}

/*!
  \params[in] n  Number of pixels along each dimension.  Should be odd
  \params[in] pixsize Size of pixels in arcsec
  \params[in] nbins Number of bins of output
  \param[out] nnonzero The number of elements of wt/bm that are filled
  \params[out] wt Returns histogrammed beam weights or order nbins
  \params[out] bm Returns histogrammed beam values.  
                   Must be pre-allocated by caller and be of length nbins
  \params[in] inverse If set, return one over the beam (for non-zero entries)

  The beam is center normalized, and binned in log sized bins.
  The first nnonzero entries contain the actual data.
 */
void beam::getBeamHist(unsigned int n, double pixsize,
		       unsigned int nbins,
		       unsigned int& nnonzero,
		       unsigned int* const wt,
		       double* const bm,
		       bool inverse) const {

  const double minval = 1e-5; //Always ignore beam values below this

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("beam", "getBeamHist", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("beam", "getBeamHist", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("beam", "getBeamHist", errstr.str(), 3);
  }
  if (nbins == 0)
    throw pofdExcept("beam", "getBeamHist", "nbins must be positive", 4);
  if (wt == NULL)
    throw pofdExcept("beam", "getBeamHist", "wt is not allocated", 5);
  if (bm == NULL)
    throw pofdExcept("beam", "getBeamHist", "bm is not allocated", 6);

  if (n == 1) {
    wt[0] = 1;
    bm[0] = 1.0;
    nnonzero = 1;
    return;
  }

  //Make factorized beam
  double* fac = new double[n];
  double sig = fwhm * pofd_coverage::fwhmfac / pixsize; //Gaussian sigma in pix
  double expfac = - 0.5 / (sig * sig);
  int ni = static_cast<int>(n);
  int no2 = ni / 2;
  double dist;
  for (int i = 0; i < ni; ++i) {
    dist = fabs(i - no2);
    fac[i] = exp(expfac * dist * dist);
  }

  double cminval = (minval > fac[0]*fac[0]) ? minval : 0.999 * fac[0]*fac[0];
  double minbinval = log2(cminval);
  double maxbinval = log2(1.001); //Assuming beam peaks at one.
  double histstep = (maxbinval - minbinval) / static_cast<double>(nbins);

  unsigned int* initwt = new unsigned int[nbins];
  double* meanbinval = new double[nbins];
  for (unsigned int i = 0; i < nbins; ++i) meanbinval[i] = 0.0;
  for (unsigned int i = 0; i < nbins; ++i) initwt[i] = 0;

  unsigned int idx;
  double fval, val;
  for (unsigned int i = 0; i < n; ++i) {
    fval = fac[i];
    for (unsigned int j = 0; j < n; ++j) {
      val = fval * fac[j];
      if (val < cminval) continue;  //Ignore
      idx = static_cast<unsigned int>((log2(val) - minbinval) / histstep);
      meanbinval[idx] += val;
      initwt[idx] += 1;
    }
  }
  delete[] fac;

  //Now cut down to those elements that are actually filled and
  // copy accross.  We use the mean value of the entries in the
  // bin as that bins coordinate.
  nnonzero = 0;
  for (unsigned int i = 0; i < nbins; ++i) wt[i] = 0;
  for (unsigned int i = 0; i < nbins; ++i) bm[i] = 0.0;
  for (unsigned int i = 0; i < nbins; ++i) {
    if (initwt[i] != 0) {
      wt[nnonzero] = initwt[i];
      bm[nnonzero] = meanbinval[i] / static_cast<double>(initwt[i]);
      ++nnonzero;
    }
  }
  delete[] initwt;
  delete[] meanbinval;
  if (nnonzero == 0)
    throw pofdExcept("beam", "getBeamHist", "No binned elements", 7);  

  if (inverse)
    for (unsigned int i = 0; i < nnonzero; ++i)
      bm[i] = 1.0 / bm[i];

}
