//doublebeam.cc

#include<sstream>

#include "../include/doublebeam.h"
#include "../include/global_settings.h"
#include "../include/pofdExcept.h"

/*!
  \param[in] FWHM1     FWHM of beam, in arcsec, band 1
  \param[in] FWHM2     FWHM of beam, in arcsec, band 2
*/
doublebeam::doublebeam(double FWHM1, double FWHM2) { 
  setFWHM(FWHM1, FWHM2); 
}

/*!
  \param[in] FWHM1    New FWHM of doublebeam, in arcsec, band 1
  \param[in] FWHM2    New FWHM of doublebeam, in arcsec, band 2
*/
void doublebeam::setFWHM(double FWHM1, double FWHM2) {
  fwhm1 = FWHM1;
  rhosq1 = pofd_coverage::rhofac / (FWHM1 * FWHM1);
  fwhm2 = FWHM2;
  rhosq2 = pofd_coverage::rhofac / (FWHM2 * FWHM2);
}

double doublebeam::getEffectiveArea1() const {
  return pofd_coverage::pi/rhosq1;
}
double doublebeam::getEffectiveArea2() const {
  return pofd_coverage::pi/rhosq2;
}

double doublebeam::getEffectiveAreaSq1() const {
  return 0.5*pofd_coverage::pi/rhosq1;
}
double doublebeam::getEffectiveAreaSq2() const {
  return 0.5*pofd_coverage::pi/rhosq2;
}

/*!
  \params[in] sel Which band to use (1 or 2)
  \params[in] n  Number of pixels along each dimension.  Should be odd
  \params[in] pixsize Size of pixels in arcsec
  \params[out] fac Returns doublebeam factor.  Must be pre-allocated by
                   caller and be of length n

  The doublebeam is center normalized.	     
 */
void doublebeam::getBeamFac(unsigned int sel, unsigned int n, double pixsize,
			    double* const fac) const {

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("doublebeam", "getBeamFac", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("doublebeam", "getBeamFac", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("doublebeam", "getBeamFac", errstr.str(), 3);
  }
  if (fac == NULL)
    throw pofdExcept("doublebeam", "getBeamFac", "fac is not allocated", 4);

  if (n == 1) {
    fac[0] = 1.0;
    return;
  }

  double fwhm;
  if (sel == 1) 
    fwhm = fwhm1;
  else if (sel == 2)
    fwhm = fwhm2;
  else {
    std::stringstream errstr;
    errstr << "Invalid band selection " << sel << "; should be 1 or 2";
    throw pofdExcept("doublebeam", "getBeamFac", errstr.str(), 5);
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
  \params[out] fac Returns beam factor.  Must be pre-allocated by
                   caller and be of length n

  The beam is center normalized.	     
 */
void doublebeam::getBeamFac1(unsigned int n, double pixsize,
			     double* const fac) const {
  return getBeamFac(1, n, pixsize, fac);
}

/*!
  \params[in] n  Number of pixels along each dimension.  Should be odd
  \params[in] pixsize Size of pixels in arcsec
  \params[out] fac Returns beam factor.  Must be pre-allocated by
                   caller and be of length n

  The beam is center normalized.	     
 */
void doublebeam::getBeamFac2(unsigned int n, double pixsize,
			     double* const fac) const {
  return getBeamFac(2, n, pixsize, fac);
}

/*!
  \params[in] sel Which band to use (1 or 2)
  \params[in] n  Number of pixels along each dimension.  Should be odd
  \params[in] pixsize Size of pixels in arcsec
  \params[out] bm Returns beam in row-major order.  Must be pre-allocated by
                   caller and be of length n * n

  The beam is center normalized.	     
 */
void doublebeam::getBeam(unsigned int sel, unsigned int n, double pixsize,
			 double* const bm) const {

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("doublebeam", "getBeam", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("doublebeam", "getBeam", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("doublebeam", "getBeam", errstr.str(), 3);
  }
  if (bm == NULL)
    throw pofdExcept("doublebeam", "getBeam", "bm is not allocated", 4);

  if (n == 1) {
    bm[0] = 1.0;
    return;
  }

  double fwhm;
  if (sel == 1) 
    fwhm = fwhm1;
  else if (sel == 2)
    fwhm = fwhm2;
  else {
    std::stringstream errstr;
    errstr << "Invalid band selection " << sel << "; should be 1 or 2";
    throw pofdExcept("doublebeam", "getBeam", errstr.str(), 5);
  }

  //Make factorized doublebeam, then multiply
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
  \params[out] bm Returns beam in row-major order.  Must be pre-allocated by
                   caller and be of length n * n

  The beam is center normalized.	     
 */
void doublebeam::getBeam1(unsigned int n, double pixsize,
			  double* const bm) const {
  return getBeam(1, n, pixsize, bm);
}


/*!
  \params[in] n  Number of pixels along each dimension.  Should be odd
  \params[in] pixsize Size of pixels in arcsec
  \params[out] bm Returns beam in row-major order.  Must be pre-allocated by
                   caller and be of length n * n

  The beam is center normalized.	     
 */
void doublebeam::getBeam2(unsigned int n, double pixsize,
			  double* const bm) const {
  return getBeam(2, n, pixsize, bm);
}

/*!
  \params[in] n  Number of pixels along each dimension.  Should be odd
  \params[in] pixsize Size of pixels in arcsec
  \params[in] nbins Maximum number of bins of output along each dimension
  \param[out] nnonzero The number of elements of wt/bm1/bm2 that actually
     have content
  \params[out] wt Returns histogrammed beam weights, must be
     pre-allocated by caller to nbins*nbins
  \params[out] bm1 Returns histogrammed beam values, band 1
                   Must be pre-allocated by caller and be of length nbins^2
  \params[out] bm2 Returns histogrammed beam values, band 2
                   Must be pre-allocated by caller and be of length nbins^2
  \params[in] inverse If set, return one over the beam (for non-zero entries)

  The beams are center normalized, and binned in log sized bins.
  The first nnonzero entries contain the actual data.
  The joint histogram is returned because that is what is needed to
  compute P(D) variables like R.
 */
void doublebeam::getBeamHist(unsigned int n, double pixsize,
			     unsigned int nbins, unsigned int& nnonzero,
			     unsigned int* const wt, double* const bm1,
			     double* const bm2, bool inverse) const {

  //Always ignore any locations where -either- beam is below this
  const double minval = 1e-6; 

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("doublebeam", "getBeamHist", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("doublebeam", "getBeamHist", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("doublebeam", "getBeamHist", errstr.str(), 3);
  }

  if (nbins == 0)
    throw pofdExcept("doublebeam", "getBeamHist", "nbins must be positive", 4);
  if (wt == NULL)
    throw pofdExcept("doublebeam", "getBeamHist", "wt is not allocated", 5);
  if (bm1 == NULL)
    throw pofdExcept("doublebeam", "getBeamHist", "bm1 is not allocated", 6);
  if (bm2 == NULL)
    throw pofdExcept("doublebeam", "getBeamHist", "bm1 is not allocated", 7);

  if (n == 1) {
    wt[0] = 1;
    bm1[0] = 1.0;
    bm2[0] = 1.0;
    nnonzero = 1;
    return;
  }

  //Make factorized beam in each band.
  double* fac1 = new double[n];
  double sig1 = fwhm1 * pofd_coverage::fwhmfac / pixsize;//Gaussian sigma in pix
  double expfac = - 0.5 / (sig1 * sig1);
  int ni = static_cast<int>(n);
  int no2 = ni / 2;
  double dist;
  for (int i = 0; i < ni; ++i) {
    dist = fabs(i - no2);
    fac1[i] = exp(expfac * dist * dist);
  }
  double* fac2 = new double[n];
  double sig2 = fwhm2 * pofd_coverage::fwhmfac / pixsize;//Gaussian sigma in pix
  expfac = - 0.5 / (sig2 * sig2);
  for (int i = 0; i < ni; ++i) {
    dist = fabs(i - no2);
    fac2[i] = exp(expfac * dist * dist);
  }

  //Compute the range we will bin over in each dimension
  //First, compute the minimum beam value we will encounter in each
  // band. We can take advantage of the symmetry to compare the corner
  // value -- the lowest -- with our cutoff.  Also -- the beam is always
  // positive.  We must skip locations where either beam falls below
  // this value because we can't store the other value in our array
  double cminval1 = (minval > fac1[0]*fac1[0]) ? minval : 
    0.999 * fac1[0]*fac1[0];
  double cminval2 = (minval > fac2[0]*fac2[0]) ? minval : 
    0.999 * fac2[0]*fac2[0];
  double minbinval1 = log2(cminval1);
  double minbinval2 = log2(cminval2);
  double maxbinval = log2(1.001); //Each beam peaks at one
  double histstep1 = (maxbinval - minbinval1) / static_cast<double>(nbins);
  double histstep2 = (maxbinval - minbinval2) / static_cast<double>(nbins);

  unsigned int ntotbins = nbins * nbins;
  unsigned int* initwt = new unsigned int[ntotbins];
  double* meanbinval1 = new double[ntotbins];
  double* meanbinval2 = new double[ntotbins];
  for (unsigned int i = 0; i < ntotbins; ++i) initwt[i] = 0;
  for (unsigned int i = 0; i < ntotbins; ++i) meanbinval1[i] = 0.0;
  for (unsigned int i = 0; i < ntotbins; ++i) meanbinval2[i] = 0.0;

  unsigned int idx1, idx2, idx; //idx is the flattened index
  double fval1, fval2, val1, val2;
  for (unsigned int i = 0; i < n; ++i) {
    fval1 = fac1[i];
    fval2 = fac2[i];
    for (unsigned int j = 0; j < n; ++j) {
      val1 = fval1 * fac1[j];
      val2 = fval2 * fac2[j];
      if (val1 < cminval1 || val2 < cminval2) continue;  //Ignore
      idx1 = static_cast<unsigned int>((log2(val1) - minbinval1) / histstep1);
      idx2 = static_cast<unsigned int>((log2(val2) - minbinval2) / histstep2);
      idx = idx1*nbins + idx2;
      initwt[idx] += 1;
      meanbinval1[idx] += val1;
      meanbinval2[idx] += val2;
    }
  }
  delete[] fac1;
  delete[] fac2;

  //Now cut down to those elements that are actually filled and
  // copy accross.  We use the mean value of the entries in the
  // bin as that bins value
  for (unsigned int i = 0; i < ntotbins; ++i) wt[i] = 0;
  for (unsigned int i = 0; i < ntotbins; ++i) bm1[i] = 0.0;
  for (unsigned int i = 0; i < ntotbins; ++i) bm2[i] = 0.0;
  nnonzero = 0;
  double d_iwt;
  for (unsigned int i = 0; i < ntotbins; ++i) {
    if (initwt[i] != 0) {
      wt[nnonzero] = initwt[i];
      d_iwt = 1.0 / static_cast<double>(initwt[i]);
      bm1[nnonzero] = meanbinval1[i] * d_iwt;
      bm2[nnonzero] = meanbinval2[i] * d_iwt;
      ++nnonzero;
    }
  }
  delete[] initwt;
  delete[] meanbinval1;
  delete[] meanbinval2;
  if (nnonzero == 0)
    throw pofdExcept("doublebeam", "getBeamHist", "No binned elements", 7);  

  if (inverse) {
    for (unsigned int i = 0; i < nnonzero; ++i) bm1[i] = 1.0 / bm1[i];
    for (unsigned int i = 0; i < nnonzero; ++i) bm2[i] = 1.0 / bm2[i];
  }

}
