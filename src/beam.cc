//beam.cc

#include<sstream>
#include<limits>

#include<fitsio.h>

#include "../include/beam.h"
#include "../include/global_settings.h"
#include "../include/pofdExcept.h"

/*!
  \param[in] FWHM     FWHM of beam, in arcsec
*/
beam::beam(double FWHM) { setFWHM(FWHM); }

/*!
  \param[in] fwhm_    New FWHM of beam, in arcsec
*/
void beam::setFWHM(double fwhm_) {
  fwhm = fwhm_;
  rhosq = pofd_coverage::rhofac / (fwhm*fwhm);
}

double beam::getEffectiveArea() const {
  return pofd_coverage::pi / rhosq;
}

double beam::getEffectiveAreaSq() const {
  return 0.5 * pofd_coverage::pi / rhosq;
}

/*!
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[out] fac Returns beam factor.  Must be pre-allocated by
                   caller and be of length n

  The beam is center normalized.  Filtering is not supported.	     
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
    throw pofdExcept("beam", "getBeamac", errstr.str(), 2);
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
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[out] bm Returns beam in row-major order.  Must be pre-allocated by
                   caller and be of length n * n

  The beam is center normalized.  Filtering not supported	     
 */
void beam::getRawBeam(unsigned int n, double pixsize,
		      double* const bm) const {

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 3);
  }
  if (bm == NULL)
    throw pofdExcept("beam", "getRawBeam", "bm is not allocated", 4);

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
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[in] oversamp Oversampling.  Must be odd
  \param[out] bm Returns beam in row-major order.  Must be pre-allocated by
                   caller and be of length n * n

  The beam is center normalized.  Note that the returned beam is
  always at the nominal sampling.  It's just that if oversampling is
  used this is generated by making a more finely sampled one and then
  repixelating.

  Filtering is not supported
 */
void beam::getRawBeam(unsigned int n, double pixsize, unsigned int oversamp,
		      double* const bm) const {

  //Quick return
  if (oversamp == 1) {
    getRawBeam(n, pixsize, bm);
    return;
  }

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 3);
  }
  if (oversamp == 0) {
    std::stringstream errstr;
    errstr << "oversampling (" << oversamp << ") should be positive";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 4);
  }
  if (oversamp % 2 == 0) {
    std::stringstream errstr;
    errstr << "oversampling (" << oversamp << ") should be odd";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 5);
  }
  if (bm == NULL)
    throw pofdExcept("beam", "getRawBeam", "bm is not allocated", 6);

  //Make factorized beam, then multiply and sum
  unsigned int ngen = oversamp * n;
  double pixgen = pixsize / static_cast<double>(oversamp);
  double* fac = new double[ngen];
  double sig = fwhm * pofd_coverage::fwhmfac / pixgen; //Gaussian sigma in pix
  double expfac = - 0.5 / (sig * sig);
  int ni = static_cast<int>(ngen);
  int no2 = ni / 2;
  double dist;
  double normfac = 1.0 / static_cast<double>(oversamp);
  for (int i = 0; i < ni; ++i) {
    dist = fabs(i - no2);
    fac[i] = normfac * exp(expfac * dist * dist);
  }

  // Zero out
  for (unsigned int i = 0; i < n * n; ++i)
    bm[i] = 0.0;

  // Form the sum
  double val;
  double* rowptr;
  unsigned int minip, maxip, minjp, maxjp;
  for (unsigned int i = 0; i < n; ++i) {
    rowptr = bm + i * n;
    minip = i * oversamp;
    maxip = minip + oversamp;
    for (unsigned int i_p = minip; i_p < maxip; ++i_p) {
      val = fac[i_p];
      for (unsigned int j = 0; j < n; ++j) {
	minjp = j * oversamp;
	maxjp = minjp + oversamp;
	for (unsigned int j_p = minjp; j_p < maxjp; ++j_p)
	  rowptr[j] += val * fac[j_p];
      }
    }
  }

  delete[] fac;
}

/*!
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[in] bmpos Beam. Must be pre-allocated by
                   caller and be of length n * n
  \param[in] filter Hi-pass filter to apply.  If null, don't apply filter

  The beam is center normalized.  Note that all spatial information is
  lost by splitting into neg/pos parts
 */
void beam::getBeam(unsigned int n, double pixsize, double* const bm, 
		   hipassFilter* const filter) const {

  // pre-filtered beam
  getRawBeam(n, pixsize, bm);

  // Apply filtering
  if (filter != NULL)
    filter->filter(pixsize, n, n, bm);
}


/*!
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[in] oversamp Oversampling.  Must be odd
  \param[in] bm Beam. Must be pre-allocated by
                   caller and be of length n * n
  \param[in] filter Hi-pass filter to apply.  If null, don't apply filter

  The beam is center normalized.  
*/
void beam::getBeam(unsigned int n, double pixsize, unsigned int oversamp,
		   double* const bm, hipassFilter* const filter) const {

  // Pre-filtered beam
  getRawBeam(n, pixsize, oversamp, bm);

  // Apply filtering
  if (filter != NULL)
    filter->filter(pixsize, n, n, bm);
}


/*!
  \param[in] outputfile Name of FITS file to write to
  \param[in] pixsize Pixel size, in arcseconds
  \param[in] nfwhm Number of FWHM out to go.
  \param[in] oversamp Oversampling to use.  Must be odd.  
  \param[in] filt High-pass filter to apply.  If NULL, no filtering is done.
  \param[in] inverse Compute inverse beam rather than beam.
*/
void beam::writeToFits(const std::string& outputfile, double pixsize, 
		       double nfwhm, unsigned int oversamp,
		       hipassFilter* const filt, bool inverse) const {
  if (nfwhm <= 0.0)
    throw pofdExcept("beam", "writeToFits", "Invalid (non-positive) nfwhm", 1);
  if (pixsize <= 0.0)
    throw pofdExcept("beam", "writeToFits", 
		     "Invalid (non-positive) pixel size", 3);
  if (oversamp % 2 == 0)
    throw pofdExcept("beam", "writeToFits", 
		     "Invalid (non-odd) oversampling", 4);

  // Figure out how many pixels we are going to use
  unsigned int npix = static_cast<unsigned int>(nfwhm * fwhm / pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;
  
  // Get
  double *bmtmp = (double*) fftw_malloc(sizeof(double) * npix * npix);
  getBeam(npix, pixsize, oversamp, bmtmp, filt); // Also filters

  // Write
  int status = 0;
  fitsfile *fp;
  
  fits_create_file(&fp, outputfile.c_str(), &status);
  if (status) {
    fits_report_error(stderr, status);
    return;
  }
  long axissize[2];
  axissize[0] = static_cast<long>(npix);
  axissize[1] = static_cast<long>(npix);
  fits_create_img(fp, DOUBLE_IMG, 2, axissize, &status);
  
  // Header
  double dtmp;
  fits_write_key(fp, TLOGICAL, const_cast<char*>("FWHM"), &inverse,
		 const_cast<char*>("Inverse beam?"), &status);
  dtmp = fwhm;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM"), &dtmp,
		 const_cast<char*>("FWHM [arcsec]"), &status);
  dtmp = pixsize;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("PIXSIZE"), &dtmp,
		 const_cast<char*>("Pixel scale [arcsec]"), &status);
  dtmp = nfwhm;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("NFWHM"), &dtmp,
		 const_cast<char*>("Number of FWHM out"), &status);
  if (filt != NULL) {
    dtmp = filt->getFiltScale();
    if (dtmp > 0)
      fits_write_key(fp, TDOUBLE, const_cast<char*>("FILTSCL"), &dtmp,
		     const_cast<char*>("Filtering scale [arcsec]"), &status);
  }
  if (oversamp > 1) {
    unsigned int utmp;
    utmp = oversamp;
    fits_write_key(fp, TUINT, const_cast<char*>("OVERSAMP"), &utmp,
		   const_cast<char*>("Oversampling"), &status);
  }
  fits_write_key(fp, TSTRING, const_cast<char*>("VERSION"),
		 const_cast<char*>(pofd_coverage::version), 
		 const_cast<char*>("pofd_coverage version"),
		 &status);
  
  fits_write_history(fp,const_cast<char*>("Beam from pofd_coverage"),
		     &status);
  fits_write_date(fp, &status);
  
  // Data.  Must transpose.  Well -- the beam is symmetric,
  // so actually we don't.  But to support possible future changes
  // in filtering, treat it as if it might not be symmetric
  double *tmpdata = new double[npix];
  long fpixel[2] = { 1, 1 };
  for ( unsigned int j = 0; j < npix; ++j ) {
    if (inverse)
      for (unsigned int i = 0; i < npix; ++i) tmpdata[i] = 
						1.0 / bmtmp[i * npix + j];
    else
      for (unsigned int i = 0; i < npix; ++i) tmpdata[i] = bmtmp[i * npix + j];
    fpixel[1] = static_cast<long>(j + 1);
    fits_write_pix(fp, TDOUBLE, fpixel, npix, tmpdata, &status);
  }
  delete[] tmpdata;
  
  //Close up and go home
  fftw_free(bmtmp);
  fits_close_file(fp, &status);
  if (status) {
    fits_report_error(stderr, status);
    return;
  }
}

////////////////////////////////////////////////

beamHist::beamHist(unsigned int NBINS) : has_data(false), inverse(false), 
					 nbins(0), fwhm(0.0), nfwhm(3.5),
					 pixsize(0.0), eff_area(0.0), 
					 oversamp(1), filtscale(0.0),
					 n_pos(0), n_neg(0) {
  if (NBINS == 0)
    throw pofdExcept("beamHist", "beamHist", "Invalid (non-positive) NBINS", 1);
  nbins = NBINS;
  wt_pos = new unsigned int[nbins];
  bm_pos = new double[nbins];
  wt_neg = new unsigned int[nbins];
  bm_neg = new double[nbins];
}

beamHist::~beamHist() {
  delete[] wt_pos;
  delete[] bm_pos;
  delete[] wt_neg;
  delete[] bm_neg;
}

/*!
  \returns Minimum/Maximum values of positive beam.  If this is the 
  inverse beam, then the minimum/maximum of the inverse positive beam
  are returned.
*/
std::pair<double, double> beamHist::getMinMaxPos() const {
  if (!has_data) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());
  if (n_pos == 0) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());

  double min, max, val;
  min = max = bm_pos[0];
  for (unsigned int i = 1; i < n_pos; ++i) {
    val = bm_pos[i];
    if (val > max) max = val; else if (val < min) min = val;
  }
  return std::make_pair(min, max);
}

/*!
  \returns Minimum/Maximum values of negative beam.  If this is the 
  inverse beam, then the minimum/maximum of the inverse negative beam
  are returned.
*/
std::pair<double, double> beamHist::getMinMaxNeg() const {
  if (!has_data) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());
  if (n_neg == 0) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());

  double min, max, val;
  min = max = bm_neg[0];
  for (unsigned int i = 1; i < n_neg; ++i) {
    val = bm_neg[i];
    if (val > max) max = val; else if (val < min) min = val;
  }
  return std::make_pair(min, max);
}

/*!
  \param[in] bm Beam we are getting the histogram for
  \param[in] num_fwhm Number of FWHM out we will go on beam
  \param[in] pixsz Pixel size (arcseconds)
  \param[in] filt High pass filter we will apply.  If NULL, don't filter.
  \param[in] inv Histogram the inverse beam
  \param[in] oversampling Oversampling of beam. Must be odd
*/
void beamHist::fill(const beam& bm, double num_fwhm, double pixsz,
		    hipassFilter* const filt, bool inv, 
		    unsigned int oversampling) {

  const double minval = 1e-5; //Always ignore beam values below this

  inverse = inv;
  pixsize = pixsz;
  fwhm = bm.getFWHM();
  oversamp = oversampling;
  nfwhm = num_fwhm;
  if (filt != NULL)
    filtscale = filt->getFiltScale();
  else
    filtscale = 0.0;

  // Get how many pixels we will go out
  unsigned int npix = static_cast<unsigned int>(nfwhm * fwhm / pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;

  // Get 2D beam
  // Temporary beam storage.  Must use fftw_malloc since we may filter
  double *bmtmp = (double*) fftw_malloc(sizeof(double) * npix * npix);
  bm.getBeam(npix, pixsize, oversamp, bmtmp, filt); // Also filters

  // Find effective area
  double val;
  val = bmtmp[0];
  eff_area = val * val;
  for (unsigned int i = 1; i < npix * npix; ++i) {
    val = fabs(bmtmp[i]);
    eff_area += val * val;
  }
  eff_area *= (pixsize * pixsize) / (3600.0 * 3600.0);
  
  // Histogram
  // Find minimum and maximum non-zero parts for neg/pos histograms
  bool has_pos = false;
  bool has_neg = false;
  double minbinval_pos, maxbinval_pos, minbinval_neg, maxbinval_neg;
  minbinval_pos = minbinval_neg = 1e100; // Will definitely never be this large
  maxbinval_pos = maxbinval_neg = -1.0;
  for (unsigned int i = 0; i < npix * npix; ++i) {
    val = bmtmp[i];
    // Ignore anything within [-minval, minval]
    if (val > minval) { // Positive part
      has_pos = true;
      if (val > maxbinval_pos) maxbinval_pos = val;
      if (val < minbinval_pos) minbinval_pos = val;
    } else if (val < -minval) { //Negative part
      val = fabs(val);
      has_neg = true;
      if (val > maxbinval_neg) maxbinval_neg = val;
      if (val < minbinval_neg) minbinval_neg = val;
    }
  }

  // Set bin size
  double histstep_pos, histstep_neg;
  const double log2outscale = 0.0014419741739063218; // log2 1.001
  if (has_pos) {
    minbinval_pos = log2(minbinval_pos) - log2outscale;
    maxbinval_pos = log2(maxbinval_pos) + log2outscale;
    histstep_pos = (maxbinval_pos - minbinval_pos) / static_cast<double>(nbins);
  } else histstep_pos = 0.0;
  if (has_neg) {
    minbinval_neg = log2(minbinval_neg) - log2outscale;
    maxbinval_neg = log2(maxbinval_neg) + log2outscale;
    histstep_neg = (maxbinval_neg - minbinval_neg) / static_cast<double>(nbins);
  } else histstep_neg = 0.0;

  // Actually histogram
  unsigned int idx, utmp;
  unsigned int *tmpwt;
  tmpwt = new unsigned int[nbins];
  double *tmphist;
  tmphist = new double[nbins];

  // Positive histogram
  n_pos = 0;
  for (unsigned int i = 0; i < nbins; ++i) wt_pos[i] = 0;
  for (unsigned int i = 0; i < nbins; ++i) bm_pos[i] = 0.0;
  double lval;
  if (has_pos) {
    for (unsigned int i = 0; i < nbins; ++i)
      tmpwt[i] = 0;
    for (unsigned int i = 0; i < nbins; ++i)
      tmphist[i] = 0.0;
    for (unsigned int i = 0; i < npix * npix; ++i) {
      val = bmtmp[i];
      lval = log2(val);
      if (lval >= minbinval_pos) {
	idx = static_cast<unsigned int>((lval - minbinval_pos) / 
					histstep_pos);
	tmphist[idx] += val;
	tmpwt[idx] += 1;
      }
    }
    for (unsigned int i = 0; i < nbins; ++i)
      if (tmpwt[i] > 0) ++n_pos;
    // Shuffle all the filled bins to the front of bm_pos
    idx = 0;
    for (unsigned int i = 0; i < nbins; ++i) {
      utmp = tmpwt[i];
      if (utmp > 0) {
	wt_pos[idx] = utmp;
	bm_pos[idx] = tmphist[i] / static_cast<double>(utmp);
	++idx;
      }
    }
    // Note we histogram in beam space, then invert.  For whatever reason,
    // this seems to work a lot better for P(D)
    if (inverse)
      for (unsigned int i = 0; i < n_pos; ++i)
	bm_pos[i] = 1.0 / bm_pos[i];
  }

  // Negative
  n_neg = 0;
  for (unsigned int i = 0; i < nbins; ++i) wt_neg[i] = 0;
  for (unsigned int i = 0; i < nbins; ++i) bm_neg[i] = 0.0;
  if (has_neg) {
    for (unsigned int i = 0; i < nbins; ++i)
      tmpwt[i] = 0;
    for (unsigned int i = 0; i < nbins; ++i)
      tmphist[i] = 0.0;
    for (unsigned int i = 0; i < npix * npix; ++i) {
      val = fabs(bmtmp[i]);
      lval = log2(val);
      if (lval >= minbinval_neg) {
	idx = static_cast<unsigned int>((lval - minbinval_neg) / 
					histstep_neg);
	tmphist[idx] += val;
	tmpwt[idx] += 1;
      }
    }
    for (unsigned int i = 0; i < nbins; ++i)
      if (tmpwt[i] > 0) ++n_neg;
    idx = 0;
    for (unsigned int i = 0; i < nbins; ++i) {
      utmp = tmpwt[i];
      if (utmp > 0) {
	wt_neg[idx] = utmp;
	bm_neg[idx] = tmphist[i] / static_cast<double>(utmp);
	++idx;
      }
    }
    if (inverse)
      for (unsigned int i = 0; i < n_neg; ++i)
	bm_neg[i] = 1.0 / bm_neg[i];
  }
  has_data = true;

  // Clean up
  fftw_free(bmtmp);
}

/*!
  \param[in] outputfile Name of FITS file to write to
*/
void beamHist::writeToFits(const std::string& outputfile) const {

  if (!has_data)
    throw pofdExcept("beamHist", "writeToFits", "Hist not filled", 1);

  // Write as binary table
  int status = 0;
  fitsfile *fp;
  
  fits_create_file(&fp, outputfile.c_str(), &status);
  if (status) {
    fits_report_error(stderr, status);
    return;
  }

  char* ttype[]= {"WEIGHT", "BEAM"};
  char* tform[] = {"1I", "1D"};
  
  int itmp;
  double dtmp;
  if (n_pos > 0) {
    fits_create_tbl(fp, BINARY_TBL, 0, 2, ttype, tform, NULL, "POSBEAM", 
		    &status );
    // Header
    itmp = inverse;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("INVERSE"), &itmp,
		   const_cast<char*>("Inverse beam?"), &status);
    dtmp = fwhm;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM"), &dtmp,
		   const_cast<char*>("FWHM [arcsec]"), &status);
    dtmp = pixsize;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("PIXSIZE"), &dtmp,
		   const_cast<char*>("Pixel scale [arcsec]"), &status);
    dtmp = nfwhm;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("NFWHM"), &dtmp,
		   const_cast<char*>("Number of FWHM out"), &status);
    if (filtscale > 0) {
      dtmp = filtscale;
      fits_write_key(fp, TDOUBLE, const_cast<char*>("FILTSCL"), &dtmp,
		     const_cast<char*>("Filtering scale [arcsec]"), &status);
    }
    if (oversamp > 1) {
      unsigned int utmp;
      utmp = oversamp;
      fits_write_key(fp, TUINT, const_cast<char*>("OVERSAMP"), &utmp,
		     const_cast<char*>("Oversampling"), &status);
    }
    fits_write_key(fp, TSTRING, const_cast<char*>("VERSION"),
		   const_cast<char*>(pofd_coverage::version), 
		   const_cast<char*>("pofd_coverage version"),
		   &status);
    fits_write_history(fp,const_cast<char*>("Beam from pofd_coverage"),
		       &status);
    fits_write_date(fp, &status);
  
    // Now do the data
    fits_insert_rows(fp, 0, n_pos, &status);
    for (unsigned int i = 0; i < n_pos; ++i) {
      itmp = static_cast<int>(wt_pos[i]);
      fits_write_col(fp, TINT, 1, i+1, 1, 1, &itmp, &status);
      dtmp= bm_pos[i];
      fits_write_col(fp, TDOUBLE, 2, i+1, 1, 1, &dtmp, &status);
    }

  }
  if (n_neg > 0) {
    fits_create_tbl(fp, BINARY_TBL, 0, 2, ttype, tform, NULL, "POSBEAM", 
		    &status );
    // Header
    itmp = inverse;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("INVERSE"), &itmp,
		   const_cast<char*>("Inverse beam?"), &status);
    dtmp = fwhm;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM"), &dtmp,
		   const_cast<char*>("FWHM [arcsec]"), &status);
    dtmp = pixsize;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("PIXSIZE"), &dtmp,
		   const_cast<char*>("Pixel scale [arcsec]"), &status);
    dtmp = nfwhm;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("NFWHM"), &dtmp,
		   const_cast<char*>("Number of FWHM out"), &status);
    if (filtscale > 0) {
      dtmp = filtscale;
      fits_write_key(fp, TDOUBLE, const_cast<char*>("FILTSCL"), &dtmp,
		     const_cast<char*>("Filtering scale [arcsec]"), &status);
    }
    if (oversamp > 1) {
      unsigned int utmp;
      utmp = oversamp;
      fits_write_key(fp, TUINT, const_cast<char*>("OVERSAMP"), &utmp,
		     const_cast<char*>("Oversampling"), &status);
    }
    fits_write_key(fp, TSTRING, const_cast<char*>("VERSION"),
		   const_cast<char*>(pofd_coverage::version), 
		   const_cast<char*>("pofd_coverage version"),
		   &status);
    fits_write_history(fp,const_cast<char*>("Beam from pofd_coverage"),
		       &status);
    fits_write_date(fp, &status);
  
    // Now do the data
    fits_insert_rows(fp, 0, n_neg, &status);
    for (unsigned int i = 0; i < n_neg; ++i) {
      itmp = static_cast<int>(wt_neg[i]);
      fits_write_col(fp, TINT, 1, i+1, 1, 1, &itmp, &status);
      dtmp= bm_neg[i];
      fits_write_col(fp, TDOUBLE, 2, i+1, 1, 1, &dtmp, &status);
    }
  }

  //Close up and go home
  fits_close_file(fp, &status);
  if (status) {
    fits_report_error(stderr, status);
    return;
  }
}

