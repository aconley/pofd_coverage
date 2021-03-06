//beam.cc

#include<sstream>
#include<limits>
#include<cstring>
#include<ctime>

#include<fitsio.h>

#include "../include/beam.h"
#include "../include/pofdExcept.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

/*!
  \param[in] FWHM     FWHM of beam, in arcsec
*/
beam::beam(double FWHM) noexcept : noise(0.0), rangen(nullptr) { setFWHM(FWHM); }

beam::~beam() { if (rangen != nullptr) delete rangen; }

/*!
  \param[in] fwhm_    New FWHM of beam, in arcsec
*/
void beam::setFWHM(double fwhm_) noexcept {
  fwhm = fwhm_;
  rhosq = pofd_coverage::rhofac / (fwhm * fwhm);
}

double beam::getEffectiveArea() const noexcept {
  return pofd_coverage::pi / rhosq;
}

double beam::getEffectiveAreaSq() const noexcept {
  return 0.5 * pofd_coverage::pi / rhosq;
}

/*!
  \param[in] n Number of values in bm
  \param[in] noiseval Gaussian noise sigma
  \param[inout] bm Data to modify, treated as 1D array of length n
*/
void beam::addNoise(unsigned int n, double noiseval, double* const bm) const {
  if (noiseval == 0) return; // Negative noise is okay
  if (rangen == nullptr) {
    unsigned long long int seed;
    seed = static_cast<unsigned long long int>(time(nullptr));
    seed += static_cast<unsigned long long int>(clock());
    rangen = new ran(seed);
  }
  for (unsigned int i = 0; i < n; ++i)
    bm[i] += noiseval * rangen->gauss();
}

/*!
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[out] fac Returns beam factor.  Must be pre-allocated by
                   caller and be of length n

  The beam is center normalized.  Filtering is not supported.
  No noise is included.
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
  if (fac == nullptr)
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
    dist = abs(i - no2);
    fac[i] = exp(expfac * dist * dist);
  }
}


/*!
  \param[in] n1 Number of pixels along dimension 1
  \param[in] n2 Number of pixels along dimension 2
  \param[in] pixsize Size of pixels in arcsec
  \param[out] bm Returns beam in row-major order.  Must be pre-allocated by
                   caller and be of length n * n

  The beam is center normalized.  Filtering not supported.
  The beam is centered as well as possible, but unless n1 and n2 are odd,
   it can't quite be centered.  No noise is included.
*/
void beam::getRawBeam(unsigned int n1, unsigned int n2, double pixsize,
                      double* const bm) const {

  //Input checks
  if (n1 == 0) {
    std::stringstream errstr;
    errstr << "n1 (" << n1 << ") should be positive";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 1);
  }
  if (n2 == 0) {
    std::stringstream errstr;
    errstr << "n2 (" << n1 << ") should be positive";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 3);
  }
  if (bm == nullptr)
    throw pofdExcept("beam", "getRawBeam", "bm is not allocated", 4);

  if ((n1 == 1) && (n2 == 1)) {
    bm[0] = 1.0;
    return;
  }

  double sig = fwhm * pofd_coverage::fwhmfac / pixsize; //Gaussian sigma in pix
  double expfac = - 0.5 / (sig * sig);

  // Include a special case for n1 and n2 the same, which can
  // be more efficient
  if (n1 == n2) {
    //Make factorized beam, then multiply
    double* fac = new double[n1];
    int ni = static_cast<int>(n1);
    int no2 = ni / 2;
    double dist;
    for (int i = 0; i < ni; ++i) {
      dist = abs(i - no2);
      fac[i] = exp(expfac * dist * dist);
    }

    double val;
    double* rowptr;
    for (unsigned int i = 0; i < n1; ++i) {
      val = fac[i];
      rowptr = bm + i * n2;
      for (unsigned int j = 0; j < n1; ++j)
        rowptr[j] = val * fac[j];
    }

    delete[] fac;
  } else {
    double* fac1 = new double[n1];
    int ni = static_cast<int>(n1);
    int no2 = ni / 2;
    double dist;
    for (int i = 0; i < ni; ++i) {
      dist = abs(i - no2);
      fac1[i] = exp(expfac * dist * dist);
    }
    double* fac2 = new double[n2];
    ni = static_cast<int>(n2);
    no2 = ni / 2;
    for (int i = 0; i < ni; ++i) {
      dist = abs(i - no2);
      fac2[i] = exp(expfac * dist * dist);
    }

    double val;
    double* rowptr;
    for (unsigned int i = 0; i < n1; ++i) {
      val = fac1[i];
      rowptr = bm + i * n2;
      for (unsigned int j = 0; j < n2; ++j)
        rowptr[j] = val * fac2[j];
    }

    delete[] fac1;
    delete[] fac2;
  }
}


/*!
  \param[in] n1  Number of pixels along first dimension
  \param[in] n2  Number of pixels along second dimension
  \param[in] pixsize Size of pixels in arcsec
  \param[in] oversamp Oversampling.  Must be odd
  \param[out] bm Returns beam in row-major order.  Must be pre-allocated by
                   caller and be of length n * n

  The beam is center normalized.  Note that the returned beam is
  always at the nominal sampling.  It's just that if oversampling is
  used this is generated by making a more finely sampled one and then
  repixelating.

  Filtering is not supported

  The beam is centered as well as possible, but unless n1 and n2 are odd,
   it can't quite be centered.  No noise is included.
*/
void beam::getRawBeam(unsigned int n1, unsigned int n2,
                      double pixsize, unsigned int oversamp,
                      double* const bm) const {

  if (oversamp == 1) {
    getRawBeam(n1, n2, pixsize, bm);
    return;
  }

  //Input checks
  if (n1 == 0) {
    std::stringstream errstr;
    errstr << "n1 (" << n1 << ") should be positive";
    throw pofdExcept("beam", "getRawBeam", errstr.str(), 1);
  }
  if (n2 == 0) {
    std::stringstream errstr;
    errstr << "n2 (" << n2 << ") should be positive";
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
  if (bm == nullptr)
    throw pofdExcept("beam", "getRawBeam", "bm is not allocated", 6);

  double pixgen = pixsize / static_cast<double>(oversamp);
  double sig = fwhm * pofd_coverage::fwhmfac / pixgen; //Gaussian sigma in pix
  double expfac = - 0.5 / (sig * sig);

  // Zero out
  std::memset(bm, 0, n1 * n2 * sizeof(double));

  //Make factorized beam, then multiply and sum
  // Allow special case of n1 == n2, which can be more efficient
  double *fac1, *fac2;
  fac1 = fac2 = nullptr;

  unsigned int ngen = oversamp * n1;

  fac1 = new double[ngen];
  int ni = static_cast<int>(ngen);
  int no2 = ni / 2;
  double dist;
  double normfac = 1.0 / static_cast<double>(oversamp);
  for (int i = 0; i < ni; ++i) {
    dist = abs(i - no2);
    fac1[i] = normfac * exp(expfac * dist * dist);
  }

  if (n1 != n2) {
    // Need a different fac2
    ngen = oversamp * n2;
    ni = static_cast<int>(ngen);
    no2 = ni / 2;
    for (int i = 0; i < ni; ++i) {
      dist = abs(i - no2);
      fac2[i] = normfac * exp(expfac * dist * dist);
    }
  } else fac2 = fac1; //Re-use -- but be careful to clean up correctly!

  // Form the sum
  double val;
  double* rowptr;
  unsigned int minip, maxip, minjp, maxjp;
  for (unsigned int i = 0; i < n1; ++i) {
    rowptr = bm + i * n2;
    minip = i * oversamp;
    maxip = minip + oversamp;
    for (unsigned int i_p = minip; i_p < maxip; ++i_p) {
      val = fac1[i_p];
      for (unsigned int j = 0; j < n2; ++j) {
        minjp = j * oversamp;
        maxjp = minjp + oversamp;
        for (unsigned int j_p = minjp; j_p < maxjp; ++j_p)
          rowptr[j] += val * fac2[j_p];
      }
    }
  }

  delete[] fac1; // Always allocated
  if (n1 != n2) delete[] fac2; // Only allocated if this is true
}

/*!
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[in] bm Beam. Must be pre-allocated by caller and be of length n * n
  \param[in] filter Hi-pass/matched filter to apply.  If null, don't apply
                     filtering

  This is only properly centered if n is odd.  Noise will be added if set.
*/
void beam::getBeam(unsigned int n, double pixsize, double* const bm,
                   const fourierFilter* const filter) const {

  // pre-filtered beam
  getRawBeam(n, n, pixsize, bm);

  // Add noise if needed
  if (noise > 0) addNoise(n * n, noise, bm);

  // Apply filtering to beam
  if (filter != nullptr)
    filter->filter(n, n, pixsize, bm);
}

/*!
  \param[in] n1 Number of pixels along dimension 1
  \param[in] n2 Number of pixels along dimension 2
  \param[in] pixsize Size of pixels in arcsec
  \param[in] bm Beam. Must be pre-allocated by caller and be of length n * n
  \param[in] filter Hi-pass/matched filter to apply.  If null, don't apply
                     filtering

  This is only properly centered if n1 and n2 are odd.  Noise will be
  added if set.
*/
void beam::getBeam(unsigned int n1, unsigned int n2,
                   double pixsize, double* const bm,
                   const fourierFilter* const filter) const {

  // pre-filtered beam
  getRawBeam(n1, n2, pixsize, bm);

  // Add noise if needed
  if (noise > 0) addNoise(n1 * n2, noise, bm);

  // Apply filtering to beam
  if (filter != nullptr)
    filter->filter(n1, n2, pixsize, bm);
}


/*!
  \param[in] n  Number of pixels along each dimension.
  \param[in] pixsize Size of pixels in arcsec
  \param[in] oversamp Oversampling.  Must be odd
  \param[in] bm Beam. Must be pre-allocated by
                   caller and be of length n * n
  \param[in] filter Hi-pass/Matched filter to apply.  If null, don't
                apply filter

  This is only properly centered if n is odd.
*/
void beam::getBeam(unsigned int n, double pixsize, unsigned int oversamp,
                   double* const bm, const fourierFilter* const filter) const {

  // Pre-filtered beam
  getRawBeam(n, n, pixsize, oversamp, bm);

  // Apply filtering
  if (filter != nullptr)
    filter->filter(n, n, pixsize, bm);
}

/*!
  \param[in] n1 Number of pixels along dimension 1
  \param[in] n2 Number of pixels along dimension 2
  \param[in] pixsize Size of pixels in arcsec
  \param[in] oversamp Oversampling.  Must be odd
  \param[in] bm Beam. Must be pre-allocated by
                   caller and be of length n * n
  \param[in] filter Hi-pass/Matched filter to apply.  If null, don't
                apply filter

  This is only properly centered if n is odd.
*/
void beam::getBeam(unsigned int n1, unsigned int n2, double pixsize,
                   unsigned int oversamp, double* const bm,
                   const fourierFilter* const filter) const {

  // Pre-filtered beam
  getRawBeam(n1, n2, pixsize, oversamp, bm);

  // Apply filtering
  if (filter != nullptr)
    filter->filter(n1, n2, pixsize, bm);
}



/*!
  \param[in] outputfile Name of FITS file to write to
  \param[in] pixsize Pixel size, in arcseconds
  \param[in] nfwhm Number of FWHM out to go.
  \param[in] oversamp Oversampling to use.  Must be odd.
  \param[in] filt High-pass filter to apply.  If nullptr, no filtering is done.
  \param[in] inverse Compute inverse beam rather than beam.
*/
void beam::writeToFits(const std::string& outputfile, double pixsize,
                       double nfwhm, unsigned int oversamp,
                       const fourierFilter* const filt, bool inverse) const {
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
  bool ubl;
  double dtmp;
  unsigned int utmp;
  fits_write_key(fp, TLOGICAL, const_cast<char*>("INVERSE"), &inverse,
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
  dtmp = noise;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("NOISE"), &dtmp,
                 const_cast<char*>("Fractional noise"), &status);

  if (filt != nullptr) {
    ubl = true;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("FILT"), &ubl,
                   const_cast<char*>("Filtered?"), &status);
    ubl = filt->isHipass();
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIPASS"), &ubl,
                   const_cast<char*>("Hipass filtered?"), &status);
    if (ubl) {
      dtmp = filt->getFiltScale();
      if (dtmp > 0)
        fits_write_key(fp, TDOUBLE, const_cast<char*>("FILTSCL"), &dtmp,
                       const_cast<char*>("Hipass filtering scale [arcsec]"),
                       &status);
      dtmp = filt->getQFactor();
      if (dtmp > 0)
        fits_write_key(fp, TDOUBLE, const_cast<char*>("QFACTOR"), &dtmp,
                       const_cast<char*>("Hipass filtering apodization"),
                       &status);
    }
    ubl = filt->isMatched();
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MATCHED"), &ubl,
                   const_cast<char*>("Match filtered?"), &status);
    if (ubl) {
      dtmp = filt->getSigInst();
      if (dtmp > 0)
        fits_write_key(fp, TDOUBLE, const_cast<char*>("MATSIGI"), &dtmp,
                       const_cast<char*>("Matched filtering sig_i [Jy]"),
                       &status);
      dtmp = filt->getSigConf();
      if (dtmp > 0)
        fits_write_key(fp, TDOUBLE, const_cast<char*>("MATSIGC"), &dtmp,
                       const_cast<char*>("Matched filtering sig_c [Jy]"),
                       &status);
    }
  } else {
    ubl = false;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("FILT"), &ubl,
                   const_cast<char*>("Filtered?"), &status);
  }
  if (oversamp > 1) {
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
  for (unsigned int j = 0; j < npix; ++j) {
    if (inverse)
      for (unsigned int i = 0; i < npix; ++i)
        tmpdata[i] = 1.0 / bmtmp[i * npix + j];
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

/*!
  \param[in] NBINS Number of bins in histogram
*/
beamHist::beamHist(unsigned int NBINS):
  has_data(false), inverse(false), nbins(0), fwhm(0.0), nfwhm(4.0),
  nfwhmkeep(4.0), pixsize(0.0), eff_area(0.0), oversamp(1),
  n_pos(0), n_neg(0) {

  if (NBINS == 0)
    throw pofdExcept("beamHist", "beamHist", "Invalid (non-positive) NBINS", 1);
  nbins = NBINS;
  wt_pos = new unsigned int[nbins];
  bm_pos = new double[nbins];
  wt_neg = new unsigned int[nbins];
  bm_neg = new double[nbins];

  minmax_pos = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN());
  minmax_neg = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN());

  isHipass = false;
  filtscale = std::numeric_limits<double>::quiet_NaN();
  qfactor = std::numeric_limits<double>::quiet_NaN();
  isMatched = false;
  matched_fwhm = std::numeric_limits<double>::quiet_NaN();
  matched_sigi = std::numeric_limits<double>::quiet_NaN();
  matched_sigc = std::numeric_limits<double>::quiet_NaN();
}

beamHist::~beamHist() {
  delete[] wt_pos;
  delete[] bm_pos;
  delete[] wt_neg;
  delete[] bm_neg;
}

/*!
  \param[in] bm Beam we are getting the histogram for
  \param[in] num_fwhm Number of FWHM out we will go on beam
  \param[in] pixsz Pixel size (arcseconds)
  \param[in] inv Histogram the inverse beam
  \param[in] oversampling Oversampling of beam. Must be odd
  \param[in] filt Fourier space filter to apply
  \param[in] num_fwhm_keep How many FWHM to keep in the histogram
              after filtering.  If 0 (the default), keeps everything.
*/
void beamHist::fill(const beam& bm, double num_fwhm, double pixsz,
                    bool inv, unsigned int oversampling,
                    const fourierFilter* const filt,
                    double num_fwhm_keep) {

  // Get how many pixels we will go out
  unsigned int npix =
    static_cast<unsigned int>(num_fwhm * bm.getFWHM() / pixsz +
                              0.9999999999);
  npix = 2 * npix + 1;
  fill(bm, npix, npix, pixsz, inv, oversampling, filt, num_fwhm_keep);
}

/*!
  \param[in] bm Beam we are getting the histogram for
  \param[in] n1 Size to build beam from, dimension 1
  \param[in] n2 Size to build beam from, dimension 2
  \param[in] pixsz Pixel size (arcseconds)
  \param[in] inv Histogram the inverse beam
  \param[in] oversampling Oversampling of beam. Must be odd
  \param[in] filt Fourier space filter to apply
  \param[in] num_fwhm_keep How many FWHM to keep in the histogram
              after filtering.  If 0 (the default), keeps everything.
              Note that if this is applied, a symmetric part of
              the beam is kept.
*/
void beamHist::fill(const beam& bm, unsigned int n1, unsigned int n2,
                    double pixsz, bool inv, unsigned int oversampling,
                    const fourierFilter* const filt,
                    double num_fwhm_keep) {

  const double minval = 1e-5; //Always ignore beam values below this

  inverse = inv;
  pixsize = pixsz;
  fwhm = bm.getFWHM();
  oversamp = oversampling;
  unsigned int minn = std::min(n1, n2);
  nfwhm = static_cast<double>(minn-1) / (2.0 * fwhm);
  eff_area = bm.getEffectiveArea();
  if (num_fwhm_keep == 0) nfwhmkeep = nfwhm;
  else nfwhmkeep = std::min(nfwhm, num_fwhm_keep);

  // Get 2D beam
  // Temporary beam storage.  Must use fftw_malloc since we may filter
  double val;
  double *bmtmp = (double*) fftw_malloc(sizeof(double) * n1 * n2);
  // Get the beam, filtering as needed
  bm.getBeam(n1, n2, pixsize, oversamp, bmtmp, filt);

  // Store information about the filtering.  This is done purely
  //  in case we want to write the beam out.
  isHipass = false;
  filtscale = std::numeric_limits<double>::quiet_NaN();
  qfactor = std::numeric_limits<double>::quiet_NaN();
  isMatched = false;
  matched_fwhm = std::numeric_limits<double>::quiet_NaN();
  matched_sigi = std::numeric_limits<double>::quiet_NaN();
  matched_sigc = std::numeric_limits<double>::quiet_NaN();
  if (filt != nullptr) {
    if (filt->isHipass()) {
      isHipass = true;
      filtscale = filt->getFiltScale();
      qfactor = filt->getQFactor();
    }
    if (filt->isMatched()) {
      isMatched = true;
      matched_fwhm = filt->getFWHM();
      matched_sigi = filt->getSigInst();
      matched_sigc = filt->getSigConf();
    }
  }

  minmax_pos = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN());
  minmax_neg = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN());

  // Set up the part we will actually use (which may mean clipping)
  unsigned int minidx1, minidx2, maxidx1, maxidx2;
  if (nfwhmkeep < nfwhm) {
    // We want to set up logical indexing into the array to only keep
    //  the part we want.  Recall the beam is centered at n1/2, n2/2 (integer)
    unsigned int nclippix =
      static_cast<unsigned int>(nfwhmkeep * fwhm / pixsize + 0.9999999999);
    unsigned int nclippix_tot = 2 * nclippix + 1;
    // Recall that the beam is centered at n1/2, n2/2
    if (nclippix_tot < n1) {
      minidx1 = n1 / 2 - nclippix; // Range is [minidx1, maxidx1)
      maxidx1 = minidx1 + nclippix_tot;
      // assert(minidx1 < n1); // unsigned, will wrap on problem
      // assert(maxidx1 <= n1);
    } else {
      minidx1 = 0;
      maxidx1 = n1;
    }
    if (nclippix_tot < n2) {
      minidx2 = n2 / 2 - nclippix;
      maxidx2 = minidx2 + nclippix_tot;
      // assert(minidx2 < n2);
      // assert(maxidx2 <= n2);
    } else {
      minidx2 = 0;
      maxidx2 = n2;
    }
  } else {
    minidx1 = minidx2 = 0;
    maxidx1 = n1;
    maxidx2 = n2;
  }

  // Histogram
  // Find minimum and maximum non-zero parts for neg/pos histograms
  bool has_pos = false;
  bool has_neg = false;
  double minbinval_pos, maxbinval_pos, minbinval_neg, maxbinval_neg;
  minbinval_pos = minbinval_neg = 1e100; // Will definitely never be this large
  maxbinval_pos = maxbinval_neg = -1.0;
  double *rowptr;
  for (unsigned int i = minidx1; i < maxidx1; ++i) {
    rowptr = bmtmp + i * n2; // Logical size npix by npix
    for (unsigned int j = minidx2; j < maxidx2; ++j) {
      val = rowptr[j];
      // Ignore anything within [-minval, minval]
      if (val > minval) { // Positive part
        has_pos = true;
        if (val > maxbinval_pos) maxbinval_pos = val;
        if (val < minbinval_pos) minbinval_pos = val;
      } else if (val < (-minval)) { //Negative part
        val = fabs(val);
        has_neg = true;
        if (val > maxbinval_neg) maxbinval_neg = val;
        if (val < minbinval_neg) minbinval_neg = val;
      }
    }
  }

  if (!(has_pos || has_neg))
    throw pofdExcept("beamHist", "fill",
                     "No pos or negative components", 1);

  if (has_pos) minmax_pos = std::make_pair(minbinval_pos, maxbinval_pos);
  if (has_pos) minmax_neg = std::make_pair(minbinval_neg, maxbinval_neg);
  if (!(has_pos || has_neg))
    throw pofdExcept("beamHist", "fill",
                     "Found no positive or negative pix", 1);

  // Set bin size
  double histstep_pos, histstep_neg;
  const double log2outscale = 0.0014419741739063218; // log2(1.001)
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

  // Do positive beam
  n_pos = 0;
  std::memset(wt_pos, 0, nbins * sizeof(unsigned int));
  std::memset(bm_pos, 0, nbins * sizeof(double));
  double lval;
  if (has_pos) {
    std::memset(tmpwt, 0, nbins * sizeof(unsigned int));
    std::memset(tmphist, 0, nbins * sizeof(double));
    for (unsigned int i = minidx1; i < maxidx1; ++i) {
      rowptr = bmtmp + i * n2;
      for (unsigned int j = minidx2; j < maxidx2; ++j) {
        val = rowptr[j];
        if (val <= minval) continue;  //skip: too close to zero or negative
        lval = log2(val);
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
  std::memset(wt_neg, 0, nbins * sizeof(unsigned int));
  std::memset(bm_neg, 0, nbins * sizeof(double));
  if (has_neg) {
    std::memset(tmpwt, 0, nbins * sizeof(unsigned int));
    std::memset(tmphist, 0, nbins * sizeof(double));
    double testval = -minval;
    for (unsigned int i = minidx1; i < maxidx1; ++i) {
      rowptr = bmtmp + i * n2;
      for (unsigned int j = minidx2; j < maxidx2; ++j) {
        val = rowptr[j];
        if (val > testval) continue;  //skip; too close to 0 or positive
        val = fabs(val); // Work with abs value
        lval = log2(val);
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

  // Empty primary HDU
  long axissize[2];
  axissize[0] = axissize[1] = 0;
  fits_create_img(fp, FLOAT_IMG, 2, axissize, &status);
  if (status) {
    fits_report_error(stderr, status);
    return;
  }

  // Primary header
  unsigned int utmp;
  int itmp;
  double dtmp;
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
  dtmp = nfwhmkeep;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("NFWHMKP"), &dtmp,
                 const_cast<char*>("Number of FWHM kept"), &status);

  if (isHipass) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT"), &itmp,
                   const_cast<char*>("Is beam hipass filtered?"), &status);
    dtmp = filtscale;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTSCL"), &dtmp,
                   const_cast<char*>("Filtering scale [arcsec]"), &status);
    dtmp = qfactor;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTQ"), &dtmp,
                   const_cast<char*>("Filtering apodization"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT"), &itmp,
                   const_cast<char*>("Is beam hipass filtered?"), &status);
  }
  if (isMatched) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT"), &itmp,
                   const_cast<char*>("Is beam match filtered?"), &status);
    dtmp = matched_fwhm;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITFWHM"), &dtmp,
                   const_cast<char*>("Matched filtering FWHM [arcsec]"),
                   &status);
    dtmp = matched_sigi;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGI"), &dtmp,
                   const_cast<char*>("Matched filtering sigi"), &status);
    dtmp = matched_sigc;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGC"), &dtmp,
                   const_cast<char*>("Matched filtering sigc"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT"), &itmp,
                   const_cast<char*>("Is beam match filtered?"), &status);
  }

  if (oversamp > 1) {
    utmp = oversamp;
    fits_write_key(fp, TUINT, const_cast<char*>("OVERSAMP"), &utmp,
                   const_cast<char*>("Oversampling"), &status);
  }

  utmp = nbins;
  fits_write_key(fp, TUINT, const_cast<char*>("NBINS"), &utmp,
                 const_cast<char*>("Number of bins"), &status);

  fits_write_key(fp, TSTRING, const_cast<char*>("VERSION"),
                 const_cast<char*>(pofd_coverage::version),
                 const_cast<char*>("pofd_coverage version"),
                 &status);
  fits_write_history(fp, const_cast<char*>("Beam from pofd_coverage"),
                     &status);
  fits_write_date(fp, &status);


  char* ttype[]= {"WEIGHT", "BEAM"};
  char* tform[] = {"1I", "1D"};

  // Pos beam
  if (n_pos > 0) {
    fits_create_tbl(fp, BINARY_TBL, 0, 2, ttype, tform, nullptr, "POSBEAM",
                    &status);
    fits_insert_rows(fp, 0, n_pos, &status);
    for (unsigned int i = 0; i < n_pos; ++i) {
      itmp = static_cast<int>(wt_pos[i]);
      fits_write_col(fp, TINT, 1, i+1, 1, 1, &itmp, &status);
      dtmp= bm_pos[i];
      fits_write_col(fp, TDOUBLE, 2, i+1, 1, 1, &dtmp, &status);
    }

  }
  if (n_neg > 0) {
    fits_create_tbl(fp, BINARY_TBL, 0, 2, ttype, tform, nullptr, "NEGBEAM",
                    &status);
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

