//doublebeam.cc

#include<sstream>
#include<limits>
#include<fitsio.h>
#include<cstring>

#include "../include/doublebeam.h"
#include "../include/pofdExcept.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

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

dblpair doublebeam::getEffectiveArea() const {
  return std::make_pair(pofd_coverage::pi/rhosq1,
			pofd_coverage::pi/rhosq2);
}

dblpair doublebeam::getEffectiveAreaSq() const {
  return std::make_pair(0.5 * pofd_coverage::pi/rhosq1,
			0.5 * pofd_coverage::pi/rhosq2);
}

/*!
  \param[in] band Which band to use (1 or 2)
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[out] fac Returns doublebeam factor.  Must be pre-allocated by
                   caller and be of length n

  The doublebeam is center normalized.	     
*/
void doublebeam::getBeamFac(unsigned int band, unsigned int n, double pixsize,
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
  if (band == 1) 
    fwhm = fwhm1;
  else if (band == 2)
    fwhm = fwhm2;
  else {
    std::stringstream errstr;
    errstr << "Invalid band selection " << band << "; should be 1 or 2";
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
  \param[in] band Which band to use (1 or 2)
  \param[in] n Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[out] bm Returns beam in row-major order.  Must be pre-allocated by
                   caller and be of length n * n

  The beam is center normalized.  Filtering not supported	     
*/
void doublebeam::getRawBeam(unsigned int band, unsigned int n, double pixsize,
			    double* const bm) const {

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("doublebeam", "getRawBeam", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("doublebeam", "getRawBeam", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("doublebeam", "getRawBeam", errstr.str(), 3);
  }
  if (bm == NULL)
    throw pofdExcept("doublebeam", "getRawBeam", "bm is not allocated", 4);

  if (n == 1) {
    bm[0] = 1.0;
    return;
  }

  double fwhm;
  if (band == 1) 
    fwhm = fwhm1;
  else if (band == 2)
    fwhm = fwhm2;
  else {
    std::stringstream errstr;
    errstr << "Invalid band selection " << band << "; should be 1 or 2";
    throw pofdExcept("doublebeam", "getBeamFac", errstr.str(), 5);
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
  \param[in] band Which band to use (1 or 2)
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
void doublebeam::getRawBeam(unsigned int band, unsigned int n, double pixsize, 
			    unsigned int oversamp, double* const bm) const {

  //Quick return
  if (oversamp == 1) {
    getRawBeam(band, n, pixsize, bm);
    return;
  }

  //Input checks
  if (n == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be positive";
    throw pofdExcept("doublebeam", "getRawBeam", errstr.str(), 1);
  }
  if (n % 2 == 0) {
    std::stringstream errstr;
    errstr << "n (" << n << ") should be odd";
    throw pofdExcept("doublebeam", "getRawBeam", errstr.str(), 2);
  }
  if (pixsize <= 0.0) {
    std::stringstream errstr;
    errstr << "pixsize (" << pixsize << ") should be positive";
    throw pofdExcept("doublebeam", "getRawBeam", errstr.str(), 3);
  }
  if (oversamp == 0) {
    std::stringstream errstr;
    errstr << "oversampling (" << oversamp << ") should be positive";
    throw pofdExcept("doublebeam", "getRawBeam", errstr.str(), 4);
  }
  if (oversamp % 2 == 0) {
    std::stringstream errstr;
    errstr << "oversampling (" << oversamp << ") should be odd";
    throw pofdExcept("doublebeam", "getRawBeam", errstr.str(), 5);
  }
  if (bm == NULL)
    throw pofdExcept("doublebeam", "getRawBeam", "bm is not allocated", 6);

  double fwhm;
  if (band == 1) 
    fwhm = fwhm1;
  else if (band == 2)
    fwhm = fwhm2;
  else {
    std::stringstream errstr;
    errstr << "Invalid band selection " << band << "; should be 1 or 2";
    throw pofdExcept("doublebeam", "getBeamFac", errstr.str(), 5);
  }

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
  \param[in] band Which band to use (1 or 2)
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[in] bm Beam. Must be pre-allocated by caller and be of length n * n
  \param[in] filter Hi-pass/Matched filter to apply.  If null, don't 
               apply filtering
*/
void doublebeam::getBeam(unsigned int band, unsigned int n, double pixsize, 
			 double* const bm, 
			 const fourierFilter* const filter) const {

  // pre-filtered beam
  getRawBeam(band, n, pixsize, bm);

  // Apply filtering
  if (filter != NULL)
    filter->filter(n, n, pixsize, bm);
}


/*!
  \param[in] band Which band to use (1 or 2)
  \param[in] n  Number of pixels along each dimension.  Should be odd
  \param[in] pixsize Size of pixels in arcsec
  \param[in] oversamp Oversampling.  Must be odd
  \param[in] bm Beam. Must be pre-allocated by
                   caller and be of length n * n
  \param[in] filter Hi-pass/matched filter to apply.  If null, don't apply 
                 filtering
*/
void doublebeam::getBeam(unsigned int band, unsigned int n, double pixsize, 
			 unsigned int oversamp, double* const bm, 
			 const fourierFilter* const filter) const {

  // Pre-filtered beam
  getRawBeam(band, n, pixsize, oversamp, bm);

  // Apply filtering
  if (filter != NULL)
    filter->filter(n, n, pixsize, bm);
}

/*!
  \param[in] outputfile Name of FITS file to write to
  \param[in] pixsize Pixel size, in arcseconds
  \param[in] nfwhm Number of FWHM out to go.
  \param[in] oversamp Oversampling to use.  Must be odd.  
  \param[in] filt1 Hi-pass/matched filter to apply, band 1, maybe band 2.  
              If null, don't apply filtering.
  \param[in] filt2 Hi-pass/matched filter to apply, band 2.  
              If null, don't apply filtering.
  \param[in] inverse Compute inverse beam rather than beam.

  If filt1 is set but filt2 is not, then filt1 is applied to both bands.
  There is no way to only filter one of the bands.
*/
void doublebeam::writeToFits(const std::string& outputfile, double pixsize, 
			     double nfwhm, unsigned int oversamp,
			     const fourierFilter* const filt1, 
			     const fourierFilter* const filt2, 
			     bool inverse) const {
  if (nfwhm <= 0.0)
    throw pofdExcept("doublebeam", "writeToFits", 
		     "Invalid (non-positive) nfwhm", 1);
  if (pixsize <= 0.0)
    throw pofdExcept("doublebeam", "writeToFits", 
		     "Invalid (non-positive) pixel size", 3);
  if (oversamp % 2 == 0)
    throw pofdExcept("doublebeam", "writeToFits", 
		     "Invalid (non-odd) oversampling", 4);

  // Figure out how many pixels we are going to use
  double fwhm = fwhm1 > fwhm2 ? fwhm1 : fwhm2;
  unsigned int npix = static_cast<unsigned int>(nfwhm * fwhm / pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;
  
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
  
  // Deal with filtering.
  const fourierFilter *f1, *f2;
  if (filt1 != NULL) {
    f1 = filt1;
    if (filt2 == NULL) f2 = filt1;
    else f2 = filt2;
  } else f1 = f2 = NULL;

  // Header for first band
  bool ubl;
  double dtmp;
  fits_write_key(fp, TLOGICAL, const_cast<char*>("INVERSE"), &inverse,
		 const_cast<char*>("Inverse beam?"), &status);
  dtmp = fwhm1;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM1"), &dtmp,
		 const_cast<char*>("FWHM band 1 [arcsec]"), &status);
  dtmp = pixsize;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("PIXSIZE"), &dtmp,
		 const_cast<char*>("Pixel scale [arcsec]"), &status);
  dtmp = nfwhm;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("NFWHM"), &dtmp,
		 const_cast<char*>("Number of FWHM out"), &status);
  if (f1 != NULL) {
    ubl = true;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("FILT"), &ubl,
		   const_cast<char*>("Beam iltered?"), &status);
    ubl = f1->isHipass();
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIPASS"), &ubl,
		   const_cast<char*>("Hipass filtered?"), &status);
    if (ubl) {
      dtmp = f1->getFiltScale();
      if (dtmp > 0)
	fits_write_key(fp, TDOUBLE, const_cast<char*>("FILTSCL"), &dtmp,
		       const_cast<char*>("Hipass filtering scale [arcsec]"), 
		       &status);
      dtmp = f1->getQFactor();
      if (dtmp > 0)
	fits_write_key(fp, TDOUBLE, const_cast<char*>("QFACTOR"), &dtmp,
		       const_cast<char*>("Hipass filtering apodization"), 
		       &status);
    }
    ubl = f1->isMatched();
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MATCHED"), &ubl,
		   const_cast<char*>("Match filtered?"), &status);
    if (ubl) {
      dtmp = f1->getFWHM();
      if (dtmp > 0)
	fits_write_key(fp, TDOUBLE, const_cast<char*>("MATFWHM"), &dtmp,
		       const_cast<char*>("Matched filtering fwhm [arcsec]"), 
		       &status);
      dtmp = f1->getSigInst();
      if (dtmp > 0)
	fits_write_key(fp, TDOUBLE, const_cast<char*>("MATSIGI"), &dtmp,
		       const_cast<char*>("Matched filtering sig_i [Jy]"), 
		       &status);
      dtmp = f1->getSigConf();
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
  
  // Get
  double *bmtmp = (double*) fftw_malloc(sizeof(double) * npix * npix);
  getBeam(1, npix, pixsize, oversamp, bmtmp, f1); // Also filters

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

  // Band 2 beam
  fits_create_img(fp, DOUBLE_IMG, 2, axissize, &status);
  fits_write_key(fp, TLOGICAL, const_cast<char*>("INVERSE"), &inverse,
		 const_cast<char*>("Inverse beam?"), &status);
  dtmp = fwhm2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM2"), &dtmp,
		 const_cast<char*>("FWHM band 2 [arcsec]"), &status);
  dtmp = pixsize;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("PIXSIZE"), &dtmp,
		 const_cast<char*>("Pixel scale [arcsec]"), &status);
  dtmp = nfwhm;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("NFWHM"), &dtmp,
		 const_cast<char*>("Number of FWHM out"), &status);

  if (f2 != NULL) {
    ubl = true;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("FILT"), &ubl,
		   const_cast<char*>("Filtered?"), &status);
    ubl = f2->isHipass();
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIPASS"), &ubl,
		   const_cast<char*>("Hipass filtered?"), &status);
    if (ubl) {
      dtmp = f2->getFiltScale();
      if (dtmp > 0)
	fits_write_key(fp, TDOUBLE, const_cast<char*>("FILTSCL"), &dtmp,
		       const_cast<char*>("Hipass filtering scale [arcsec]"), 
		       &status);
      dtmp = f2->getQFactor();
      if (dtmp > 0)
	fits_write_key(fp, TDOUBLE, const_cast<char*>("QFACTOR"), &dtmp,
		       const_cast<char*>("Hipass filtering apodization"), 
		       &status);
    }
    ubl = f2->isMatched();
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MATCHED"), &ubl,
		   const_cast<char*>("Band 1 match filtered?"), &status);
    if (ubl) {
      dtmp = f2->getFWHM();
      if (dtmp > 0)
	fits_write_key(fp, TDOUBLE, const_cast<char*>("MATFWHM"), &dtmp,
		       const_cast<char*>("Matched filtering fwhm [arcsec]"), 
		       &status);
      dtmp = f2->getSigInst();
      if (dtmp > 0)
	fits_write_key(fp, TDOUBLE, const_cast<char*>("MATSIGI"), &dtmp,
		       const_cast<char*>("Matched filtering1 sig_i [Jy]"), 
		       &status);
      dtmp = f2->getSigConf();
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
  getBeam(2, npix, pixsize, oversamp, bmtmp, f2); // Also filters
  for ( unsigned int j = 0; j < npix; ++j ) {
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

/////////////////////////////////////////////////////


/*!
  \param[in] NBINS Number of bins in histogram along each dimension
*/
doublebeamHist::doublebeamHist(unsigned int NBINS): 
  has_data(false), inverse(false), nbins(0), fwhm1(0.0), fwhm2(0.0),
  nfwhm(4.0), nfwhmkeep(4.0), pixsize(0.0), eff_area1(0.0), eff_area2(0.0), 
  oversamp(1) {

  if (NBINS == 0)
    throw pofdExcept("doublebeamHist", "doublebeamHist", 
		     "Invalid (non-positive) NBINS", 1);
  nbins = NBINS;
  for (unsigned int i = 0; i < 4; ++i) n[i] = 0;
  for (unsigned int i = 0; i < 4; ++i) wt[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) bm1[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) bm2[i] = NULL;

  dblpair NaN = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			       std::numeric_limits<double>::quiet_NaN());
  for (unsigned int i = 0; i < 4; ++i) minmax1[i] = NaN;
  for (unsigned int i = 0; i < 4; ++i) minmax2[i] = NaN;

  isHipass = isMatched = std::make_pair(false, false);
  filtscale = qfactor = matched_fwhm = matched_sigi = matched_sigc = NaN;
}

doublebeamHist::~doublebeamHist() {
  for (unsigned int i = 0; i < 4; ++i)
    if (wt[i] != NULL) delete[] wt[i];
  for (unsigned int i = 0; i < 4; ++i)
    if (bm1[i] != NULL) delete[] bm1[i];
  for (unsigned int i = 0; i < 4; ++i)
    if (bm2[i] != NULL) delete[] bm2[i];
}

/*!
  \param[in] bm Beam we are getting the histogram for
  \param[in] num_fwhm Number of FWHM out we will go on beam
  \param[in] pixsz Pixel size (arcseconds)
  \param[in] inv Histogram the inverse beam
  \param[in] oversampling Oversampling of beam. Must be odd
  \param[in] filt1 Fourier space filter to apply in band 1, and maybe band 2
  \param[in] filt2 Fourier space filter to apply in band 2
  \param[in] num_fwhm_keep How many FWHM to keep in the histogram
              after filtering.  If 0 (the default), keeps everything.

  If filt1 is set but not filt2, then filt1 is applied to both bands.
  There is no way to filter only one of the bands.
*/
void doublebeamHist::fill(const doublebeam& bm, double num_fwhm, double pixsz,
			  bool inv, unsigned int oversampling,
			  const fourierFilter* const filt1,
			  const fourierFilter* const filt2,
			  double num_fwhm_keep) {

  //Always ignore beam values below this
  // Use a lower value than in the 1D version because either beam can
  // fail this test
  const double minval = 5e-7; 

  // Clean out old arrays if present
  for (unsigned int i = 0; i < 4; ++i) n[i] = 0;
  for (unsigned int i = 0; i < 4; ++i)
    if (wt[i] != NULL) { delete[] wt[i]; wt[i] = NULL; }
  for (unsigned int i = 0; i < 4; ++i)
    if (bm1[i] != NULL) { delete[] bm1[i]; bm1[i] = NULL; }
  for (unsigned int i = 0; i < 4; ++i)
    if (bm2[i] != NULL) { delete[] bm2[i]; bm2[i] = NULL; }
  dblpair nan = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			       std::numeric_limits<double>::quiet_NaN());
  for (unsigned int i = 0; i < 4; ++i) minmax1[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) minmax2[i] = nan;

  dblpair tmp;
  inverse = inv;
  pixsize = pixsz;
  tmp = bm.getFWHM();
  fwhm1 = tmp.first;
  fwhm2 = tmp.second;
  oversamp = oversampling;
  nfwhm = num_fwhm;
  tmp = bm.getEffectiveArea();
  eff_area1 = tmp.first;
  eff_area2 = tmp.second;
  if (num_fwhm_keep == 0) nfwhmkeep = nfwhm;
  else nfwhmkeep = std::min(num_fwhm, num_fwhm_keep);

  // Get how many pixels we will go out
  double fwhm = std::max(fwhm1, fwhm2);
  unsigned int npix = static_cast<unsigned int>(nfwhm * fwhm / pixsize + 
						0.9999999999);
  npix = 2 * npix + 1;

  // Figure out filtering; we will apply filt1 to both if it is set
  //  but f2 is not.
  const fourierFilter *f1, *f2;
  if (filt1 != NULL) {
    f1 = filt1;
    if (filt2 != NULL) f2 = filt2; else f2 = filt1;
  } else if (filt2 != NULL) {
    f1 = NULL;
    f2 = filt2;
  } else f1 = f2 = NULL;

  // Get 2D beams
  // Temporary beam storage.  Must use fftw_malloc since we may filter
  double *bmtmp1 = (double*) fftw_malloc(sizeof(double) * npix * npix);
  double *bmtmp2 = (double*) fftw_malloc(sizeof(double) * npix * npix);
  // Get the beams
  bm.getBeam(1, npix, pixsize, oversamp, bmtmp1, f1); 
  bm.getBeam(2, npix, pixsize, oversamp, bmtmp2, f2); 

  unsigned int minidx, maxidx;
  if (nfwhmkeep < num_fwhm) {
    // We want to set up logical indexing into the array to only keep
    //  the part we want.
    unsigned int nclippix = 
      static_cast<unsigned int>(num_fwhm_keep * fwhm / pixsize + 0.9999999999);
    nclippix = 2 * nclippix + 1;
    if (nclippix < npix) {
      minidx = (npix - nclippix) / 2;
      maxidx = npix - minidx;
    } else {
      minidx = 0;
      maxidx = npix;
    }
  } else {
    minidx = 0;
    maxidx = npix;
  }
  
  // Store information about the filtering.  This is done purely
  //  in case we want to write the beam out.
  // Just reset all variables to initial state for simplicity first
  isHipass = std::make_pair(false, false);
  filtscale = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			     std::numeric_limits<double>::quiet_NaN());
  qfactor = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			   std::numeric_limits<double>::quiet_NaN());
  isMatched = std::make_pair(false, false);
  matched_fwhm = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN());
  matched_sigi = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN());
  matched_sigc = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN());
  if (f1 != NULL) {
    if (f1->isHipass()) {
      isHipass.first = true;
      filtscale.first = f1->getFiltScale();
      qfactor.first = f1->getQFactor();
    } 
    if (f1->isMatched()) {
      isMatched.first = true;
      matched_fwhm.first = f1->getFWHM();
      matched_sigi.first = f1->getSigInst();
      matched_sigc.first = f1->getSigConf();
    } 
  }
  if (f2 != NULL) {
    if (f2->isHipass()) {
      isHipass.second = true;
      filtscale.second = f2->getFiltScale();
      qfactor.second = f2->getQFactor();
    }
    if (f2->isMatched()) {
      isMatched.second = true;
      matched_fwhm.second = f2->getFWHM();
      matched_sigi.second = f2->getSigInst();
      matched_sigc.second = f2->getSigConf();
    }
  }

  // Histogram

  // Count up number of elements in each sign component, 
  // and do min/max values in each sign component.
  // The max beam value should certainly always be much smaller than 1e100
  // Also make a temporary nbins by nbins array that says what sign
  // component each thing will be.
  unsigned int ninbm[4];
  double minbinval1[4], maxbinval1[4];
  double minbinval2[4], maxbinval2[4];
  unsigned int npix2 = npix * npix;
  unsigned char* comparr;
  double *rowptr1, *rowptr2;
  comparr = new unsigned char[npix2];
  std::memset(ninbm, 0, 4 * sizeof(unsigned int));
  for (unsigned int i = 0; i < 4; ++i) minbinval1[i] = 1e100;
  for (unsigned int i = 0; i < 4; ++i) maxbinval1[i] = -1.0;
  for (unsigned int i = 0; i < 4; ++i) minbinval2[i] = 1e100;
  for (unsigned int i = 0; i < 4; ++i) maxbinval2[i] = -1.0;
  unsigned int comp;
  double val1, val2, fval1, fval2;
  for (unsigned int i = minidx; i < maxidx; ++i) {
    rowptr1 = bmtmp1 + i * npix;
    rowptr2 = bmtmp2 + i * npix;
    for (unsigned int j = minidx; j < maxidx; ++j) {
      val1 = rowptr1[j];
      val2 = rowptr2[j];
      fval1 = fabs(val1);
      fval2 = fabs(val2);
      // Ignore anything within [-minval, minval]
      if ((fval1 <= minval) || (fval2 <= minval)) {
	comparr[i*npix + j] = 4; //Invalid value -- will be skipped
	continue;
      }
      
      // Determine the sign component
      if (val1 > 0) {
	if (val2 > 0) comp = 0; else comp = 1;
      } else {
	if (val2 > 0) comp = 2; else comp = 3;
      }
      ++ninbm[comp];
      comparr[i*npix + j] = comp;
      
      // Min/max bit
      if (fval1 > maxbinval1[comp]) maxbinval1[comp] = fval1;
      else if (fval1 < minbinval1[comp]) minbinval1[comp] = fval1;
      if (fval2 > maxbinval2[comp]) maxbinval2[comp] = fval2;
      else if (fval2 < minbinval2[comp]) minbinval2[comp] = fval2;
    }
  }

  for (unsigned int i = 0; i < 4; ++i)
    if (ninbm[i] > 0) {
      minmax1[i] = std::make_pair(minbinval1[i], maxbinval1[i]);
      minmax2[i] = std::make_pair(minbinval2[i], maxbinval2[i]);
    }

  // Set bin sizes for each
  double ihiststep1[4], ihiststep2[4];
  const double log2outscale = 0.0014419741739063218; // log2(1.001)
  double dnbins = static_cast<double>(nbins);
  for (unsigned int i = 0; i < 4; ++i) 
    if (ninbm[i] > 0) {
      minbinval1[i] = log2(minbinval1[i]) - log2outscale;
      maxbinval1[i] = log2(maxbinval1[i]) + log2outscale;
      minbinval2[i] = log2(minbinval2[i]) - log2outscale;
      maxbinval2[i] = log2(maxbinval2[i]) + log2outscale;
      ihiststep1[i] = dnbins / (maxbinval1[i] - minbinval1[i]);
      ihiststep2[i] = dnbins / (maxbinval2[i] - minbinval2[i]);
    } else {
      ihiststep1[i] = 1.0;
      ihiststep2[i] = 1.0;
    }

  // Actually histogram, one component at a time
  unsigned int nbins2 = nbins * nbins;
  unsigned int curr_n, idx1, idx2, totidx, utmp;
  unsigned int *tmpwt, *wtptr;
  double wtval;
  double *tmphist1, *tmphist2, *bm1ptr, *bm2ptr;
  
  tmpwt = new unsigned int[nbins2];
  tmphist1 = new double[nbins2];
  tmphist2 = new double[nbins2];

  double min1, min2, istp1, istp2;
  for (unsigned int compidx = 0; compidx < 4; ++compidx) {
    curr_n = ninbm[compidx]; // Number in input beam, not histogram
    if (curr_n > 0) {
      // Zero temporary arrays
      std::memset(tmpwt, 0, nbins2 * sizeof(unsigned int));
      std::memset(tmphist1, 0, nbins2 * sizeof(double));
      std::memset(tmphist2, 0, nbins2 * sizeof(double));

      // Now loop over pixels
      min1 = minbinval1[compidx];
      min2 = minbinval2[compidx];
      istp1 = ihiststep1[compidx];
      istp2 = ihiststep2[compidx];
      for (unsigned int i = minidx; i < maxidx; ++i) {
	rowptr1 = bmtmp1 + i * npix;
	rowptr2 = bmtmp2 + i * npix;
	for (unsigned int j = minidx; j < maxidx; ++j) {
	  // Skip if not in this component
	  //  Also skips too-close-to-zeros
	  if (comparr[i * npix + j] != compidx) continue; 
	  fval1 = fabs(rowptr1[j]);
	  fval2 = fabs(rowptr2[j]);
	  idx1 = static_cast<unsigned int>((log2(fval1) - min1) * istp1);
	  idx2 = static_cast<unsigned int>((log2(fval2) - min2) * istp2);
	  totidx = idx1 * nbins + idx2;
	  tmpwt[totidx] += 1;
	  tmphist1[totidx] += fval1;
	  tmphist2[totidx] += fval2;
	}
      }

      // Count the number of histogram bins.
      unsigned int n_nonzero = 0;
      for (unsigned int i = 0; i < nbins2; ++i)
	if (tmpwt[i] > 0) ++n_nonzero; // Guaranteed to be > 0

      // Now we copy the bins to the final array.  Not necessary to zero
      wtptr = new unsigned int[n_nonzero];
      bm1ptr = new double[n_nonzero];
      bm2ptr = new double[n_nonzero];
      unsigned int binidx = 0;
      for (unsigned int i = 0; i < nbins2; ++i) {
	utmp = tmpwt[i];
	if (utmp > 0) {
	  wtptr[binidx] = utmp;
	  wtval = 1.0 / static_cast<double>(utmp);
	  bm1ptr[binidx] = tmphist1[i] * wtval;
	  bm2ptr[binidx] = tmphist2[i] * wtval;
	  ++binidx;
	}
      }

      if (inverse) {
	for (unsigned int i = 0; i < n_nonzero; ++i)
	  bm1ptr[i] = 1.0 / bm1ptr[i];
	for (unsigned int i = 0; i < n_nonzero; ++i)
	  bm2ptr[i] = 1.0 / bm2ptr[i];
      }

      n[compidx] = n_nonzero;
      wt[compidx] = wtptr;
      bm1[compidx] = bm1ptr;
      bm2[compidx] = bm2ptr;
    }
  }
      
  has_data = true;

  // Clean up
  fftw_free(bmtmp1);
  fftw_free(bmtmp2);
  delete[] comparr;
  delete[] tmpwt;
  delete[] tmphist1;
  delete[] tmphist2;
}

/*!
  \param[in] outputfile Name of FITS file to write to
*/
void doublebeamHist::writeToFits(const std::string& outputfile) const {

  if (!has_data)
    throw pofdExcept("doublebeamHist", "writeToFits", "Hist not filled", 1);

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
  dtmp = fwhm1;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM1"), &dtmp,
		 const_cast<char*>("FWHM1 [arcsec]"), &status);
  dtmp = fwhm2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM2"), &dtmp,
		 const_cast<char*>("FWHM1 [arcsec]"), &status);
  dtmp = pixsize;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("PIXSIZE"), &dtmp,
		 const_cast<char*>("Pixel scale [arcsec]"), &status);
  dtmp = nfwhm;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("NFWHM"), &dtmp,
		 const_cast<char*>("Number of FWHM out"), &status);
  dtmp = nfwhmkeep;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("NFWHMKP"), &dtmp,
		 const_cast<char*>("Number of FWHM kept"), &status);  

  if (isHipass.first) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT1"), &itmp,
		   const_cast<char*>("Is beam1 hipass filtered?"), &status);
    dtmp = filtscale.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTSCL1"), &dtmp,
		   const_cast<char*>("Filtering scale1 [arcsec]"), &status);
    dtmp = qfactor.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTQ1"), &dtmp,
		   const_cast<char*>("Filtering apodization1"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT1"), &itmp,
		   const_cast<char*>("Is beam1 hipass filtered?"), &status);
  }
  if (isHipass.second) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT2"), &itmp,
		   const_cast<char*>("Is beam2 hipass filtered?"), &status);
    dtmp = filtscale.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTSCL2"), &dtmp,
		   const_cast<char*>("Filtering scale2 [arcsec]"), &status);
    dtmp = qfactor.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTQ2"), &dtmp,
		   const_cast<char*>("Filtering apodization2"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT2"), &itmp,
		   const_cast<char*>("Is beam2 hipass filtered?"), &status);
  }
  if (isMatched.first) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT1"), &itmp,
		   const_cast<char*>("Is beam1 match filtered?"), &status);
    dtmp = matched_fwhm.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITFWHM"), &dtmp,
		   const_cast<char*>("Matched filtering1 FWHM [arcsec]"), 
		   &status);
    dtmp = matched_sigi.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGI1"), &dtmp,
		   const_cast<char*>("Matched filtering1 sigi"), &status);
    dtmp = matched_sigc.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGC1"), &dtmp,
		   const_cast<char*>("Matched filtering1 sigc"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT1"), &itmp,
		   const_cast<char*>("Is beam1 match filtered?"), &status);
  }
  if (isMatched.second) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT2"), &itmp,
		   const_cast<char*>("Is beam2 match filtered?"), &status);
    dtmp = matched_fwhm.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITFWHM"), &dtmp,
		   const_cast<char*>("Matched filtering2 FWHM [arcsec]"), 
		   &status);
    dtmp = matched_sigi.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGI2"), &dtmp,
		   const_cast<char*>("Matched filtering2 sigi"), &status);
    dtmp = matched_sigc.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGC2"), &dtmp,
		   const_cast<char*>("Matched filtering2 sigc"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT2"), &itmp,
		   const_cast<char*>("Is beam2 match filtered?"), &status);
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
  

  char* ttype[]= {"WEIGHT", "BEAM1", "BEAM2"};
  char* tform[] = {"1I", "1D", "1D"};

  char* label[4] = {"POSPOS", "POSNEG", "NEGPOS", "NEGNEG"};

  for (unsigned int idx=0; idx < 4; ++idx)
    if (n[idx] > 0) {
      fits_create_tbl(fp, BINARY_TBL, 0, 3, ttype, tform, NULL, label[idx],
		      &status );
      fits_insert_rows(fp, 0, n[idx], &status);
      for (unsigned int i = 0; i < n[idx]; ++i) {
	itmp = static_cast<int>(wt[idx][i]);
	fits_write_col(fp, TINT, 1, i+1, 1, 1, &itmp, &status);
	dtmp= bm1[idx][i];
	fits_write_col(fp, TDOUBLE, 2, i+1, 1, 1, &dtmp, &status);
	dtmp= bm2[idx][i];
	fits_write_col(fp, TDOUBLE, 3, i+1, 1, 1, &dtmp, &status);
      }
    }

  //Close up and go home
  fits_close_file(fp, &status);
  if (status) {
    fits_report_error(stderr, status);
    return;
  }
}
