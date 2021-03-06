#include<sstream>
#include<cmath>
#include<limits>

#include "../include/fourierFilter.h"
#include "../include/pofdExcept.h"

double NaN = std::numeric_limits<double>::quiet_NaN();

/*!
  \param[in] pixsize Size of pixels, in arcseconds
  \param[in] FWHM FWHM of Gaussian beam, in arcseconds
  \param[in] sigi Instrument noise, in Jy
  \param[in] sigc Confusion noise, in Jy
  \param[in] quickfft If set, use FFTW_ESTIMATE for the plans.  Otherwise
                    use FFTW_MEASURE.  Set this if you are only planning
                    on calling this once.
  \param[in] fixedsize Only allow this to be resized the first time.

  This version of the constructor sets up only matched filtering
*/
fourierFilter::fourierFilter(double pixsize, double FWHM, double sigi,
                             double sigc, bool quickfft, bool fixedsize):
  initialized(false), doHipass(false), doMatched(true),
  allowResize(!fixedsize), nx(0), ny(0), pixscale(pixsize), filtscale(NaN),
  qfactor(NaN),  fwhm(FWHM), sig_inst(sigi), sig_conf(sigc), nyhalf(0),
  plan(nullptr), plan_inv(nullptr), filt_fft(nullptr), map_fft(nullptr) {

  // Check inputs
  if (fwhm <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) FWHM", 3);
  if (pixsize <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) pixsize", 4);
  if (fwhm / pixsize < 2.0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Insufficiently sampled FWHM", 5);
  if (sigi <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) sigi", 6);
  if (sigc <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) sigc", 7);


  // Set plan style.  We don't support WISDOM here because our image
  //  sizes are unlikely to be regular enough to be prepared.
  if (quickfft) fftw_plan_style = FFTW_ESTIMATE; 
  else fftw_plan_style = FFTW_MEASURE;
}

/*!
  \param[in] pixsize Size of pixels, in arcseconds
  \param[in] fscale Filtering scale, in arcseconds
  \param[in] q Gaussian sigma of apodization as a fraction of fscale
  \param[in] quickfft If set, use FFTW_ESTIMATE for the plans.  Otherwise
                    use FFTW_MEASURE.  Set this if you are only planning
                    on calling this once.
  \param[in] fixedsize Only allow this to be resized the first time.

  This version of the constructor sets up only hipass filtering
*/
fourierFilter::fourierFilter(double pixsize, double fscale, double q,
                             bool quickfft, bool fixedsize):
  initialized(false), doHipass(true), doMatched(false), allowResize(!fixedsize),
  nx(0), ny(0), pixscale(pixsize), filtscale(fscale), qfactor(q),
  fwhm(NaN), sig_inst(NaN), sig_conf(NaN), nyhalf(0),
  plan(nullptr), plan_inv(nullptr), filt_fft(nullptr), map_fft(nullptr) {

  // Check inputs
  if (filtscale <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) filtscale", 3);
  if (pixsize <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) pixsize", 4);
  if (q < 0.0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) q", 5);


  // Set plan style.  We don't support WISDOM here because our image
  //  sizes are unlikely to be regular enough to be prepared.
  if (quickfft) fftw_plan_style = FFTW_ESTIMATE; 
  else fftw_plan_style = FFTW_MEASURE;
}


/*!
  \param[in] pixsize Size of pixels, in arcseconds
  \param[in] FWHM FWHM of Gaussian beam, in arcseconds
  \param[in] sigi Instrument noise, in Jy
  \param[in] sigc Confusion noise, in Jy
  \param[in] fscale Filtering scale, in arcseconds
  \param[in] q Gaussian sigma of apodization as a fraction of fscale
  \param[in] quickfft If set, use FFTW_ESTIMATE for the plans.  Otherwise
                    use FFTW_MEASURE.  Set this if you are only planning
                    on calling this once.
  \param[in] fixedsize Only allow this to be resized the first time.

  This version of the constructor sets up both hipass and matched filtering
*/
fourierFilter::fourierFilter(double pixsize, double FWHM, double sigi,
                             double sigc, double fscale, double q,
                             bool quickfft, bool fixedsize):
  initialized(false), doHipass(true), doMatched(true), allowResize(!fixedsize),
  nx(0), ny(0), pixscale(pixsize), filtscale(fscale), qfactor(q), 
  fwhm(FWHM), sig_inst(sigi), sig_conf(sigc), nyhalf(0),
  plan(nullptr), plan_inv(nullptr), filt_fft(nullptr), map_fft(nullptr) {

  // Check inputs
  if (fwhm <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) FWHM", 3);
  if (pixsize <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) pixsize", 4);
  if (fwhm / pixsize < 2.0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Insufficiently sampled FWHM", 5);
  if (sigi <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) sigi", 6);
  if (sigc <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) sigc", 7);
  if (filtscale <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) filtscale", 8);
  if (pixsize <= 0)
    throw pofdExcept("fourierFilter", "fourierFilter", 
                     "Invalid (non-positive) pixsize", 9);

  // Set plan style.  We don't support WISDOM here because our image
  //  sizes are unlikely to be regular enough to be prepared.
  if (quickfft) fftw_plan_style = FFTW_ESTIMATE; 
  else fftw_plan_style = FFTW_MEASURE;
}

fourierFilter::~fourierFilter() {
  if (filt_fft != nullptr) fftw_free(filt_fft);
  if (map_fft != nullptr) fftw_free(map_fft);
  if (plan != nullptr) fftw_destroy_plan(plan);
  if (plan_inv != nullptr) fftw_destroy_plan(plan_inv);
}

/*!
  Set up the beam for fourierFilter::setup_matched
*/
void fourierFilter::setup_beam(double* const bm) const {
  double bmsig = pofd_coverage::fwhmfac * fwhm / pixscale;
  double expfac = -0.5 / (bmsig * bmsig);
  
  // Do the beam in quadrants (because it wraps)
  double distsq, xdsq;
  double *rowptr;
  for (unsigned int i = 0; i < nx / 2; ++i) {
    rowptr = bm + i * ny;
    xdsq = static_cast<double>(i * i);
    for (unsigned int j = 0; j < ny / 2; ++j) {
      distsq = xdsq + j * j;
      rowptr[j] = exp(expfac * distsq);
    }
    for (unsigned int j = ny / 2; j < ny; ++j) {
      distsq = xdsq + static_cast<double>((ny - j) * (ny - j));
      rowptr[j] = exp(expfac * distsq);
    }
  }
  for (unsigned int i = nx / 2; i < nx; ++i) {
    rowptr = bm + i * ny;
    xdsq = static_cast<double>((nx - i) * (nx - i));
    for (unsigned int j = 0; j < ny / 2; ++j) {
      distsq = xdsq + j * j;
      rowptr[j] = exp(expfac * distsq);
    }
    for (unsigned int j = ny / 2; j < ny; ++j) {
      distsq = xdsq + static_cast<double>((ny - j) * (ny - j));
      rowptr[j] = exp(expfac * distsq);
    }
  }
}

/*!
  \param[in] rl Real variable (size nx by ny) to plan with.
  \param[in] im Complex variable (size nx by ny/2+1) to plan with

  Note that both of these will be destroyed on output for any
  style of planning except FFTW_ESTIMATE
*/
void fourierFilter::setup_plans(double* const rl,
                                fftw_complex* const im) const {

  unsigned intx = static_cast<int>(nx);
  unsigned inty = static_cast<int>(ny);
  if (plan != nullptr) fftw_destroy_plan(plan);
  plan = fftw_plan_dft_r2c_2d(intx, inty, rl, im, fftw_plan_style);
  if (plan == nullptr) {
    std::stringstream str;
    str << "Plan creation failed for forward transform of size: " << 
      nx << " by " << ny;
    throw pofdExcept("fourierFilter", "setup_plans", str.str(), 1);
  }
  if (plan_inv != nullptr) fftw_destroy_plan(plan_inv);

  plan_inv = fftw_plan_dft_c2r_2d(intx, inty, im, rl, fftw_plan_style);
  if (plan_inv == nullptr) {
    std::stringstream str;
    str << "Inverse plan creation failed for forward transform of size: " << 
      nx << " by " << ny;
    throw pofdExcept("fourierFilter", "setup_plans", str.str(), 2);
  }
  initialized = false;
}

/*! 
  Set up the filter for matched filtering.  Only touches mutable variables.
*/
void fourierFilter::setup_matched() const {
  if (initialized) return; 

  // Allocate internal arrays.  They can't have been touched
  filt_fft = (fftw_complex*) fftw_malloc(nx * nyhalf * sizeof(fftw_complex));
  map_fft = (fftw_complex*) fftw_malloc(nx * nyhalf * sizeof(fftw_complex));

  // We will need some auxilliary storage.  A real space array
  //  to hold the beam.  We will allocate that on the heap because 
  //  this function is only called once.  Do this first so we can use
  //  it to plan
  double *bm;
  bm = (double*) fftw_malloc(nx * ny * sizeof(double));
  
  // Set up plans
  setup_plans(bm, map_fft);

  // Now set up the beam
  setup_beam(bm);
  
  // Get the variance of the PSF using the two-pass algorithm
  double mn, nxny;
  nxny = static_cast<double>(nx * ny);
  mn = bm[0];
  for (unsigned int i = 1; i < nx * ny; ++i)
    mn += bm[i];
  mn /= nxny;
  double sum1 = 0.0, sum2 = 0.0, val;
  for (unsigned int i = 0; i < nx * ny; ++i) {
    val = bm[i] - mn;
    sum1 += val * val;
    sum2 += val;
  }
  double var_bm = (sum1 - sum2 * sum2 / nxny) / (nxny - 1.0);
  
  // FFT the beam into map_fft
  fftw_execute_dft_r2c(plan, bm, map_fft);

  // Now build the (unscaled, conjugated) filter in filt_fft, which is
  //   conj(map_fft) / (nx * ny * sigi^2 + (sig_c/sig_psf)^2 map_fft^2
  double inst_term = sig_inst * sig_inst * nxny;
  double conf_scale = sig_conf * sig_conf / var_bm;
  double rmap, imap, pnoise;
  for (unsigned int i = 0; i < nx * nyhalf; ++i) {
    rmap = map_fft[i][0]; // This is the FFTed beam
    imap = map_fft[i][1];
    pnoise = inst_term + conf_scale * (rmap * rmap + imap * imap);
    filt_fft[i][0] = rmap / pnoise;
    filt_fft[i][1] = -imap / pnoise;
  }

  // Now we compute the scaling by applying the filter to the original
  //  beam (modifying map_fft) and then transforming back into bm.
  //  Since we are per-beam normalized, and setup_beam produces a map
  //  with a peak value of unity, the scaling is just one over
  //  the maximum.  And the maximum should be where we put it -- in the 0
  //  pixel.
  double rfilt, ifilt;
  for (unsigned int i = 0; i < nx * nyhalf; ++i) {
    rmap = map_fft[i][0]; // Still the FFTed beam -- until we modify
    rfilt = filt_fft[i][0]; // it here
    imap = map_fft[i][1];
    ifilt = filt_fft[i][1];
    map_fft[i][0] = rmap * rfilt - imap * ifilt;
    map_fft[i][1] = imap * rfilt + rmap * ifilt;
  }
  fftw_execute_dft_c2r(plan_inv, map_fft, bm);

  // The matched filter shouldn't move the maximum, which was in the
  //  0th pixel.
  double maxval = bm[0];

  if (maxval <= 0.0)
    throw pofdExcept("fourierFilter", "setup",
                     "Invalid recovered maximum value", 3);

  double scalefac = 1.0 / maxval;
  for (unsigned int i = 0; i < nx * nyhalf; ++i) {
    filt_fft[i][0] *= scalefac;
    filt_fft[i][1] *= scalefac;
  }
    
  // Done and clean up
  fftw_free(bm);

  initialized = true;
}

/*! 
  Set up the filter for hipass filtering.  Only touches mutable variables.
*/
void fourierFilter::setup_hipass() const {
  if (initialized) return;
  
  // This one is quite simple -- we just need to allocate
  //  map_fft because hipass filtering doesn't require setting up filt_fft,
  //  and set up the plans
  // Note we can do this even if setup_hipass has been called
  //  because of the quick return on initialized above.
  map_fft = (fftw_complex*) fftw_malloc(nx * nyhalf * sizeof(fftw_complex));

  // Need some temporary planning storage -- even FFTW_ESTIMATE
  // doesn't work correctly without something to plan through
  double *tmp;
  tmp = (double*) fftw_malloc(nx * ny * sizeof(double));
  setup_plans(tmp, map_fft);
  fftw_free(tmp);

  initialized = true;
}

/*!
  \param[in] NX New size of transform, dimension 1
  \param[in] NY New size of transform, dimension 2
  
  \returns true if resizing was done.
  Set up all filtering
*/
bool fourierFilter::setup(unsigned int NX, unsigned int NY) const {

  if (NX < 2)
    throw pofdExcept("fourierFilter", "setup", "Invalid NX (<2)", 1);
  if (NY < 2)
    throw pofdExcept("fourierFilter", "setup", "Invalid NY (<2)", 2);

  // Figure out if the size has changed
  bool resized = false;
  if ((NX != nx) || (NY != ny)) {
    if ((!allowResize) && (nx != 0))
      throw pofdExcept("fourierFilter", "filter",
                       "Can only resize the first time this is called", 1);
    // Have to resize; clean up and mark as uninitialized.
    if (filt_fft != nullptr) { fftw_free(filt_fft); filt_fft = nullptr; }
    if (map_fft != nullptr) { fftw_free(map_fft); map_fft = nullptr; }
    if (plan != nullptr) {
      fftw_destroy_plan(plan);
      plan = nullptr;
    }
    if (plan_inv != nullptr) {
      fftw_destroy_plan(plan_inv);
      plan_inv = nullptr;
    }
    initialized = false;
    resized = true;
    nx = NX;
    ny = NY;
    nyhalf = ny / 2 + 1;
  }

  if (doMatched) setup_matched();
  if (doHipass) setup_hipass();
  return resized;
}

/*!
  \param[inout] data Data to mean subtract; modified on output
  \returns The mean of the data

  Mean subtracts the input data
*/
double fourierFilter::meanSub(double* const data) const {
  double tot;
  tot = data[0];
  for (unsigned int i = 1; i < nx * ny; ++i)
    tot += data[i];
  double meanval = tot / static_cast<double>(nx * ny);
  for (unsigned int i = 0; i < nx * ny; ++i) data[i] -= meanval;
  return meanval;
}


/*!
  \param[in] n1 Dimension one extent of data
  \param[in] n2 Dimension two extent of data
  \param[in] pixsize Pixel scale of data, in arcsec
  \param[inout] data Data to filter.  Modified on output.  This must be
                      allocated with fftw_malloc. 

  n1, n2, and pixscale are just used to check against what this was set up 
  with.  If they don't match, the code throws an exception.

  The data must have the same pixel size as the filter was set up with,
  but the other variables do not need to match.

  The order is matched, then hipass.  The matched filter maintains
  peak (per-beam) normalization for isolated objects, but the high-pass
  filter does not.  In general, the high-pass filter will mildly suppress
  the peak values of sources more than it does white noise -- which means 
  that (unsurprisingly), by throwing out low frequency information some 
  signal-to-noise will be lost.

  This will resize the filter internally if needed, but that is an expensive 
  operation.  So the caller should do their best to call this
  with the same n1, n2 every time to avoid setup overheads.
*/
void fourierFilter::filter(unsigned int n1, unsigned int n2, double pixsize,
                           double* const data) const {
  const double nsig = 5.0; //Number of sigma out we go out in highpass Gaussian

  double reldiff = fabs(pixsize - pixscale) / pixscale;
  if (reldiff > 1e-4) {
    std::stringstream errstr;
    errstr << "Pixel size (" << pixsize << ") doesn't match what this filter"
           << " was set up with (" << pixscale << ")";
    throw pofdExcept("fourierFilter", "filter", errstr.str(), 1);
  }

  // Handle the special case of no matched filtering, and hipass filtering
  // on scales larger than the map.  That's just mean subtracting.
  // Note this can be done without calling setup!
  double n1n2 = static_cast<double>(n1 * n2);
  if (doHipass && (!doMatched)) {
    double filtscale_pix = filtscale / pixscale;
    if (filtscale_pix * filtscale_pix > n1n2) {
      meanSub(data);
      return;
    }
  }

  // If we got this far, we need to allocate some stuff, maybe 
  //  set up the matched filter.  This will set nx=n1, ny=n2.
  // If this is true already and we have done this before, this will
  //  quick return.
  setup(n1, n2);
  double nxny = static_cast<double>(nx * ny);

  // Forward transform the data into map_fft
  // We use the advanced 'array execute' interface in case we are 
  // re-using an old plan which may have been set up on a different
  // data array.
  fftw_execute_dft_r2c(plan, data, map_fft);

  if (doMatched) {
    // Multiply the filter into the map.  Note that we actually store 
    // the conjugated fft of the filter so that this is cross-correlation 
    // rather than convolution, which is what we want
    double rmap, rfilt, imap, ifilt;
    for (unsigned int i = 0; i < nx * nyhalf; ++i) {
      rmap = map_fft[i][0];
      rfilt = filt_fft[i][0];
      imap = map_fft[i][1];
      ifilt = filt_fft[i][1];
      map_fft[i][0] = rmap * rfilt - imap * ifilt;
      map_fft[i][1] = imap * rfilt + rmap * ifilt;
    }
  }

  if (doHipass) {
    // We do this on the fly, since it's more efficient.
    // We are high-pass filtering, so we want to remove components below 
    //  a certain wavenumber.
    // The filter is 1 outside a radius filtscale (in pixels),
    //  and inside that there is a Gaussian lip down to
    //  of sigma qmult * filtscale_pix out to nsig sigmas, then 0
    // The mean is always set to zero.
    // It is more convenient to work with the index rather than carrying
    //  along all the scaling factors.  This works out such that the
    //  filter is 1 for all i^2 + j^2 > nx * ny / filtscale_pix^2 
    unsigned int rowidx, iwrap;
    double kx2, set1dist2, dist2;
    double filtscale_pix = filtscale / pixscale;
  
    // If i^2 + j^2 > than this, filter is unity
    set1dist2 = nxny / (filtscale_pix * filtscale_pix);

    if (qfactor > 0) {
      // Outside set1dist, the filter is 1
      // Inside set0dist, the filter is 0
      // Between them it's Gaussian.  We work with the squares 
      //  when possible.
      double qpixfac, set0dist2, set1dist, reldist, sigfac, expmult;
      qpixfac = 1.0 - qfactor * nsig; // How far down we go before setting to 0
      if (qpixfac <= 0) set0dist2 = 0.0; // Gaussian all the way down
      else set0dist2 = qpixfac * qpixfac * set1dist2;
      set1dist = sqrt(set1dist2);
      sigfac = -0.5 / (qfactor * qfactor * set1dist2); // -1/2*sigma^2
      for (unsigned int i = 0; i < nx; ++i) {
        if (i <= nx / 2) iwrap = i; else iwrap = nx - i;
        kx2 = static_cast<double>(iwrap * iwrap);
        if (kx2 > set1dist2) continue; // Filter is always 1 for this i
        // If that didn't already continue, we have to check j by j
        rowidx = i * nyhalf;
        for (unsigned int j = 0; j < nyhalf; ++j) {
          dist2 = kx2 + static_cast<double>(j * j); // Distance2 from 0
          if (dist2 > set1dist2) {} //Filter is 1
          else if (dist2 <= set0dist2) { // Filter is 0
            map_fft[rowidx + j][0] = 0.0;
            map_fft[rowidx + j][1] = 0.0;
          } else {
            // On the Gaussian.  This is the messiest case
            reldist = set1dist - sqrt(dist2); // distance from edge of 1 region
            expmult = exp(sigfac * reldist * reldist);
            map_fft[rowidx + j][0] *= expmult;
            map_fft[rowidx + j][1] *= expmult;
          }
        }
      }
    } else {
      // No gaussian apodization
      // 1 outside this radius^2, 0 inside
      for (unsigned int i = 0; i < nx; ++i) {
        if (i <= nx / 2) iwrap = i; else iwrap = nx - i;
        kx2 = static_cast<double>(iwrap * iwrap);
        if (kx2 > set1dist2) continue;
        rowidx = i * nyhalf;
        for (unsigned int j = 0; j < nyhalf; ++j) {
          dist2 = kx2 + static_cast<double>(j * j); // Distance from 0
          if (dist2 <= set1dist2) { // Filter is 0
            map_fft[rowidx + j][0] = 0.0;
            map_fft[rowidx + j][1] = 0.0;
          }
        }
      }
    }
    // Since FFTW does unscaled transforms, we need to fix the scaling. 
    // However, if we also did a matched filter, we don't want to do 
    //  it again.
    if (!doMatched) {
      double scalfac = 1.0 / nxny;
      for (unsigned int i = 1; i < nx * nyhalf; ++i) {
        map_fft[i][0] *= scalfac;
        map_fft[i][1] *= scalfac;
      }
    }
  }

  // Always 0 mean
  map_fft[0][0] = 0.0;
  map_fft[0][1] = 0.0;

  // Transform back -- same use of array execute versions
  // We don't have to rescale because that's already included
  //  in the filtering operations above
  fftw_execute_dft_c2r(plan_inv, map_fft, data);
}
