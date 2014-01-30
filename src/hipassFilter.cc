#include<sstream>
#include<cmath>
#include<limits>

#include "../include/hipassFilter.h"
#include "../include/pofdExcept.h"

/*!
  \param[in] f Radius outside which filter is unity (in pixels)
  \param[in] q Sigma as a fraction of filtscale if apodizing with Gaussian
  \param[in] quickfft If set, use FFTW_ESTIMATE for the plans.  Otherwise
                    use FFTW_MEASURE.  Set this if you are only planning
                    on calling this once.
 */
hipassFilter::hipassFilter(double f, double q, bool quickfft) : 
  filtscale(f), qfactor(q), nx(0), ny(0), nyhalf(0), 
  transdata(NULL), nxplan(0), nyplan(0), plan(NULL), plan_inv(NULL) {

  // Set plan style.  We don't support WISDOM here because our image
  //  sizes are unlikely to be regular enough to be prepared.
  if (quickfft) fftw_plan_style = FFTW_ESTIMATE; 
  else fftw_plan_style = FFTW_MEASURE;

}


hipassFilter::~hipassFilter() {
  if (transdata != NULL) fftw_free(transdata);
  if (plan != NULL) fftw_destroy_plan(plan);
  if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
}

double hipassFilter::meanSub(unsigned int N, double* const data) const {
  if (N == 0) return std::numeric_limits<double>::quiet_NaN();
  double tot;
  tot = data[0];
  for (unsigned int i = 1; i < N; ++i)
    tot += data[i];
  double meanval = tot / static_cast<double>(N);
  for (unsigned int i = 0; i < N; ++i) data[i] -= meanval;
  return meanval;
}

/*!
  \param[in] pixscale Pixel scale in arcseconds.
  \param[in] NX size of input data array along 1st dimension
  \param[in] NY size of input data array along 2nd dimension
  \param[inout] data Data to filter.  Modified on output.  This must be
                      allocated with fftw_malloc.
  In wavenumber space, the filter is unity for wavenumbers larger
  than the filterscale, followed by a Gaussian decay with 
  sigma=qfactor * filtscale.
 */
void hipassFilter::filter(double pixscale, unsigned int NX, unsigned int NY, 
			  double* const data) {
			  

  const double nsig = 5.0; //!< Number of sigma out we in for Gaussian bit

  // Input checks
  if (qfactor < 0) throw pofdExcept("hipassFilter", "filter",
				    "Invalid (negative) qfactor", 1);
  if (NX * NY == 0) return; //Nothing to do
  if (filtscale < 0.0) return; // Do nothing
  if (filtscale == 0) {
    // Mean sub only
    meanSub(NX * NY, data);
    return;
  }

  double filtscale_pix = filtscale / pixscale;

  // More detailed check if we will be doing any filtering.
  // This follows from the fact that the filter is 1 (i.e., it does nothing)
  // for all radius^2 (in pix) > nx * ny / filtscale_pix^2.  So if the
  // smallest non-zero value (1) is larger than this, all we are doing is
  // setting the filter at frequency 0 to 0.
  double nxny = static_cast<double>(NX * NY);
  if (filtscale_pix * filtscale_pix > nxny) {
    // Just mean sub
    meanSub(NX * NY, data);
    return;
  }

  // Allocate space if needed
  if (transdata == NULL) {
    // First alloc
    nx = NX;
    ny = NY;
    nyhalf = ny / 2 + 1;
    transdata = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * nyhalf);
  } else if (NX > nx || NY > ny) {
    fftw_free(transdata);
    nx = NX;
    ny = NY;
    nyhalf = ny / 2 + 1;
    transdata = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * nyhalf);
  }

  // Deal with plans
  // We use FFTW measure because frequently we are doing this many times
  if (plan == NULL) { // New
    unsigned intx = static_cast<int>(NX);
    unsigned inty = static_cast<int>(NY);
    plan = fftw_plan_dft_r2c_2d(intx, inty, data, transdata, fftw_plan_style);
    plan_inv = fftw_plan_dft_c2r_2d(intx, inty, transdata, data,
				    fftw_plan_style);
    nxplan = NX;
    nyplan = NY;
  } else if ((NX != nxplan) || (NY != nyplan)) {
    fftw_destroy_plan(plan);
    fftw_destroy_plan(plan_inv);
    unsigned intx = static_cast<int>(NX);
    unsigned inty = static_cast<int>(NY);
    plan = fftw_plan_dft_r2c_2d(intx, inty, data, transdata, fftw_plan_style);
    plan_inv = fftw_plan_dft_c2r_2d(intx, inty, transdata, data,
				    fftw_plan_style);
    nxplan = NX;
    nyplan = NY;
  }
  if (plan == NULL) {
    std::stringstream str;
    str << "Plan creation failed for forward transform of size: " << 
      NX << " by " << NY;
    throw pofdExcept("hipassFilter", "filter", str.str(), 2);
  }
  if (plan_inv == NULL) {
    std::stringstream str;
    str << "Reverse plan creation failed for forward transform of size: " << 
      NX << " by " << NY;
    throw pofdExcept("hipassFilter", "filter", str.str(), 3);
  }

  // Forward transform the data into transdata
  // We use the advanced 'array execute' interface in case we are 
  // re-using an old plan which may have been set up on a different
  // input array (but with the same dimensions)
  fftw_execute_dft_r2c(plan, data, transdata);
  
  // Now, apply the filtering.  We are high-pass filtering, so we
  // want to remove components below a certain wavenumber.
  // The filter is 1 outside a radius filtscale (in pixels),
  // and inside that there is a Gaussian lip down to
  // of sigma qmult * filtscale_pix out to nsig sigmas, then 0
  // The mean is always set to zero.
  // It is more convenient to work with the index rather than carrying
  // along all the scaling factors.  This works out such that the
  // filter is 1 for all i^2 + j^2 > nx * ny / filtscale_pix^2 

  unsigned int rowidx, iwrap;
  double kx2, set1dist2, dist2;
  
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
      for (unsigned int j = 0; j <= nyhalf; ++j) {
	dist2 = kx2 + static_cast<double>(j * j); // Distance2 from 0
	if (dist2 > set1dist2) {} //Filter is 1
	else if (dist2 <= set0dist2) { // Filter is 0
	  transdata[rowidx + j][0] = 0.0;
	  transdata[rowidx + j][1] = 0.0;
	} else {
	  // On the Gaussian.  This is the messiest case
	  reldist = set1dist - sqrt(dist2); // distance from edge of 1 region
	  expmult = exp(sigfac * reldist * reldist);
	  transdata[rowidx + j][0] *= expmult;
	  transdata[rowidx + j][1] *= expmult;
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
      for (unsigned int j = 0; j <= nyhalf; ++j) {
	dist2 = kx2 + static_cast<double>(j * j); // Distance from 0
	if (dist2 <= set1dist2) { // Filter is 0
	  transdata[rowidx + j][0] = 0.0;
	  transdata[rowidx + j][1] = 0.0;
	}
      }
    }
  }
  // Always 0 mean
  transdata[0][0] = 0.0;
  transdata[0][1] = 0.0;
  
  // Transform back -- same use of array execute versions
  fftw_execute_dft_c2r(plan_inv, transdata, data);

  // Rescale because forward followed by inverse FFTW transforms
  // scale array up by it's size
  double scalfac = 1.0 / (static_cast<double>(NX) * static_cast<double>(NY));
  for (unsigned int i = 0; i < NX * NY; ++i)
    data[i] *= scalfac;
}
