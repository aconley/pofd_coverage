#include<sstream>
#include<cmath>

#include "../include/hipassFilter.h"
#include "../include/pofdExcept.h"

/*!
  \param[in] f Radius where filter is unity (in pixels)
  \param[in] q Sigma as a fraction of filtscale
 */
hipassFilter::hipassFilter(double f, double q) : filtscale(f), qfactor(q),
						 nx(0), ny(0), nyhalf(0), 
						 transdata(NULL), nxplan(0), 
						 nyplan(0), plan(NULL), 
						 plan_inv(NULL){}

hipassFilter::~hipassFilter() {
  if (transdata != NULL) fftw_free(transdata);
  if (plan != NULL) fftw_destroy_plan(plan);
  if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
}

/*!
  \param[in] NX size of input data array along 1st dimension
  \param[in] NY size of input data array along 2nd dimension
  \param[inout] data Data to filter.  Modified on output.  This must be
                      allocated with fftw_malloc.

  The filter is unity within filtscale, followed by a Gaussian decay
  with sigma qfactor * filtscale;
 */
void hipassFilter::filter(unsigned int NX, unsigned int NY, 
			  double* const data) {
			  

  const double nsig = 5.0; //!< Number of sigma out we go for Gaussian bit

  // Input checks
  if (NX * NY == 0) return; //Nothing to do
  if (filtscale <= 0.0) return;
  if (qfactor < 0) throw pofdExcept("hipassFilter", "filter",
				    "Invalid (negative) qfactor", 1);

  // If maximum radius is less than filtscale, we aren't going to do
  // anything.  So check for that.
  if (NX * NX + NY*NY < filtscale * filtscale) return;

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
  if (plan == NULL) { // New
    unsigned intx = static_cast<int>(NX);
    unsigned inty = static_cast<int>(NY);
    plan = fftw_plan_dft_r2c_2d(intx, inty, data, transdata, FFTW_ESTIMATE);
    plan_inv = fftw_plan_dft_c2r_2d(intx, inty, transdata, data,
				    FFTW_ESTIMATE);
    nxplan = NX;
    nyplan = NY;
  } else if ((NX != nxplan) || (NY != nyplan)) {
    fftw_destroy_plan(plan);
    fftw_destroy_plan(plan_inv);
    unsigned intx = static_cast<int>(NX);
    unsigned inty = static_cast<int>(NY);
    plan = fftw_plan_dft_r2c_2d(intx, inty, data, transdata, FFTW_ESTIMATE);
    plan_inv = fftw_plan_dft_c2r_2d(intx, inty, transdata, data,
				    FFTW_ESTIMATE);
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
  
  // Now, apply the filtering.  The filter is 1 within a radius
  // of filtscale (which should be in pixels), then a Gaussian lip
  // of sigma qmult * filtscale out to nsig sigmas, then 0
  unsigned int rowidx, iwrap;
  double kx2, cutdist2, dist2, set0dist2, sigfac, nsig2, expmult;
  sigfac = 1.0 / (filtscale * filtscale * qfactor * qfactor); // For computing nsig
  cutdist2 = filtscale * filtscale; // Inside this filter is 1
  set0dist2 = (1.0 + nsig * qfactor) * filtscale; // Outside this it's 0
  for (unsigned int i = 0; i < nx; ++i) {
    // Get k_x, minding the wrapping
    if (i <= nx / 2) iwrap = i; else iwrap = nx - i;
    kx2 = static_cast<double>(iwrap * iwrap);
    rowidx = i * nyhalf;
    if (kx2 > set0dist2) { 
      // Outside 0 radius, so will be zero for all y values
      for (unsigned int j = 0; j <= nyhalf; ++j) {
	transdata[rowidx + j][0] = 0.0;
	transdata[rowidx + j][1] = 0.0;
      }
    } else {
      // Have to check each y value
      for (unsigned int j = 0; j <= nyhalf; ++j) {
	dist2 = kx2 + static_cast<double>(j * j); // Distance from 0
	if (dist2 <= cutdist2) continue; // Filter is 1
	if (dist2 >= set0dist2) { // Filter is 0
	  transdata[rowidx + j][0] = 0.0;
	  transdata[rowidx + j][1] = 0.0;
	} else { // On Gaussian bit
	  nsig2 = dist2 - filtscale;
	  nsig2 = nsig2 * nsig2 * sigfac;
	  expmult = exp(-0.5 * nsig2);
	  transdata[rowidx + j][0] *= expmult;
	  transdata[rowidx + j][1] *= expmult;
	}
      }
    }
  }
  
  // Transform back -- same use of array execute versions
  fftw_execute_dft_c2r(plan_inv, transdata, data);

  // Rescale because forward followed by inverse FFTW transforms
  // scale array up by it's size
  double scalfac = 1.0 / (static_cast<double>(NX) * static_cast<double>(NY));
  for (unsigned int i = 0; i < NX * NY; ++i)
    data[i] *= scalfac;
}
