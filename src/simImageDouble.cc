#include<ctime>
#include<cmath>
#include<limits>
#include<cstring>

#include<fitsio.h>
#include<fftw3.h> // For fftw_malloc

#include "../include/global_settings.h"
#include "../include/simImageDouble.h"
#include "../include/pofdExcept.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

/*!
  \param[in] N1 Dimension of simulated image along 1st axis
  \param[in] N2 Dimension of simulated image along 2nd axis
  \param[in] PIXSIZE Pixel size in arcseconds.
  \param[in] FWHM1 The FWHM of the simulated image in arcsec, before
              filtering or smoothing, band 1
  \param[in] FWHM2 The FWHM of the simulated image in arcsec, before
              filtering or smoothing, band 2
  \param[in] SIGI1 Gaussian instrumental noise in the units the model is in,
              band 1. This is the noise before additional smoothing or 
	      filtering.
  \param[in] SIGI2 Gaussian instrumental noise in the units the model is in,
              band 2. This is the noise before additional smoothing or 
	      filtering.
  \param[in] ESMOOTH1 Additional Gaussian smoothing FWHM, band 1.  If non-zero,
              additional smoothing is applied.
  \param[in] ESMOOTH2 Additional Gaussian smoothing FWHM, band 2.  If non-zero,
              additional smoothing is applied.
  \param[in] OVERSAMPLE Oversampling of simulated image.  Must be odd.
              1 means no additional oversampling.
  \param[in] NBINS Number of bins, if binning is applied.
  \param[in] powerspecfile File containing power spectrum.  If set, the
              source positions are generated using this P(k).  If not set, they
	      are uniformly distributed across the image.
*/
simImageDouble::simImageDouble(unsigned int N1, unsigned int N2, double PIXSIZE,
			       double FWHM1, double FWHM2, double SIGI1, 
			       double SIGI2, double ESMOOTH1, double ESMOOTH2,
			       unsigned int OVERSAMPLE, unsigned int NBINS,
			       const std::string& powerspecfile) {

  if (N1 == 0)
    throw pofdExcept("simImageDouble", "simImageDouble", 
		     "Invalid (non-positive) N1", 1);
  if (N2 == 0)
    throw pofdExcept("simImageDouble", "simImageDouble", 
		     "Invalid (non-positive) N2", 2);
  if (PIXSIZE <= 0.0)
    throw pofdExcept("simImageDouble", "simImageDouble", 
		     "Invalid (non-positive) PIXSIZE", 3);
  if (FWHM1 <= 0.0)
    throw pofdExcept("simImageDouble", "simImageDouble", 
		     "Invalid (non-positive) FWHM1", 4);
  if (FWHM2 <= 0.0)
    throw pofdExcept("simImageDouble", "simImageDouble", 
		     "Invalid (non-positive) FWHM2", 5);
  if (OVERSAMPLE % 2 == 0)
    throw pofdExcept("simImageDouble", "simImageDouble", 
		     "Invalid (even) OVERSAMPLE", 6);
  
  n1 = N1;
  n2 = N2;
  oversample = OVERSAMPLE;
  ngen1 = n1 * oversample;
  ngen2 = n2 * oversample;
  data1 = (double*) fftw_malloc(sizeof(double) * n1 * n2);
  data2 = (double*) fftw_malloc(sizeof(double) * n1 * n2);
  work = (double*) fftw_malloc(sizeof(double) * ngen1 * ngen2);
  if (oversample > 1) {
    gen_1 = (double*) fftw_malloc(sizeof(double) * ngen1 * ngen2);
    gen_2 = (double*) fftw_malloc(sizeof(double) * ngen1 * ngen2);
  } else
    gen_1 = gen_2 = NULL;
  pixsize = PIXSIZE;
  pixsize_gen = PIXSIZE / static_cast<double>(oversample);

  fwhm1 = FWHM1;
  fwhm2 = FWHM2;
  sigi1 = SIGI1;
  sigi2 = SIGI2;
  is_full = false;
  esmooth1 = ESMOOTH1;
  esmooth2 = ESMOOTH2;

  is_binned = false;
  nbins = NBINS;
  bin_sparcity = 1;
  bincent01 = 0.0;
  bindelta1 = 0.0;
  bincent02 = 0.0;
  bindelta2 = 0.0;
  binval = NULL;

  sigi_final_computed = false;
  sigi_final_ntrials = 0;
  sigi_final1 = sigi_final2 = 0.0;

  //Set RNG seed
  unsigned long long int seed;
  seed = static_cast<unsigned long long int>(time(NULL));
  seed += static_cast<unsigned long long int>(clock());
  rangen.set_seed(seed);

  // Filtering
  isHipass = isMatched = std::make_pair(false, false);
  dblpair NaNpr =  std::make_pair(std::numeric_limits<double>::quiet_NaN(),
				  std::numeric_limits<double>::quiet_NaN());
  filtscale = qfactor = matched_fwhm = matched_sigi = matched_sigc = NaNpr;

  // Position generator if needed; note we generate in non-oversampled
  //  space then use oversampling to get positions
  use_clustered_pos = false;
  posgen = NULL;
  if (!powerspecfile.empty()) {
    use_clustered_pos = true;
    posgen = new positionGeneratorClustered(n1, n2, pixsize, powerspecfile);
  }
 
  //Set up array to hold 1D beams in each band (center normalized)
  //We have to decide how far out to go.  Note that the beams are set
  // up in oversampled space
  const double nfwhm = 4.5;
  ngauss1 = static_cast<unsigned int>(nfwhm * fwhm1 / pixsize_gen + 0.99999999);
  ngauss1 = 2 * ngauss1 + 1;
  gauss1 = new double[ngauss1];
  ngauss2 = static_cast<unsigned int>(nfwhm * fwhm2 / pixsize_gen + 0.99999999);
  ngauss2 = 2 * ngauss2 + 1;
  gauss2 = new double[ngauss2];
  doublebeam bm(fwhm1, fwhm2);
  bm.getBeamFac(1, ngauss1, pixsize_gen, gauss1);
  bm.getBeamFac(2, ngauss2, pixsize_gen, gauss2);

  //Set up additional smoothing 
  // This is done in un-oversampled space
  if (esmooth1 > 0 || esmooth2 > 0) {
    doublebeam ebm(esmooth1, esmooth2);
    if (esmooth1 > 0) {
      ngauss_add1 = static_cast<unsigned int>(nfwhm * esmooth1 / pixsize +
					      0.99999999);
      ngauss_add1 = 2 * ngauss_add1 + 1;
      gauss_add1 = new double[ngauss_add1];
      ebm.getBeamFac(1, ngauss_add1, pixsize, gauss_add1);
    } else {
      ngauss_add1 = 0;
      gauss_add1 = NULL;
    }
    if (esmooth2 > 0) {
      ngauss_add2 = static_cast<unsigned int>(nfwhm * esmooth2 / pixsize +
					      0.99999999);
      ngauss_add2 = 2 * ngauss_add2 + 1;
      gauss_add2 = new double[ngauss_add2];
      ebm.getBeamFac(2, ngauss_add2, pixsize, gauss_add2);
    } else {
      ngauss_add2 = 0;
      gauss_add2 = NULL;
    }
  } else {
    ngauss_add1 = ngauss_add2 = 0;
    gauss_add1 = gauss_add2 = NULL;
  }
}

simImageDouble::~simImageDouble() {
  fftw_free(data1);
  fftw_free(data2);
  fftw_free(work);
  if (gen_1 != NULL) fftw_free(gen_1);
  if (gen_2 != NULL) fftw_free(gen_2);
  if (posgen != NULL) delete posgen;
  delete[] gauss1;
  delete[] gauss2;
  if (ngauss_add1 > 0) delete[] gauss_add1;
  if (ngauss_add2 > 0) delete[] gauss_add2;
  if (binval != NULL) delete[] binval;
}

bool simImageDouble::isValid() const {
  if (n1 == 0) return false;
  if (n2 == 0) return false;
  if (oversample == 0) return false;
  if (pixsize <= 0.0) return false;
  if (fwhm1 <= 0.0) return false;
  if (fwhm2 <= 0.0) return false;
  if (sigi1 < 0.0) return false;
  if (sigi2 < 0.0) return false;
  if (esmooth1 < 0.0) return false;
  if (esmooth2 < 0.0) return false;
  return true;
}


/*!
  Downsample the internal image (inarr) into the output one (outarr)
  This rebinning preserves the mean flux per pixel
*/
void simImageDouble::downSample(unsigned int ni1, unsigned int ni2, 
				double* const inarr, unsigned int no1, 
				unsigned int no2, double* const outarr) {
  unsigned int osamp = ni1 / no1;
  if (osamp == 0)
    pofdExcept("simImageDouble", "downSample", "osamp is invalid (0)", 1);
  if (no2 * osamp != ni2)
    pofdExcept("simImageDouble", "downSample", "Dimension 1 doesn't match", 2);
  if (osamp == 1) {
    for (unsigned int i = 0; i < no1 * no2; ++i) outarr[i] = inarr[i];
    return;
  }

  double sum;
  unsigned int minxidx, maxxidx, minyidx, maxyidx;
  double *in_rowptr, *out_rowptr;
  double corrfac = 1.0 / static_cast<double>(osamp * osamp);
  for (unsigned int i = 0; i < no1; ++i) {
    minxidx = i * osamp;
    maxxidx = minxidx + osamp;
    out_rowptr =  outarr + i * no2;
    for (unsigned int j = 0; j < no2; ++j) {
      //This is the range in inarr we will sum over
      minyidx = j * osamp;
      maxyidx = minyidx + osamp;

      sum = 0.0;
      for (unsigned int k = minxidx; k < maxxidx; ++k) {
	in_rowptr = inarr + k * ni2;
	for (unsigned int h = minyidx; h < maxyidx; ++h)
	  sum += in_rowptr[h];
      }

      out_rowptr[j] = sum * corrfac;
    }
  }
}


void simImageDouble::convolveInner(unsigned int n, const double* const arr,
				   unsigned int ni1, unsigned int ni2,
				   double* const inarr,
				   double* const outarr) const {
  //inarr, outarr must be same size (ni1 by ni2)

  //Here we take major advantage of the fact that a Gaussian beam factorizes
  //to do this as two 1D convolutions, along rows then columns

  //This uses the work array as an internal working array

  //Do the column convolution first -- that is, convolve
  // along the second index of the input array, store in work
  //The first and last n items have to be handled specially 
  // -- or, at least, it's faster to do so
  //We do edge wrapping in the convolution to avoid apodizing
  // down the edges
  //There is one tricky index game here -- we store the results of
  // the column convolution in column-major rather than row-major order
  // to speed re-accessing them in the second convolution

  unsigned minidx, maxidx;
  double *rowptr, *subarrptr, currval;
  unsigned int nover2 = n/2;
  for (unsigned int i = 0; i < ni1; ++i) {
    rowptr = inarr + i * ni2;
    maxidx = n;
    for (unsigned int j = 0; j < nover2; ++j) {
      //Two steps -- wrapped part, unwrapped part
      minidx = nover2 - j; //Minimum index of Gauss array used in unwrapped part

      //wrapped part
      subarrptr = rowptr + j + ni2 - nover2;
      currval = subarrptr[0] * arr[0];
      for (unsigned int k = 1; k < minidx; ++k)
	currval += subarrptr[k] * arr[k];

      //unwrapped part
      subarrptr = rowptr + j - nover2;
      for (unsigned int k = minidx; k < maxidx; ++k)
	currval += subarrptr[k] * arr[k];
      work[j * ni1 + i] = currval;
    }
    for (unsigned int j = nover2; j < ni2 - nover2; ++j) {
      subarrptr = rowptr + j - nover2;
      currval = subarrptr[0] * arr[0];
      for (unsigned int k = 1; k < n; ++k)
	currval += subarrptr[k] * arr[k];
      work[j * ni1 + i] = currval;
    }
    minidx = 0;
    //Same wrapped/unwrapped part
    for (unsigned int j = ni2 - nover2; j < ni2; ++j) {
      maxidx = nover2 + (ni2 - j); //Limit on top end of array

      subarrptr = rowptr + j - nover2;
      currval = subarrptr[minidx] * arr[minidx];
      for (unsigned int k = minidx + 1; k < maxidx; ++k)
	currval += subarrptr[k] * arr[k];

      subarrptr = rowptr - maxidx;
      for (unsigned int k = maxidx; k < n; ++k)
	currval += subarrptr[k] * arr[k];

      work[j * ni1 + i] = currval;
    }
  }

  //Now convolve along the first index in the array, store that
  // into output array
  // Basically the same thing over again
  double *workptr;
  for (unsigned int j = 0; j < ni2; ++j) {
    maxidx = n;
    workptr = work + j * ni1;
    for (unsigned int i = 0; i < nover2; ++i) {
      minidx = nover2 - i;
      subarrptr = workptr + i + ni1 - nover2;
      currval = subarrptr[0] * arr[0];
      for (unsigned int k = 1; k < minidx; ++k)
	currval += subarrptr[k] * arr[k];
      subarrptr = workptr + i - nover2;
      for (unsigned int k = minidx; k < maxidx; ++k)
	currval += subarrptr[k] * arr[k];
      outarr[i * ni2 + j] = currval;
    }
    for (unsigned int i = nover2; i < ni1 - nover2; ++i) {
      subarrptr = workptr + i - nover2;
      currval = subarrptr[0] * arr[0];
      for (unsigned int k = 1; k < n; ++k)
	currval += subarrptr[k] * arr[k];
      outarr[i * ni2 + j] = currval;
    }
    minidx = 0;
    for (unsigned int i = ni1 - nover2; i < ni1; ++i) {
      maxidx = nover2 + (ni1 - i); //Limit on top end of array
      subarrptr = workptr + i - nover2;
      currval = subarrptr[minidx] * arr[minidx];
      for (unsigned int k = minidx + 1; k < maxidx; ++k)
	currval += subarrptr[k] * arr[k];
      subarrptr = workptr - maxidx;
      for (unsigned int k = maxidx; k < n; ++k)
	currval += subarrptr[k] * arr[k];
      outarr[i * ni2 + j] = currval;
    }
  }
}

void simImageDouble::convolveWithBeam() {
  if (!is_full)
    throw pofdExcept("simImageDouble", "convolveWithBeam",
		     "Trying to convolve empty image", 1);
  if (oversample > 1) {
    //First, convolve generated image in band 1, using work as a temporary
    convolveInner(ngauss1, gauss1, ngen1, ngen2, gen_1, gen_1);
    //Now downsample to data1
    downSample(ngen1, ngen2, gen_1, n1, n2, data1);
    //Same for band 2
    convolveInner(ngauss2, gauss2, ngen1, ngen2, gen_2, gen_2);
    downSample(ngen1, ngen2, gen_2, n1, n2, data2);
  } else {
    //Simpler -- image was generated in data1/data2, so just convolve to there
    convolveInner(ngauss1, gauss1, n1, n2, data1, data1);
    convolveInner(ngauss2, gauss2, n1, n2, data2, data2);
  }
}

void simImageDouble::convolveWithAdd() {
  //This can only be done to data1, data2
  if (! is_full )
    throw pofdExcept("simImageDouble", "convolveWithAdd",
		     "Trying to convolve empty image", 1);

  if (esmooth1 > 0) {
    //From data to data, no oversampling
    convolveInner(ngauss_add1, gauss_add1, n1, n2, data1, data1);

    //Normalization step; the idea is to keep the peak flux values
    // the same for sources.  This depends both on the beam size
    // and the extra amount of smoothing, and is derived from the
    // relation between the peak value of a Gaussian beam and it's
    // area
    const double prefac = 4*std::log(2)/pofd_coverage::pi;
    double normval = prefac * (fwhm1*fwhm1 + esmooth1*esmooth1) * 
      pixsize*pixsize / (fwhm1*fwhm1*esmooth1*esmooth1);
    for (unsigned int i = 0; i < n1*n2; ++i)
      data1[i] *= normval;
  }

  if (esmooth2 > 0) {
    //Same, band 2
    convolveInner(ngauss_add2, gauss_add2, n1, n2, data2, data2);

    const double prefac = 4*std::log(2)/pofd_coverage::pi;
    double normval = prefac * (fwhm2*fwhm2 + esmooth2*esmooth2) * 
      pixsize*pixsize / (fwhm2*fwhm2*esmooth2*esmooth2);
    for (unsigned int i = 0; i < n1*n2; ++i)
      data2[i] *= normval;
  }

}


double 
simImageDouble::getFinalNoiseHelper(unsigned int ntrials,
				    double* const data, double sigi,
				    double fwhm, double esmooth, 
				    unsigned int ngauss_add,
				    const double* const gauss_add,
				    const fourierFilter* const filt) const {
  // Compute esmooth prefactor if needed
  const double prefac = 4 * std::log(2) / pofd_coverage::pi;
  double normval = 1.0;
  if (esmooth > 0) {
    double fwhmsq = fwhm * fwhm;
    double esmoothsq = esmooth * esmooth;
    normval = prefac * (fwhmsq + esmoothsq) * pixsize * pixsize / 
      (fwhmsq + esmoothsq);
  }

  double var, var1, var2, mn, dtmp;
  var = 0.0;
  for (unsigned int idx = 0; idx < ntrials; ++idx) {
    // Fill with random noise
    for (unsigned int i = 0; i < n1 * n2; ++i)
      data[i] = sigi * rangen.gauss();

    // Apply additional Gaussian smoothing if needed
    if (esmooth > 0) {
      convolveInner(ngauss_add, gauss_add, n1, n2, data, data);
      // Then renormalize
      for (unsigned int i = 0; i < n1 * n2; ++i)
	data[i] *= normval;
    }

    // Now filtering
    filt->filter(n1, n2, pixsize, data);

    // Measure using two pass algorithm
    mn = data[0];
    for (unsigned int i = 1; i < n1 * n2; ++i)
      mn += data[i];
    mn /= static_cast<double>(n1 * n2);
    var1 = var2 = 0.0;
    for (unsigned int i = 0; i < n1 * n2; ++i) {
      dtmp = data[i] - mn;
      var1 += dtmp * dtmp;
      var2 += dtmp;
    }
    var1 -= var2 * var2 / static_cast<double>(n1 * n2);
    var += var1 / static_cast<double>(n1 * n2 - 1);
  }
  
  if (ntrials > 1)
    var /= static_cast<double>(ntrials);
  return sqrt(var);
}

/*
  \param[in] ntrials Number of trials to do if measuring filtered noise.
  \param[in] filt1 Fourier space filter for band1, maybe band 2
  \param[in] filt2 Fourier space filter for band2

  \returns Instrument noise including effects of filtering, etc.
  
  This can include both filtering and additional Gaussian smoothing.
  Note that it is not necessary for the data array to be filled
  for this function to work, but it can take additional storage
  equal to the size of the image.  

  This is not particularly cheap the first time it is called.
  After that it uses the previously computed value if possible.

  If filt1 is set but filt2 is not set, filt1 is used for both
*/
std::pair<double, double> 
simImageDouble::getFinalNoise(unsigned int ntrials,
			      const fourierFilter* const filt1,
			      const fourierFilter* const filt2) const {

  // esmooth normalization factor input, when analyticity is possible
  const double noise_prefac = sqrt(2 * std::log(2) / pofd_coverage::pi);

  bool recompute = true;
  if (sigi_final_computed) {
    // Decide if we can re-use the previous value
    if ((filt1 == NULL) && (filt2 == NULL)) recompute = false; // Ignore ntrials -- not needed
    else if (ntrials <= sigi_final_ntrials) recompute = false;
  }
  if (!recompute)
    return std::make_pair(sigi_final1, sigi_final2);

  // Figure out if we need temporary storage for sims
  double *tmpdata;
  if ((filt1 != NULL) && (filt2 != NULL)) {
    if (ntrials == 0)
      throw pofdExcept("simImageDouble", "getFinalNoise",
		       "Invalid (non-positive) ntrials", 1);
    tmpdata = (double*) fftw_malloc(sizeof(double) * n1 * n2);
  } else tmpdata = NULL;

  const fourierFilter *f1, *f2;
  if (filt1 != NULL) {
    f1 = filt1;
    if (filt2 == NULL) f2 = filt2; else f2 = filt1;
  } else if (filt2 != NULL) {
    f1 = NULL;
    f2 = filt2;
  } else f1 = f2 = NULL;

  // First, sigi1
  if (sigi1 == 0)
    sigi_final1 = 0.0;
  else if (f1 == NULL) {
    if (esmooth1 <= 0.0) 
      sigi_final1 = sigi1;
    else {
      // Only extra Gaussian smoothing.  Can be done analytically
      double fwhm1sq = fwhm1 * fwhm1;
      sigi_final1 = noise_prefac * sigi1 * pixsize * 
	(fwhm1sq + esmooth1 * esmooth1) / (esmooth1 * fwhm1sq);
    }
  } else {
    // Filtering, and maybe smoothing.  
    // We have to do this by making a simulated image, filling it with 
    // noise, etc.  This is why this function may not be cheap,
    // and why we support multiple trials
    sigi_final1 = getFinalNoiseHelper(ntrials, tmpdata, sigi1, fwhm1, esmooth1, 
				      ngauss_add1, gauss_add1, f1);
  }
  if (sigi2 == 0)
    sigi_final2 = 0.0;
  else if (f2 == NULL) {
    if (esmooth2 <= 0.0) 
      sigi_final2 = sigi2;
    else {
      double fwhm2sq = fwhm2 * fwhm2;
      sigi_final2 = noise_prefac * sigi2 * pixsize * 
	(fwhm2sq + esmooth2 * esmooth2) / (esmooth2 * fwhm2sq);
    }
  } else {
    sigi_final2 = getFinalNoiseHelper(ntrials, tmpdata, sigi2, fwhm2, esmooth2, 
				      ngauss_add2, gauss_add2, f2);
  }
  if (tmpdata != NULL) fftw_free(tmpdata);

  sigi_final_computed = true;
  sigi_final_ntrials = ntrials;
  return std::make_pair(sigi_final1, sigi_final2);
}

double simImageDouble::getArea() const {
  double val = pixsize / 3600.0;
  return val * val * n1 * n2;
}

/*!
  Generates a simulated image
  \params[in] model Base number counts model
  \params[in] n0 Number of sources per area to generate
  \params[in] meansub Do mean subtraction
  \param[in] filt1 Fourier space filter to apply to band 1, maybe band 2.  
             If NULL, no filtering.
  \param[in] filt2 Fourier space filter to apply to band 2	     
  \params[in] bin Create binned image data
  \params[in] sparsebin Only take every this many pixels in binned image.
                         Does nothing if no binning.

  If filt1 is set but filt2 is not, then filt1 is applied to both.
  There is no way to filter only the first one.
*/
void simImageDouble::realize(const numberCountsDouble& model,
			     double n0, bool meansub, 
			     const fourierFilter* const filt1, 
			     const fourierFilter* const filt2,
			     bool bin, unsigned int sparsebin) {

  if (!isValid())
    throw pofdExcept("simImageDouble", "realize",
		     "Trying to realize image with invalid parameters", 1);
  if (!model.isValid())
    throw pofdExcept("simImageDouble", "realize",
		     "Trying to realize model with invalid parameters", 2);
  if (n0 <= 0.0)
    throw pofdExcept("simImageDouble", "realize", "n0 must be positive", 3);

  double area = getArea();
  unsigned int nsrcs = static_cast<unsigned int>(area * n0);

  // Set up the position generator if it will be needed
  if (use_clustered_pos) posgen->generate(rangen);

  //Inject sources
  std::pair<double, double> src;
  if (oversample > 1) {
    //Generate in oversampled gen_1, gen_2
    std::memset(gen_1, 0, ngen1 * ngen2 * sizeof(double));
    std::memset(gen_2, 0, ngen1 * ngen2 * sizeof(double));
    if (nsrcs > 0) {
      unsigned int idx1, idx2, combidx;
      if (use_clustered_pos) {
	// Clustered positions
	std::pair<unsigned int, unsigned int> pos;
	for (unsigned int i = 0; i < nsrcs; ++i) {
	  pos = posgen->getPosition(rangen, oversample);
	  combidx = pos.first * ngen2 + pos.second;
	  src = model.genSource(rangen.doub(), rangen.gauss());
	  gen_1[combidx] += src.first;
	  gen_2[combidx] += src.second;
	} 
      } else {
	for (unsigned int i = 0; i < nsrcs; ++i) {
	  idx1 = static_cast<unsigned int>(rangen.doub() * ngen1);
	  idx2 = static_cast<unsigned int>(rangen.doub() * ngen2);
	  combidx = idx1 * ngen2 + idx2;
	  src = model.genSource(rangen.doub(), rangen.gauss());
	  gen_1[combidx] += src.first;
	  gen_2[combidx] += src.second;
	}
      }
    }
  } else {
    //Generate in data1, data2
    std::memset(data1, 0, n1 * n2 * sizeof(double));
    std::memset(data2, 0, n1 * n2 * sizeof(double));

    //Inject sources
    if (nsrcs > 0) {
      unsigned int idx1, idx2, combidx;
      if (use_clustered_pos) {
	// Clustered positions
	std::pair<unsigned int, unsigned int> pos;
	for (unsigned int i = 0; i < nsrcs; ++i) {
	  pos = posgen->getPosition(rangen);
	  combidx = pos.first * n2 + pos.second;
	  src = model.genSource(rangen.doub(), rangen.gauss());
	  data1[combidx] += src.first;
	  data2[combidx] += src.second;
	} 
      } else {
	for (unsigned int i = 0; i < nsrcs; ++i) {
	  idx1 = static_cast<unsigned int>(rangen.doub() * n1);
	  idx2 = static_cast<unsigned int>(rangen.doub() * n2);
	  combidx = idx1 * n2 + idx2;
	  src = model.genSource(rangen.doub(), rangen.gauss());
	  data1[combidx] += src.first;
	  data2[combidx] += src.second;
	}
      }
    }
  }
  is_full = true;
  convolveWithBeam();

  //Add instrument noise
  if (sigi1 > 0.0)
    for (unsigned int i = 0; i < n1*n2; ++i)
      data1[i] += sigi1 * rangen.gauss();
  if (sigi2 > 0.0)
    for (unsigned int i = 0; i < n1*n2; ++i)
      data2[i] += sigi2 * rangen.gauss();

  //Extra smoothing, if set
  if ((esmooth1 > 0) || (esmooth2 > 0))
    convolveWithAdd();

  // Apply filtering.  Note this is done after adding noise and
  // downsampling to the final resolution (if oversampling is used).
  // Filtering will always result in mean subtraction, so no need to
  // do that twice.
  isHipass = isMatched = std::make_pair(false, false);
  dblpair NaNpr = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
				 std::numeric_limits<double>::quiet_NaN());
  filtscale = qfactor = matched_fwhm = matched_sigi = matched_sigc = NaNpr;
  const fourierFilter *f1, *f2;
  if (filt1 != NULL) {
    f1 = filt1;
    if (filt2 == NULL) f2 = f1; else f2 = filt2;
  } else if (filt2 != NULL) {
    f1 = NULL;
    f2 = filt2;
  } else f1 = f2 = NULL;
  // Deal with band 1
  if (f1 != NULL) {
    f1->filter(n1, n2, pixsize, data1);
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
  // And band 2
  if (f2 != NULL) {
    f2->filter(n1, n2, pixsize, data2);
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
  if (meansub && ((f1 == NULL) || (f2 == NULL))) meanSubtract();

  //bin
  is_binned = false;
  if (bin) applyBinning(sparsebin);

}

std::pair<double,double> simImageDouble::meanSubtract() {
  std::pair<double,double> mn;
  mn = getMean();
  double mn1 = mn.first;
  for (unsigned int i = 0; i < n1*n2; ++i)
    data1[i] -= mn1;
  double mn2 = mn.second;
  for (unsigned int i = 0; i < n1*n2; ++i)
    data2[i] -= mn2;
  if (is_binned) {
    bincent01 -= mn1;
    bincent02 -= mn2;
  }
  return mn;
}

void simImageDouble::getMinMax(double& min1, double& max1, double& min2,
			       double& max2) const {
  if (! is_full ) {
    min1 = std::numeric_limits<double>::quiet_NaN();
    max1 = std::numeric_limits<double>::quiet_NaN();
    min2 = std::numeric_limits<double>::quiet_NaN();
    max2 = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  min1 = data1[0];
  max1 = data1[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    if (data1[i] < min1) min1 = data1[i];
  for (unsigned int i = 1; i < n1*n2; ++i)
    if (data1[i] > max1) max1 = data1[i];
  min2 = data2[0];
  max2 = data2[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    if (data2[i] < min2) min2 = data2[i];
  for (unsigned int i = 1; i < n1*n2; ++i)
    if (data2[i] > max2) max2 = data2[i];
}

std::pair<double, double> 
simImageDouble::getMean() const {
  if (!is_full)
    throw pofdExcept("simImageDouble", "getMean",
		     "Trying to get means of empty images", 1);
  double norm = 1.0 / static_cast<double>(n1 * n2);
  double mn1 = data1[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    mn1 += data1[i];
  mn1 *= norm;
  double mn2 = data2[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    mn2 += data2[i];
  mn2 *= norm;
  return std::make_pair(mn1, mn2);
}

void simImageDouble::getMeanAndVar(double& mn1, double& var1,
				   double& mn2, double& var2) const {
  if (!is_full)
    throw pofdExcept("simImageDouble", "getMeanAndVar",
		     "Trying to get mean and vars of empty images", 1);
  //Use corrected two pass algorithm
  double norm = 1.0/static_cast<double>(n1 * n2);
  mn1 = data1[0];
  for (unsigned int i = 1; i < n1 * n2; ++i)
    mn1 += data1[i];
  mn1 *= norm;

  double sum1, sum2, tmp;
  tmp = data1[0]-mn1;
  sum1 = tmp * tmp;
  sum2 = tmp;
  for (unsigned int i = 1; i < n1 * n2; ++i) {
    tmp = data1[i] - mn1;
    sum1 += tmp * tmp;
    sum2 += tmp;
  }
  var1 = (sum1 - norm * sum2 * sum2)/static_cast<double>(n1 * n2 - 1);

  mn2 = data2[0];
  for (unsigned int i = 1; i < n1 * n2; ++i)
    mn2 += data2[i];
  mn2 *= norm;
  tmp = data2[0] - mn2;
  sum1 = tmp * tmp;
  sum2 = tmp;
  for (unsigned int i = 1; i < n1*n2; ++i) {
    tmp = data2[i] - mn2;
    sum1 += tmp * tmp;
    sum2 += tmp;
  }
  var2 = (sum1 - norm * sum2 * sum2) / static_cast<double>(n1 * n2 - 1);
}

std::pair<double,double> simImageDouble::getBeamSum() const {
  if (ngauss1 == 0)
    throw pofdExcept("simImageDouble", "getBeamSum",
		     "No beam pixels, band 1", 1);
  if (ngauss2 == 0)
    throw pofdExcept("simImageDouble", "getBeamSum",
		     "No beam pixels, band 2", 2);
  
  double sum1D_1 = gauss1[0];
  for (unsigned int i = 1; i < ngauss1; ++i)
    sum1D_1 += gauss1[i];

  double sum1D_2 = gauss2[0];
  for (unsigned int i = 1; i < ngauss2; ++i)
    sum1D_2 += gauss2[i];

  return std::make_pair(sum1D_1 * sum1D_1, sum1D_2 * sum1D_2);
}


std::pair<double,double> simImageDouble::getBeamSumSq() const {
  if (ngauss1 == 0)
    throw pofdExcept("simImageDouble", "getBeamSumSq",
		     "No beam pixels, band 1", 1);
  if (ngauss2 == 0)
    throw pofdExcept("simImageDouble", "getBeamSumSq",
		     "No beam pixels, band 2", 2);
  
  double tmp = gauss1[0];
  double sum1D_1 = tmp * tmp;
  for (unsigned int i = 1; i < ngauss1; ++i) {
    tmp = gauss1[i];
    sum1D_1 += tmp * tmp;
  }

  tmp = gauss2[0];
  double sum1D_2 = tmp * tmp;
  for (unsigned int i = 1; i < ngauss2; ++i) {
    tmp = gauss2[i];
    sum1D_2 += tmp * tmp;
  }

  return std::make_pair(sum1D_1 * sum1D_1, sum1D_2 * sum1D_2);
}

/*!
  \param[in] sparsebin Only take every this many pixels when binning. 1 means
                        fully sampled (0 is also interpreted to mean the same
			thing)

  Keeps the original, unbinned image around as well
 */
void simImageDouble::applyBinning(unsigned int sparsebin) {
  if (!is_full)
    throw pofdExcept("simImageDouble", "applyBinning",
		     "Trying to bin empty image", 1);

  if (nbins == 0) throw pofdExcept("simImageDouble", "applyBinning",
				   "Trying to bin with no bins", 2);
  if (sparsebin >= n1 * n2) 
    throw pofdExcept("simImageDouble", "applyBinning",
		     "Sparse binning factor larger than simulated image", 3);

  //Only allocated the first time we bin
  //There is no way to change nbins after initialization
  if (binval == NULL)
    binval = new unsigned int[nbins * nbins];
  std::memset(binval, 0, nbins * nbins * sizeof(unsigned int));

  //First, we need the minimum and maximum
  double minval1, maxval1, minval2, maxval2;
  getMinMax(minval1, maxval1, minval2, maxval2);
  
  //We want to put the max and min in the center of the top and
  // bottom bin
  bincent01 = minval1;
  if (nbins == 1)
    bindelta1 = 2*(maxval1-minval1);
  else
    bindelta1 = (maxval1-minval1)/static_cast<double>(nbins-1);
  bincent02 = minval2;
  if (nbins == 1)
    bindelta2 = 2*(maxval2-minval2);
  else
    bindelta2 = (maxval2-minval2)/static_cast<double>(nbins-1);

  //And... bin
  double ibindelta1 = 1.0/bindelta1;
  double ibindelta2 = 1.0/bindelta2;
  unsigned int idx1, idx2;
  if (sparsebin > 1) bin_sparcity = sparsebin; else bin_sparcity = 1;
  for (unsigned int i = 0; i < n1 * n2; i += bin_sparcity) {
    idx1 = static_cast<unsigned int>((data1[i] - bincent01) * 
				     ibindelta1 + 0.5);
    idx2 = static_cast<unsigned int>((data2[i] - bincent02) * 
				     ibindelta2 + 0.5);
    binval[idx1 * nbins + idx2] += 1;
  }
  is_binned = true;
}

/*!
  Writes a single band to an already open fits file, creating the extension

  \param[in] outfile File to write to
  \param[in] idx Which band to write (1 or 2)
  \returns 0 on success, an error code (!=0) for anything else
*/
int simImageDouble::writeFits(fitsfile* fp, unsigned int idx) const {

  if (idx < 1 || idx > 2)
    throw pofdExcept("simImageDouble", "writeFits",
		     "Invalid index", 1);

  int status = 0;

  //Make image array
  long axissize[2];
  axissize[0] = static_cast<long>(n1);
  axissize[1] = static_cast<long>(n2);
  fits_create_img(fp, DOUBLE_IMG, 2, axissize, &status);
  
  // Write header
  int itmp = 0;
  double tmpval;

  // Band
  fits_write_key(fp, TUINT, const_cast<char*>("BAND"),
		 &idx, const_cast<char*>("Which band is this image?"),
		 &status);

  // Model params
  fits_write_key(fp, TSTRING, const_cast<char*>("MODEL"),
		 const_cast<char*>("Spline-Log Normal"), 
		 const_cast<char*>("Model type"),
		 &status);

  //Sim params
  tmpval = fwhm1;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM1"), &tmpval, 
		 const_cast<char*>("Beam fwhm, band 1 [arcsec]"), 
		 &status);
  tmpval = fwhm2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM2"), &tmpval, 
		 const_cast<char*>("Beam fwhm, band 2 [arcsec]"), 
		 &status);
  if ((esmooth1 > 0) && (esmooth2 > 0)) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &itmp,
		   const_cast<char*>("Additional smoothing applied"), 
		   &status);
    tmpval = esmooth1;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("ESMOOTH1"), &tmpval, 
		   const_cast<char*>("Extra smoothing fwhm, band 1 [arcsec]"), 
		   &status);
    tmpval = esmooth2;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("ESMOOTH2"), &tmpval, 
		   const_cast<char*>("Extra smoothing fwhm, band 2 [arcsec]"), 
		   &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &itmp,
		   const_cast<char*>("Additional smoothing applied"), 
		   &status);
  }

  tmpval = sigi1;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGI_1"), &tmpval, 
		 const_cast<char*>("Raw instrument noise, band 1"), 
		 &status);
  tmpval = sigi2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGI_2"), &tmpval, 
		 const_cast<char*>("Raw instrument noise, band 2"), 
		 &status);
  if (sigi_final_computed) {
    tmpval = sigi_final1;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGIFNL1"), &tmpval, 
		   const_cast<char*>("Final instrument noise, band 1"), 
		   &status);
    tmpval = sigi_final2;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGIFNL2"), &tmpval, 
		   const_cast<char*>("Final instrument noise, band 2"), 
		   &status);
  }
  
  if (oversample > 1) {
    unsigned int utmp = oversample;
    fits_write_key(fp, TUINT, const_cast<char*>("OVERSMPL"), &utmp, 
		    const_cast<char*>("Oversampling factor"), 
		    &status);
  }

  if (use_clustered_pos) itmp = 1; else itmp = 0;
  fits_write_key(fp, TLOGICAL, const_cast<char*>("CLUSTPOS"), &itmp,
		 const_cast<char*>("Use clustered positions"), &status);

  if (isHipass.first) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT1"), &itmp,
		   const_cast<char*>("Is band1 hipass filtered?"), 
		   &status);
    tmpval = filtscale.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTSCL1"), &tmpval,
		   const_cast<char*>("Filtering scale1 [arcsec]"), &status);
    tmpval = qfactor.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTQ1"), &tmpval,
		   const_cast<char*>("Filtering apodization1"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT1"), &itmp,
		   const_cast<char*>("Is band1 hipass filtered?"), &status);
  }
  if (isHipass.second) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT2"), &itmp,
		   const_cast<char*>("Is band2 hipass filtered?"), 
		   &status);
    tmpval = filtscale.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTSCL2"), &tmpval,
		   const_cast<char*>("Filtering scale2 [arcsec]"), &status);
    tmpval = qfactor.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTQ2"), &tmpval,
		   const_cast<char*>("Filtering apodization2"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT2"), &itmp,
		   const_cast<char*>("Is band2 hipass filtered?"), &status);
  }
  if (isMatched.first) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT1"), &itmp,
		   const_cast<char*>("Is band1 match filtered?"), &status);
    tmpval = matched_fwhm.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITFWHM1"), &tmpval,
		   const_cast<char*>("Matched filtering FWHM1 [arcsec]"), 
		   &status);
    tmpval = matched_sigi.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGI1"), &tmpval,
		   const_cast<char*>("Matched filtering sigi1"), &status);
    tmpval = matched_sigc.first;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGC1"), &tmpval,
		   const_cast<char*>("Matched filtering sigc1"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT1"), &itmp,
		   const_cast<char*>("Is band1 match filtered?"), &status);
  }
  if (isMatched.second) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT2"), &itmp,
		   const_cast<char*>("Is band1 match filtered?"), &status);
    tmpval = matched_fwhm.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITFWHM2"), &tmpval,
		   const_cast<char*>("Matched filtering FWHM2 [arcsec]"), 
		   &status);
    tmpval = matched_sigi.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGI2"), &tmpval,
		   const_cast<char*>("Matched filtering sigi2"), &status);
    tmpval = matched_sigc.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FITSIGC2"), &tmpval,
		   const_cast<char*>("Matched filtering sigc2"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT2"), &itmp,
		   const_cast<char*>("Is band2 match filtered?"), &status);
  }

  fits_write_key(fp, TSTRING, const_cast<char*>("VERSION"),
		 const_cast<char*>(pofd_coverage::version), 
		 const_cast<char*>("pofd_coverage version"),
		 &status);
  fits_write_history(fp, 
		     const_cast<char*>("Simulated image from pofd_coverage"),
		     &status);
  fits_write_date(fp, &status);

  // Astrometry
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE1"),
		 const_cast<char*>("RA---TAN"),
		 const_cast<char*>("WCS: Projection type axis 1"),&status);
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE2"),
		 const_cast<char*>("DEC--TAN"),
		 const_cast<char*>("WCS: Projection type axis 2"),&status);
  tmpval = n1/2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRPIX1"), &tmpval, 
		 const_cast<char*>("Ref pix of axis 1"), &status);
  tmpval = n2/2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRPIX2"), &tmpval, 
		 const_cast<char*>("Ref pix of axis 2"), &status);
  tmpval = 90.0; //Arbitrary
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRVAL1"), &tmpval, 
		 const_cast<char*>("val at ref pix axis 1"), &status);
  tmpval = 10.0; //Arbitrary
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRVAL2"), &tmpval, 
		 const_cast<char*>("val at ref pix axis 2"), &status);
  tmpval = - pixsize/3600.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD1_1"), &tmpval, 
		 const_cast<char*>("Pixel scale axis 1,1"), &status);
  tmpval = 0.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD1_2"), &tmpval, 
		 const_cast<char*>("Pixel scale axis 1,2"), &status);
  tmpval = 0.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD2_1"), &tmpval, 
		 const_cast<char*>("Pixel scale axis 2,1"), &status);
  tmpval = pixsize/3600.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD2_2"), &tmpval, 
		 const_cast<char*>("Pixel scale axis 2,2"), &status);
  tmpval = 2000.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("EPOCH"), &tmpval, 
		 const_cast<char*>("WCS: Epoch of celestial pointing"), 
		 &status);
  fits_write_key(fp, TDOUBLE, const_cast<char*>("EQUINOX"), &tmpval, 
		 const_cast<char*>("WCS: Equinox of celestial pointing"), 
		 &status);

  // History, date
  fits_write_history(fp, 
		     const_cast<char*>("Simulated image from pofd_coverage"),
		     &status);
  fits_write_date(fp, &status);

  //Do data writing.  We have to make a transposed copy of the
  // data to do this, which is irritating as hell
  double *tmpdata = new double[n1];
  long fpixel[2] = { 1, 1 };
  for ( unsigned int j = 0; j < n2; ++j ) {
    if (idx == 1)
      for (unsigned int i = 0; i < n1; ++i) tmpdata[i] = data1[i * n2 + j];
    else
      for (unsigned int i = 0; i < n1; ++i) tmpdata[i] = data2[i * n2 + j];
    fpixel[1] = static_cast<long>(j+1);
    fits_write_pix(fp, TDOUBLE, fpixel, n1, tmpdata, &status);
  }
  delete[] tmpdata;
  
  return status;
}

/*!
  \param[in] outputfile1 File to write to.  
  \returns 0 on success, an error code (!=0) for anything else

  Written as two image extensions.
*/
int simImageDouble::writeToFits(const std::string& outputfile) const {

  if (!is_full)
    throw pofdExcept("simImageDouble", "writeToFits",
		     "Trying to write image without realizing", 1);

  int status = 0;
  fitsfile *fp;
  fits_create_file(&fp, outputfile.c_str(), &status);
  if (status) {
    fits_report_error(stderr, status);
    return status;
  }

  status = writeFits(fp, 1);
  if (status)
    throw pofdExcept("simImageDouble", "writeToFits",
		     "Failure writing band 1 map", 2);

  status = writeFits(fp, 2);
  if (status) {
    throw pofdExcept("simImageDouble", "writeToFits",
		     "Failure writing band 2 map", 3);
    return status;
  }

  // Close up
  fits_close_file(fp, &status);
  if (status)
    fits_report_error(stderr, status);
  return status;
}

/*!
  \param[in] outfile File to write to
*/
int simImageDouble::writeProbImageToFits(const std::string& outfile) const {
  if (!use_clustered_pos)
    throw pofdExcept("simImageDouble", "writeProbImageToFits",
		     "No probability image to write", 1);
  return posgen->writeProbToFits(outfile);
}

/*!
  \param[in] obj_id Open HDF5 handle to write to.

  Writes information about the position generator to the open handle,
  creating a group called "PositionGenerator"
*/
void simImageDouble::writePositionGeneratorToHDF5Handle(hid_t obj_id) const {

  if (H5Iget_ref(obj_id) < 0)
    throw pofdExcept("simImageDouble", "writePositionGeneratorToHDF5Handle",
		     "Given non-open obj_id to write to", 1);

  hid_t group_id;
  group_id = H5Gcreate(obj_id, "PositionGenerator", H5P_DEFAULT, H5P_DEFAULT, 
		      H5P_DEFAULT);
  if (H5Iget_ref(group_id) < 0)
    throw pofdExcept("simImageDouble", "writePositionGeneratorToHDF5",
		     "Failed to create HDF5 model group", 2);

  hsize_t adims;
  hid_t mems_id, att_id;
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  if (use_clustered_pos) {
    const char postype[] = "Clustered";
    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, strlen(postype)); 
    att_id = H5Acreate1(group_id, "PositionType", datatype,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, datatype, postype);
    H5Aclose(att_id);
    posgen->writeToHDF5Handle(group_id);
  } else {
    const char postype[] = "Uniform";
    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, strlen(postype)); 
    att_id = H5Acreate1(group_id, "PositionType", datatype,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, datatype, postype);
    H5Aclose(att_id);
  }
  H5Sclose(mems_id);
  H5Gclose(group_id); 
}
