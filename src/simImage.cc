#include<ctime>
#include<cmath>
#include<limits>
#include<cstring>
#include<fstream>

#include<fitsio.h>
#include<fftw3.h> // For fftw_malloc

#include "../include/global_settings.h"
#include "../include/simImage.h"
#include "../include/beam.h"
#include "../include/pofdExcept.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

/*!
  \param[in] N1 Dimension of simulated image along 1st axis
  \param[in] N2 Dimension of simulated image along 2nd axis
  \param[in] PIXSIZE Pixel size in arcseconds.
  \param[in] FWHM The FWHM of the simulated image in arcsec, before
              filtering or smoothing
  \param[in] SIGI Gaussian instrumental noise in the units the model is in.
              This is the noise before additional smoothing or filtering.
  \param[in] ESMOOTH Additional Gaussian smoothing FWHM.  If non-zero,
              additional smoothing is applied.
  \param[in] OVERSAMPLE Oversampling of simulated image.  Must be odd.
              1 means no additional oversampling.
  \param[in] NBINS Number of bins, if binning is applied.
  \param[in] powerspecfile File containing power spectrum.  If set, the
              source positions are generated using this P(k).  If not set, they
              are uniformly distributed across the image.
*/
simImage::simImage(unsigned int N1, unsigned int N2, double PIXSIZE,
                   double FWHM, double SIGI, double ESMOOTH, 
                   unsigned int OVERSAMPLE, unsigned int NBINS, 
                   const std::string& powerspecfile) {

  if (N1 == 0)
    throw pofdExcept("simImage", "simImage", "Invalid (non-positive) N1", 1);
  if (N2 == 0)
    throw pofdExcept("simImage", "simImage", "Invalid (non-positive) N2", 2);
  if (PIXSIZE <= 0.0)
    throw pofdExcept("simImage", "simImage", 
                     "Invalid (non-positive) PIXSIZE", 3);
  if (FWHM <= 0.0)
    throw pofdExcept("simImage", "simImage", "Invalid (non-positive) FWHM", 4);
  if (OVERSAMPLE % 2 == 0)
    throw pofdExcept("simImage", "simImage", "Invalid (even) OVERSAMPLE", 5);
  
  n1 = N1;
  n2 = N2;
  oversample = OVERSAMPLE;
  ngen1 = n1 * oversample;
  ngen2 = n2 * oversample;
  data = (double*) fftw_malloc(sizeof(double) * n1 * n2);
  work = (double*) fftw_malloc(sizeof(double) * ngen1 * ngen2);
  if (oversample > 1)
    gen_image = (double*) fftw_malloc(sizeof(double) * ngen1 * ngen2);
  else
    gen_image = nullptr;
  pixsize = PIXSIZE;
  pixsize_gen = PIXSIZE / static_cast<double>(oversample);

  fwhm = FWHM;
  sigi = SIGI;
  is_full = false;
  esmooth = ESMOOTH;
  is_binned = false;
  nbins = NBINS;
  bin_sparcity = 1;
  bincent0 = 0.0;
  bindelta = 0.0;
  binval = nullptr;

  sigi_final_computed = false;
  sigi_final_ntrials = 0;
  sigi_final = 0.0;

  //Set RNG seed
  unsigned long long int seed;
  seed = static_cast<unsigned long long int>(time(nullptr));
  seed += static_cast<unsigned long long int>(clock());
  rangen.set_seed(seed);
  
  // Filtering variables
  isHipass = false;
  filtscale = qfactor = std::numeric_limits<double>::quiet_NaN();
  isMatched = false;
  matched_fwhm = std::numeric_limits<double>::quiet_NaN();
  matched_sigi = matched_sigc = std::numeric_limits<double>::quiet_NaN();

  // Position generator if needed; note we generate in non-oversampled
  //  space but use interpolation to get to interpolated space
  use_clustered_pos = false;
  posgen = nullptr;
  if (!powerspecfile.empty()) {
    use_clustered_pos = true;
    posgen = new positionGeneratorClustered(n1, n2, pixsize, powerspecfile);
  }

  // Set up array to hold 1D beams (center normalized)
  // Note that the beam is set up in oversampled space
  const double nfwhm = 4.5;
  ngauss = static_cast<unsigned int>(nfwhm * fwhm / pixsize_gen + 0.99999999);
  ngauss = 2 * ngauss + 1;
  gauss = new double[ngauss];
  beam bm(fwhm);
  bm.getBeamFac(ngauss, pixsize_gen, gauss);

  //Set up additional smoothing 
  // This is done in un-oversampled space
  if (esmooth > 0) {
    ngauss_add = static_cast<unsigned int>(nfwhm * esmooth / pixsize +
                                           0.99999999);
    ngauss_add = 2 * ngauss_add + 1;
    gauss_add = new double[ngauss_add];
    beam ebm(esmooth);
    ebm.getBeamFac(ngauss_add, pixsize, gauss_add);
  } else {
    ngauss_add = 0;
    gauss_add = nullptr;
  }
}

simImage::~simImage() {
  fftw_free(data);
  fftw_free(work);
  if (gen_image != nullptr) fftw_free(gen_image);
  if (posgen != nullptr) delete posgen;
  delete[] gauss;
  if (gauss_add != nullptr) delete[] gauss_add;
  if (binval != nullptr) delete[] binval;
}

bool simImage::isValid() const {
  if (n1 == 0) return false;
  if (n2 == 0) return false;
  if (pixsize <= 0.0) return false;
  if (oversample == 0) return false;
  if (fwhm <= 0.0) return false;
  if (sigi < 0.0) return false;
  if (esmooth < 0.0) return false;
  return true;
}

/*!
  Downsample the internal image into data.
  This rebinning preserves the mean flux per pixel
  
 */
void simImage::downSample(unsigned int ni1, unsigned int ni2, 
                          double* const inarr, unsigned int no1, 
                          unsigned int no2, double* const outarr) {
  unsigned int osamp = ni1 / no1;
  if (osamp == 0)
    pofdExcept("simImage", "downSample", "osamp is invalid (0)", 1);
  if (no2 * osamp != ni2)
    pofdExcept("simImage", "downSample", "Dimension 1 doesn't match", 2);
  if (osamp == 1) {
    for (unsigned int i = 0; i < no1 * no2; ++i) outarr[i] = inarr[i];
    return;
  }

  double sum;
  unsigned int minxidx, maxxidx, minyidx, maxyidx;
  double *in_rowptr, *out_rowptr;
  double corrfac = 1.0 / static_cast<double>(osamp*osamp);
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

void simImage::convolveInner(unsigned int n, const double* const arr,
                             unsigned int ni1, unsigned int ni2,
                             double* const inarr,
                             double* const outarr) const {
  //inarr, outarr must be same size (ni1 by ni2)
  //This makes use of the work array to hold the intermediate product
  //In use, the input array may be data (if we aren't oversampling) or 
  // gen_image (if we are oversampling)

  //Here we take advantage of the fact that a Gaussian beam factorizes
  //to do this as two 1D convolutions, along rows then columns.  This is
  // a large speedup

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
  // either into data (if we aren't oversampling) or gen_image (if we are).
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

void simImage::convolveWithBeamInPlace(unsigned int n1, unsigned int n2,
                                       double* const inout) {
  convolveInner(ngauss, gauss, n1, n2, inout, inout);
}

void simImage::convolveWithBeam(unsigned int ni1, unsigned int ni2,
                                double* const input, unsigned int no1,
                                unsigned int no2, double* const output) {
  //This is the general version -- it can handle downsampling if needed
  if (ni1 != no1 || ni2 != no2) {
    //From input to work2 via work
    convolveInner(ngauss, gauss, ni1, ni2, input, gen_image);
    //From gen_work2 to output
    downSample(ni1, ni2, gen_image, no1, no2, output);
  } else {
    //Simpler case -- from input to output via work
    convolveInner(ngauss, gauss, ni1, ni2, input, output);
  }
}

void simImage::convolveWithAdd() {
  //This can only be done to data
  if (!is_full)
    throw pofdExcept("simImage", "convolveWithAdd",
                     "Trying to convolve empty image", 1);

  //From data to data, no oversampling
  convolveInner(ngauss_add, gauss_add, n1, n2, data, data);

  //Normalization step; the idea is to keep the peak flux values
  // the same for sources.  This depends both on the beam size
  // and the extra amount of smoothing, and is derived from the
  // relation between the peak value of a Gaussian beam and it's
  // area
  const double prefac = 4 * std::log(2)/pofd_coverage::pi;
  double fwhm2 = fwhm * fwhm;
  double esmooth2 = esmooth * esmooth;
  double normval = prefac * (fwhm2 + esmooth2) * pixsize * pixsize / 
    (fwhm2 + esmooth2);
  for (unsigned int i = 0; i < n1 * n2; ++i)
    data[i] *= normval;
}

/*
  \param[in] ntrials Number of trials to do if measuring filtered noise.
  \param[in] filt fourier space filter to apply
  \returns Instrument noise including effects of filtering, etc.
  
  This can include both filtering and additional Gaussian smoothing.
  Note that it is not necessary for the data array to be filled
  for this function to work, but it can take additional storage
  equal to the size of the image.  

  This is not particularly cheap the first time it is called.
  After that it uses the previously computed value if possible.  
*/
double simImage::getFinalNoise(unsigned int ntrials,
                               const fourierFilter* const filt) const {

  // esmooth normalization factor input
  const double prefac = 4 * std::log(2) / pofd_coverage::pi;
  const double noise_prefac = sqrt(2 * std::log(2) / pofd_coverage::pi);

  // First decide if we can re-use the previous value
  bool recompute = true;
  if (sigi_final_computed) {
    if (filt == nullptr) recompute = false; // Ignore ntrials -- not needed
    else if (ntrials <= sigi_final_ntrials) recompute = false;
  }
  if (ntrials == 0)
    throw pofdExcept("simImage", "getFinalNoise",
                     "Invalid (non-positive) ntrials", 1);

  // If necessary, do the computation
  if (recompute) {
    if (sigi == 0) {
      sigi_final_computed = true;
      sigi_final_ntrials = ntrials; // Well, not really...
      sigi_final = 0.0;
    } else if (filt == nullptr) {
      if (esmooth <= 0.0) {
        sigi_final_computed = true;
        sigi_final_ntrials = ntrials;
        sigi_final = sigi;
      } else {
        // Only extra Gaussian smoothing.  Can be done analytically
        double fwhm2 = fwhm * fwhm;
        sigi_final_computed = true;
        sigi_final_ntrials = ntrials;
        sigi_final = noise_prefac * sigi * pixsize * 
          (fwhm2 + esmooth * esmooth) / (esmooth * fwhm2);
      }
    } else {
      // Filtering, and maybe smoothing.  
      // We have to do this by making a simulated image, filling it with 
      // noise, etc.  This is why this function may not be cheap,
      // and why we support multiple trials
      double* tmpdata = (double*) fftw_malloc(sizeof(double) * n1 * n2);

      // Compute esmooth prefactor if needed
      double normval = 1.0;
      if (esmooth > 0) {
        double fwhm2 = fwhm * fwhm;
        double esmooth2 = esmooth * esmooth;
        normval = prefac * (fwhm2 + esmooth2) * pixsize * pixsize / 
          (fwhm2 + esmooth2);
      }

      double var, var1, var2, mn, dtmp;
      var = 0.0;
      for (unsigned int idx = 0; idx < ntrials; ++idx) {
        // Fill with random noise
        for (unsigned int i = 0; i < n1 * n2; ++i)
          tmpdata[i] = sigi * rangen.gauss();

        // Apply additional Gaussian smoothing if needed
        if (esmooth > 0) {
          convolveInner(ngauss_add, gauss_add, n1, n2, tmpdata, tmpdata);
          // Then renormalize
          for (unsigned int i = 0; i < n1 * n2; ++i)
            tmpdata[i] *= normval;
        }

        // Now filtering
        filt->filter(n1, n2, pixsize, tmpdata);

        // Measure using two pass algorithm
        mn = tmpdata[0];
        for (unsigned int i = 1; i < n1 * n2; ++i)
          mn += tmpdata[i];
        mn /= static_cast<double>(n1 * n2);
        var1 = var2 = 0.0;
        for (unsigned int i = 0; i < n1 * n2; ++i) {
          dtmp = tmpdata[i] - mn;
          var1 += dtmp * dtmp;
          var2 += dtmp;
        }
        var1 -= var2 * var2 / static_cast<double>(n1 * n2);
        var += var1 / static_cast<double>(n1 * n2 - 1);
      }
      if (ntrials > 1)
        var /= static_cast<double>(ntrials);

      fftw_free(tmpdata);
      sigi_final_computed = true;
      sigi_final_ntrials = ntrials;
      sigi_final = sqrt(var);
    }
  }

  return sigi_final;
}


double simImage::getArea() const {
  double val = pixsize / 3600.0;
  return val * val * n1 * n2;
}

/*!
  Generates a simulated image
  \param[in] model Base number counts model
  \param[in] n0 Number of sources per area to generate
  \param[in] meansub Do mean subtraction.  Note that filtering
             automatically results in mean subtraction.
  \param[in] Fourier space filter to apply.  If nullptr, no filtering.
  \param[in] bin Create binned image data
  \param[in] sparsebin Only take every this many pixels in binned image.
                         Does nothing if no binning.
*/
void simImage::realize(const numberCounts& model, double n0,
                       bool meansub, const fourierFilter* const filt,
                       bool bin, unsigned int sparsebin) {

  if (!isValid())
    throw pofdExcept("simImage", "realize",
                     "Trying to realize image with invalid parameters", 1);
  if (!model.isValid())
    throw pofdExcept("simImage", "realize",
                     "Trying to realize model with invalid parameters", 2);
  if (n0 <= 0.0)
    throw pofdExcept("simImage", "realize", "n0 must be positive", 3);

  double area = getArea();
  unsigned int nsrcs = static_cast<unsigned int>(area * n0);
  
  // Set up the position generator if it will be needed
  if (use_clustered_pos) posgen->generate(rangen);

  //Inject sources
  if (oversample > 1) {
    //Generate in oversampled gen_image
    std::memset(gen_image, 0, ngen1 * ngen2 * sizeof(double));
    if (nsrcs > 0) {
      unsigned int idx1, idx2;
      if (use_clustered_pos) {
        // Clustered positions
        std::pair<unsigned int, unsigned int> pos;
        for (unsigned int i = 0; i < nsrcs; ++i) {
          pos = posgen->getPosition(rangen, oversample);
          idx1 = pos.first;
          idx2 = pos.second;
          gen_image[idx1 * ngen2 + idx2] += model.genSource(rangen.doub());
        }
      } else {
        // Uniform distribution -- easy!
        for (unsigned int i = 0; i < nsrcs; ++i) {
          idx1 = static_cast<unsigned int>(rangen.doub() * ngen1);
          idx2 = static_cast<unsigned int>(rangen.doub() * ngen2);
          gen_image[idx1 * ngen2 + idx2] += model.genSource(rangen.doub());
        }
      }

      // This also moves the data from gen_image to the data array
      // and downsamples
      convolveWithBeam(ngen1, ngen2, gen_image, n1, n2, data);
    }
  } else {
    //Generate in data
    std::memset(data, 0, n1 * n2 * sizeof(double));

    //Inject sources, much like above except no downsampling
    if (nsrcs > 0) {
      unsigned int idx1, idx2;
      if (use_clustered_pos) {
        std::pair<unsigned int, unsigned int> pos;
        for (unsigned int i = 0; i < nsrcs; ++i) {
          pos = posgen->getPosition(rangen);
          idx1 = pos.first;
          idx2 = pos.second;
          data[idx1 * ngen2 + idx2] += model.genSource(rangen.doub());
        }
      } else {
        for (unsigned int i = 0; i < nsrcs; ++i) {
          idx1 = static_cast<unsigned int>(rangen.doub() * n1);
          idx2 = static_cast<unsigned int>(rangen.doub() * n2);
          data[idx1 * n2 + idx2] += model.genSource(rangen.doub());
        }
      }
    }
    convolveWithBeamInPlace(n1, n2, data);
  }
  is_full = true;

  //Add instrument noise
  if (sigi > 0.0)
    for (unsigned int i = 0; i < n1 * n2; ++i)
      data[i] += sigi * rangen.gauss();

  //Extra smoothing, if set
  if (esmooth > 0.0 && ((nsrcs > 0) || (sigi > 0.0)))
    convolveWithAdd();

  // Apply filtering.  Note this is done after adding noise and
  //  downsampling to the final resolution (if oversampling is used).
  // Filtering will always result in mean subtraction, so no need to do
  //  that twice
  isHipass = isMatched = false;
  filtscale = qfactor = matched_fwhm = matched_sigi = matched_sigc =
        std::numeric_limits<double>::quiet_NaN();
  if (filt != nullptr) {
    filt->filter(n1, n2, pixsize, data);
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
  } else if (meansub) meanSubtract();

  is_binned = false;
  if (bin) applyBinning(sparsebin);
}



/*!
  \param[in] sparsebin Only take every this many pixels when binning. 1 means
                        fully sampled (0 is also interpreted to mean the same
                        thing)

  Keeps the original, unbinned image around as well
*/
void simImage::applyBinning(unsigned int sparsebin) {
  if (!is_full)
    throw pofdExcept("simImage", "applyBinning",
                     "Trying to bin empty image", 1);

  if (nbins == 0) throw pofdExcept("simImage", "applyBinning",
                                   "Trying to bin with no bins", 2);

  if (sparsebin >= n1 * n2) 
    throw pofdExcept("simImage", "applyBinning",
                     "Sparse binning factor larger than simulated image", 3);

  //Only allocated the first time we bin
  //There is no way to change nbins after initialization
  if (binval == nullptr)
    binval = new unsigned int[nbins];
  std::memset(binval, 0, nbins * sizeof(unsigned int));

  //First, we need the minimum and maximum
  double minval, maxval;
  getMinMax(minval, maxval);
  
  //We want to put the max and min in the center of the top and
  // bottom bin
  bincent0 = minval;
  if (nbins == 1)
    bindelta = 2 * (maxval - minval);
  else
    bindelta = (maxval - minval) / static_cast<double>(nbins - 1);

  //And... bin
  double ibindelta = 1.0 / bindelta;
  unsigned int idx;
  if (sparsebin > 1) bin_sparcity = sparsebin; else bin_sparcity = 1;
  for (unsigned int i = 0; i < n1*n2; i += bin_sparcity) {
    idx = static_cast<unsigned int>((data[i] - bincent0) * ibindelta + 0.5);
    binval[idx] += 1;
  } 
  is_binned = true;
}

double simImage::meanSubtract() {
  if (!is_full)
    throw pofdExcept("simImage", "meanSubtract", "No data", 1);
  double mn = getMean();
  for (unsigned int i = 0; i < n1*n2; ++i)
    data[i] -= mn;
  if (is_binned)
    bincent0 -= mn;
  return mn;

}

void simImage::getMinMax(double& min, double& max) const {
  if (! is_full )
    throw pofdExcept("simImage", "getMinMax",
                     "Trying to getMinMax of empty image", 1);

  min = data[0];
  max = data[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    if (data[i] < min) min = data[i];
  for (unsigned int i = 1; i < n1*n2; ++i)
    if (data[i] > max) max = data[i];
}

double simImage::getMean() const {
  if (!is_full)
    throw pofdExcept("simImage", "getMean",
                     "Trying to get mean of empty image", 1);
  double norm = 1.0/static_cast<double>(n1*n2);
  double mn = data[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    mn += data[i];
  mn *= norm;
  return mn;
}

void simImage::getMeanAndVar(double& mn, double& var) const {
  if (!is_full)
    throw pofdExcept("simImage", "getMeanAndVar",
                     "Trying to get mean and var of empty image", 1);
  //Use corrected two pass algorithm
  double norm = 1.0/static_cast<double>(n1 * n2);
  mn = data[0];
  for (unsigned int i = 1; i < n1 * n2; ++i)
    mn += data[i];
  mn *= norm;

  double sum1, sum2, tmp;
  tmp = data[0]-mn;
  sum1 = tmp * tmp;
  sum2 = tmp;
  for (unsigned int i = 1; i < n1 * n2; ++i) {
    tmp = data[i] - mn;
    sum1 += tmp * tmp;
    sum2 += tmp;
  }
  var = (sum1 - norm * sum2 * sum2) / static_cast<double>(n1 * n2 - 1);
}

double simImage::getBeamSum() const {
  if (ngauss == 0)
    throw pofdExcept("simImage", "getBeamSum",
                     "No beam pixels", 1);
  
  double sum1D = gauss[0];
  for (unsigned int i = 1; i < ngauss; ++i)
    sum1D += gauss[i];
  if (oversample > 1) sum1D /= static_cast<double>(oversample);
  return sum1D * sum1D;
}


double simImage::getBeamSumSq() const {
  if (ngauss == 0)
    throw pofdExcept("simImage", "getBeamSumSq",
                     "No beam pixels", 1);
  
  double tmp = gauss[0];
  double sum1D = tmp*tmp;
  for (unsigned int i = 1; i < ngauss; ++i) {
    tmp = gauss[i];
    sum1D += tmp * tmp;
  }
  if (oversample > 1) sum1D /= static_cast<double>(oversample);
  return sum1D * sum1D;
}

double simImage::getEsmooth() const {
  if (esmooth <= 0.0) return std::numeric_limits<double>::quiet_NaN();
  else return esmooth;
}

/*!
  \param[in] outputfile File to write to
  \returns 0 on success, an error code (!=0) for anything else

  Note this doesn't throw exceptions in keeping with the CFITSIO
  error handling strategy
*/
int simImage::writeToFits(const std::string& outputfile) const {

  if (!is_full)
    throw pofdExcept("simImage", "writeToFits",
                     "Trying to write image without realizing", 1);

  //Make the fits file
  int status = 0;
  fitsfile *fp;

  fits_create_file(&fp, outputfile.c_str(), &status);
  if (status) {
    fits_report_error(stderr, status);
    return status;
  }

  //Now make image array
  long axissize[2];
  axissize[0] = static_cast<long>(n1);
  axissize[1] = static_cast<long>(n2);
  fits_create_img(fp, DOUBLE_IMG, 2, axissize, &status);
  
  //Model params
  fits_write_key(fp, TSTRING, const_cast<char*>("MODEL"),
                 const_cast<char*>("Spline"), 
                 const_cast<char*>("Model type"),
                 &status);

  //Sim params
  int itmp = 0;
  double dtmp = fwhm;

  fits_write_key(fp, TDOUBLE, const_cast<char*>("FWHM"), &dtmp, 
                 const_cast<char*>("Beam fwhm [arcsec]"), 
                 &status);
  if (esmooth > 0) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &itmp,
                   const_cast<char*>("Additional smoothing applied"), 
                   &status);
    dtmp = esmooth;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("ESMOOTH"), &dtmp, 
                   const_cast<char*>("Extra smoothing fwhm [arcsec]"), 
                   &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &itmp,
                   const_cast<char*>("Additional smoothing applied"), 
                   &status);
  }

  dtmp = sigi;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGI"), &dtmp, 
                 const_cast<char*>("Raw instrument noise"), 
                 &status);

  if (sigi_final_computed) {
    dtmp = sigi_final;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGIFNL"), &dtmp, 
                   const_cast<char*>("Final instrument noise"), 
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

  if (isHipass) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT"), &itmp,
                   const_cast<char*>("Is map hipass filtered?"), &status);
    dtmp = filtscale;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTSCL"), &dtmp,
                   const_cast<char*>("Filtering scale [arcsec]"), &status);
    dtmp = qfactor;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("FLTQ"), &dtmp,
                   const_cast<char*>("Filtering apodization"), &status);
  } else {
    itmp = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("HIFLT"), &itmp,
                   const_cast<char*>("Is map hipass filtered?"), &status);
  }
  if (isMatched) {
    itmp = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("MTCHFLT"), &itmp,
                   const_cast<char*>("Is map match filtered?"), &status);
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
                   const_cast<char*>("Is map match filtered?"), &status);
  }

  fits_write_key(fp, TSTRING, const_cast<char*>("VERSION"),
                 const_cast<char*>(pofd_coverage::version), 
                 const_cast<char*>("pofd_coverage version"),
                 &status);

  fits_write_history(fp,const_cast<char*>("Simulated image from pofd_coverage"),
                     &status);
  fits_write_date(fp, &status);



  // Astrometry
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE1"),
                 const_cast<char*>("RA---TAN"),
                 const_cast<char*>("WCS: Projection type axis 1"),&status);
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE2"),
                 const_cast<char*>("DEC--TAN"),
                 const_cast<char*>("WCS: Projection type axis 2"),&status);
  dtmp = n1/2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRPIX1"), &dtmp, 
                 const_cast<char*>("Ref pix of axis 1"), &status);
  dtmp = n2/2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRPIX2"), &dtmp, 
                 const_cast<char*>("Ref pix of axis 2"), &status);
  dtmp = 90.0; //Arbitrary
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRVAL1"), &dtmp, 
                 const_cast<char*>("val at ref pix axis 1"), &status);
  dtmp = 10.0; //Arbitrary
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRVAL2"), &dtmp, 
                 const_cast<char*>("val at ref pix axis 2"), &status);
  dtmp = - pixsize/3600.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD1_1"), &dtmp, 
                 const_cast<char*>("Pixel scale axis 1,1"), &status);
  dtmp = 0.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD1_2"), &dtmp, 
                 const_cast<char*>("Pixel scale axis 1,2"), &status);
  dtmp = 0.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD2_1"), &dtmp, 
                 const_cast<char*>("Pixel scale axis 2,1"), &status);
  dtmp = pixsize/3600.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD2_2"), &dtmp, 
                 const_cast<char*>("Pixel scale axis 2,2"), &status);
  dtmp = 2000.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("EPOCH"), &dtmp, 
                 const_cast<char*>("WCS: Epoch of celestial pointing"), 
                 &status);
  fits_write_key(fp, TDOUBLE, const_cast<char*>("EQUINOX"), &dtmp, 
                 const_cast<char*>("WCS: Equinox of celestial pointing"), 
                 &status);

  //Do data writing.  We have to make a transposed copy of the
  // data to do this, which is irritating as hell
  double *tmpdata = new double[n1];
  long fpixel[2] = { 1, 1 };
  for ( unsigned int j = 0; j < n2; ++j ) {
    for (unsigned int i = 0; i < n1; ++i) tmpdata[i] = data[i * n2 + j];
    fpixel[1] = static_cast<long>(j + 1);
    fits_write_pix(fp, TDOUBLE, fpixel, n1, tmpdata, &status);
  }
  delete[] tmpdata;

  //Close up and go home
  fits_close_file(fp, &status);
  if (status) {
    fits_report_error(stderr, status);
    return status;
  }

  return status;
}

/*!
  \param[in] outfile File to write to
*/
int simImage::writeProbImageToFits(const std::string& outfile) const {
  if (!use_clustered_pos)
    throw pofdExcept("simImage", "writeProbImageToFits",
                     "No probability image to write", 1);
  return posgen->writeProbToFits(outfile);
}

/*!
  \param[in] obj_id Open HDF5 handle to write to.

  Writes information about the position generator to the open handle,
  creating a group called "PositionGenerator"
*/
void simImage::writePositionGeneratorToHDF5Handle(hid_t obj_id) const {

  if (H5Iget_ref(obj_id) < 0)
    throw pofdExcept("simImage", "writePositionGeneratorToHDF5Handle",
                     "Given non-open obj_id to write to", 1);

  hid_t group_id;
  group_id = H5Gcreate(obj_id, "PositionGenerator", H5P_DEFAULT, H5P_DEFAULT, 
                      H5P_DEFAULT);
  if (H5Iget_ref(group_id) < 0)
    throw pofdExcept("simImage", "writePositionGeneratorToHDF5",
                     "Failed to create HDF5 model group", 2);

  hsize_t adims;
  hid_t mems_id, att_id;
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);
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
