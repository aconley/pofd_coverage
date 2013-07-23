#include<ctime>
#include<cmath>
#include<limits>

#include<iostream> //For debuggin

#include<fitsio.h>

#include "../include/global_settings.h"
#include "../include/simImageDouble.h"
#include "../include/pofdExcept.h"

simImageDouble::simImageDouble(unsigned int N1, unsigned int N2, double PIXSIZE,
			       double FWHM1, double FWHM2, double SIGI1, 
			       double SIGI2, double ESMOOTH1, double ESMOOTH2,
			       unsigned int OVERSAMPLE, unsigned int NBINS) {
  n1 = N1;
  n2 = N2;
  oversample = OVERSAMPLE;
  ngen1 = n1 * oversample;
  ngen2 = n2 * oversample;
  data1 = new double[n1 * n2];
  data2 = new double[n1 * n2];
  work  = new double[ngen1 * ngen2];
  if (oversample > 1) {
    gen_1 = new double[ngen1 * ngen2];
    gen_2 = new double[ngen1 * ngen2];
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
  smooth_applied = false;
  is_binned = false;
  nbins = NBINS;
  bincent01 = 0.0;
  bindelta1 = 0.0;
  bincent02 = 0.0;
  bindelta2 = 0.0;
  binval = NULL;

  //Set RNG seed
  unsigned long long int seed;
  seed = static_cast<unsigned long long int>(time(NULL));
  seed += static_cast<unsigned long long int>(clock());
  rangen.set_seed(seed);
  
  //Set up array to hold 1D beams in each band (center normalized)
  //We have to decide how far out to go.  Note that the beams are set
  // up in oversampled space
  const double nfwhm = 3.5;
  ngauss1 = static_cast<unsigned int>(nfwhm * fwhm1 / pixsize_gen + 0.99999999);
  ngauss1 = 2 * ngauss1 + 1;
  gauss1 = new double[ngauss1];
  ngauss2 = static_cast<unsigned int>(nfwhm * fwhm2 / pixsize_gen + 0.99999999);
  ngauss2 = 2 * ngauss2 + 1;
  gauss2 = new double[ngauss2];
  doublebeam bm(fwhm1, fwhm2);
  bm.getBeamFac1(ngauss1, pixsize_gen, gauss1);
  bm.getBeamFac2(ngauss2, pixsize_gen, gauss2);

  //Set up additional smoothing 
  // This is done in un-oversampled space
  if (esmooth1 > 0 || esmooth2 > 0) {
    doublebeam ebm(esmooth1, esmooth2);
    if (esmooth1 > 0) {
      ngauss_add1 = static_cast<unsigned int>(nfwhm * esmooth1 / pixsize +
					      0.99999999);
      ngauss_add1 = 2 * ngauss_add1 + 1;
      gauss_add1 = new double[ngauss_add1];
      ebm.getBeamFac1(ngauss_add1, pixsize, gauss_add1);
    } else {
      ngauss_add1 = 0;
      gauss_add1 = NULL;
    }
    if (esmooth2 > 0) {
      ngauss_add2 = static_cast<unsigned int>(nfwhm * esmooth2 / pixsize +
					      0.99999999);
      ngauss_add2 = 2 * ngauss_add2 + 1;
      gauss_add2 = new double[ngauss_add2];
      ebm.getBeamFac2(ngauss_add2, pixsize, gauss_add2);
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
  delete[] data1;
  delete[] data2;
  delete[] work;
  if (gen_1 != NULL) delete[] gen_1;
  if (gen_2 != NULL) delete[] gen_2;
  delete[] gauss1;
  delete[] gauss2;
  if (ngauss_add1 > 0) delete[] gauss_add1;
  if (ngauss_add2 > 0) delete[] gauss_add2;
  if (binval != NULL) delete[] binval;
}

bool simImageDouble::isValid() const {
  if (n1 == 0) return false;
  if (n2 == 0) return false;
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
    pofdExcept("simImageDouble","downSample","osamp is invalid (0)",1);
  if (no2 * osamp != ni2)
    pofdExcept("simImageDouble","downSample","Dimension 1 doesn't match",2);
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
				   double* const outarr) {
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
    throw pofdExcept("simImageDouble","convolveWithBeam",
		     "Trying to convolve empty image",1);
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
  smooth_applied = false;
}

void simImageDouble::convolveWithAdd() {
  //This can only be done to data1, data2
  if (! is_full )
    throw pofdExcept("simImageDouble","convolveWithAdd",
		     "Trying to convolve empty image",1);

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
    smooth_applied = true;
  }

  if (esmooth2 > 0) {
    //Same, band 2
    convolveInner(ngauss_add2, gauss_add2, n1, n2, data2, data2);

    const double prefac = 4*std::log(2)/pofd_coverage::pi;
    double normval = prefac * (fwhm2*fwhm2 + esmooth2*esmooth2) * 
      pixsize*pixsize / (fwhm2*fwhm2*esmooth2*esmooth2);
    for (unsigned int i = 0; i < n1*n2; ++i)
      data2[i] *= normval;
    smooth_applied = true;
  }

}


std::pair<double,double> simImageDouble::getNoise() const {
  if (smooth_applied) return getSmoothedNoiseEstimate(); else
    return std::make_pair(sigi1, sigi2);
}

//Will return estimated smoothed noise, even if current image is not smoothed
std::pair<double,double> simImageDouble::getSmoothedNoiseEstimate() const {
  const double prefac = sqrt(2*std::log(2)/pofd_coverage::pi);
  double s1, s2;

  if ( sigi1 == 0.0 ) s1 = 0.0;
  else if (esmooth1 <= 0.0) s1 = sigi1;
  else s1 = prefac*sigi1*pixsize*(fwhm1*fwhm1+esmooth1*esmooth1) / 
	 ( esmooth1 * fwhm1*fwhm1 );

  if ( sigi2 == 0.0 ) s2 = 0.0;
  else if (esmooth2 <= 0.0) s2 = sigi2;
  else s2 = prefac*sigi2*pixsize*(fwhm2*fwhm2+esmooth2*esmooth2) / 
	 ( esmooth2 * fwhm2*fwhm2 );

  return std::make_pair(s1, s2);

}

double simImageDouble::getArea() const {
  return pixsize*pixsize*n1*n2/(3600.0*3600.0);
}

/*!
  Generates a simulated image
  \params[in] model Base number counts model
  \params[in] n0 Number of sources per area to generate
  \params[in] extra_smooth Apply additional Gaussian smoothing
  \params[in] meansub Do mean subtraction
  \params[in] bin Create binned image data
 */
void simImageDouble::realize(const numberCountsDouble& model,
			     double n0, bool extra_smooth, bool meansub, 
			     bool bin) {

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

  //Inject sources
  std::pair<double, double> src;
  if (oversample > 1) {
    //Generate in oversampled gen_1, gen_2
    for (unsigned int i = 0; i < ngen1 * ngen2; ++i) gen_1[i] = 0.0;
    for (unsigned int i = 0; i < ngen1 * ngen2; ++i) gen_2[i] = 0.0;
    if (nsrcs > 0) {
      unsigned int idx1, idx2, combidx;
      for (unsigned int i = 0; i < nsrcs; ++i) {
	idx1 = static_cast<unsigned int>(rangen.doub() * ngen1);
	idx2 = static_cast<unsigned int>(rangen.doub() * ngen2);
	src = model.genSource(rangen.doub(), rangen.gauss());
	combidx = idx1 * ngen2 + idx2;
	gen_1[combidx] += src.first;
	gen_2[combidx] += src.second;
      }
    }
  } else {
    //Generate in data1, data2
    for (unsigned int i = 0; i < n1 * n2; ++i) data1[i] = 0.0;
    for (unsigned int i = 0; i < n1 * n2; ++i) data2[i] = 0.0;

    //Inject sources
    if (nsrcs > 0) {
      unsigned int idx1, idx2, combidx;
      for (unsigned int i = 0; i < nsrcs; ++i) {
	idx1 = static_cast<unsigned int>(rangen.doub() * n1);
	idx2 = static_cast<unsigned int>(rangen.doub() * n2);
	src = model.genSource(rangen.doub(), rangen.gauss());
	combidx = idx1 * ngen2 + idx2;
	data1[combidx] += src.first;
	data2[combidx] += src.second;
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
  if ( extra_smooth && ( (nsrcs > 0) || (sigi1 > 0.0 || sigi2 > 0.0) ) ) {
    if (ngauss_add1 == 0)
       throw pofdExcept("simImageDouble","realize",
			"Trying to apply extra smoothing without setting esmooth1",2);
    if (ngauss_add2 == 0)
       throw pofdExcept("simImageDouble","realize",
			"Trying to apply extra smoothing without setting esmooth2",2);
    convolveWithAdd();
  }


  //Mean subtract
  if (meansub) meanSubtract();

  //bin
  if (bin) applyBinning();

}

std::pair<double,double> simImageDouble::meanSubtract() {
  double mn1, mn2;
  getMean(mn1, mn2);
  for (unsigned int i = 0; i < n1*n2; ++i)
    data1[i] -= mn1;
  for (unsigned int i = 0; i < n1*n2; ++i)
    data2[i] -= mn2;
  std::pair<double,double> ret;
  ret.first = mn1;
  ret.second = mn2;
  if (is_binned) {
    bincent01 -= mn1;
    bincent02 -= mn2;
  }
  return ret;
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

void simImageDouble::getMean(double& mn1, double& mn2) const {
  if (!is_full)
    throw pofdExcept("simImageDouble","getMean",
		     "Trying to get means of empty images",1);
  double norm = 1.0/static_cast<double>(n1 * n2);
  mn1 = data1[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    mn1 += data1[i];
  mn1 *= norm;
  mn2 = data2[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    mn2 += data2[i];
  mn2 *= norm;
}

void simImageDouble::getMeanAndVar(double& mn1, double& var1,
				   double& mn2, double& var2) const {
  if (!is_full)
    throw pofdExcept("simImageDouble","getMeanAndVar",
		     "Trying to get mean and vars of empty images",1);
  //Use corrected two pass algorithm
  double norm = 1.0/static_cast<double>(n1*n2);
  mn1 = data1[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    mn1 += data1[i];
  mn1 *= norm;

  double sum1, sum2, tmp;
  tmp = data1[0]-mn1;
  sum1 = tmp*tmp;
  sum2 = tmp;
  for (unsigned int i = 1; i < n1*n2; ++i) {
    tmp = data1[i]-mn1;
    sum1 += tmp*tmp;
    sum2 += tmp;
  }
  var1 = (sum1 - norm*sum2*sum2)/static_cast<double>(n1*n2-1);

  mn2 = data2[0];
  for (unsigned int i = 1; i < n1*n2; ++i)
    mn2 += data2[i];
  mn2 *= norm;
  tmp = data2[0]-mn2;
  sum1 = tmp*tmp;
  sum2 = tmp;
  for (unsigned int i = 1; i < n1*n2; ++i) {
    tmp = data2[i]-mn2;
    sum1 += tmp*tmp;
    sum2 += tmp;
  }
  var2 = (sum1 - norm*sum2*sum2)/static_cast<double>(n1*n2-1);
}

std::pair<double,double> simImageDouble::getBeamSum() const {
  if (ngauss1 == 0)
    throw pofdExcept("simImageDouble","getBeamSum",
		     "No beam pixels, band 1",1);
  if (ngauss2 == 0)
    throw pofdExcept("simImageDouble","getBeamSum",
		     "No beam pixels, band 2",2);
  
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
    throw pofdExcept("simImageDouble","getBeamSumSq",
		     "No beam pixels, band 1",1);
  if (ngauss2 == 0)
    throw pofdExcept("simImageDouble","getBeamSumSq",
		     "No beam pixels, band 2",2);
  
  double tmp = gauss1[0];
  double sum1D_1 = tmp*tmp;
  for (unsigned int i = 1; i < ngauss1; ++i) {
    tmp = gauss1[i];
    sum1D_1 += tmp*tmp;
  }

  double sum1D_2 = tmp*tmp;
  for (unsigned int i = 1; i < ngauss2; ++i) {
    tmp = gauss2[i];
    sum1D_2 += tmp*tmp;
  }

  return std::make_pair(sum1D_1 * sum1D_1, sum1D_2 * sum1D_2);
}

/*!
  Keeps the original, unbinned image around as well
 */
void simImageDouble::applyBinning() {
  if (!is_full)
    throw pofdExcept("simImageDouble","applyBinning",
		     "Trying to bin empty image",1);

  if (nbins == 0) throw pofdExcept("simImageDouble","applyBinning",
				   "Trying to bin with no bins",2);

  //Only allocated the first time we bin
  //There is no way to change nbins after initialization
  if (binval == NULL)
    binval = new unsigned int[nbins*nbins];
  for (unsigned int i = 0; i < nbins*nbins; ++i)
    binval[i] = 0;

  //First, we need the minimum and maximum
  double minval1, maxval1, minval2, maxval2;
  getMinMax(minval1,maxval1,minval2,maxval2);
  
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
  for (unsigned int i = 0; i < n1*n2; ++i) {
    idx1 = static_cast<unsigned int>( (data1[i]-bincent01)*ibindelta1 + 0.5 );
    idx2 = static_cast<unsigned int>( (data2[i]-bincent02)*ibindelta2 + 0.5 );
    binval[idx1*nbins+idx2] += 1;
  }
  is_binned = true;
}



/*!
  \param[in] outputfile File to write to
  \returns 0 on success, an error code (!=0) for anything else

  Note this doesn't throw exceptions in keeping with the CFITSIO
  error handling strategy
*/
int simImageDouble::writeToFits(const std::string& outputfile) const {

  //Make the fits file
  int status = 0;
  fitsfile *fp;

  fits_create_file(&fp, outputfile.c_str(), &status);

  if (status) {
    fits_report_error(stderr,status);
    return status;
  }

  //Stuff for the primary header
  int tval = 1;
  fits_write_key(fp, TLOGICAL, const_cast<char*>("SIMPLE"),
		 &tval, const_cast<char*>("Primary Header"),&status);
  tval = 8;
  fits_write_key(fp, TINT, const_cast<char*>("BITPIX"),
		 &tval, NULL,&status);
  tval = 0;
  fits_write_key(fp, TINT, const_cast<char*>("NAXIS"),
		 &tval, NULL,&status);
  tval = 1;
  fits_write_key(fp, TLOGICAL, const_cast<char*>("Extend"),
		 &tval, const_cast<char*>("Extensions may be present"),
		 &status);
  fits_write_history( fp, const_cast<char*>("Simulated image from pofd_coverage"),
		      &status);
  fits_write_date(fp, &status);

  //Model params
  double tmpval;
  fits_write_key(fp, TSTRING, const_cast<char*>("MODEL"),
		 const_cast<char*>("Spline-LogNormal"), 
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
  if (smooth_applied) {
    int tmpi = 1;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &tmpi,
		   const_cast<char*>("Additional smoothing applied"), 
		   &status);
    tmpval = esmooth1;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("ESMOOTH1"), &tmpval, 
		   const_cast<char*>("Extra smoothing fwhm, band 1 [arcsec]"), 
		   &status);
    tmpval = esmooth2;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("ESMOOTH1"), &tmpval, 
		   const_cast<char*>("Extra smoothing fwhm, band 2 [arcsec]"), 
		   &status);
  } else {
    int tmpi = 0;
    fits_write_key(fp, TLOGICAL, const_cast<char*>("ADDSMTH"), &tmpi,
		   const_cast<char*>("Additional smoothing applied"), 
		   &status);
  }

  tmpval = sigi1;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGI_1"), &tmpval, 
		 const_cast<char*>("Instrument noise, band 1"), 
		 &status);
  tmpval = sigi2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGI_2"), &tmpval, 
		 const_cast<char*>("Instrument noise, band 1"), 
		 &status);
  
  if (smooth_applied) {
    double tmpval2;
    std::pair<double, double> ns;
    ns = getSmoothedNoiseEstimate();
    tmpval = ns.first;
    tmpval2 = ns.second;
    fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGISM1"), &tmpval, 
		   const_cast<char*>("Smoothed instrument noise, band 1"), 
		   &status);
    fits_write_key(fp, TDOUBLE, const_cast<char*>("SIGISM2"), &tmpval2, 
		   const_cast<char*>("Smoothed instrument noise, band 2"), 
		   &status);
  }
  fits_write_key(fp, TSTRING, const_cast<char*>("VERSION"),
		 const_cast<char*>(pofd_coverage::version), 
		 const_cast<char*>("pofd_coverage version"),
		 &status);
  fits_write_history(fp, 
		     const_cast<char*>("Simulated image from pofd_coverage"),
		     &status);
  fits_write_date(fp, &status);

  //Now make two image arrays, one for each band
  long axissize[2];
  axissize[0] = static_cast<long>(n1);
  axissize[1] = static_cast<long>(n2);
  
  fits_create_img(fp, DOUBLE_IMG, 2, axissize, &status);
  
  //Header stuff for this image
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE1"),
		 const_cast<char*>("RA---TAN"),
		 const_cast<char*>("WCS: Projection type axis 1"),&status);
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE2"),
		 const_cast<char*>("DEC--TAN"),
		 const_cast<char*>("WCS: Projection type axis 2"),&status);
  tmpval = n1/2;
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX1"), &tmpval, 
		 const_cast<char*>("Ref pix of axis 1"), &status);
  tmpval = n2/2;
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX2"), &tmpval, 
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

  //Do data writing.  We have to make a transposed copy of the
  // data to do this, which is irritating as hell
  double *tmpdata = new double[ n1 ];
  long fpixel[2] = { 1, 1 };
  for ( unsigned int j = 0; j < n2; ++j ) {
    for (unsigned int i = 0; i < n1; ++i) tmpdata[i] = data1[i*n2+j];
    fpixel[1] = static_cast<long>(j+1);
    fits_write_pix(fp, TDOUBLE, fpixel, n1, tmpdata, &status);
  }

  //Second data array, band 2
  fits_create_img(fp, DOUBLE_IMG, 2, axissize, &status);
  
  //Header stuff for this image
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE1"),
		 const_cast<char*>("RA---TAN"),
		 const_cast<char*>("WCS: Projection type axis 1"),&status);
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE2"),
		 const_cast<char*>("DEC--TAN"),
		 const_cast<char*>("WCS: Projection type axis 2"),&status);
  tmpval = n1/2;
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX1"), &tmpval, 
		 const_cast<char*>("Ref pix of axis 1"), &status);
  tmpval = n2/2;
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX2"), &tmpval, 
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

  //Data
  for ( unsigned int j = 0; j < n2; ++j ) {
    for (unsigned int i = 0; i < n1; ++i) tmpdata[i] = data2[i*n2+j];
    fpixel[1] = static_cast<long>(j+1);
    fits_write_pix(fp, TDOUBLE, fpixel, n1, tmpdata, &status);
  }
  delete[] tmpdata;

  fits_close_file(fp, &status);

  if (status) {
    fits_report_error(stderr,status);
    return status;
  }
  return status;
}
