//PD.cc
#include<limits>

#include<fitsio.h>

#include "../include/PD.h"
#include "../include/pofdExcept.h"

const double PD::lowsigval = 3.0;

PD::PD(unsigned int N, FFTFLOAT MINFLUX, FFTFLOAT DFLUX, bool LOG) : 
  n(N), capacity(N), logflat(LOG), minflux(MINFLUX), dflux(DFLUX) {
  if (capacity == 0) pd_ = NULL; else
    pd_ = (FFTFLOAT *) FFTWMALLOC(sizeof(FFTFLOAT) * capacity);
}

PD::~PD() {
  if (pd_ != NULL) FFTWFREE(pd_);
}

/*!
  Generally doesn't preserve data
 */
void PD::resize(unsigned int N) {
  //Doesn't actually resize arrays if it can avoid it
  unsigned int newcap = N;
  if ( newcap > capacity ) {
    if (pd_ != NULL) FFTWFREE(pd_);
    if (newcap > 0) 
      pd_ = (FFTFLOAT *) FFTWMALLOC(sizeof(FFTFLOAT) * newcap);
    else pd_ = NULL;
    capacity = newcap;
  }
  n = N;
}

/*!
  Tries to preserve data
 */
void PD::shrink() {
  unsigned int newcap = n;
  if ( newcap < capacity ) {
    if (newcap > 0) {
      FFTFLOAT* tmp = (FFTFLOAT*) FFTWMALLOC(sizeof(FFTFLOAT) * newcap);
      for (unsigned int i = 0; i < newcap; ++i)
	tmp[i] = pd_[i];
      if (pd_ != NULL) FFTWFREE(pd_);
      pd_ = tmp;
    } else {
      if (pd_ != NULL) FFTWFREE(pd_);
      pd_ = NULL;
    }
    capacity = newcap;
  }
}

/*!
  Generally doesn't preserve data
 */
void PD::strict_resize(unsigned int N) {
  unsigned int newcap = n;
  if ( newcap != capacity ) {
    if (pd_ != NULL) FFTWFREE(pd_);
    if (newcap > 0) pd_ = (FFTFLOAT*) FFTWMALLOC(sizeof(FFTFLOAT) * newcap);
    else pd_ = NULL;
    capacity = newcap;
  }
  n = N;
}

double PD::getTotal() const {
  if (n == 0)
    return std::numeric_limits<double>::quiet_NaN();
  double retval;
  if (logflat) {
    retval = exp2(pd_[0]);
    for (unsigned int i = 1; i < n; ++i)
      retval += exp2(pd_[i]);
  } else {
    retval = pd_[0];
    for (unsigned int i = 1; i < n; ++i)
      retval += pd_[i];
  }
  return retval;
}

double PD::getIntegral() const {
  if (n == 0)
    return std::numeric_limits<double>::quiet_NaN();
  
  double tot;
  if (logflat) {
    tot = 0.5 * exp2(pd_[0]);
    for (unsigned int i = 1; i < n-1; ++i)
      tot += exp2(pd_[i]);
    tot += 0.5 * exp2(pd_[n-1]);
  } else {
    tot = 0.5 * pd_[0];
    for (unsigned int i = 1; i < n-1; ++i)
      tot += pd_[i];
    tot += 0.5 * pd_[n-1];
  }
  return tot*dflux;
}

/*!
  Normalize the P(D), using the trapezoidal rule
  to integrate
 */
void PD::normalize() {
  if (n == 0)
    throw pofdExcept("PD","normalize",
		     "No information present to normalize",1);
  double tot = getIntegral();
  unsigned int sz = n;
  if (logflat) {
    FFTFLOAT lgtot = static_cast<FFTFLOAT>(log2(tot));
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] -= lgtot;
  } else {
    FFTFLOAT itot = static_cast<FFTFLOAT>(1.0 / tot);
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] *= itot;
  }
}

void PD::applyLog(bool nocheck) {
  if (logflat) return;
  unsigned int sz = n;
  if (nocheck) {
    for (unsigned int i = 0; i < sz; ++i) 
      pd_[i] = log2(pd_[i]);
  } else {
    FFTFLOAT val;
    for (unsigned int i = 0; i < sz; ++i) {
      val = pd_[i];
      if (val > 0.0) pd_[i] = log2(val);
      else pd_[i] = pofd_coverage::smalllogval;
    }
  }
  logflat = true;
}


void PD::deLog() {
  if (!logflat) return;
  unsigned int sz = n;
  for (unsigned int i = 0; i < sz; ++i)
    pd_[i] = exp2(pd_[i]);
  logflat = false;
}


/*
  \param[in] donorm Do not assume P(D) is normalized
 */
void PD::edgeFix(bool donorm) {
  //Compute mean and stdev
  if (n < 3) return; //No point

  if (logflat)
    throw pofdExcept("PD", "edgeFix",
		     "Not supported for logged PDs", 1);

  //Get mean and vars
  double mn, var;
  getMeanAndVar(mn, var, donorm);
  if (std::isnan(mn) || std::isinf(mn) ||
      std::isnan(var) || std::isinf(var))
    throw pofdExcept("PD", "edgeFix",
		     "Problem with means/vars", 2);
  
  double istdev = 1.0 / sqrt(var);

  //Figure out what indexes these represent in x and y
  double maxfluxfix, ddflux, dminflux;
  int maxidx;
  ddflux = static_cast<double>(dflux);
  dminflux = static_cast<double>(minflux);
  maxfluxfix = mn - PD::lowsigval * sqrt(var);
  maxidx = static_cast<int>((maxfluxfix - dminflux) / ddflux);
  maxfluxfix = dminflux + maxidx * ddflux;
  
  double tval, stepfac, subfac, pdval, preconst;
  if (maxidx > 1) {
    pdval = static_cast<double>(pd_[maxidx]);
    tval = (maxfluxfix - mn) * istdev;
    preconst = pdval * exp(0.5 * tval * tval);
    subfac = (dminflux - mn) * istdev;
    stepfac = ddflux * istdev;
    for (int i = 0; i < maxidx; ++i) {
      tval = subfac + i * stepfac;
      pd_[i] = static_cast<FFTFLOAT>(preconst * exp(-0.5 * tval * tval));
    }
  }
}


/*!
  \param[out] mean mean of P(D)
  \param[in] donorm Do not assume that P(D) is normalized.
 */
void PD::getMean(double& mean, bool donorm) const {
  if (n == 0) {
    mean = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  
  //We use the trapezoidal rule here for the integrals
  // so it isn't quite a simple sum
  mean = 0.0;

  if (logflat) {
    //i=0 has weight 0
    for (unsigned int i = 1; i < n-1; ++i)
      mean += static_cast<double>(i)*exp2(pd_[i]);
    mean += 0.5*static_cast<double>(n-1)*exp2(pd_[n-1]);
  } else {
    for (unsigned int i = 1; i < n-1; ++i)
      mean += static_cast<double>(i)*pd_[i];
    mean += 0.5*static_cast<double>(n-1)*pd_[n-1];
  }
  //Add on step sizes for each integral,
  // which is both area and step size in x,y
  //dflux twice -- once for the conversion from i to flux, once for
  // the step size
  mean *= dflux*dflux;

  if (donorm) mean /= getIntegral();

  mean += minflux;
}

/////////////////// STOP ////////////////////

/*!
  \param[out] mean mean value
  \param[out] var  variance
  \param[in] donorm Do not assume that P(D) is normalized.
 */
void PD::getMeanAndVar(double& mean, double& var,
		       bool donorm) const {
  if (n == 0) {
    mean = std::numeric_limits<double>::quiet_NaN();
    var = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  
  double normfac = 1.0;
  if (donorm) normfac = 1.0/getIntegral();
  
  //We do this as a <x^2> - <x>^2 calculation
  //Why not just call getMeans?  To avoid calling getIntegral
  // twice.  After this, mean1 and mean2 will be the actual
  // means/dflux - minflux.
  var = mean = 0.0;
  if (logflat) {
    //i=0 has weight 0
    for (unsigned int i = 1; i < n-1; ++i)
      mean += i*exp2(pd_[i]);
    mean += 0.5*(n-1)*exp2(pd_[n-1]);
    for (unsigned int i = 1; i < n-1; ++i)
      var += i*i*exp2(pd_[i]);
    var += 0.5*(n-1)*(n-1)*exp2(pd_[n-1]);
  } else {
    for (unsigned int i = 1; i < n-1; ++i)
      mean += i*pd_[i];
    mean += 0.5*(n-1)*pd_[n-1];
    for (unsigned int i = 1; i < n-1; ++i)
      var += i*i*pd_[i];
    var += 0.5*(n-1)*(n-1)*pd_[n-1];
  }

  //Multiply flux bits in
  mean *= dflux*dflux;
  var  *= dflux*dflux*dflux;

  //Apply normalization
  if (donorm) {
    mean *= normfac;
    var  *= normfac;
  }
  
  //Form real var
  var = var - mean*mean;
  //Correct mean for offset
  mean += minflux;
}

PD& PD::operator=(const PD& other) {
  if ( this == &other ) return *this; //Self-copy
  resize(other.n);
  minflux = other.minflux;
  dflux   = other.dflux;
  if (n > 0)
    for (unsigned int i = 0; i < n; ++i)
      pd_[i] = other.pd_[i];
  logflat = other.logflat;
  return *this;
}

void PD::fill(unsigned int N, FFTFLOAT MINFLUX, FFTFLOAT DFLUX,
	      const FFTFLOAT* const PD, bool LOG) {
  logflat = LOG;
  resize(N);
  minflux = MINFLUX;
  dflux = DFLUX;
  if (n > 0)
    for (unsigned int i = 0; i < n; ++i)
      pd_[i] = PD[i];
}

FFTFLOAT PD::getPDVal(double x) const {
  if (pd_ == NULL) return std::numeric_limits<FFTFLOAT>::quiet_NaN();

  //look up the effective indexes
  FFTFLOAT fx = static_cast<FFTFLOAT>(x);
  int idx = static_cast<int>((fx - minflux) / dflux);

  FFTFLOAT maxflux = minflux + static_cast<FFTFLOAT>(n-1) * dflux;

  if (fx < minflux) return pd_[0];
  if (fx > (maxflux - dflux)) return pd_[n-1]; 
  FFTFLOAT p0 = pd_[idx];
  FFTFLOAT dx = fx - (minflux + static_cast<FFTFLOAT>(idx) * dflux);
  return p0 + dx * (pd_[idx+1] - p0) / dflux;
}


std::ostream& PD::writeToStream(std::ostream& os) const {
  os << n << " " << minflux << " " << dflux << std::endl;
  if (n > 0) {
    os << pd_[0];
    for (unsigned int i = 1; i < n-1; ++i)
      os << " " << pd_[i];
    os << std::endl;
  }
  return os;
}

/*!
  \param[in] outputfile File to write to
  \returns 0 on success, an error code (!=0) for anything else
 */
int PD::writeToFits(const std::string& outputfile) const {

  //Make the fits file
  int status = 0;
  fitsfile *fp;

  fits_create_file(&fp, outputfile.c_str(), &status);

  if (status) {
    fits_report_error(stderr,status);
    throw pofdExcept("PD","writeToFits",
		     "Error creating FITS output file",1);
  }

  long axissize[1];
  axissize[0] = static_cast<long>(n);
  fits_create_img(fp, FITSIOIMG, 1, axissize, &status);
  
  //Add "WCS" info to hdr
  float crpix = 1;
  FFTFLOAT tmpval;
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE1"),
		 const_cast<char*>("FLUX"),
		 const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		 const_cast<char*>("Ref pix of axis 1"), &status );
  tmpval = minflux;
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CRVAL1"), &tmpval, 
		 const_cast<char*>("val at ref pix"), &status );
  tmpval = dflux;
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CDELT1"), &tmpval,
		 const_cast<char*>("delta along axis 1"), &status );

  int lg = static_cast<int>(logflat);
  fits_write_key(fp, TLOGICAL, const_cast<char*>("LOG"),&lg,
		 const_cast<char*>("Is log P(D) stored?"), &status);

  //Do data writing.  
  long fpixel[1] = { 1 };
  fits_write_pix(fp, TFITSFLOAT, fpixel, n, pd_, &status);

  fits_close_file(fp, &status);

  if (status) {
    fits_report_error(stderr,status);
    throw pofdExcept("PD","writeToFits",
		     "Error doing FITS write",2);
  }
  return status;
}

double PD::getLogLike(const simImage& data) const {
  if (pd_ == NULL) throw pofdExcept("PD","getLogLike",
                                    "pd not filled before likelihood calc",1);
  unsigned int ndata = data.getN1()*data.getN2();
  if (ndata == 0) throw pofdExcept("PD","getLogLike",
				   "No data present",2);

  if (data.isBinned()) return getLogLikeBinned(data);
  else return getLogLikeUnbinned(data);
}

double PD::getLogLikeBinned(const simImage& data) const {

  if (!data.isBinned())
    throw pofdExcept("PD","getLogLikeBinned",
		     "Data is not binned",1);

  //Quantities for edge test
  double maxflux = static_cast<double>(minflux + static_cast<double>(n-1) * 
				       dflux);

  int idx; //!< Index look up
  const unsigned int* bins;
  unsigned int nbins, ninbin;
  double cflux, bincent0, bindelta, loglike;
  double t, delt, maxfluxval;
  loglike = 0.0;
  bins = data.getBinnedData();
  bincent0 = data.getBinCent0();
  bindelta = data.getBinDelta();
  nbins = data.getNBins();
  double idflux = 1.0 / static_cast<double>(dflux);
  maxfluxval = maxflux - static_cast<double>(dflux);

  if (logflat) {
    for (unsigned int i = 0; i < nbins; ++i) {
      ninbin = bins[i];
      if (ninbin == 0) continue;
      cflux = bincent0 + static_cast<double>(i) * bindelta;
      //Get effective index
      if (cflux < minflux) 
	loglike += static_cast<double>(ninbin) * static_cast<double>(pd_[0]);
      else if (cflux > maxfluxval) 
	loglike += static_cast<double>(ninbin) * static_cast<double>(pd_[n-1]);
      else {
        //Not off edge
	delt = (cflux - minflux) * idflux;
	idx = static_cast<int>(delt);
	t   = delt - static_cast<double>(idx);
	loglike += ((1.0-t)*static_cast<double>(pd_[idx]) + 
		    t * static_cast<double>(pd_[idx+1])) *
	  static_cast<double>(ninbin);
      }
    }
  } else {
    //Not stored as log -- inefficient, but supported
    //Note that it would be insane to do this multiplicatively,
    // then take the log.  Also, it's better to interpolate
    // in log space than interpolate, then log
    for (unsigned int i = 0; i < nbins; ++i) {
      ninbin = bins[i];
      if (ninbin == 0) continue;
      cflux = bincent0 + static_cast<double>(i) * bindelta;
      delt = (cflux - minflux) * idflux;
      idx = static_cast<int>(delt);
      if (cflux < minflux) 
	loglike += static_cast<double>(ninbin) * log2(pd_[0]);
      else if ( cflux > maxfluxval ) 
	loglike += static_cast<double>(ninbin) * log2(pd_[n-1]);
      else {
	t   = delt - static_cast<double>(idx);
	loglike += ((1.0-t)*log2(pd_[idx]) + t*log2(pd_[idx+1])) *
	  static_cast<double>(ninbin);
      }
    }
  }
  //This has been base 2 -- convert back to base e
  return pofd_coverage::log2toe * loglike;
}

double PD::getLogLikeUnbinned(const simImage& data) const {

  unsigned int ndata = data.getN1()*data.getN2();

  //Quantities for edge test
  double minfluxi = static_cast<double>(minflux);
  double dfluxi = static_cast<double>(dflux);
  double maxflux = minfluxi + static_cast<double>(n-1) * dfluxi;

  int idx; //!< Index look up
  const double* flux;
  double cflux, loglike, t, delt, maxfluxval;
  loglike = 0.0;
  flux = data.getData();
  double idflux = 1.0 / dfluxi;
  maxfluxval = maxflux - dfluxi;

  if (logflat) {
    //Stored as log2 P(D)
    for (unsigned int i = 0; i < ndata; ++i) {
      cflux = flux[i];
      //Get effective index
      if (cflux <= minfluxi) loglike += static_cast<double>(pd_[0]);
      else if (cflux >= maxfluxval) 
	loglike += static_cast<double>(pd_[n-1]);
      else {
        //Not off edge
	delt = (cflux - minfluxi) * idflux;
	idx = static_cast<int>(delt);
	t   = delt - static_cast<double>(idx);
	loglike += (1.0-t)*pd_[idx] + t*pd_[idx+1];
      }
    }
  } else {
    //Not stored as log -- inefficient, but supported
    //Note that it would be insane to do this multiplicatively,
    // then take the log.  Also, it's better to interpolate
    // in log space than interpolate, then log
    for (unsigned int i = 0; i < ndata; ++i) {
      cflux = flux[i]; 
      if (cflux <= minfluxi) loglike += log2(pd_[0]);
      else if (cflux >= maxfluxval) loglike += log2(pd_[n-1]);
      else {
	delt = (cflux - minfluxi) * idflux;
	idx = static_cast<int>(delt);
	t   = delt - static_cast<double>(idx);
	loglike += (1.0-t)*log2(pd_[idx]) + t*log2(pd_[idx+1]);
      }
    }
  }
  //This has been base 2 -- convert back to base e
  return pofd_coverage::log2toe * loglike;
}


std::ostream& operator<<(std::ostream& os, const PD& b) {
  b.writeToStream(os);
  return os;
}
