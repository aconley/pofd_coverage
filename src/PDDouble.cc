//PDDouble.cc
#include<limits>

#include<fitsio.h>

#include "../include/PDDouble.h"
#include "../include/pofdExcept.h"

const double PDDouble::lowsigval = 3.0;

PDDouble::PDDouble(unsigned int N1, FFTFLOAT MINFLUX1, FFTFLOAT DFLUX1,
		   unsigned int N2, FFTFLOAT MINFLUX2, FFTFLOAT DFLUX2,
		   bool LOG) : n1(N1), n2(N2), capacity(N1*N2),
			       logflat(LOG),
			       minflux1(MINFLUX1), dflux1(DFLUX1),
			       minflux2(MINFLUX2), dflux2(DFLUX2) {
  if (capacity == 0) pd_ = NULL; else
    pd_ = (FFTFLOAT *) FFTWMALLOC(sizeof(FFTFLOAT)*capacity);
}

PDDouble::~PDDouble() {
  if (pd_ != NULL) FFTWFREE(pd_);
}

/*!
  Generally doesn't preserve data
 */
void PDDouble::resize(unsigned int N1, unsigned int N2) {
  //Doesn't actually resize arrays if it can avoid it
  unsigned int newcap = N1*N2;
  if ( newcap > capacity ) {
    if (pd_ != NULL) FFTWFREE(pd_);
    if (newcap > 0) pd_ = (FFTFLOAT *) FFTWMALLOC( sizeof(FFTFLOAT)*newcap );
    else pd_ = NULL;
    capacity = newcap;
  }
  n1 = N1;
  n2 = N2;
}

/*!
  Tries to preserve data
 */
void PDDouble::shrink() {
  unsigned int newcap = n1*n2;
  if ( newcap < capacity ) {
    if (newcap > 0) {
      FFTFLOAT* tmp = (FFTFLOAT*) FFTWMALLOC( sizeof(FFTFLOAT)*newcap );
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
void PDDouble::strict_resize(unsigned int N1, unsigned int N2) {
  unsigned int newcap = N1*N2;
  if (newcap != capacity) {
    if (pd_ != NULL) FFTWFREE(pd_);
    if (newcap > 0) pd_ = (FFTFLOAT*) FFTWMALLOC( sizeof(FFTFLOAT)*newcap );
    else pd_ = NULL;
    capacity = newcap;
  }
  n1 = N1;
  n2 = N2;
}

double PDDouble::getTotal() const {
  if ( (n1 == 0) || (n2 == 0) )
    return std::numeric_limits<double>::quiet_NaN();
  double retval;
  unsigned int sz = n1*n2;
  if (logflat) {
    retval = exp2(pd_[0]);
    for (unsigned int i = 1; i < sz; ++i)
      retval += exp2(pd_[i]);
  } else {
    retval = pd_[0];
    for (unsigned int i = 1; i < sz; ++i)
      retval += pd_[i];
  }
  return retval;
}

double PDDouble::getIntegral() const {
  if ( (n1 == 0) || (n2 == 0) )
    return std::numeric_limits<double>::quiet_NaN();
  
  double tot;
  FFTFLOAT *rowptr;
  if (logflat) {
    tot = 0.5 * exp2(pd_[0]);
    for (unsigned int j = 1; j < n2-1; ++j)
      tot += exp2(pd_[j]);
    tot += 0.5 * exp2(pd_[n2-1]);
    tot *= 0.5;
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + i * n2;
      tot += 0.5 * exp2(rowptr[0]);
      for (unsigned int j = 1; j < n2-1; ++j)
	tot += exp2(rowptr[j]);
      tot += 0.5 * exp2(rowptr[n2-1]);
    }
    rowptr = pd_ + (n1-1)*n2;
    tot += 0.25 * exp2(rowptr[0]);
    for (unsigned int j = 1; j < n2-1; ++j)
      tot += 0.5 * exp2(rowptr[j]);
    tot += 0.25 * exp2(rowptr[n2-1]);
  } else {
    tot = 0.5 * pd_[0];
    for (unsigned int j = 1; j < n2-1; ++j)
      tot += pd_[j];
    tot += 0.5 * pd_[n2-1];
    tot *= 0.5;
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + i*n2;
      tot += 0.5 * rowptr[0];
      for (unsigned int j = 1; j < n2-1; ++j)
	tot += rowptr[j];
      tot += 0.5 * rowptr[n2-1];
    }
    rowptr = pd_ + (n1-1)*n2;
    tot += 0.25 * rowptr[0];
    for (unsigned int j = 1; j < n2-1; ++j)
      tot += 0.5 * rowptr[j];
    tot += 0.25 * rowptr[n2-1];
  }
  return tot * dflux1 * dflux2;
}

/*!
  Normalize the P(D), using the trapezoidal rule
  to integrate
 */
void PDDouble::normalize() {
  if ((n1 == 0) || (n2 == 0))
    throw pofdExcept("PDDouble", "normalize",
		     "No information present to normalize", 1);
  //Note, because of the 0.5 edge pieces we don't just use
  // getTotal
  double tot = getIntegral();
  unsigned int sz = n1 * n2;
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

void PDDouble::applyLog(bool nocheck) {
  if (logflat) return;
  unsigned int sz = n1*n2;
  if (nocheck)
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] = log2(pd_[i]);
  else {
    FFTFLOAT val;
    for (unsigned int i = 0; i < sz; ++i) {
      val = pd_[i];
      if (val <= 0.0) pd_[i] = pofd_coverage::smalllogval; 
      else pd_[i] = log2(val);
    }
  }
  logflat = true;
}


void PDDouble::deLog() {
  if (!logflat) return;
  unsigned int sz = n1 * n2;
  for (unsigned int i = 0; i < sz; ++i)
    pd_[i] = exp2(pd_[i]);
  logflat = false;
}


/*
  \param[in] donorm Do not assume P(D) is normalized
 */
void PDDouble::edgeFix(bool donorm) {
  //Compute mean and stdev in each row and column.
  if (n1 < 3 || n2 < 3) return; //No point

  if (logflat)
    throw pofdExcept("PDDouble", "edgeFix",
		     "Not supported for logged PDs", 1);

  //Get mean and vars
  double mn1, mn2, var1, var2;
  getMeansAndVars(mn1,mn2,var1,var2,donorm);
  if (std::isnan(mn1) || std::isinf(mn1) ||
      std::isnan(var1) || std::isinf(var1) ||
      std::isnan(mn2) || std::isinf(mn2) ||
      std::isnan(var2) || std::isinf(var2))
    throw pofdExcept("PDDouble", "edgeFix",
		     "Problem with means/vars", 2);
  
  double istdev1 = 1.0 / sqrt(var1);
  double istdev2 = 1.0 / sqrt(var2);

  //Figure out what indexes these represent in x and y
  double maxfluxfix1, maxfluxfix2, ddflux1, ddflux2, dminflux1, dminflux2;
  int maxidx1, maxidx2;
  ddflux1 = static_cast<double>(dflux1);
  dminflux1 = static_cast<double>(minflux1);
  ddflux2 = static_cast<double>(dflux2);
  dminflux2 = static_cast<double>(minflux2);
  maxfluxfix1 = mn1 - PDDouble::lowsigval*sqrt(var1);
  maxfluxfix2 = mn2 - PDDouble::lowsigval*sqrt(var2);
  maxidx1 = static_cast<int>((maxfluxfix1 - dminflux1) / ddflux1);
  maxidx2 = static_cast<int>((maxfluxfix2 - dminflux2) / ddflux2);
  maxfluxfix1 = dminflux1 + maxidx1 * ddflux1;
  maxfluxfix2 = dminflux2 + maxidx2 * ddflux2;
  
  //Do edges now
  FFTFLOAT *rowptr;
  double tval, stepfac, subfac, pdval, preconst;
  if (maxidx1 > 1) {
    int minidx2 = (maxidx2 > 0) ? maxidx2 : 0;
    for (unsigned int j =  minidx2; j < n2; ++j) {
      pdval = static_cast<double>(pd_[maxidx1*n2+j]);
      tval = (maxfluxfix1 - mn1) * istdev1;
      preconst = pdval * exp(0.5 * tval * tval);
      subfac = (dminflux1 - mn1) * istdev1;
      stepfac = ddflux1 * istdev1;
      for (int i = 0; i < maxidx1; ++i) {
	tval = subfac + i * stepfac;
	pd_[i*n2+j] = static_cast<FFTFLOAT>(preconst * exp(-0.5 * tval * tval));
      }
    }
  }
  if (maxidx2 > 1) {
    int minidx1 = (maxidx1 > 0) ? maxidx1 : 0;
    for (unsigned int i =  minidx1; i < n1; ++i) {
      rowptr = pd_ + i*n2;
      pdval = static_cast<double>(rowptr[maxidx2]);
      tval = (maxfluxfix2 - mn2) * istdev2;
      preconst = pdval * exp(0.5 * tval * tval);
      subfac = (dminflux2 - mn2) * istdev2;
      stepfac = ddflux2 * istdev2;
      for (int j = 0; j < maxidx2; ++j) {
	tval = subfac + j * stepfac;
	rowptr[j] = static_cast<FFTFLOAT>(preconst * exp(-0.5 * tval * tval));
      }
    }
  }

  //Corner is tricky.
  // We will extrapolate in from both sides (if available)
  // and take the geometric mean
  if (maxidx1 > 0 && maxidx2 > 0) {
    for (int i = 0; i < maxidx1; ++i) {
      rowptr = pd_ + i*n2;
      pdval = static_cast<double>(rowptr[maxidx2]);
      tval = (maxfluxfix2 - mn2) * istdev2;
      preconst = pdval * exp(0.5 * tval * tval);
      subfac = (dminflux2 - mn2) * istdev2;
      stepfac = ddflux2 * istdev2;
      for (int j = 0; j < maxidx2; ++j) {
	tval = subfac + j * stepfac;
	rowptr[j] = static_cast<FFTFLOAT>(preconst * exp(-0.5 * tval * tval));
      }
    }
    double pdval2;
    for (int j =  0; j < maxidx2; ++j) {
      pdval = static_cast<double>(pd_[maxidx1*n2+j]);
      tval = (maxfluxfix1 - mn1) * istdev1;
      preconst = pdval * exp(0.5 * tval * tval);
      subfac = (dminflux1 - mn1) * istdev1;
      stepfac = ddflux1 * istdev1;
      for (int i = 0; i < maxidx1; ++i) {
	tval = subfac + i * stepfac;
	pdval2 = static_cast<double>(pd_[i*n2+j]);
	pd_[i*n2+j] = static_cast<FFTFLOAT>(sqrt(pdval2 * preconst * 
						 exp(-0.5*tval*tval)));
      }
    }
  }
}


/*!
  \param[out] mean1 mean along axis 1
  \param[out] mean2 mean along axis 2
  \param[in] donorm Do not assume that P(D) is normalized.
 */
void PDDouble::getMeans(double& mean1, double& mean2,
			bool donorm) const {
  if ( (n1 == 0) || (n2 == 0) ) {
    mean1 = std::numeric_limits<double>::quiet_NaN();
    mean2 = std::numeric_limits<double>::quiet_NaN();
    return;
  }


  //We use the trapezoidal rule here for the integrals
  // so it isn't quite a simple sum
  FFTFLOAT *rowptr, pval;
  double xsum;
  mean1 = mean2 = 0.0;

  if (!logflat) {
    //Integrate over lowest x value.  Note that x = 0
    // here (modulo the minflux, which we add later anyways)
    // so there is no mean1 contribution
    for (unsigned int j = 1; j < n2-1; ++j)
      mean2 += j*pd_[j]; 
    mean2 += 0.5*(n2-1)*pd_[n2-1]; //(y=0 at pd_[0])
    mean2 *= 0.5; //End bit of trap in x, so 1/2 factor

    //Now main body of trap
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + n2*i;
      xsum = 0.5 * rowptr[0]; //xsum will be the row sum, mult by x later
      //mean2 += 0.5*0*rowptr[0] obviously not needed
      for (unsigned int j = 1; j < n2-1; ++j) {
	pval = rowptr[j];
	xsum += pval;
	mean2 += static_cast<double>(j) * pval;
      }      
      xsum += 0.5 * rowptr[n2-1];
      mean1 += xsum * static_cast<double>(i); //Multiply in x value
      mean2 += 0.5 * (n2-1) * rowptr[n2-1];
    }

    //Endpiece, all multiplied by 1/2 since last x bit
    rowptr = pd_ + (n1-1)*n2;
    xsum = 0.5 * rowptr[0];
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = rowptr[j];
      xsum += pval;
      mean2 += 0.5 * static_cast<double>(j) * pval;
    }
    xsum += 0.5 * rowptr[n2-1];
    mean1 += 0.5 * xsum * (n1 - 1);
    mean2 += 0.25 * (n2 - 1) * rowptr[n2-1];
  } else {
    for (unsigned int j = 1; j < n2-1; ++j)
      mean2 += j * exp2(pd_[j]);
    mean2 += 0.5 * (n2 - 1) * exp2(pd_[n2-1]);
    mean2 *= 0.5; 
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + n2*i;
      xsum = 0.5 * exp2(rowptr[0]); 
      for (unsigned int j = 1; j < n2-1; ++j) {
	pval = exp2(rowptr[j]);
	xsum += pval;
	mean2 += static_cast<double>(j) * pval;
      }      
      pval = exp2(rowptr[n2-1]);
      xsum += 0.5 * pval;
      mean1 += xsum * static_cast<double>(i);
      mean2 += 0.5 * (n2 - 1) * pval;
    }
    rowptr = pd_ + (n1 - 1) * n2;
    xsum = 0.5*exp2(rowptr[0]);
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = exp2(rowptr[j]);
      xsum += pval;
      mean2 += 0.5 * static_cast<double>(j) * pval;
    }
    pval = exp2(rowptr[n2 - 1]);
    xsum += 0.5 * pval;
    mean1 += 0.5 * xsum * (n1-1);
    mean2 += 0.25 * (n2 - 1) * pval;
  }

  //Add on step sizes for each integral,
  // which is both area and step size in x,y
  mean1 *= dflux1 * dflux1 * dflux2;
  mean2 *= dflux1 * dflux2 * dflux2;

  if (donorm) {
    double inorm = 1.0 / getIntegral();
    mean1 *= inorm;
    mean2 *= inorm;
  }
  mean1 += minflux1;
  mean2 += minflux2;
}

/*!
  \param[out] mean1 mean along axis 1
  \param[out] mean2 mean along axis 2
  \param[out] var1  variance along axis 1
  \param[out] var2  variance along axis 2
  \param[in] donorm Do not assume that P(D) is normalized.
 */
void PDDouble::getMeansAndVars(double& mean1, double& mean2,
			       double& var1, double& var2,
			       bool donorm) const {
  if ( (n1 == 0) || (n2 == 0) ) {
    mean1 = std::numeric_limits<double>::quiet_NaN();
    mean2 = std::numeric_limits<double>::quiet_NaN();
    var1  = std::numeric_limits<double>::quiet_NaN();
    var2  = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  double normfac = 1.0;
  if (donorm) normfac = 1.0/getIntegral();

  //First compute means
  //Why not just call getMeans?  To avoid calling getIntegral
  // twice.  After this, mean1 and mean2 will be the actual
  // means/dflux - minflux.
  FFTFLOAT pval, *rowptr;
  double xsum;
  mean1 = mean2 = 0.0;
  if (!logflat) {
    for (unsigned int j = 1; j < n2-1; ++j)
      mean2 += j*pd_[j]; 
    mean2 += 0.5*(n2-1)*pd_[n2-1]; 
    mean2 *= 0.5; 
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + n2*i;
      xsum = 0.5*rowptr[0];
      for (unsigned int j = 1; j < n2-1; ++j) {
	pval = rowptr[j];
	xsum += pval;
	mean2 += static_cast<double>(j) * pval;
      }      
      xsum += 0.5*rowptr[n2-1];
      mean1 += xsum * static_cast<double>(i); 
      mean2 += 0.5*(n2-1)*rowptr[n2-1];
    }
    rowptr = pd_ + (n1-1)*n2;
    xsum = 0.5*rowptr[0];
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = rowptr[j];
      xsum += pval;
      mean2 += 0.5*static_cast<double>(j) * pval;
    }
    xsum += 0.5*rowptr[n2-1];
    mean1 += 0.5*xsum*(n1-1);
    mean2 += 0.25*(n2-1)*rowptr[n2-1];
  } else {
    for (unsigned int j = 1; j < n2-1; ++j)
      mean2 += j*exp2(pd_[j]);
    mean2 += 0.5*(n2-1)*exp2(pd_[n2-1]);
    mean2 *= 0.5; 
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + n2*i;
      xsum = 0.5*exp2(rowptr[0]); 
      for (unsigned int j = 1; j < n2-1; ++j) {
	pval = exp2(rowptr[j]);
	xsum += pval;
	mean2 += static_cast<double>(j) * pval;
      }      
      pval = exp2(rowptr[n2-1]);
      xsum += 0.5*pval;
      mean1 += xsum * static_cast<double>(i);
      mean2 += 0.5*(n2-1)*pval;
    }
    rowptr = pd_ + (n1-1)*n2;
    xsum = 0.5*exp2(rowptr[0]);
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = exp2(rowptr[j]);
      xsum += pval;
      mean2 += 0.5*static_cast<double>(j) * pval;
    }
    pval = exp2(rowptr[n2-1]);
    xsum += 0.5*pval;
    mean1 += 0.5*xsum*(n1-1);
    mean2 += 0.25*(n2-1)*pval;
  }
  mean1 *= dflux1*dflux2;
  mean2 *= dflux1*dflux2;

  if (donorm) {
    mean1 *= normfac;
    mean2 *= normfac;
  }
  
  //Now variances, pretty much the same calculation
  // Recall mean1, mean2 are means/dflux - minflux
  var1 = var2 = 0.0;
  double deltax, deltay;
  if (!logflat) {
    //i=0 bit, 1/2 factor for end
    pval = pd_[0];
    xsum = 0.5*pval;
    var2 = -0.5*mean2*mean2*pval; // deltay = -mean2
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = pd_[j];
      xsum += pval;
      deltay = j - mean2;
      var2 += pval*deltay*deltay;
    }
    pval = pd_[n2-1];
    xsum += 0.5*pval;
    //1/2 factor for end bit, x-<mean> was -mean1 here
    var1 = 0.5 * mean1 * mean1 * xsum; 
    deltay = (n2-1) - mean2;
    var2 += 0.5*deltay*deltay*pval;
    var2 *= 0.5; //1/2 factor for end
    
    //Now core bit
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + i*n2;
      pval = rowptr[0];
      deltax = i - mean1;
      xsum = 0.5*pval;
      var2 += 0.5*mean2*mean2*pval; // deltay = -mean2
      for (unsigned int j = 1; j < n2-1; ++j) {
	pval = rowptr[j];
	xsum += pval;
	deltay = j - mean2;
	var2 += deltay*deltay*pval;
      }
      pval = rowptr[n2-1];
      xsum += 0.5*pval;
      var1 += xsum*deltax*deltax;
      deltay = n2-1-mean2;
      var2 += 0.5*deltay*deltay*pval;
    }

    //And now the top end in x
    rowptr = pd_ + (n1-1)*n2;
    pval = rowptr[0];
    deltax = n1-1-mean1;
    xsum = 0.5*pval;
    var2 += 0.25*mean2*mean2*pval;
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = rowptr[j];
      xsum += pval;
      deltay = j - mean2;
      var2 += 0.5*deltay*deltay*pval;
    }
    pval = rowptr[n2-1];
    xsum += 0.5*pval;
    var1 += 0.5*xsum*deltax*deltax; //1/2 for end
    deltay = n2-1-mean2;
    var2 += 0.25*deltay*deltay*pval;  //1/4 is 1/2*1/2 for double end
  } else {
    pval = exp2(pd_[0]);
    xsum = 0.5*pval;
    var2 = -0.5*mean2*mean2*pval; // deltay = -mean2
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = exp2(pd_[j]);
      xsum += pval;
      deltay = j - mean2;
      var2 += pval*deltay*deltay;
    }
    pval = exp2(pd_[n2-1]);
    xsum += 0.5*pval;
    var1 = 0.5 * mean1 * mean1 * xsum; 
    deltay = (n2-1) - mean2;
    var2 += 0.5*deltay*deltay*pval;
    var2 *= 0.5;
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + i*n2;
      pval = exp2(rowptr[0]);
      deltax = i - mean1;
      xsum = 0.5*pval;
      var2 += 0.5*mean2*mean2*pval; // deltay = -mean2
      for (unsigned int j = 1; j < n2-1; ++j) {
	pval = exp2(rowptr[j]);
	xsum += pval;
	deltay = j - mean2;
	var2 += deltay*deltay*pval;
      }
      pval = exp2(rowptr[n2-1]);
      xsum += 0.5*pval;
      var1 += xsum*deltax*deltax;
      deltay = n2-1-mean2;
      var2 += 0.5*deltay*deltay*pval;
    }
    rowptr = pd_ + (n1-1)*n2;
    pval = exp2(rowptr[0]);
    deltax = n1-1-mean1;
    xsum = 0.5*pval;
    var2 += 0.25*mean2*mean2*pval;
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = exp2(rowptr[j]);
      xsum += pval;
      deltay = j - mean2;
      var2 += 0.5*deltay*deltay*pval;
    }
    pval = exp2(rowptr[n2-1]);
    xsum += 0.5*pval;
    var1 += 0.5*xsum*deltax*deltax; //1/2 for end
    deltay = n2-1-mean2;
    var2 += 0.25*deltay*deltay*pval;  //1/4 is 1/2*1/2 for double end
  }

  //Integral dx dy
  var1 *= dflux1*dflux2;
  var2 *= dflux1*dflux2;

  if (donorm) {
    var1 *= normfac;
    var2 *= normfac;
  }

  //Put in step size in x/y bit
  mean1 *= dflux1;
  mean2 *= dflux2;
  var1  *= dflux1*dflux1;
  var2  *= dflux2*dflux2;

  mean1 += minflux1;
  mean2 += minflux2;
}



PDDouble& PDDouble::operator=(const PDDouble& other) {
  if ( this == &other ) return *this; //Self-copy
  resize(other.n1,other.n2);
  minflux1 = other.minflux1;
  dflux1 = other.dflux1;
  minflux2 = other.minflux2;
  dflux2 = other.dflux2;
  unsigned int sz = n1 * n2;
  if (sz > 0) {
    if (other.pd_ == NULL)
      throw pofdExcept("PDDouble", "operator=", 
			"Copying from uninitialized other", 1);
    if (pd_ == NULL)
      throw pofdExcept("PDDouble", "operator=", 
			"initialization of this space failed", 2);
      
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] = other.pd_[i];
  }
  logflat = other.logflat;
  return *this;
}

void PDDouble::fill(unsigned int N1, FFTFLOAT MINFLUX1, FFTFLOAT DFLUX1,
		    unsigned int N2, FFTFLOAT MINFLUX2, FFTFLOAT DFLUX2,
		    const FFTFLOAT* const PD, bool LOG) {

  logflat = LOG;
  resize(N1,N2);
  minflux1 = MINFLUX1;
  dflux1 = DFLUX1;
  minflux2 = MINFLUX2;
  dflux2 = DFLUX2;
  unsigned int sz = n1*n2;
  if (sz > 0) {
    if (PD == NULL)
      throw pofdExcept("PDDouble", "fill", 
			"PD is not initialized", 1);
    if (pd_ == NULL)
      throw pofdExcept("PDDouble", "fill", 
			"Initialization of this space failed", 2);

    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] = PD[i];
  }
}

/*!
  This doesn't interpolate the same way as getLogLike, returning
  the interpolation in whatever way the PD is stored.
 */
FFTFLOAT PDDouble::getPDVal(double x, double y, bool logval) const {
  if (pd_ == NULL) return std::numeric_limits<FFTFLOAT>::quiet_NaN();

  //look up the effective indexes
  FFTFLOAT fx = static_cast<FFTFLOAT>(x);
  FFTFLOAT fy = static_cast<FFTFLOAT>(y);
  int idx1 = static_cast<int>((fx - minflux1) / dflux1);
  int idx2 = static_cast<int>((fy - minflux2) / dflux2);
  int n2idx1 = n2 * idx1;

  FFTFLOAT maxflux1 = minflux1 + static_cast<FFTFLOAT>(n1-1) * dflux1;
  FFTFLOAT maxflux2 = minflux2 + static_cast<FFTFLOAT>(n2-1) * dflux2;

  unsigned int n2n1 = n2 * n1;
  FFTFLOAT interp_val;
  //Check to see if we are off the edge
  if (fx < minflux1) {
    if (fy < minflux2) interp_val = pd_[0];
    else if (fy > maxflux2) interp_val = pd_[n2 - 1];
    else interp_val = pd_[idx2];
  } else if (fx > maxflux1) {
    if (fy < minflux2) interp_val = pd_[n2n1 - n1];
    else if (fy > maxflux2) interp_val = pd_[n2n1 - 1];
    else interp_val =  pd_[n2n1-n1+idx2];
  } else if (fy < minflux2) {
    interp_val =  pd_[n2idx1];
  } else if (fy > maxflux2) {
    interp_val = pd_[n2idx1 + n2 - 1];
  } else {
    //Actual interpolation
    FFTFLOAT u, t, omu, omt;
    t = (x - minflux1) / dflux1 - static_cast<double>(idx1);
    u = (y - minflux2) / dflux2 - static_cast<double>(idx2);
    omu = 1.0 - u; omt = 1.0 - t;
    
    unsigned int baseidx = n2idx1+idx2;
    interp_val = omt*(omu * pd_[baseidx] + u * pd_[baseidx+1]) +
      t * (omu * pd_[baseidx+n2] + u * pd_[baseidx+n2+1]);
  }
  if (!logval) {
    if (logflat) return exp2(interp_val); else return interp_val;
  } else {
    //Note we return ln, not log2
    if (logflat) return pofd_coverage::log2toe*interp_val; 
    else return log(interp_val);
  }
}


std::ostream& PDDouble::writeToStream(std::ostream& os) const {
  os << n1 << " " << minflux1 << " " << dflux1 << std::endl;
  os << n2 << " " << minflux2 << " " << dflux2 << std::endl;
  FFTFLOAT *rowptr;
  if (n1 * n2 > 0) {
    for (unsigned int i = 0; i < n1; ++i) {
      rowptr = pd_ + i*n2;
      os << rowptr[0];
      for (unsigned int j = 1; j < n2; ++j)
	os << " " << rowptr[j];
      os << std::endl;
    }
  }
  return os;
}

/*!
  \param[in] outputfile File to write to
  \returns 0 on success, an error code (!=0) for anything else
 */
int PDDouble::writeToFits( const std::string& outputfile ) const {

  //Make the fits file
  int status = 0;
  fitsfile *fp;

  fits_create_file(&fp, outputfile.c_str(), &status);

  if (status) {
    fits_report_error(stderr, status);
    throw pofdExcept("PDDouble","writeToFits",
		       "Error creating FITS output file",1);
  }

  long axissize[2];
  axissize[0] = static_cast<long>(n1);
  axissize[1] = static_cast<long>(n2);
  
  fits_create_img(fp, FITSIOIMG, 2, axissize, &status);
  
  //Add "WCS" info to hdr
  float crpix = 1;
  FFTFLOAT tmpval;
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE1"),
		 const_cast<char*>("FLUX1"),
		 const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		 const_cast<char*>("Ref pix of axis 1"), &status);
  tmpval = minflux1;
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CRVAL1"), &tmpval, 
		 const_cast<char*>("val at ref pix"), &status);
  tmpval = dflux1;
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CDELT1"), &tmpval,
		 const_cast<char*>("delta along axis 1"), &status);

  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE2"),
		 const_cast<char*>("FLUX2"),
		 const_cast<char*>("Type of Data axis 2"), &status);
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX2"), &crpix, 
		 const_cast<char*>("Ref pix of axis 2"), &status);
  tmpval = minflux2;
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CRVAL2"), &tmpval, 
		 const_cast<char*>("val at ref pix"), &status);
  tmpval = dflux2;
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CDELT2"), &tmpval,
		 const_cast<char*>("delta along axis 2"), &status);
  tmpval = dflux1;
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CD1_1"), &tmpval,
		 const_cast<char*>("WCS matrix element 1 1"), &status);
  tmpval = 0.0;
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CD1_2"), &tmpval, 
		 const_cast<char*>("WCS matrix element 1 2"),
		 &status);
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CD2_1"), &tmpval, 
		 const_cast<char*>("WCS matrix element 2 1"),
		 &status);
  tmpval = dflux2;
  fits_write_key(fp, TFITSFLOAT, const_cast<char*>("CD2_2"), &tmpval, 
		 const_cast<char*>("WCS matrix element 2 2"),
		 &status);

  int lg = static_cast<int>(logflat);
  fits_write_key(fp, TLOGICAL, const_cast<char*>("LOG"),&lg,
		 const_cast<char*>("Is log P(D) stored?"), &status);

  //Do data writing.  We have to make a transposed copy of the
  // data to do this, which is irritating as hell
  FFTFLOAT *tmpdata = new FFTFLOAT[ n1 ];
  long fpixel[2] = { 1, 1 };
  for ( unsigned int j = 0; j < n2; ++j ) {
    for (unsigned int i = 0; i < n1; ++i) tmpdata[i] = pd_[i * n2 + j];
    fpixel[1] = static_cast<long>(j+1);
    fits_write_pix(fp, TFITSFLOAT, fpixel, n1, tmpdata, &status);
  }
  delete[] tmpdata;

  fits_close_file(fp,&status);

  if (status) {
    fits_report_error(stderr,status);
    throw pofdExcept("PDDouble","writeToFits",
		       "Error doing FITS write",2);
  }
  return status;
}

double PDDouble::getLogLike(const simImageDouble& data) const {
  if (pd_ == NULL) throw pofdExcept("PDDouble", "getLogLike",
				    "pd not filled before likelihood calc", 1);
  unsigned int ndata = data.getN1() * data.getN2();
  if (ndata == 0) throw pofdExcept("PDDouble", "getLogLike",
				   "No data present", 2);

  if (data.isBinned()) return getLogLikeBinned(data);
  else return getLogLikeUnbinned(data);
}

double PDDouble::getLogLikeUnbinned(const simImageDouble& data) const {
  unsigned int ndata = data.getN1() * data.getN2();

  //Quantities for edge test
  double dflux1i = static_cast<double>(dflux1); 
  double minflux1i = static_cast<double>(minflux1);
  double maxflux1 = minflux1i + static_cast<double>(n1-1)*dflux1i;
  double dflux2i = static_cast<double>(dflux2); 
  double minflux2i = static_cast<double>(minflux1); 
  double maxflux2 = minflux2i + static_cast<double>(n2-1)*dflux2i;

  int idx1, idx2, n2idx1; //!< Index look up
  unsigned int n2n1, baseidx;
  n2n1 = n2*n1;

  const double* flux1;
  const double* flux2;
  double cflux1, cflux2, loglike, delt1, delt2;
  FFTFLOAT interp_val;
  double u, t, omu, omt;
  double idflux1 = 1.0 / dflux1i;
  double idflux2 = 1.0 / dflux2i;

  loglike = 0.0;
  flux1 = data.getData1();
  flux2 = data.getData2();

  //Do interpolation.  Note we don't call getPDval, since
  // we do the interpolation here always in log space no matter
  // how it is stored internally, and because it's more efficient
  // to do it in house.
  if (logflat) {
    //Internal information is stored as log2 of P(D)
    for (unsigned int i = 0; i < ndata; ++i) {
      cflux1 = flux1[i]; cflux2 = flux2[i];
      //Get effective indices
      delt1 = (cflux1 - minflux1i) * idflux1;
      delt2 = (cflux2 - minflux2i) * idflux2;
      idx1 = static_cast<int>(delt1);
      idx2 = static_cast<int>(delt2);
      n2idx1 = n2 * idx1;
      if (cflux1 <= minflux1i) {
	if (cflux2 <= minflux2i) interp_val = pd_[0];
        else if (cflux2 >= maxflux2) interp_val = pd_[n2-1];
        else interp_val = pd_[idx2];
      } else if (cflux1 >= maxflux1) {
        if (cflux2 <= minflux2i) interp_val = pd_[n2n1-n1];
        else if (cflux2 >= maxflux2) interp_val = pd_[n2n1-1];
        else interp_val = pd_[n2n1-n1+idx2];
      } else if (cflux2 <= minflux2i) {
        interp_val = pd_[n2idx1];
      } else if (cflux2 >= maxflux2) {
        interp_val = pd_[n2idx1+n2-1];
      } else {
        //Not off edge
        t = delt1 - static_cast<double>(idx1);
        u = delt2 - static_cast<double>(idx2);
        omu = 1.0 - u; omt = 1.0 - t;
        baseidx = n2idx1 + idx2;
        interp_val = omt * (omu * pd_[baseidx] + u * pd_[baseidx+1]) +
          t * (omu * pd_[baseidx+n2] + u * pd_[baseidx+n2+1]);
      }
      loglike += interp_val;
    }
  } else {
    //Not stored as log2 -- inefficient, but supported
    //Note that it would be insane to do this multiplicatively,
    // then take the log.  Also, it's better to interpolate
    // in log space than interpolate, then log
    for (unsigned int i = 0; i < ndata; ++i) {
      cflux1 = flux1[i]; cflux2 = flux2[i];
      delt1 = (cflux1 - minflux1i) * idflux1;
      delt2 = (cflux2 - minflux2i) * idflux2;
      idx1 = static_cast<int>(delt1);
      idx2 = static_cast<int>(delt2);
      n2idx1 = n2 * idx1;
      
      if (cflux1 < minflux1i) {
        if (cflux2 < minflux2i) interp_val = log2(pd_[0]);
        else if (cflux2 > maxflux2) interp_val = log2(pd_[n2-1]);
        else interp_val = log2(pd_[idx2]);
      } else if (cflux1 > maxflux1) {
        if (cflux2 < minflux2i) interp_val = log2(pd_[n2n1-n1]);
        else if (cflux2 > maxflux2) interp_val = log2(pd_[n2n1-1]);
        else interp_val = log2(pd_[n2n1-n1+idx2]);
      } else if (cflux2 < minflux2i) {
        interp_val = log2(pd_[n2idx1]);
      } else if (cflux2 > maxflux2) {
        interp_val = log2(pd_[n2idx1+n2-1]);
      } else {
        //Not off edge
        t = delt1 - static_cast<double>(idx1);
        u = delt2 - static_cast<double>(idx2);
        omu = 1.0 - u; omt = 1.0-t;
        
        baseidx = n2idx1 + idx2;
        interp_val = omt * (omu*log2(pd_[baseidx]) + u*log2(pd_[baseidx+1])) +
          t * (omu*log2(pd_[baseidx+n2]) + u*log2(pd_[baseidx+n2+1]));
      }
      loglike += interp_val;
    }
  }
  //This has been base 2 -- convert back to base e
  return pofd_coverage::log2toe*loglike;
}


double PDDouble::getLogLikeBinned(const simImageDouble& data) const {
  if (!data.isBinned())
    throw pofdExcept("PDDouble","getLogLikeBinned",
		     "Data is not binned",1);

  //Quantities for edge test
  double dflux1i = static_cast<double>(dflux1); 
  double minflux1i = static_cast<double>(minflux1);
  double maxflux1i = minflux1i + static_cast<double>(n1-1)*dflux1i;
  double dflux2i = static_cast<double>(dflux2); 
  double minflux2i = static_cast<double>(minflux1); 
  double maxflux2i = minflux2i + static_cast<double>(n2-1)*dflux2i;

  int idx1, idx2, n2idx1; //!< Index look up
  unsigned int baseidx, n2n1;
  n2n1 = n2 * n1;

  const unsigned int *bins, *binptr;
  unsigned int ninbin, nbins;
  double cflux1, cflux2, bincent01, bincent02, bindelta1, bindelta2;
  double loglike, idflux1, idflux2;
  FFTFLOAT interp_val;
  double u, t, omu, omt;

  loglike = 0.0;
  nbins = data.getNBins();
  bincent01 = data.getBinCent01();
  bincent02 = data.getBinCent02();
  bindelta1 = data.getBinDelta1();
  bindelta2 = data.getBinDelta2();
  bins = data.getBinnedData();
  idflux1 = 1.0 / dflux1i;
  idflux2 = 1.0 / dflux2i;

  //Do interpolation.  Note we don't call getPDval, since
  // we do the interpolation here always in log space no matter
  // how it is stored internally, and because it's more efficient
  // to do it in house.
  if (logflat) {
    //Internal information is stored as log2 of P(D)
    for (unsigned int i = 0; i < nbins; ++i) {
      cflux1 = bincent01 + static_cast<double>(i) * bindelta1;
      idx1 = static_cast<int>((cflux1 - minflux1i) * idflux1);
      n2idx1 = n2 * idx1;
      binptr = bins + i*nbins;
      for (unsigned int j = 0; j < nbins; ++j) {
	ninbin = binptr[j];
	if (ninbin == 0) continue;  //skip calculation
	cflux2 = bincent02 + static_cast<double>(j)*bindelta2;
	idx2 = static_cast<int>((cflux2 - minflux2i) * idflux2);
	if (cflux1 < minflux1i) {
	  if (cflux2 < minflux2i) interp_val = pd_[0];
	  else if (cflux2 > maxflux2i) interp_val = pd_[n2-1];
	  else interp_val = pd_[idx2];
	} else if (cflux1 > maxflux1i) {
	  if (cflux2 < minflux2i) interp_val = pd_[n2n1-n1];
	  else if (cflux2 > maxflux2i) interp_val = pd_[n2n1-1];
	  else interp_val = pd_[n2n1-n1+idx2];
	} else if (cflux2 < minflux2i) {
	  interp_val = pd_[n2idx1];
	} else if (cflux2 > maxflux2i) {
	  interp_val = pd_[n2idx1+n2-1];
	} else {
	  //Not off edge
	  t = (cflux1 - minflux1i)*idflux1 - static_cast<double>(idx1);
	  u = (cflux2 - minflux2i)*idflux2 - static_cast<double>(idx2);
	  omu = 1.0 - u; omt = 1.0 - t;
	  baseidx = n2idx1 + idx2;
	  interp_val = omt*(omu*pd_[baseidx] + u*pd_[baseidx+1]) +
	    t*(omu*pd_[baseidx+n2] + u*pd_[baseidx+n2+1]);
	}
	loglike += static_cast<double>(ninbin)*interp_val;
      }
    }
  } else {
    //Not stored as log2 -- inefficient, but supported
    //Note that it would be insane to do this multiplicatively,
    // then take the log.  Also, it's better to interpolate
    // in log space than interpolate, then log
    for (unsigned int i = 0; i < nbins; ++i) {
      cflux1 = bincent01 + static_cast<double>(i)*bindelta1;
      idx1 = static_cast<int>((cflux1 - minflux1i) * idflux1);
      n2idx1 = n2 * idx1;
      binptr = bins + i * nbins;
      for (unsigned int j = 0; j < nbins; ++j) {
	ninbin = binptr[j];
	if (ninbin == 0) continue;  //skip calculation
	cflux2 = bincent02 + static_cast<double>(j) * bindelta2;
	idx2 = static_cast<int>((cflux2 - minflux2i) * idflux2);
	
	if (cflux1 < minflux1i) {
	  if (cflux2 < minflux2i) interp_val = log2(pd_[0]);
	  else if (cflux2 > maxflux2i) interp_val = log2(pd_[n2-1]);
	  else interp_val = log2(pd_[idx2]);
	} else if (cflux1 > maxflux1i) {
	  if (cflux2 < minflux2i) interp_val = log2(pd_[n2n1-n1]);
	  else if (cflux2 > maxflux2i) interp_val = log2(pd_[n2n1-1]);
	  else interp_val = log2(pd_[n2n1-n1+idx2]);
	} else if (cflux2 < minflux2i) {
	  interp_val = log2(pd_[n2idx1]);
	} else if (cflux2 > maxflux2i) {
	  interp_val = log2(pd_[n2idx1+n2-1]);
	} else {
	  //Not off edge
	  t = (cflux1 - minflux1i) * idflux1 - idx1;
	  u = (cflux2 - minflux2i) * idflux2 - idx2;
	  omu = 1.0 - u; omt = 1.0-t;
	  
	  baseidx = n2idx1 + idx2;
	  interp_val = omt*(omu*log2(pd_[baseidx]) + u*log2(pd_[baseidx+1])) +
	    t*(omu*log2(pd_[baseidx+n2]) + u*log2(pd_[baseidx+n2+1]));
	}
	loglike += static_cast<double>(ninbin)*interp_val;
      }
    }
  }
  //This has been base 2 -- convert back to base e
  return pofd_coverage::log2toe*loglike;
}


std::ostream& operator<<(std::ostream& os, const PDDouble& b) {
  b.writeToStream(os);
  return os;
}
