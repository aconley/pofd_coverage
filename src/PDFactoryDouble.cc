#include<sstream>
#include<cmath>
#include<cstring>
#include<limits>
#include<iomanip>

#include "hdf5.h"

#include "../include/global_settings.h"
#include "../include/PDFactoryDouble.h"
#include "../include/pofdExcept.h"

// Helper function for RFlux setup
// Returns dflux
double initRFluxInternal(unsigned int n, double minflux, double maxflux,
			 double* const rflux, double& minflux_realized) {
  // Assume n >= 1 (checked by initRFlux)
  double inm1 = 1.0 / static_cast<double>(n - 1);
  double dflux = (maxflux - minflux) * inm1;
  if (minflux >= 0.0) {
    // Easy case
    rflux[0] = minflux;
    for (unsigned int i = 1; i < n; ++i)
      rflux[i] = static_cast<double>(i) * dflux + minflux;
    minflux_realized = minflux;
  } else {
    // Here the complication is that we would really like to have
    // RFlux = 0 included in the array.  Also, we are wrapping 
    // negative fluxes around to the top of the array.
    // We do this by tweaking minflux slightly
    dflux = maxflux / (n - floor(-minflux / dflux) - 1.0);
    // Figure out what index we go up to with positive fills
    unsigned int maxpos = static_cast<unsigned int>(maxflux / dflux + 0.999999);
    rflux[0] = 0.0;
    for (unsigned int i = 1; i < maxpos + 1; ++i) // Pos Rflux
      rflux[i] = static_cast<double>(i) * dflux;

    // Figure out new minimum flux
    double wrapval = - static_cast<double>(n) * dflux;
    for (unsigned int i = maxpos + 1; i < n; ++i) // Wrapped neg Rflux
      rflux[i] = static_cast<double>(i) * dflux + wrapval;
    minflux_realized = rflux[maxpos + 1];
  }
  return dflux;
}

const double PDFactoryDouble::lowEdgeRMult=1e-9;
//Control of how we do the edge integrals -- linear or log?
const bool PDFactoryDouble::use_edge_log_x = true;
const bool PDFactoryDouble::use_edge_log_y = false;

/*!
  \param[in] nedge Size of edge integrals
 */
PDFactoryDouble::PDFactoryDouble(unsigned int nedge) {
  init(nedge);
}

/*
  \param[in] wisfile Wisdom file filename
  \param[in] nedge Size of edge integrals
 */
PDFactoryDouble::PDFactoryDouble(const std::string& wisfile,
				 unsigned int nedge) {
  init(nedge);
  addWisdom(wisfile);
}

PDFactoryDouble::~PDFactoryDouble() {
  if (RFlux1 != nullptr) fftw_free(RFlux1);
  if (RFlux2 != nullptr) fftw_free(RFlux2);

  if (rvals != nullptr) fftw_free(rvals);
  if (rsum != nullptr) fftw_free(rsum);
  if (rtrans != nullptr) fftw_free(rtrans);
  if (pofd != nullptr) fftw_free(pofd);
  if (pval != nullptr) fftw_free(pval);

  if (REdgeFlux1 != nullptr) fftw_free(REdgeFlux1);
  if (REdgeFlux2 != nullptr) fftw_free(REdgeFlux2);
  if (REdgeWork != nullptr) fftw_free(REdgeWork);

  if (plan != nullptr) fftw_destroy_plan(plan); 
  if (plan_inv != nullptr) fftw_destroy_plan(plan_inv);
}

/*!
  \param[in] NEDGE Number of edge integral steps
 */
void PDFactoryDouble::init(unsigned int NEDGE) {
  currsize = 0;

#ifdef TIMING
  resetTime();
#endif

  rvars_allocated = false;
  RFlux1 = RFlux2 = nullptr;
  minfluxR_1 = minfluxR_2 = 0.0;
  rinitialized = false;
  rdflux = false;
  rvals = nullptr;
  rsum = nullptr;
  rtrans = nullptr;
  pofd = nullptr;
  pval = nullptr;

  nedge = NEDGE;
  edgevars_allocated = false;
  REdgeFlux1 = nullptr;
  REdgeFlux2 = nullptr;
  REdgeWork  = nullptr;

  dflux1 = dflux2 = 0.0;

  plans_valid = false;
  plan = plan_inv = nullptr;

  verbose = false;
  has_wisdom = false;
  fftw_plan_style = FFTW_MEASURE;

  sigma1 = sigma2 = std::numeric_limits<double>::quiet_NaN();
  varnoi1 = varnoi2 = std::numeric_limits<double>::quiet_NaN();
  mn1 = mn2 = sg1 = sg2 = std::numeric_limits<double>::quiet_NaN();
  initialized = false;
}

/*!
  \param[in] nedg Edge integration size
 */
void PDFactoryDouble::setNEdge(unsigned int nedg) {
  if (nedg == nedge) return;
  nedge = nedg;
  if (edgevars_allocated) {
    freeEdgevars();
    allocateEdgevars();
  }
}


#ifdef TIMING
void PDFactoryDouble::resetTime() {
  RTime = p0Time = fftTime = posTime = copyTime = normTime = 0;
  meanTime = logTime = 0;
}

void PDFactoryDouble::summarizeTime(unsigned int nindent) const {
  std::string prestring(nindent,' ');
    
  std::cout << "R time: " << prestring 
	    << 1.0*RTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "p0 time: " << prestring 
	    << 1.0*p0Time/CLOCKS_PER_SEC << std::endl;
  std::cout << "fft time: " << prestring 
	    << 1.0*fftTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "pos time: " << prestring 
	    << 1.0*posTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "copy time: " << prestring 
	    << 1.0*copyTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "norm time: " << prestring 
	    << 1.0*normTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "mean time: " << prestring 
	    << 1.0*meanTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "log time: " << prestring 
	    << 1.0*logTime/CLOCKS_PER_SEC << std::endl;
}
#endif


/*
  Doesn't force resize if there is already enough room available

  \param[in] NSIZE new size
  \returns True if a resize was needed
*/
//I don't check for 0 NSIZE because that should never happen
// in practice, and it isn't worth the idiot-proofing
// rtrans always gets nulled in here, then refilled when you
//  call initPD
bool PDFactoryDouble::resize(unsigned int NSIZE) {
  if (NSIZE == currsize) return false;
  freeRvars(); // Also sets initialized variables
  currsize = NSIZE;
  return true;
}

void PDFactoryDouble::allocateRvars() {
  if (rvars_allocated) return;
  if (currsize == 0)
    throw pofdExcept("PDFactoryDouble","allocate_rvars",
		     "Invalid (0) currsize",1);
  RFlux1 = (double*) fftw_malloc(sizeof(double) * currsize);
  RFlux2 = (double*) fftw_malloc(sizeof(double) * currsize);
  rsum = (double*) fftw_malloc(sizeof(double) * currsize);
  unsigned int fsize = currsize * currsize;
  rvals = (double*) fftw_malloc(sizeof(double) * fsize);
  pofd  = (double*) fftw_malloc(sizeof(double) * fsize);
  fsize = currsize * (currsize / 2 + 1);
  rtrans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fsize);
  pval = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fsize);

  plans_valid = false;
  rvars_allocated = true;
  rinitialized = false;
  initialized = false;
}

void PDFactoryDouble::freeRvars() {
  if (RFlux1 != nullptr) { fftw_free(RFlux1); RFlux1=nullptr; }
  if (RFlux2 != nullptr) { fftw_free(RFlux2); RFlux2=nullptr; }
  if (rvals != nullptr) { fftw_free(rvals); rvals=nullptr; }
  if (rsum != nullptr) { fftw_free(rsum); rsum=nullptr; }
  if (rtrans != nullptr) { fftw_free(rtrans); rtrans=nullptr; }
  if (pval != nullptr) { fftw_free(pval); pval = nullptr; }
  if (pofd != nullptr) { fftw_free(pofd); pofd = nullptr; }
  plans_valid = false;
  rvars_allocated = false;
  rinitialized = false;
  initialized = false;
}

void PDFactoryDouble::allocateEdgevars() {
  if (edgevars_allocated) return;
  if (nedge > 0) {
    REdgeFlux1 = (double*) fftw_malloc(sizeof(double) * nedge);
    REdgeFlux2 = (double*) fftw_malloc(sizeof(double) * nedge);
    REdgeWork  = (double*) fftw_malloc(sizeof(double) * nedge * nedge);
    edgevars_allocated = true;
  } else {
    REdgeWork = REdgeFlux1 = REdgeFlux2 = nullptr;
    edgevars_allocated = false;
  }
  rinitialized = false;
  initialized = false; 
}

void PDFactoryDouble::freeEdgevars() {
  if (REdgeFlux1 != nullptr) { fftw_free(REdgeFlux1); REdgeFlux1 = nullptr; }
  if (REdgeFlux2 != nullptr) { fftw_free(REdgeFlux2); REdgeFlux2 = nullptr; }
  if (REdgeWork != nullptr) { fftw_free(REdgeWork); REdgeWork = nullptr; }
  edgevars_allocated = false;
  rinitialized = false;
  initialized = false;
}

/*!
  Frees all internal memory
*/
void PDFactoryDouble::free() {
  freeRvars();
  freeEdgevars();
}

/*!
  \param[in] filename Name of wisdom file
*/
bool PDFactoryDouble::addWisdom(const std::string& filename) {
  if (fftw_import_wisdom_from_filename(filename.c_str()) == 0) {
    std::stringstream str;
    str << "Error reading wisdom file: " << filename;
    throw pofdExcept("PDFactoryDouble", "addWisdom", str.str(), 2);
  }
  fftw_plan_style = FFTW_WISDOM_ONLY;
  has_wisdom = true;
  wisdom_file = filename;
  if (plan != nullptr) {
    fftw_destroy_plan(plan); 
    plan = nullptr;
  }
  if (plan_inv != nullptr) {
    fftw_destroy_plan(plan_inv);
    plan_inv = nullptr;
  }
  // Doesn't affect R initialization status, but the plans are definitely
  // no longer valid
  plans_valid = false;
  initialized = false;
  return true;
}

/*!
  \param[in] n Number of elements
  \param[in] minflux1 Minimum flux band 1 to use in R
  \param[in] maxflux1 Maximum flux band 1 to use in R
  \param[in] minflux2 Minimum flux band 2 to use in R
  \param[in] maxflux2 Maximum flux band 2 to use in R

  Sets up RFlux1 and Rflux2.  They may not quite get to minflux/maxflux
  if minflux is negative.  
*/
void PDFactoryDouble::initRFlux(unsigned int n, double minflux1, 
				double maxflux1, double minflux2,
				double maxflux2) {
  // Make sure there is room
  resize(n);

  if (n == 0)
    throw pofdExcept("PDFactoryDouble", "initR", "Invalid (0) n", 1);
  if (n == 1) {
    dflux1 = dflux2 = 0.1;
    minfluxR_1 = RFlux1[0] = minflux1;
    minfluxR_2 = RFlux2[0] = minflux2;
    return;
  }

  if (maxflux1 < minflux1) std::swap(minflux1, maxflux1);
  if (maxflux2 < minflux2) std::swap(minflux2, maxflux2);

  dflux1 = initRFluxInternal(n, minflux1, maxflux1, RFlux1, minfluxR_1);
  dflux2 = initRFluxInternal(n, minflux2, maxflux2, RFlux2, minfluxR_2);

  rinitialized = false;
  initialized = false;
}

/*!
  \param[in] n        Size of transform 
  \param[in] minflux1 Minimum flux to use in R, band 1
  \param[in] maxflux1 Maximum flux to use in R, band 1
  \param[in] minflux2 Minimum flux to use in R, band 2
  \param[in] maxflux2 Maximum flux to use in R, band 2
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] bm      Histogrammed inverse beam
  \param[in] setEdge Set the edge of R using edge integration
  \param[in] muldflux Multiply R by dflux1 * dflux2

  dflux1, dflux2 is also set, are are RFlux1, RFlux2 and the edge vars
*/
void PDFactoryDouble::initR(unsigned int n, double minflux1, double maxflux1, 
			    double minflux2, double maxflux2, 
			    const numberCountsDouble& model,
			    const doublebeamHist& bm, bool setEdge,
			    bool muldflux) {
  
  // This version is much more complex than the 1D case because of the
  // edge bits
  if (!rvars_allocated) allocateRvars();

  // Fill in RFlux values
  initRFlux(n, minflux1, maxflux1, minflux2, maxflux2);

  //Now fill in R.  The edges (e.g., between f=0 and the first element)
  // require special care.  This first call will fill zero values into 
  // the lower edges, which we will later overwrite if we are doing setEdge
  //Note that we do -not- multiply by dflux yet because of the edge stuff.
  double *rptr;  

#ifdef TIMING
  starttime = std::clock();
#endif
  model.getR(n, RFlux1, n, RFlux2, bm, rvals);
#ifdef TIMING
  RTime += std::clock() - starttime;
#endif

  if (setEdge) {
    //Now fill in the lower edges by doing integrals
    // and setting to the mean value inside that.
    //Use the trapezoidal rule in either log or linear flux
    // space

    if (nedge == 0) 
      throw pofdExcept("PDFactoryDouble", "initR",
		       "Invalid (0) value of nedge with setEdge set", 1);

    //Edge bits
    //Minimum values in integral; maximum are dflux1, dflux2
    double minedge1 = dflux1 *  PDFactoryDouble::lowEdgeRMult;
    double minedge2 = dflux2 *  PDFactoryDouble::lowEdgeRMult;
    double inedgem1 = 1.0 / static_cast<double>(nedge - 1);
    double dinterpfluxedge1, dinterpfluxedge2;
    double iRxnorm = 0.0, iRynorm = 0.0, iR00norm = 0.0;
    
    if (!edgevars_allocated) allocateEdgevars();
    //Make sure edge variables successfully allocated
    if (REdgeFlux1 == nullptr)
      throw pofdExcept("PDFactoryDouble", "initR",
		       "R edge flux 1 was not allocated", 2);
    if (REdgeFlux2 == nullptr)
      throw pofdExcept("PDFactoryDouble", "initR",
		       "R edge flux 2 was not allocated", 3);
    if (REdgeWork == nullptr)
      throw pofdExcept("PDFactoryDouble", "initR",
		       "R edge work was not allocated", 4);

    if (use_edge_log_x) {
      dinterpfluxedge1 = -log(PDFactoryDouble::lowEdgeRMult) * inedgem1;
      for (unsigned int i = 0; i < nedge; ++i)
	REdgeFlux1[i] = minedge1*exp(static_cast<double>(i)*dinterpfluxedge1);
    } else {
      dinterpfluxedge1 = (dflux1-minedge1)*inedgem1;
      for (unsigned int i = 0; i < nedge; ++i)
	REdgeFlux1[i] = minedge1 + static_cast<double>(i)*dinterpfluxedge1;
    }
    if (use_edge_log_y) {
      dinterpfluxedge2 = -log(PDFactoryDouble::lowEdgeRMult)*inedgem1;
      for (unsigned int i = 0; i < nedge; ++i)
	REdgeFlux2[i] = minedge2*exp(static_cast<double>(i)*dinterpfluxedge2);
    } else {
      dinterpfluxedge2 = (dflux2-minedge2)*inedgem1;
      for (unsigned int i = 0; i < nedge; ++i)
	REdgeFlux2[i] = minedge2 + static_cast<double>(i)*dinterpfluxedge2;
    }
    iRxnorm  = dinterpfluxedge1 / (dflux1 - minedge1);
    iRynorm  = dinterpfluxedge2 / (dflux2 - minedge2);
    iR00norm = dinterpfluxedge1 * dinterpfluxedge2/
      ((dflux1 - minedge1) * (dflux2 - minedge2));

    //First, do r[0,0]
    double scriptr;
#ifdef TIMING
    starttime = std::clock();
#endif
    model.getR(nedge, REdgeFlux1, nedge, REdgeFlux2, bm, REdgeWork);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif

    //Do y integral first, store in REdgeWork[0,*]
    if (use_edge_log_y) {
      for (unsigned int i = 0; i < nedge; ++i) {
	rptr = REdgeWork + i*nedge; //row pointer
	scriptr = 0.5*(REdgeFlux2[0] * rptr[0] + 
		       REdgeFlux2[nedge-1] * rptr[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += REdgeFlux2[j] * rptr[j];
	REdgeWork[i] = scriptr;
      }
    } else {
      for (unsigned int i = 0; i < nedge; ++i) {
	rptr = REdgeWork + i*nedge; 
	scriptr = 0.5*(rptr[0] + rptr[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += rptr[j];
	REdgeWork[i] = scriptr;
      }
    }
    //Now X integral, put in integral step size and area
    // of bin, store in R[0,0]
    if (use_edge_log_x) {
      scriptr = 0.5*(REdgeFlux1[0] * REdgeWork[0]+
		     REdgeFlux1[nedge-1] * REdgeWork[nedge-1]);
      for (unsigned int i = 1; i < nedge-1; ++i)
	scriptr += REdgeFlux1[i] * REdgeWork[i];
      rvals[0] = scriptr * iR00norm;
    } else {
      scriptr = 0.5*(REdgeWork[0] + REdgeWork[nedge-1]);
      for (unsigned int i = 1; i < nedge-1; ++i)
	scriptr += REdgeWork[i];
      rvals[0] = scriptr * iR00norm;
    }
    
    //Now do Rx = R[0,y], integral along x
    double fixed_value;
    for (unsigned int j = 1; j < n; ++j) {
      fixed_value = RFlux2[j];
      //REdgeWork is more than big enough
#ifdef TIMING
      starttime = std::clock();
#endif
      model.getR(nedge, REdgeFlux1, 1, &fixed_value, bm, REdgeWork);
#ifdef TIMING
      RTime += std::clock() - starttime;
#endif
      if (use_edge_log_x) {
	scriptr = 0.5*(REdgeFlux1[0]*REdgeWork[0]+
		       REdgeFlux1[nedge-1]*REdgeWork[nedge-1]);
	for (unsigned int i = 1; i < nedge-1; ++i)
	  scriptr += REdgeFlux1[i]*REdgeWork[i];
      } else {
	scriptr = 0.5*(REdgeWork[0]+REdgeWork[nedge-1]);
	for (unsigned int i = 1; i < nedge-1; ++i)
	  scriptr += REdgeWork[i];
      }
      rvals[j] = scriptr*iRxnorm;
    }
    
    //And Ry = R[x,0]
    for (unsigned int i = 1; i < n; ++i) {
      fixed_value = RFlux1[i];
      //REdgeWork is more than big enough
#ifdef TIMING
      starttime = std::clock();
#endif
      model.getR(1, &fixed_value, nedge, REdgeFlux2, bm, REdgeWork);
#ifdef TIMING
      RTime += std::clock() - starttime;
#endif
      if (use_edge_log_y) {
	scriptr = 0.5*(REdgeFlux2[0]*REdgeWork[0]+
		       REdgeFlux2[nedge-1]*REdgeWork[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += REdgeFlux2[j]*REdgeWork[j];
      } else {
	scriptr = 0.5*(REdgeWork[0]+REdgeWork[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += REdgeWork[j];
      }
      rvals[i*n] = scriptr*iRynorm;
    }
  }

  //Multiply R by dflux1 * dflux2
  if (muldflux) {
    double fluxfac = dflux1 * dflux2;
    for (unsigned int i = 0; i < n * n; ++i)
      rvals[i] *= fluxfac;
    rdflux = true;
  } else rdflux = false;
  
  initialized = false;
  rinitialized = true;
}

/*!
  \param[in] model Number counts model
  \param[in] bm Inverse beam histogram
  \returns Estimate of min/max flux values where R is non-zero in each band.
*/
std::pair<dblpair, dblpair> 
PDFactoryDouble::getMinMaxR(const numberCountsDouble& model, 
			    const doublebeamHist& bm) const {
  const double safetyfac = 1.01;

  double minflux1, maxflux1, minflux2, maxflux2;
  minflux1 = maxflux1 = minflux2 = maxflux2 = 0.0;
  
  dblpair maxf = model.getMaxFluxEstimate();
  double cval;
  if (bm.hasSign(0)) {
    // pos/pos bit.  Affects maximum values in bands 1 and 2
    cval = safetyfac * maxf.first * bm.getMinMax1(0).second;
    if (cval > maxflux1) maxflux1 = cval;
    cval = safetyfac * maxf.second * bm.getMinMax2(0).second;
    if (cval > maxflux2) maxflux2 = cval;
  }
  if (bm.hasSign(1)) {
    //pos/neg bit.  Affects maximum in band 1, minimum in band 2
    cval = safetyfac * maxf.first * bm.getMinMax1(1).second;
    if (cval > maxflux1) maxflux1 = cval;
    cval = -safetyfac * maxf.second * bm.getMinMax2(1).second;
    if (cval < minflux2) minflux2 = cval;
  }
  if (bm.hasSign(2)) {
    //neg/pos bit.  Affects minimum in band 1, maximum in band 2
    cval = -safetyfac * maxf.first * bm.getMinMax1(2).second;
    if (cval < minflux1) minflux1 = cval;
    cval = safetyfac * maxf.second * bm.getMinMax2(2).second;
    if (cval > maxflux2) maxflux2 = cval;
  }
  if (bm.hasSign(3)) {
    //neg/neg bit.  Affects minimum in bands 1 and 2
    cval = -safetyfac * maxf.first * bm.getMinMax1(3).second;
    if (cval < minflux1) minflux1 = cval;
    cval = -safetyfac * maxf.second * bm.getMinMax2(3).second;
    if (cval < minflux2) minflux2 = cval;
  }

  return std::make_pair(std::make_pair(minflux1, maxflux1),
			std::make_pair(minflux2, maxflux2));
}

/*!
  \param[in] n Number of elements to use
  \param[in] model Number counts model
  \param[in] bm Inverse histogrammed beam
  \param[in] minmaxR A pair of pairs of min/max R ranges (see getMinMaxR)
  \returns A pair of pairs containing the mean and variance for each dimension

  Computes mean and variance of the P(D) using R along each dimension.
  This also fills R. The instrument noise is not included
  in the returned moments.  
*/
std::pair<dblpair, dblpair>
  PDFactoryDouble::getRMoments(unsigned int n, const numberCountsDouble& model, 
			       const doublebeamHist& bm, 
			       std::pair<dblpair, dblpair>& minmaxR,
			       bool setEdge) {

  // Recall that the formulae for the mean and central 2nd moment
  // along each axis (not including instrumental noise) are
  //  <x> = \int x R dx dy
  //  <y> = \int y R dx dy
  //  <(x - <x>)^2> = \int x^2 R dx dy
  //  <(y - <y>)^2> = \int y^2 R dx dy

  if (n == 0)
    throw pofdExcept("PDFactoryDouble", "getRMoments",
		     "Invalid n (== 0)", 1);

  // Now fill in, set up R, Rflux, etc. in all their glory.
  // This will also resize as needed.
  initR(n, minmaxR.first.first, minmaxR.first.second,
	minmaxR.second.first, minmaxR.second.second,
	model, bm, setEdge, false);

  // And... trap rule the integrals.
  // Always move along j since array access is
  //  faster that way (rvals is row-major order)
  // This is a simple calculation, somewhat tedious to write out.
  // Note that R is already multiplied by dflux1 * dflux2 by initR

  double cf1, cf2, cR, mean1, mean2, var1, var2, prod1, prod2;
  double *rowptr;
  
  // First row, has a 0.5 factor
  rowptr = rvals;
  cf1 = RFlux1[0];
  cR = 0.5 * rowptr[0]; // This handles the 0.5 factor
  prod1 = 0.5 * cf1 * cR; // Extra 0.5 for first col
  mean1 = prod1;
  var1 = cf1 * prod1;
  cf2 = RFlux2[0];
  prod2 = 0.5 * cf2 * cR;
  mean2 = prod2;
  var2 = cf2 * prod2;
  for (unsigned int j = 1; j < n-1; ++j) {
    cR = 0.5 * rowptr[j]; // Still 0.5 for first row
    prod1 = cf1 * cR;
    mean1 += prod1;
    var1 += cf1 * prod1;
    cf2 = RFlux2[j];
    prod2 = cf2 * cR;
    mean2 += prod2;
    var2 += cf2 * prod2;
  }
  cR = 0.5 * rowptr[n-1];
  prod1 = 0.5 * cf1 * cR;
  mean1 += prod1;
  var1 += cf1 * prod1;
  cf2 = RFlux2[n-1];
  prod2 = 0.5 * cf2 * cR;
  mean2 += prod2;
  var2 += cf2 * prod2;

  // Do middle rows
  // Note that, while R is allocated to rsize by rsize, we only fill up
  // to n by n
  for (unsigned int i = 1; i < n-1; ++i) {
    rowptr = rvals + i * n; 
    cf1 = RFlux1[i];
    cR = rowptr[0]; // No 0.5 any more -- we are in the middle
    prod1 = 0.5 * cf1 * cR; // 0.5 for first col
    mean1 += prod1;
    var1 += cf1 * prod1;
    cf2 = RFlux2[0];
    prod2 = 0.5 * cf2 * cR;
    mean2 += prod2;
    var2 += cf2 * prod2;
    for (unsigned int j = 1; j < n-1; ++j) {
      cR = rowptr[j];
      prod1 = cf1 * cR;
      mean1 += prod1;
      var1 += cf1 * prod1;
      cf2 = RFlux2[j];
      prod2 = cf2 * cR;
      mean2 += prod2;
      var2 += cf2 * prod2;
    }
    cR = rowptr[n-1];
    prod1 = 0.5 * cf1 * cR;
    mean1 += prod1;
    var1 += cf1 * prod1;
    cf2 = RFlux2[n-1];
    prod2 = 0.5 * cf2 * cR;
    mean2 += prod2;
    var2 += cf2 * prod2;
  }

  // And final row, has an extra 0.5 factor again
  rowptr = rvals + (n - 1) * n;
  cf1 = RFlux1[n-1];
  cR = 0.5 * rowptr[0]; 
  prod1 = 0.5 * cf1 * cR;
  mean1 += prod1;
  var1 += cf1 * prod1;
  cf2 = RFlux2[0];
  prod2 = 0.5 * cf2 * cR;
  mean2 += prod2;
  var2 += cf2 * prod2;
  for (unsigned int j = 1; j < n-1; ++j) {
    cR = 0.5 * rowptr[j];
    prod1 = cf1 * cR;
    mean1 += prod1;
    var1 += cf1 * prod1;
    cf2 = RFlux2[j];
    prod2 = cf2 * cR;
    mean2 += prod2;
    var2 += cf2 * prod2;
  }
  cR = 0.5 * rowptr[n-1];
  prod1 = 0.5 * cf1 * cR;
  mean1 += prod1;
  var1 += cf1 * prod1;
  cf2 = RFlux2[n-1];
  prod2 = 0.5 * cf2 * cR;
  mean2 += prod2;
  var2 += cf2 * prod2;

  // Apply dflux factors
  double dfact = dflux1 * dflux2;
  mean1 *= dfact;
  mean2 *= dfact;
  var1 *= dfact;
  var2 *= dfact;

  return std::make_pair(std::make_pair(mean1, var1),
			std::make_pair(mean2, var2));
}
  
/*!
  \param[in] n0 Value of n0 in use
  \param[in] n Number of elements
  \param[out] pd Holds P(D) on output, normalized, mean subtracted,
                 and with positivity enforced.

  This does the unwrapping of the internal P(D) into pd.
*/
// This should only ever be called by getPD, so we don't really
// check the inputs
void PDFactoryDouble::unwrapPD(double n0, unsigned int n, PDDouble& pd) const {
  // This is similar to the 1D case -- we want to unwrap the PD where
  //  the negative and positive bits overlap.  But using the minimum
  //  is unstable due to numeric 'noise' so we instead look for the PD
  //  being down by a certain amount relative to peak.
  const double peakfrac = 1e-12; // How far down from the peak we want

  // Our acceptance testing is a bit complicated.
  // If the minimum point is more than nsig1 away from the expected mean (0)
  //  then we just accept it.  If it is more than nsig2 away, we make sure
  //  that the min/max ratio of the P(D) along that axis is more than
  //  maxminratio.  If it is less than nsig2, we just flat out reject.
  const double nsig1 = 4.0; 
  const double nsig2 = 2.0;
  const double minmaxratio = 1e-5;

  // Recall that varnoi1, varnoi2 were computed for the max n0, not the 
  // current one.  So fix that internally
  double n0ratio = n0 / max_n0;
  double curr_sigma1 = sqrt(n0ratio * varnoi1 + sigma1 * sigma1);
  double curr_sigma2 = sqrt(n0ratio * varnoi2 + sigma2 * sigma2);

  //Enforce positivity
#ifdef TIMING
  starttime = std::clock();
#endif
  for (unsigned int idx = 0; idx < n * n; ++idx)
    if (pofd[idx] < 0) pofd[idx] = 0.0;
#ifdef TIMING
  posTime += std::clock() - starttime;
#endif

  // Figure out the indices of the split point in each dimension
  // Rather than find the 2D index, which would be noisier, 
  //  intergrate out each axis to find the minimum, storing in rsum.
  // Start from the top and move down -- if there is a tie for
  //  some reason we prefer the index be high because that is more
  //  likely to be right in practice
  // Start with min along the first index.  We also find the max, but
  // don't keep track of its index.
#ifdef TIMING
  starttime = std::clock();
#endif
  // Sum into rsum, making this a 1D problem
  int nm1 = static_cast<int>(n - 1);
  double cval, *rowptr;
  for (unsigned int i = 0; i < n; ++i) {
    rowptr = pofd + i * n;
    cval = 0.5 * rowptr[0];
    for (unsigned int j = 1; j < n - 1; ++j)
      cval += rowptr[j];
    cval += 0.5 * rowptr[n-1];
    rsum[i] = cval;
  }
  int mdx = nm1; // Curr min index
  double minval, maxval;
  minval = maxval = rsum[nm1];
  for (int i = n - 2; i >= 0; --i) {
    cval = rsum[i];
    if (cval < minval) {
      minval = cval;
      mdx = i;
    } else if (cval > maxval) maxval = cval;
  }
  unsigned int splitidx1 = static_cast<unsigned int>(mdx);
  double splitval1 = minval;
  if (minval / maxval < peakfrac) {
    // There should therefore be a more robust point that isn't so
    //  far down -- look for it from the top again
    // There is such a point -- look for it again from the top
    double targval = peakfrac * maxval;
    for (int i = n - 1; i >= 0; --i)
      if (rsum[i] <= targval) {
	splitidx1 = static_cast<unsigned int>(i);
	splitval1 = rsum[i];
	break;
      }
  }

  // Sanity check stuff
  double fwrap_plus = 
    static_cast<double>(splitidx1) * dflux1; // Wrap in pos flux
  double fwrap_minus 
    = static_cast<double>(n - splitidx1) * dflux1; // Abs neg wrap
  double cs1, cs2;
  cs1 = nsig1 * curr_sigma1;
  cs2 = nsig2 * curr_sigma1;

  // If fwrap_plus/minus are within cs2 of the peak, we will always
  //  fail the tests below.  So, in that case, we try to adjust the
  //  split point out.  This can work because the minimum of the P(D)
  //  is very flat in some cases, and so it doesn't cost us to try it.
  //  We do the same in the 1D code (PDFactory::unwrap)
  if (fwrap_plus < cs2) {
    unsigned int splitidx1_trial = 
      static_cast<unsigned int>((cs2 - fwrap_plus) / dflux1) + 1 + splitidx1;
    // Can only go up to n
    if (splitidx1_trial >= n) splitidx1_trial = n - 1;
    double splitval1_trial = rsum[splitidx1_trial];
    if ((splitval1_trial / maxval) < minmaxratio) {
      // Worth doing as an attempt to save things
      splitidx1 = splitidx1_trial;
      splitval1 = splitval1_trial;
    }
    fwrap_plus = static_cast<double>(splitidx1) * dflux1; // Wrap in pos flux
    fwrap_minus = static_cast<double>(n - splitidx1) * dflux1;
  } else if (fwrap_minus < cs2) {
    unsigned int splitidx1_delta = 
      static_cast<unsigned int>((cs2 - fwrap_minus) / dflux1) + 1;
    if (splitidx1_delta > splitidx1) splitidx1_delta = splitidx1;
    unsigned int splitidx1_trial = splitidx1 - splitidx1_delta;
    double splitval1_trial = rsum[splitidx1_trial];
    if ((splitval1_trial / maxval) < minmaxratio) {
      splitidx1 = splitidx1_trial;
      splitval1 = splitval1_trial;
    }
    fwrap_plus = static_cast<double>(splitidx1) * dflux1; // Wrap in pos flux
    fwrap_minus = static_cast<double>(n - splitidx1) * dflux1;
  }

  // Now the checks
  if ((fwrap_plus < cs1) || (fwrap_minus < cs1)) {
    // Worth further investigation
    if (fwrap_plus < cs2) {
      std::stringstream errstr;
      errstr << "Top wrapping problem dim 1; wrapping point at "
	     << fwrap_plus << " which is only " << fwrap_plus / curr_sigma1
	     << " sigma away from expected (0) mean with sigma "
	     << curr_sigma1 << " at n0: " << n0;
      throw pofdExcept("PDFactoryDouble", "unwrapPD", errstr.str(), 1);
    }
    if (fwrap_minus < cs2) {
      std::stringstream errstr;
      errstr << "Bottom wrapping problem dim 1; wrapping point at "
	     << -fwrap_minus << " which is only " << fwrap_minus / curr_sigma1
	     << " sigma away from expected (0) mean, with sigma "
	     << curr_sigma1 << " at n0: " << n0;
      throw pofdExcept("PDFactoryDouble", "unwrapPD", errstr.str(), 2);
    } 
    // Min/max ratio test
    if (splitval1 / maxval > minmaxratio) {
      std::stringstream errstr;
      errstr << "Dim 1 wrapping problem with wrapping fluxes: "
	     << fwrap_plus << " and " << -fwrap_minus << " with min/max ratio: "
	     << splitval1 / maxval << " and sigma: " << curr_sigma1
	     << " with n0: " << n0;
      throw pofdExcept("PDFactoryDouble", "unwrapPD", errstr.str(), 3);
    }
  }
  
  // Now second dimension.  This one is slower due to stride issues,
  //  but otherwise is the same. Doing the sum in this order is much faster
  for (unsigned int j = 0; j < n; ++j) rsum[j] = 0.5 * pofd[j]; // i=0
  for (unsigned int i = 1; i < n - 1; ++i) { // center
    rowptr = pofd + i * n;
    for (unsigned int j = 0; j < n; ++j)
      rsum[j] += rowptr[j];
  }
  rowptr = pofd + nm1 * n; // i = n-1
  for (unsigned int j = 0; j < n; ++j)
    rsum[j] += 0.5 * rowptr[j];
  mdx = nm1; // Curr min index
  minval = maxval = rsum[nm1];
  for (int i = n - 2; i >= 0; --i) {
    cval = rsum[i];
    if (cval < minval) {
      minval = cval;
      mdx = i;
    } else if (cval > maxval) maxval = cval;
  }
  unsigned int splitidx2 = static_cast<unsigned int>(mdx);
  double splitval2 = minval;
  if (minval / maxval < peakfrac) {
    double targval = peakfrac * maxval;
    for (int i = n - 1; i >= 0; --i)
      if (rsum[i] <= targval) {
	splitidx2 = static_cast<unsigned int>(i);
	splitval2 = rsum[i];
	break;
      }
  }

  // Same sanity checks/tweaks
  fwrap_plus = static_cast<double>(splitidx2) * dflux2;
  fwrap_minus = static_cast<double>(n - splitidx2) * dflux2;
  cs1 = nsig1 * curr_sigma2;
  cs2 = nsig2 * curr_sigma2;
  if (fwrap_plus < cs2) {
    unsigned int splitidx2_trial = 
      static_cast<unsigned int>((cs2 - fwrap_plus) / dflux2) + 1 + splitidx2;
    if (splitidx2_trial >= n) splitidx2_trial = n - 1;
    double splitval2_trial = rsum[splitidx2_trial];
    if ((splitval2_trial / maxval) < minmaxratio) {
      splitidx2 = splitidx2_trial;
      splitval2 = splitval2_trial;
    }
    fwrap_plus = static_cast<double>(splitidx2) * dflux2;
    fwrap_minus = static_cast<double>(n - splitidx2) * dflux2;
  } else if (fwrap_minus < cs2) {
    unsigned int splitidx2_delta = 
      static_cast<unsigned int>((cs2 - fwrap_minus) / dflux2) + 1;
    if (splitidx2_delta > splitidx2) splitidx2_delta = splitidx2;
    unsigned int splitidx2_trial = splitidx2 - splitidx2_delta;
    double splitval2_trial = rsum[splitidx2_trial];
    if ((splitval2_trial / maxval) < minmaxratio) {
      splitidx2 = splitidx2_trial;
      splitval2 = splitval2_trial;
    }
    fwrap_plus = static_cast<double>(splitidx2) * dflux2;
    fwrap_minus = static_cast<double>(n - splitidx2) * dflux2;
  }

  if ((fwrap_plus < cs1) || (fwrap_minus < cs1)) {
    // Worth further investigation
    if (fwrap_plus < cs2) {
      std::stringstream errstr;
      errstr << "Top wrapping problem dim 2; wrapping point at "
	     << fwrap_plus << " which is only " << fwrap_plus / curr_sigma2
	     << " sigma away from expected (0) mean with sigma "
	     << curr_sigma2 << " at n0: " << n0;
      throw pofdExcept("PDFactoryDouble", "unwrapPD", errstr.str(), 4);
    }
    if (fwrap_minus < cs2) {
      std::stringstream errstr;
      errstr << "Bottom wrapping problem dim 2; wrapping point at "
	     << -fwrap_minus << " which is only " << fwrap_minus / curr_sigma2
	     << " sigma away from expected (0) mean, with sigma "
	     << curr_sigma2 << " at n0: " << n0;
      throw pofdExcept("PDFactoryDouble", "unwrapPD", errstr.str(), 5);
    } 
    // Min/max ratio test
    if (splitval2 / maxval > minmaxratio) {
      std::stringstream errstr;
      errstr << "Dim 2 wrapping problem with wrapping fluxes: "
	     << fwrap_plus << " and " << -fwrap_minus << " with min/max ratio: "
	     << splitval2 / maxval << " and sigma: " << curr_sigma2
	     << " with n0: " << n0;
      throw pofdExcept("PDFactoryDouble", "unwrapPD", errstr.str(), 6);
    }
  }

  // Now the actual copying, which is an exercise in index gymnastics
  pd.resize(n, n);
  size_t size_double = sizeof(double); // Size of double
  std::memset(pd.pd_, 0, n * n * size_double);

  size_t rowsz;
  double *ptr_curr, *ptr_out; // convenience pointers
  // We start with the neg, neg bit -- that is stuff >= both
  // splitidx1 and splitidx2, which ends up in the 0, 0 part of the output
  ptr_out = pd.pd_;
  ptr_curr = pofd + splitidx1 * n + splitidx2;
  rowsz = (n - splitidx2) * size_double;
  for (unsigned int i = 0; i < n - splitidx1; ++i)
    std::memcpy(ptr_out + i * n, ptr_curr + i * n, rowsz);

  // Next pos neg
  ptr_out = pd.pd_ + (n - splitidx1) * n;
  ptr_curr = pofd + splitidx2;
  for (unsigned int i = 0; i < splitidx1; ++i) 
    std::memcpy(ptr_out + i * n, ptr_curr + i * n, rowsz);

  // pos, pos
  ptr_out = pd.pd_ + (n - splitidx1) * n + (n - splitidx2);
  ptr_curr = pofd;
  rowsz = splitidx2 * size_double;
  for (unsigned int i = 0; i < splitidx1; ++i)
    std::memcpy(ptr_out + i * n, ptr_curr + i * n, rowsz);
  
  // neg, pos
  ptr_out = pd.pd_ + (n - splitidx2);
  ptr_curr = pofd + splitidx1 * n;
  for (unsigned int i = 0; i < n - splitidx1; ++i)
    std::memcpy(ptr_out + i * n, ptr_curr + i * n, rowsz);

  pd.logflat = false;
  pd.minflux1 = 0.0; pd.dflux1 = dflux1;
  pd.minflux2 = 0.0; pd.dflux2 = dflux2;

#ifdef TIMING
  copyTime += std::clock() - starttime;
#endif

  // Normalize
#ifdef TIMING
  starttime = std::clock();
#endif
  pd.normalize();
#ifdef TIMING
  normTime += std::clock() - starttime;
#endif

  // Mean subtract axes
#ifdef TIMING
  starttime = std::clock();
#endif
  double tmn1, tmn2;
  pd.getMeans(tmn1, tmn2, false);
  if (std::isinf(tmn1) || std::isnan(tmn1) || std::isinf(tmn2) ||
       std::isnan(tmn2)) {
    std::stringstream str;
    str << "Un-shift amounts not finite band1: " << tmn1 << " "
        << tmn2 << std::endl;
    str << "At length: " << n << " with noise: " << sigma1 << " "
        << sigma2;
    throw pofdExcept("PDFactoryDouble", "unwrapPD", str.str(), 7);
  }
  pd.minflux1 = -tmn1;
  pd.minflux2 = -tmn2;
#ifdef TIMING
  meanTime += std::clock() - starttime;
#endif

}

/*!
  Gets ready for P(D) computation by preparing R and transforming it
 
  \param[in] n Size of transform (square)
  \param[in] inst_sigma1 Instrument sigma, band 1
  \param[in] inst_sigma2 Instrument sigma, band 2
  \param[in] maxflux1 Desired maximum flux for R, dimension 1
  \param[in] maxflux2 Desired maximum flux for R, dimension 2
  \param[in] maxn0   Maximum m0 supported (number of sources per area)
  \param[in] model    number counts model to use for fill.  Params must be set
  \param[in] bm       Beam
  \param[in] setEdge  Use integral of mean values at the edges
*/
void PDFactoryDouble::initPD(unsigned int n, 
			     double inst_sigma1, double inst_sigma2, 
			     double maxflux1, double maxflux2, 
			     double maxn0, const numberCountsDouble& model,
			     const doublebeamHist& bm, bool setEdge) {

  if (n == 0)
    throw pofdExcept("PDFactoryDouble", "initPD",
		     "Invalid (non-positive) n",1);  

  if (inst_sigma1 < 0.0)
    throw pofdExcept("PDFactoryDouble", "initPD",
		     "Invalid (negative) inst_sigma1",2);
  if (inst_sigma2 < 0.0)
    throw pofdExcept("PDFactoryDouble", "initPD",
		     "Invalid (negative) inst_sigma2",3);
  if (maxflux1 <= 0.0)
    throw pofdExcept("PDFactoryDouble", "initPD",
		     "Invalid (non-positive) maxflux1",4);
  if (maxflux2 <= 0.0)
    throw pofdExcept("PDFactoryDouble", "initPD",
		     "Invalid (non-positive) maxflux2",5);
  if (!model.isValid())
    throw pofdExcept("PDFactoryDouble", "initPD", "model not valid", 6);  
  if (maxn0 <= 0.0)
    throw pofdExcept("PDFactoryDouble", "initPD", 
		     "Invalid (non-positive) n0", 7);

  //Must allocate space before we plan -- including the R variables
  resize(n); 
  allocateRvars();
  
  //Make the plans, or keep the old ones if possible
  // Do this before doing R fill, as plan construction can
  // mess with the input/output variables
  //Note that the forward transform dumps into out_part,
  // but the backwards one comes from pval.  The idea
  // is that out_part holds the working bit.  These are
  // convolved together into pval.  This is inefficient
  // if there is only one sign present, but the special
  // case doesn't seem worth the effort
  //If we resized, we must make the new plans because the
  // addresses changed
  //We will have to use the advanced interfact to point
  // specifically at the rvals subindex we are using on the
  // forward plan, but the backwards plan is fine
  int intn = static_cast<int>(n);
  if (!plans_valid) {
    if (plan != nullptr) fftw_destroy_plan(plan);
    plan = fftw_plan_dft_r2c_2d(intn, intn, rvals, rtrans,
				fftw_plan_style);
    if (plan == nullptr) {
      std::stringstream str;
      str << "Plan creation failed for forward transform of size: " << intn
	  << " by " << intn;
      if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			  << " that size";
      throw pofdExcept("PDFactoryDouble", "initPD", str.str(), 14);
    }

    if (plan_inv != nullptr) fftw_destroy_plan(plan_inv);
    plan_inv = fftw_plan_dft_c2r_2d(intn, intn, pval, pofd,
				    fftw_plan_style);
    if (plan_inv == nullptr) {
      std::stringstream str;
      str << "Plan creation failed for inverse transform of size: " << intn
	  << " by " << intn;
      if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			  << " that size";
      throw pofdExcept("PDFactoryDouble", "initPD", str.str(), 15);
    }
    plans_valid = true;
  }

  //This will cause R wrapping problems, so check maxn0 relative to
  // the model base n0 value
  base_n0 = model.getBaseN0();
  if (maxn0 < base_n0) {
    std::stringstream errstr;
    errstr << "maxn0 (" << maxn0 << ") must be greater than model base N0 ("
	   << base_n0 << ")";
    throw pofdExcept("PDFactoryDouble", "initPD", errstr.str(), 8);
  }

  // Figure out the min/max nonzero R values in each band.
  std::pair<dblpair, dblpair> minmaxR = getMinMaxR(model, bm);

  double n0ratio = maxn0 / base_n0; 

  //Estimate the mean and standard deviation of the resulting P(D) using R.
  //Since the usage model for pofd_coverage is to call initPD once and then
  // re-use it a bunch of times, we can afford to compute R twice to get
  // a better estimate.  Note we compute these for the maximum n0
  std::pair<dblpair, dblpair> moments;
  moments = getRMoments(n, model, bm, minmaxR, setEdge);
  mn1 = n0ratio * moments.first.first;
  varnoi1 = n0ratio * moments.first.second;
  sg1 = sqrt(varnoi1 + inst_sigma1 * inst_sigma1);
  mn2 = n0ratio * moments.second.first;
  varnoi2 = n0ratio * moments.second.second;
  sg2 = sqrt(varnoi2 + inst_sigma2 * inst_sigma2);

  // Now that we have that estimate, figure out what range we will
  // as R to be computed over for the actual computation we will use
  // to form P(D).  This is rather messy, as it happens.
  // We want to ensure there is enough padding at the top.
  double maxfluxR_1, maxfluxR_2;
  maxfluxR_1 = maxflux1 + pofd_coverage::n_zero_pad * sg1;
  maxfluxR_2 = maxflux2 + pofd_coverage::n_zero_pad * sg2;

  //Get final R for base model, multiplying by dflux1 * dflux2.  
  initR(n, minfluxR_1, maxfluxR_1, minfluxR_2, maxfluxR_2, 
	model, bm, setEdge, true);
  
  //Decide if we will shift, and if so by how much
  // The idea is to shift the mean to zero -- but we only
  // do the shift if the step is larger than one actual step size
  // because otherwise we can't represent it well.
  doshift1 = fabs(mn1) >= dflux1;
  doshift2 = fabs(mn2) >= dflux2;
  if (doshift1) shift1 = - mn1; else shift1 = 0.0;
  if (doshift2) shift2 = - mn2; else shift2 = 0.0;

  if (verbose) {
    std::cout << " For max_n0: " << maxn0 << std::endl;
    std::cout << "  Initial mean estimate band1: " << mn1 << " band2: "
	      << mn2 << std::endl;
    std::cout << "  Initial stdev estimate band1: " << sg1 << " band2: "
	      << sg2 << std::endl;
    if (doshift1)
      std::cout << "  Additional shift applied in band 1: " 
		<< shift1 << std::endl;
    else std::cout << " Not applying additional shift in band 1" << std::endl;
    if (doshift2)
      std::cout << "  Additional shift applied in band 2: " 
		<< shift2 << std::endl;
    else std::cout << "  Not applying additional shift in band 2" << std::endl;
  }

  //Compute forward transform of this r value, store in rtrans
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute_dft_r2c(plan, rvals, rtrans);
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif
  
  max_n0 = maxn0;
  sigma1 = inst_sigma1;
  sigma2 = inst_sigma2;
  initialized = true;
}

/*!
  Calculates P(D) for a dataset based on already computed R (by initPD)
 
  \param[in] n0 Desired number of sources per sq deg
  \param[out] pd Holds P(D1, D2) on output
  \param[in] setLog If true, pd is log(P(D1,D2)) on output; convenient
              for likelihood evaluation.

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux often won't quite match the target values.

  You must call initPD first.
*/
void PDFactoryDouble::getPD(double n0, PDDouble& pd, bool setLog) {

  // The basic idea is to compute the P(D) from the previously filled
  // R values, adding in noise and all that fun stuff, filling pd
  // for output

  if (! initialized )
    throw pofdExcept("PDFactoryDouble","getPD",
		     "Must call initPD first",1);
  if (n0 > max_n0) {
    std::stringstream errstr("");
    errstr << "N_0 " << n0
	   << " larger than maximum prepared value " << max_n0
	   << std::endl;
    errstr << "initPD should have been called with at least " << n0;
    throw pofdExcept("PDFactoryDouble","getPD",errstr.str(),2);
  }
  double n0ratio = n0 / base_n0;

  //Output array from 2D FFT is n * (n/2+1)
  unsigned int n = currsize;
  unsigned int ncplx = n/2 + 1;
      
  //Calculate p(omega) = exp( r(omega1,omega2) - r(0,0) ),
  // which is what we will transform back into pofd.
  // There are some complications because of shifts and all that.
  //The 2D output real FFT format makes this a bit tricky.  The output
  // array is n by (n/2+1).  The first dimension has both
  // positive and negative frequencies, the second dimension has
  // only positive dimensions.
  //In particular, if i is the index over the first dimension,
  // the frequencies along the first dimension are:
  //  f1 = i/delta1*n        for i = [0,n/2]
  //     = - (n-i)/delta1*n  for i = [n/2+1,n-1]
  // where delta1 is dflux1.
  // And if j is the second dimension index, then
  //  f2 = j/delta2*n  
  // and delta2 = dflux2.  We work in w instead of f (2 pi f)
  // and actually compute 
  //  exp( r(omega1,omega2) - r(0,0) - i*shift1*omega1 - i*shift2*omega2
  //       - 1/2 sigma1^2 omega1^2 - 1/2 sigma2^2 * omega2^2)
  
  double r0 = n0ratio * rtrans[0][0]; //r[0,0] is pure real
  double iflux1 = pofd_coverage::two_pi / (n * dflux1);
  double iflux2 = pofd_coverage::two_pi / (n * dflux2);
  
  if (verbose) std::cout << "  Computing p(w1,w2)" << std::endl;

#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif

  // It's possible to compute p more efficiently if there are no shifts,
  // but for simplicity just act as if there are anyways.
  // Note that the shifts were computed for max_n0, not the current
  // or base n0
  double curr_shift1 = n0 * shift1 / max_n0;
  double curr_shift2 = n0 * shift2 / max_n0;
  double sigfac1 = 0.5 * sigma1 * sigma1;
  double sigfac2 = 0.5 * sigma2 * sigma2;  
  
  //First, Pos freq
  double expfac, rval, ival;
  fftw_complex *row_current_out; //Output out variable
  fftw_complex *r_input_rowptr; //Row pointer into out_part (i.e., input)
  double sigprod1, meanprod1, didx2, w1, w2;
  for (unsigned int idx1 = 0; idx1 < ncplx; ++idx1) {
    r_input_rowptr = rtrans + idx1 * ncplx; //Input
    row_current_out = pval + idx1 * ncplx; //Output
    w1 = iflux1 * static_cast<double>(idx1);
    meanprod1 = curr_shift1 * w1;
    sigprod1  = sigfac1 * w1 * w1;
    for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
      didx2 = static_cast<double>(idx2);
      w2 = iflux2 * didx2;
      rval = n0ratio * r_input_rowptr[idx2][0] - r0 
	- sigprod1 - sigfac2 * w2 * w2;
      ival = n0ratio * r_input_rowptr[idx2][1] - meanprod1 - curr_shift2 * w2;
      expfac = exp(rval);
      row_current_out[idx2][0] = expfac*cos(ival);
      row_current_out[idx2][1] = expfac*sin(ival);
    }
  }
  //Now, Neg freq
  for (unsigned int idx1 = ncplx; idx1 < n; ++idx1) {
    r_input_rowptr = rtrans + idx1 * ncplx; //Input
    row_current_out = pval + idx1 * ncplx; //Output
    w1 = - iflux1 * static_cast<double>(n - idx1);
    meanprod1 = curr_shift1 * w1;
    sigprod1  = sigfac1 * w1 * w1;
    for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
      didx2 = static_cast<double>(idx2);
      w2 = iflux2 * didx2;
      rval = n0ratio * r_input_rowptr[idx2][0] - 
	r0 - sigprod1 - sigfac2 * w2 * w2;
      ival = n0ratio * r_input_rowptr[idx2][1] - meanprod1 - curr_shift2 * w2;
      expfac = exp(rval);
      row_current_out[idx2][0] = expfac*cos(ival);
      row_current_out[idx2][1] = expfac*sin(ival);
    }
  }

  //p(0,0) is special
  pval[0][0] = 1.0;
  pval[0][1] = 0.0;

#ifdef TIMING
  p0Time += std::clock() - starttime;
#endif

  
  //Now transform back
  if (verbose) std::cout << " Reverse transform" << std::endl;
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute(plan_inv); //overwrites pofd with reverse transform of pval
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif

  // Copy into output variable, also normalizing, mean subtracting, 
  // making positive
  unwrapPD(n0, n, pd);

  //Turn PD to log for more efficient log computation of likelihood
#ifdef TIMING
  starttime = std::clock();
#endif
  if (setLog) pd.applyLog(false);
#ifdef TIMING
  logTime += std::clock() - starttime;
#endif
}


/*!
  Writes out current R
 
  \param[in] filename File to write to

  You must call initPD first, or bad things will probably happen.
*/
void PDFactoryDouble::writeRToFile(const std::string& filename) const {
  if (!rinitialized)
    throw pofdExcept("PDFactoryDouble", "writeRToFile",
		     "Must initialize R first", 1);

  std::ofstream ofs(filename.c_str());
  if (!ofs)
    throw pofdExcept("PDFactoryDouble", "writeRToFile",
		     "Couldn't open output file", 2);

  ofs << currsize << std::endl;
  if (rdflux) {
    double ifluxfac = 1.0 / (dflux1 * dflux2);
    for (unsigned int i = 0; i < currsize; ++i)
      for (unsigned int j = 0; j < currsize; ++j)
	ofs << RFlux1[i] << " " << RFlux2[j] << " " << 
	  rvals[i * currsize + j] * ifluxfac << std::endl;
  } else {
    for (unsigned int i = 0; i < currsize; ++i)
      for (unsigned int j = 0; j < currsize; ++j)
	ofs << RFlux1[i] << " " << RFlux2[j] << " " << 
	  rvals[i * currsize + j] << std::endl;
  }
  ofs.close();
}

/*!
  Writes out current R to a HDF5 file
 
  \param[in] filename File to write to

  You must call initR or initPD first, or bad things will probably happen.
*/
void PDFactoryDouble::writeRToHDF5(const std::string& filename) const {
  if (!rinitialized)
    throw pofdExcept("PDFactoryDouble", "writeRToHDF5",
		     "Must call initPD or initR first", 1);

  hid_t file_id;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);

  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw pofdExcept("PDFactoryDouble", "writeToHDF5",
		     "Failed to open HDF5 file to write", 2);
  }

  // Write it as one dataset -- Rflux1, Rflux2, R. 
  hsize_t adims;
  hid_t mems_id, att_id, dat_id;
  
  // First, some properties
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  att_id = H5Acreate2(file_id, "dflux1", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux1);
  H5Aclose(att_id);
  att_id = H5Acreate2(file_id, "dflux2", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux2);
  H5Aclose(att_id);
  att_id = H5Acreate2(file_id, "N0", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &base_n0);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  // Rfluxes
  adims = currsize;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  dat_id = H5Dcreate2(file_id, "RFlux1", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, RFlux1);
  H5Dclose(dat_id);
  dat_id = H5Dcreate2(file_id, "RFlux2", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, RFlux2);
  H5Dclose(dat_id);
  H5Sclose(mems_id);

  // R -- which we may need to copy to remove the dflux
  hsize_t dims_steps[2] = {currsize, currsize};
  mems_id = H5Screate_simple(2, dims_steps, nullptr);
  dat_id = H5Dcreate2(file_id, "R", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (rdflux) {
    double* tmp = new double[currsize * currsize];
    double idflux = 1.0 / (dflux1 * dflux2);
    for (unsigned int i = 0; i < currsize * currsize; ++i) 
      tmp[i] = rvals[i] * idflux;
    H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,  H5P_DEFAULT, tmp);
    delete[] tmp;
  } else
    H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,  H5P_DEFAULT, rvals);
  H5Dclose(dat_id);
  H5Sclose(mems_id);

  // Done
  H5Fclose(file_id);
}
