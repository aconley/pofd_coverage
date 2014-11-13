#include<fstream>
#include<sstream>
#include<limits>

#include<fitsio.h>

#include "../include/positionGenerator.h"
#include "../include/pofdExcept.h"
#include "../include/utility.h"

////////////////// Powerspectrum ////////////////
/*!
  \param[in] filename File to read k, P(k) from
*/
powerSpectrum::powerSpectrum(const std::string& filename) {
  // Just read the file, pass it to init
  std::ifstream ifs(filename.c_str());
  if (!ifs) {
    std::stringstream errstr;
    errstr << "Unable to open input model file: " << filename;
    throw pofdExcept("numberCounts", "numberCounts", errstr.str(), 1);
  }

  //Do the read into temporary vectors, then copy
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double currval;
  std::vector<double> kin, pkin;
  const unsigned int nreq = 2;
  while (!ifs.eof()) {
    std::getline(ifs, line);
    if (line[0] == '#' || line[0] == '%') continue; //Skip comments
    
    //Parse into words, stipping spaces
    utility::stringwords(line, words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#' or words[0][0] == '%') continue; //Comment line
    if (words.size() < nreq) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currval;
    kin.push_back(currval);
    str.str(words[1]); str.clear(); str >> currval;
    pkin.push_back(currval);
  }
  ifs.close();

  init(kin, pkin);
}

/*!
  \param[in] k k vector in inverse arcmin.  Assumed sorted
  \param[in] pk P(k)
*/
powerSpectrum::powerSpectrum(const std::vector<double>& k,
			     const std::vector<double>& pk) {
  init(k, pk);
}

powerSpectrum::~powerSpectrum() {
  delete[] logk;
  delete[] logpk;
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline_logpk);
}


/*!
  \param[in] k k vector in inverse arcmin.  Assumed sorted
  \param[in] pk P(k)
*/
void powerSpectrum::init(const std::vector<double>& k,
			 const std::vector<double>& pk) {
  nk = k.size();
  if (nk == 0)
    throw pofdExcept("powerSpectrum", "powerSpectrum",
		     "No elements in k", 1);
  if (nk != pk.size())
    throw pofdExcept("powerSpectrum", "powerSpectrum",
		     "P(k) not same length as k", 2);
  if (nk < 3)
    throw pofdExcept("powerSpectrum", "powerSpectrum",
		     "Need at least 3 elements in k, P(k)", 3);

  logk = new double[nk];
  mink = k[0]; // We assume k is sorted
  maxk = k[nk-1];
  if (mink <= 0.0) throw pofdExcept("powerSpectrum", "powerSpectrum",
				    "Invalid (non-positive) k", 4);

  for (unsigned int i = 0; i < nk; ++i)
    logk[i] = log(k[i]);

  logpk = new double[nk];
  for (unsigned int i = 0; i < nk; ++i)
    if (k[i] <= 0.0)
      throw pofdExcept("powerSpectrum", "powerSpectrum",
		       "Invalid (non-positive) P(k)", 5);
  for (unsigned int i = 0; i < nk; ++i)
    logpk[i] = log(pk[i]);

  acc = gsl_interp_accel_alloc();
  spline_logpk = gsl_spline_alloc(gsl_interp_cspline,
				  static_cast<size_t>(nk));
  gsl_spline_init(spline_logpk, logk, logpk, static_cast<size_t>(nk));
}

double powerSpectrum::getPk(double k) const {
  // Taking advantage of the fact that mink must be positive
  if (k < mink || k > maxk) return 0;
  return exp(gsl_spline_eval(spline_logpk, log(k), acc));
}

double powerSpectrum::getLogPk(double k) const {
  if (k < mink || k > maxk) std::numeric_limits<double>::quiet_NaN();
  return gsl_spline_eval(spline_logpk, log(k), acc);
}

/*!
  \param[in] objid HDF5 object ID to write to.  Must already be open
*/
void powerSpectrum::writeToHDF5Handle(hid_t obj_id) const {
  if (H5Iget_ref(obj_id) < 0)
    throw pofdExcept("powerSpectrum", "writeToHDF5Handle",
		     "Given non-open obj_id to write to", 1);

  hsize_t adims;
  hid_t mems_id, att_id, dat_id;

  // Single item attributes
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);

  att_id = H5Acreate2(obj_id, "Mink", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &mink);
  H5Aclose(att_id);  
  att_id = H5Acreate2(obj_id, "Maxk", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &maxk);
  H5Aclose(att_id);  
  H5Sclose(mems_id);

  // Knot positions and values as data
  double* tmp = new double[nk];
  adims = nk;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  dat_id = H5Dcreate2(obj_id, "k", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for (unsigned int i = 0; i < nk; ++i)
    tmp[i] = exp(logk[i]);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
	   H5P_DEFAULT, tmp);
  H5Dclose(dat_id);

  dat_id = H5Dcreate2(obj_id, "Pk", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for (unsigned int i = 0; i < nk; ++i)
    tmp[i] = exp(logpk[i]);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
	   H5P_DEFAULT, tmp);
  H5Dclose(dat_id);
  H5Sclose(mems_id);
  delete[] tmp;
}


/////////////////// positionGeneratorClustered ///////////////////
positionGeneratorClustered::
positionGeneratorClustered(unsigned int NX, unsigned int NY, 
			   double PIXSIZE, const std::string& powerfile) :
  nx(NX), ny(NY), pixsize(PIXSIZE), powspec(powerfile), probarr(nullptr),
  probarr_trans(nullptr) {

  if (nx == 0)
    throw pofdExcept("positionGeneratorClustered", 
		     "positionGeneratorClustered", 
		     "invalid (0) x dimension", 1);
  if (ny == 0)
    throw pofdExcept("positionGeneratorClustered", 
		     "positionGeneratorClustered", 
		     "invalid (0) y dimension", 1);

  if (pixsize <= 0.0)
    throw pofdExcept("positionGeneratorClustered", 
		     "positionGeneratorClustered",
		     "Invalid (non-positive) pixel size", 1);

  // Now, generate the scaling factor
  // First we generate the k values, and store them in scl
  double tail_invk = 60.0 / pixsize;
  unsigned int nkx = nx / 2 + 1;
  double kxfac = tail_invk / static_cast<double>(nx);
  unsigned int nky = ny / 2 + 1;
  double kyfac = tail_invk / static_cast<double>(ny);
  double kxfacsq = kxfac * kxfac;
  double kyfacsq = kyfac * kyfac;

  unsigned int xtmp, ytmp;
  double kxval, kyval;
  scl = new double[nx * ny];
  for (unsigned int i = 0; i < nkx; ++i) {
    kxval = kxfacsq * static_cast<double>(i * i);
    for (unsigned int j = 0; j < nky; ++j) {
      kyval = kyfacsq * static_cast<double>(j * j);
      scl[ny * i + j] = sqrt(kxval + kyval);
    }
    for (unsigned int j = nky; j < ny; ++j) {
      ytmp = ny - j;
      kyval = kyfacsq * static_cast<double>(ytmp * ytmp);
      scl[ny * i + j] = sqrt(kxval + kyval);
    }
  }
  for (unsigned int i = nkx; i < nx; ++i) {
    xtmp = nx - i;
    kxval = kxfacsq * static_cast<double>(xtmp * xtmp);
    for (unsigned int j = 0; j < nky; ++j) {
      kyval = kyfacsq * static_cast<double>(j * j);
      scl[ny * i + j] = sqrt(kxval + kyval);
    }
    for (unsigned int j = nky; j < ny; ++j) {
      ytmp = ny - j;
      kyval = kyfacsq * static_cast<double>(ytmp * ytmp);
      scl[ny * i + j] = sqrt(kxval + kyval);
    }
  }

  // Now build the scaling factor, which is sqrt(P(k))
  scl[0] = 0.0; //Mean 0
  for (unsigned int i = 1; i < nx * ny; ++i)
    scl[i] = sqrt(powspec.getPk(scl[i]));
}

positionGeneratorClustered::~positionGeneratorClustered() {
  delete[] scl;
  if (probarr != nullptr) fftw_free(probarr);
  if (probarr_trans != nullptr) fftw_free(probarr_trans);
  if (plan != nullptr) fftw_destroy_plan(plan); 
  if (plan_inv != nullptr) fftw_destroy_plan(plan_inv);
}

/*!
  \param[in] rangen Random number generator
*/
void positionGeneratorClustered::generate(ran& rangen) {

  if (probarr == nullptr) {
    probarr = (double*) fftw_malloc(sizeof(double) * nx * ny);
    nyhalf = ny / 2 + 1;
    probarr_trans = (fftw_complex*) 
      fftw_malloc(sizeof(fftw_complex) * nx * nyhalf);
    unsigned intx = static_cast<int>(nx);
    unsigned inty = static_cast<int>(ny);
    plan = fftw_plan_dft_r2c_2d(intx, inty, probarr, probarr_trans,
				FFTW_ESTIMATE);
    if (plan == nullptr) {
      std::stringstream str;
      str << "Plan creation failed for forward transform of size: " << 
	  intx << " by " << inty;
      throw pofdExcept("positionGeneratorClustered", "generate",
		       str.str(), 1);
    }
    plan_inv = fftw_plan_dft_c2r_2d(intx, inty, probarr_trans, probarr,
				    FFTW_ESTIMATE);
    if (plan_inv == nullptr) {
      std::stringstream str;
      str << "Inverse plan creation failed for forward transform of size: " << 
	  intx << " by " << inty;
      throw pofdExcept("positionGeneratorClustered", "generate",
		       str.str(), 2);
    }
  }

  // Fill
  for (unsigned int i = 0; i < nx * ny; ++i)
    probarr[i] = rangen.doub();

  // Transform
  fftw_execute(plan);

  // Multiply by the scaling factor
  double val;
  for (unsigned int i = 0; i < nx; ++i)
    for (unsigned int j = 0; j < nyhalf; ++j) {
      val = scl[i * ny + j];
      probarr_trans[i * nyhalf + j][0] *= val;
      probarr_trans[i * nyhalf + j][1] *= val;
    }

  // Transform back
  fftw_execute(plan_inv);

  // Renormalize
  double minval, maxval, currval;
  minval = maxval = probarr[0];
  for (unsigned int i = 0; i < nx * ny; ++i) {
    currval = probarr[i];
    if (currval > maxval) maxval = currval;
    if (currval < minval) minval = currval;
  }

  double rng = maxval - minval;
  if (rng == 0)
    for (unsigned int i = 0; i < nx * ny; ++i)
      probarr[i] -= minval - 1.0; //Make it all one
  else {
    double irng = 1.0 / rng;
    for (unsigned int i = 0; i < nx * ny; ++i)
      probarr[i] = (probarr[i] - minval) * irng; 
  }
    
}

/*!
  \param[inout] rangen Random number generator to use
  \returns Positions in the range [0, nx), [0, ny)

  This version does not have oversampling -- no interpolation
  is used.
*/
std::pair<unsigned int, unsigned int> 
positionGeneratorClustered::getPosition(ran& rangen) const {
  unsigned int idx1, idx2;
  double zval;

  do {
    idx1 = static_cast<unsigned int>(rangen.doub() * nx);
    idx2 = static_cast<unsigned int>(rangen.doub() * ny);
    zval = rangen.doub();
  } while (zval > probarr[idx1 * ny + idx2]);

  return std::make_pair(idx1, idx2);
}

/*!
  \param[in] x x position, in range [0, nx)
  \param[in] y y position, in range [0, ny)
  \returns Interpolated value of probarr

  Uses bilinear interpolation
*/
double positionGeneratorClustered::interpolate(double x, double y) const {
  // Check if out of bounds
  if (x < 0 || x >= (nx-1) || (y < 0) || y >= (ny-1)) return 0.0;

  unsigned int idx1 = static_cast<unsigned int>(x);
  unsigned int idx2 = static_cast<unsigned int>(y);
  const double *prob_base = probarr + idx1 * ny + idx2;

  double u, t, omu, omt;
  t = x - floor(x);
  u = y - floor(y);
  omu = 1.0 - u; omt = 1.0 - t;  
  return omt * (omu * prob_base[0] + u * prob_base[1]) +
    t * (omu * prob_base[ny] + u * prob_base[ny + 1]);  
}

/*!
  \param[inout] rangen Random number generator to use
  \param[in] oversamp Oversampling to use in position generation
  \returns Positions in the range [0, oversamp * nx), [0, oversamp * ny)

  Bilinear interpolation is used on the probability image to generate
  sub-positions.

  This gets called a lot, so the inputs are not checked for speed.
*/
std::pair<unsigned int, unsigned int> 
positionGeneratorClustered::getPosition(ran& rangen, 
					unsigned int oversamp) const {
  double x, y, zval;
  double onx = static_cast<double>(oversamp * nx);
  double ony = static_cast<double>(oversamp * ny);
  double idoversamp = 1.0 / static_cast<double>(oversamp);
  
  do {
    x = floor(onx * rangen.doub()); // [0, oversamp * nx), integral
    y = floor(ony * rangen.doub()); // [0, oversamp * ny), integral
    zval = rangen.doub();
  } while (zval > interpolate(x * idoversamp, y * idoversamp));

  return std::make_pair(static_cast<unsigned int>(x),
			static_cast<unsigned int>(y));
}
 
/*!
  \param[in] outfile Name of FITS file to write to
*/
int
positionGeneratorClustered::writeProbToFits(const std::string& outfile) const {

  if (probarr == nullptr)
    throw pofdExcept("positionGeneratorClustered", "writeProbToFits",
		     "probarr not prepared", 1);

  int status = 0;
  fitsfile *fp;

  fits_create_file(&fp, outfile.c_str(), &status);
  if (status) {
    fits_report_error(stderr, status);
    return status;
  }

  long axissize[2];
  axissize[0] = static_cast<long>(nx);
  axissize[1] = static_cast<long>(ny);
  fits_create_img(fp, DOUBLE_IMG, 2, axissize, &status);

  fits_write_history(fp, const_cast<char*>("Probability image"), &status);
  fits_write_date(fp, &status);

  // Astrometry
  double tmpval;
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE1"),
		 const_cast<char*>("RA---TAN"),
		 const_cast<char*>("WCS: Projection type axis 1"), &status);
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE2"),
		 const_cast<char*>("DEC--TAN"),
		 const_cast<char*>("WCS: Projection type axis 2"), &status);
  tmpval = nx / 2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRPIX1"), &tmpval, 
		 const_cast<char*>("Ref pix of axis 1"), &status);
  tmpval = ny / 2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRPIX2"), &tmpval, 
		 const_cast<char*>("Ref pix of axis 2"), &status);
  tmpval = 90.0; //Arbitrary
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRVAL1"), &tmpval, 
		 const_cast<char*>("val at ref pix axis 1"), &status);
  tmpval = 10.0; //Arbitrary
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRVAL2"), &tmpval, 
		 const_cast<char*>("val at ref pix axis 2"), &status);
  tmpval = - pixsize / 3600.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD1_1"), &tmpval, 
		 const_cast<char*>("Pixel scale axis 1,1"), &status);
  tmpval = 0.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD1_2"), &tmpval, 
		 const_cast<char*>("Pixel scale axis 1,2"), &status);
  tmpval = 0.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD2_1"), &tmpval, 
		 const_cast<char*>("Pixel scale axis 2,1"), &status);
  tmpval = pixsize / 3600.0;
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
  double *tmpdata = new double[nx];
  long fpixel[2] = { 1, 1 };
  for ( unsigned int j = 0; j < ny; ++j ) {
    for (unsigned int i = 0; i < nx; ++i) tmpdata[i] = probarr[i * ny + j];
    fpixel[1] = static_cast<long>(j + 1);
    fits_write_pix(fp, TDOUBLE, fpixel, nx, tmpdata, &status);
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
  \param[in] obj_id Open HDF5 handle to write to.

  Writes summary of power spectrum properties
*/
void positionGeneratorClustered::writeToHDF5Handle(hid_t obj_id) const {
  if (H5Iget_ref(obj_id) < 0)
    throw pofdExcept("positionGeneratorClustered", "writeToHDF5Handle",
		     "Given non-open obj_id to write to", 1);

  // Single item attributes
  hsize_t adims;
  hid_t mems_id, att_id;

  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);

  att_id = H5Acreate2(obj_id, "nx", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nx);
  H5Aclose(att_id);  
  att_id = H5Acreate2(obj_id, "ny", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &ny);
  H5Aclose(att_id);  
  att_id = H5Acreate2(obj_id, "pixsize", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &pixsize);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  powspec.writeToHDF5Handle(obj_id);
}
