#include<iostream>

#include<getopt.h>
#include<hdf5.h>

#include "../include/global_settings.h"
#include "../include/utility.h"
#include "../include/beam.h"
#include "../include/doublebeam.h"
#include "../include/numberCounts.h"
#include "../include/numberCountsDouble.h"
#include "../include/pofdExcept.h"

//All sub-parses have to have the same options to avoid
// getting warnings.  So we give them all the same long_options,
// but then only process the appropriate ones, ignoring the rest
//Yes, this is a bit complicated and error prone, but such is life
static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"double", no_argument, 0, 'd'},
  {"filterscale", required_argument, 0, 'F'},
  {"matched", no_argument, 0, 'm'},
  {"nbins", required_argument, 0, 'n'},
  {"nfwhm", required_argument, 0, 'N'},
  {"nkeep", required_argument, 0, '1'},
  {"oversamp", required_argument, 0, 'o'},
  {"qfactor", required_argument, 0, 'q'},
  {"sigc", required_argument, 0, '3'},
  {"sigi", required_argument, 0, '2'},
  {"sigi1", required_argument, 0, '4'},
  {"sigi2", required_argument, 0, '5'},
  {"verbose", no_argument, 0, 'v'},
  {"version", no_argument, 0, 'V'}, 
  {0,0,0,0}
};

char optstring[] = "hdF:mn:N:1:o:q:2:3:4:5:vV";

//One-D version
int getRSingle(int argc, char** argv) {

  std::string modelfile; //Init file (having model we want)
  std::string outfile; //File to write to
  bool verbose, matched;
  double n0, nfwhm, pixsize, minflux, maxflux, fwhm, nkeep;
  unsigned int nflux, nbins, oversamp;
  double filterscale, qfactor, sigi, sigc;

  // Defaults
  verbose = false;
  nfwhm = 40.0;
  nbins = 120;
  filterscale = 0.0;
  qfactor = 0.2;
  matched = false;
  sigi = 0.002;
  sigc = 0.006;
  oversamp = 1;
  nkeep = 0; // Means keep all

  int c;
  int option_index = 0;
  optind = 1; //!< Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			    &option_index)) != -1) 
    switch(c) {
    case 'F':
      filterscale = atof(optarg);
      break;
    case 'm':
      matched = true;
      break;
    case 'n':
      nbins = atoi(optarg);
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case '1':
      nkeep = atof(optarg);
      break;
    case 'o':
      oversamp = atoi(optarg);
      break;
    case 'q':
      qfactor = atof(optarg);
      break;
    case '2':
      sigi = atof(optarg);
      break;
    case '3':
      sigc = atof(optarg);
      break;
    case 'v':
      verbose = true;
      break;
    }

  if (optind >= argc - 7) {
    std::cerr << "Required arguments missing; run with --help" << std::endl;
    return 1;
  }
  modelfile = std::string(argv[optind]);
  n0 = atof(argv[optind+1]);
  fwhm = atof(argv[optind+2]);
  pixsize = atof(argv[optind+3]);
  minflux = atof(argv[optind+4]);
  maxflux = atof(argv[optind+5]);
  nflux = static_cast<unsigned int>(atoi(argv[optind+6]));
  outfile = std::string(argv[optind+7]);

  if (nflux == 0) {
    std::cout << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (nbins == 0) {
    std::cout << "Invalid (non-positive) number of beam histogram bins"
	      << std::endl;
    return 1;
  }
  if (fwhm <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM" << std::endl;
    return 1;
  }
  if (nfwhm <= 0) {
    std::cout << "Invalid (non-positive) number of beam FWHMs"
	      << std::endl;
    return 1;
  }
  if (oversamp % 2 == 0) {
    std::cout << "Invalid (non-odd) oversampling " << oversamp << std::endl;
    return 1;
  }
  if (maxflux <= 0.0) {
    std::cout << "Invalid (non-positive) maxflux " << maxflux << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (non-positive) n0 " << n0 << std::endl;
    return 1;
  }
  if (nkeep < 0) {
    std::cout << "Invalid (negative) nkeep " << nkeep << std::endl;
    return 1;
  }
  if (filterscale < 0.0) {
    std::cout << "Invalid (negative) filter scale " << filterscale << std::endl;
    return 1;
  }
  if (qfactor < 0.0) {
    std::cout << "Invalid (negative) high-pass filter q factor" 
	      << qfactor << std::endl;
    return 1;
  }
  if (matched) {
    if (sigi <= 0.0) {
      std::cout << "Invalid (non-positive) sigi for matched filter"
		<< std::endl;
      return 1;
    }
    if (sigc <= 0.0) {
      std::cout << "Invalid (non-positive) sigc for matched filter"
		<< std::endl;
      return 1;
    }
  }

  double *R = NULL;
  double *flux = NULL;
  try {
    numberCounts model(modelfile);
    beam bm(fwhm);
    beamHist inv_bmhist(nbins);

    // Set up filter
    fourierFilter *filt = NULL;
    if (filterscale > 0) {
      if (matched) {
	filt = new fourierFilter(pixsize, fwhm, sigi, sigc,
				 filterscale, qfactor, true, true);
      } else
	filt = new fourierFilter(pixsize, filterscale, qfactor, true, true);
    } else if (matched)
      filt = new fourierFilter(pixsize, fwhm, sigi, sigc, true, true);

    // Fill beam
    inv_bmhist.fill(bm, nfwhm, pixsize, true, oversamp, filt, nkeep);

    if (filt != NULL) delete filt;

    if (n0 == 0)
      n0 = model.getBaseN0();

    if (verbose) {
      printf("   Beam fwhm:          %0.2f\n", bm.getFWHM());
      printf("   Beam area:          %0.3e\n", bm.getEffectiveArea());
      printf("   Pixel size:         %0.2f\n", pixsize);
      printf("   Flux per area:      %0.2f\n",
	     model.getBaseFluxPerArea());
      printf("   Base N0:            %0.4e\n", model.getBaseN0());
      printf("   N0:                 %0.4e\n", n0);
      if (filterscale > 0.0) {
	printf("   filter scale:       %0.1f\n", filterscale);
	printf("   filter q:           %0.2f\n", qfactor);
      }
      if (matched) {
	printf("   matched fwhm:       %0.1f\n", fwhm);
	printf("   matched sigi:       %0.4f\n", sigi);
	printf("   matched sigc:       %0.4f\n", sigc);
      }
      if (oversamp != 1)
	printf("   oversamp:           %u\n", oversamp);
    }

    // Get R
    double dflux;
    if (nflux > 1) dflux = (maxflux - minflux) / static_cast<double>(nflux - 1);
    else dflux = 0.0;
    flux = new double[nflux];
    for (unsigned int i = 0; i < nflux; ++i)
      flux[i] = static_cast<double>(i) * dflux + minflux;
    R = new double[nflux];

    model.getR(nflux, flux, inv_bmhist, R);
    
    // Adjust for N0
    double n0fac = n0 / model.getBaseN0();
    if (n0fac != 1)
      for (unsigned int i = 0; i < nflux; ++i)
	R[i] *= n0fac;

    // Write.
    utility::outfiletype oft = utility::getOutputFileType(outfile);
    if (oft == utility::HDF5 || oft == utility::UNKNOWN) {
      if (verbose) std::cout << "Writing as HDF5" << std::endl;
      hid_t file_id;
      file_id = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			  H5P_DEFAULT);
      if (H5Iget_ref(file_id) < 0) {
	H5Fclose(file_id);
	throw pofdExcept("pofd_coverge_getR", "pofd_coverage_getR",
			 "Failed to open HDF5 file to write", 1);
      }
      hsize_t adims;
      hid_t mems_id, att_id, dat_id;
      
      // Properties
      adims = 1;
      mems_id = H5Screate_simple(1, &adims, NULL);
      att_id = H5Acreate2(file_id, "dflux", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux);
      H5Aclose(att_id);
      att_id = H5Acreate2(file_id, "N0", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(att_id, H5T_NATIVE_DOUBLE, &n0);
      H5Aclose(att_id);
      H5Sclose(mems_id);

      // Rflux
      adims = nflux;
      mems_id = H5Screate_simple(1, &adims, NULL);
      dat_id = H5Dcreate2(file_id, "RFlux", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, flux);
      H5Dclose(dat_id);

      dat_id = H5Dcreate2(file_id, "R", H5T_NATIVE_DOUBLE, mems_id,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, R);
      H5Dclose(dat_id);
      H5Sclose(mems_id);

      H5Fclose(file_id);
    } else if (oft == utility::TXT) {
      if (verbose) std::cout << "Writing as text" << std::endl;
      FILE *fp;
      fp = fopen(outfile.c_str(), "w");
      if (!fp) {
	std::cerr << "Failed to open output file" << std::endl;
	return 128;
      }
      fprintf(fp, "#%-11s   %-12s\n", "Flux", "R");
      for (unsigned int i = 0; i < nflux; ++i) 
	fprintf(fp, "%12.6e   %15.9e\n", flux[i], R[i]);
      fclose(fp);
    } else if (oft == utility::FITS) {
      std::cerr << "Output to FITS is not supported." << std::endl;
      return 256;
    }

    delete[] flux;    
    delete[] R;
  } catch (const pofdExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (R != NULL) delete[] R;
    if (flux != NULL) delete[] flux;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (flux != NULL) delete[] flux;
    if (R != NULL) delete[] R;
    return 16;
  }
  return 0;
}

//Two-D version
int getRDouble(int argc, char** argv) {

  std::string modelfile; //Init file (having model we want)
  std::string outfile; //File to write to
  bool verbose, matched;
  double minflux1, maxflux1, minflux2, maxflux2;
  double n0, nfwhm, pixsize, fwhm1, fwhm2, nkeep;
  unsigned int nflux1, nflux2, nbins, oversamp;
  double filterscale, qfactor, sigi1, sigi2, sigc;

  // Defaults
  verbose = false;
  nfwhm = 40.0;
  nbins = 150;
  filterscale = 0.0;
  qfactor = 0.2;
  matched = false;
  sigc = 0.006;
  sigi1 = 0.002;
  sigi2 = 0.002;
  oversamp = 1;
  nkeep = 0.0;

  int c;
  int option_index = 0;
  optind = 1; //!< Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			    &option_index)) != -1) 
    switch(c) {
    case 'F':
      filterscale = atof(optarg);
      break;
    case 'm':
      matched = true;
      break;
    case 'n':
      nbins = atoi(optarg);
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case '1':
      nkeep = atof(optarg);
      break;
    case 'o':
      oversamp = atoi(optarg);
      break;
    case 'q':
      qfactor = atof(optarg);
      break;
    case '3':
      sigc = atof(optarg);
      break;
    case '4':
      sigi1 = atof(optarg);
      break;
    case '5':
      sigi2 = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc - 11) {
    std::cerr << "Required arguments missing; run with --help" << std::endl;
    return 1;
  }
  modelfile = std::string(argv[optind]);
  n0 = atof(argv[optind + 1]);
  fwhm1 = atof(argv[optind + 2]);
  fwhm2 = atof(argv[optind + 3]);
  pixsize = atof(argv[optind + 4]);
  minflux1 = atof(argv[optind + 5]);
  maxflux1 = atof(argv[optind + 6]);
  nflux1 = static_cast<unsigned int>(atoi(argv[optind + 7]));
  minflux2 = atof(argv[optind + 8]);
  maxflux2 = atof(argv[optind + 9]);
  nflux2 = static_cast<unsigned int>(atoi(argv[optind + 10]));
  outfile = std::string(argv[optind + 11]);

  if (nflux1 == 0) {
    std::cout << "Error -- number of fluxes requested is zero in dim 1."
	      << std::endl;
    return 1;
  }
  if (nflux2 == 0) {
    std::cout << "Error -- number of fluxes requested is zero in dim 2."
	      << std::endl;
    return 1;
  }
  if (nbins == 0) {
    std::cout << "Invalid (non-positive) number of beam histogram bins"
	      << std::endl;
    return 1;
  }
  if (fwhm1 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM1" << fwhm1 << std::endl;
    return 1;
  }
  if (fwhm2 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM2" << fwhm2 << std::endl;
    return 1;
  }
  if (nfwhm <= 0) {
    std::cout << "Invalid (non-positive) number of beam FWHMs"
	      << std::endl;
    return 1;
  }
  if (oversamp % 2 == 0) {
    std::cout << "Invalid (non-odd) oversampling " << oversamp << std::endl;
    return 1;
  }
  if (minflux1 > maxflux1) std::swap(minflux1, maxflux1);
  if (maxflux1 <= 0.0) {
    std::cout << "Invalid (non-positive) maxflux1 " << maxflux1 << std::endl;
    return 1;
  }
  if (minflux2 > maxflux2) std::swap(minflux1, maxflux1);
  if (maxflux2 <= 0.0) {
    std::cout << "Invalid (non-positive) maxflux1 " << maxflux2 << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (non-positive) n0 " << n0 << std::endl;
    return 1;
  }
  if (nkeep < 0) {
    std::cout << "Invalid (negative) nkeep " << nkeep << std::endl;
    return 1;
  }
  if (filterscale < 0.0) {
    std::cout << "Invalid (negative) filter scale " << filterscale << std::endl;
    return 1;
  }
  if (qfactor < 0.0) {
    std::cout << "Invalid (negative) high-pass filter q factor" 
	      << qfactor << std::endl;
    return 1;
  }
  if (matched) {
    if (sigi1 <= 0.0) {
      std::cout << "Invalid (non-positive) sigi1 for matched filter"
		<< std::endl;
      return 1;
    }
    if (sigi2 <= 0.0) {
      std::cout << "Invalid (non-positive) sigi2 for matched filter"
		<< std::endl;
      return 1;
    }
    if (sigc <= 0.0) {
      std::cout << "Invalid (non-positive) sigc for matched filter"
		<< std::endl;
      return 1;
    }
  }

  double *R = NULL;
  double *flux1 = NULL;
  double *flux2 = NULL;
  try {
    numberCountsDouble model(modelfile);
    doublebeam bm(fwhm1, fwhm2);
    doublebeamHist inv_bmhist(nbins);

    // Set up filter
    fourierFilter *filt1 = NULL, *filt2 = NULL;
    if (filterscale > 0) {
      if (matched) {
	filt1 = new fourierFilter(pixsize, fwhm1, sigi1, sigc,
				  filterscale, qfactor, true, true);
	filt2 = new fourierFilter(pixsize, fwhm2, sigi2, sigc,
				  filterscale, qfactor, true, true);
      } else
	filt1 = new fourierFilter(pixsize, filterscale, 0.1, true);
    } else if (matched) {
      filt1 = new fourierFilter(pixsize, fwhm1, sigi1, sigc, true);
      filt2 = new fourierFilter(pixsize, fwhm2, sigi2, sigc, true);
    }
    
    // Fill beam
    inv_bmhist.fill(bm, nfwhm, pixsize, true, oversamp, filt1, filt2, nkeep);

    if (filt1 != NULL) delete filt1;
    if (filt2 != NULL) delete filt2;

    if (n0 == 0)
      n0 = model.getBaseN0();

    if (verbose) {
      std::pair<double, double> dpr;
      dpr = bm.getFWHM();
      printf("   Beam fwhm1:         %0.2f\n", dpr.first);
      printf("   Beam fwhm2:         %0.2f\n", dpr.second);
      dpr = bm.getEffectiveArea();
      printf("   Beam area1:         %0.3e\n", dpr.first);
      printf("   Beam area2:         %0.3e\n", dpr.second);
      printf("   Pixel size:         %0.2f\n", pixsize);
      printf("   Flux per area1:     %0.2f\n",
	     model.getBaseFluxPerArea1());
      printf("   Flux per area2:     %0.2f\n",
	     model.getBaseFluxPerArea2());
      printf("   Base N0:            %0.4e\n", model.getBaseN0());
      printf("   N0:                 %0.4e\n", n0);
      if (filterscale > 0.0) {
	printf("   filter scale:       %0.1f\n", filterscale);
	printf("   filter q:           %0.2f\n", qfactor);
      }
      if (matched) {
	printf("   matched fwhm1:      %0.1f\n", fwhm1);
	printf("   matched fwhm2:      %0.1f\n", fwhm2);
	printf("   matched sigi1:      %0.4f\n", sigi1);
	printf("   matched sigi2:      %0.4f\n", sigi2);
	printf("   matched sigc:       %0.4f\n", sigc);
      }
      if (oversamp != 1)
	printf("   oversamp:           %u\n", oversamp);
    }

    // Set up fluxes
    double dflux1, dflux2;
    if (nflux1 > 1) 
      dflux1 = (maxflux1 - minflux1) / static_cast<double>(nflux1 - 1);
    else
      dflux1 = 0.0;
    if (nflux2 > 1)
      dflux2 = (maxflux2 - minflux2)/ static_cast<double>(nflux2 - 1);
    else
      dflux2 = 0.0;

    flux1 = new double[nflux1];
    for (unsigned int i = 0; i < nflux1; ++i)
      flux1[i] = dflux1 * static_cast<double>(i) + minflux1;
    flux2 = new double[nflux2];
    for (unsigned int i = 0; i < nflux2; ++i)
      flux2[i] = dflux2 * static_cast<double>(i) + minflux2;

    // Get R
    R = new double[nflux1 * nflux2];

    model.getR(nflux1, flux1, nflux2, flux2, inv_bmhist, R);

    // Adjust for N0
    double n0fac = n0 / model.getBaseN0();
    if (fabs(n0fac - 1.0) > 1e-5)
      for (unsigned int i = 0; i < nflux1 * nflux2; ++i)
	R[i] *= n0fac;

    // Write
    utility::outfiletype oft = utility::getOutputFileType(outfile);
    if (oft == utility::HDF5 || oft == utility::UNKNOWN) {
      if (verbose) std::cout << "Writing as HDF5" << std::endl;
      hid_t file_id;
      file_id = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			  H5P_DEFAULT);
      if (H5Iget_ref(file_id) < 0) {
	H5Fclose(file_id);
	throw pofdExcept("pofd_coverge_getR", "pofd_coverage_getR",
			 "Failed to open HDF5 file to write", 1);
      }
      hsize_t adims;
      hid_t mems_id, att_id, dat_id;
      
      // Properties
      adims = 1;
      mems_id = H5Screate_simple(1, &adims, NULL);
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
      H5Awrite(att_id, H5T_NATIVE_DOUBLE, &n0);
      H5Aclose(att_id);
      H5Sclose(mems_id);

      // Rfluxes
      adims = nflux1;
      mems_id = H5Screate_simple(1, &adims, NULL);
      dat_id = H5Dcreate2(file_id, "RFlux1", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, flux1);
      H5Dclose(dat_id);
      H5Sclose(mems_id);
      adims = nflux2;
      mems_id = H5Screate_simple(1, &adims, NULL);
      dat_id = H5Dcreate2(file_id, "RFlux2", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, flux2);
      H5Dclose(dat_id);
      H5Sclose(mems_id);

      // R, which is 2D
      hsize_t dims_steps[2] = {nflux1, nflux2};
      mems_id = H5Screate_simple(2, dims_steps, NULL);
      dat_id = H5Dcreate2(file_id, "R", H5T_NATIVE_DOUBLE, mems_id,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, R);			  
      H5Dclose(dat_id);
      H5Sclose(mems_id);

      H5Fclose(file_id);
    } else if (oft == utility::TXT) {
      if (verbose) std::cout << "Writing as text" << std::endl;
      FILE *fp;
      fp = fopen(outfile.c_str(), "w");
      if (!fp) {
	std::cerr << "Failed to open output file" << std::endl;
	return 128;
      }
      fprintf(fp,"#%4u %4u\n", nflux1, nflux2);
      fprintf(fp,"#minflux1: %12.6e dflux1: %12.6e\n", 0., dflux1);
      fprintf(fp,"#minflux2: %12.6e dflux2: %12.6e\n", 0., dflux2);
      for (unsigned int i = 0; i < nflux1; ++i) {
	for (unsigned int j = 0; j < nflux2 - 1; ++j)
	  fprintf(fp,"%13.7e ",R[nflux2 * i + j]);
	fprintf(fp,"%13.7e\n",R[nflux2 * i + nflux2 - 1]);
      }
      fclose(fp);
    } else if (oft == utility::FITS) {
      std::cerr << "Writing to FITS is not supported" << std::endl;
      return 256;
    }

    delete[] R;
    delete[] flux1;
    delete[] flux2;
  } catch (const pofdExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (R != NULL) delete[] R;
    if (flux1 != NULL) delete[] flux1;
    if (flux2 != NULL) delete[] flux2;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (R != NULL) delete[] R;
    if (flux1 != NULL) delete[] flux1;
    if (flux2 != NULL) delete[] flux2;
    return 16;
  }
  return 0;
}


////////////////////////////////////////

int main(int argc, char** argv) {
  bool twod;

  twod = false;

  //Only interested in a) displaying help and b) figuring out
  // if this is 1D or 2D c) displaying the version number
  int c;
  int option_index = 0;
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'h' :
      std::cout << "NAME" << std::endl;
      std::cout << "\tpofd_coverage_getR -- get R for a number counts model."
		<< "  Both" << std::endl;
      std::cout << "\tone-dimensional and two-dimensional models are supported."
		<< std::endl;
      std::cout << std::endl;
      std::cout << "SYNOPSIS" << std::endl;
      std::cout << "\tEither" << std::endl;
      std::cout << std::endl;
      std::cout << "\t pofd_coverage_getR [options] modelfile n0 fwhm pixsize"
		<< " minflux maxflux" << std::endl;
      std::cout << "\t\tnflux outfile" << std::endl;
      std::cout << std::endl;
      std::cout << "\tfor the 1D case or" << std::endl;
      std::cout << std::endl;
      std::cout << "\t pofd_coverage_getR -d [options] modelfile n0 fwhm1 fwhm2 "
		<< "pixsize" << std::endl;
      std::cout << "\t\tminflux1 maxflux1 nflux1 minflux2 maxflux2 nflux2"
		<< " outfile" << std::endl;
      std::cout << std::endl;
      std::cout << "\tfor the 2D case." << std::endl;
      std::cout << std::endl;
      std::cout << "DESCRIPTION" << std::endl;
      std::cout << "\tEvaluates R for the model in initfile using the P(D)"
		<< " formalism and" << std::endl;
      std::cout << "\twrites it to outfile.  The 1D model is a log-space "
		<< "spline" << std::endl;
      std::cout << "\tmodel for the number counts, and the 2D model is the 1D" 
		<< " spline" << std::endl;
      std::cout << "\tmodel times a log-normal color function for the second"
		<< " band," << std::endl;
      std::cout << "\twith the log-space variance and mean color stored as"
		<< " splines" << std::endl;
      std::cout << "\tin the flux of the first band." << std::endl;
      std::cout << std::endl;
      std::cout << "\tmodelfile is a text file specifying the model; the exact"
		<< " details" << std::endl;
      std::cout << "\t(given below) depend on whether the 1D or 2D case is"
		<< " being used." << std::endl;
      std::cout << std::endl;
      std::cout << "\tFor the 1D case, modelfile is a text file giving the "
		<< "positions" << std::endl;
      std::cout << "\tof the spline knots and their values in the format"
		<< " knotflux value." << std::endl;
      std::cout << "\tAdditional elements on each line are ignored."
		<< std::endl;
      std::cout << "\tFor the 2D case, modelfile is a text file giving the "
		<< "positions" << std::endl;
      std::cout << "\tof the knot points and their values, followed by the "
		<< "sigma" << std::endl;
      std::cout << "\tknot positions and their values, then likewise for the "
		<< "colour" << std::endl;
      std::cout << "\toffset.  The format is three numbers on the first line, "
		<< "giving" << std::endl;
      std::cout << "\tthe number of number count knots, sigma knots, and "
		<< "offset knots," << std::endl;
      std::cout << "\tfollowed by a number of lines again with the format"
		<< std::endl;
      std::cout << "\tknotpos value.  The sigmas and offsets are in log space."
		<< std::endl;
      std::cout << std::endl;
      std::cout << "\tn0 is the number of sources per square degree.  If set"
		<< " to" << std::endl;
      std::cout << "\tzero, then the value from the input model is used "
		<< "directly." << std::endl;
      std::cout << "\tOtherwise, the model is scaled to match this value."
		<< std::endl;
      std::cout << std::endl;
      std::cout << "\tfwhm is the beam FWHM (in arcsec) in the 1D case, and"
		<< std::endl;
      std::cout << "\tfwhm1, fwhm2 are the FWHM values for each band in the 2D"
		<< "case." << std::endl;
      std::cout << "\tpixscale is the pixel scale (in arcsec)." << std::endl;
      std::cout << std::endl;
      std::cout << "\tminflux, maxflux give the minimum and maximum"
		<< " flux density" << std::endl;
      std::cout << "\tto evaluate R for in the 1D case.  For the 2D case this" 
		<< std::endl;
      std::cout << "\tis extended to the minimum and maximum in each band."
		<< " nflux" << std::endl;
      std::cout << "\tis the number of fluxes to generate; in the 2D case along"
		<< " each" << std::endl;
      std::cout << "\tdimension." << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "\tThe format of the output (text, hdf5) is set by"
		<< " the" << std::endl;
      std::cout << "\textension of outfile.  The default is hdf5." << std::endl;
      std::cout << std::endl;
      std::cout << "OPTIONS" << std::endl;
      std::cout << "\t-d, --double" << std::endl;
      std::cout << "\t\tUse the 2D model." << std::endl;
      std::cout << "\t-F, --filterscale VALUE" << std::endl;
      std::cout << "\t\tHigh-pass filter scale, in arcsec.  Zero means no"
		<< " filtering" << std::endl;
      std::cout << "\t\tis applied (def: 0)" << std::endl;
      std::cout << "\t-m, --matched" << std::endl;
      std::cout << "\t\tApply matched filtering to the beam, with a FWHM"
		<< " matching the" << std::endl;
      std::cout << "\t\tbeam (each band in the 2D case).  Off by default." 
		<< std::endl;
      std::cout << "\t-n, --nbins VALUE" << std::endl;
      std::cout << "\t\tNumber of beam histogram bins (def: 120)" << std::endl;
      std::cout << "\t-N, --nfwhm VALUE" << std::endl;
      std::cout << "\t\tNumber of FWHM to go out in beam representation. "
		<< "(def: 40.0)" << std::endl;
      std::cout << "\t--nkeep VALUE" << std::endl;
      std::cout << "\t\tNumber of FWHM out to keep after filtering in beam"
		<< std::endl;
      std::cout << "\t\trepresentation.  The default is to keep all of it."
		<< std::endl;
      std::cout << "\t-o, --oversamp VALUE" << std::endl;
      std::cout << "\t\tOversampling of pixels used to generate beam. One means"
		<< std::endl;
      std::cout << "\t\tno oversampling.  Must be odd (def: 1)" << std::endl;
      std::cout << "\t-q, --qfactor VALUE" << std::endl;
      std::cout << "\t\tHigh-pass filter apodization sigma as fraction of"
		<< std::endl;
      std::cout << "\t\tfiltscale. (def: 0.2)." << std::endl;
      std::cout << "\t--sigc VALUE" << std::endl;
      std::cout << "\t\tConfusion noise for matched filtering, in Jy. (Def:"
		<< " 0.006)" << std::endl;
      std::cout << "\t-V, --version" << std::endl;
      std::cout << "\t\tOutput the version number and exit." << std::endl;
      std::cout << "ONE-DIMENSIONAL OPTIONS" << std::endl;
      std::cout << "\t--sigi VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, in Jy. (Def:"
		<< " 0.002)" << std::endl;
      std::cout << "TWO-DIMENSIONAL OPTIONS" << std::endl;
      std::cout << "\t--sigi1 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 1, in Jy."
		<< std::endl;
      std::cout << "\t\t(Def: 0.002)" << std::endl;
      std::cout << "\t--sigi2 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 2, in Jy."
		<< std::endl;
      std::cout << "\t\t(Def: 0.002)" << std::endl;
      std::cout << std::endl;
      return 0;
      break;
    case 'd' :
      twod = true;
      break;
    case 'V' :
      std::cout << "pofd_coverage version number: " << pofd_coverage::version
		<< std::endl;
      return 0;
      break;
    }

  if (! twod)
    return getRSingle(argc,argv);
  else
    return getRDouble(argc,argv);
}
