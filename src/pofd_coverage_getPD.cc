#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cstdlib>

#include<getopt.h>

#include "../include/global_settings.h"
#include "../include/utility.h"
#include "../include/PDFactory.h"
#include "../include/beam.h"
#include "../include/numberCounts.h"
#include "../include/PDFactoryDouble.h"
#include "../include/doublebeam.h"
#include "../include/numberCountsDouble.h"
#include "../include/pofdExcept.h"

static struct option long_options[] = {
  {"double", no_argument, 0, 'd'},
  {"edgeinterp", required_argument, 0, 'e'},
  {"help", no_argument, 0, 'h'},
  {"filterscale", required_argument, 0, 'F'},
  {"log", no_argument, 0, 'l'},
  {"matched", no_argument, 0, 'm'},
  {"nflux", required_argument, 0, 'n'},
  {"nbins", required_argument, 0, '0'},
  {"nfwhm", required_argument, 0, 'N'},
  {"nkeep", required_argument, 0, '1'},
  {"oversamp", required_argument, 0, 'o'},
  {"qfactor", required_argument, 0, 'q'},
  {"rfile", required_argument, 0, 'r'},
  {"sigma", required_argument, 0, 's'},
  {"sigma1", required_argument, 0, '3'},
  {"sigma2", required_argument, 0, '4'},
  {"sigc", required_argument, 0, '6'},
  {"sigc1", required_argument, 0, '!'},
  {"sigc2", required_argument, 0, '@'},
  {"sigi", required_argument, 0, '5'},
  {"sigi1", required_argument, 0, '7'},
  {"sigi2", required_argument, 0, '8'},
  {"verbose", no_argument, 0, 'v'},
  {"version", no_argument, 0, 'V'},
  {"wisdom", required_argument, 0, 'w'},
  {0, 0, 0, 0}
};
char optstring[] = "dhe:F:lmn:N:1:0:o:q:r:s:3:4:5:6:!:@:7:8:vVw:";

///////////////////////////////

int getPDSingle(int argc, char **argv) {

  std::string modelfile; //Knot parameters
  double n0; //Number of sources
  double maxflux, nkeep, fwhm, nfwhm, pixsize;
  double sigma; //Instrument noise
  double filterscale, qfactor, sigi, sigc; // Filtering parameters
  std::string outputfile; //Ouput pofd option
  unsigned int nflux, nbins;
  bool has_wisdom, verbose, return_log, matched, write_r;
  std::string wisdom_file, r_file;
  unsigned int oversample;

  //Defaults
  sigma               = 0.0;
  has_wisdom          = false;
  nflux               = 131072;
  nbins               = 120;
  nfwhm               = 4.5;
  nkeep               = 0.0;
  verbose             = false;
  return_log          = false;
  write_r             = false;
  oversample          = 1;
  filterscale         = 0.0;
  qfactor             = 0.2;
  matched             = false;
  sigi                = 0.0; // Means: use sigma
  sigc                = 0.006;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'F' :
      filterscale = atof(optarg);
      break;
    case 'l' :
      return_log = true;
      break;
    case 'm':
      matched = true;
      break;
    case 'n' :
      nflux = static_cast<unsigned int>(atoi(optarg));
      break;
    case '0':
      nbins = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case '1':
      nkeep = atof(optarg);
      break;
    case 'o':
      oversample = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'q':
      qfactor = atof(optarg);
      break;
    case 'r':
      write_r = true;
      r_file = std::string(optarg);
      break;
    case 's' :
      sigma = atof(optarg);
      break;
    case '5':
      sigi = atof(optarg);
      break;
    case '6':
      sigc = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    case 'w' :
      has_wisdom = true;
      wisdom_file = std::string(optarg);
      break;
    }

  if (optind >= argc - 5) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string(argv[optind]);
  n0         = atof(argv[optind + 1]);
  fwhm       = atof(argv[optind + 2]);
  pixsize    = atof(argv[optind + 3]);
  maxflux    = atof(argv[optind + 4]);
  outputfile = std::string(argv[optind + 5]);

  if (matched && (sigi == 0)) sigi = sigma;

  //Input tests
  if (nflux == 0) {
    std::cout << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (sigma < 0.0) {
    std::cout << "Invalid noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (negative) number of sources per area"
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
  if (pixsize >= fwhm / 2.0) {
    std::cout << "Insufficient (FWHM/2) beam sampling based on pixel size"
	      << std::endl;
    return 1;
  }
  if (oversample % 2 == 0) {
    std::cout << "Invalid (non-odd) oversampling " << oversample << std::endl;
    return 1;
  }
  if (nkeep < 0.0) {
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

  try {
    numberCounts model(modelfile);
    beam bm(fwhm);
    PDFactory pfactory;
    PD pd;

    if (verbose) pfactory.setVerbose();

    bool success;
    if (has_wisdom) {
      if (verbose)
	std::cout << "Reading in wisdom file: " << wisdom_file 
		  << std::endl;
      success = pfactory.addWisdom(wisdom_file);
      if (!success) {
	std::cout << "Error reading wisdom file: " << wisdom_file << std::endl;
	return 4;
      }
    }

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
      printf("   sigma:              %0.4f\n", sigma);
      if (filterscale > 0.0) {
	printf("   filter scale:       %0.1f\n", filterscale);
	printf("   filter q:           %0.2f\n", qfactor);
      }
      if (matched) {
	printf("   matched fwhm:       %0.1f\n", fwhm);
	printf("   matched sigi:       %0.4f\n", sigi);
	printf("   matched sigc:       %0.4f\n", sigc);
      }
      if (oversample != 1)
	printf("   oversamp:           %u\n", oversample);
      if (return_log) 
	printf("   Returning log(P(D)) rather than P(D)\n");
    }

    fourierFilter *filt = NULL;
    if (filterscale > 0) {
      if (matched) {
	filt = new fourierFilter(pixsize, fwhm, sigi, sigc,
				 filterscale, qfactor, true);
      } else
	filt = new fourierFilter(pixsize, filterscale, qfactor, true);
    } else if (matched)
      filt = new fourierFilter(pixsize, fwhm, sigi, sigc, true);


    // Get histogrammed inverse beam
    beamHist inv_bmhist(nbins);
    inv_bmhist.fill(bm, nfwhm, pixsize, true, oversample, filt, nkeep);

    if (filt != NULL) delete filt;

    //Get P(D)
    if (verbose) std::cout << "Getting P(D) with transform length: " 
			   << nflux << " and max flux: " 
			   << maxflux << std::endl;
    double base_n0 = model.getBaseN0();
    pfactory.initPD(nflux, sigma, maxflux, base_n0 > n0 ? base_n0 : n0, 
		    model, inv_bmhist);

    if (write_r) {
      if (verbose) std::cout << "  Writing R to " << r_file;
      utility::outfiletype oftype = utility::getOutputFileType(r_file);
      switch (oftype) {
      case utility::UNKNOWN:
      case utility::HDF5:
	pfactory.writeRToHDF5(r_file);
	break;
      case utility::FITS:
	throw pofdExcept("pofd_coverage_getPD", "getPDSingle",
			 "Don't support writing R to FITS", 1);
	break;
      case utility::TXT:
	pfactory.writeRToFile(r_file);
	break;
      }
    }
   
    pfactory.getPD(n0, pd, return_log);
    
    //Write it
    if (verbose) std::cout << "Writing P(D) to " << outputfile << std::endl;
    pd.writeToFile(outputfile);

#ifdef TIMING
    std::cout << "Timing results:" << std::endl;
    pfactory.summarizeTime(2);
#endif
  } catch ( const pofdExcept& ex ) {
    std::cout << "Error encountered" << std::endl;
    std::cout << ex << std::endl;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cout << "Bad allocation error: " << ba.what() << std::endl;
    return 16;
  } 

  return 0;
}

///////////////////////////////////////

int getPDDouble(int argc, char **argv) {

  std::string modelfile; //Knot parameters
  double n0; //Number of sources
  double maxflux1, maxflux2; //Maximum fluxes requested
  double fwhm1, fwhm2, nfwhm, pixsize, nkeep;
  double filterscale, qfactor; //Highpass filter params
  double sigi, sigi1, sigi2, sigc, sigc1, sigc2; // Matched filter params
  double sigma1, sigma2; //Instrument noise
  std::string outputfile; //Ouput pofd option
  unsigned int nflux, nbins, oversample;
  bool has_wisdom, verbose, return_log, matched, write_r, single_filt;
  std::string wisdom_file, r_file;

  //Defaults
  sigma1              = 0.0;
  sigma2              = 0.0;
  has_wisdom          = false;
  nflux               = 2048;
  nbins               = 150;
  nfwhm               = 4.5;
  nkeep               = 0.0;
  filterscale         = 0.0;
  qfactor             = 0.2;
  matched             = false;
  single_filt         = true;
  sigc                = 0.006;
  sigc1               = 0.0; // Use sigc if not set
  sigc2               = 0.0; // Use sigc if not set
  sigi                = 0.002;
  sigi1               = 0.0; // Use sigi if not set
  sigi2               = 0.0; // Use sigi if not set
  verbose             = false;
  return_log          = false;
  write_r             = false;
  oversample          = 1;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'F':
      filterscale = atof(optarg);
      break;
    case 'l':
      return_log = true;
      break;
    case 'm':
      matched = true;
      break;
    case 'n':
      nflux = static_cast<unsigned int>(atoi(optarg));
      break;
    case '0':
      nbins = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case '1':
      nkeep = atof(optarg);
      break;
    case 'o':
      oversample = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'q':
      qfactor = atof(optarg);
      break;
    case 'r':
      write_r = true;
      r_file = std::string(optarg);
      break;
    case '3':
      sigma1 = atof(optarg);
      break;
    case '4':
      sigma2 = atof(optarg);
      break;
    case '6':
      sigc = atof(optarg);
      break;
    case '!':
      sigc1 = atof(optarg);
      single_filt = false;
      break;
    case '@':
      sigc2 = atof(optarg);
      single_filt = false;
      break;
    case '5':
      sigi = atof(optarg);
      break;
    case '7':
      sigi1 = atof(optarg);
      single_filt = false;
      break;
    case '8':
      sigi2 = atof(optarg);
      single_filt = false;
      break;
    case 'v':
      verbose = true;
      break;
    case 'w':
      has_wisdom = true;
      wisdom_file = std::string(optarg);
      break;
    }

  if (optind >= argc - 7) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string(argv[optind]);
  n0         = atof(argv[optind + 1]);
  fwhm1      = atof(argv[optind + 2]);
  fwhm2      = atof(argv[optind + 3]);
  pixsize    = atof(argv[optind + 4]);
  maxflux1   = atof(argv[optind + 5]);
  maxflux2   = atof(argv[optind + 6]);
  outputfile = std::string(argv[optind + 7]);

  //Input tests
  if (nflux == 0) {
    std::cout << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (sigma1 < 0.0) {
    std::cout << "Invalid noise level (band1): must be >= 0.0" << std::endl;
    return 1;
  }
  if (sigma2 < 0.0) {
    std::cout << "Invalid noise level (band2): must be >= 0.0" << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (negative) number of sources per area"
	      << std::endl;
    return 1;
  }
  if (nbins == 0) {
    std::cout << "Invalid (non-positive) number of beam histogram bins"
	      << std::endl;
    return 1;
  }
  if (fwhm1 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM (band 1)" << std::endl;
    return 1;
  }
  if (fwhm2 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM (band 2)" << std::endl;
    return 1;
  }
  if (nfwhm <= 0) {
    std::cout << "Invalid (non-positive) number of beam FWHMs"
	      << std::endl;
    return 1;
  }
  if (pixsize >= fwhm1 / 2.0 || pixsize >= fwhm2 / 2.0) {
    std::cout << "Insufficient (FWHM/2) beam sampling based on pixel size"
	      << std::endl;
    return 1;
  }
  if (oversample % 2 == 0) {
    std::cout << "Invalid (non-odd) oversampling " << oversample << std::endl;
    return 1;
  }
  if (nkeep < 0.0) {
    std::cout << "Invalid (negative) nkeep " << nkeep << std::endl;
    return 1;
  }
  if (filterscale < 0.0) {
    std::cout << "Invalid (negative) filter scale: " << filterscale
	      << std::endl;
    return 1;
  }
  if (qfactor < 0.0) {
    std::cout << "Invalid (negative) high-pass filter q factor" 
	      << qfactor << std::endl;
    return 1;
  }

  if (matched) {
    if (single_filt) {
      if (sigc <= 0.0) {
	std::cout << "Invalid sigma_confusion for single matched filter: "
		  << sigc << std::endl;
	return 1;
      }
      if (sigi <= 0.0) {
	std::cout << "Invalid sigma_instrument for single matched filter: "
		  << sigi << std::endl;
	return 1;
      }
    } else {
      // Different filters for each band.  More complex
      if (sigi1 == 0) sigi1 = sigi;
      if (sigi2 == 0) sigi2 = sigi;
      if (sigc1 == 0) sigc1 = sigc;
      if (sigc2 == 0) sigc2 = sigc;
      if (sigc1 <= 0.0) {
	std::cout << "Invalid sigma_confusion1 for double matched filters: "
		  << sigc1 << std::endl;
	return 1;
      }
      if (sigc2 <= 0.0) {
	std::cout << "Invalid sigma_confusion2 for double matched filters: "
		  << sigc2 << std::endl;
	return 1;
      }
      if (sigi1 <= 0.0) {
	std::cout << "Invalid sigma_instrument1 for double matched filter: "
		  << sigi1 << std::endl;
	return 1;
      }
      if (sigi2 <= 0.0) {
	std::cout << "Invalid sigma_instrument2 for double matched filter: "
		  << sigi2 << std::endl;
	return 1;
      }
    }
  }

  try {
    numberCountsDouble model(modelfile);
    doublebeam bm(fwhm1, fwhm2);
    PDFactoryDouble pfactory;
    PDDouble pd;

    if (verbose) pfactory.setVerbose();

    bool success;
    if (has_wisdom) {
      if (verbose)
	std::cout << "Reading in wisdom file: " << wisdom_file 
		  << std::endl;
      success = pfactory.addWisdom(wisdom_file);
      if (!success) {
	std::cout << "Error reading wisdom file: " << wisdom_file << std::endl;
	return 4;
      }
    }
    
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
      printf("   sigma1:             %0.4f\n", sigma1);
      printf("   sigma2:             %0.4f\n", sigma2);
      if (filterscale > 0.0) {
	printf("   filter scale:       %0.1f\n", filterscale);
	printf("   filter q:           %0.2f\n", qfactor);
      }
      if (matched) {
	dpr = bm.getFWHM();
	if (single_filt) {
	  printf("   matched fwhm:       %0.1f\n", dpr.first);
	  printf("   matched sigi:       %0.4f\n", sigi);
	  printf("   matched sigc:       %0.4f\n", sigc);
	} else {
	  printf("   matched fwhm1:      %0.1f\n", dpr.first);
	  printf("   matched fwhm2:      %0.1f\n", dpr.second);
	  printf("   matched sigi1:      %0.4f\n", sigi1);
	  printf("   matched sigi2:      %0.4f\n", sigi2);
	  printf("   matched sigc1:      %0.4f\n", sigc1);
	  printf("   matched sigc2:      %0.4f\n", sigc2);
	}
      }
      if (oversample != 1)
	printf("  oversamp:            %u\n", oversample);
      if (return_log) 
	printf("  Returning log(P(D)) rather than P(D)\n");
    }

    // Set up filter
    fourierFilter *filt1 = NULL, *filt2 = NULL;
    if (filterscale > 0) {
      if (matched) {
	// Both highpass and matched
	if (single_filt) {
	  filt1 = new fourierFilter(pixsize, fwhm1, sigi, sigc,
				    filterscale, qfactor, true);
	} else {
	  filt1 = new fourierFilter(pixsize, fwhm1, sigi1, sigc1,
				    filterscale, qfactor, true);
	  filt2 = new fourierFilter(pixsize, fwhm2, sigi2, sigc2,
				    filterscale, qfactor, true);
	}
      } else // Just highpass, can always use one filter
	filt1 = new fourierFilter(pixsize, filterscale, qfactor, true);
    } else if (matched) {
      if (single_filt) {
	filt1 = new fourierFilter(pixsize, fwhm1, sigi, sigc, true);
      } else {
	filt1 = new fourierFilter(pixsize, fwhm1, sigi1, sigc1, true);
	filt2 = new fourierFilter(pixsize, fwhm2, sigi2, sigc2, true);
      }
    }

    // Get histogrammed inverse beam
    doublebeamHist inv_bmhist(nbins);
    inv_bmhist.fill(bm, nfwhm, pixsize, true, oversample, filt1, filt2, nkeep);

    if (filt1 != NULL) delete filt1; // Don't need these any more
    if (filt2 != NULL) delete filt2; 

    //Get P(D)
    if (verbose) std::cout << "Getting P(D) with transform length: " 
			   << nflux << " and max fluxes: " 
			   << maxflux1 << " " << maxflux2 << std::endl;
    double base_n0 = model.getBaseN0();
    pfactory.initPD(nflux, sigma1, sigma2, maxflux1, maxflux2, 
		    base_n0 > n0 ? base_n0 : n0, model, inv_bmhist, true);

    if (write_r) {
      if (verbose) std::cout << "  Writing R to " << r_file;
      utility::outfiletype oftype = utility::getOutputFileType(r_file);
      switch (oftype) {
      case utility::UNKNOWN:
      case utility::HDF5:
	pfactory.writeRToHDF5(r_file);
	break;
      case utility::FITS:
	throw pofdExcept("pofd_coverage_getPD", "getPDouble",
			 "Don't support writing R to FITS", 1);
	break;
      case utility::TXT:
	pfactory.writeRToFile(r_file);
	break;
      }
    }
    pfactory.getPD(n0, pd, return_log);
    
    //Write it
    if (verbose) std::cout << "Writing P(D) to " << outputfile << std::endl;
    pd.writeToFile(outputfile);

#ifdef TIMING
    std::cout << "Timing results:" << std::endl;
    pfactory.summarizeTime(2);
#endif
  } catch ( const pofdExcept& ex ) {
    std::cout << "Error encountered" << std::endl;
    std::cout << ex << std::endl;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cout << "Bad allocation error: " << ba.what() << std::endl;
    return 16;
  } 

  return 0;
}


////////////////////////////////////
int main( int argc, char** argv ) {
  bool twod;

  twod = false;

  //Only interested in a) displaying help and b) figuring out
  // if this is 1D or 2D c) displaying the version number
  int c;
  int option_index = 0;
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cout << "NAME" << std::endl;
      std::cout << "\tpofd_coverage_getPD -- get the P(D) for a spline"
		<< std::endl;
      std::cout << "\t type model with a Gaussian beam (1D) or the same type"
		<< " of model" << std::endl;
      std::cout << "\t paired with a log-normal color model in flux_2/flux_1."
		<< std::endl;
      std::cout << std::endl;
      std::cout << "SYNOPSIS" << std::endl;
      std::cout << "\t One-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_getPD [options] modelfile n0 fwhm pixsize"
		<< std::endl;
      std::cout << "\t                      maxflux outfile" << std::endl;
      std::cout << std::endl;
      std::cout << "\t Two-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_getPD [options] -d modelfile n0 fwhm1"
		<< " fwhm2" << std::endl;
      std::cout << "\t                      pixsize maxflux1 maxflux2 outfile" 
		<< std::endl; 
      std::cout << std::endl;
      std::cout << "DESCRIPTION" << std::endl;
      std::cout << "\tEvaluates P(D) for the specified model and writes it to" 
		<< std::endl;
      std::cout << "\toutfile.  The number counts model in 1D is a spline"
		<< std::endl;
      std::cout << "\tmodel specified by modelfile, and by the number of"
		<< std::endl;
      std::cout << "\tsources per unit area n0." << std::endl;
      std::cout << std::endl;
      std::cout << "\tIn the 2D case the model is the 1D model in band 1 times"
		<< " a" << std::endl;
      std::cout << "\tLog-Normal distribution in flux2/flux1.  The mu and sigma"
		<< " Log-Normal" << std::endl;
      std::cout << "\tmodel parameters are stored as splines as a function of"
		<< " the" << std::endl;
      std::cout << "\tflux in the first band." << std::endl;
      std::cout << std::endl;
      std::cout << "\tmodelfile should be a text file.  For 1D it consists of"
		<< " nknots" << std::endl;
      std::cout << "\tlines of the form:" << std::endl << std::endl;
      std::cout << "\t\tflux_density n" << std::endl << std::endl;
      std::cout << "\twhere flux_density gives the positions of the knots"
		<< " in Jy" << std::endl;
      std::cout << "\tand n is the log10 differential number counts in"
		<< " deg^-2 Jy^-1" << std::endl;
      std::cout << "\tat the corresponding flux density.  Additional entries" 
		<< " on each"<< std::endl;
      std::cout << "\tline are ignored, and # denotes a comment line."
		<< std::endl;
      std::cout << std::endl;
      std::cout << "\tIn the 2D case the file should start with a line giving"
		<< " the" << std::endl;
      std::cout << "\tnumber of knots in the band 1 model, the number of"
		<< " knots in" << std::endl;
      std::cout << "\tthe sigma spline, and then the number in the mu spline."
		<< " This" << std::endl;
      std::cout << "\tshould be followed by nknots + nspline + nmu lines"
		<< std::endl;
      std::cout << "\tof the same form as the 1D model, with the first nknots"
		<< std::endl;
      std::cout << "\tspecifying the band 1 model as in the 1D case, and the"
		<< std::endl;
      std::cout << "\tfollowing lines giving the knot positions and values"
		<< " for" << std::endl;
      std::cout << "\tof the sigma and mu splines." << std::endl;
      std::cout << std::endl;
      std::cout << "\tfwhm is the beam FWHM in arcsec.  The beam is assumed "
		<< "Gaussian." << std::endl;
      std::cout << "\tThe pixel size, in arcsec, is specified by pixsize." 
		<< std::endl;
      std::cout << "\tIn the 2D case, fwhm1 and fwhm2 are the values for each"
		<< " band." << std::endl;
      std::cout << std::endl;
      std::cout << "\tmaxflux is the maximum flux density generated.  In"
		<< " general" << std::endl;
      std::cout << "\tthe maxflux values will not quite be realized.  Again,"
		<< std::endl;
      std::cout << "\tin the 2D case maxflux1 and maxflux2 are the values in"
		<< " each" << std::endl;
      std::cout << "\tband.  Pixsize has the same meaning; a single size must" 
		<< std::endl;
      std::cout << "\tbe used for both bands." << std::endl;
      std::cout << std::endl;
      std::cout << "\tThe format of the output (text, fits, hdf5) is set by"
		<< " the" << std::endl;
      std::cout << "\textension of outfile." << std::endl;
      std::cout << std::endl;
      std::cout << "OPTIONS" << std::endl;
      std::cout << "\t-d, --twod" << std::endl;
      std::cout << "\t\tIf set, the two-dimensional model is used."
		<< std::endl;
      std::cout << "\t-h --help" << std::endl;
      std::cout << "\t\tPrint this message and exit." << std::endl;
      std::cout << "\t-F, --filtscale VALUE" << std::endl;
      std::cout << "\t\tRadius of high-pass filter in arcseconds. If zero,"
		<< std::endl;
      std::cout << "\t\tno filtering is applied (def: 0)." << std::endl;
      std::cout << "\t-l, --log" << std::endl;
      std::cout << "\t\tReturn the log P(D) rather than the P(D)."
		<< std::endl;
      std::cout << "\t-m, --matched" << std::endl;
      std::cout << "\t\tApply matched filtering to the beam, with a FWHM"
		<< " matching the" << std::endl;
      std::cout << "\t\tbeam (each beam in the 2d case)."
		<< std::endl;
      std::cout << "\t-n, --nflux value" << std::endl;
      std::cout << "\t\tThe number of requested fluxes along each dimension."
		<< std::endl;
      std::cout << "\t\tAlso sets the transform size used. (def: 131072 in 1D,)"
		<< std::endl;
      std::cout << "\t\tand 2048 in 2D)." << std::endl;
      std::cout << "\t--nbins value" << std::endl;
      std::cout << "\t\tNumber of bins to use in histogrammed beam. (def: 120)"
		<< std::endl;
      std::cout << "\t-N, --nfwhm value" << std::endl;
      std::cout << "\t\tNumber of beam FWHM out to go when computing beam."
		<< "(def: 3.5)" << std::endl;
      std::cerr << "\t--nkeep VALUE" << std::endl;
      std::cerr << "\t\tNumber of FWHM out to keep after filtering in beam"
		<< std::endl;
      std::cerr << "\t\trepresentation.  The default is to keep all of it."
		<< std::endl;
      std::cout << "\t-o, --oversample VALUE" << std::endl;
      std::cout << "\t\tAmount to oversample the beam; must be odd integer."
		<< " (def: 1)" << std::endl;
      std::cout << "\t-q, --qfactor VALUE" << std::endl;
      std::cout << "\t\tHigh-pass filter apodization sigma as fraction of"
		<< std::endl;
      std::cout << "\t\tfiltscale. (def: 0.2)." << std::endl;
      std::cout << "\t-r, --rfile FILENAME" << std::endl;
      std::cout << "\t\tWrite the R used to this file as text." << std::endl;
      std::cout << "\t--sigc VALUE" << std::endl;
      std::cout << "\t\tConfusion noise for matched filtering, in Jy. (Def:"
		<< " 0.006)" << std::endl;
      std::cout << "\t--sigi VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, in Jy"
		<< std::endl;
      std::cout << "\t\t(Def: 0.002)." << std::endl;
      std::cout << "\t-v, --verbose" << std::endl;
      std::cout << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cout << "\t-V, --version" << std::endl;
      std::cout << "\t\tOutput version number and exit" << std::endl;
      std::cout << "\t-w, --wisdom wisdomfile" << std::endl;
      std::cout << "\t\tName of wisdom file (prepared with fftw-wisdom)." 
		<< std::endl;
      std::cout << "ONE-D MODEL OPTIONS" << std::endl;
      std::cout << "\t-s, --sigma VALUE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise (def: 0)" << std::endl;
      std::cout << "TWO-D MODEL OPTIONS" << std::endl;
      std::cout << "\t--sigma1 NOISE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise, band 1 (def: 0)." 
		<< std::endl;
      std::cout << "\t--sigma2 NOISE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise, band 2 (def: 0)." 
		<< std::endl;
      std::cout << "\t--sigc1 VALUE" << std::endl;
      std::cout << "\t\tConfusion noise, band 1 for matched filtering, in Jy."
		<< std::endl;
      std::cout << "\t\t(Def: sigc)" << std::endl;
      std::cout << "\t--sigc2 VALUE" << std::endl;
      std::cout << "\t\tConfusion noise, band 2 for matched filtering, in Jy."
		<< std::endl;
      std::cout << "\t\t(Def: sigc)" << std::endl;
      std::cout << "\t--sigi1 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 1, in "
		<< std::endl;
      std::cout << "\t\tJy. (Def: sigi)" << std::endl;
      std::cout << "\t--sigi2 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 2, in "
		<< std::endl;
      std::cout << "\t\tJy. (Def: sigi)." << std::endl;
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

  if (!twod)
    return getPDSingle(argc, argv);
  else
    return getPDDouble(argc, argv);
}
