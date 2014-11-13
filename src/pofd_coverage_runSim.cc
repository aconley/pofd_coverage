#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cstdlib>

#include<getopt.h>

#include "../include/simManager.h"
#include "../include/simManagerDouble.h"
#include "../include/pofdExcept.h"
#include "../include/global_settings.h"
#include "../include/utility.h"

static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"bin", no_argument, 0, 'b'},
  {"double", no_argument, 0, 'd'},
  {"esmooth", required_argument, 0, 'e'},
  {"extra_smooth1", required_argument, 0, '!'},
  {"extra_smooth2", required_argument, 0, '@'},
  {"fftsize", required_argument, 0, 'f'},
  {"filtscale", required_argument, 0, 'F'},
  {"matched", no_argument, 0, 'm'},
  {"nbeambins", required_argument, 0, '5'},
  {"nbins", required_argument, 0, '1'},
  {"nfwhm", required_argument, 0, '6'},
  {"nlike", required_argument, 0, '2'},
  {"nolike", no_argument, 0, 'N'},
  {"nsims", required_argument, 0, 'n'},
  {"n0initrange", required_argument, 0, '3'},
  {"n0rangefrac", required_argument, 0, '4'},
  {"oversample", required_argument, 0, 'o'},
  {"powspec", required_argument, 0, 'p'},
  {"qfactor", required_argument, 0, 'q'},
  {"sigma", required_argument, 0, 's'},
  {"sigma1", required_argument, 0, '#'},
  {"sigma2", required_argument, 0, '$'},
  {"sigc", required_argument, 0, ','},
  {"sigc1", required_argument, 0, '.'},
  {"sigc2", required_argument, 0, '/'},
  {"sigma_matched", required_argument, 0, '7'},
  {"sigma_matched1", required_argument, 0, '9'},
  {"sigma_matched2", required_argument, 0, '0'},
  {"sparcity", required_argument, 0, '%'},
  {"seed", required_argument, 0, 'S'},
  {"verbose", no_argument, 0, 'v'},
  {"version", no_argument, 0, 'V'},
  {"wisdom", required_argument, 0, 'w'},
  {0, 0, 0, 0}
};

char optstring[] = "hbde:!:@:f:F:1:mn:5:6:N2:3:4:o:p:q:s:#:$:,:.:/:7:9:0:%:S:vVw:";

int runSimSingle(int argc, char **argv) {

  std::string modelfile; //Base model file
  double n0; //Model params
  double fwhm; //Calculation params (req)
  double esmooth; //Extra smoothing amount
  double sigma; // Instrument noise
  unsigned int nsims, nlike, n1, n2, fftsize, nbins;
  unsigned int oversample, sparcity, nbeambins;
  double pixsize, n0rangefrac, n0initrange, nfwhm;
  double filtscale, qfactor, sigmc, sigmi; //Filtering parameters
  std::string outputfile; //Ouput pofd option
  std::string powerspecfile; //Power spectrum file
  bool verbose, has_wisdom, has_user_seed, use_binning, map_like, matched;
  unsigned long long int seed;
  std::string wisdom_file;

  //Defaults
  fftsize             = 131072;
  nsims               = 100;
  nlike               = 500;
  n0rangefrac         = 0.1;
  n0initrange         = 0.15;
  sigma               = 0.002;
  verbose             = false;
  has_wisdom          = false;
  esmooth             = 0.0;
  has_user_seed       = false;
  seed                = 1024;
  oversample          = 1;
  filtscale           = 0.0;
  qfactor             = 0.2;
  matched             = false;
  sigmi               = 0.0;  // Means: use sigma
  sigmc               = 0.006;
  sparcity            = 1;
  nbins               = 1000;
  use_binning         = false;
  map_like            = true;
  powerspecfile       = "";
  nfwhm               = 15.0;
  nbeambins           = 100;

  int c;
  int option_index = 0;
  optind = 1; //Rewind 
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'b':
      use_binning = true;
      break;
    case 'e':
      esmooth = atof(optarg);
      break;
    case 'f':
      fftsize = atoi(optarg);
      break;
    case 'F':
      filtscale = atof(optarg);
      break;
    case 'm':
      matched = true;
      break;
    case '5':
      nbeambins = atoi(optarg);
      break;
    case '1':
      nbins = atoi(optarg);
      break;
    case '6':
      nfwhm = atof(optarg);
      break;
    case 'n':
      nsims = atoi(optarg);
      break;
    case 'N':
      map_like = false;
      break;
    case '2':
      nlike = atoi(optarg);
      break;
    case '3':
      n0initrange = atof(optarg);
      break;
    case '4':
      n0rangefrac = atof(optarg);
      break;
    case 'o':
      oversample = atoi(optarg);
      break;
    case 'p':
      powerspecfile = std::string(optarg);
      break;
    case 'q':
      qfactor = atof(optarg);
      break;
    case 's':
      sigma = atof(optarg);
      break;
    case ',':
      sigmc = atof(optarg);
      break;
    case '7':
      sigmi = atof(optarg);
      break;
    case '%':
      sparcity = atoi(optarg);
      break;
    case 'S':
      has_user_seed = true;
      seed = static_cast<unsigned long long int>( atoi(optarg) );
      break;
    case 'v':
      verbose = true;
      break;
    case 'w':
      has_wisdom = true;
      wisdom_file = std::string( optarg );
      break;
    }

  if (optind >= argc-6) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  n0         = atof(argv[optind]);
  modelfile  = std::string(argv[optind + 1]);
  fwhm       = atof(argv[optind + 2]);
  pixsize    = atof(argv[optind + 3]);
  n1         = atoi(argv[optind + 4]);
  n2         = atoi(argv[optind + 5]);
  outputfile = std::string(argv[optind + 6]);

  if (matched && sigmi == 0) sigmi = sigma;

  if (sigma < 0.0) {
    std::cout << "Invalid instrument noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (n0initrange <= 0.0) {
    std::cout << "Invalid n0initrange: must be > 0" << std::endl;
    return 1;
  }
  if (n0initrange >= 1.0) {
    std::cout << "Invalid n0initrange: must be < 1" << std::endl;
    return 1;
  }
  if (fwhm <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM" << std::endl;
    return 1;
  }
  if (nfwhm <= 0.0) {
    std::cout << "Invalid (non-positive) NFWHM: " << nfwhm << std::endl;
    return 1;
  }
  if (nbeambins == 0) {
    std::cout << "Invalid number of beam bins (0)" << std::endl;
    return 1;
  }
  if (nbeambins > 10000) {
    std::cout << "Invalid number of beam bins -- too large: " 
	      << nbeambins << std::endl;
    return 1;
  }
  if (esmooth < 0.0) {
    std::cout << "Invalid (negative) additional smoothing" << std::endl;
    return 1;
  }
  if (use_binning && nbins == 0) {
    std::cout << "Invalid (zero) number of bins" << std::endl;
    return 1;
  }
  if (oversample == 0) {
    std::cout << "Invalid (non-positive) oversampling" << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (negative) n0" << std::endl;
    return 1;
  }
  if (pixsize <= 0.0) {
    std::cout << "Invalid (non-positve) pixsize" << std::endl;
    return 1;
  }
  if (n1 * n2 == 0) {
    std::cout << "Simulated map has zero size" << std::endl;
    return 1;
  }
  if (sparcity > n1 * n2) {
    std::cout << "Sampling sparcity less frequent than map size" << std::endl;
    return 1;
  }
  if (map_like && n0rangefrac <= 0.0) {
    std::cout << "Invalid n0rangefrac: must be > 0" << std::endl;
    return 1;
  }
  if (map_like && n0rangefrac >= 1.0) {
    std::cout << "Invalid n0rangefrac: must be < 1" << std::endl;
    return 1;
  }
  if (filtscale < 0.0) {
    std::cout << "Invalid (negative) filter scale: " << filtscale << std::endl;
    return 1;
  }
  if (qfactor < 0.0) {
    std::cout << "Invalid (negative) high-pass filter q factor" 
	      << qfactor << std::endl;
    return 1;
  }
  if (matched) {
    if (sigmi == 0) {
      std::cout << "Invalid (non-positive) sigma instrument for matched "
		<< "filtering" << std::endl;
      return 1;
    }
    if (sigmc <= 0.0) {
      std::cout << "Invalid (non-positive) sigma confusion for matched filter"
		<< std::endl;
      return 1;
    }
  }

  try {
    double base_n0;
    if (n0 == 0 || verbose) {
      // Use base model for n0
      numberCounts model(modelfile);
      base_n0 = model.getBaseN0();
      if (n0 == 0)
	n0 = base_n0;
    }

    if (verbose) {
      double area = n1*n2*std::pow(pixsize/3600.0,2);
      printf("   base model file:    %s\n", modelfile.c_str());
      printf("   base n0:            %0.3e\n", base_n0);
      printf("   nsims:              %u\n", nsims);
      printf("   n0initrange         %0.3f\n", n0initrange);
      printf("   Beam fwhm:          %0.2f\n", fwhm);
      printf("   Nfwhm:              %0.2f\n", nfwhm);
      printf("   N Beam Bins:        %u\n", nbeambins);
      printf("   Pixel size:         %0.1f\n", pixsize);
      printf("   Number of pixels:   %u x %u\n", n1, n2);
      printf("   Area:               %0.2f\n", area);
      printf("   N0:                 %0.3e\n", n0);
      printf("   sigma:              %0.4f\n", sigma);
      printf("   fftsize:            %u\n", fftsize);
      if (esmooth > 0) 
	printf("   esmooth:            %0.2f\n",esmooth);
      if (oversample > 1)
	printf("   oversampling:       %u\n", oversample);
      if (filtscale > 0) {
	printf("   filtering scale:    %0.1f\n", filtscale);
	printf("   filtering q:        %0.2f\n", qfactor);
      }
      if (matched > 0) {
	printf("   matched fwhm:       %0.1f\n", fwhm);
	printf("   matched sigi:       %0.4f\n", sigmi);
	printf("   matched sigc:       %0.4f\n", sigmc);
      }
      if (sparcity > 1)
	printf("   sparcity:           %u\n", sparcity);
      if (!powerspecfile.empty())
	printf("   clustering P(k):    %s\n", powerspecfile.c_str());
      if (use_binning)
	printf("   nbins:              %u\n", nbins);
      if (map_like) {
	printf("   nlike:              %u\n", nlike);
	printf("   n0rangefrac         %0.3f\n",n0rangefrac);
      }
    }

    simManager sim(modelfile, nsims, n0initrange, map_like, nlike, 
		   n0rangefrac, fftsize, n1, n2, pixsize, fwhm, nfwhm,
		   sigma, filtscale, qfactor, matched, sigmi, sigmc, nbeambins, 
		   n0, esmooth, oversample, powerspecfile, sparcity, 
		   use_binning, nbins);
    if (has_wisdom) sim.addWisdom(wisdom_file);
    if (has_user_seed) sim.setSeed(seed);

    sim.doSims(verbose);

    //Write it
    if (verbose) std::cout << "Writing simulation results to " << outputfile 
			   << std::endl;
    utility::outfiletype oft = utility::getOutputFileType(outputfile);
    if (oft == utility::HDF5 || oft == utility::UNKNOWN)
      sim.writeToHDF5(outputfile);
    else if (oft == utility::FITS) {
      int status = sim.writeToFits(outputfile);
      if (status != 0) return status;
    } else if (oft == utility::TXT)
      throw pofdExcept("pofd_coverage_runSim", "runSimSingle",
		       "Output to text not supported", 1);
  } catch (const pofdExcept& ex) {
    std::cout << "Error encountered" << std::endl;
    std::cout << ex << std::endl;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cout << "Bad allocation error: " << ba.what() << std::endl;
    return 16;
  } 

  return 0;
}

////////////////////////////////////////

int runSimDouble(int argc, char **argv) {

  std::string modelfile; //Base model file
  double n0; //Model params
  double fwhm1, fwhm2; //Calculation params (req)
  double esmooth1, esmooth2; //Extra smoothing amount
  double sigma1, sigma2; //Instrument noise, unsmoothed
  unsigned int nsims, nlike, n1, n2, fftsize, nbins, oversample;
  unsigned int nbeambins, sparcity;
  bool verbose, has_wisdom, has_user_seed, use_binning, map_like;
  double pixsize, n0rangefrac, n0initrange, nfwhm;

  // Filtering params
  bool matched, single_filt;
  double filterscale, qfactor; // Hipass
  double sigm1, sigm2, sigc1, sigc2, sigm, sigc; // Matched

  std::string outputfile; //Ouput pofd option
  std::string powerspecfile; // Power spectrum file

  unsigned long long int seed;
  std::string wisdom_file;

  //Defaults
  fftsize             = 4096;
  nsims               = 100;
  nlike               = 500;
  n0rangefrac         = 0.1;
  n0initrange         = 0.15;
  sigma1              = 0.002;
  sigma2              = 0.002;
  verbose             = false;
  has_wisdom          = false;
  esmooth1            = 0.0;
  esmooth2            = 0.0;
  has_user_seed       = false;
  seed                = 1024;
  oversample          = 1;
  filterscale         = 0.0;
  qfactor             = 0.2;
  matched             = false;
  single_filt         = true;
  sigm                = 0.0; // Use sigma 1 if needed
  sigm1               = 0.0; // Means use sigma1 if used
  sigm2               = 0.0; // Means use sigma2
  sigc                = 0.006;
  sigc1               = 0.0; // Use sigc if not set
  sigc2               = 0.0; // ditto
  sparcity            = 1;
  nbins               = 1000;
  use_binning         = false;
  map_like            = true;
  powerspecfile       = "";
  nbeambins           = 150;
  nfwhm               = 15.0;

  int c;
  int option_index = 0;
  optind = 1; //Rewind 
  while ((c = getopt_long(argc,argv,optstring,long_options,
			    &option_index)) != -1) 
    switch(c) {
    case 'b':
      use_binning = true;
      break;
    case '!':
      esmooth1 = atof(optarg);
      break;
    case '@':
      esmooth2 = atof(optarg);
      break;
    case 'f':
      fftsize = atoi(optarg);
      break;
    case 'F':
      filterscale = atof(optarg);
      break;
    case '1':
      nbins = atoi(optarg);
      break;
    case '5':
      nbeambins = atoi(optarg);
      break;
    case '6':
      nfwhm = atof(optarg);
      break;
    case 'm':
      matched = true;
      break;
    case 'n':
      nsims = atoi(optarg);
      break;
    case 'N':
      map_like = false;
      break;
    case '2':
      nlike = atoi(optarg);
      break;
    case '3':
      n0initrange = atof(optarg);
      break;
    case '4':
      n0rangefrac = atof(optarg);
      break;
    case 'o':
      oversample = atoi(optarg);
      break;
    case 'p':
      powerspecfile = std::string(optarg);
      break;
    case 'q':
      qfactor = atof(optarg);
      break;
    case '#':
      sigma1 = atof(optarg);
      break;
    case '$':
      sigma2 = atof(optarg);
      break;
    case ',':
      sigc = atof(optarg);
      break;
    case '.':
      sigc1 = atof(optarg);
      single_filt = false;
      break;
    case '/':
      sigc2 = atof(optarg);
      single_filt = false;
      break;
    case '7':
      sigm = atof(optarg);
      break;
    case '9':
      sigm1 = atof(optarg);
      single_filt = false;
      break;
    case '0':
      sigm2 = atof(optarg);
      single_filt = false;
      break;
    case '%':
      sparcity = atoi(optarg);
      break;
    case 'S':
      has_user_seed = true;
      seed = static_cast<unsigned long long int>( atoi(optarg) );
      break;
    case 'v':
      verbose = true;
      break;
    case 'w':
      has_wisdom = true;
      wisdom_file = std::string( optarg );
      break;
    }

  if (optind >= argc - 7) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  n0         = atof(argv[optind]);
  modelfile  = std::string(argv[optind + 1]);
  fwhm1      = atof(argv[optind + 2]);
  fwhm2      = atof(argv[optind + 3]);
  pixsize    = atof(argv[optind + 4]);
  n1         = atoi(argv[optind + 5]);
  n2         = atoi(argv[optind + 6]);
  outputfile = std::string(argv[optind + 7]);

  if (sigma1 < 0.0) {
    std::cout << "Invalid instrument noise level, band 1: must be >= 0.0 "
	      << "but is: " << sigma1 << std::endl;
    return 1;
  }
  if (sigma1 < 0.0) {
    std::cout << "Invalid instrument noise level, band 2: must be >= 0.0 "
	      << "but is: " << sigma2 << std::endl;
    return 1;
  }
  if (n0initrange <= 0.0) {
    std::cout << "Invalid n0initrange: must be > 0" << std::endl;
    return 1;
  }
  if (n0initrange >= 1.0) {
    std::cout << "Invalid n0initrange: must be < 1" << std::endl;
    return 1;
  }
  if (fwhm1 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM, band 1: " << fwhm1 << std::endl;
    return 1;
  }
  if (fwhm2 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM, band 2: " << fwhm2 << std::endl;
    return 1;
  }
  if (esmooth1 < 0.0) {
    std::cout << "Invalid (negative) additional smoothing, band 1:" 
	      << esmooth1 << std::endl;
    return 1;
  }
  if (esmooth2 < 0.0) {
    std::cout << "Invalid (negative) additional smoothing, band 2:" 
	      << esmooth2 << std::endl;
    return 1;
  }
  if (use_binning && nbins == 0) {
    std::cout << "Invalid (zero) number of bins" << std::endl;
    return 1;
  }
  if (oversample == 0) {
    std::cout << "Invalid (non-positive) oversampling" << std::endl;
    return 1;
  }
  if (nfwhm <= 0.0) {
    std::cout << "Invalid (non-positive) NFWHM: " << nfwhm << std::endl;
    return 1;
  }
  if (nbeambins == 0) {
    std::cout << "Invalid number of beam bins (0)" << std::endl;
    return 1;
  }
  if (nbeambins > 10000) {
    std::cout << "Invalid number of beam bins -- too large: " 
	      << nbeambins << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (negative) n0" << std::endl;
    return 1;
  }
  if (pixsize <= 0.0) {
    std::cout << "Invalid (non-positve) pixsize" << std::endl;
    return 1;
  }
  if (n1 * n2 == 0) {
    std::cout << "Simulated image has zero size" << std::endl;
    return 1;
  }
  if (sparcity > n1 * n2) {
    std::cout << "Sampling sparcity less frequent than map size" << std::endl;
    return 1;
  }
  if (map_like && n0rangefrac <= 0.0) {
    std::cout << "Invalid n0rangefrac: must be > 0" << std::endl;
    return 1;
  }
  if (map_like && n0rangefrac >= 1.0) {
    std::cout << "Invalid n0rangefrac: must be < 1" << std::endl;
    return 1;
  }
  if (filterscale < 0.0) {
    std::cout << "Invalid (negative) filter scale: " << filterscale << std::endl;
    return 1;
  }
  if (qfactor < 0.0) {
    std::cout << "Invalid (negative) high-pass filter q factor" 
	      << qfactor << std::endl;
    return 1;
  }


  // Set up filtering parameters.  Rather complex
  // if matched filtering.  The idea is to use
  // the same filter in both bands unless the caller
  // has specified one of the single band variables
  //  (sigma_matched1, sigc1, sigma_matched2, sigc2)
  if (matched) {
    if (single_filt) {
      // Same filter for both bands.  Set into band 1 variables
      if (sigm == 0) sigm1 = sigma1; else sigm1 = sigm;
      sigc1 = sigc;
      if (sigc1 <= 0.0) {
	std::cout << "Invalid sigma_confusion for single matched filter: "
		  << sigc1 << std::endl;
	return 1;
      }
      if (sigm1 <= 0.0) {
	std::cout << "Invalid sigma_instrument for single matched filter: "
		  << sigm1 << std::endl;
	return 1;
      }
    } else {
      // Different filters for each band.  More complex
      if (sigm1 == 0) sigm1 = sigma1;
      if (sigm2 == 0) sigm2 = sigma2;
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
      if (sigm1 <= 0.0) {
	std::cout << "Invalid sigma_instrument1 for double matched filter: "
		  << sigm1 << std::endl;
	return 1;
      }
      if (sigm2 <= 0.0) {
	std::cout << "Invalid sigma_instrument2 for double matched filter: "
		  << sigm2 << std::endl;
	return 1;
      }
    }
  }

  // Main execution block
  try {
    double base_n0;
    if (n0 == 0 || verbose) {
      // Use base model for n0
      numberCountsDouble model(modelfile);
      base_n0 = model.getBaseN0();
      if (n0 == 0)
	n0 = base_n0;
    }

    if (verbose) {
      double area = n1*n2*std::pow(pixsize/3600.0,2);
      printf("   base model file:    %s\n", modelfile.c_str());
      printf("   base n0:            %0.3e\n", base_n0);
      printf("   nsims:              %u\n",nsims);
      printf("   n0initrange         %0.3f\n",n0initrange);
      printf("   Beam fwhm1:         %0.2f\n",fwhm1);
      printf("   Beam fwhm2:         %0.2f\n",fwhm2);
      printf("   Nfwhm:              %0.2f\n", nfwhm);
      printf("   N Beam Bins:        %u\n", nbeambins);
      printf("   Pixel size:         %0.1f\n",pixsize);
      printf("   Number of pixels:   %u x %u\n",n1,n2);
      printf("   Area:               %0.2f\n",area);
      printf("   N0:                 %0.3e\n",n0);
      printf("   sigma1:             %0.4f\n",sigma1);
      printf("   sigma2:             %0.4f\n",sigma2);
      printf("   fftsize:            %u by %u\n", fftsize, fftsize);
      if (esmooth1 > 0) 
	printf("   esmooth1:           %0.2f\n",esmooth1);
      if (esmooth2 > 0) 
	printf("   esmooth2:           %0.2f\n",esmooth2);
      if (oversample > 1)
	printf("   oversampling:       %u\n",oversample);
      if (filterscale > 0) {
	printf("   filtering scale:    %0.1f\n", filterscale);
	printf("   filtering q:        %0.2f\n", qfactor);
      }
      if (matched > 0) {
	if (single_filt) {
	  printf("   matched fwhm:       %0.1f\n", fwhm1);
	  printf("   matched sigi:       %0.4f\n", sigm1);
	  printf("   matched sigc:       %0.4f\n", sigc1);
	} else {
	  printf("   matched fwhm1:      %0.1f\n", fwhm1);
	  printf("   matched fwhm2:      %0.1f\n", fwhm2);
	  printf("   matched sigm1:      %0.4f\n", sigm1);
	  printf("   matched sigm2:      %0.4f\n", sigm2);
	  printf("   matched sigc1:      %0.4f\n", sigc1);
	  printf("   matched sigc2:      %0.4f\n", sigc2);
	}
      }
      if (sparcity > 1)
	printf("   sparcity:           %u\n", sparcity);
      if (!powerspecfile.empty())
	printf("   clustering P(k):    %s\n", powerspecfile.c_str());
      if (use_binning)
	printf("   nbins:              %u\n",nbins);
      if (map_like) {
	printf("   nlike:              %u\n",nlike);
	printf("   n0rangefrac         %0.3f\n",n0rangefrac);
      }
    }

    simManagerDouble sim(modelfile, nsims, n0initrange, map_like, nlike, 
			 n0rangefrac, fftsize, n1, n2, pixsize, fwhm1, fwhm2, 
			 nfwhm, sigma1, sigma2, single_filt, filterscale,
			 qfactor, matched, sigm1, sigm2, sigc1, sigc2,
			 nbeambins, n0, esmooth1, esmooth2, oversample, 
			 powerspecfile, sparcity, use_binning, nbins);
    if (has_wisdom) sim.addWisdom(wisdom_file);
    if (has_user_seed) sim.setSeed(seed);

    // Main loop
    sim.doSims(verbose);

    //Write it
    if (verbose) std::cout << "Writing simulation results to " << outputfile 
			   << std::endl;
    utility::outfiletype oft = utility::getOutputFileType(outputfile);
    if (oft == utility::HDF5 || oft == utility::UNKNOWN)
      sim.writeToHDF5(outputfile);
    else if (oft == utility::FITS) {
      int status = sim.writeToFits(outputfile);
      if (status != 0) return status;
    } else if (oft == utility::TXT)
      throw pofdExcept("pofd_coverage_runSim", "runSimSingle",
		       "Output to text not supported", 1);
  } catch (const pofdExcept& ex) {
    std::cout << "Error encountered" << std::endl;
    std::cout << ex << std::endl;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cout << "Bad allocation error: " << ba.what() << std::endl;
    return 16;
  } 

  return 0;
}


////////////////////////////////////////



int main(int argc, char **argv) {
  bool twod;
  
  twod = false;
  
  //Only interested in a) displaying help and b) figuring out
  // if this is 1D or 2D c) displaying the version number
  int c;
  int option_index = 0;
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h':
      std::cout << "NAME" << std::endl;
      std::cout << "\tpofd_coverage_runSim -- make a set of "
		<< "simulated images for" << std::endl;
      std::cout << "\t a spline model with a"
		<< " Gaussian beam (1D) and measure" << std::endl;
      std::cout << "\t the number of objects per sq deg from them."
		<< std::endl;
      std::cout << "\t The 2D case uses the same 1D model paired with a"
		<< " log-normal" << std::endl;
      std::cout << "\t color model in flux_2/flux_1." << std::endl;
      std::cout << std::endl;
      std::cout << "SYNOPSIS" << std::endl;
      std::cout << "\t One-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_runSim [options] n0 modelfile "
		<< "fwhm pixsize n1 n2" << std::endl;
      std::cout << "\t    outputfile" << std::endl;      
      std::cout << std::endl;
      std::cout << "\t Two-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_runSim [options] -d n0 modelfile "
		<< "fwhm1 fwhm2 pixsize" << std::endl;
      std::cout << "\t    n1 n2 outputfile" << std::endl;
      std::cout << "DESCRIPTION" << std::endl;
      std::cout << "\tCreates simulated images for a given model, computes the"
		<< " P(D)" << std::endl;
      std::cout << "\tfor a range of input n0 values, and finds the best"
		<< " fitting value" << std::endl;
      std::cout << "\tof n0 for each image.  Optionally, this is then followed"
		<< " by" << std::endl;
      std::cout << "\tmapping out a range of n0 values around the best fit and"
		<< " storing" << std::endl;
      std::cout << "\tthe resulting log likelihood.  The results are written to"
		<< " outputfile" << std::endl;
      std::cout << "\twith the file type (hdf5 or fits) controlled by the "
		<< "extension." << std::endl;
      std::cout << "\tHDF5 is the default. If n0 is zero, then the value from "
		<< " the" << std::endl;
      std::cout << "\traw model file is adopted." << std::endl;
      std::cout << std::endl;
      std::cout << "\tThe number counts model in 1D is a spline" << std::endl;
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
      std::cout << "\tfwhm is the beam FWHM in arcsec , and the beam is "
		<< "assumed" << std::endl;
      std::cout << "\tGaussian.  In 2D, fwhm1 and fwhm2 are the FWHM values for"
		<< std::endl;
      std::cout << "\tband." << std::endl;
      std::cout << std::endl;
      std::cout << "\tPixsize gives the pixel size (in arcsec), n1 and n2 give"
		<< " the" << std::endl;
      std::cout << "\tnumber of pixels along each dimension in the simulated "
		<< "data." << std::endl;
      std::cout << std::endl;
      std::cout << "OPTIONS" << std::endl;
      std::cout << "\t-h --help" << std::endl;
      std::cout << "\t\tPrint this message and exit." << std::endl;
      std::cout << "\t-b, --bin" << std::endl;
      std::cout << "\t\tBin the simulated image for the likelihood calculation."
		<< std::endl;
      std::cout << "\t-F, --filtscale VALUE" << std::endl;
      std::cout << "\t\tRadius of high-pass filter in arcseconds. If zero,"
		<< std::endl;
      std::cout << "\t\tno filtering is applied (def: 0)." << std::endl;
      std::cout << "\t-f, --fftsize FFTSIZE" << std::endl;
      std::cout << "\t\tSize of FFT to use when computing P(D) (def: 131072 in"
		<< std::endl;
      std::cout << "\t\tone dimension and 4096 in two dimensions.)" 
		<< std::endl;
      std::cout << "\t-m, --matched" << std::endl;
      std::cout << "\t\tApply matched filtering to the beam, with a FWHM"
		<< " matching the" << std::endl;
      std::cout << "\t\tbeam.  In the two-d case, a different filter is "
		<< "applied in" << std::endl;
      std::cout << "\t\teach band matching the properties of that band.  Off"
		<< " by default." << std::endl;
      std::cout << "\t--nbeambins NBINS" << std::endl;
      std::cout << "\t\tNumber of histogram bins to use for beams. (def: 100"
		<< std::endl;
      std::cout << "\t\tin the 1D case, 150 for the 2D case.)" << std::endl;
      std::cout << "\t--nbins NBINS" << std::endl;
      std::cout << "\t\tNumber of bins to use if binning simulated image."
		<< std::endl;
      std::cout << "\t--nfwhm NFWHM" << std::endl;
      std::cout << "\t\tNumber of FWHM kept after filtering. (def: 15.0)"
		<< std::endl;
      std::cout << "\t-n, --nsims NSIMS" << std::endl;
      std::cout << "\t\tThe number of simulations to do (def: 100)." 
		<< std::endl;
      std::cout << "\t-N, --nolike" << std::endl;
      std::cout << "\t\tDo not map out the likelihood values." << std::endl;
      std::cout << "\t--nlike NLIKE" << std::endl;
      std::cout << "\t\tThe number of likelihoods to compute for each sim if" 
		<< std::endl;
      std::cout << "\t\t --nolike is not set (def: 500)." << std::endl;
      std::cout << "\t--n0initrange N0INITRANGE" << std::endl;
      std::cout << "\t\tFractional range used to bracket likelihood in "
		<< "initial" << std::endl;
      std::cout << "\t\tsearch for maximum likelihood (def: 0.15)."
		<< std::endl;
      std::cout << "\t--n0rangefrac N0RANGEFRAC" << std::endl;
      std::cout << "\t\tThe fractional range in n0 to explore in each "
		<< "direction" << std::endl;
      std::cout << "\t\tif doing likelihood map (def: 0.1)."
		<< std::endl;
      std::cout << "\t-o, --oversample VALUE" << std::endl;
      std::cout << "\t\tAmount of oversampling to use (integral) when " 
		<< "generating" << std::endl;
      std::cout << "\t\timage.  The data is then down-binned to the specified"
		<< "size." << std::endl;
      std::cout << "\t\tThe default is to apply no oversampling." << std::endl;
      std::cout << "\t-p, --powerspec FILENAME" << std::endl;
      std::cout << "\t\tName of text file giving k, P(k) (in 1/arcmin and "
		<< "Jy/sr)" << std::endl;
      std::cout << "\t\tfor on-sky source clustering to include in simulation."
		<< " If" << std::endl;
      std::cout << "\t\tnot specified, the sources are uniformly distributed."
		<< std::endl;
      std::cout << "\t-q, --qfactor VALUE" << std::endl;
      std::cout << "\t\tHigh-pass filter apodization sigma as fraction of"
		<< std::endl;
      std::cout << "\t\tfiltscale. (def: 0.2)." << std::endl;
      std::cout << "\t-S, --seed SEED" << std::endl;
      std::cout << "\t\tSet user specified seed, otherwise taken from time."
		<< std::endl;
      std::cout << "\t--sparcity SPARCITY" << std::endl;
      std::cout << "\t\tOnly sample the simulated maps every this many pixels"
		<< std::endl;
      std::cout << "\t\twhen computing the likelihood." << std::endl;
      std::cout << "\t-v, --verbose" << std::endl;
      std::cout << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cout << "\t-V, --version" << std::endl;
      std::cout << "\t\tOutput version number and exit" << std::endl;
      std::cout << "\t-w, --wisdom wisdomfile" << std::endl;
      std::cout << "\t\tName of wisdom file (prepared with fftw-wisdom)." 
		<< std::endl;
      std::cout << "ONE-DIMENSIONAL OPTIONS" << std::endl;
      std::cout << "\t-e, --esmooth ESMOOTH" << std::endl;
      std::cout << "\t\tExtra smoothing FWHM (in arcsec)" << std::endl;
      std::cout << "\t-s, --sigma noise" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise (assumed Gaussian).  This"
		<< " is" << std::endl;
      std::cout << "\t\tthe level before any extra smoothing (def: 0.002)."
		<< std::endl;
      std::cout << "\t--sigc VALUE" << std::endl;
      std::cout << "\t\tConfusion noise for matched filtering, in Jy. (Def:"
		<< " 0.006)" << std::endl;
      std::cout << "\t--sigma_matched VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, in Jy. (Def:"
		<< " sigma)" << std::endl;
      std::cout << "TWO-DIMENSIONAL OPTIONS" << std::endl;
      std::cout << "\t---esmooth1 ESMOOTH" << std::endl;
      std::cout << "\t\tExtra smoothing FWHM (in arcsec), band 1" << std::endl;
      std::cout << "\t---esmooth2 ESMOOTH" << std::endl;
      std::cout << "\t\tExtra smoothing FWHM (in arcsec), band 2" << std::endl;
      std::cout << "\t--sigma1 noise" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise (assumed Gaussian), band1."
		<< " This" << std::endl;
      std::cout << "\t\tis the level before any extra smoothing (def: 0.002)."
		<< std::endl;
      std::cout << "\t--sigma2 noise" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise (assumed Gaussian), band2."
		<< " This" << std::endl;
      std::cout << "\t\tis the level before any extra smoothing (def: 0.002)."
		<< std::endl;
      std::cout << "\t--sigc1 VALUE" << std::endl;
      std::cout << "\t\tConfusion noise for matched filtering, in Jy, band 1."
		<< std::endl;
      std::cout << "\t\t(Def: 0.006)" << std::endl;
      std::cout << "\t--sigc2 VALUE" << std::endl;
      std::cout << "\t\tConfusion noise for matched filtering, in Jy, band 2."
		<< std::endl;
      std::cout << "\t\t(Def: 0.006)" << std::endl;
      std::cout << "\t--sigma_matched1 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 1, in Jy."
		<< std::endl;
      std::cout << "\t\t(Def: sigma1)" << std::endl;
      std::cout << "\t--sigma_matched2 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 2, in Jy."
		<< std::endl;
      std::cout << "\t\t(Def: sigma2)" << std::endl;
      return 0;
      break;
    case 'd':
      twod = true;
      break;
    case 'V':
      std::cout << "pofd_coverage version number: " << pofd_coverage::version 
		<< std::endl;
      return 0;
      break;
    }

  if (!twod)
    return runSimSingle(argc, argv);
  else 
    return runSimDouble(argc, argv);
}
