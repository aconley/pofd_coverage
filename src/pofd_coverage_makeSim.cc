#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cstdlib>

#include<getopt.h>

#include "../include/pofdExcept.h"
#include "../include/global_settings.h"
#include "../include/simImage.h"
#include "../include/simImageDouble.h"

//Set up global option index that can be used for both single and double case
static struct option long_options[] = {
  {"double", no_argument, 0, 'd'},
  {"extra_smooth", required_argument, 0, 'e'},
  {"extra_smooth1", required_argument, 0, '1'},
  {"extra_smooth2", required_argument, 0, '2'},
  {"filtscale", required_argument, 0, 'F'},
  {"help", no_argument, 0, 'h'},
  {"matched", no_argument, 0, 'm'},
  {"nfinaltrials", required_argument, 0, 'n'},
  {"oversample", required_argument, 0, 'o'},
  {"powerspec", required_argument, 0, 'p'},
  {"probimfile", required_argument, 0, 'P'},
  {"qfactor", required_argument, 0, 'q'},
  {"seed", required_argument, 0, 'S'},
  {"sigma", required_argument, 0, 's'},
  {"sigma1", required_argument, 0, '3'},
  {"sigma2", required_argument, 0, '4'},
  {"sigc", required_argument, 0, '6'},
  {"sigc1", required_argument, 0, '.'},
  {"sigc2", required_argument, 0, '/'},
  {"sigma_matched", required_argument, 0, '5'},
  {"sigma_matched1", required_argument, 0, '7'},
  {"sigma_matched2", required_argument, 0, '8'},
  {"verbose", no_argument, 0, 'v'},
  {"version", no_argument, 0, 'V'},
  {0,0,0,0}
};
char optstring[] = "de:1:2:F:hmn:o:p:P:q:S:s:3:4:5:6:.:/:7:8:vV";


int makeSimSingle(int argc, char **argv) {

  unsigned int n1, n2;
  double n0, pixsize, sigma, fwhm;
  double filterscale, qfactor, sigi, sigc; // Filtering params
  double extra_smooth; //Additional smoothing
  std::string modelfile, outputfile, powspecfile, probimfile;
  unsigned long long int user_seed;
  bool verbose, have_user_seed, matched;
  unsigned int oversample, nfinal;

  //Defaults
  extra_smooth        = 0.0;
  sigma               = 0.0;
  filterscale         = 0.0;
  qfactor             = 0.2;
  matched             = false;
  sigi                = 0.0; // Means: use sigma
  sigc                = 0.006;
  verbose             = false;
  user_seed           = 0;
  have_user_seed      = false;
  oversample          = 1;
  powspecfile         = "";
  probimfile          = "";
  nfinal              = 0;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc, argv, optstring, long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'e':
      extra_smooth = atof(optarg);
      break;
    case 'F':
      filterscale = atof(optarg);
      break;
    case 'm':
      matched = true;
      break;
    case 'n':
      nfinal = atoi(optarg);
      break;
    case 'o':
      oversample = atoi(optarg);
      break;
    case 'p':
      powspecfile = std::string(optarg);
      break;
    case 'P':
      probimfile = std::string(optarg);
      break;
    case 'q':
      qfactor = atof(optarg);
      break;
    case 'S':
      have_user_seed = true;
      user_seed = static_cast<unsigned long long int>( atoi(optarg) );
      break;
    case 's':
      sigma = atof(optarg);
      break;
    case '5':
      sigi = atof(optarg);
      break;
    case '6':
      sigc = atof(optarg);
      break;
    case 'v':
      verbose = true;
      break;
    }

  if (optind >= argc - 6) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string(argv[optind]);
  n0         = atof(argv[optind + 1]);
  fwhm       = atof(argv[optind + 2]);
  pixsize    = atof(argv[optind + 3]);
  n1         = atoi(argv[optind + 4]);
  n2         = atoi(argv[optind + 5]);
  outputfile = std::string(argv[optind + 6]);

  if (matched && (sigi == 0)) sigi = sigma;

  if (n0 < 0.0) {
    std::cout << "Invalid (negative) n0: " << n0 << std::endl;
    return 1;
  }
  if (std::isnan(n0) || std::isinf(n0)) {
    std::cout << "Invalid (non-finite) n0: " << n0 << std::endl;
    return 1;
  }
  if (sigma < 0.0) {
    std::cout << "Invalid noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (fwhm < 0.0) {
    std::cout << "Invalid (non-positive) FWHM" << std::endl;
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
  if (extra_smooth < 0.0) {
    std::cout << "Invalid (non-positive) extra smoothing FWHM" << std::endl;
    return 1;
  }
  if (oversample == 0) {
    std::cout << "Invalid (non-positive) oversampling" << std::endl;
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
  if ((!probimfile.empty()) && powspecfile.empty()) {
    std::cout << "Specified output of probim, but didn't provide clustering."
	      << std::endl;
    return 1;
  }

  try {
    numberCounts model(modelfile);
    if (n0 == 0) {
      n0 = model.getBaseN0();
      if (verbose)
	std::cout << "Using model n0: " << n0 << std::endl;
    } else if (verbose)
      std::cout << "Base model n0: " << model.getBaseN0()
		<< " Your value: " << n0 << std::endl;

    fourierFilter *filt = NULL;
    if (filterscale > 0) {
      if (matched) {
	filt = new fourierFilter(pixsize, fwhm, sigi, sigc,
				 filterscale, qfactor, true);
      } else
	filt = new fourierFilter(pixsize, filterscale, qfactor, true);
    } else if (matched)
	filt = new fourierFilter(pixsize, fwhm, sigi, sigc, true);

    simImage im(n1, n2, pixsize, fwhm, sigma, extra_smooth,
		oversample, 1000, powspecfile);
    if (have_user_seed) im.setSeed(user_seed);

    // Get final noise if specified
    if (verbose && (nfinal > 0) && (filt != NULL))
      std::cout << "Estimating final (filtered) noise with "
		<< nfinal << " trials" << std::endl;
    if (nfinal != 0) im.getFinalNoise(nfinal, filt);

    // Generate with mean subtraction
    im.realize(model, n0, true, filt, false); 

    if (filt != NULL) delete filt;

    // Write it
    if (verbose) std::cout << "Writing simulated image to " << outputfile
			   << std::endl;
    int status = im.writeToFits(outputfile);
    if (status != 0) return status;

    // Write prob image if needed
    if (!probimfile.empty()) {
      if (verbose)
	std::cout << "Outputting probability map to: " << probimfile
		  << std::endl;
      status = im.writeProbImageToFits(probimfile);
      if (status != 0) return status;
    }

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

///////////////////////////////////

int makeSimDouble(int argc, char **argv) {

  unsigned int n1, n2;
  double n0, pixsize, sigma1, sigma2, fwhm1, fwhm2;
  // Filtering params
  bool matched, single_filt;
  double filterscale, qfactor; // Hipass
  double sigm1, sigm2, sigc1, sigc2, sigm, sigc; // Matched

  double extra_smooth1, extra_smooth2; //Additional smoothing
  std::string modelfile, outputfile, powspecfile, probimfile; 
  unsigned long long int user_seed;
  bool verbose, have_user_seed;
  unsigned int oversample, nfinal;

  //Defaults
  extra_smooth1       = 0.0;
  extra_smooth2       = 0.0;
  sigma1              = 0.0;
  sigma2              = 0.0;
  filterscale         = 0.0;
  qfactor             = 0.2;
  matched             = false;
  single_filt         = true;
  sigc                = 0.006; // Value if using single filt
  sigc1               = 0.0; // Value if multiple filts; 0 means don't use
  sigc2               = 0.0; // Value if multiple filts; 0 means don't use
  sigm                = 0; // Value if using single filt; if 0, use sigi1
  sigm1               = 0.0; // Value if multiple filts; 0 means use sigma1 if mult
  sigm2               = 0.0; // Value if multiple filts; 0 means use sigma2 if mult
  verbose             = false;
  user_seed           = 0;
  have_user_seed      = false;
  oversample          = 1;
  powspecfile         = "";
  probimfile          = "";
  nfinal              = 0;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc, argv, optstring, long_options,
			  &option_index)) != -1) 
    switch(c) {
    case '1':
      extra_smooth1 = atof(optarg);
      break;
    case '2':
      extra_smooth2 = atof(optarg);
      break;
    case 'F':
      filterscale = atof(optarg);
      break;
    case 'm':
      matched = true;
      break;
    case 'n':
      nfinal = atoi(optarg);
      break;
    case 'o':
      oversample = atoi(optarg);
      break;
    case 'p':
      powspecfile = std::string(optarg);
      break;
    case 'P':
      probimfile = std::string(optarg);
      break;
    case 'q':
      qfactor = atof(optarg);
      break;
    case 'S':
      have_user_seed = true;
      user_seed = static_cast<unsigned long long int>( atoi(optarg) );
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
    case '.':
      sigc1 = atof(optarg);
      single_filt = false;
      break;
    case '/':
      sigc2 = atof(optarg);
      single_filt = false;
      break;
    case '5':
      sigm = atof(optarg);
      break;
    case '7':
      sigm1 = atof(optarg);
      single_filt = false;
      break;
    case '8':
      sigm2 = atof(optarg);
      single_filt = false;
      break;
    case 'v':
      verbose = true;
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
  n1         = atoi(argv[optind + 5]);
  n2         = atoi(argv[optind + 6]);
  outputfile = std::string(argv[optind + 7]);


  if (n0 < 0.0) {
    std::cout << "Invalid (negative) n0: " << n0 << std::endl;
    return 1;
  }
  if (std::isnan(n0) || std::isinf(n0)) {
    std::cout << "Invalid (non-finite) n0: " << n0 << std::endl;
    return 1;
  }
  if (sigma1 < 0.0) {
    std::cout << "Invalid band 1 noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (sigma2 < 0.0) {
    std::cout << "Invalid band 2 noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (fwhm1 < 0.0) {
    std::cout << "Invalid (non-positive) FWHM band 1" << std::endl;
    return 1;
  }
  if (fwhm2 < 0.0) {
    std::cout << "Invalid (non-positive) FWHM band 1" << std::endl;
    return 1;
  }
  if (extra_smooth1 < 0.0) {
    std::cout << "Invalid (non-positive) extra smoothing FWHM band 1" 
	      << std::endl;
    return 1;
  }
  if (extra_smooth2 < 0.0) {
    std::cout << "Invalid (non-positive) extra smoothing FWHM band 2" 
	      << std::endl;
    return 1;
  }
  if (oversample == 0) {
    std::cout << "Invalid (non-positive) oversampling" << std::endl;
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
  if ((!probimfile.empty()) && powspecfile.empty()) {
    std::cout << "Specified output of probim, but didn't provide clustering."
	      << std::endl;
    return 1;
  }

  // Set up filtering parameters.  Rather complex
  // if matched filtering.  The idea is to use
  // the same filter in both bands unless the caller
  // has specified one of the single band variables
  //  (sigma_matched1, sigc1, sigma_matched2, sigc2)
  if (matched) {
    if (single_filt) {
      // Same filter for both bands
      if (sigm == 0) sigm = sigma1;
      if (sigc <= 0.0) {
	std::cout << "Invalid sigma_confusion for single matched filter: "
		  << sigc << std::endl;
	return 1;
      }
      if (sigm <= 0.0) {
	std::cout << "Invalid sigma_instrument for single matched filter: "
		  << sigm << std::endl;
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
    numberCountsDouble model(modelfile);
    if (n0 == 0) {
      n0 = model.getBaseN0();
      if (verbose)
	std::cout << "Using model n0: " << n0 << std::endl;
    } else if (verbose)
      std::cout << "Base model n0: " << model.getBaseN0()
		<< " Your value: " << n0 << std::endl;

    // Set up filtering
    fourierFilter *filt1 = NULL, *filt2 = NULL;
    if (filterscale > 0) {
      if (matched) {
	if (single_filt) {
	  if (verbose)
	    std::cout << "Using single matched + highpass filter" 
		      << std::endl;
	  filt1 = new fourierFilter(pixsize, fwhm1, sigm, sigc,
				    filterscale, qfactor, true);
	} else {
	  if (verbose)
	    std::cout << "Using double matched + highpass filter" 
		      << std::endl;
	  filt1 = new fourierFilter(pixsize, fwhm1, sigm1, sigc1,
				    filterscale, qfactor, true);
	  filt2 = new fourierFilter(pixsize, fwhm2, sigm2, sigc2,
				    filterscale, qfactor, true);
	} 
      } else {
	if (verbose)
	  std::cout << "Using single highpass filter" << std::endl;
	filt1 = new fourierFilter(pixsize, filterscale, qfactor, true);
      }
    } else if (matched) {
      if (single_filt) {
	if (verbose)
	    std::cout << "Using single matched filter" << std::endl;
	filt1 = new fourierFilter(pixsize, fwhm1, sigm, sigc, true);
      } else {
	if (verbose)
	    std::cout << "Using double matched filter" << std::endl;
	filt1 = new fourierFilter(pixsize, fwhm1, sigm1, sigc1, true);
	filt2 = new fourierFilter(pixsize, fwhm2, sigm2, sigc2, true);
      }
    }

    simImageDouble dim(n1, n2, pixsize, fwhm1, fwhm2, sigma1, sigma2, 
		       extra_smooth1, extra_smooth2, oversample, 
		       1000, powspecfile);
    if (have_user_seed) dim.setSeed(user_seed);
    
    // Get final noise if specified
    if (verbose && (nfinal > 0) && (filt1 != NULL || filt2 != NULL))
      std::cout << "Estimating final (filtered) noise with "
		<< nfinal << " trials" << std::endl;
    if (nfinal != 0) dim.getFinalNoise(nfinal, filt1, filt2);

    // Generate with mean subtraction
    dim.realize(model, n0, true, filt1, filt2, false);

    if (filt1 != NULL) delete filt1;
    if (filt2 != NULL) delete filt2;

    //Write it
    if (verbose) std::cout << "Writing simulated images to " << outputfile
			   << std::endl;
    int status = dim.writeToFits(outputfile);
    if (status != 0) return status;

    // Write prob image if needed
    if (!probimfile.empty()) {
      if (verbose)
	std::cout << "Outputting probability map to: " << probimfile
		  << std::endl;
      status = dim.writeProbImageToFits(probimfile);
      if (status != 0) return status;
    }

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

///////////////////////////////////

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
    case 'h':
      std::cout << "NAME" << std::endl;
      std::cout << "\tpofd_coverage_makeSim -- make simulated images for"
		<< " a" << std::endl;
      std::cout << "\t spline type model with Gaussian beams (1D) or the"
		<< std::endl;
      std::cout << "\t same type of model paired with a log-normal color model"
		<< std::endl;
      std::cout << "\t in flux_2 / flux_1." << std::endl;
      std::cout << std::endl;
      std::cout << "SYNOPSIS" << std::endl;
      std::cout << "\t One-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_makeSim [options] modelfile n0 fwhm "
		<< "pixsize" << std::endl;
      std::cout << "\t    n1 n2 outputfile" << std::endl;
      std::cout << std::endl;
      std::cout << "\t Two-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_makeSim [options] -d modelfile n0 fwhm1 "
		<< "fwhm2" << std::endl;
      std::cout << "\t    pixsize n1 n2 outputfile" << std::endl;
      std::cout << std::endl;
      std::cout << "DESCRIPTION" << std::endl;
      std::cout << "\tCreates a simulated image for a given model, and writes"
		<< " them" << std::endl;
      std::cout << "\tto outfile.  The number counts model in 1D is a spline"
		<< std::endl;
      std::cout << "\ta model specified by modelfile, and by the number of"
		<< std::endl;
      std::cout << "\tsources per unit area n0.  If you set n0 to zero, then" 
		<< std::endl;
      std::cout << "\tthe raw model from modelfile is used without adjustment."
		<< std::endl;
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
      std::cout << "\tfwhm is the beam FWHM in arcsec (1D), and fwhm1, fwhm2"
		<< " the values" << std::endl;
      std::cout << "\tin each band in the 2D case.  The beams are assumed "
		<< "Gaussian. " << std::endl;
      std::cout << "\tPixsize gives the pixel size (in arcsec), while n1 and "
		<< "n2 are the" << std::endl;
      std::cout << "\tnumber of pixels along each dimension." << std::endl;
      std::cout << std::endl;
      std::cout << "\tThe output file format is FITS.  For the 2D case, the"
		<< " two" << std::endl;
      std::cout << "\tbands are written to the 0th and 1st extension,"
		<< " respectively." << std::endl;
      std::cout << std::endl;
      std::cout << "\tA variety of filtering options are supported: highpass "
		<< " filtering" << std::endl;
      std::cout << "\twith an apodized edge, matched filtering, and the " 
		<< "combination" << std::endl;
      std::cout << "\tof both.  This is relatively simple in the one-band case."
		<< " In" << std::endl;
      std::cout << "\tthe two band case, the default is to use the same filter"
		<< " in" << std::endl;
      std::cout << "\tboth bands for simplicity.  With matched filtering, this"
		<< " is" << std::endl;
      std::cout << "\tbased on the band 1 beam.  To force different filters in"
		<< " each" << std::endl;
      std::cout << "\tband (which uses the FWHM for each band), the user should"
		<< " specify" << std::endl;
      std::cout << "\tone or more of --sigc1, --sigc2, --sigma_matched1, "
		<< "--sigma_matched2," << std::endl;
      std::cout << "\tand if this isn't desired stick to --sigc, --sigma_matched."
		<< std::endl;
      std::cout << std::endl;
      std::cout << "OPTIONS" << std::endl;
      std::cout << "\t-h --help" << std::endl;
      std::cout << "\t\tPrint this message and exit." << std::endl;
      std::cout << "\t-d, --double" << std::endl;
      std::cout << "\t\tUse the 2D model instead of the 1D one." << std::endl;
      std::cout << "\t-F, --filtscale VALUE" << std::endl;
      std::cout << "\t\tRadius of high-pass filter in arcseconds. If zero,"
		<< std::endl;
      std::cout << "\t\tno filtering is applied (def: 0)." << std::endl;
      std::cout << "\t-m, --matched" << std::endl;
      std::cout << "\t\tApply matched filtering to the beam, with a FWHM"
		<< " matching the" << std::endl;
      std::cout << "\t\tbeam (in each band for the 2D case).  Off by default." 
		<< std::endl;
      std::cout << "\t-n, --nfinaltrials VALUE" << std::endl;
      std::cout << "\t\tThe number of trials to do when estimating the final"
		<< " (filtered)" << std::endl;
      std::cout << "\t\tinstrument noise value.  This is only of interest if"
		<< " filtering" << std::endl;
      std::cout << "\t\tis being applied -- if not, the final noise are just "
		<< "the" << std::endl;
      std::cout << "\t\tinput instrument noise values.  A value of 0 (the" 
		<< " default)" << std::endl;
      std::cout << "\t\tmeans that the final noise values are not computed."
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
      std::cout << "\t-P, --probimfile VALUE" << std::endl;
      std::cerr << "\t\tWrite probability image (for clustering) to this file"
		<< std::endl;
      std::cerr << "\t\tin FITS format.  Must have specified --powerspec."
		<< std::endl;
      std::cout << "\t-q, --qfactor VALUE" << std::endl;
      std::cout << "\t\tHigh-pass filter apodization sigma as fraction of"
		<< std::endl;
      std::cout << "\t\tfiltscale. (def: 0.2)." << std::endl;
      std::cout << "\t-S, --seed SEED" << std::endl;
      std::cout << "\t\tUse this seed for the random number generator." 
		<< std::endl;
      std::cout << "\t--sigc VALUE" << std::endl;
      std::cout << "\t\tConfusion noise for matched filtering, in Jy.  This"
		<< std::endl;
      std::cout << "\t\tis the value applied if the same filter is applied in"
		<< " both" << std::endl;
      std::cout << "\t\tbands in the 2D case, but can be overridden on a per-band"
		<< std::endl;
      std::cout << "\t\tbasis by sigc1, sigc2 (Def: 0.006)." << std::endl;
      std::cout << "\t--sigma_matched VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, in Jy.  This"
		<< std::endl;
      std::cout << "\t\tis the value applied if the same filter is applied in"
		<< "both" << std::endl;
      std::cout << "\t\tbands in the 2D case, but can be overridden on a per-band"
		<< std::endl;
      std::cout << "\t\tbasis by sigma_matched1, sigma_mached2.  If not set, the "
		<< "actual" << std::endl;
      std::cout << "\t\tinstrument noise is used (in band 1 for the 2D case)." 
		<< std::endl;
      std::cout << "\t-v, --verbose" << std::endl;
      std::cout << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cout << "\t-V, --version" << std::endl;
      std::cout << "\t\tOutput version number and exit" << std::endl;
      std::cout << "ONE-D MODEL OPTIONS" << std::endl;
      std::cout << "\t-e, --extra_smooth FWHM" << std::endl;
      std::cout << "\t\tApply additional smoothing with a Gaussian of this"
		<< " FWHM" << std::endl;
      std::cout << "\t\t(in arcseconds); this is applied after noise is added."
		<< std::endl;
      std::cout << "\t-s, --sigma NOISE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise (def: 0)." << std::endl;
      std::cout << "TWO-D MODEL OPTIONS" << std::endl;
      std::cout << "\t--extra_smooth1 FWHM" << std::endl;
      std::cout << "\t\tApply additional smoothing in band 1 with a Gaussian of"
		<< std::endl;
      std::cout << "\t\tthis FWHM (in arcseconds); this is applied after noise "
		<< "is added."<< std::endl;
      std::cout << "\t--extra_smooth2 FWHM" << std::endl;
      std::cout << "\t\tApply additional smoothing in band 2 with a Gaussian of"
		<< std::endl;
      std::cout << "\t\tthis FWHM (in arcseconds); this is applied after noise "
		<< "is added."<< std::endl;
      std::cout << "\t--sigma1 NOISE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise, band 1 (def: 0)." 
		<< std::endl;
      std::cout << "\t--sigma2 NOISE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise, band 2 (def: 0)." 
		<< std::endl;
      std::cout << "\t--sigc1 VALUE" << std::endl;
      std::cout << "\t\tConfusion noise for matched filtering, band 1, in Jy."
		<< std::endl;
      std::cout << "\t\tIf this is set, then different filters are used in each"
		<< " band." << std::endl;
      std::cout << "\t--sigc2 VALUE" << std::endl;
      std::cout << "\t\tConfusion noise for matched filtering, band 2, in Jy."
		<< std::endl;
      std::cout << "\t\tIf this is set, then different filters are used in each"
		<< " band." << std::endl;
      std::cout << "\t\t(Def: The instrument noise in band 1)." << std::endl;
      std::cout << "\t--sigma_matched1 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 1, in Jy."
		<< std::endl;
      std::cout << "\t\tIf this is set, then different filters are used in each"
		<< " band." << std::endl;
      std::cout << "\t\tDefaults to sigma1 if a value is required." << std::endl;
      std::cout << "\t--sigma_matched2 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 2, in Jy."
		<< std::endl;
      std::cout << "\t\tIf this is set, then different filters are used in each"
		<< " band." << std::endl;
      std::cout << "\t\tDefaults to sigma2 if a value is required." << std::endl;
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
    return makeSimSingle(argc,argv);
  else 
    return makeSimDouble(argc,argv);

}



