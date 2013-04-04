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
  {"double",no_argument,0,'d'},
  {"extra_smooth",required_argument,0,'e'},
  {"extra_smooth1",required_argument,0,'1'},
  {"extra_smooth2",required_argument,0,'2'},
  {"help", no_argument, 0, 'h'},
  {"oversample", required_argument, 0, 'o'},
  {"seed", required_argument, 0,'S'},
  {"sigma",required_argument,0,'s'},
  {"sigma1",required_argument,0,'3'},
  {"sigma2",required_argument,0,'4'},
  {"verbose",no_argument,0,'v'},
  {"version",no_argument,0,'V'},
  {0,0,0,0}
};
char optstring[] = "de:1:2:ho:S:s:3:4:vV";


int makeSimSingle(int argc, char **argv) {

  unsigned int n1, n2;
  double n0, pixsize, sigma, fwhm;
  double extra_smooth; //Additional smoothing
  std::string modelfile, outputfile; 
  unsigned long long int user_seed;
  bool verbose, do_extra_smooth, have_user_seed;
  unsigned int oversample;

  //Defaults
  do_extra_smooth     = false;
  extra_smooth        = 0.0;
  sigma               = 0.0;
  verbose             = false;
  user_seed           = 0;
  have_user_seed      = false;
  oversample          = 1;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'e' :
      do_extra_smooth = true;
      extra_smooth = atof(optarg);
      break;
    case 'o':
      oversample = atoi(optarg);
      break;
    case 'S' :
      have_user_seed = true;
      user_seed = static_cast<unsigned long long int>( atoi(optarg) );
      break;
    case 's' :
      sigma = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc-6 ) {
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

  if (n0 <= 0.0) {
    std::cout << "Invalid (non-positive) n0: " << n0 << std::endl;
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
  if (do_extra_smooth && extra_smooth < 0.0) {
    std::cout << "Invalid (non-positive) extra smoothing FWHM" << std::endl;
    return 1;
  }
  if (oversample == 0) {
    std::cout << "Invalid (non-positive) oversampling" << std::endl;
    return 1;
  }
  

  try {
    numberCounts model(modelfile);
    if (verbose)
      std::cout << "Base model n0: " << model.getBaseN0()
		<< " Your value: " << n0 << std::endl;

    simImage dim(n1, n2, pixsize, fwhm, sigma, extra_smooth,
		 oversample);
    if (have_user_seed) dim.setSeed( user_seed );
    dim.realize(model, n0, do_extra_smooth, true, false); //Do mean subtract

    //Write it
    if (verbose) std::cout << "Writing simulated image to " << outputfile 
			   << std::endl;
    int status = dim.writeToFits(outputfile);
    if (status != 0) return status;
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
  double extra_smooth1, extra_smooth2; //Additional smoothing
  std::string modelfile, outputfile; 
  unsigned long long int user_seed;
  bool verbose, do_extra_smooth, have_user_seed;
  unsigned int oversample;

  //Defaults
  do_extra_smooth     = false;
  extra_smooth1       = 0.0;
  extra_smooth2       = 0.0;
  sigma1              = 0.0;
  sigma2              = 0.0;
  verbose             = false;
  user_seed           = 0;
  have_user_seed      = false;
  oversample          = 1;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case '1' :
      do_extra_smooth = true;
      extra_smooth1 = atof(optarg);
      break;
    case '2' :
      do_extra_smooth = true;
      extra_smooth2 = atof(optarg);
      break;
    case 'o':
      oversample = atoi(optarg);
      break;
    case 'S' :
      have_user_seed = true;
      user_seed = static_cast<unsigned long long int>( atoi(optarg) );
      break;
    case '3' :
      sigma1 = atof(optarg);
      break;
    case '4' :
      sigma2 = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc-7 ) {
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

  if (n0 <= 0.0) {
    std::cout << "Invalid (non-positive) n0: " << n0 << std::endl;
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
  if (do_extra_smooth && extra_smooth1 < 0.0) {
    std::cout << "Invalid (non-positive) extra smoothing FWHM band 1" 
	      << std::endl;
    return 1;
  }
  if (do_extra_smooth && extra_smooth2 < 0.0) {
    std::cout << "Invalid (non-positive) extra smoothing FWHM band 2" 
	      << std::endl;
    return 1;
  }
  if (oversample == 0) {
    std::cout << "Invalid (non-positive) oversampling" << std::endl;
    return 1;
  }
  

  try {
    numberCountsDouble model(modelfile);
    if (verbose)
      std::cout << "Base model n0: " << model.getBaseN0()
		<< " Your value: " << n0 << std::endl;

    simImageDouble dim(n1, n2, pixsize, fwhm1, fwhm2, sigma1, sigma2, 
		       extra_smooth1, extra_smooth2, oversample);
    if (have_user_seed) dim.setSeed( user_seed );
    dim.realize(model, n0, do_extra_smooth, true, false); //Do mean subtract

    //Write it
    if (verbose) std::cout << "Writing simulated image to " << outputfile 
			   << std::endl;
    int status = dim.writeToFits(outputfile);
    if (status != 0) return status;
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
    case 'h' :
      std::cout << "NAME" << std::endl;
      std::cout << "\tpofd_coverage_makeSim -- make simulated images for"
		<< " a broken" << std::endl;
      std::cout << "\t power law type model with Gaussian beams (1D) or the"
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
      std::cout << "\tto outfile.  The number counts model in 1D is a broken "
		<< "power" << std::endl;
      std::cout << "\ta law model specified by modelfile, and by the number of"
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
      std::cout << "\tfwhm is the beam FWHM in arcsec (1D), and fwhm1, fwhm2"
		<< " the values" << std::endl;
      std::cout << "\tin each band in the 2D case.  The beams are assumed "
		<< "Gaussian. " << std::endl;
      std::cout << "\tPixsize gives the pixel size (in arcsec), while n1 and "
		<< "n2 are the" << std::endl;
      std::cout << "\tnumber of pixels along each dimension." << std::endl;
      std::cout << std::endl;
      std::cout << "OPTIONS" << std::endl;
      std::cout << "\t-h --help" << std::endl;
      std::cout << "\t\tPrint this message and exit." << std::endl;
      std::cout << "\t-d, --double" << std::endl;
      std::cout << "\t\tUse the 2D model instead of the 1D one." << std::endl;
      std::cout << "\t-o, --oversample VALUE" << std::endl;
      std::cout << "\t\tAmount of oversampling to use (integral) when " 
		<< "generating" << std::endl;
      std::cout << "\t\timage.  The data is then down-binned to the specified"
		<< "size." << std::endl;
      std::cout << "\t\tThe default is to apply no oversampling." << std::endl;
      std::cout << "\t--S, --seed SEED" << std::endl;
      std::cout << "\t\tUse this seed for the random number generator." 
		<< std::endl;
      std::cout << "\t-v, --verbose" << std::endl;
      std::cout << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cout << "\t-V, --version" << std::endl;
      std::cout << "\t\tOutput version number and exit" << std::endl;
      std::cout << "\tONE-D MODEL OPTIONS" << std::endl;
      std::cout << "\t-e, --extra_smooth FWHM" << std::endl;
      std::cout << "\t\tApply additional smoothing with a Gaussian of this"
		<< " FWHM" << std::endl;
      std::cout << "\t\t(in arcseconds); this is applied after noise is added."
		<< std::endl;
      std::cout << "\t--sigma NOISE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise (def: 0)." << std::endl;
      std::cout << "\tTWO-D MODEL OPTIONS" << std::endl;
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
    return makeSimSingle(argc,argv);
  else 
    return makeSimDouble(argc,argv);

}



